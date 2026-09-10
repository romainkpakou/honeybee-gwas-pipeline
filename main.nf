#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    honeybee-gwas-pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Auteur      : Romain KPAKOU
    Description : Pipeline WGS et GWAS pour Apis mellifera mellifera
    Version     : 1.0.0
    GitHub      : https://github.com/romainkpakou/honeybee-gwas-pipeline

    ÉTAPES :
    1.  QC reads          FastQC · fastp · MultiQC
    2.  Alignement        BWA-MEM2 · SAMtools
    3.  BAM               Picard MarkDuplicates
    4.  Variant calling   GATK HaplotypeCaller (GVCF)
    5.  Génotypage joint  GATK GenomicsDBImport · GenotypeGVCFs
    6.  Filtrage          GATK VariantFiltration · bcftools
    7.  Annotation        SnpEff (base construite depuis le GFF3)
    8.  Génétique pop.    PLINK2 · ADMIXTURE · vcftools
    9.  GWAS              GEMMA LMM · PLINK2
    10. Visualisation     R (Manhattan · QQ · PCA · Admixture · LD · FST)
    11. Rapport           R Markdown HTML

    NOTE DSL2 v26 :
    En Nextflow DSL2 version 26+, TOUT le code exécutable doit être
    à l'intérieur d'un bloc workflow, process ou function.
    Les statements globaux (if, log.info, appels de fonction,
    workflow.onComplete) sont INTERDITS au niveau du script.
    Seuls les blocs suivants sont autorisés au niveau global :
      - nextflow.enable.dsl = 2
      - include { ... } from '...'
      - def maFonction() { ... }
      - workflow { ... }
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// ── Import des modules ────────────────────────────────────────────────────────
// Les includes sont les seuls statements autorisés au niveau global en DSL2
include { FASTQC                   } from './modules/fastqc'
include { FASTP                    } from './modules/fastp'
include { MULTIQC as MULTIQC_QC    } from './modules/multiqc'
include { MULTIQC as MULTIQC_FINAL } from './modules/multiqc'
include { BWA_MEM2_INDEX           } from './modules/bwa_mem2'
include { BWA_MEM2_ALIGN           } from './modules/bwa_mem2'
include { SAMTOOLS_SORT            } from './modules/samtools'
include { SAMTOOLS_INDEX           } from './modules/samtools'
include { SAMTOOLS_FLAGSTAT        } from './modules/samtools'
include { PICARD_MARKDUPLICATES    } from './modules/picard'
include { GATK_HAPLOTYPECALLER     } from './modules/gatk'
include { GATK_DICT                } from './modules/gatk'
include { SAMTOOLS_FAIDX           } from './modules/samtools'
include { GATK_GENOMICSDBIMPORT    } from './modules/gatk'
include { GATK_GENOTYPEGVCFS       } from './modules/gatk'
include { GATK_VARIANTFILTRATION   } from './modules/gatk'
include { BCFTOOLS_FILTER          } from './modules/bcftools'
include { BCFTOOLS_STATS           } from './modules/bcftools'
include { SNPEFF_BUILD             } from './modules/snpeff'
include { SNPEFF_ANNOTATE          } from './modules/snpeff'
include { SNPEFF_COMPRESS          } from './modules/snpeff'
include { PLINK2_QC                } from './modules/plink2'
include { PLINK2_PCA               } from './modules/plink2'
include { PLINK2_GWAS              } from './modules/plink2'
include { BUILD_PHENOTYPE          } from './modules/phenotype'
include { ADMIXTURE_RUN            } from './modules/admixture'
include { VCFTOOLS_FST             } from './modules/vcftools'
include { VCFTOOLS_LD              } from './modules/vcftools'
include { VCFTOOLS_PI              } from './modules/vcftools'
include { GEMMA_KINSHIP            } from './modules/gemma'
include { GEMMA_LMM                } from './modules/gemma'
include { PLOT_MANHATTAN           } from './modules/r_plots'
include { PLOT_QQ                  } from './modules/r_plots'
include { PLOT_PCA                 } from './modules/r_plots'
include { PLOT_ADMIXTURE           } from './modules/r_plots'
include { PLOT_LD_DECAY            } from './modules/r_plots'
include { PLOT_FST                 } from './modules/r_plots'
include { GWAS_REPORT              } from './modules/report'

// ── Fonctions utilitaires ─────────────────────────────────────────────────────
// Les définitions de fonctions sont autorisées au niveau global en DSL2

// Parse le samplesheet CSV et crée un canal Nextflow
// Chaque élément = [meta, [fastq_1, fastq_2]]
// meta = map Groovy : {id, population, sex}
def parseSamplesheet(csv) {
    channel
        .fromPath(csv)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [
                id         : row.sample,
                population : row.population ?: 'unknown',
                sex        : row.sex        ?: 'unknown'
            ]
            def fq1 = file(row.fastq_1, checkIfExists: true)
            def fq2 = row.fastq_2
                ? file(row.fastq_2, checkIfExists: true)
                : null
            // Retourne [meta, [fq1, fq2]] si paired-end
            // ou [meta, [fq1]] si single-end
            fq2 ? [meta, [fq1, fq2]] : [meta, [fq1]]
        }
}

// Détecte si des phénotypes sont disponibles pour le GWAS.
// Deux sources possibles :
//   1. --phenotype_file <fichier dédié>
//   2. une colonne 'phenotype' dans le samplesheet avec >=1 valeur non vide
// Lecture synchrone du CSV (exécutée dans le workflow, sur le nœud principal).
def hasPhenotypes(csv) {
    if (params.phenotype_file) return true
    def f = file(csv)
    if (!f.exists()) return false
    def lines = f.readLines().findAll { row -> row?.trim() }
    if (lines.size() < 2) return false
    def header = lines[0].split(',').collect { col -> col.trim().toLowerCase() }
    def idx = header.indexOf('phenotype')
    if (idx < 0) return false
    return lines.drop(1).any { row ->
        def cols = row.split(',', -1)
        idx < cols.size() &&
            !(cols[idx].trim().toLowerCase() in ['', 'na', 'nan', '-9', 'none', 'null', '.'])
    }
}

// Valide les paramètres obligatoires
// Appelée DANS le workflow — pas au niveau global
def validateParams() {
    if (!params.input) {
        error "ERREUR : --input est obligatoire.\nExemple : --input samplesheet.csv"
    }
    if (!params.genome) {
        error "ERREUR : --genome est obligatoire.\nExemple : --genome data/reference/Amel_HAv3.1.fa"
    }
    if (!file(params.input).exists()) {
        error "ERREUR : Samplesheet introuvable : ${params.input}"
    }
    if (!file(params.genome).exists()) {
        error "ERREUR : Génome introuvable : ${params.genome}"
    }
}

// ── Workflow principal ────────────────────────────────────────────────────────
workflow {

    // ── Validation et message de démarrage ────────────────────────────────────
    // En DSL2 v26, TOUT statement exécutable doit être ici
    validateParams()

    // GWAS activé si un fichier phénotype OU une colonne 'phenotype' est fourni
    run_gwas = hasPhenotypes(params.input)

    // Annotation SnpEff activée si un GFF3 valide est fourni
    run_snpeff = params.gff ? file(params.gff).exists() : false

    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║        honeybee-gwas-pipeline v1.0.0
    ║   WGS & GWAS — Apis mellifera mellifera
    ╚══════════════════════════════════════════════════════════════════╝
      Samplesheet  : ${params.input}
      Génome       : ${params.genome}
      Phénotypes   : ${params.phenotype_file ?: (run_gwas ? "colonne 'phenotype' du samplesheet" : 'non fourni — GWAS ignoré')}
      GWAS         : ${run_gwas ? 'ACTIVÉ' : 'désactivé'}
      Annotation   : ${run_snpeff ? "SnpEff (${params.gff})" : 'désactivée (pas de GFF)'}
      Sortie       : ${params.outdir}
      Modèle GWAS  : ${params.gwas_model}
      MAF          : ${params.maf}
      ADMIXTURE K  : ${params.admixture_k}
      Profil       : ${workflow.profile}
    """.stripIndent()

    // ── Canaux d'entrée ───────────────────────────────────────────────────────

    // Canal des reads : un élément [meta, reads] par échantillon
    ch_reads = parseSamplesheet(params.input)

    // Canal du génome : valeur unique partagée par tous les process
    ch_genome = channel.value(file(params.genome))

    // Indexer le génome pour GATK (fai + dict requis)
    SAMTOOLS_FAIDX(ch_genome)
    GATK_DICT(ch_genome)
    ch_fai  = SAMTOOLS_FAIDX.out.fai
    ch_dict = GATK_DICT.out.dict

    // Canal de la source des phénotypes : fichier dédié si fourni, sinon le
    // samplesheet lui-même (colonne 'phenotype'). BUILD_PHENOTYPE convertit
    // cette source aux formats attendus par GEMMA et PLINK2.
    ch_pheno_source = channel.value(file(params.phenotype_file ?: params.input))

    // Canal des valeurs K pour ADMIXTURE
    // "2,3,4,5" → channel émettant 2, 3, 4, 5 séquentiellement
    // Le mot-clé 'each' dans ADMIXTURE_RUN lancera un job par valeur de K
    ch_k_values = channel
        .of(params.admixture_k.split(',').collect { k -> k.trim() as Integer })
        .flatten()

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 1 — Contrôle qualité des reads
    //
    // FASTQC    : rapport HTML qualité par échantillon (avant trimming)
    // FASTP     : trimming adaptateurs + filtrage qualité + rapport JSON
    // MULTIQC   : agrège tous les rapports en un seul HTML interactif
    //
    // Flux de données :
    //   ch_reads → FASTQC (parallèle pour chaque échantillon)
    //   ch_reads → FASTP  (parallèle pour chaque échantillon)
    //   FASTQC.zip + FASTP.json → collect() → MULTIQC
    // ─────────────────────────────────────────────────────────────────────────
    FASTQC(ch_reads)
    FASTP(ch_reads)

    // collect() attend que TOUS les échantillons soient traités
    // avant de lancer MultiQC (qui a besoin de tous les rapports)
    ch_qc_reports = FASTQC.out.zip.map { meta_file -> meta_file[1] }
        .mix(FASTP.out.json.map { meta_file -> meta_file[1] })
        .collect()
    MULTIQC_QC(ch_qc_reports, 'qc')

    // Les reads filtrés par fastp alimentent l'alignement
    ch_trimmed = FASTP.out.reads

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 2 — Alignement sur le génome de référence Amel_HAv3.1
    //
    // BWA_MEM2_INDEX   : indexe le génome UNE SEULE FOIS
    //                    (Nextflow ne le relancera pas si déjà fait)
    // BWA_MEM2_ALIGN   : aligne les reads de chaque échantillon en parallèle
    //                    Ajoute le Read Group (obligatoire pour GATK)
    //                    Pipe BWA → SAMtools pour éviter le SAM intermédiaire
    // SAMTOOLS_SORT    : trie le BAM par coordonnées génomiques (requis GATK)
    // SAMTOOLS_INDEX   : crée l'index .bai (accès rapide aux régions)
    // SAMTOOLS_FLAGSTAT: % reads alignés, dupliqués, etc. → MultiQC
    // ─────────────────────────────────────────────────────────────────────────
    BWA_MEM2_INDEX(ch_genome)

    // BWA_MEM2_ALIGN reçoit : [meta, reads] + index (tuple genome + fichiers index)
    BWA_MEM2_ALIGN(ch_trimmed, BWA_MEM2_INDEX.out.index)

    SAMTOOLS_SORT(BWA_MEM2_ALIGN.out.bam)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.bam)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 3 — Traitement du BAM : suppression des duplicats PCR
    //
    // PICARD_MARKDUPLICATES : identifie les fragments PCR dupliqués
    //   REMOVE_DUPLICATES=false : marque sans supprimer
    //   GATK ignore automatiquement les reads marqués DUPLICATE
    //   Produit : BAM dédupliqué + index .bai + métriques de duplication
    //
    // ch_dedup_bam : canal [meta, bam, bai] pour GATK HaplotypeCaller
    // join() réunit le BAM et son index dans le même tuple par meta.id
    // ─────────────────────────────────────────────────────────────────────────
    PICARD_MARKDUPLICATES(SAMTOOLS_SORT.out.bam)

    ch_dedup_bam = PICARD_MARKDUPLICATES.out.bam
        .join(PICARD_MARKDUPLICATES.out.bai)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 4 — Variant calling individuel (mode GVCF)
    //
    // GATK_HAPLOTYPECALLER : détecte les variants dans chaque échantillon
    //   Mode -ERC GVCF : produit un GVCF (Genomic VCF) qui contient
    //   la confiance de génotypage à CHAQUE position du génome,
    //   pas seulement aux positions variantes.
    //   Indispensable pour le génotypage joint de l'étape suivante.
    // ─────────────────────────────────────────────────────────────────────────
    GATK_HAPLOTYPECALLER(ch_dedup_bam, ch_genome, ch_fai, ch_dict)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 5 — Génotypage joint de tous les échantillons
    //
    // GATK_GENOMICSDBIMPORT : consolide tous les GVCFs en une base de données
    //   collect() attend que TOUS les GVCFs soient produits avant de démarrer
    // GATK_GENOTYPEGVCFS    : génotypage joint → VCF multi-échantillons final
    //   Plus puissant qu'un génotypage individuel : exploite l'info de
    //   toute la cohorte pour appeler les variants rares et corriger
    //   les erreurs de génotypage individuels
    // ─────────────────────────────────────────────────────────────────────────
    ch_all_gvcfs = GATK_HAPLOTYPECALLER.out.gvcf.map { meta_file -> meta_file[1] }
        .mix(GATK_HAPLOTYPECALLER.out.tbi.map { meta_file -> meta_file[1] })
        .collect()

    GATK_GENOMICSDBIMPORT(ch_all_gvcfs, ch_genome)
    GATK_GENOTYPEGVCFS(GATK_GENOMICSDBIMPORT.out.db, ch_genome, ch_fai, ch_dict)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 6 — Filtrage des variants
    //
    // GATK_VARIANTFILTRATION : applique les filtres GATK Best Practices
    //   Filtre les SNPs de mauvaise qualité (QD, FS, MQ, SOR)
    //   Marque les variants FILTER=PASS ou FILTER=nom_du_filtre
    // BCFTOOLS_FILTER        : garde uniquement PASS, bialléliques, SNPs
    // BCFTOOLS_STATS         : calcule Ts/Tv, nb SNPs, distribution MAF
    //   Ts/Tv attendu ~2.0 pour un génome de bonne qualité
    // ─────────────────────────────────────────────────────────────────────────
    GATK_VARIANTFILTRATION(GATK_GENOTYPEGVCFS.out.vcf, ch_genome, ch_fai, ch_dict)
    BCFTOOLS_FILTER(GATK_VARIANTFILTRATION.out.vcf)
    BCFTOOLS_STATS(BCFTOOLS_FILTER.out.vcf)
    ch_filtered_vcf = BCFTOOLS_FILTER.out.vcf

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 7 — Annotation fonctionnelle des variants (SnpEff)
    //
    // SNPEFF_BUILD    : construit une base SnpEff locale (génome + GFF3 NCBI)
    //                   → assemblage identique à celui du variant calling
    // SNPEFF_ANNOTATE : ajoute le champ ANN= (effet + impact) à chaque variant
    //                   + rapport HTML/CSV agrégé par MultiQC
    // SNPEFF_COMPRESS : recompresse (bgzip) et indexe (tabix) le VCF annoté
    //
    // Activé uniquement si params.gff pointe vers un fichier existant.
    // ─────────────────────────────────────────────────────────────────────────
    ch_snpeff_csv = channel.empty()
    if (run_snpeff) {
        SNPEFF_BUILD(ch_genome, channel.value(file(params.gff)))
        SNPEFF_ANNOTATE(ch_filtered_vcf, SNPEFF_BUILD.out.db)
        SNPEFF_COMPRESS(SNPEFF_ANNOTATE.out.vcf)
        ch_snpeff_csv = SNPEFF_ANNOTATE.out.csv
    }

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 8 — Génétique des populations
    //
    // PLINK2_QC     : filtres MAF/geno/mind/HWE + élagage LD
    //                 Produit les fichiers binaires PLINK (.bed/.bim/.fam)
    // PLINK2_PCA    : ACP génomique → structure de population
    //                 Produit .eigenvec (coordonnées) et .eigenval (variance)
    // ADMIXTURE_RUN : proportions d'ascendance pour K=2..5 en parallèle
    //                 Le mot-clé 'each' dans le module lance un job par K
    // VCFTOOLS_FST  : différenciation génétique FST entre populations
    //                 en fenêtres glissantes de 50 kb
    // VCFTOOLS_LD   : déclin du LD (r²) en fonction de la distance physique
    //                 Chez l'abeille : portée du LD ~10-50 kb
    // VCFTOOLS_PI   : diversité nucléotidique π par fenêtres de 50 kb
    // ─────────────────────────────────────────────────────────────────────────
    PLINK2_QC(ch_filtered_vcf)
    ch_plink    = PLINK2_QC.out.plink_files   // SNPs élagués LD  → PCA, ADMIXTURE, kinship
    ch_plink_qc = PLINK2_QC.out.qc_files      // SNPs QC non élagués → tests d'association

    PLINK2_PCA(ch_plink)

    // ch_k_values émet 2, 3, 4, 5 → 4 jobs ADMIXTURE en parallèle
    ADMIXTURE_RUN(ch_plink, ch_k_values)

    VCFTOOLS_FST(ch_filtered_vcf)
    VCFTOOLS_LD(ch_filtered_vcf)
    VCFTOOLS_PI(ch_filtered_vcf)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 9 — GWAS (si colonne 'phenotype' du samplesheet ou --phenotype_file)
    //
    // BUILD_PHENOTYPE : convertit la source des phénotypes aux formats
    //                 GEMMA (-p) et PLINK2 (--pheno), ordre du .fam respecté
    // GEMMA_KINSHIP : calcule la matrice de parenté centrée N×N
    //                 Capture la structure de population ET la parenté cryptique
    // GEMMA_LMM     : GWAS par modèle mixte linéaire (LMM)
    //                 Corrige la stratification via la matrice de parenté
    //                 Calcule 3 tests : Wald, LRT, Score
    //                 Plus robuste que la régression simple sur les PC
    // PLINK2_GWAS   : association linéaire/logistique complémentaire
    //                 Plus rapide, utile pour validation croisée
    // ─────────────────────────────────────────────────────────────────────────
    if (run_gwas) {
        // Construit phenotype.gemma.txt + phenotype.plink.tsv alignés sur le .fam
        BUILD_PHENOTYPE(ch_plink, ch_pheno_source)
        ch_pheno_gemma = BUILD_PHENOTYPE.out.gemma
        ch_pheno_plink = BUILD_PHENOTYPE.out.plink

        // Kinship sur les SNPs élagués LD
        GEMMA_KINSHIP(ch_plink, ch_pheno_gemma)
        // LMM sur le jeu QC complet, corrigé par la kinship
        GEMMA_LMM(
            ch_plink_qc,
            GEMMA_KINSHIP.out.kinship,
            ch_pheno_gemma
        )
        // Le module GEMMA expose 'annotated' (résultats avec lambda GC)
        ch_gwas_results = GEMMA_LMM.out.annotated

        // Association PLINK2 complémentaire sur le même jeu QC
        PLINK2_GWAS(ch_plink_qc, ch_pheno_plink)
    } else {
        log.warn "Aucun phénotype fourni (--phenotype_file ou colonne 'phenotype') — étapes GWAS ignorées."
        ch_gwas_results = channel.empty()
    }

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 10 — Visualisation (figures publication-ready)
    //
    // PLOT_PCA       : ACP colorée par population (PC1 vs PC2/PC3)
    //                  Reçoit eigenvec ET eigenval pour afficher % variance
    // PLOT_ADMIXTURE : barplot des proportions d'ascendance pour K=2..5
    //                  collect() attend tous les fichiers .Q avant de tracer
    // PLOT_LD_DECAY  : courbe r² en fonction de la distance (kb)
    // PLOT_FST       : Manhattan FST par fenêtres chromosomiques
    // PLOT_MANHATTAN : résultats GWAS (seulement si phénotype fourni)
    // PLOT_QQ        : QQ plot avec lambda GC (contrôle inflation)
    // ─────────────────────────────────────────────────────────────────────────
    PLOT_PCA(
        PLINK2_PCA.out.eigenvec,
        PLINK2_PCA.out.eigenval
    )
    PLOT_ADMIXTURE(
        ADMIXTURE_RUN.out.q_files.collect(),
        PLINK2_QC.out.plink_files.map { bed_bim_fam -> bed_bim_fam[2] }.first()
    )
    PLOT_LD_DECAY(VCFTOOLS_LD.out.ld)
    PLOT_FST(VCFTOOLS_FST.out.fst)

    if (run_gwas) {
        PLOT_MANHATTAN(ch_gwas_results)
        PLOT_QQ(ch_gwas_results)
    }

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 11 — Rapport scientifique et MultiQC final
    //
    // GWAS_REPORT    : rapport HTML via R Markdown (rocker/tidyverse)
    //                  Intègre toutes les figures et statistiques clés
    // MULTIQC_FINAL  : agrège les stats alignement + déduplication + variants
    //
    // Note : .map { it[1] } extrait le fichier du tuple [meta, fichier]
    //        car MultiQC ne veut pas les métadonnées
    // ─────────────────────────────────────────────────────────────────────────
    // ch_gwas_results = Channel.empty() si run_gwas est faux → .mix inoffensif
    ch_report_inputs = BCFTOOLS_STATS.out.stats
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { meta_file -> meta_file[1] })
        .mix(PICARD_MARKDUPLICATES.out.metrics.map { meta_file -> meta_file[1] })
        .mix(PLINK2_PCA.out.eigenvec)
        .mix(ADMIXTURE_RUN.out.q_files.flatten())
        .mix(VCFTOOLS_LD.out.ld)
        .mix(VCFTOOLS_FST.out.fst)
        .mix(ch_gwas_results)
        .collect()

    // Rapport R Markdown HTML — conteneur rocker/tidyverse (rmarkdown + pandoc)
    if (run_gwas) {
        GWAS_REPORT(ch_report_inputs, file("${projectDir}/report/gwas_report.Rmd"))
    }

    // MultiQC final agrège alignement + déduplication + stats variants
    // + résumé SnpEff (ch_snpeff_csv est vide si l'annotation est désactivée)
    ch_final_multiqc = BCFTOOLS_STATS.out.stats
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { meta_file -> meta_file[1] })
        .mix(PICARD_MARKDUPLICATES.out.metrics.map { meta_file -> meta_file[1] })
        .mix(ch_snpeff_csv)
        .collect()
    MULTIQC_FINAL(ch_final_multiqc, 'final')

    // ── Handlers de fin (DSL2 v26 : doivent être DANS le workflow) ────────────
         
  
   // fin du workflow
}
