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
    7.  Annotation        SnpEff (Apis_mellifera)
    8.  Génétique pop.    PLINK2 · ADMIXTURE · vcftools
    9.  GWAS              GEMMA LMM · PLINK2
    10. Visualisation     R (Manhattan · QQ · PCA · Admixture · LD · FST)
    11. Rapport           R Markdown HTML/PDF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// ── Import des modules ────────────────────────────────────────────────────────
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
include { GATK_GENOMICSDBIMPORT    } from './modules/gatk'
include { GATK_GENOTYPEGVCFS       } from './modules/gatk'
include { GATK_VARIANTFILTRATION   } from './modules/gatk'
include { BCFTOOLS_FILTER          } from './modules/bcftools'
include { BCFTOOLS_STATS           } from './modules/bcftools'
include { SNPEFF_ANNOTATE          } from './modules/snpeff'
include { PLINK2_QC                } from './modules/plink2'
include { PLINK2_PCA               } from './modules/plink2'
include { PLINK2_GWAS              } from './modules/plink2'
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

// ── Fonction : parser le samplesheet CSV ─────────────────────────────────────
// Lit le samplesheet ligne par ligne et crée un canal Nextflow
// Chaque élément = [meta, [fastq_1, fastq_2]]
// meta = map contenant id, population, sex
def parseSamplesheet(csv) {
    Channel
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
            fq2 ? [meta, [fq1, fq2]] : [meta, [fq1]]
        }
}

// ── Fonction : validation des paramètres obligatoires ─────────────────────────
// Vérifie que les fichiers requis existent avant de démarrer
// En DSL2, les fonctions peuvent être définies en dehors du workflow
def validateParams() {
    if (!params.input) {
        error "ERREUR : --input est obligatoire. Exemple : --input samplesheet.csv"
    }
    if (!params.genome) {
        error "ERREUR : --genome est obligatoire. Exemple : --genome data/reference/Amel_HAv3.1.fa"
    }
    if (!file(params.input).exists()) {
        error "ERREUR : Samplesheet introuvable : ${params.input}"
    }
    if (!file(params.genome).exists()) {
        error "ERREUR : Génome introuvable : ${params.genome}"
    }
}

// ── Workflow principal ────────────────────────────────────────────────────────
// En DSL2, TOUT le code exécutable doit être dans un bloc workflow,
// process ou function. Jamais au niveau global du script.
workflow {

    // Validation des paramètres obligatoires
    validateParams()

    // Message de démarrage — dans le workflow car DSL2 l'exige
    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║        honeybee-gwas-pipeline v${manifest.version}
    ║   WGS & GWAS — Apis mellifera mellifera
    ╚══════════════════════════════════════════════════════════════════╝
      Samplesheet  : ${params.input}
      Génome       : ${params.genome}
      Phénotypes   : ${params.phenotype_file ?: 'non fourni — GWAS ignoré'}
      Sortie       : ${params.outdir}
      Modèle GWAS  : ${params.gwas_model}
      MAF          : ${params.maf}
      ADMIXTURE K  : ${params.admixture_k}
      Profil       : ${workflow.profile}
    """.stripIndent()

    // ── Canaux d'entrée ───────────────────────────────────────────────────────
    // Canal des reads : un élément par échantillon
    ch_reads = parseSamplesheet(params.input)

    // Canal du génome : valeur unique partagée par tous les process
    ch_genome = Channel.value(file(params.genome))

    // Canal des phénotypes : vide si pas fourni (GWAS ignoré)
    ch_pheno = params.phenotype_file
        ? Channel.value(file(params.phenotype_file))
        : Channel.empty()

    // Canal des valeurs K pour ADMIXTURE
    // params.admixture_k = "2,3,4,5" → Channel.of(2, 3, 4, 5)
    ch_k_values = Channel
        .of(params.admixture_k.split(',').collect { it.trim() as Integer })
        .flatten()

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 1 — Contrôle qualité des reads
    // FastQC analyse la qualité brute des reads
    // fastp filtre les adaptateurs et les reads de mauvaise qualité
    // MultiQC agrège tous les rapports QC en un seul HTML interactif
    // ─────────────────────────────────────────────────────────────────────────
    FASTQC(ch_reads)
    FASTP(ch_reads)

    ch_qc_reports = FASTQC.out.zip
        .mix(FASTP.out.json)
        .collect()
    MULTIQC_QC(ch_qc_reports, 'qc')

    // Reads filtrés → étape alignement
    ch_trimmed = FASTP.out.reads

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 2 — Alignement sur le génome de référence Amel_HAv3.1
    // BWA_MEM2_INDEX : indexe le génome une seule fois pour tous les échantillons
    // BWA_MEM2_ALIGN : aligne les reads de chaque échantillon en parallèle
    // SAMTOOLS_SORT  : trie le BAM par coordonnées génomiques (requis par GATK)
    // SAMTOOLS_INDEX : crée l'index .bai pour accès rapide aux régions
    // SAMTOOLS_FLAGSTAT : calcule les statistiques d'alignement
    // ─────────────────────────────────────────────────────────────────────────
    BWA_MEM2_INDEX(ch_genome)
    BWA_MEM2_ALIGN(ch_trimmed, BWA_MEM2_INDEX.out.index)
    SAMTOOLS_SORT(BWA_MEM2_ALIGN.out.bam)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.bam)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 3 — Traitement du BAM
    // Picard MarkDuplicates identifie les duplicats PCR
    // Les duplicats sont marqués (REMOVE_DUPLICATES=false)
    // GATK les ignore automatiquement lors du variant calling
    // ─────────────────────────────────────────────────────────────────────────
    PICARD_MARKDUPLICATES(SAMTOOLS_SORT.out.bam)

    // BAM dédupliqué + index réunis pour GATK HaplotypeCaller
    ch_dedup_bam = PICARD_MARKDUPLICATES.out.bam
        .join(PICARD_MARKDUPLICATES.out.bai)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 4 — Variant calling individuel (mode GVCF)
    // HaplotypeCaller produit un GVCF par échantillon
    // Le GVCF contient la confiance de génotypage à chaque position
    // du génome — pas seulement aux variants détectés
    // ─────────────────────────────────────────────────────────────────────────
    GATK_HAPLOTYPECALLER(ch_dedup_bam, ch_genome)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 5 — Génotypage joint de tous les échantillons
    // GenomicsDBImport consolide tous les GVCFs en une base de données
    // GenotypeGVCFs produit le VCF multi-échantillons final
    // Le génotypage joint exploite l'info de toute la cohorte
    // ─────────────────────────────────────────────────────────────────────────
    ch_all_gvcfs = GATK_HAPLOTYPECALLER.out.gvcf
        .mix(GATK_HAPLOTYPECALLER.out.tbi)
        .collect()

    GATK_GENOMICSDBIMPORT(ch_all_gvcfs, ch_genome)
    GATK_GENOTYPEGVCFS(GATK_GENOMICSDBIMPORT.out.db, ch_genome)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 6 — Filtrage des variants
    // VariantFiltration applique les filtres GATK Best Practices (hard-filter)
    // bcftools_filter garde uniquement les variants PASS bialléliques SNPs
    // bcftools_stats calcule Ts/Tv, nombre de SNPs, distribution MAF
    // ─────────────────────────────────────────────────────────────────────────
    GATK_VARIANTFILTRATION(GATK_GENOTYPEGVCFS.out.vcf, ch_genome)
    BCFTOOLS_FILTER(GATK_VARIANTFILTRATION.out.vcf)
    BCFTOOLS_STATS(BCFTOOLS_FILTER.out.vcf)
    ch_filtered_vcf = BCFTOOLS_FILTER.out.vcf

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 7 — Annotation fonctionnelle des variants
    // SnpEff prédit l'effet de chaque SNP sur les gènes d'Apis mellifera
    // Base de données Apis_mellifera construite sur Amel_HAv3.1
    // ─────────────────────────────────────────────────────────────────────────
    SNPEFF_ANNOTATE(ch_filtered_vcf)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 8 — Génétique des populations
    // PLINK2_QC      : filtres MAF/geno/mind/HWE + élagage LD
    // PLINK2_PCA     : analyse en composantes principales
    // ADMIXTURE_RUN  : structure de population K=2..5 en parallèle
    // VCFTOOLS_FST   : différenciation génétique entre populations
    // VCFTOOLS_LD    : déclin du déséquilibre de liaison
    // VCFTOOLS_PI    : diversité nucléotidique π
    // ─────────────────────────────────────────────────────────────────────────
    PLINK2_QC(ch_filtered_vcf)
    ch_plink = PLINK2_QC.out.plink_files

    PLINK2_PCA(ch_plink)
    ADMIXTURE_RUN(ch_plink, ch_k_values)
    VCFTOOLS_FST(ch_filtered_vcf)
    VCFTOOLS_LD(ch_filtered_vcf)
    VCFTOOLS_PI(ch_filtered_vcf)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 9 — GWAS (uniquement si fichier phénotype fourni)
    // GEMMA_KINSHIP : calcule la matrice de parenté génomique N×N
    // GEMMA_LMM     : GWAS par modèle mixte linéaire
    //                 Corrige la stratification via la matrice de parenté
    // PLINK2_GWAS   : association complémentaire linéaire/logistique
    // ─────────────────────────────────────────────────────────────────────────
    if (params.phenotype_file) {
        GEMMA_KINSHIP(ch_plink)
        GEMMA_LMM(
            ch_plink,
            GEMMA_KINSHIP.out.kinship,
            ch_pheno
        )
        ch_gwas_results = GEMMA_LMM.out.annotated
        PLINK2_GWAS(ch_plink, ch_pheno)
    } else {
        log.warn "Pas de fichier phénotype fourni — étapes GWAS ignorées."
        ch_gwas_results = Channel.empty()
    }

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 10 — Visualisation (figures publication-ready)
    // PLOT_PCA       : structure de population (PC1 vs PC2/PC3)
    // PLOT_ADMIXTURE : barplot proportions d'ascendance (K=2..5)
    // PLOT_LD_DECAY  : courbe de déclin du LD
    // PLOT_FST       : différenciation génétique par fenêtres chromosomiques
    // PLOT_MANHATTAN : résultats GWAS — seulement si phénotype fourni
    // PLOT_QQ        : contrôle inflation génomique (lambda GC)
    // ─────────────────────────────────────────────────────────────────────────
    PLOT_PCA(
        PLINK2_PCA.out.eigenvec,
        PLINK2_PCA.out.eigenval
    )
    PLOT_ADMIXTURE(
        ADMIXTURE_RUN.out.q_files.collect(),
        ch_plink.map { bed, bim, fam -> fam }
    )
    PLOT_LD_DECAY(VCFTOOLS_LD.out.ld)
    PLOT_FST(VCFTOOLS_FST.out.fst)

    if (params.phenotype_file) {
        PLOT_MANHATTAN(ch_gwas_results)
        PLOT_QQ(ch_gwas_results)
    }

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 11 — Rapport scientifique final
    // Collecte tous les résultats intermédiaires
    // Génère un rapport HTML/PDF complet via R Markdown
    // ─────────────────────────────────────────────────────────────────────────
    ch_report_inputs = BCFTOOLS_STATS.out.stats
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { meta, f -> f })
        .mix(PICARD_MARKDUPLICATES.out.metrics.map { meta, f -> f })
        .mix(PLINK2_PCA.out.eigenvec)
        .mix(ADMIXTURE_RUN.out.q_files.flatten())
        .mix(VCFTOOLS_LD.out.ld)
        .mix(VCFTOOLS_FST.out.fst)
        .collect()

    if (params.phenotype_file) {
        ch_report_inputs = ch_report_inputs.mix(ch_gwas_results)
    }

    GWAS_REPORT(ch_report_inputs.collect())

    // MultiQC final — agrège stats alignement, déduplication et variants
    ch_final_multiqc = BCFTOOLS_STATS.out.stats
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { meta, f -> f })
        .mix(PICARD_MARKDUPLICATES.out.metrics.map { meta, f -> f })
        .collect()
    MULTIQC_FINAL(ch_final_multiqc, 'final')
}

// ── Messages de fin ───────────────────────────────────────────────────────────
// workflow.onComplete et workflow.onError sont des handlers spéciaux
// qui peuvent être définis en dehors du bloc workflow
workflow.onComplete {
    def status = workflow.success ? "SUCCÈS" : "ÉCHEC"
    log.info """
    ════════════════════════════════════════════════════════
    Pipeline terminé — ${status}
    ────────────────────────────────────────────────────────
    Durée      : ${workflow.duration}
    Résultats  : ${params.outdir}/
    Rapport    : ${params.outdir}/06_report/gwas_report.html
    MultiQC    : ${params.outdir}/pipeline_info/report.html
    ════════════════════════════════════════════════════════
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline échoué — trace : ${params.outdir}/pipeline_info/trace.txt"
}
