#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    honeybee-gwas-pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Author      : Romain KPAKOU
    Description : End-to-end WGS and GWAS pipeline for Apis mellifera mellifera
    Version     : 1.0.0
    GitHub      : https://github.com/romainkpakou/honeybee-gwas-pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ÉTAPES DU PIPELINE
    ──────────────────
    1.  QC reads          FastQC, fastp, MultiQC
    2.  Alignement        BWA-MEM2, SAMtools
    3.  BAM processing    Picard MarkDuplicates
    4.  Variant calling   GATK HaplotypeCaller (GVCF)
    5.  Génotypage joint  GATK GenomicsDBImport, GenotypeGVCFs
    6.  Filtrage          GATK VariantFiltration, bcftools
    7.  Annotation        SnpEff
    8.  Génétique pop.    PLINK2 QC+PCA, ADMIXTURE, vcftools FST+LD
    9.  GWAS              GEMMA kinship + LMM, PLINK2
    10. Visualisation     R (Manhattan, QQ, PCA, Admixture, LD, FST)
    11. Rapport           R Markdown HTML/PDF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// ── Aide ──────────────────────────────────────────────────────────────────────
if (params.help) {
    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║        honeybee-gwas-pipeline v${manifest.version}
    ║   WGS & GWAS — Apis mellifera mellifera
    ╚══════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf \\
            --input samplesheet.csv \\
            --genome data/reference/Amel_HAv3.1.fa \\
            --phenotype_file data/phenotypes/phenotypes.txt \\
            --outdir results \\
            -profile docker

    Arguments obligatoires:
        --input           Samplesheet CSV (sample,fastq_1,fastq_2,sex,population)
        --genome          Génome de référence FASTA (Amel_HAv3.1)

    Arguments optionnels:
        --phenotype_file  Fichier phénotypes pour GWAS
        --outdir          Répertoire de sortie (défaut: results)
        --maf             Seuil MAF (défaut: 0.05)
        --gwas_model      Modèle GWAS: lmm ou linear (défaut: lmm)
        --admixture_k     Valeurs K ADMIXTURE (défaut: 2,3,4,5)
        --help            Afficher cette aide

    Profils:
        -profile docker       Docker (local)
        -profile singularity  Singularity (HPC)
        -profile slurm        SLURM + Singularity
        -profile test         Données de test
    """.stripIndent()
    System.exit(0)
}

// ── Import des modules ────────────────────────────────────────────────────────
include { FASTQC                    } from './modules/fastqc'
include { FASTP                     } from './modules/fastp'
include { MULTIQC as MULTIQC_QC     } from './modules/multiqc'
include { MULTIQC as MULTIQC_FINAL  } from './modules/multiqc'
include { BWA_MEM2_INDEX            } from './modules/bwa_mem2'
include { BWA_MEM2_ALIGN            } from './modules/bwa_mem2'
include { SAMTOOLS_SORT             } from './modules/samtools'
include { SAMTOOLS_INDEX            } from './modules/samtools'
include { SAMTOOLS_FLAGSTAT         } from './modules/samtools'
include { PICARD_MARKDUPLICATES     } from './modules/picard'
include { GATK_HAPLOTYPECALLER      } from './modules/gatk'
include { GATK_GENOMICSDBIMPORT     } from './modules/gatk'
include { GATK_GENOTYPEGVCFS        } from './modules/gatk'
include { GATK_VARIANTFILTRATION    } from './modules/gatk'
include { BCFTOOLS_FILTER           } from './modules/bcftools'
include { BCFTOOLS_STATS            } from './modules/bcftools'
include { SNPEFF_ANNOTATE           } from './modules/snpeff'
include { PLINK2_QC                 } from './modules/plink2'
include { PLINK2_PCA                } from './modules/plink2'
include { PLINK2_GWAS               } from './modules/plink2'
include { ADMIXTURE_RUN             } from './modules/admixture'
include { VCFTOOLS_FST              } from './modules/vcftools'
include { VCFTOOLS_LD               } from './modules/vcftools'
include { VCFTOOLS_PI               } from './modules/vcftools'
include { GEMMA_KINSHIP             } from './modules/gemma'
include { GEMMA_LMM                 } from './modules/gemma'
include { PLOT_MANHATTAN            } from './modules/r_plots'
include { PLOT_QQ                   } from './modules/r_plots'
include { PLOT_PCA                  } from './modules/r_plots'
include { PLOT_ADMIXTURE            } from './modules/r_plots'
include { PLOT_LD_DECAY             } from './modules/r_plots'
include { PLOT_FST                  } from './modules/r_plots'
include { GWAS_REPORT               } from './modules/report'

// ── Fonction : parser le samplesheet ─────────────────────────────────────────
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
            def fq2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
            fq2 ? [meta, [fq1, fq2]] : [meta, [fq1]]
        }
}

// ── Log de démarrage ──────────────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════════════╗
║        honeybee-gwas-pipeline v${manifest.version}
║   WGS & GWAS — Apis mellifera mellifera
╚══════════════════════════════════════════════════════════════════╝
  Samplesheet     : ${params.input}
  Génome          : ${params.genome}
  Phénotypes      : ${params.phenotype_file ?: 'non fourni — GWAS ignoré'}
  Sortie          : ${params.outdir}
  Modèle GWAS     : ${params.gwas_model}
  MAF             : ${params.maf}
  ADMIXTURE K     : ${params.admixture_k}
  Profil          : ${workflow.profile}
""".stripIndent()

// ── Workflow principal ────────────────────────────────────────────────────────
workflow {

    // ── Canaux d'entrée ───────────────────────────────────────────────────────
    ch_reads  = parseSamplesheet(params.input)
    ch_genome = Channel.value(file(params.genome))
    ch_pheno  = params.phenotype_file
        ? Channel.value(file(params.phenotype_file))
        : Channel.empty()

    // Valeurs de K pour ADMIXTURE
    ch_k_values = Channel.of(
        params.admixture_k.split(',').collect { it.trim() as Integer }
    ).flatten()

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 1 — Contrôle qualité des reads
    // ─────────────────────────────────────────────────────────────────────────
    FASTQC(ch_reads)
    FASTP(ch_reads)

    // Rapport MultiQC QC
    ch_qc_reports = FASTQC.out.zip
        .mix(FASTP.out.json)
        .collect()
    MULTIQC_QC(ch_qc_reports, 'qc')

    // Reads filtrés → étape suivante
    ch_trimmed = FASTP.out.reads

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 2 — Alignement
    // ─────────────────────────────────────────────────────────────────────────
    BWA_MEM2_INDEX(ch_genome)
    BWA_MEM2_ALIGN(ch_trimmed, BWA_MEM2_INDEX.out.index)

    SAMTOOLS_SORT(BWA_MEM2_ALIGN.out.bam)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.bam)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 3 — Traitement BAM
    // ─────────────────────────────────────────────────────────────────────────
    PICARD_MARKDUPLICATES(SAMTOOLS_SORT.out.bam)

    // BAM + index pour GATK
    ch_dedup_bam = PICARD_MARKDUPLICATES.out.bam
        .join(PICARD_MARKDUPLICATES.out.bai)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 4 — Variant calling (GVCF par échantillon)
    // ─────────────────────────────────────────────────────────────────────────
    GATK_HAPLOTYPECALLER(ch_dedup_bam, ch_genome)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 5 — Génotypage joint
    // ─────────────────────────────────────────────────────────────────────────
    // Collecter tous les GVCFs avant GenomicsDBImport
    ch_all_gvcfs = GATK_HAPLOTYPECALLER.out.gvcf
        .mix(GATK_HAPLOTYPECALLER.out.tbi)
        .collect()

    GATK_GENOMICSDBIMPORT(ch_all_gvcfs, ch_genome)
    GATK_GENOTYPEGVCFS(GATK_GENOMICSDBIMPORT.out.db, ch_genome)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 6 — Filtrage des variants
    // ─────────────────────────────────────────────────────────────────────────
    GATK_VARIANTFILTRATION(GATK_GENOTYPEGVCFS.out.vcf, ch_genome)
    BCFTOOLS_FILTER(GATK_VARIANTFILTRATION.out.vcf)
    BCFTOOLS_STATS(BCFTOOLS_FILTER.out.vcf)
    ch_filtered_vcf = BCFTOOLS_FILTER.out.vcf

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 7 — Annotation fonctionnelle
    // ─────────────────────────────────────────────────────────────────────────
    SNPEFF_ANNOTATE(ch_filtered_vcf)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 8 — Génétique des populations
    // ─────────────────────────────────────────────────────────────────────────
    PLINK2_QC(ch_filtered_vcf)
    ch_plink = PLINK2_QC.out.plink_files

    PLINK2_PCA(ch_plink)
    ADMIXTURE_RUN(ch_plink, ch_k_values)
    VCFTOOLS_FST(ch_filtered_vcf)
    VCFTOOLS_LD(ch_filtered_vcf)
    VCFTOOLS_PI(ch_filtered_vcf)

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 9 — GWAS
    // ─────────────────────────────────────────────────────────────────────────
    if (params.phenotype_file) {

        // GEMMA : matrice de parenté puis LMM
        GEMMA_KINSHIP(ch_plink)
        GEMMA_LMM(
            ch_plink,
            GEMMA_KINSHIP.out.kinship,
            ch_pheno
        )
        ch_gwas_results = GEMMA_LMM.out.annotated

        // PLINK2 en complément
        PLINK2_GWAS(ch_plink, ch_pheno)

    } else {
        log.warn "Pas de fichier phénotype fourni — étapes GWAS ignorées."
        ch_gwas_results = Channel.empty()
    }

    // ─────────────────────────────────────────────────────────────────────────
    // ÉTAPE 10 — Visualisation
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
    // ÉTAPE 11 — Rapport final
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

    // Rapport MultiQC final
    ch_final_multiqc = BCFTOOLS_STATS.out.stats
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { meta, f -> f })
        .mix(PICARD_MARKDUPLICATES.out.metrics.map { meta, f -> f })
        .collect()
    MULTIQC_FINAL(ch_final_multiqc, 'final')
}

// ── Messages de fin ───────────────────────────────────────────────────────────
workflow.onComplete {
    def status = workflow.success ? "SUCCÈS" : "ÉCHEC"
    log.info """
    ════════════════════════════════════════════════════
    Pipeline terminé — ${status}
    ────────────────────────────────────────────────────
    Durée       : ${workflow.duration}
    Résultats   : ${params.outdir}/
    Rapport     : ${params.outdir}/06_report/gwas_report.html
    MultiQC     : ${params.outdir}/pipeline_info/report.html
    ════════════════════════════════════════════════════
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline échoué — consultez la trace : ${params.outdir}/pipeline_info/trace.txt"
}
