/*
    MODULE : GATK4
    Outil   : GATK v4.5.0.0
    Rôle    : Variant calling suivant les GATK Best Practices
    Docker  : broadinstitute/gatk:4.5.0.0

    4 PROCESS DANS CE MODULE :

    1. GATK_HAPLOTYPECALLER
       Détecte les variants dans chaque échantillon individuellement.
       Mode GVCF : produit un fichier intermédiaire contenant la
       confiance de génotypage à CHAQUE position du génome, pas
       seulement aux variants détectés. Indispensable pour le
       génotypage joint.

    2. GATK_GENOMICSDBIMPORT
       Consolide tous les GVCFs individuels dans une base de données
       génomique optimisée. Plus efficace que CombineGVCFs pour
       de grandes cohortes.

    3. GATK_GENOTYPEGVCFS
       Génotypage joint de tous les échantillons simultanément.
       Plus puissant qu'un génotypage individuel car exploite
       l'information de toute la cohorte pour appeler les variants
       rares et corriger les erreurs de génotypage.

    4. GATK_VARIANTFILTRATION
       Filtre les variants de mauvaise qualité selon les critères
       GATK Best Practices (hard-filter) :
       QD  < 2.0  : variant dans une région à faible couverture
       FS  > 60.0 : biais de brin important
       MQ  < 40.0 : mauvaise qualité d'alignement
       SOR > 3.0  : autre mesure de biais de brin
*/

// ── Process 1 : HaplotypeCaller (GVCF mode) ───────────────────────────────────
process GATK_HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/03_variants/gvcf", mode: 'copy'

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path genome

    output:
    tuple val(meta), path("${meta.id}.g.vcf.gz"),     emit: gvcf
    tuple val(meta), path("${meta.id}.g.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                               emit: versions

    script:
    """
    # Indexer le génome si nécessaire
    if [ ! -f ${genome}.fai ]; then
        samtools faidx ${genome}
    fi

    # Créer le dictionnaire si nécessaire
    if [ ! -f \${genome%.fa}.dict ] && [ ! -f \${genome%.fasta}.dict ]; then
        gatk CreateSequenceDictionary -R ${genome}
    fi

    mkdir -p tmp

    gatk ${params.gatk_java} HaplotypeCaller \\
        -R ${genome} \\
        -I ${bam} \\
        -O ${meta.id}.g.vcf.gz \\
        -ERC GVCF \\
        --sample-name ${meta.id} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --tmp-dir ./tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | head -1 | sed 's/The Genome Analysis Toolkit (GATK) v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.g.vcf.gz
    touch ${meta.id}.g.vcf.gz.tbi
    touch versions.yml
    """
}

// ── Process 2 : GenomicsDBImport ──────────────────────────────────────────────
process GATK_GENOMICSDBIMPORT {
    tag "all_samples"
    label 'process_high_memory'

    publishDir "${params.outdir}/03_variants/genomicsdb", mode: 'copy'

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    path(gvcfs)
    path genome

    output:
    path "genomicsdb/", emit: db
    path "versions.yml", emit: versions

    script:
    // Construire la liste des arguments -V pour chaque GVCF
    def vcf_args = gvcfs
        .findAll { it.name.endsWith('.g.vcf.gz') }
        .collect { "-V ${it}" }
        .join(" \\\n        ")

    """
    # Récupérer la liste des chromosomes depuis le génome
    samtools faidx ${genome}
    cut -f1 ${genome}.fai > chromosomes.list

    mkdir -p tmp

    gatk ${params.gatk_java} GenomicsDBImport \\
        ${vcf_args} \\
        --genomicsdb-workspace-path genomicsdb \\
        -L chromosomes.list \\
        --reader-threads ${task.cpus} \\
        --batch-size 50 \\
        --tmp-dir ./tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | head -1 | sed 's/The Genome Analysis Toolkit (GATK) v//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p genomicsdb
    touch versions.yml
    """
}

// ── Process 3 : GenotypeGVCFs ─────────────────────────────────────────────────
process GATK_GENOTYPEGVCFS {
    tag "joint_genotyping"
    label 'process_high'

    publishDir "${params.outdir}/03_variants/vcf", mode: 'copy'

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    path db
    path genome

    output:
    tuple path("cohort.vcf.gz"), path("cohort.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                      emit: versions

    script:
    """
    mkdir -p tmp

    gatk ${params.gatk_java} GenotypeGVCFs \\
        -R ${genome} \\
        -V gendb://${db} \\
        -O cohort.vcf.gz \\
        --tmp-dir ./tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | head -1 | sed 's/The Genome Analysis Toolkit (GATK) v//')
    END_VERSIONS
    """

    stub:
    """
    touch cohort.vcf.gz cohort.vcf.gz.tbi
    touch versions.yml
    """
}

// ── Process 4 : VariantFiltration ────────────────────────────────────────────
process GATK_VARIANTFILTRATION {
    tag "variant_filtration"
    label 'process_medium'

    publishDir "${params.outdir}/03_variants/filtered", mode: 'copy'

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple path(vcf), path(tbi)
    path genome

    output:
    tuple path("cohort.filtered.vcf.gz"),
          path("cohort.filtered.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                       emit: versions

    script:
    """
    mkdir -p tmp

    # Extraire uniquement les SNPs
    gatk ${params.gatk_java} SelectVariants \\
        -R ${genome} \\
        -V ${vcf} \\
        --select-type-to-include SNP \\
        -O snps_only.vcf.gz

    # Appliquer les filtres hard-filter GATK Best Practices
    gatk ${params.gatk_java} VariantFiltration \\
        -R ${genome} \\
        -V snps_only.vcf.gz \\
        --filter-expression "QD < ${params.snp_qd}" \\
            --filter-name "QD_filter" \\
        --filter-expression "FS > ${params.snp_fs}" \\
            --filter-name "FS_filter" \\
        --filter-expression "MQ < ${params.snp_mq}" \\
            --filter-name "MQ_filter" \\
        --filter-expression "MQRankSum < ${params.snp_mqranksum}" \\
            --filter-name "MQRankSum_filter" \\
        --filter-expression "ReadPosRankSum < ${params.snp_readposrank}" \\
            --filter-name "ReadPosRankSum_filter" \\
        --filter-expression "SOR > ${params.snp_sor}" \\
            --filter-name "SOR_filter" \\
        -O cohort.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | head -1 | sed 's/The Genome Analysis Toolkit (GATK) v//')
    END_VERSIONS
    """

    stub:
    """
    touch cohort.filtered.vcf.gz cohort.filtered.vcf.gz.tbi
    touch versions.yml
    """
}
