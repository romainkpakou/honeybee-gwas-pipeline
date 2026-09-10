
/*
    MODULE : MultiQC
    Outil   : MultiQC v1.21
    Rôle    : Agrège les rapports QC de tous les échantillons
              en un seul rapport HTML interactif
    Input   : Collection de fichiers rapports (FastQC zip, fastp JSON,
              Picard metrics, SAMtools flagstat, bcftools stats...)
    Output  : Rapport HTML unique + répertoire de données
    Docker  : quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0

    NOTE : Ce module est appelé deux fois dans le pipeline :
      1. Après QC des reads (FastQC + fastp)
      2. En fin de pipeline (tous les outils)
    Le paramètre 'stage' distingue les deux rapports.
*/

process MULTIQC {
    tag "${stage}"
    label 'process_low'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    input:
    path(reports)
    val(stage)

    output:
    path "multiqc_${stage}.html",   emit: report
    path "multiqc_${stage}_data/",  emit: data
    path "versions.yml",            emit: versions

    script:
    """
    multiqc \\
        --title "honeybee-gwas-pipeline — ${stage} QC" \\
        --filename multiqc_${stage}.html \\
        --force \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_${stage}.html
    mkdir -p multiqc_${stage}_data
    touch versions.yml
    """
}
