/*
    MODULE : FastQC
    Outil   : FastQC v0.12.1
    Rôle    : Contrôle qualité des reads FASTQ bruts
    Input   : [meta, reads] — meta = map d'infos sur l'échantillon
                              reads = liste de fichiers FASTQ
    Output  : Rapport HTML + archive ZIP par échantillon
    Docker  : quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
*/

process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/01_qc/fastqc/${meta.id}", mode: 'copy'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    path "versions.yml",             emit: versions

    script:
    def read_list = reads.join(' ')
    """
    fastqc \\
        --outdir . \\
        --threads ${task.cpus} \\
        --format fastq \\
        ${read_list}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_R1_fastqc.html ${meta.id}_R1_fastqc.zip
    touch ${meta.id}_R2_fastqc.html ${meta.id}_R2_fastqc.zip
    touch versions.yml
    """
}

