
/*
    MODULE : fastp
    Outil   : fastp v0.23.4
    Rôle    : Trimming des adaptateurs + filtrage qualité des reads
    Input   : [meta, reads]
    Output  : Reads filtrés (FASTQ gzippés) + rapport JSON + rapport HTML
    Docker  : quay.io/biocontainers/fastp:0.23.4--hadf994f_2
*/

process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/01_qc/fastp", mode: "copy"

    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"),             emit: json
    tuple val(meta), path("*.html"),             emit: html
    path "versions.yml",                          emit: versions

    script:
    def prefix = meta.id
    def paired = reads.size() == 2

    def input_args = paired
        ? "--in1 ${reads[0]} --in2 ${reads[1]} --out1 ${prefix}_R1_trimmed.fastq.gz --out2 ${prefix}_R2_trimmed.fastq.gz"
        : "--in1 ${reads[0]} --out1 ${prefix}_trimmed.fastq.gz"

    """
    fastp \\
        ${input_args} \\
        --thread ${task.cpus} \\
        --qualified_quality_phred ${params.min_quality} \\
        --length_required ${params.min_length} \\
        --detect_adapter_for_pe \\
        --correction \\
        --overrepresentation_analysis \\
        --json ${prefix}_fastp.json \\
        --html ${prefix}_fastp.html \\
        --report_title "${prefix} — fastp QC report"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | head -1 | sed 's/fastp //')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    def paired = reads.size() == 2
    """
    ${paired ? "touch ${prefix}_R1_trimmed.fastq.gz ${prefix}_R2_trimmed.fastq.gz" : "touch ${prefix}_trimmed.fastq.gz"}
    touch ${prefix}_fastp.json ${prefix}_fastp.html
    touch versions.yml
    """
}
