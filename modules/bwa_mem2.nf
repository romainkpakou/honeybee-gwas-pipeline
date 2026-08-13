/*
    MODULE : BWA-MEM2
    Outil   : bwa-mem2 v2.2.1
    Rôle    : Indexation du génome de référence et alignement des reads
    Génome  : Amel_HAv3.1 (GCF_003254395.2) — 16 chromosomes Apis mellifera
    Docker  : quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5

    DEUX PROCESS DANS CE MODULE :
      1. BWA_MEM2_INDEX — indexe le génome (lancé une seule fois)
      2. BWA_MEM2_ALIGN — aligne les reads (lancé pour chaque échantillon)

    FORMAT DE SORTIE : BAM non trié (sera trié par SAMtools dans l'étape suivante)

    READ GROUP (-R) : métadonnées obligatoires pour GATK
      ID  = identifiant unique du run
      SM  = nom de l'échantillon (utilisé par GATK pour nommer les génotypes)
      PL  = plateforme de séquençage (ILLUMINA)
      LB  = librairie (pour détecter les duplicats inter-librairies)
      PU  = unité de plateforme (flowcell + lane)
*/

// ── Process 1 : Indexation du génome ─────────────────────────────────────────
process BWA_MEM2_INDEX {
    tag "Amel_HAv3.1"
    label 'process_high_memory'

    publishDir "${params.outdir}/02_alignment/index", mode: 'copy'

    container 'quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5'

    input:
    path genome

    output:
    tuple path(genome), path("${genome}.*"), emit: index
    path "versions.yml",                      emit: versions

    script:
    """
    bwa-mem2 index ${genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch ${genome}.0123
    touch ${genome}.amb
    touch ${genome}.ann
    touch ${genome}.bwt.2bit.64
    touch ${genome}.pac
    touch versions.yml
    """
}

// ── Process 2 : Alignement des reads ─────────────────────────────────────────
process BWA_MEM2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/02_alignment/bam", mode: 'copy'

    container 'quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5'

    input:
    tuple val(meta), path(reads)
    tuple path(genome), path(index)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path "versions.yml",                      emit: versions

    script:
    // Read Group : obligatoire pour GATK
    // Sans Read Group, GATK HaplotypeCaller refuse de tourner
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}_lib1\\tPU:${meta.id}"
    def reads_input = reads.size() == 2
        ? "${reads[0]} ${reads[1]}"
        : "${reads[0]}"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        -M \\
        ${genome} \\
        ${reads_input} | \\
    samtools view \\
        -@ ${task.cpus} \\
        -bS \\
        -o ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bam
    touch versions.yml
    """
}
