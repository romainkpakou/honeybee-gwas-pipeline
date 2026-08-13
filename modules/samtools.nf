nextflow.enable.dsl = 2

/*
    MODULE : SAMtools
    Outil   : SAMtools v1.19.2
    Rôle    : Tri, indexation et statistiques des fichiers BAM
    Docker  : quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1

    TROIS PROCESS DANS CE MODULE :

    1. SAMTOOLS_SORT
       Trie le BAM par coordonnées chromosomiques.
       GATK exige un BAM trié — sans ça HaplotypeCaller refuse de tourner.
       Option -m 2G : alloue 2 Go de RAM par thread de tri.

    2. SAMTOOLS_INDEX
       Crée le fichier d'index .bai (Binary Alignment Index).
       Permet à GATK et aux outils de visualisation (IGV) d'accéder
       directement à n'importe quelle région du génome sans lire
       tout le fichier BAM.

    3. SAMTOOLS_FLAGSTAT
       Calcule les statistiques d'alignement :
       - % reads alignés
       - % reads en paires correctement alignées
       - % reads dupliqués
       - % reads avec le compagnon aligné sur un autre chromosome
       Collecté par MultiQC pour comparaison entre échantillons.
*/

// ── Process 1 : Tri du BAM ────────────────────────────────────────────────────
process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/02_alignment/bam", mode: 'copy'

    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam
    path "versions.yml",                             emit: versions

    script:
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -m 2G \\
        -o ${meta.id}.sorted.bam \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    touch versions.yml
    """
}

// ── Process 2 : Indexation du BAM ────────────────────────────────────────────
process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/02_alignment/bam", mode: 'copy'

    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam}.bai"), emit: bai
    path "versions.yml",                  emit: versions

    script:
    """
    samtools index -@ ${task.cpus} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bai
    touch versions.yml
    """
}

// ── Process 3 : Statistiques d'alignement ────────────────────────────────────
process SAMTOOLS_FLAGSTAT {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/02_alignment/flagstat", mode: 'copy'

    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.flagstat"), emit: flagstat
    path "versions.yml",                           emit: versions

    script:
    """
    samtools flagstat \\
        -@ ${task.cpus} \\
        ${bam} > ${meta.id}.flagstat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.flagstat
    touch versions.yml
    """
}
