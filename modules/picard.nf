/*
    MODULE : Picard MarkDuplicates
    Outil   : Picard v3.1.1
    Rôle    : Identification et marquage des duplicats PCR
    Input   : BAM trié et indexé
    Output  : BAM dédupliqué + index + fichier de métriques
    Docker  : quay.io/biocontainers/picard:3.1.1--hdfd78af_0

    PRINCIPE :
    Picard compare les coordonnées de début et de fin de chaque paire
    de reads. Si deux paires ont exactement les mêmes coordonnées,
    elles sont considérées comme des duplicats PCR.
    Seule la paire avec le meilleur score de qualité est conservée
    comme "primaire" — les autres sont marquées DUPLICATE dans le BAM.
    GATK ignore automatiquement les reads marqués DUPLICATE lors
    du variant calling.

    REMOVE_DUPLICATES=false : on marque sans supprimer pour garder
    la possibilité d'analyser les taux de duplication.

    MÉTRIQUES IMPORTANTES dans le fichier .dup_metrics.txt :
    - ESTIMATED_LIBRARY_SIZE : taille estimée de la librairie
    - PERCENT_DUPLICATION    : % de reads dupliqués (idéal < 20%)
    - READ_PAIRS_EXAMINED    : nombre total de paires analysées
*/

process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_high_memory'

    publishDir "${params.outdir}/02_alignment/dedup", mode: 'copy'

    container 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.dedup.bam"),      emit: bam
    tuple val(meta), path("${meta.id}.dedup.bam.bai"),  emit: bai
    tuple val(meta), path("${meta.id}.dup_metrics.txt"), emit: metrics
    path "versions.yml",                                 emit: versions

    script:
    """
    # Créer le répertoire temporaire
    mkdir -p tmp

    picard ${params.gatk_java} MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${meta.id}.dedup.bam \\
        METRICS_FILE=${meta.id}.dup_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=./tmp

    # Picard crée l'index avec extension .bai
    # On le renomme au format attendu par GATK : .bam.bai
    mv ${meta.id}.dedup.bai ${meta.id}.dedup.bam.bai 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -i version | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.dedup.bam
    touch ${meta.id}.dedup.bam.bai
    touch ${meta.id}.dup_metrics.txt
    touch versions.yml
    """
}
