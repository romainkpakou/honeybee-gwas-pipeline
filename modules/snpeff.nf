
/*
    MODULE : SnpEff
    Outil   : SnpEff v5.2
    Rôle    : Annotation fonctionnelle des variants
    Docker  : quay.io/biocontainers/snpeff:5.2--hdfd78af_0

    PRINCIPE :
    Pour chaque variant du VCF, SnpEff consulte l'annotation du génome
    et prédit l'effet biologique : gène / exon / intron touché, changement
    d'acide aminé (missense), codon stop prématuré (stop_gained), etc.

    IMPACT des variants :
      HIGH     : perte de fonction probable (stop_gained, frameshift)
      MODERATE : effet modéré (missense_variant)
      LOW      : effet faible (synonymous_variant)
      MODIFIER : régions non-codantes (intron, intergénique)

    APPROCHE HORS-LIGNE :
    Plutôt que `snpEff download Apis_mellifera` (nécessite internet, et la
    base publique peut reposer sur un assemblage différent), on CONSTRUIT
    la base localement à partir du génome de référence et de son annotation
    GFF3 (NCBI). La base correspond ainsi exactement à l'assemblage utilisé
    pour l'alignement et le variant calling.

    Les flags -noCheckCds / -noCheckProtein évitent d'avoir besoin des
    fichiers CDS/protéines de référence (non disponibles hors-ligne). Les
    prédictions d'effet restent valides ; seule la validation croisée des
    traductions est ignorée.

    TROIS PROCESS :
      1. SNPEFF_BUILD    — construit la base SnpEff (génome + GFF3)
      2. SNPEFF_ANNOTATE — annote le VCF, produit le rapport HTML/CSV
      3. SNPEFF_COMPRESS — recompresse et indexe le VCF annoté (bgzip + tabix)
*/

// ── Process 1 : Construction de la base SnpEff ────────────────────────────────
process SNPEFF_BUILD {
    tag "snpeff_build:${genome.baseName}"
    label 'process_medium'

    publishDir "${params.outdir}/03_variants/annotated", mode: 'copy', pattern: 'snpEff.config'

    container 'quay.io/biocontainers/snpeff:5.2--hdfd78af_0'

    input:
    path genome
    path gff

    output:
    tuple val(db), path("snpeff_data"), path("snpEff.config"), emit: db
    path "versions.yml",                                       emit: versions

    script:
    // Nom de la base = nom du FASTA sans extension (.fa/.fasta/.fna[.gz])
    db = genome.name.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
    def avail_mem = task.memory ? task.memory.toGiga() : 6
    """
    mkdir -p snpeff_data/${db}

    # SnpEff attend le génome (sequences.fa) et l'annotation (genes.gff)
    # dans data/<db>/ — décompressés si besoin
    case "${genome}" in
        *.gz) zcat ${genome} > snpeff_data/${db}/sequences.fa ;;
        *)    cp -L ${genome} snpeff_data/${db}/sequences.fa ;;
    esac
    case "${gff}" in
        *.gz) zcat ${gff} > snpeff_data/${db}/genes.gff ;;
        *)    cp -L ${gff} snpeff_data/${db}/genes.gff ;;
    esac

    # Configuration locale minimale (data.dir relatif au répertoire de travail)
    {
        echo "data.dir = ./snpeff_data/"
        echo "${db}.genome : ${params.genome_species} (${db})"
    } > snpEff.config

    snpEff -Xmx${avail_mem}g build -gff3 -v \\
        -noCheckCds -noCheckProtein \\
        -c snpEff.config \\
        ${db}

    test -s snpeff_data/${db}/snpEffectPredictor.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(snpEff -version 2>&1 | grep -oiE '[0-9]+\\.[0-9]+[a-z]?' | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p snpeff_data/${genome.baseName}
    touch snpeff_data/${genome.baseName}/snpEffectPredictor.bin
    touch snpEff.config versions.yml
    """
}

// ── Process 2 : Annotation du VCF ─────────────────────────────────────────────
process SNPEFF_ANNOTATE {
    tag "snpeff_annotate"
    label 'process_medium'

    // Le VCF brut n'est pas publié : SNPEFF_COMPRESS publie la version bgzip+tabix
    publishDir "${params.outdir}/03_variants/annotated", mode: 'copy',
               pattern: "{snpEff_summary.*,versions.yml}"

    container 'quay.io/biocontainers/snpeff:5.2--hdfd78af_0'

    input:
    tuple path(vcf), path(tbi)
    tuple val(db), path(snpeff_data), path(config)

    output:
    path "cohort.annotated.vcf",    emit: vcf
    path "snpEff_summary.html",     emit: report
    path "snpEff_summary.csv",      emit: csv
    path "snpEff_summary.genes.txt", emit: genes
    path "versions.yml",            emit: versions

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 6
    """
    # snpEff lit directement le VCF gzippé
    snpEff -Xmx${avail_mem}g ann -v \\
        -c ${config} \\
        -stats snpEff_summary.html \\
        -csvStats snpEff_summary.csv \\
        ${db} ${vcf} > cohort.annotated.vcf

    echo "=== Distribution des effets (les plus fréquents) ==="
    grep -v '^#' cohort.annotated.vcf \\
        | grep -oE 'ANN=[^[:space:]]+' \\
        | grep -oE '\\|(missense_variant|synonymous_variant|stop_gained|stop_lost|start_lost|frameshift_variant|splice_[a-z_]+variant|intron_variant|intergenic_region|[35]_prime_UTR_variant|upstream_gene_variant|downstream_gene_variant)\\|' \\
        | sort | uniq -c | sort -rn | head -15 || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(snpEff -version 2>&1 | grep -oiE '[0-9]+\\.[0-9]+[a-z]?' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch cohort.annotated.vcf snpEff_summary.html snpEff_summary.csv snpEff_summary.genes.txt versions.yml
    """
}

// ── Process 3 : Recompression + indexation du VCF annoté ──────────────────────
process SNPEFF_COMPRESS {
    tag "snpeff_compress"
    label 'process_low'

    publishDir "${params.outdir}/03_variants/annotated", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    path vcf

    output:
    tuple path("cohort.annotated.vcf.gz"),
          path("cohort.annotated.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools view -Oz -o cohort.annotated.vcf.gz ${vcf}
    bcftools index --tbi cohort.annotated.vcf.gz
    """

    stub:
    """
    touch cohort.annotated.vcf.gz cohort.annotated.vcf.gz.tbi
    """
}
