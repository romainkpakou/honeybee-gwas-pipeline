nextflow.enable.dsl = 2

/*
    MODULE : SnpEff
    Outil   : SnpEff v5.2
    Rôle    : Annotation fonctionnelle des variants
    Base de données : Apis_mellifera (construite sur Amel_HAv3.1)
    Docker  : quay.io/biocontainers/snpeff:5.2--hdfd78af_0

    PRINCIPE :
    Pour chaque SNP dans le VCF, SnpEff consulte l'annotation du
    génome Amel_HAv3.1 et prédit l'effet biologique du variant :
    - Est-il dans un gène ? Dans un exon ? Dans un intron ?
    - Change-t-il un acide aminé (missense) ?
    - Crée-t-il un codon stop prématuré (stop_gained) ?
    - Est-il synonyme (même acide aminé) ?

    IMPACT des variants :
    HIGH     : perte de fonction probable (stop_gained, frameshift)
    MODERATE : effet modéré (missense_variant)
    LOW      : effet faible (synonymous_variant)
    MODIFIER : régions non-codantes (intron, intergenic)

    UTILITÉ POUR LE GWAS :
    Quand on identifie des SNPs significatifs en GWAS, l'annotation
    SnpEff permet de comprendre si le variant impacte un gène connu,
    et quel est son effet probable sur la protéine correspondante.
    C'est l'étape qui donne du sens biologique aux résultats GWAS.

    RAPPORT HTML :
    SnpEff génère un rapport HTML résumant la distribution des effets
    dans tout le VCF — collecté par MultiQC pour le rapport final.
*/

process SNPEFF_ANNOTATE {
    tag "snpeff"
    label 'process_medium'

    publishDir "${params.outdir}/03_variants/annotated", mode: 'copy'

    container 'quay.io/biocontainers/snpeff:5.2--hdfd78af_0'

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("cohort.annotated.vcf.gz"),
          path("cohort.annotated.vcf.gz.tbi"), emit: vcf
    path "snpEff_summary.html",                emit: report
    path "snpEff_genes.txt",                   emit: genes
    path "versions.yml",                        emit: versions

    script:
    """
    # Télécharger la base de données Apis mellifera si absente
    snpEff download Apis_mellifera -v

    # Annoter les variants
    snpEff -Xmx8g \\
        -v \\
        -stats snpEff_summary.html \\
        -csvStats snpEff_genes.txt \\
        Apis_mellifera \\
        ${vcf} | \\
    bgzip -c > cohort.annotated.vcf.gz

    # Indexer le VCF annoté
    tabix -p vcf cohort.annotated.vcf.gz

    # Résumé des effets dans les logs
    echo "=== Distribution des effets SnpEff ==="
    bcftools stats cohort.annotated.vcf.gz | grep "^SN" | head -10

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(snpEff -version 2>&1 | head -2 | tail -1)
    END_VERSIONS
    """

    stub:
    """
    touch cohort.annotated.vcf.gz cohort.annotated.vcf.gz.tbi
    touch snpEff_summary.html snpEff_genes.txt
    touch versions.yml
    """
}
