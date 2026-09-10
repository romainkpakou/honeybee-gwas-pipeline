
/*
    MODULE : bcftools
    Outil   : bcftools v1.19
    Rôle    : Filtrage final du VCF et statistiques
    Docker  : quay.io/biocontainers/bcftools:1.19--h8b25389_1

    DEUX PROCESS DANS CE MODULE :

    1. BCFTOOLS_FILTER
       Après GATK VariantFiltration, le VCF contient TOUS les variants
       avec une colonne FILTER indiquant s'ils passent ou non les filtres.
       bcftools filter extrait uniquement les variants PASS bialléliques.

       Filtres appliqués :
       --apply-filters PASS    : garde uniquement les variants PASS
       --min-alleles 2         : exclut les monomorphes (pas de variant)
       --max-alleles 2         : garde uniquement les variants bialléliques
       --types snps            : garde uniquement les SNPs (pas les indels)

       Pourquoi biallélique uniquement ?
       PLINK2 et GEMMA travaillent sur des variants bialléliques.
       Les variants multialléliques compliquent les analyses GWAS
       et sont rares chez l'abeille.

    2. BCFTOOLS_STATS
       Calcule des statistiques complètes sur le VCF filtré :
       - Nombre total de SNPs
       - Ratio Ts/Tv (transitions/transversions)
         Attendu ~2.0 pour un génome de bonne qualité
         Si Ts/Tv < 1.5 : trop de faux positifs
       - Distribution des fréquences alléliques (MAF spectrum)
       - Taux de génotypage manquant
       Collecté par MultiQC pour le rapport final.
*/

// ── Process 1 : Filtrage PASS + SNPs bialléliques ─────────────────────────────
process BCFTOOLS_FILTER {
    tag "bcftools_filter"
    label 'process_medium'

    publishDir "${params.outdir}/03_variants/filtered", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("cohort.pass.vcf.gz"),
          path("cohort.pass.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                  emit: versions

    script:
    """
    # Garder uniquement les variants PASS, bialléliques et SNPs
    bcftools view \\
        --apply-filters PASS \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        --types snps \\
        -Oz -o cohort.pass.vcf.gz \\
        ${vcf}

    # Indexer le VCF filtré
    bcftools index --tbi cohort.pass.vcf.gz

    # Afficher un résumé rapide
    echo "=== Variants après filtrage ===" 
    bcftools stats cohort.pass.vcf.gz | grep "^SN"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch cohort.pass.vcf.gz cohort.pass.vcf.gz.tbi
    touch versions.yml
    """
}

// ── Process 2 : Statistiques du VCF ──────────────────────────────────────────
process BCFTOOLS_STATS {
    tag "bcftools_stats"
    label 'process_low'

    publishDir "${params.outdir}/03_variants/stats", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_1'

    input:
    tuple path(vcf), path(tbi)

    output:
    path "cohort.stats.txt", emit: stats
    path "versions.yml",     emit: versions

    script:
    """
    bcftools stats ${vcf} > cohort.stats.txt

    # Afficher les statistiques clés dans les logs
    echo "=== Statistiques clés ==="
    grep "^SN" cohort.stats.txt | head -20
    echo ""
    echo "=== Ratio Ts/Tv ==="
    grep "Ts/Tv" cohort.stats.txt | head -5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch cohort.stats.txt
    touch versions.yml
    """
}
