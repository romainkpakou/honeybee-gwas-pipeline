/*
    MODULE : vcftools
    Outil   : vcftools v0.1.16
    Rôle    : Statistiques de génétique des populations
    Docker  : quay.io/biocontainers/vcftools:0.1.16--pl5321h9a82719_6

    TROIS ANALYSES DANS CE MODULE :

    1. FST DE WEIR & COCKERHAM (1984)
       Mesure la différenciation génétique entre populations.
       FST = 0.0 : populations identiques (pas de différenciation)
       FST = 0.15 : différenciation modérée
       FST = 1.0 : populations complètement fixées pour des allèles différents

       Interprétation pour l'abeille noire :
       FST élevé entre AMM et ligustica/carnica = bonne différenciation
       = les populations AMM sont génétiquement distinctes
       FST faible entre populations AMM = hybridation possible

       Calculé en fenêtres glissantes (50 kb, pas 10 kb) pour
       détecter des régions sous sélection (îlots de différenciation).

    2. DÉCLIN DU LD (Linkage Disequilibrium decay)
       Mesure la corrélation (r²) entre paires de SNPs en fonction
       de leur distance physique sur le chromosome.
       Chez l'abeille (taux de recombinaison élevé ~19 cM/Mb),
       le LD décline rapidement — typiquement r² < 0.2 à 10-50 kb.
       À comparer avec l'humain : r² < 0.2 à ~250 kb.
       Ce paramètre est crucial pour définir la taille de fenêtre
       optimale dans un GWAS et interpréter les régions candidates.

    3. DIVERSITÉ NUCLÉOTIDIQUE (π — pi)
       Mesure la diversité génétique intra-population.
       = probabilité que deux allèles tirés au hasard soient différents.
       Calculée en fenêtres glissantes de 50 kb.
       Permet de détecter des régions de faible diversité = balayages
       sélectifs (selective sweeps) — signe de sélection positive.
*/

// ── Process 1 : FST entre populations ────────────────────────────────────────
process VCFTOOLS_FST {
    tag "vcftools_fst"
    label 'process_medium'

    publishDir "${params.outdir}/04_population/fst", mode: 'copy'

    container 'quay.io/biocontainers/vcftools:0.1.16--pl5321h9a82719_6'

    input:
    tuple path(vcf), path(tbi)

    output:
    path "*.windowed.weir.fst", emit: fst
    path "versions.yml",        emit: versions

    script:
    """
    # Extraire les listes d'individus par population
    # à partir des identifiants d'échantillons
    # Format attendu : POPULATION_PAYS_NUMERO (ex: AMM_FR_001)
    bcftools query -l ${vcf} | grep "^AMM" > pop_AMM.txt || true
    bcftools query -l ${vcf} | grep -v "^AMM" > pop_OTHER.txt || true

    echo "Population AMM : \$(wc -l < pop_AMM.txt) individus"
    echo "Autres populations : \$(wc -l < pop_OTHER.txt) individus"

    # Calculer le FST en fenêtres glissantes si les deux populations existent
    if [ -s pop_AMM.txt ] && [ -s pop_OTHER.txt ]; then
        vcftools \\
            --gzvcf ${vcf} \\
            --weir-fst-pop pop_AMM.txt \\
            --weir-fst-pop pop_OTHER.txt \\
            --fst-window-size 50000 \\
            --fst-window-step 10000 \\
            --out AMM_vs_OTHER

        echo "=== FST moyen AMM vs autres ==="
        awk 'NR>1 {sum+=\$5; count++} END {print "FST moyen =", sum/count}' \\
            AMM_vs_OTHER.windowed.weir.fst
    else
        echo "Populations insuffisantes pour le FST — création d'un fichier vide"
        echo "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIR_AND_COCKERHAM_FST" > \\
            AMM_vs_OTHER.windowed.weir.fst
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(vcftools --version 2>&1 | head -1 | sed 's/VCFtools - //')
    END_VERSIONS
    """

    stub:
    """
    touch AMM_vs_OTHER.windowed.weir.fst
    touch versions.yml
    """
}

// ── Process 2 : Déclin du LD ──────────────────────────────────────────────────
process VCFTOOLS_LD {
    tag "vcftools_ld"
    label 'process_medium'

    publishDir "${params.outdir}/04_population/ld", mode: 'copy'

    container 'quay.io/biocontainers/vcftools:0.1.16--pl5321h9a82719_6'

    input:
    tuple path(vcf), path(tbi)

    output:
    path "cohort_ld.geno.ld", emit: ld
    path "versions.yml",      emit: versions

    script:
    """
    # Calculer r² entre toutes les paires de SNPs
    # dans une fenêtre de 500 kb
    # --min-r2 0.001 : ignore les corrélations négligeables
    #                  pour réduire la taille du fichier de sortie
    vcftools \\
        --gzvcf ${vcf} \\
        --ld-window-bp 500000 \\
        --min-r2 0.001 \\
        --geno-r2 \\
        --out cohort_ld

    echo "=== Aperçu du déclin du LD ==="
    head -5 cohort_ld.geno.ld

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(vcftools --version 2>&1 | head -1 | sed 's/VCFtools - //')
    END_VERSIONS
    """

    stub:
    """
    touch cohort_ld.geno.ld
    touch versions.yml
    """
}

// ── Process 3 : Diversité nucléotidique (π) ──────────────────────────────────
process VCFTOOLS_PI {
    tag "vcftools_pi"
    label 'process_medium'

    publishDir "${params.outdir}/04_population/diversity", mode: 'copy'

    container 'quay.io/biocontainers/vcftools:0.1.16--pl5321h9a82719_6'

    input:
    tuple path(vcf), path(tbi)

    output:
    path "cohort_diversity.windowed.pi", emit: pi
    path "versions.yml",                 emit: versions

    script:
    """
    # Calculer la diversité nucléotidique π en fenêtres de 50 kb
    vcftools \\
        --gzvcf ${vcf} \\
        --window-pi 50000 \\
        --window-pi-step 10000 \\
        --out cohort_diversity

    echo "=== Diversité nucléotidique moyenne ==="
    awk 'NR>1 {sum+=\$5; count++} END {print "π moyen =", sum/count}' \\
        cohort_diversity.windowed.pi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(vcftools --version 2>&1 | head -1 | sed 's/VCFtools - //')
    END_VERSIONS
    """

    stub:
    """
    touch cohort_diversity.windowed.pi
    touch versions.yml
    """
}
