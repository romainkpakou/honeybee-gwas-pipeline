
/*
    MODULE : GEMMA
    Outil   : GEMMA v0.98.5
    Rôle    : GWAS par modèle mixte linéaire (LMM)
    Docker  : quay.io/biocontainers/gemma:0.98.5--h38cc83e_1

    DEUX PROCESS DANS CE MODULE :

    1. GEMMA_KINSHIP — Matrice de parenté génomique (kinship)
       Calcule la matrice de parenté centrée (centered relatedness
       matrix) entre tous les individus à partir de l'ensemble
       des SNPs génotypés.
       Cette matrice capture :
         - La structure de population (différences entre groupes)
         - La parenté cryptique (individus apparentés non déclarés)
         - La consanguinité (inbreeding)

       Deux types de matrices disponibles :
         -gk 1 : matrice de parenté centrée (recommandée)
         -gk 2 : matrice de parenté standardisée

       Le fichier .cXX.txt contient une matrice N×N
       (N = nombre d'individus) de coefficients de parenté.

    2. GEMMA_LMM — Linear Mixed Model GWAS
       Pour chaque SNP, ajuste un modèle :
         y = Xβ + Zu + ε
       où :
         y = vecteur des phénotypes
         X = génotypes au SNP testé
         β = effet du SNP (ce qu'on cherche)
         Z = matrice de conception
         u = effets aléatoires (parenté) ~ N(0, σ²_g × K)
         ε = résidus ~ N(0, σ²_e × I)

       Calcule trois statistiques de test :
         Wald test    : le plus utilisé en pratique
         LRT (likelihood ratio test) : plus précis mais plus lent
         Score test   : rapide, bonne approximation

       SEUILS DE SIGNIFICATIVITÉ :
         Bonferroni strict : p < 0.05 / N_SNPs
         Suggestif        : p < 1e-5
         Exploratoire     : p < 1e-4

       FACTEUR D'INFLATION GÉNOMIQUE (λ) :
         Calculé sur les p-values de l'ensemble des SNPs.
         λ proche de 1.0 = pas de stratification résiduelle
         λ > 1.1 = inflation possible = vérifier la correction
*/

// ── Process 1 : Matrice de parenté (kinship) ──────────────────────────────────
process GEMMA_KINSHIP {
    tag "gemma_kinship"
    label 'process_high'

    publishDir "${params.outdir}/05_gwas/kinship", mode: 'copy'

    container 'quay.io/biocontainers/gemma:0.98.5--h38cc83e_1'

    input:
    tuple path(bed), path(bim), path(fam)
    path phenotype

    output:
    path "output/honeybee.cXX.txt", emit: kinship
    path "output/honeybee.log.txt", emit: log
    path "versions.yml",            emit: versions

    script:
    // Préfixe PLINK déduit du .bed (kinship calculée sur les SNPs élagués LD)
    def prefix = bed.name - ~/\.bed$/
    """
    # GEMMA calcule la matrice de parenté centrée (-gk 1)
    # Input  : fichiers PLINK binaires (bed/bim/fam) + phénotypes (-p)
    # -p : GEMMA exclut du calcul les individus au phénotype manquant
    # Output : matrice N×N dans output/honeybee.cXX.txt
    gemma \\
        -bfile ${prefix} \\
        -p ${phenotype} \\
        -gk 1 \\
        -o honeybee \\
        -outdir output

    echo "=== Matrice de parenté générée ==="
    echo "Dimensions : \$(wc -l < output/honeybee.cXX.txt) individus"
    echo "Aperçu des 3 premières valeurs de parenté :"
    head -1 output/honeybee.cXX.txt | cut -d' ' -f1-3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gemma: \$(gemma -v 2>&1 | grep -im1 version || echo 0.98.5)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output
    touch output/honeybee.cXX.txt
    touch output/honeybee.log.txt
    touch versions.yml
    """
}

// ── Process 2 : LMM GWAS ──────────────────────────────────────────────────────
process GEMMA_LMM {
    tag "gemma_lmm"
    label 'process_high'

    publishDir "${params.outdir}/05_gwas/gemma", mode: 'copy'

    container 'quay.io/biocontainers/gemma:0.98.5--h38cc83e_1'

    input:
    tuple path(bed), path(bim), path(fam)
    path kinship
    path phenotype

    output:
    path "output/honeybee.assoc.txt",           emit: results
    path "output/honeybee.assoc.annotated.txt", emit: annotated
    path "output/honeybee.log.txt",             emit: log
    path "versions.yml",                         emit: versions

    script:
    // Préfixe PLINK déduit du .bed — l'association est testée sur le jeu QC
    // (non élagué LD) tandis que la kinship vient du jeu élagué.
    def prefix = bed.name - ~/\.bed$/
    """
    # ── Lancer le LMM GWAS ────────────────────────────────────────────────────
    # -p     : fichier de phénotypes (une valeur par individu, ordre du .fam)
    # -lmm 4 : calcule les trois tests (Wald + LRT + Score) simultanément
    # -n 1   : utilise le premier phénotype du fichier
    # -k     : chemin vers la matrice de parenté calculée précédemment
    gemma \\
        -bfile ${prefix} \\
        -p ${phenotype} \\
        -k ${kinship} \\
        -lmm 4 \\
        -n 1 \\
        -o honeybee \\
        -outdir output

    # ── Post-traitement des résultats (awk pur — pas de dépendance Python) ────
    # Le lambda GC est recalculé par les scripts R (Manhattan / QQ). Ici on se
    # limite à annoter chaque SNP avec les seuils de significativité et à
    # afficher un résumé dans les logs Nextflow.
    ASSOC=output/honeybee.assoc.txt
    ANNOT=output/honeybee.assoc.annotated.txt

    N_SNPS=\$(tail -n +2 "\$ASSOC" | wc -l)
    P_COL=\$(head -1 "\$ASSOC" | tr '\\t' '\\n' | grep -nx 'p_wald' | cut -d: -f1)
    BONF=\$(awk -v n="\$N_SNPS" 'BEGIN { if (n > 0) printf "%.6e", 0.05 / n; else print "NA" }')

    echo "SNPs testés               : \$N_SNPS"
    echo "Seuil Bonferroni (0.05/N) : \$BONF"
    echo "Seuil suggestif           : 1e-5"

    awk -v OFS='\\t' -v bonf="\$BONF" -v pc="\$P_COL" '
        NR == 1 { print \$0, "bonferroni_threshold", "suggestive_threshold"; next }
        {
            print \$0, bonf, 1e-5
            if (\$pc != "" && \$pc != "NA" && (\$pc + 0) < (bonf + 0)) nsig++
            if (\$pc != "" && \$pc != "NA" && (\$pc + 0) < 1e-5)        nsug++
        }
        END {
            print "SNPs significatifs (Bonferroni) : " nsig+0 > "/dev/stderr"
            print "SNPs suggestifs (p < 1e-5)      : " nsug+0 > "/dev/stderr"
        }
    ' "\$ASSOC" > "\$ANNOT"

    echo "=== Top 10 SNPs (p_wald croissant) ==="
    head -1 "\$ASSOC"
    tail -n +2 "\$ASSOC" | sort -g -k"\$P_COL","\$P_COL" | head -10

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gemma: \$(gemma -v 2>&1 | grep -im1 version || echo 0.98.5)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output
    touch output/honeybee.assoc.txt
    touch output/honeybee.assoc.annotated.txt
    touch output/honeybee.log.txt
    touch versions.yml
    """
}
