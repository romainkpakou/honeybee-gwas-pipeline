nextflow.enable.dsl = 2

/*
    MODULE : GEMMA
    Outil   : GEMMA v0.98.5
    Rôle    : GWAS par modèle mixte linéaire (LMM)
    Docker  : quay.io/biocontainers/gemma:0.98.5--hdcf5f25_4

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

    container 'quay.io/biocontainers/gemma:0.98.5--hdcf5f25_4'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    path "output/honeybee.cXX.txt", emit: kinship
    path "output/honeybee.log.txt", emit: log
    path "versions.yml",            emit: versions

    script:
    """
    # GEMMA calcule la matrice de parenté centrée (-gk 1)
    # Input  : fichiers PLINK binaires (bed/bim/fam)
    # Output : matrice N×N dans output/honeybee.cXX.txt
    gemma \\
        -bfile honeybee.pruned \\
        -gk 1 \\
        -o honeybee \\
        -outdir output

    echo "=== Matrice de parenté générée ==="
    echo "Dimensions : \$(wc -l < output/honeybee.cXX.txt) individus"
    echo "Aperçu des 3 premières valeurs de parenté :"
    head -1 output/honeybee.cXX.txt | cut -d' ' -f1-3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gemma: \$(gemma -v 2>&1 | head -1)
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

    container 'quay.io/biocontainers/gemma:0.98.5--hdcf5f25_4'

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
    """
    # ── Lancer le LMM GWAS ────────────────────────────────────────────────────
    # -lmm 4 : calcule les trois tests (Wald + LRT + Score) simultanément
    # -n 1   : utilise le premier phénotype du fichier
    # -k     : chemin vers la matrice de parenté calculée précédemment
    gemma \\
        -bfile honeybee.pruned \\
        -k ${kinship} \\
        -lmm 4 \\
        -n 1 \\
        -o honeybee \\
        -outdir output

    # ── Post-traitement des résultats ─────────────────────────────────────────
    python3 - <<'PYEOF'
import pandas as pd
import numpy as np
import scipy.stats as stats

# Charger les résultats GEMMA
df = pd.read_csv('output/honeybee.assoc.txt', sep='\t')
n_snps = len(df)

print(f"Nombre total de SNPs testés : {n_snps:,}")

# ── Calcul du facteur d'inflation génomique (lambda GC) ──────────────────────
# Lambda mesure si les p-values sont globalement inflées
# Lambda = 1.0 : pas d'inflation (idéal)
# Lambda > 1.1 : inflation probable = problème de stratification
chisq_obs  = stats.chi2.ppf(1 - df['p_wald'].dropna(), df=1)
lambda_gc  = np.median(chisq_obs) / stats.chi2.ppf(0.5, df=1)
print(f"Facteur d'inflation génomique λ = {lambda_gc:.4f}")
if lambda_gc > 1.1:
    print("ATTENTION : λ > 1.1 — vérifier la correction pour stratification")
else:
    print("λ acceptable — pas d'inflation détectée")

# ── Seuils de significativité ─────────────────────────────────────────────────
bonferroni = 0.05 / n_snps
suggestif  = 1e-5

df['significant_bonferroni'] = df['p_wald'] < bonferroni
df['significant_suggestive'] = df['p_wald'] < suggestif

n_sig = df['significant_bonferroni'].sum()
n_sug = df['significant_suggestive'].sum()

print(f"Seuil Bonferroni (p < {bonferroni:.2e}) : {n_sig} SNPs significatifs")
print(f"Seuil suggestif  (p < {suggestif:.0e}) : {n_sug} SNPs suggestifs")

# ── Top 10 SNPs les plus significatifs ───────────────────────────────────────
print("\\n=== Top 10 SNPs les plus significatifs ===")
top10 = df.nsmallest(10, 'p_wald')[['chr','ps','rs','af','beta','se','p_wald']]
print(top10.to_string(index=False))

# ── Sauvegarder les résultats annotés ─────────────────────────────────────────
df['lambda_gc']             = lambda_gc
df['bonferroni_threshold']  = bonferroni
df['suggestive_threshold']  = suggestif
df.to_csv('output/honeybee.assoc.annotated.txt', sep='\t', index=False)
print("\\nRésultats annotés sauvegardés.")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gemma: \$(gemma -v 2>&1 | head -1)
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
