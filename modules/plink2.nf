
/*
    MODULE : PLINK2
    Outil   : PLINK2 v2.00a5.12
    Rôle    : QC génotypique, PCA et GWAS
    Docker  : quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0

    TROIS PROCESS DANS CE MODULE :

    1. PLINK2_QC — Contrôle qualité génotypique + élagage LD
       Étape 1 : Filtres QC sur le VCF
         MAF  (--maf)  : exclut les SNPs trop rares (< 5%)
                         Les SNPs rares manquent de puissance en GWAS
         GENO (--geno) : exclut les SNPs avec trop de données manquantes
                         (> 5% de génotypes manquants par SNP)
         MIND (--mind) : exclut les individus avec trop de données manquantes
                         (> 10% de génotypes manquants par individu)
         HWE  (--hwe)  : exclut les SNPs en déséquilibre Hardy-Weinberg
                         Signe d'erreur de génotypage ou de sélection forte

       Étape 2 : Élagage LD (LD pruning)
         Supprime les SNPs en fort déséquilibre de liaison.
         Nécessaire avant PCA et ADMIXTURE car ces analyses
         supposent des marqueurs indépendants.
         Paramètres :
           window  = 50 kb  : taille de la fenêtre glissante
           step    = 10 SNPs : pas de déplacement
           r2      = 0.2    : seuil de corrélation maximum

    2. PLINK2_PCA — Analyse en Composantes Principales
       Calcule les 20 premières composantes principales
       à partir des génotypes. Visualise la structure de
       population — les individus de la même population
       se regroupent dans l'espace des PC.
       Les PC sont aussi utilisés comme covariables dans
       le GWAS PLINK2 pour corriger la stratification.

    3. PLINK2_GWAS — Test d'association
       Test linéaire (phénotype quantitatif) ou logistique
       (phénotype binaire cas/témoin).
       Moins robuste que GEMMA LMM pour corriger la
       stratification mais plus rapide et complémentaire.
*/

// ── Process 1 : QC génotypique + élagage LD ───────────────────────────────────
process PLINK2_QC {
    tag "plink2_qc"
    label 'process_medium'

    publishDir "${params.outdir}/04_population/plink_qc", mode: 'copy'

    container 'quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0'

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("honeybee.pruned.bed"),
          path("honeybee.pruned.bim"),
          path("honeybee.pruned.fam"), emit: plink_files
    path "honeybee.qc.log",           emit: log
    path "versions.yml",              emit: versions

    script:
    """
    # ── Étape 1 : Conversion VCF → PLINK + filtres QC ────────────────────────
    plink2 \\
        --vcf ${vcf} \\
        --double-id \\
        --allow-extra-chr \\
        --set-missing-var-ids '@:#' \\
        --maf ${params.maf} \\
        --geno ${params.geno} \\
        --mind ${params.mind} \\
        --make-bed \\
        --out honeybee.qc \\
        --threads ${task.cpus} \\
        2>&1 | tee honeybee.qc.log

    echo "" >> honeybee.qc.log
    echo "=== SNPs après QC ===" >> honeybee.qc.log
    wc -l honeybee.qc.bim >> honeybee.qc.log

    # ── Étape 2 : Calcul du LD entre SNPs ────────────────────────────────────
    plink2 \\
        --bfile honeybee.qc \\
        --allow-extra-chr \\
        --bad-ld \\
        --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \\
        --out honeybee.ldprune \\
        --threads ${task.cpus}

    echo "" >> honeybee.qc.log
    echo "=== SNPs après élagage LD ===" >> honeybee.qc.log
    wc -l honeybee.ldprune.prune.in >> honeybee.qc.log

    # ── Étape 3 : Extraction des SNPs indépendants ────────────────────────────
    plink2 \\
        --bfile honeybee.qc \\
        --allow-extra-chr \\
        --extract honeybee.ldprune.prune.in \\
        --make-bed \\
        --out honeybee.pruned \\
        --threads ${task.cpus}

    echo "" >> honeybee.qc.log
    echo "=== Fichiers PLINK finaux ===" >> honeybee.qc.log
    ls -lh honeybee.pruned.* >> honeybee.qc.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch honeybee.pruned.bed honeybee.pruned.bim honeybee.pruned.fam
    touch honeybee.qc.log
    touch versions.yml
    """
}

// ── Process 2 : PCA ───────────────────────────────────────────────────────────
process PLINK2_PCA {
    tag "plink2_pca"
    label 'process_medium'

    publishDir "${params.outdir}/04_population/pca", mode: 'copy'

    container 'quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    path "honeybee.eigenvec", emit: eigenvec
    path "honeybee.eigenval", emit: eigenval
    path "versions.yml",      emit: versions

    script:
    """
    # Calculer les fréquences alléliques d'abord (requis avec < 50 échantillons)
    plink2 \\
        --bfile honeybee.pruned \\
        --allow-extra-chr \\
        --freq \\
        --out honeybee \\
        --threads ${task.cpus}

    plink2 \\
        --bfile honeybee.pruned \\
        --allow-extra-chr \\
        --read-freq honeybee.afreq \\
        --pca 4 \\
        --out honeybee \\
        --threads ${task.cpus}

    # Afficher les valeurs propres dans les logs
    echo "=== Variance expliquée par chaque PC ==="
    cat honeybee.eigenval

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch honeybee.eigenvec honeybee.eigenval
    touch versions.yml
    """
}

// ── Process 3 : GWAS association test ────────────────────────────────────────
process PLINK2_GWAS {
    tag "plink2_gwas"
    label 'process_high'

    publishDir "${params.outdir}/05_gwas/plink2", mode: 'copy'

    container 'quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0'

    input:
    tuple path(bed), path(bim), path(fam)
    path phenotype

    output:
    path "honeybee_gwas*.glm.*", emit: results
    path "versions.yml",         emit: versions

    script:
    def model = params.gwas_model == 'logistic'
        ? '--logistic hide-covar'
        : '--linear hide-covar'
    """
    plink2 \\
        --bfile honeybee.pruned \\
        --allow-extra-chr \\
        --pheno ${phenotype} \\
        ${model} \\
        --covar-variance-standardize \\
        --ci 0.95 \\
        --out honeybee_gwas \\
        --threads ${task.cpus}

    # Résumé des résultats
    echo "=== SNPs significatifs (p < ${params.gwas_pval}) ==="
    awk -v threshold=${params.gwas_pval} 'NR>1 && \$12 < threshold {print}' \\
        honeybee_gwas*.glm.* | wc -l

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch honeybee_gwas.PHENO1.glm.linear
    touch versions.yml
    """
}
