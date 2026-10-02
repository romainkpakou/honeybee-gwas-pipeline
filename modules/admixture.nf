
/*
    MODULE : ADMIXTURE
    Outil   : ADMIXTURE v1.3.0
    Rôle    : Estimation de la structure de population
              et des proportions d'ascendance
    Docker  : quay.io/biocontainers/admixture:1.3.0--0

    PRINCIPE :
    ADMIXTURE modélise chaque individu comme un mélange de K
    populations ancestrales. Il estime simultanément :
      - Q : matrice des proportions d'ascendance (individus × K)
      - P : matrice des fréquences alléliques ancestrales (SNPs × K)

    CROSS-VALIDATION :
    Le flag --cv calcule l'erreur de validation croisée pour chaque K.
    On teste plusieurs valeurs de K (par défaut 2,3,4,5) et on choisit
    le K qui minimise le CV error — c'est le nombre de populations
    ancestrales qui explique le mieux les données.

    INTERPRÉTATION POUR L'ABEILLE NOIRE :
    K=2 : AMM vs autres sous-espèces
    K=3 : AMM + ligustica + carnica (les 3 sous-espèces principales)
    K=4+ : structure fine au sein des populations AMM

    FICHIERS DE SORTIE :
    honeybee.pruned.2.Q → proportions d'ascendance pour K=2
    honeybee.pruned.3.Q → proportions d'ascendance pour K=3
    honeybee.pruned.2.P → fréquences alléliques ancestrales K=2
    cv_K<k>.txt         → erreur CV de chaque K
    admixture_cv_summary.tsv + best_K.txt → K retenu (ADMIXTURE_BEST_K)
*/

process ADMIXTURE_RUN {
    tag "admixture_K${k}"
    label 'process_high'

    publishDir "${params.outdir}/04_population/admixture", mode: 'copy'

    container 'quay.io/biocontainers/admixture:1.3.0--0'

    input:
    tuple path(bed), path(bim), path(fam)
    each k

    output:
    path "*.Q",          emit: q_files
    path "*.P",          emit: p_files
    path "cv_K${k}.txt", emit: cv_error
    path "versions.yml", emit: versions

    script:
    """
    # ADMIXTURE nécessite que le BED soit dans le répertoire courant
    # et que le BIM utilise des identifiants de chromosomes numériques
    # On crée une copie avec les chromosomes renommés si nécessaire

    # Vérifier le format des chromosomes dans le BIM
    # Amel_HAv3.1 utilise des identifiants NC_XXXXXXX
    # ADMIXTURE requiert des entiers ou X/Y
    awk '{
        gsub(/NC_001566\\.1/, "MT");
        gsub(/NC_037638\\.1/, "1");
        gsub(/NC_037639\\.1/, "2");
        gsub(/NC_037640\\.1/, "3");
        gsub(/NC_037641\\.1/, "4");
        gsub(/NC_037642\\.1/, "5");
        gsub(/NC_037643\\.1/, "6");
        gsub(/NC_037644\\.1/, "7");
        gsub(/NC_037645\\.1/, "8");
        gsub(/NC_037646\\.1/, "9");
        gsub(/NC_037647\\.1/, "10");
        gsub(/NC_037648\\.1/, "11");
        gsub(/NC_037649\\.1/, "12");
        gsub(/NC_037650\\.1/, "13");
        gsub(/NC_037651\\.1/, "14");
        gsub(/NC_037652\\.1/, "15");
        gsub(/NC_037653\\.1/, "16");
        print
    }' ${bim} > honeybee.pruned.bim.tmp
    mv honeybee.pruned.bim.tmp honeybee_admix.bim
    cp ${bed} honeybee_admix.bed
    cp ${fam} honeybee_admix.fam

    # Lancer ADMIXTURE avec cross-validation
    admixture --cv honeybee_admix.bed ${k} -j${task.cpus} | tee admixture_K${k}.log

    # Renommer les fichiers de sortie
    mv honeybee_admix.${k}.Q honeybee.pruned.${k}.Q
    mv honeybee_admix.${k}.P honeybee.pruned.${k}.P

    # Extraire et sauvegarder le CV error
    CV=\$(grep "CV error" admixture_K${k}.log | awk '{print \$NF}')
    echo "K=${k}\tCV_error=\${CV}" > cv_K${k}.txt
    echo "K=${k} CV error = \${CV}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: \$(admixture | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch honeybee.pruned.${k}.Q
    touch honeybee.pruned.${k}.P
    touch cv_K${k}.txt
    touch versions.yml
    """
}

// ── Sélection du meilleur K (erreur de validation croisée minimale) ──────────
// Agrège les cv_K*.txt de tous les jobs ADMIXTURE_RUN : un tableau récapitulatif
// et le K retenu, au lieu d'une lecture manuelle des logs. Le K minimisant
// l'erreur CV est un critère statistique ; l'interprétation biologique
// (sous-espèces attendues) reste à faire, surtout si les erreurs sont proches.
process ADMIXTURE_BEST_K {
    tag "admixture_best_k"
    label 'process_low'

    publishDir "${params.outdir}/04_population/admixture", mode: 'copy'

    container 'quay.io/biocontainers/admixture:1.3.0--0'

    input:
    path cv_files

    output:
    path "admixture_cv_summary.tsv", emit: summary
    path "best_K.txt",               emit: best_k

    script:
    """
    # cv_K*.txt : "K=<k>\\tCV_error=<valeur>" — un fichier par job
    printf "K\\tCV_error\\n" > admixture_cv_summary.tsv
    cat ${cv_files} \\
        | sed -e 's/^K=//' -e 's/CV_error=//' \\
        | awk -F'\\t' '\$2 != ""' \\
        | sort -t\$'\\t' -k1,1n >> admixture_cv_summary.tsv

    if [ "\$(tail -n +2 admixture_cv_summary.tsv | wc -l)" -eq 0 ]; then
        echo "ERREUR : aucune erreur CV exploitable dans les logs ADMIXTURE." >&2
        exit 1
    fi

    tail -n +2 admixture_cv_summary.tsv | sort -t\$'\\t' -k2,2g | head -1 | cut -f1 > best_K.txt
    echo "=== Erreurs de validation croisée ==="
    cat admixture_cv_summary.tsv
    echo "K retenu (CV minimale) : \$(cat best_K.txt)"
    """

    stub:
    """
    printf "K\\tCV_error\\n2\\t0.5\\n" > admixture_cv_summary.tsv
    echo 2 > best_K.txt
    """
}
