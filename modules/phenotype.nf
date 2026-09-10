
/*
    MODULE : BUILD_PHENOTYPE
    Outil   : Python 3 (stdlib)
    Rôle    : Construire les fichiers de phénotypes pour le GWAS

    GEMMA et PLINK2 attendent des formats différents ET une correspondance
    STRICTE avec l'ordre des individus du fichier .fam :

      - GEMMA  (-p) : une valeur numérique par ligne, dans l'ordre EXACT
                      du .fam, 'NA' pour une valeur manquante.
      - PLINK2 (--pheno) : fichier tabulé avec en-tête '#FID  IID  PHENO1'.

    Source des phénotypes (auto-détectée) :
      1. Le samplesheet lui-même s'il contient une colonne 'phenotype'
         (en-tête CSV : sample,fastq_1,fastq_2,...,phenotype)
      2. Un fichier dédié passé via --phenotype_file :
           - 2 colonnes : <sample> <valeur>
           - 3 colonnes : <FID> <IID> <valeur>
         séparateur virgule, tabulation ou espace.

    Les individus absents de la source reçoivent 'NA' et seront ignorés
    par GEMMA (et exclus du calcul de la kinship).
*/

process BUILD_PHENOTYPE {
    tag "build_phenotype"
    label 'process_low'

    publishDir "${params.outdir}/05_gwas/phenotype", mode: 'copy'

    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple path(bed), path(bim), path(fam)
    path pheno_source

    output:
    path "phenotype.gemma.txt", emit: gemma
    path "phenotype.plink.tsv", emit: plink
    path "phenotype_summary.txt", emit: summary
    path "versions.yml",        emit: versions

    script:
    """
    export PHENO_SOURCE="${pheno_source}"
    export FAM_FILE="${fam}"

    python3 - <<'PYEOF'
    import os, sys, csv, io

    src_path = os.environ["PHENO_SOURCE"]
    fam_path = os.environ["FAM_FILE"]

    MISSING = {"", "na", "nan", "-9", "none", "null", "."}

    def clean(v):
        v = (v or "").strip()
        return "NA" if v.lower() in MISSING else v

    # ── Charger la table sample -> phénotype depuis la source ────────────────
    pheno = {}
    with open(src_path, newline="") as fh:
        text = fh.read()

    first = text.splitlines()[0] if text.splitlines() else ""
    if "," in first:
        delim = ","
    elif "\\t" in first:
        delim = "\\t"
    else:
        delim = None  # espace(s)

    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        rows.append(line.split(delim) if delim else line.split())

    header = [h.strip().lstrip("#").lower() for h in rows[0]]

    if "phenotype" in header and "sample" in header:
        # Cas 1 : samplesheet avec colonne 'phenotype'
        s_i, p_i = header.index("sample"), header.index("phenotype")
        for r in rows[1:]:
            if len(r) > max(s_i, p_i):
                pheno[r[s_i].strip()] = clean(r[p_i])
    else:
        # Cas 2 : fichier dédié, en-tête optionnel
        data = rows
        if header and header[0] in ("sample", "fid", "iid", "id"):
            data = rows[1:]
        for r in data:
            r = [c.strip() for c in r if c.strip() != ""]
            if len(r) >= 3:
                pheno[r[1]] = clean(r[2])   # FID IID valeur -> clé = IID
                pheno.setdefault(r[0], clean(r[2]))
            elif len(r) == 2:
                pheno[r[0]] = clean(r[1])

    # ── Parcourir le .fam et produire les deux formats ──────────────────────
    n_total = n_pheno = 0
    with open(fam_path) as fh, \\
         open("phenotype.gemma.txt", "w") as g, \\
         open("phenotype.plink.tsv", "w") as p:
        p.write("#FID\\tIID\\tPHENO1\\n")
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            fid, iid = parts[0], parts[1]
            val = pheno.get(iid, pheno.get(fid, "NA"))
            n_total += 1
            if val != "NA":
                n_pheno += 1
            g.write(val + "\\n")
            p.write(f"{fid}\\t{iid}\\t{val}\\n")

    with open("phenotype_summary.txt", "w") as s:
        s.write(f"Individus dans le .fam        : {n_total}\\n")
        s.write(f"Individus avec phénotype      : {n_pheno}\\n")
        s.write(f"Individus sans phénotype (NA) : {n_total - n_pheno}\\n")
        s.write(f"Source                        : {src_path}\\n")

    print(open("phenotype_summary.txt").read())

    if n_pheno < 3:
        sys.exit("ERREUR : moins de 3 individus phénotypés — GWAS impossible.")
    PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch phenotype.gemma.txt phenotype.plink.tsv phenotype_summary.txt versions.yml
    """
}
