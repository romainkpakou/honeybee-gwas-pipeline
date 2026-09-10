#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# Téléchargement des données publiques WGS — Apis mellifera mellifera
# Author  : Romain KPAKOU
# Dataset : Wallberg et al. 2019 (PRJNA473480)
# Génome  : Amel_HAv3.1 (GCF_003254395.2)
# Usage   : bash bin/download_data.sh [nombre_échantillons]
# Défaut  : 10 échantillons
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

N=${1:-10}
REF="data/reference"
SMP="data/samples"

echo "============================================================"
echo "  honeybee-gwas-pipeline — Téléchargement des données"
echo "  Génome : Amel_HAv3.1"
echo "  Dataset : PRJNA473480 (Wallberg et al. 2019)"
echo "  Échantillons : ${N}"
echo "============================================================"

mkdir -p ${REF} ${SMP}

# ── 1. Génome de référence ────────────────────────────────────────────────────
echo ""
echo "[1/3] Téléchargement du génome Amel_HAv3.1..."
if [ ! -f "${REF}/Amel_HAv3.1.fa" ]; then
    wget -q --show-progress \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz" \
        -O "${REF}/Amel_HAv3.1.fa.gz"
    gunzip "${REF}/Amel_HAv3.1.fa.gz"
    echo "    Génome téléchargé : ${REF}/Amel_HAv3.1.fa"
else
    echo "    Génome déjà présent. Étape ignorée."
fi

# ── 2. Annotation SnpEff ──────────────────────────────────────────────────────
echo ""
echo "[2/3] Téléchargement de l'annotation GFF..."
if [ ! -f "${REF}/Amel_HAv3.1.gff.gz" ]; then
    wget -q --show-progress \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz" \
        -O "${REF}/Amel_HAv3.1.gff.gz"
    echo "    Annotation téléchargée."
else
    echo "    Annotation déjà présente. Étape ignorée."
fi

# ── 3. Données WGS SRA ────────────────────────────────────────────────────────
echo ""
echo "[3/3] Téléchargement de ${N} échantillons WGS depuis PRJNA473480..."

# Accessions SRA — Apis mellifera mellifera Europe occidentale
SAMPLES=(
    "SRR7190186"  # AMM_FR_001 — France
    "SRR7190187"  # AMM_FR_002 — France
    "SRR7190188"  # AMM_FR_003 — France
    "SRR7190189"  # AMM_UK_001 — Royaume-Uni
    "SRR7190190"  # AMM_UK_002 — Royaume-Uni
    "SRR7190191"  # AMM_NO_001 — Norvège
    "SRR7190192"  # AMM_NO_002 — Norvège
    "SRR7190193"  # AMM_ES_001 — Espagne
    "SRR7190194"  # AMM_ES_002 — Espagne
    "SRR7190195"  # AMM_DE_001 — Allemagne
)

COUNTRIES=("FR" "FR" "FR" "UK" "UK" "NO" "NO" "ES" "ES" "DE")

for i in $(seq 0 $((N-1))); do
    SRR=${SAMPLES[$i]}
    COUNTRY=${COUNTRIES[$i]}
    SAMPLE_NAME="AMM_${COUNTRY}_$(printf '%03d' $((i+1)))"

    echo "    Téléchargement ${SRR} (${SAMPLE_NAME})..."

    if [ ! -f "${SMP}/${SRR}_1.fastq.gz" ]; then
        prefetch ${SRR} --output-directory ${SMP}
        fasterq-dump ${SMP}/${SRR}/${SRR}.sra \
            --outdir ${SMP} \
            --split-files \
            --threads 4
        gzip ${SMP}/${SRR}_1.fastq ${SMP}/${SRR}_2.fastq 2>/dev/null || true
        rm -rf ${SMP}/${SRR}/
        echo "    ${SRR} → OK"
    else
        echo "    ${SRR} déjà présent. Étape ignorée."
    fi
done

# ── 4. Génération du samplesheet ──────────────────────────────────────────────
echo ""
echo "Génération du samplesheet..."
SHEET="samplesheet.csv"
# La colonne 'phenotype' est laissée vide : la renseigner (trait quantitatif
# ou binaire 0/1) pour activer l'étape GWAS. Vide => GWAS ignoré.
echo "sample,fastq_1,fastq_2,sex,population,phenotype" > ${SHEET}

for i in $(seq 0 $((N-1))); do
    SRR=${SAMPLES[$i]}
    COUNTRY=${COUNTRIES[$i]}
    SAMPLE_NAME="AMM_${COUNTRY}_$(printf '%03d' $((i+1)))"
    echo "${SAMPLE_NAME},${SMP}/${SRR}_1.fastq.gz,${SMP}/${SRR}_2.fastq.gz,unknown,AMM," >> ${SHEET}
done

echo ""
echo "============================================================"
echo "  Téléchargement terminé !"
echo ""
echo "  Génome    : ${REF}/Amel_HAv3.1.fa"
echo "  Échantillons : ${SMP}/"
echo "  Samplesheet  : ${SHEET}"
echo ""
echo "  Lancer le pipeline :"
echo "  nextflow run main.nf \\"
echo "      -params-file conf/params.yml \\"
echo "      -profile docker"
echo "============================================================"
