# honeybee-gwas-pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/container-Docker-blue.svg)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Pipeline WGS et GWAS de bout en bout pour la génomique des populations
d'*Apis mellifera mellifera* (abeille noire)**

> Auteur : Romain KPAKOU | Master 2 Bioinformatique, Nantes Université
> GitHub : [github.com/romainkpakou](https://github.com/romainkpakou)

---

## Vue d'ensemble

`honeybee-gwas-pipeline` est un pipeline Nextflow DSL2 entièrement automatisé
pour l'analyse WGS et les études d'association pangénomique (GWAS) appliquées
à *Apis mellifera mellifera*, la sous-espèce native d'Europe sous pression
de conservation. Il suit les GATK Best Practices pour le variant calling et
implémente des méthodes modernes de génomique des populations.

---

## Étapes du pipeline

FASTQ bruts
↓
ÉTAPE 1 — QC FastQC · fastp · MultiQC
↓
ÉTAPE 2 — Alignement BWA-MEM2 · SAMtools
↓
ÉTAPE 3 — BAM Picard MarkDuplicates
↓
ÉTAPE 4-5 — Variants GATK HaplotypeCaller · GenomicsDBImport · GenotypeGVCFs
↓
ÉTAPE 6 — Filtrage GATK VariantFiltration · bcftools
↓
ÉTAPE 7 — Annotation SnpEff (base construite localement depuis le GFF3)
↓
ÉTAPE 8 — Pop. gen. PLINK2 · ADMIXTURE · vcftools (FST · LD · π)
↓
ÉTAPE 9 — GWAS GEMMA LMM · PLINK2
↓
ÉTAPE 10 — Figures Manhattan · QQ · PCA · Admixture · FST · LD decay
↓
ÉTAPE 11 — Rapport R Markdown HTML


---

## Prérequis

| Outil | Version | Installation |
|---|---|---|
| Nextflow | ≥ 23.04.0 | `curl -s https://get.nextflow.io \| bash` |
| Docker | any | [docs.docker.com](https://docs.docker.com) |
| Java | ≥ 11 | `sudo apt install default-jdk` |

Tous les outils bioinformatiques sont téléchargés automatiquement
via Docker — aucune installation manuelle requise.

---

## Démarrage rapide

### 1. Télécharger les données publiques

```bash
# Génome Amel_HAv3.1 + 10 échantillons WGS (PRJNA473480)
bash bin/download_data.sh 10
```

### 2. Configurer les paramètres

```bash
# Éditer conf/params.yml selon vos besoins
nano conf/params.yml
```

### 3. Lancer le pipeline

```bash
# Exécution locale avec Docker
nextflow run main.nf \
    -params-file conf/params.yml \
    -profile docker

# Cluster HPC avec SLURM + Singularity
nextflow run main.nf \
    -params-file conf/params.yml \
    -profile slurm

# Test rapide (données minimales)
nextflow run main.nf -profile test,docker
```

---

## Format du samplesheet

```csv
sample,fastq_1,fastq_2,sex,population,phenotype
AMM_FR_001,data/samples/SRR7190186_1.fastq.gz,data/samples/SRR7190186_2.fastq.gz,unknown,AMM,62.5
AMM_UK_001,data/samples/SRR7190189_1.fastq.gz,data/samples/SRR7190189_2.fastq.gz,unknown,AMM,45.0
```

La colonne **`phenotype`** est optionnelle :

- **présente avec au moins une valeur** → l'étape GWAS (GEMMA LMM + PLINK2 +
  Manhattan/QQ) est activée automatiquement ; valeur vide ou `NA` = individu
  exclu de l'association.
- **absente** (ou `--phenotype_file` non fourni) → le GWAS est ignoré, le reste
  du pipeline s'exécute normalement.

Alternative : `--phenotype_file` pointant vers un fichier `sample,valeur` ou
`FID IID valeur`.

---

## Résultats

results/
├── 01_qc/ FastQC, fastp, MultiQC
├── 02_alignment/ BAM triés et dédupliqués
├── 03_variants/ VCF filtré + VCF annoté SnpEff + rapport d'effets
├── 04_population/ PCA, ADMIXTURE, FST, LD decay
├── 05_gwas/ GEMMA LMM, PLINK2, Manhattan + QQ plots
├── 06_report/ Rapport HTML complet (R Markdown)
└── pipeline_info/ Nextflow report, timeline, trace, DAG

---

## Paramètres principaux

| Paramètre | Défaut | Description |
|---|---|---|
| `--maf` | 0.05 | Seuil fréquence allélique mineure |
| `--geno` | 0.05 | Taux max de génotypage manquant par SNP |
| `--mind` | 0.10 | Taux max de génotypage manquant par individu |
| `--hwe` | 1e-6 | Seuil Hardy-Weinberg |
| `--ld_r2` | 0.2 | Seuil r² pour l'élagage LD |
| `--admixture_k` | `2,3,4,5` | Valeurs de K à tester |
| `--gff` | `data/reference/Amel_HAv3.1.gff.gz` | Annotation GFF3 → active SnpEff (vide = désactivé) |
| `--phenotype_file` | `null` | Fichier phénotypes (sinon colonne `phenotype` du samplesheet) |
| `--gwas_model` | `lmm` | Modèle GWAS : `lmm` (GEMMA) ou `logistic` (trait binaire) |
| `--gwas_pval` | 1e-6 | Seuil de significativité GWAS |

Tous les paramètres sont modifiables dans `conf/params.yml`.

---

## Données utilisées

| Ressource | Accession | Référence |
|---|---|---|
| Données WGS | PRJNA473480 | Wallberg et al. (2019) *Nat Ecol Evol* |
| Génome référence | GCF_003254395.2 | Amel_HAv3.1 — 16 chromosomes |
| Base SnpEff | Apis_mellifera | Cingolani et al. (2012) |

---

## Citation
Romain KPAKOU (2026). honeybee-gwas-pipeline: End-to-end WGS/GWAS pipeline
for Apis mellifera mellifera population genomics. v1.0.0.
https://github.com/romainkpakou/honeybee-gwas-pipeline

---

## Licence

MIT — voir [LICENSE](LICENSE)

---

## Contact

**Romain KPAKOU**
Master 2 Bioinformatique, Biostatistique & Biologie Computationnelle
Nantes Université
kpakouromain@gmail.com | [github.com/romainkpakou](https://github.com/romainkpakou)
