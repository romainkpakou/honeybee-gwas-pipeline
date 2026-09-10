# Changelog

Format : [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/) ·
Versionnage : [SemVer](https://semver.org/lang/fr/).

## [1.0.0] - 2026-09-10

Première version **fonctionnelle de bout en bout** : les étapes GWAS et
annotation, présentes mais jamais exécutées dans l'ébauche 0.9.0, sont
opérationnelles et couvertes par l'intégration continue.

### Ajouté
- **Branche GWAS** opérationnelle : `BUILD_PHENOTYPE` (formats GEMMA `-p` et
  PLINK2 `--pheno`, alignés sur l'ordre du `.fam`) → `GEMMA_KINSHIP` /
  `GEMMA_LMM` → `PLINK2_GWAS` → Manhattan / QQ. Activation automatique par une
  colonne `phenotype` du samplesheet ou `--phenotype_file`.
- **Annotation SnpEff hors-ligne** (étape 7) : `SNPEFF_BUILD` construit la base
  depuis le génome + GFF3 NCBI (assemblage identique au variant calling),
  `SNPEFF_ANNOTATE` produit le VCF annoté + rapport HTML/CSV (agrégé par
  MultiQC), `SNPEFF_COMPRESS` bgzip + tabix. Activée par `--gff`.
- **`GWAS_REPORT`** : rapport HTML R Markdown (QC, PCA, ADMIXTURE, Manhattan,
  top SNPs) via `rocker/tidyverse` (rmarkdown + pandoc).
- **Intégration continue** (`.github/workflows/ci.yml`) : `nextflow lint` +
  run `-stub` complet du workflow sur fixtures versionnées (`test/stub/`),
  sans conteneur ni donnée réelle.
- **Séparation démo / production** : `-profile test` (5 génomes, filtres
  relâchés) ; `conf/params.yml` = modèle de production commenté ;
  `docs/running_a_real_cohort.md` = guide de passage à l'échelle.
- `process.resourceLimits` : plafonne cpus / mémoire / temps de toutes les
  tâches aux valeurs `max_*`.
- Paramètres exposés : `pca_components`, `plink_bad_ld`, `hwe_filter`, `gff`,
  `genome_species` (auparavant codés en dur).
- `PLINK2_QC` expose le jeu de SNPs QC non élagué : l'association est testée
  dessus, la kinship reste calculée sur le jeu élagué LD.

### Corrigé
- `GEMMA_KINSHIP` / `GEMMA_LMM` ne recevaient aucun phénotype (colonne 6 du
  `.fam` à `-9`) → passage explicite via `-p`.
- Image Docker GEMMA obsolète (`0.98.5--hdcf5f25_4`) → `0.98.5--h38cc83e_1`.
- `GEMMA_LMM` : post-traitement réécrit en awk (pandas/numpy/scipy absents du
  conteneur) ; demande mémoire 32 GB → 8 GB.
- `PLINK2_GWAS` : `--glm allow-no-covars` (l'ancien `--covar-variance-standardize`
  échouait sans covariables) ; colonne P repérée par en-tête.
- `PLOT_MANHATTAN` / rapport : la colonne `chr` de GEMMA (accession NCBI
  `NC_037638.1`) est mappée vers un index de chromosome séquentiel.
- `BWA_MEM2_INDEX` : ajout d'un plafond mémoire explicite (OOM sous le profil
  de test).
- `nextflow lint` : 0 avertissement (`it` implicite → paramètres nommés,
  `Channel` → `channel`).

## [0.9.0] - 2026-08-14

Ébauche : pipeline WGS complet et fonctionnel. Les modules GWAS, SnpEff et
rapport existaient mais n'avaient jamais tourné de bout en bout.

### Ajouté
- Pipeline WGS : FastQC, fastp, BWA-MEM2, SAMtools, Picard, GATK4
- Génotypage joint : GATK GenomicsDBImport + GenotypeGVCFs
- Filtrage : GATK VariantFiltration + bcftools
- Génétique des populations : PLINK2 QC + PCA, ADMIXTURE, vcftools (FST, LD, π)
- Modules GWAS (GEMMA, PLINK2), visualisations R, rapport R Markdown — non testés
- Profils Docker, Singularity, SLURM, test
- Script de téléchargement des données SRA (PRJNA473480)
- Documentation (README, tools_and_versions.md, params.yml)
