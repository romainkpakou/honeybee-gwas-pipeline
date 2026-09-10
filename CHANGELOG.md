# Changelog

## [Non publié]

### Ajouté
- `BUILD_PHENOTYPE` : nouveau module qui construit les fichiers de phénotypes
  aux formats GEMMA (`-p`) et PLINK2 (`--pheno`), alignés sur l'ordre du `.fam`
- Activation automatique du GWAS via une colonne `phenotype` dans le samplesheet
  (en plus de `--phenotype_file`)
- `PLINK2_QC` expose désormais le jeu de SNPs QC non élagué (`qc_files`), utilisé
  pour les tests d'association (la kinship reste calculée sur le jeu élagué LD)
- `GWAS_REPORT` réactivé : rapport HTML R Markdown (résumé QC, PCA, ADMIXTURE,
  Manhattan, top SNPs), généré si le GWAS est activé
- **Annotation SnpEff réactivée** (étape 7), en mode hors-ligne : `SNPEFF_BUILD`
  construit la base à partir du génome + GFF3 NCBI (assemblage identique au
  variant calling), `SNPEFF_ANNOTATE` produit le VCF annoté + rapport HTML/CSV
  (agrégé par MultiQC), `SNPEFF_COMPRESS` bgzip + tabix. Activée par `--gff`.

### Corrigé
- `GWAS_REPORT` : conteneur `quay.io/biocontainers/r-base` (sans pandoc) →
  `rocker/tidyverse:4.3.1` (rmarkdown + pandoc inclus) ; suppression du
  `install.packages` à l'exécution ; template `.Rmd` passé en entrée du process ;
  dépendances `kableExtra`/`gridExtra` retirées du template
- Branche GWAS jamais fonctionnelle : `GEMMA_KINSHIP`/`GEMMA_LMM` ne recevaient
  aucun phénotype (colonne 6 du `.fam` à `-9`) → passage explicite via `-p`
- Image Docker GEMMA obsolète (`0.98.5--hdcf5f25_4` retirée de quay.io)
  → `0.98.5--h38cc83e_1`
- `GEMMA_LMM` : post-traitement réécrit en awk (dépendance pandas/numpy/scipy
  absente du conteneur biocontainers)
- `GEMMA_LMM` demandait 32 GB de RAM (irréaliste) → 8 GB
- `PLINK2_GWAS` : `--glm allow-no-covars` (l'ancien `--covar-variance-standardize`
  échouait sans fichier de covariables) + détection de la colonne P par en-tête
- `PLOT_MANHATTAN` : la colonne `chr` de GEMMA (accession NCBI `NC_037638.1`)
  est mappée vers un index de chromosome séquentiel pour l'axe X

## [1.0.0] - 2026-08-14

### Ajouté
- Pipeline WGS complet : FastQC, fastp, BWA-MEM2, SAMtools, Picard, GATK4
- Génotypage joint : GATK GenomicsDBImport + GenotypeGVCFs
- Filtrage variants : GATK VariantFiltration + bcftools
- Annotation fonctionnelle : SnpEff (base de données Apis_mellifera)
- Génétique des populations : PLINK2 QC + PCA, ADMIXTURE (K=2..5)
- Statistiques populationnelles : vcftools FST, LD decay, diversité π
- GWAS : GEMMA kinship + LMM, PLINK2 association
- Visualisations R publication-ready : Manhattan, QQ, PCA, Admixture, FST, LD
- Rapport R Markdown HTML/PDF automatique
- Profils Docker, Singularity, SLURM, test
- Script de téléchargement automatique des données SRA (PRJNA473480)
- Documentation complète (README, tools_and_versions.md, params.yml)
