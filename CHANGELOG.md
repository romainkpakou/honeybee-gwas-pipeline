# Changelog

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
