# Tools & Dependencies — honeybee-gwas-pipeline
Author  : Romain KPAKOU
Updated : August 2026

## Workflow manager

### Nextflow DSL2 (v26.04.6)
**Rôle :** Orchestrateur du pipeline. Nextflow permet de chaîner
les étapes bioinformatiques en modules indépendants, de gérer
automatiquement la parallélisation, les reprises en cas d'erreur,
et le déploiement sur HPC ou cloud. DSL2 permet la modularité :
chaque outil est un module réutilisable.
**Pourquoi DSL2 et pas Snakemake ?** Nextflow est plus adapté
aux pipelines HPC/cloud et est le standard dans les pipelines
NGS professionnels (nf-core).

---

## Système & Environnement

### Java OpenJDK 25.0.2
**Rôle :** Requis par Nextflow et GATK4 qui sont des applications
Java. Sans Java, aucun des deux ne peut tourner.

### Conda/Mamba (v2.5.0)
**Rôle :** Gestionnaire d'environnements et de paquets.
Mamba est une réimplémentation de Conda en C++ — il est
5 à 10 fois plus rapide pour résoudre les dépendances.
Utilisé pour créer l'environnement `honeybee-gwas` qui
isole tous les outils du pipeline.

### Docker (v29.1.3)
**Rôle :** Conteneurisation. Chaque outil du pipeline tourne
dans son propre conteneur Docker — cela garantit la
reproductibilité exacte des résultats, quelle que soit
la machine utilisée. Les conteneurs sont automatiquement
téléchargés depuis quay.io/biocontainers.

### Git (v2.43.0)
**Rôle :** Versioning du code. Toutes les modifications du
pipeline sont tracées. Permet la collaboration et le
déploiement depuis GitHub directement avec Nextflow.

---

## Contrôle qualité des reads

### FastQC
**Rôle :** Analyse la qualité des reads FASTQ bruts.
Génère un rapport HTML avec les distributions de qualité
par position, le contenu en GC, les séquences sur-représentées
et les adaptateurs détectés.
**Quand :** Avant ET après le trimming pour comparer.

### fastp
**Rôle :** Trimming et filtrage des reads. Supprime les
adaptateurs Illumina, filtre les reads de mauvaise qualité
(Phred < 20), élimine les reads trop courts (< 50 bp).
Plus rapide que Trimmomatic et génère son propre rapport QC.

### MultiQC
**Rôle :** Agrège les rapports de tous les outils QC
(FastQC, fastp, Picard, SAMtools flagstat, bcftools stats)
en un seul rapport HTML interactif. Indispensable pour
comparer la qualité de dizaines d'échantillons d'un coup.

---

## Alignement

### BWA-MEM2
**Rôle :** Aligne les reads Illumina paired-end sur le
génome de référence Amel_HAv3.1. Version optimisée de
BWA-MEM, 3 à 4x plus rapide grâce aux instructions SIMD.
Produit un fichier SAM converti immédiatement en BAM.
**Génome référence :** Amel_HAv3.1 (GCF_003254395.2)
16 chromosomes d'Apis mellifera.

### SAMtools (v1.19.2)
**Rôle :** Manipulation des fichiers BAM/SAM.
- `samtools sort` : trie le BAM par coordonnées génomiques
  (requis par GATK)
- `samtools index` : crée l'index .bai (accès rapide)
- `samtools flagstat` : statistiques d'alignement
  (% reads alignés, reads dupliqués, etc.)

---

## Traitement des BAM

### Picard MarkDuplicates
**Rôle :** Identifie et marque les duplicats PCR — reads
qui proviennent du même fragment d'ADN original amplifié
plusieurs fois. Les duplicats biaiseraient le comptage
des variants (sur-représentation artificielle d'allèles).
GATK les ignore automatiquement lors du variant calling.

---

## Variant calling

### GATK4 HaplotypeCaller
**Rôle :** Détecte les SNPs et indels dans chaque échantillon
individuellement. Mode GVCF : produit un fichier intermédiaire
contenant la confiance de génotypage à chaque position,
pas seulement aux variants détectés. Indispensable pour
le génotypage joint.
**Suit les GATK Best Practices** : standard de référence
en génomique.

### GATK4 GenomicsDBImport + GenotypeGVCFs
**Rôle :** Génotypage joint de tous les échantillons
simultanément. GenomicsDBImport consolide tous les GVCFs
dans une base de données génomique. GenotypeGVCFs produit
le VCF final multi-échantillons. Le génotypage joint est
plus puissant qu'un génotypage individuel car il exploite
l'information de tous les échantillons pour appeler les
variants rares.

### GATK4 VariantFiltration
**Rôle :** Filtre les variants de mauvaise qualité selon
les critères GATK Best Practices :
- QD (Quality by Depth) < 2.0
- FS (FisherStrand) > 60.0
- MQ (MappingQuality) < 40.0
- SOR (StrandOddsRatio) > 3.0

### bcftools (v1.19)
**Rôle :** Manipulation et filtrage des VCF.
- `bcftools filter` : garde uniquement les variants PASS
  bialléliques (un seul allèle alternatif)
- `bcftools stats` : statistiques complètes du VCF
  (nb SNPs, transitions/transversions, distribution des
  fréquences alléliques)

### SnpEff (v5.2)
**Rôle :** Annotation fonctionnelle des variants. Pour chaque
SNP, prédit l'effet biologique : synonyme, non-synonyme,
stop-gain, intronique, intergénique... Indispensable pour
identifier les variants candidats GWAS.
**Base de données :** construite localement par le process
`SNPEFF_BUILD` à partir du génome de référence et du GFF3 NCBI
(`snpEff build -gff3`). Aucun téléchargement à l'exécution, et
l'assemblage correspond exactement à celui du variant calling.
Les options `-noCheckCds -noCheckProtein` évitent d'avoir besoin
des fichiers CDS/protéines de référence (indisponibles hors-ligne).

---

## Génétique des populations

### PLINK2 (v2.00a5.12)
**Rôle :** Outil central de la génomique des populations
et du GWAS. Utilisé pour :
1. **QC génotypique** : filtre MAF < 5%, taux de génotypage
   manquant, déséquilibre Hardy-Weinberg
2. **Élagage LD** (LD pruning) : supprime les SNPs en fort
   déséquilibre de liaison pour obtenir un set indépendant
3. **PCA** : analyse en composantes principales pour
   visualiser la structure de population
4. **GWAS** : tests d'association linéaire/logistique

### ADMIXTURE (v1.3.0)
**Rôle :** Estime les proportions d'ascendance de chaque
individu dans K populations ancestrales. Teste plusieurs
valeurs de K et sélectionne la meilleure par
validation croisée (cross-validation error).
Produit les matrices Q (proportions d'ascendance)
et P (fréquences alléliques ancestrales).

### vcftools (v0.1.16)
**Rôle :** Calcul de statistiques de génétique des
populations directement depuis le VCF :
- **FST de Weir & Cockerham** : mesure de différenciation
  génétique entre populations (0=identiques, 1=fixés)
- **Déclin du LD** : r² en fonction de la distance physique
  — caractéristique de chaque population
- **Diversité nucléotidique (π)** : diversité intra-population

### GEMMA (v0.98.5)
**Rôle :** GWAS par modèle mixte linéaire (LMM). Corrige
la stratification de population en intégrant une matrice
de parenté génomique (kinship). Plus puissant que la
régression sur les PC car capture la parenté cryptique.
Calcule trois statistiques : Wald, rapport de vraisemblance,
score — robustesse maximale des résultats.

---

## Données

### SRA Toolkit (v3.4.1)
**Rôle :** Téléchargement des données WGS publiques depuis
NCBI SRA. `prefetch` télécharge les fichiers SRA,
`fasterq-dump` les convertit en FASTQ gzippés.
**Dataset utilisé :** PRJNA473480 (Wallberg et al. 2019)
WGS d'Apis mellifera mellifera populations européennes.

---

## Visualisation & Rapport

### R + ggplot2
**Rôle :** Génération des figures publication-ready :
Manhattan plot, QQ plot avec λ GC, PCA plot, barplot
ADMIXTURE, heatmap FST, courbe de déclin LD.

### R Markdown / knitr
**Rôle :** Rapport scientifique HTML/PDF automatique
intégrant toutes les figures et statistiques clés.
Reproductible : le rapport se regénère automatiquement
à chaque run du pipeline.

### qqman (R package)
**Rôle :** Package R spécialisé pour les Manhattan plots
et QQ plots de GWAS. Calcule automatiquement le facteur
d'inflation génomique λ (lambda GC) qui mesure la
présence de biais systématique dans les p-values GWAS.
Un λ proche de 1.0 indique une absence de stratification.
