# Lancer le pipeline sur une cohorte réelle

Ce document décrit comment passer de la **démonstration** (5 échantillons,
`-profile test`) à une **analyse GWAS réelle** (≥ 20–30 individus).

---

## 1. Préparer les données

```bash
# Génome + annotation + N échantillons WGS depuis PRJNA473480
bash bin/download_data.sh 30
```

Le script écrit un `samplesheet.csv` de base. **Ajouter une colonne
`phenotype`** (trait quantitatif ou binaire 0/1) pour activer le GWAS :

```csv
sample,fastq_1,fastq_2,sex,population,phenotype
AMM_FR_001,data/samples/SRR7190186_1.fastq.gz,data/samples/SRR7190186_2.fastq.gz,unknown,AMM,62.5
AMM_ES_004,data/samples/SRR7190194_1.fastq.gz,data/samples/SRR7190194_2.fastq.gz,unknown,AMM,48.1
...
```

- `population` sert au calcul du **FST** (comparaison AMM vs autres) et à la
  coloration des figures PCA / ADMIXTURE.
- `phenotype` vide ou `NA` → l'individu est exclu de l'association (mais garde
  sa place dans le variant calling et la génétique des populations).

---

## 2. Configurer

Copier et éditer le modèle de production :

```bash
cp conf/params.yml conf/my_run.yml
```

> `conf/demo10.yml` est un exemple intermédiaire prêt à l'emploi : 10 génomes
> A. m. mellifera, paramètres adaptés à un petit effectif (`pca_components: 8`,
> `plink_bad_ld: true`), sans phénotype.

Points d'attention pour une vraie cohorte :

| Paramètre | Démo (5 ind.) | Production | Remarque |
|---|---|---|---|
| `maf` | 0.01 | **0.05** | les SNPs rares n'ont pas de puissance en GWAS |
| `geno` | 0.5 | **0.05** | max 5 % de génotypes manquants par SNP |
| `mind` | 0.5 | **0.10** | exclut les individus mal génotypés |
| `pca_components` | 4 | **20** | doit rester `< nombre d'individus` |
| `plink_bad_ld` | true | **false** | `true` seulement si `< 50` individus |
| `admixture_k` | "2,3" | **"2,3,4,5"** | l'abeille : 3–4 lignées (A, M, C, O) |
| `hwe_filter` | false | false | garder `false` en population structurée |

---

## 3. Lancer

```bash
# Machine unique avec Docker
nextflow run main.nf -params-file conf/my_run.yml -profile docker -resume

# Cluster SLURM + Singularity
nextflow run main.nf -params-file conf/my_run.yml -profile slurm -resume
```

Adapter les plafonds de ressources dans `conf/my_run.yml` (`max_cpus`,
`max_memory`, `max_time`) : ils bornent automatiquement toutes les tâches
via `process.resourceLimits`.

---

## 4. Ordres de grandeur (30 génomes d'abeille, couverture ~10×)

| Ressource | Estimation |
|---|---|
| Téléchargement SRA | ~150–250 Go |
| Espace disque `work/` | ~400–600 Go (nettoyable avec `nextflow clean -f`) |
| Étape la plus longue | `GATK_HAPLOTYPECALLER` (~1–3 h/individu, parallélisable) |
| Durée totale (8 cœurs) | ~1–2 jours ; quelques heures sur un cluster |
| RAM pic | `BWA_MEM2_ALIGN` (~10 Go, index du génome) |

---

## 5. Résultats clés

| Fichier | Contenu |
|---|---|
| `results/03_variants/annotated/cohort.annotated.vcf.gz` | VCF + effets SnpEff (`ANN=`) |
| `results/04_population/pca/honeybee.eigenvec` | coordonnées PCA |
| `results/04_population/admixture/*.Q` | proportions d'ascendance |
| `results/05_gwas/gemma/output/honeybee.assoc.txt` | statistiques d'association (Wald / LRT / Score) |
| `results/05_gwas/plots/manhattan_plot.png` | Manhattan plot |
| `results/06_report/gwas_report.html` | rapport de synthèse |

Pour interpréter un pic GWAS : croiser la position du SNP avec le champ `ANN=`
du VCF annoté (`bcftools view -H cohort.annotated.vcf.gz <chr:pos>`).
