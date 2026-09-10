# Exemple de sorties

Rapports produits par le **profil de démonstration** :

```bash
nextflow run main.nf -profile test,docker
```

> ⚠️ **Démonstration technique, pas une analyse biologique.**
> Jeu de 5 génomes *Apis mellifera mellifera* (PRJNA473480, ~10×), filtres QC
> volontairement relâchés. Les statistiques de génétique des populations et le
> GWAS n'ont donc **pas** de valeur biologique à cette échelle. Ce dossier
> illustre uniquement que le pipeline s'exécute et produit les livrables
> attendus. Pour une analyse réelle : [`../running_a_real_cohort.md`](../running_a_real_cohort.md).

Seuls les rapports **réellement informatifs à n = 5** sont inclus ici (le
contrôle qualité et l'annotation fonctionnelle portent sur de vraies données) :

| Fichier | Contenu |
|---|---|
| `multiqc_report.html` | MultiQC : Samtools flagstat, Picard MarkDuplicates, bcftools stats, résumé SnpEff — 5 génomes |
| `snpeff_summary.html` | Rapport SnpEff : distribution des effets, impacts, Ts/Tv |
| `snpeff_summary.csv` | Mêmes chiffres, format tableau ([rendu directement par GitHub](snpeff_summary.csv)) |

**Visualiser les HTML** (GitHub affiche le code source, pas la page rendue) :
soit via [htmlpreview.github.io](https://htmlpreview.github.io/?https://github.com/romainkpakou/honeybee-gwas-pipeline/blob/main/docs/example_output/snpeff_summary.html),
soit en activant **GitHub Pages** (Settings → Pages → *Deploy from branch* `main` `/docs`)
puis `https://romainkpakou.github.io/honeybee-gwas-pipeline/example_output/multiqc_report.html`.

## Chiffres clés de l'annotation SnpEff

| Métrique | Valeur |
|---|---|
| Variants annotés | 227 / 227 |
| Effets prédits | 2 264 |
| `missense_variant` | 88 |
| `synonymous_variant` | 13 |
| `stop_gained` | 5 |
| Ratio missense / silent | 5,87 |

Base SnpEff construite localement depuis le GFF3 NCBI (Release 104,
`GCF_003254395.2`) — assemblage identique à celui du variant calling.
