#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Manhattan plot — publication-ready
# Author  : Romain KPAKOU
# Input   : fichier de résultats GEMMA (.assoc.txt)
# Output  : manhattan_plot.png + manhattan_plot.pdf
#
# STRUCTURE DU MANHATTAN PLOT :
# Axe X : position génomique (chromosomes 1-16 de l'abeille)
# Axe Y : -log10(p-value) du test de Wald (GEMMA)
# Ligne rouge    : seuil Bonferroni (p < 0.05 / N_SNPs)
# Ligne orange   : seuil suggestif (p < 1e-5)
# Points au-dessus de la ligne rouge = SNPs significatifs
#
# INTERPRÉTATION :
# Un pic net au-dessus du seuil Bonferroni indique une association
# forte entre cette région génomique et le phénotype testé.
# Chez l'abeille (fort taux de recombinaison), les pics sont
# généralement étroits — bonne résolution de cartographie.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
opt_list <- list(
    make_option("--input",      type="character", help="Fichier résultats GEMMA"),
    make_option("--output",     type="character", default="manhattan_plot"),
    make_option("--threshold",  type="double",    default=5e-8),
    make_option("--suggestive", type="double",    default=1e-5),
    make_option("--title",      type="character",
                default="GWAS Manhattan Plot — Apis mellifera mellifera")
)
opt <- parse_args(OptionParser(option_list=opt_list))

cat("Chargement des résultats GWAS :", opt$input, "\n")

# ── Chargement et préparation des données ────────────────────────────────────
df <- read.table(opt$input, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Standardiser les noms de colonnes (GEMMA vs PLINK2)
if ("p_wald"  %in% names(df)) df$P   <- df$p_wald
if ("ps"      %in% names(df)) df$BP  <- df$ps
if ("chr"     %in% names(df)) df$CHR <- as.integer(df$chr)

# Nettoyer les données
df <- df[!is.na(df$P) & df$P > 0 & !is.na(df$CHR), ]
df$LOG10P <- -log10(df$P)

n_snps <- nrow(df)
bonferroni <- 0.05 / n_snps
cat(sprintf("SNPs analysés : %d\n", n_snps))
cat(sprintf("Seuil Bonferroni : p < %.2e\n", bonferroni))
cat(sprintf("SNPs significatifs : %d\n", sum(df$P < bonferroni)))

# ── Calcul des positions cumulatives ─────────────────────────────────────────
# Nécessaire pour placer les chromosomes bout à bout sur l'axe X
chr_lengths <- df %>%
    group_by(CHR) %>%
    summarise(max_bp = max(BP), .groups='drop') %>%
    arrange(CHR) %>%
    mutate(bp_offset = cumsum(lag(max_bp, default=0)))

df <- df %>%
    left_join(chr_lengths, by="CHR") %>%
    mutate(BP_cum = BP + bp_offset)

# Position centrale de chaque chromosome pour les labels de l'axe X
axis_df <- df %>%
    group_by(CHR) %>%
    summarise(center = (max(BP_cum) + min(BP_cum)) / 2, .groups='drop')

# ── Palette de couleurs alternée par chromosome ───────────────────────────────
chr_colours <- rep(c("#2C3E50", "#E67E22"), length(unique(df$CHR)))

# ── Construction du Manhattan plot ───────────────────────────────────────────
p <- ggplot(df, aes(x=BP_cum, y=LOG10P, colour=as.factor(CHR))) +
    # Points SNPs
    geom_point(size=0.6, alpha=0.7) +
    # Couleurs alternées par chromosome
    scale_colour_manual(values=chr_colours) +
    # Axe X : numéros de chromosomes centrés
    scale_x_continuous(
        labels = axis_df$CHR,
        breaks = axis_df$center,
        expand = c(0.01, 0)
    ) +
    scale_y_continuous(expand=c(0, 0.3)) +
    # Seuil Bonferroni (rouge)
    geom_hline(
        yintercept = -log10(bonferroni),
        colour = "red", linetype = "dashed", linewidth = 0.7
    ) +
    # Seuil suggestif (orange)
    geom_hline(
        yintercept = -log10(opt$suggestive),
        colour = "orange", linetype = "dashed", linewidth = 0.5
    ) +
    # Annotations des seuils
    annotate("text",
        x = max(df$BP_cum) * 0.98,
        y = -log10(bonferroni) + 0.4,
        label = sprintf("Bonferroni (p=%.0e)", bonferroni),
        colour = "red", size = 3, hjust = 1
    ) +
    annotate("text",
        x = max(df$BP_cum) * 0.98,
        y = -log10(opt$suggestive) + 0.4,
        label = sprintf("Suggestif (p=%.0e)", opt$suggestive),
        colour = "orange", size = 3, hjust = 1
    ) +
    # Labels
    labs(
        title    = opt$title,
        subtitle = sprintf(
            "n = %s SNPs | Bonferroni p < %.0e | %d SNPs significatifs",
            format(n_snps, big.mark=","),
            bonferroni,
            sum(df$P < bonferroni)
        ),
        x = "Chromosome",
        y = expression(-log[10](p))
    ) +
    # Thème publication
    theme_classic(base_size=12) +
    theme(
        legend.position   = "none",
        panel.grid        = element_blank(),
        axis.text.x       = element_text(size=9),
        axis.text.y       = element_text(size=10),
        plot.title        = element_text(face="bold", size=14),
        plot.subtitle     = element_text(size=10, colour="grey40"),
        plot.margin       = margin(10, 20, 10, 10)
    )

# ── Sauvegarde ────────────────────────────────────────────────────────────────
ggsave(paste0(opt$output, ".png"), p, width=14, height=6, dpi=300, bg="white")
ggsave(paste0(opt$output, ".pdf"), p, width=14, height=6)
cat("Figures sauvegardées :", paste0(opt$output, ".png / .pdf"), "\n")
