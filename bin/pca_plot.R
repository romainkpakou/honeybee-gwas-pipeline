#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# PCA plot — Structure de population
# Author : Romain KPAKOU
#
# INTERPRÉTATION :
# Chaque point = un individu
# Les individus de la même population se regroupent
# PC1 capture la plus grande source de variation (souvent la
# différence entre sous-espèces)
# PC2 capture la deuxième source de variation (souvent la
# structure géographique au sein des sous-espèces)
# Les ellipses représentent l'intervalle de confiance à 95%
# pour chaque population
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
})

opt_list <- list(
    make_option("--eigenvec", type="character", help="Fichier .eigenvec PLINK2"),
    make_option("--eigenval", type="character", help="Fichier .eigenval PLINK2"),
    make_option("--output",   type="character", default="pca_plot")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# ── Chargement ────────────────────────────────────────────────────────────────
df <- read.table(opt$eigenvec, header=FALSE, stringsAsFactors=FALSE)
colnames(df)[1:2] <- c("FID", "IID")
n_pc <- ncol(df) - 2
colnames(df)[3:ncol(df)] <- paste0("PC", 1:n_pc)

# Extraire la population depuis le nom de l'échantillon
# Format : POPULATION_PAYS_NUMERO (ex: AMM_FR_001)
df$Population <- sub("_[^_]+_[^_]+$", "", df$IID)
df$Country    <- sub("^[^_]+_([^_]+)_.*$", "\\1", df$IID)

# Variance expliquée par chaque PC
if (!is.null(opt$eigenval)) {
    eigenval <- scan(opt$eigenval, quiet=TRUE)
    pct_var  <- round(eigenval / sum(eigenval) * 100, 1)
    x_label  <- sprintf("PC1 (%.1f%% variance)", pct_var[1])
    y_label  <- sprintf("PC2 (%.1f%% variance)", pct_var[2])
} else {
    x_label <- "PC1"
    y_label <- "PC2"
}

cat(sprintf("Individus : %d\n", nrow(df)))
cat(sprintf("Populations : %s\n", paste(unique(df$Population), collapse=", ")))

# ── PC1 vs PC2 ────────────────────────────────────────────────────────────────
p1 <- ggplot(df, aes(x=PC1, y=PC2, colour=Population, shape=Country)) +
    geom_point(size=3, alpha=0.85) +
    stat_ellipse(aes(group=Population), level=0.95,
                 linetype="dashed", linewidth=0.6) +
    labs(
        title    = "Structure de population — PCA (PC1 vs PC2)",
        subtitle = "Apis mellifera mellifera | honeybee-gwas-pipeline",
        x        = x_label,
        y        = y_label,
        colour   = "Population",
        shape    = "Pays"
    ) +
    theme_classic(base_size=12) +
    theme(
        legend.position = "right",
        plot.title      = element_text(face="bold", size=14),
        plot.subtitle   = element_text(size=10, colour="grey40")
    )

# ── PC1 vs PC3 ────────────────────────────────────────────────────────────────
p2 <- ggplot(df, aes(x=PC1, y=PC3, colour=Population, shape=Country)) +
    geom_point(size=3, alpha=0.85) +
    stat_ellipse(aes(group=Population), level=0.95,
                 linetype="dashed", linewidth=0.6) +
    labs(
        title    = "Structure de population — PCA (PC1 vs PC3)",
        subtitle = "Apis mellifera mellifera | honeybee-gwas-pipeline",
        x        = x_label,
        y        = ifelse(!is.null(opt$eigenval),
                          sprintf("PC3 (%.1f%% variance)", pct_var[3]),
                          "PC3"),
        colour   = "Population",
        shape    = "Pays"
    ) +
    theme_classic(base_size=12) +
    theme(
        legend.position = "right",
        plot.title      = element_text(face="bold", size=14),
        plot.subtitle   = element_text(size=10, colour="grey40")
    )

# ── Sauvegarde ────────────────────────────────────────────────────────────────
ggsave(paste0(opt$output, "_PC1_PC2.png"), p1, width=10, height=8,
       dpi=300, bg="white")
ggsave(paste0(opt$output, "_PC1_PC3.png"), p2, width=10, height=8,
       dpi=300, bg="white")
ggsave(paste0(opt$output, "_PC1_PC2.pdf"), p1, width=10, height=8)
cat("PCA plots sauvegardés\n")
