#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Admixture barplot — Proportions d'ascendance
# Author : Romain KPAKOU
#
# INTERPRÉTATION :
# Chaque barre = un individu
# Chaque couleur = une population ancestrale
# Une barre d'une seule couleur = individu "pur" pour cette population
# Une barre mélangée = individu hybride entre populations
# Pour l'abeille noire : on cherche des barres uniformes dans la
# population AMM = bonne pureté génétique
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(tidyr)
    library(dplyr)
})

opt_list <- list(
    make_option("--qdir",   type="character", default=".",
                help="Répertoire contenant les fichiers .Q"),
    make_option("--fam",    type="character", help="Fichier .fam PLINK"),
    make_option("--output", type="character", default="admixture_plot")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# ── Charger les fichiers Q pour chaque K ─────────────────────────────────────
q_files <- list.files(opt$qdir, pattern="\\.Q$", full.names=TRUE)
q_files <- sort(q_files)

if (length(q_files) == 0) {
    cat("Aucun fichier .Q trouvé dans", opt$qdir, "\n")
    quit(status=0)
}

# ── Charger le FAM pour les identifiants ─────────────────────────────────────
if (!is.null(opt$fam) && file.exists(opt$fam)) {
    fam <- read.table(opt$fam, header=FALSE)
    sample_ids <- fam$V2
} else {
    sample_ids <- NULL
}

# ── Palette de couleurs ───────────────────────────────────────────────────────
palette_colors <- c("#E74C3C","#3498DB","#2ECC71","#F39C12",
                    "#9B59B6","#1ABC9C","#E67E22","#34495E")

# ── Créer un barplot pour chaque K ───────────────────────────────────────────
plots <- list()

for (f in q_files) {
    k <- as.integer(gsub(".*\\.(\\d+)\\.Q$", "\\1", basename(f)))
    q <- read.table(f, header=FALSE)
    colnames(q) <- paste0("Anc", 1:ncol(q))

    if (!is.null(sample_ids) && nrow(q) == length(sample_ids)) {
        q$Sample <- sample_ids
        q$Population <- sub("_[^_]+_[^_]+$", "", sample_ids)
    } else {
        q$Sample <- paste0("S", 1:nrow(q))
        q$Population <- "Unknown"
    }

    # Trier les individus par population puis par proportion de Anc1
    q <- q %>%
        arrange(Population, desc(Anc1)) %>%
        mutate(Order = row_number())

    q_long <- pivot_longer(q,
        cols      = starts_with("Anc"),
        names_to  = "Ancestry",
        values_to = "Proportion"
    )

    p <- ggplot(q_long, aes(x=Order, y=Proportion, fill=Ancestry)) +
        geom_bar(stat="identity", width=1) +
        scale_fill_manual(values=palette_colors[1:k]) +
        scale_y_continuous(expand=c(0, 0)) +
        scale_x_continuous(expand=c(0, 0)) +
        labs(
            title = sprintf("K = %d", k),
            x     = "Individus",
            y     = "Proportion d'ascendance"
        ) +
        theme_classic(base_size=10) +
        theme(
            axis.text.x     = element_blank(),
            axis.ticks.x    = element_blank(),
            legend.position = "right",
            plot.title      = element_text(face="bold", size=12)
        )

    plots[[length(plots) + 1]] <- p
}

# ── Sauvegarder tous les K dans une seule figure ──────────────────────────────
if (requireNamespace("gridExtra", quietly=TRUE)) {
    g <- gridExtra::arrangeGrob(grobs=plots, ncol=1)
    ggsave(paste0(opt$output, ".png"), g,
           width=14, height=3*length(plots), dpi=300, bg="white")
    ggsave(paste0(opt$output, ".pdf"), g,
           width=14, height=3*length(plots))
} else {
    ggsave(paste0(opt$output, ".png"), plots[[1]],
           width=14, height=4, dpi=300, bg="white")
}
cat("Admixture plot sauvegardé\n")
