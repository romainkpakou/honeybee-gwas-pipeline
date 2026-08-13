#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# LD decay plot — Déclin du déséquilibre de liaison
# Author : Romain KPAKOU
#
# INTERPRÉTATION :
# Axe X : distance physique entre paires de SNPs (kb)
# Axe Y : corrélation moyenne (r²) entre paires de SNPs
# La courbe montre à quelle vitesse le LD décline avec la distance.
# Chez Apis mellifera, le LD décline très rapidement (~10-50 kb)
# car le taux de recombinaison est très élevé (~19 cM/Mb).
# La ligne pointillée rouge indique r² = 0.2 (seuil commun)
# La distance à laquelle la courbe croise ce seuil définit la
# "portée du LD" — information clé pour interpréter les pics GWAS.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
})

opt_list <- list(
    make_option("--input",  type="character", help="Fichier .geno.ld vcftools"),
    make_option("--output", type="character", default="ld_decay"),
    make_option("--max_kb", type="integer",   default=500,
                help="Distance maximale à afficher (kb)")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# ── Chargement ────────────────────────────────────────────────────────────────
cat("Chargement des données LD...\n")
df <- read.table(opt$input, header=TRUE, stringsAsFactors=FALSE)
colnames(df) <- c("CHR", "POS1", "POS2", "N_INDV", "R2")

# Distance en kb
df$DIST_KB <- abs(df$POS2 - df$POS1) / 1000
df <- df[df$DIST_KB <= opt$max_kb & !is.na(df$R2), ]

cat(sprintf("Paires de SNPs analysées : %s\n", format(nrow(df), big.mark=",")))

# ── Binning par distance ──────────────────────────────────────────────────────
# On divise la distance en intervalles de 5 kb et on calcule
# le r² moyen dans chaque intervalle
df$BIN_KB <- cut(df$DIST_KB,
                 breaks=seq(0, opt$max_kb, 5),
                 include.lowest=TRUE)

ld_mean <- df %>%
    group_by(BIN_KB) %>%
    summarise(
        mean_r2 = mean(R2, na.rm=TRUE),
        mid_kb  = mean(DIST_KB, na.rm=TRUE),
        n       = n(),
        .groups = "drop"
    ) %>%
    filter(!is.na(mid_kb))

# Distance où r² croise 0.2
ld_threshold <- ld_mean %>%
    filter(mean_r2 <= 0.2) %>%
    slice_min(mid_kb, n=1)

if (nrow(ld_threshold) > 0) {
    cat(sprintf("Portée du LD (r² < 0.2) : ~%.0f kb\n",
                ld_threshold$mid_kb))
}

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(ld_mean, aes(x=mid_kb, y=mean_r2)) +
    geom_line(colour="#2C3E50", linewidth=1.2) +
    geom_hline(yintercept=0.2,
               linetype="dashed", colour="red", linewidth=0.7) +
    annotate("text", x=opt$max_kb*0.95, y=0.22,
             label="r² = 0.2", colour="red", size=3.5, hjust=1) +
    {if (nrow(ld_threshold) > 0)
        geom_vline(xintercept=ld_threshold$mid_kb,
                   linetype="dotted", colour="blue", linewidth=0.6)
    } +
    labs(
        title    = "Déclin du LD — Apis mellifera mellifera",
        subtitle = sprintf(
            "Portée du LD (r² < 0.2) : ~%.0f kb | Taux recombinaison élevé (~19 cM/Mb)",
            ifelse(nrow(ld_threshold) > 0, ld_threshold$mid_kb, NA)
        ),
        x = "Distance physique (kb)",
        y = expression(Moyenne~r^2)
    ) +
    scale_x_continuous(expand=c(0, 5)) +
    scale_y_continuous(limits=c(0, 1), expand=c(0, 0)) +
    theme_classic(base_size=12) +
    theme(
        plot.title    = element_text(face="bold", size=14),
        plot.subtitle = element_text(size=10, colour="grey40")
    )

ggsave(paste0(opt$output, ".png"), p, width=10, height=6, dpi=300, bg="white")
ggsave(paste0(opt$output, ".pdf"), p, width=10, height=6)
cat("LD decay plot sauvegardé\n")
