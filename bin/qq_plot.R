#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# QQ plot avec facteur d'inflation génomique (lambda GC)
# Author : Romain KPAKOU
#
# STRUCTURE DU QQ PLOT :
# Axe X : p-values attendues sous H0 (-log10 de la distribution uniforme)
# Axe Y : p-values observées (-log10)
# Diagonale rouge : la droite y=x (attendu si aucune association)
# Zone grise : intervalle de confiance à 95%
#
# INTERPRÉTATION :
# Si les points suivent la diagonale jusqu'à la queue de distribution
# puis s'en écartent brusquement vers le haut → signal GWAS réel
# Si les points s'écartent dès le début → inflation = problème de
# stratification non corrigé → revoir la correction par kinship
#
# LAMBDA GC (facteur d'inflation génomique) :
# λ = médiane(χ² observé) / médiane(χ² attendu sous H0)
# λ = 1.00 : aucune inflation (parfait)
# λ = 1.05 : légère inflation acceptable
# λ > 1.10 : inflation à investiguer
# λ >> 1.0 : stratification non corrigée (problème sérieux)
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
})

opt_list <- list(
    make_option("--input",  type="character", help="Fichier résultats GEMMA"),
    make_option("--output", type="character", default="qq_plot"),
    make_option("--title",  type="character",
                default="QQ Plot — Apis mellifera mellifera GWAS")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# ── Chargement ────────────────────────────────────────────────────────────────
df <- read.table(opt$input, header=TRUE, sep="\t", stringsAsFactors=FALSE)
if ("p_wald" %in% names(df)) df$P <- df$p_wald
df <- df[!is.na(df$P) & df$P > 0, ]
n  <- nrow(df)

# ── Calcul du lambda GC ───────────────────────────────────────────────────────
chisq  <- qchisq(1 - df$P, df=1)
lambda <- round(median(chisq, na.rm=TRUE) / qchisq(0.5, df=1), 4)
cat(sprintf("Lambda GC = %.4f\n", lambda))
cat(sprintf("Interprétation : %s\n",
    ifelse(lambda < 1.05, "Pas d'inflation détectée",
    ifelse(lambda < 1.10, "Légère inflation acceptable",
    "Inflation significative — vérifier la correction"))))

# ── Données QQ ────────────────────────────────────────────────────────────────
observed <- sort(-log10(df$P))
expected <- -log10(ppoints(n))

# Intervalle de confiance à 95%
ci_lo <- -log10(qbeta(0.975, 1:n, n:1))
ci_hi <- -log10(qbeta(0.025, 1:n, n:1))

qq_df <- data.frame(
    expected = expected,
    observed = observed,
    ci_lo    = sort(ci_lo),
    ci_hi    = sort(ci_hi)
)

# ── Construction du QQ plot ───────────────────────────────────────────────────
p <- ggplot(qq_df, aes(x=expected, y=observed)) +
    # Intervalle de confiance 95%
    geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi),
                fill="grey80", alpha=0.7) +
    # Points observés
    geom_point(size=0.8, alpha=0.8, colour="#2C3E50") +
    # Diagonale (attendu sous H0)
    geom_abline(intercept=0, slope=1,
                colour="red", linewidth=0.8) +
    # Annotation lambda
    annotate("text",
        x = max(expected) * 0.05,
        y = max(observed) * 0.95,
        label = sprintf("λ = %.4f", lambda),
        size = 5, hjust=0, fontface="bold",
        colour = ifelse(lambda > 1.1, "red", "black")
    ) +
    labs(
        title    = opt$title,
        subtitle = sprintf(
            "n = %s SNPs | λ (inflation génomique) = %.4f",
            format(n, big.mark=","), lambda
        ),
        x = expression(Expected~~-log[10](p)),
        y = expression(Observed~~-log[10](p))
    ) +
    theme_classic(base_size=12) +
    theme(
        plot.title    = element_text(face="bold", size=14),
        plot.subtitle = element_text(size=10, colour="grey40")
    )

# ── Sauvegarde ────────────────────────────────────────────────────────────────
ggsave(paste0(opt$output, ".png"), p, width=8, height=8, dpi=300, bg="white")
ggsave(paste0(opt$output, ".pdf"), p, width=8, height=8)
cat("QQ plot sauvegardé\n")
