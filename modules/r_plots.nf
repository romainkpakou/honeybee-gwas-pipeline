/*
    MODULE : r_plots
    Outil   : R + ggplot2
    Rôle    : Génération des figures publication-ready
    Docker  : quay.io/biocontainers/r-base:4.3.3

    IMPORTANT DSL2 :
    Le code R inline dans les scripts Nextflow doit être protégé
    avec des guillemets simples (heredoc << 'REOF') pour éviter
    que Nextflow n'interprète les variables R ($var) comme des
    variables Groovy/Nextflow.

    6 PROCESS :
    1. PLOT_MANHATTAN  — Manhattan plot GWAS
    2. PLOT_QQ         — QQ plot avec lambda GC
    3. PLOT_PCA        — Structure de population (PC1 vs PC2/PC3)
    4. PLOT_ADMIXTURE  — Barplot proportions d'ascendance
    5. PLOT_LD_DECAY   — Courbe de déclin du LD
    6. PLOT_FST        — Différenciation génétique FST
*/

// ── Process 1 : Manhattan plot ────────────────────────────────────────────────
process PLOT_MANHATTAN {
    tag "manhattan"
    label 'process_low'

    publishDir "${params.outdir}/05_gwas/plots", mode: 'copy'

    container 'quay.io/biocontainers/r-base:4.3.3'

    input:
    path gwas_results

    output:
    path "manhattan_plot.png", emit: png
    path "manhattan_plot.pdf", emit: pdf
    path "versions.yml",       emit: versions

    script:
    // On passe les paramètres Nextflow comme variables shell
    // puis R les lit via Sys.getenv() — évite les conflits $ Groovy vs R
    def threshold = params.gwas_pval
    """
    export GWAS_FILE="${gwas_results}"
    export THRESHOLD="${threshold}"

    Rscript - << 'REOF'
    suppressPackageStartupMessages({
        library(ggplot2)
        library(dplyr)
    })

    # Lire les paramètres depuis les variables d'environnement
    gwas_file <- Sys.getenv("GWAS_FILE")
    threshold <- as.numeric(Sys.getenv("THRESHOLD", "1e-6"))

    cat("Chargement :", gwas_file, "\\n")
    df <- read.table(gwas_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE)

    # Standardiser les noms de colonnes GEMMA → standard
    if ("p_wald" %in% names(df)) df\$P   <- df\$p_wald
    if ("ps"     %in% names(df)) df\$BP  <- df\$ps
    if ("chr"    %in% names(df)) df\$CHR <- as.integer(df\$chr)

    df <- df[!is.na(df\$P) & df\$P > 0 & !is.na(df\$CHR), ]
    df\$LOG10P <- -log10(df\$P)

    n_snps     <- nrow(df)
    bonferroni <- 0.05 / n_snps
    cat(sprintf("SNPs : %d | Bonferroni : %.2e | Sig : %d\\n",
        n_snps, bonferroni, sum(df\$P < bonferroni)))

    # Positions cumulatives pour l'axe X
    chr_lengths <- df %>%
        group_by(CHR) %>%
        summarise(max_bp = max(BP), .groups='drop') %>%
        arrange(CHR) %>%
        mutate(bp_offset = cumsum(lag(max_bp, default=0)))

    df <- df %>%
        left_join(chr_lengths, by="CHR") %>%
        mutate(BP_cum = BP + bp_offset)

    axis_df <- df %>%
        group_by(CHR) %>%
        summarise(center = (max(BP_cum) + min(BP_cum)) / 2, .groups='drop')

    chr_colours <- rep(c("#2C3E50","#E67E22"), length(unique(df\$CHR)))

    p <- ggplot(df, aes(x=BP_cum, y=LOG10P, colour=as.factor(CHR))) +
        geom_point(size=0.6, alpha=0.7) +
        scale_colour_manual(values=chr_colours) +
        scale_x_continuous(labels=axis_df\$CHR, breaks=axis_df\$center,
                           expand=c(0.01,0)) +
        scale_y_continuous(expand=c(0, 0.3)) +
        geom_hline(yintercept=-log10(bonferroni),
                   colour="red", linetype="dashed", linewidth=0.7) +
        geom_hline(yintercept=-log10(1e-5),
                   colour="orange", linetype="dashed", linewidth=0.5) +
        labs(
            title    = "GWAS Manhattan Plot — Apis mellifera mellifera",
            subtitle = sprintf("n = %s SNPs | Bonferroni p < %.0e | %d sig.",
                format(n_snps, big.mark=","), bonferroni,
                sum(df\$P < bonferroni)),
            x = "Chromosome",
            y = expression(-log[10](p))
        ) +
        theme_classic(base_size=12) +
        theme(legend.position="none",
              plot.title=element_text(face="bold", size=14),
              plot.subtitle=element_text(size=10, colour="grey40"))

    ggsave("manhattan_plot.png", p, width=14, height=6, dpi=300, bg="white")
    ggsave("manhattan_plot.pdf", p, width=14, height=6)
    cat("Manhattan plot sauvegardé\\n")
REOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch manhattan_plot.png manhattan_plot.pdf versions.yml
    """
}

// ── Process 2 : QQ plot ───────────────────────────────────────────────────────
process PLOT_QQ {
    tag "qq_plot"
    label 'process_low'

    publishDir "${params.outdir}/05_gwas/plots", mode: 'copy'

    container 'quay.io/biocontainers/r-base:4.3.3'

    input:
    path gwas_results

    output:
    path "qq_plot.png",  emit: png
    path "qq_plot.pdf",  emit: pdf
    path "versions.yml", emit: versions

    script:
    """
    export GWAS_FILE="${gwas_results}"

    Rscript - << 'REOF'
    suppressPackageStartupMessages({ library(ggplot2) })

    gwas_file <- Sys.getenv("GWAS_FILE")
    df <- read.table(gwas_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE)
    if ("p_wald" %in% names(df)) df\$P <- df\$p_wald
    df <- df[!is.na(df\$P) & df\$P > 0, ]
    n  <- nrow(df)

    # Facteur d'inflation génomique lambda GC
    chisq  <- qchisq(1 - df\$P, df=1)
    lambda <- round(median(chisq, na.rm=TRUE) / qchisq(0.5, df=1), 4)
    cat(sprintf("Lambda GC = %.4f\\n", lambda))

    observed <- sort(-log10(df\$P))
    expected <- -log10(ppoints(n))
    ci_lo    <- -log10(qbeta(0.975, 1:n, n:1))
    ci_hi    <- -log10(qbeta(0.025, 1:n, n:1))

    qq_df <- data.frame(
        expected = expected,
        observed = observed,
        ci_lo    = sort(ci_lo),
        ci_hi    = sort(ci_hi)
    )

    p <- ggplot(qq_df, aes(x=expected, y=observed)) +
        geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi), fill="grey80", alpha=0.7) +
        geom_point(size=0.8, alpha=0.8, colour="#2C3E50") +
        geom_abline(intercept=0, slope=1, colour="red", linewidth=0.8) +
        annotate("text", x=min(expected)+0.5, y=max(observed)*0.95,
                 label=sprintf("lambda = %.4f", lambda),
                 size=5, hjust=0, fontface="bold",
                 colour=ifelse(lambda > 1.1, "red", "black")) +
        labs(
            title    = "QQ Plot — Apis mellifera mellifera GWAS",
            subtitle = sprintf("n = %s SNPs | lambda GC = %.4f",
                format(n, big.mark=","), lambda),
            x = expression(Expected~~-log[10](p)),
            y = expression(Observed~~-log[10](p))
        ) +
        theme_classic(base_size=12) +
        theme(plot.title=element_text(face="bold", size=14),
              plot.subtitle=element_text(size=10, colour="grey40"))

    ggsave("qq_plot.png", p, width=8, height=8, dpi=300, bg="white")
    ggsave("qq_plot.pdf", p, width=8, height=8)
    cat("QQ plot sauvegardé\\n")
REOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch qq_plot.png qq_plot.pdf versions.yml
    """
}

// ── Process 3 : PCA plot ──────────────────────────────────────────────────────
process PLOT_PCA {
    tag "pca_plot"
    label 'process_low'

    publishDir "${params.outdir}/04_population/plots", mode: 'copy'

    container 'quay.io/biocontainers/r-base:4.3.3'

    input:
    path eigenvec
    path eigenval

    output:
    path "pca_plot_PC1_PC2.png", emit: pc12_png
    path "pca_plot_PC1_PC3.png", emit: pc13_png
    path "pca_plot_PC1_PC2.pdf", emit: pc12_pdf
    path "versions.yml",         emit: versions

    script:
    """
    export EIGENVEC_FILE="${eigenvec}"
    export EIGENVAL_FILE="${eigenval}"

    Rscript - << 'REOF'
    suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

    eigenvec_file <- Sys.getenv("EIGENVEC_FILE")
    eigenval_file <- Sys.getenv("EIGENVAL_FILE")

    df <- read.table(eigenvec_file, header=FALSE, stringsAsFactors=FALSE)
    colnames(df)[1:2] <- c("FID","IID")
    n_pc <- ncol(df) - 2
    colnames(df)[3:ncol(df)] <- paste0("PC", 1:n_pc)
    df\$Population <- sub("_[^_]+_[^_]+\$", "", df\$IID)

    eigenval <- scan(eigenval_file, quiet=TRUE)
    pct_var  <- round(eigenval / sum(eigenval) * 100, 1)

    make_pca_plot <- function(pc_x, pc_y, pct_x, pct_y) {
        ggplot(df, aes_string(x=pc_x, y=pc_y, colour="Population")) +
            geom_point(size=3, alpha=0.85) +
            stat_ellipse(level=0.95, linetype="dashed", linewidth=0.6) +
            labs(
                title    = sprintf("PCA — %s vs %s", pc_x, pc_y),
                subtitle = "Apis mellifera mellifera | honeybee-gwas-pipeline",
                x        = sprintf("%s (%.1f%% variance)", pc_x, pct_x),
                y        = sprintf("%s (%.1f%% variance)", pc_y, pct_y)
            ) +
            theme_classic(base_size=12) +
            theme(plot.title=element_text(face="bold"),
                  legend.position="right")
    }

    p1 <- make_pca_plot("PC1","PC2", pct_var[1], pct_var[2])
    p2 <- make_pca_plot("PC1","PC3", pct_var[1], pct_var[3])

    ggsave("pca_plot_PC1_PC2.png", p1, width=10, height=8, dpi=300, bg="white")
    ggsave("pca_plot_PC1_PC3.png", p2, width=10, height=8, dpi=300, bg="white")
    ggsave("pca_plot_PC1_PC2.pdf", p1, width=10, height=8)
    cat("PCA plots sauvegardés\\n")
REOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch pca_plot_PC1_PC2.png pca_plot_PC1_PC3.png pca_plot_PC1_PC2.pdf
    touch versions.yml
    """
}

// ── Process 4 : Admixture barplot ─────────────────────────────────────────────
process PLOT_ADMIXTURE {
    tag "admixture_plot"
    label 'process_low'

    publishDir "${params.outdir}/04_population/plots", mode: 'copy'

    container 'quay.io/biocontainers/r-base:4.3.3'

    input:
    path q_files
    path fam

    output:
    path "admixture_plot.png", emit: png
    path "admixture_plot.pdf", emit: pdf
    path "versions.yml",       emit: versions

    script:
    """
    export FAM_FILE="${fam}"

    Rscript - << 'REOF'
    suppressPackageStartupMessages({
        library(ggplot2); library(dplyr); library(tidyr)
    })

    fam_file <- Sys.getenv("FAM_FILE")
    q_files  <- list.files(".", pattern="\\\\.Q\$", full.names=TRUE)
    q_files  <- sort(q_files)

    if (length(q_files) == 0) {
        cat("Aucun fichier .Q trouvé\\n")
        file.create("admixture_plot.png")
        file.create("admixture_plot.pdf")
        quit(status=0)
    }

    # Charger les identifiants depuis le FAM
    if (file.exists(fam_file)) {
        fam <- read.table(fam_file, header=FALSE)
        sample_ids <- fam\$V2
    } else {
        sample_ids <- NULL
    }

    palette_colors <- c("#E74C3C","#3498DB","#2ECC71","#F39C12",
                        "#9B59B6","#1ABC9C","#E67E22","#34495E")

    plots <- lapply(q_files, function(f) {
        k <- as.integer(gsub(".*\\\\.(\\\\d+)\\\\.Q\$", "\\\\1", basename(f)))
        q <- read.table(f, header=FALSE)
        colnames(q) <- paste0("Anc", 1:ncol(q))

        if (!is.null(sample_ids) && nrow(q) == length(sample_ids)) {
            q\$Sample <- sample_ids
        } else {
            q\$Sample <- paste0("S", 1:nrow(q))
        }

        q <- q %>% arrange(desc(Anc1)) %>% mutate(Order = row_number())
        q_long <- pivot_longer(q, starts_with("Anc"),
                               names_to="Ancestry", values_to="Proportion")

        ggplot(q_long, aes(x=Order, y=Proportion, fill=Ancestry)) +
            geom_bar(stat="identity", width=1) +
            scale_fill_manual(values=palette_colors[1:k]) +
            scale_y_continuous(expand=c(0,0)) +
            scale_x_continuous(expand=c(0,0)) +
            labs(title=sprintf("K = %d", k), x="Individus",
                 y="Proportion d'ascendance") +
            theme_classic(base_size=10) +
            theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  plot.title=element_text(face="bold"))
    })

    if (requireNamespace("gridExtra", quietly=TRUE)) {
        g <- gridExtra::arrangeGrob(grobs=plots, ncol=1)
        ggsave("admixture_plot.png", g, width=14,
               height=3*length(plots), dpi=300, bg="white")
        ggsave("admixture_plot.pdf", g, width=14, height=3*length(plots))
    } else {
        ggsave("admixture_plot.png", plots[[1]], width=14, height=4,
               dpi=300, bg="white")
        ggsave("admixture_plot.pdf", plots[[1]], width=14, height=4)
    }
    cat("Admixture plot sauvegardé\\n")
REOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch admixture_plot.png admixture_plot.pdf versions.yml
    """
}

// ── Process 5 : LD decay plot ─────────────────────────────────────────────────
process PLOT_LD_DECAY {
    tag "ld_decay"
    label 'process_low'

    publishDir "${params.outdir}/04_population/plots", mode: 'copy'

    container 'quay.io/biocontainers/r-base:4.3.3'

    input:
    path ld_file

    output:
    path "ld_decay.png",  emit: png
    path "ld_decay.pdf",  emit: pdf
    path "versions.yml",  emit: versions

    script:
    """
    export LD_FILE="${ld_file}"

    Rscript - << 'REOF'
    suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

    ld_file <- Sys.getenv("LD_FILE")
    cat("Chargement LD :", ld_file, "\\n")

    df <- read.table(ld_file, header=TRUE, stringsAsFactors=FALSE)
    colnames(df) <- c("CHR","POS1","POS2","N_INDV","R2")
    df\$DIST_KB <- abs(df\$POS2 - df\$POS1) / 1000
    df <- df[df\$DIST_KB <= 500 & !is.na(df\$R2), ]

    df\$BIN_KB <- cut(df\$DIST_KB, breaks=seq(0,500,5), include.lowest=TRUE)
    ld_mean <- df %>%
        group_by(BIN_KB) %>%
        summarise(mean_r2=mean(R2, na.rm=TRUE),
                  mid_kb=mean(DIST_KB, na.rm=TRUE), .groups='drop') %>%
        filter(!is.na(mid_kb))

    ld_thresh <- ld_mean %>% filter(mean_r2 <= 0.2) %>%
        slice_min(mid_kb, n=1)
    if (nrow(ld_thresh) > 0)
        cat(sprintf("Portee LD (r2 < 0.2) : ~%.0f kb\\n", ld_thresh\$mid_kb))

    p <- ggplot(ld_mean, aes(x=mid_kb, y=mean_r2)) +
        geom_line(colour="#2C3E50", linewidth=1.2) +
        geom_hline(yintercept=0.2, linetype="dashed",
                   colour="red", linewidth=0.7) +
        annotate("text", x=480, y=0.22, label="r2 = 0.2",
                 colour="red", size=3.5, hjust=1) +
        labs(
            title    = "Declin du LD — Apis mellifera mellifera",
            subtitle = "Taux de recombinaison eleve (~19 cM/Mb)",
            x        = "Distance physique (kb)",
            y        = "Moyenne r2"
        ) +
        scale_x_continuous(expand=c(0,5)) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        theme_classic(base_size=12) +
        theme(plot.title=element_text(face="bold", size=14),
              plot.subtitle=element_text(size=10, colour="grey40"))

    ggsave("ld_decay.png", p, width=10, height=6, dpi=300, bg="white")
    ggsave("ld_decay.pdf", p, width=10, height=6)
    cat("LD decay plot sauvegarde\\n")
REOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch ld_decay.png ld_decay.pdf versions.yml
    """
}

// ── Process 6 : FST plot ──────────────────────────────────────────────────────
process PLOT_FST {
    tag "fst_plot"
    label 'process_low'

    publishDir "${params.outdir}/04_population/plots", mode: 'copy'

    container 'quay.io/biocontainers/r-base:4.3.3'

    input:
    path fst_file

    output:
    path "fst_plot.png",  emit: png
    path "fst_plot.pdf",  emit: pdf
    path "versions.yml",  emit: versions

    script:
    """
    export FST_FILE="${fst_file}"

    Rscript - << 'REOF'
    suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

    fst_file <- Sys.getenv("FST_FILE")
    cat("Chargement FST :", fst_file, "\\n")

    df <- read.table(fst_file, header=TRUE, stringsAsFactors=FALSE)
    colnames(df) <- c("CHROM","BIN_START","BIN_END","N_VARIANTS","FST_VAL")
    df <- df[!is.na(df\$FST_VAL) & df\$FST_VAL >= 0, ]
    df\$CHR <- suppressWarnings(as.integer(gsub("[^0-9]","",df\$CHROM)))
    df <- df[!is.na(df\$CHR), ]

    if (nrow(df) == 0) {
        cat("Pas assez de donnees FST — figure vide generee\\n")
        png("fst_plot.png"); plot.new(); dev.off()
        pdf("fst_plot.pdf"); plot.new(); dev.off()
        quit(status=0)
    }

    fst_mean      <- mean(df\$FST_VAL, na.rm=TRUE)
    fst_sd        <- sd(df\$FST_VAL,   na.rm=TRUE)
    fst_threshold <- fst_mean + 3 * fst_sd

    cat(sprintf("FST moyen : %.4f | Seuil outlier : %.4f\\n",
                fst_mean, fst_threshold))

    p <- ggplot(df, aes(x=BIN_START/1e6, y=FST_VAL,
                        colour=as.factor(CHR))) +
        geom_point(size=0.5, alpha=0.6) +
        geom_hline(yintercept=fst_threshold,
                   colour="red", linetype="dashed", linewidth=0.6) +
        facet_wrap(~CHR, scales="free_x", nrow=4) +
        scale_colour_manual(values=rep(c("#2C3E50","#E67E22"), 16)) +
        labs(
            title    = "Differentiation genetique FST — AMM vs autres",
            subtitle = sprintf("FST moyen = %.4f | Seuil outlier = %.4f",
                fst_mean, fst_threshold),
            x = "Position (Mb)",
            y = "FST (Weir & Cockerham)"
        ) +
        theme_classic(base_size=10) +
        theme(legend.position="none",
              strip.text=element_text(face="bold"),
              plot.title=element_text(face="bold", size=12),
              plot.subtitle=element_text(size=9, colour="grey40"))

    ggsave("fst_plot.png", p, width=16, height=10, dpi=300, bg="white")
    ggsave("fst_plot.pdf", p, width=16, height=10)
    cat("FST plot sauvegarde\\n")
REOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch fst_plot.png fst_plot.pdf versions.yml
    """
}
