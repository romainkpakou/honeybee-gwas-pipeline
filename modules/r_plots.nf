/*
    MODULE : r_plots
    Outil   : R + ggplot2
    Rôle    : Génération des figures publication-ready
    Docker  : quay.io/biocontainers/r-base:4.3.3

    6 PROCESS DANS CE MODULE :
    1. PLOT_MANHATTAN  — Manhattan plot GWAS
    2. PLOT_QQ         — QQ plot avec lambda GC
    3. PLOT_PCA        — Structure de population (PC1 vs PC2/PC3)
    4. PLOT_ADMIXTURE  — Barplot proportions d'ascendance
    5. PLOT_LD_DECAY   — Courbe de déclin du LD
    6. PLOT_FST        — Visualisation de la différenciation FST
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
    """
    # Installer les packages R nécessaires
    Rscript -e "
    pkgs <- c('ggplot2','dplyr','optparse')
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
    if (length(missing) > 0)
        install.packages(missing, repos='https://cran.r-project.org', quiet=TRUE)
    "

    Rscript ${projectDir}/bin/manhattan_plot.R \\
        --input ${gwas_results} \\
        --output manhattan_plot \\
        --threshold ${params.gwas_pval}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
        ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """

    stub:
    """
    touch manhattan_plot.png manhattan_plot.pdf
    touch versions.yml
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
    Rscript -e "
    pkgs <- c('ggplot2','dplyr','optparse')
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
    if (length(missing) > 0)
        install.packages(missing, repos='https://cran.r-project.org', quiet=TRUE)
    "

    Rscript ${projectDir}/bin/qq_plot.R \\
        --input ${gwas_results} \\
        --output qq_plot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch qq_plot.png qq_plot.pdf
    touch versions.yml
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
    Rscript -e "
    pkgs <- c('ggplot2','dplyr','optparse')
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
    if (length(missing) > 0)
        install.packages(missing, repos='https://cran.r-project.org', quiet=TRUE)
    "

    Rscript ${projectDir}/bin/pca_plot.R \\
        --eigenvec ${eigenvec} \\
        --eigenval ${eigenval} \\
        --output pca_plot

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
    Rscript -e "
    pkgs <- c('ggplot2','dplyr','tidyr','optparse','gridExtra')
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
    if (length(missing) > 0)
        install.packages(missing, repos='https://cran.r-project.org', quiet=TRUE)
    "

    Rscript ${projectDir}/bin/admixture_plot.R \\
        --qdir . \\
        --fam ${fam} \\
        --output admixture_plot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch admixture_plot.png admixture_plot.pdf
    touch versions.yml
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
    Rscript -e "
    pkgs <- c('ggplot2','dplyr','optparse')
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
    if (length(missing) > 0)
        install.packages(missing, repos='https://cran.r-project.org', quiet=TRUE)
    "

    Rscript ${projectDir}/bin/ld_decay.R \\
        --input ${ld_file} \\
        --output ld_decay

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch ld_decay.png ld_decay.pdf
    touch versions.yml
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
    Rscript -e "
    pkgs <- c('ggplot2','dplyr','optparse')
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
    if (length(missing) > 0)
        install.packages(missing, repos='https://cran.r-project.org', quiet=TRUE)
    "

    # Script FST inline (visualisation des fenêtres FST par chromosome)
    Rscript - <<'REOF'
    library(ggplot2)
    library(dplyr)

    fst_file <- "${fst_file}"
    df <- read.table(fst_file, header=TRUE, stringsAsFactors=FALSE)
    colnames(df) <- c("CHROM","BIN_START","BIN_END","N_VARIANTS","FST")
    df <- df[!is.na(df$FST) & df$FST >= 0, ]
    df$CHR <- as.integer(gsub("[^0-9]", "", df$CHROM))
    df <- df[!is.na(df$CHR), ]

    fst_mean <- mean(df$FST, na.rm=TRUE)
    fst_sd   <- sd(df$FST, na.rm=TRUE)
    threshold_outlier <- fst_mean + 3 * fst_sd

    cat(sprintf("FST moyen : %.4f\\n", fst_mean))
    cat(sprintf("Seuil outlier (mean + 3SD) : %.4f\\n", threshold_outlier))

    p <- ggplot(df, aes(x=BIN_START/1e6, y=FST,
                        colour=as.factor(CHR))) +
        geom_point(size=0.5, alpha=0.6) +
        geom_hline(yintercept=threshold_outlier,
                   colour="red", linetype="dashed", linewidth=0.6) +
        facet_wrap(~CHR, scales="free_x", nrow=4) +
        scale_colour_manual(
            values=rep(c("#2C3E50","#E67E22"), 16)
        ) +
        labs(
            title    = "Différenciation génétique (FST) — AMM vs autres",
            subtitle = sprintf(
                "FST moyen = %.4f | Seuil outlier (mean+3SD) = %.4f",
                fst_mean, threshold_outlier
            ),
            x = "Position (Mb)",
            y = "FST (Weir & Cockerham)"
        ) +
        theme_classic(base_size=10) +
        theme(
            legend.position = "none",
            strip.text      = element_text(face="bold"),
            plot.title      = element_text(face="bold", size=12),
            plot.subtitle   = element_text(size=9, colour="grey40")
        )

    ggsave("fst_plot.png", p, width=16, height=10, dpi=300, bg="white")
    ggsave("fst_plot.pdf", p, width=16, height=10)
    cat("FST plot sauvegardé\\n")
REOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch fst_plot.png fst_plot.pdf
    touch versions.yml
    """
}
