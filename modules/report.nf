
/*
    MODULE : GWAS Report
    Outil   : R Markdown + knitr
    Rôle    : Génération du rapport scientifique HTML/PDF
    Docker  : quay.io/biocontainers/r-base:4.3.3

    NOTE : manifest.version n'est pas accessible dans les modules.
    On passe la version via une variable shell depuis main.nf.
*/

process GWAS_REPORT {
    tag "gwas_report"
    label 'process_medium'

    publishDir "${params.outdir}/06_report", mode: 'copy'

    container 'quay.io/biocontainers/r-base:4.3.3'

    input:
    path(all_results)

    output:
    path "gwas_report.html", emit: html
    path "gwas_report.pdf",  emit: pdf,  optional: true
    path "versions.yml",     emit: versions

    script:
    // Passer les params Nextflow comme variables d'environnement
    // pour les rendre accessibles dans le code R
    def maf_val   = params.maf
    def gwas_mod  = params.gwas_model
    """
    export PIPELINE_VERSION="1.0.0"
    export RUN_DATE=\$(date '+%Y-%m-%d %H:%M')
    export MAF_THRESHOLD="${maf_val}"
    export GWAS_MODEL="${gwas_mod}"
    export RESULTS_DIR="."

    # Installer les packages R dans un répertoire local accessible en écriture
    mkdir -p /tmp/Rlibs
    Rscript -e "
    lib_path <- '/tmp/Rlibs'
    .libPaths(c(lib_path, .libPaths()))
    pkgs <- c('rmarkdown','knitr','ggplot2','dplyr','tidyr','kableExtra')
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
    if (length(missing) > 0)
        install.packages(missing, repos='https://cran.r-project.org',
                         lib=lib_path, quiet=TRUE)
    "
    export R_LIBS_USER=/tmp/Rlibs

    # Copier le template R Markdown
    cp ${projectDir}/report/gwas_report.Rmd .

    # Rendre le rapport HTML
    Rscript - << 'REOF'
    # Charger les packages depuis le répertoire local
    .libPaths(c('/tmp/Rlibs', .libPaths()))
    pipeline_version <- Sys.getenv("PIPELINE_VERSION", "1.0.0")
    run_date         <- Sys.getenv("RUN_DATE")
    maf_threshold    <- as.numeric(Sys.getenv("MAF_THRESHOLD", "0.05"))
    gwas_model       <- Sys.getenv("GWAS_MODEL", "lmm")
    results_dir      <- Sys.getenv("RESULTS_DIR", ".")

    rmarkdown::render(
        'gwas_report.Rmd',
        output_format = rmarkdown::html_document(
            toc            = TRUE,
            toc_float      = TRUE,
            toc_depth      = 3,
            code_folding   = 'hide',
            theme          = 'flatly',
            highlight      = 'tango',
            self_contained = TRUE
        ),
        output_file = 'gwas_report.html',
        params = list(
            pipeline_version = pipeline_version,
            run_date         = run_date,
            results_dir      = results_dir,
            maf_threshold    = maf_threshold,
            gwas_model       = gwas_model
        )
    )
    cat("Rapport HTML genere : gwas_report.html\\n")
REOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
        rmarkdown: \$(Rscript -e "cat(as.character(packageVersion('rmarkdown')))")
    END_VERSIONS
    """

    stub:
    """
    touch gwas_report.html gwas_report.pdf versions.yml
    """
}
