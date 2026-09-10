
/*
    MODULE : GWAS Report
    Outil   : R Markdown + knitr + pandoc
    Rôle    : Génération du rapport scientifique HTML
    Docker  : rocker/tidyverse:4.3.1  (rmarkdown + pandoc + tidyverse inclus)

    NOTE : manifest.version n'est pas accessible dans les modules.
    On passe la version et les paramètres via des variables d'environnement.

    Le template .Rmd est fourni en entrée du process (pas de ${projectDir}
    dans le script) pour une exécution hermétique.
*/

process GWAS_REPORT {
    tag "gwas_report"
    label 'process_medium'

    publishDir "${params.outdir}/06_report", mode: 'copy'

    container 'rocker/tidyverse:4.3.1'

    input:
    path(all_results)
    path(rmd)

    output:
    path "gwas_report.html", emit: html
    path "versions.yml",     emit: versions

    script:
    def maf_val  = params.maf
    def gwas_mod = params.gwas_model
    """
    export PIPELINE_VERSION="1.0.0"
    export RUN_DATE=\$(date '+%Y-%m-%d %H:%M')
    export MAF_THRESHOLD="${maf_val}"
    export GWAS_MODEL="${gwas_mod}"
    export RESULTS_DIR="."

    # rocker/tidyverse fournit déjà rmarkdown, knitr, ggplot2, dplyr, tidyr
    # et pandoc — aucune installation à l'exécution.
    Rscript - << 'REOF'
    pipeline_version <- Sys.getenv("PIPELINE_VERSION", "1.0.0")
    run_date         <- Sys.getenv("RUN_DATE")
    maf_threshold    <- as.numeric(Sys.getenv("MAF_THRESHOLD", "0.05"))
    gwas_model       <- Sys.getenv("GWAS_MODEL", "lmm")
    results_dir      <- Sys.getenv("RESULTS_DIR", ".")

    rmarkdown::render(
        "${rmd}",
        output_format = rmarkdown::html_document(
            toc            = TRUE,
            toc_float      = TRUE,
            toc_depth      = 3,
            code_folding   = "hide",
            theme          = "flatly",
            highlight      = "tango",
            self_contained = TRUE
        ),
        output_file = file.path(getwd(), "gwas_report.html"),
        knit_root_dir = getwd(),
        params = list(
            pipeline_version = pipeline_version,
            run_date         = run_date,
            results_dir      = results_dir,
            maf_threshold    = maf_threshold,
            gwas_model       = gwas_model
        )
    )
    cat("Rapport HTML généré : gwas_report.html\\n")
REOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
        rmarkdown: \$(Rscript -e "cat(as.character(packageVersion('rmarkdown')))")
        pandoc: \$(pandoc --version | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch gwas_report.html versions.yml
    """
}
