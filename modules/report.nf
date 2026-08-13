/*
    MODULE : GWAS Report
    Outil   : R Markdown + knitr
    Rôle    : Génération du rapport scientifique complet HTML/PDF
    Docker  : quay.io/biocontainers/r-base:4.3.3

    CONTENU DU RAPPORT :
    1. Résumé du pipeline (paramètres, versions des outils)
    2. Statistiques QC des reads (depuis MultiQC)
    3. Statistiques d'alignement (flagstat)
    4. Statistiques du VCF (Ts/Tv, nombre de SNPs, MAF spectrum)
    5. Structure de population (PCA + Admixture)
    6. Déclin du LD
    7. Différenciation FST
    8. Résultats GWAS (Manhattan plot + QQ plot + table des top SNPs)
    9. Versions de tous les outils utilisés

    FORMAT :
    HTML  : rapport interactif avec table des matières flottante
    PDF   : version imprimable pour soumission/archivage
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
    path "gwas_report.pdf",  emit: pdf,    optional: true
    path "versions.yml",     emit: versions

    script:
    """
    # Installer les packages R nécessaires
    Rscript -e "
    pkgs <- c('rmarkdown','knitr','ggplot2','dplyr','tidyr',
              'kableExtra','DT','plotly')
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
    if (length(missing) > 0)
        install.packages(missing,
                         repos='https://cran.r-project.org',
                         quiet=TRUE)
    "

    # Copier le template R Markdown depuis le projet
    cp ${projectDir}/report/gwas_report.Rmd .

    # Rendre le rapport HTML
    Rscript -e "
    rmarkdown::render(
        'gwas_report.Rmd',
        output_format = rmarkdown::html_document(
            toc            = TRUE,
            toc_float      = TRUE,
            toc_depth      = 3,
            code_folding   = 'hide',
            theme          = 'flatly',
            highlight      = 'tango',
            self_contained = TRUE,
            fig_width      = 10,
            fig_height     = 6
        ),
        output_file = 'gwas_report.html',
        params = list(
            pipeline_version = '${manifest.version}',
            run_date         = format(Sys.time(), '%Y-%m-%d %H:%M'),
            results_dir      = '.',
            n_samples        = '${params.input}',
            maf_threshold    = ${params.maf},
            gwas_model       = '${params.gwas_model}'
        )
    )
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | awk '{print \$3}')
        rmarkdown: \$(Rscript -e "cat(as.character(packageVersion('rmarkdown')))")
        knitr: \$(Rscript -e "cat(as.character(packageVersion('knitr')))")
    END_VERSIONS
    """

    stub:
    """
    touch gwas_report.html gwas_report.pdf
    touch versions.yml
    """
}
