#' Execute the drake RNA-seq plan
#'
#' This function wraps the \code{"drake::plan()"} and copies the R markdown report to the
#' user's working directory.
#'
#' @param metadata_id Synapse ID to clean metadata file with sample identifiers in a
#' column and variables of interest as column names. There cannot be any missing values.
#' @param metadata_version Optionally, include Synapse file version number. If omitted,
#' current version will be downloaded.
#' @param counts_id Synapse ID to counts data frame with identifiers to the metadata as
#' column names and gene ids in a column.
#' @param counts_version Optionally, include Synapse file version number.
#' @param gene_id_input Column name of the gene ids in the counts_id file.
#' @param sample_id_input Column name of the sample ids in the metadata_id file.
#' @param factor_input Vector of factor variables. Variables must be present
#' in the metadata as column names.
#' @param continuous_input Vector of continuous variables. Variables must be
#' present in the metadata as column names.
#' @param biomart_id Synapse ID to biomart object.
#' @param biomart_version Optionally, include Synapse file version number.
#' @param x_var_for_plot Variable to separate groups for boxplot.
#' @inheritParams plot_sexcheck
#' @inheritParams get_biomart
#' @inheritParams filter_genes
#' @export
rnaseq_plan <- function(metadata_id, metadata_version, counts_id,
                        counts_version, gene_id_input, sample_id_input,
                        factor_input, continuous_input,
                        biomart_id, biomart_version, host, filters,
                        organism, conditions, cpm_threshold = 1,
                        conditions_threshold = 0.5,
                        x_var_for_plot, sex_var){

  # Copies markdown to user's working directory
  if (!file.exists("sageseqr-report.Rmd")) {
    fs::file_copy(system.file("sageseqr-report.Rmd", package = "sageseqr"),
                  new_path = getwd())
  }

  drake::drake_plan(
    import_metadata = get_data(!!metadata_id,
                               !!metadata_version),
    import_counts = get_data(!!counts_id,
                             !!counts_version),
    counts = tibble::column_to_rownames(import_counts,
                                        var = !!gene_id_input),
    clean_md = clean_covariates(md = import_metadata,
                                factors = !!factor_input,
                                continuous = !!continuous_input,
                                sample_identifier = !!sample_id_input),
    biomart_results = get_biomart(count_df = counts,
                                synid = !!biomart_id,
                                version = !!biomart_version,
                                filters = !!filters,
                                host = !!host,
                                organism = !!organism),
    filtered_counts = filter_genes(clean_metadata = clean_md,
                                   count_df = counts,
                                   conditions = !!conditions,
                                   cpm_threshold = 1,
                                   conditions_threshold = 0.5),
    cqn_counts = cqn(filtered_counts,
                     biomart_results),
    boxplots = boxplot_vars(md = clean_md,
                            include_vars = !!continuous_input,
                            x_var = !!x_var_for_plot),
    sex_plot = conditional_plot_sexcheck(clean_md,
                                         counts,
                                         biomart_results,
                                         !!sex_var),
    correlation_plot = get_association_statistics(clean_md),
    report = rmarkdown::render(
      drake::knitr_in("sageseqr-report.Rmd"),
      output_file = drake::file_out(
        !!glue::glue(getwd(), "/sageseqr-report.html")
        ),
      output_dir = "."
      )
    )
}
