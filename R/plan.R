#' Execute the drake RNA-seq plan
#'
#' This function wraps the \code{"drake::plan()"} and copies the R markdown
#' report to the user's working directory.
#'
#' @param metadata_id Synapse ID to clean metadata file with sample identifiers
#' in a column and variables of interest as column names. There cannot be any
#' missing values.
#' @param metadata_version Optionally, include Synapse file version number. If
#' omitted, current version will be downloaded.
#' @param counts_id Synapse ID to counts data frame with identifiers to the
#' metadata as column names and gene ids in a column.
#' @param counts_version Optionally, include Synapse file version number.
#' @param gene_id_input Column name of the gene ids in the counts_id file.
#' @param sample_id_input Column name of the sample ids in the metadata_id file.
#' @param factor_input Vector of factor variables. Variables must be present
#' in the metadata as column names.
#' @param continuous_input Vector of continuous variables. Variables must be
#' present in the metadata as column names.
#' @param biomart_id Synapse ID to biomart object.
#' @param biomart_version Optionally, include Synapse file version number.
#' @param primary_variable Baseline variable for model selection and variable to
#'  stratify groups in the boxplot.
#' @param report_name Name of output markdown file.
#' @param de_contrasts A list of variables to specify which groups to compare.
#' @inheritParams plot_sexcheck
#' @inheritParams get_biomart
#' @inheritParams filter_genes
#' @inheritParams identify_outliers
#' @inheritParams store_results
#' @inheritParams prepare_results
#' @inheritParams differential_expression
#' @param force_model Optional. A vector of variables to include in the differential
#' expression model if you want to skip the output of the stepwise regression.
#' @param gene_list Optional. A vector of genes to label in the volcano plot.
#' @param skip_model If TRUE, does not run regression model.
#' @export
rnaseq_plan <- function(metadata_id, metadata_version, counts_id,
                        counts_version, gene_id_input, sample_id_input,
                        factor_input, continuous_input,
                        biomart_id, biomart_version, host, filters,
                        organism, conditions, cpm_threshold = 1,
                        conditions_threshold = 0.5,
                        primary_variable, de_contrasts, log_fold_threshold,
                        p_value_threshold, sex_var, color, shape, size,
                        report_name, skip_model, parent_id,
                        rownames, config_file, force_model, gene_list) {
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
                                   cpm_threshold = !!cpm_threshold,
                                   conditions_threshold = !!conditions_threshold),
    biotypes = summarize_biotypes(filtered_counts, biomart_results),
    cqn_counts = cqn(filtered_counts,
                     biomart_results),
    gene_coexpression = plot_coexpression(cqn_counts),
    boxplots = boxplot_vars(md = clean_md,
                            include_vars = !!continuous_input,
                            x_var = !!primary_variable),
    sex_plot = conditional_plot_sexcheck(clean_md,
                                         counts,
                                         biomart_results,
                                         !!sex_var),
    sex_plot_pca = plot_sexcheck_pca(
      clean_md,
      counts,
      biomart_results,
      !!sex_var),
    correlation_plot = get_association_statistics(clean_md),
    significant_covariates_plot = run_pca_and_plot_correlations(cqn_counts$E,
                                                                clean_md),
    outliers = identify_outliers(filtered_counts, clean_md, !!color, !!shape,
                                 !!size),
     model = stepwise_regression(
       clean_md,
       primary_variable = !!primary_variable,
       cqn_counts = cqn_counts,
       skip = !!skip_model
        ),
    de_1 = target(
      command = wrap_de(
        conditions = !!de_contrasts,
        filtered_counts = filtered_counts,
        cqn_counts = cqn_counts$E,
        md = clean_md,
        biomart_results = biomart_results,
        p_value_threshold = !!p_value_threshold,
        log_fold_threshold = !!log_fold_threshold,
        model_variables = model$variables_in_model
      ),
      trigger = trigger(condition = is.na(!!force_model))
    ),
    de_2 = target(
      command = wrap_de(
        conditions = !!de_contrasts,
        filtered_counts = filtered_counts,
        cqn_counts = cqn_counts$E,
        md = clean_md,
        biomart_results = biomart_results,
        p_value_threshold = !!p_value_threshold,
        log_fold_threshold = !!log_fold_threshold,
        model_variables = !!force_model
      ),
      trigger = trigger(condition = !is.na(!!force_model))
    ),
    report = rmarkdown::render(
      drake::knitr_in("sageseqr-report.Rmd"),
      output_file = drake::file_out(
        !!glue::glue("{getwd()}/{report_name}.html")
        ),
      output_dir = "."
      ),
    document_inputs = provenance_helper(
      !!metadata_id, !!counts_id,
      !!metadata_version, !!counts_version,
      !!biomart_id, !!biomart_version
    ),
    plot_de_volcano = target(
      plot_volcano(
        data,
        p_value_threshold = !!p_value_threshold,
        log_fold_threshold = !!log_fold_threshold,
        gene_list = !!gene_list
      ),
      transform = map(data = c(de_1, de_2))
    ),
    de_store = c(names(de_1), names(de_2)),
    Synapse = target(
      store_results(
        parent_id = !!parent_id,
        cqn_counts = cqn_counts$counts,
        clean_md = clean_md,
        filtered_counts = filtered_counts,
        biomart_results = biomart_results,
        de_results = data,
        rownames = !!rownames,
        syn_names = list("Covariates", "Filtered counts (greater than 1cpm)",
                   "BioMart query results", "Normalized counts (CQN)",
                   glue::glue("Differential Expression ({de_store})")),
        data_names = list("clean_md", "filtered_counts", "biomart_results",
                        "cqn_counts", de_store),
        inputs = document_inputs,
        activity_provenance = "Analyze RNA-seq data with sageseqr pkg",
        config_file = !!config_file,
        report_name = !!report_name
      ),
      transform = map(data = c(de_1, de_2)))
    )
  }
