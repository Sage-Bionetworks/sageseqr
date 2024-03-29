library(config)
library(sageseqr)
library(targets)
# build plan
list(
  tar_target(
    import_metadata,
    get_data(
      get("metadata")$synID,
      get("metadata")$version
    )
  ),
  tar_target(
    import_counts,
    get_data(
      get("counts")$synID,
      get("counts")$version
    )
  ),
  tar_target(
    counts,
    tibble::column_to_rownames(
      import_counts,
      var = get("counts")$`gene id`
    )
  ),
  tar_target(
    clean_md,
    clean_covariates(
      md = import_metadata,
      factors = get("factors"),
      continuous = get("continuous"),
      sample_identifier = get("metadata")$`sample id`
    )
  ),
  tar_target(
    check,
    check_mismatch(
      md = clean_md,
      count_df = counts
    )
  ),
  tar_target(
    biomart_results,
    get_biomart(
      count_df = counts,
      synid = get("biomart")$synID,
      version = get("biomart")$version,
      filters = get("biomart")$filters,
      host = get("biomart")$host,
      organism = get("biomart")$organism,
      isexon = get("biomart")$`exon only`,
      custom = get("biomart")$`custom build`,
      cores = get("cores"),
      gtfID = get("biomart")$gtfID,
      gtfVersion = get("biomart")$gtfVersion,
      fastaID = get("biomart")$fastaID,
      fastaVersion = get("biomart")$fastaVersion
    )
  ),
  tar_target(
    filtered_counts,
    filter_genes(
      clean_metadata = clean_md,
      count_df = counts,
      conditions = get("conditions"),
      cpm_threshold = get("cpm threshold"),
      conditions_threshold = get("percent threshold")
    )
  ),
  tar_target(
    biotypes,
    summarize_biotypes(
      filtered_counts,
      biomart_results
    )
  ),
  tar_target(
    cqn_counts,
    cqn(
      filtered_counts,
      biomart_results
    )
  ),
  tar_target(
    gene_coexpression,
    plot_coexpression(
      cqn_counts
    )
  ),
  tar_target(
    boxplots,
    boxplot_vars(
      md = clean_md,
      include_vars = get("continuous"),
      x_var = get("x_var")
    )
  ),
  tar_target(
    sex_plot,
    conditional_plot_sexcheck(
      clean_md,
      counts,
      biomart_results,
      get("sex check")
    )
  ),
  tar_target(
    sex_plot_pca,
    plot_sexcheck_pca(
      clean_md,
      counts,
      biomart_results,
      get("sex check")
    )
  ),
  tar_target(
    correlation_plot,
    get_association_statistics(
      clean_md
    )
  ),
  tar_target(
    significant_covariates_plot,
    run_pca_and_plot_correlations(
      cqn_counts$E,clean_md
    )
  ),
  tar_target(
    outliers,
    identify_outliers(
      filtered_counts,
      clean_md,
      get("dimensions")$color,
      get("dimensions")$shape,
      get("dimensions")$size
    )
  ),
  tar_target(
    dropped,
    dropped_genes(
      filtered_counts = filtered_counts,
      cqn_counts = cqn_counts$E
    )
  ),
  tar_target(
    model,
    stepwise_regression(
      clean_md,
      add_model = get("force null model"),
      primary_variable = get("x_var"),
      cqn_counts = cqn_counts,
      skip = get("skip model"),
      random_effect = get("random_effect")
    )
  ),
  tar_target(
    selected_model,
    if(is.null(get("force model with"))) {
      model$variables_in_model
    } else {
      get("force model with")
    }
  ),
  tar_target(
    residualized_counts,
    purrr::map(
      get("de contrasts"),
      function(x) compute_residuals(
        clean_md,
        filtered_counts,
        dropped,
        random_effect = get("random_effect"),
        cqn_counts = cqn_counts$E,
        primary_variable = x$primary,
        is_num = x$is_numeric_int,
        num_var = x$numeric,
        model_variables = selected_model,
        cores = get("cores")
        )
      )
    ),
  tar_target(
    de,
    wrap_de(
      conditions = get("de contrasts"),
      filtered_counts,
      cqn_counts$E,
      random_effect = get("random_effect"),
      clean_md,
      dropped,
      biomart_results,
      p_value_threshold = get("de p-value threshold"),
      fold_change_threshold = get("de FC"),
      model_variables = selected_model,
      cores = get("cores")
    )
  ),
  tar_target(
    get_gene_list,
    purrr::map(de, function(x) x$differential_expression)
  ),
  tar_target(
    plot_de_volcano,
    purrr::map(
      get_gene_list,
      function(x) plot_volcano(
        x,
        p_value_threshold = get("de p-value threshold"),
        fold_change_threshold = get("de FC"),
        gene_list = get_gene_list
      )
    )
  ),
  tarchetypes::tar_render(
    report,
    "sageseqr-report.Rmd",
    output_file = glue::glue("{getwd()}/{config::get('report')}.html")
  ),
  tar_target(
    document_inputs,
    provenance_helper(
      get("metadata")$synID,
      get("counts")$synID,
      get("metadata")$version,
      get("counts")$version,
      get("biomart")$synID,
      get("biomart")$version
    )
  ),
  tar_target(
    Synapse,
    store_results(
      parent_id = get("store output"),
      cqn_counts = cqn_counts$E,
      clean_md = clean_md,
      filtered_counts = filtered_counts,
      biomart_results = biomart_results,
      de_results = de,
      residualized_counts = residualized_counts,
      report = report,
      rownames = list(
        config::get("metadata")$`sample id`,
        config::get("counts")$`gene id`,
        config::get("biomart")$filters,
        config::get("counts")$`gene id`
      ),
      syn_names = as.list(
        c(
          "Covariates",
          "Filtered counts (greater than 1cpm)",
          "BioMart query results",
          "Normalized counts (CQN)",
          glue::glue("Residualized counts ({names(residualized_counts)})"),
          glue::glue("Differential Expression ({names(de)})")
          )
        ),
      data_names = as.list(
        c(
          "clean_md",
          "filtered_counts",
          "biomart_results",
          "cqn_counts",
          names(residualized_counts),
          names(de)
          )
        ),
      inputs = document_inputs,
      activity_provenance = "Analyze RNA-seq data with sageseqr pkg",
      config_file = "config.yml",
      report_name = get("report")
    )
  )
)
