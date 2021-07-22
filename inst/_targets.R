library(sageseqr)

# Login to Synapse. Make a Synapse account and use synapser to login:
# https://r-docs.synapse.org/articles/manageSynapseCredentials.html
synapser::synLogin()

#Set the cofiguration profile here:
Sys.setenv(R_CONFIG_ACTIVE = "INSERT_CONFIG_PROFILE_HERE")

if (!file.exists("sageseqr-report.Rmd")) {
  fs::file_copy(system.file("sageseqr-report.Rmd", package = "sageseqr"),
                new_path = getwd())
}

#Build Plan
list(
  targets::tar_target(
    import_metadata,
    get_data(!!config::get("metadata")$synID,
             !!config::get("metadata")$version)
  ),
  targets::tar_target(
    import_counts,
    get_data(!!config::get("counts")$synID,
             !!config::get("counts")$version)
  ),
  targets::tar_target(
    counts,
    tibble::column_to_rownames(import_counts,
                               var = !!config::get("counts")$`gene id`)
  ),
  targets::tar_target(
    clean_md,
    clean_covariates(md = import_metadata,
                     factors = !!config::get("factors"),
                     continuous = !!config::get("continuous"),
                     sample_identifier = !!config::get("metadata")$`sample id`)
  ),
  targets::tar_target(
    biomart_results,
    get_biomart(count_df = counts,
                synid = !!config::get("biomart")$synID,
                version = !!config::get("biomart")$version,
                filters = !!config::get("biomart")$filters,
                host = !!config::get("biomart")$host,
                organism = !!config::get("biomart")$organism)
  ),
  targets::tar_target(
    filtered_counts,
    filter_genes(clean_metadata = clean_md,
                 count_df = counts,
                 conditions = !!config::get("conditions"),
                 cpm_threshold = !!config::get("cpm threshold"),
                 conditions_threshold = !!config::get("percent threshold"))
  ),
  targets::tar_target(
    biotypes,
    summarize_biotypes(filtered_counts, biomart_results)
  ),
  targets::tar_target(
    cqn_counts,
    cqn(filtered_counts, biomart_results)
  ),
  targets::tar_target(
    gene_coexpression,
    plot_coexpression(cqn_counts)
  ),

  targets::tar_target(
    boxplots,
    boxplot_vars(md = clean_md,
                 include_vars = !!config::get("continuous"),
                 x_var = !!config::get("x_var"))
  ),
  targets::tar_target(
    sex_plot,
    conditional_plot_sexcheck(clean_md,
                              counts,
                              biomart_results,
                              !!config::get("sex check"))
  ),
  targets::tar_target(sex_plot_pca,
                      plot_sexcheck_pca(
                        clean_md,
                        counts,
                        biomart_results,
                        !!config::get("sex check"))
  ),
  targets::tar_target(
    correlation_plot,
    get_association_statistics(clean_md)
  ),
  targets::tar_target(
    significant_covariates_plot,
    run_pca_and_plot_correlations(cqn_counts$E,clean_md)
  ),
  targets::tar_target(
    outliers,
    identify_outliers(filtered_counts,
                      clean_md,
                      !!config::get("dimensions")$color,
                      !!config::get("dimensions")$shape,
                      !!config::get("dimensions")$size
                     )
  ),
  targets::tar_target(
    model,
    stepwise_regression(
      clean_md,
      primary_variable = !!config::get("x_var"),
      cqn_counts = cqn_counts,
      skip = !!config::get("skip model")
    )
  ),
  targets::tar_target(
    report,
      tarchetypes::tar_render(name="sageseqr-report.Rmd",
                              path=glue::glue("{getwd()}",
                                              '/inst/sageseqr-report.Rmd')
      )
  ),
  targets::tar_target(
    document_inputs,
    provenance_helper(
      !!config::get("metadata")$synID, !!config::get("counts")$synID,
      !!config::get("metadata")$version, !!config::get("counts")$version,
      !!config::get("biomart")$synID, !!config::get("biomart")$version
    )
  ),
  targets::tar_target(
    Synapse,
    store_results(
      parent_id = !!config::get("store output"),
      cqn_counts = cqn_counts$counts,
      clean_md = clean_md,
      filtered_counts = filtered_counts,
      biomart_results = biomart_results,
      rownames = !!list(
        config::get("metadata")$`sample id`,
        config::get("counts")$`gene id`,
        config::get("biomart")$filters,
        config::get("counts")$`gene id`
      ),
      syn_names = list("Covariates", "Filtered counts (greater than 1cpm)",
                       "BioMart query results", "Normalized counts (CQN)"),
      data_names = list("clean_md", "filtered_counts", "biomart_results",
                        "cqn_counts"),
      inputs = document_inputs,
      activity_provenance = "Analyze RNA-seq data with sageseqr pkg",
      config_file = "config.yml",
      report_name = !!config::get("report")
    )
  )

)
