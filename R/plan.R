plan <- drake_plan(
  import_metadata = get_data(config::get("metadata")$synID,
                config::get("metadata")$version),
  import_counts = get_data(config::get("counts")$synID,
                    config::get("counts")$version),
  counts = tibble::column_to_rownames(import_counts, var = config::get("counts")$`gene id`),
  clean_md = clean_covariates(md = import_metadata, factors = config::get("factors"),
                              continuous = config::get("continuous")),
  #covar_correlation = CovariateAnalysis::getAssociationStatistics(clean_md, PVAL = 0.05),
  biomart_results = get_biomart(count_df = counts,
                                gene_id = config::get("biomart")$`gene id`,
                                synid = config::get("biomart")$synID,
                                version = config::get("biomart")$version,
                                filters = config::get("biomart")$filters,
                                host = config::get("biomart")$host,
                                organism = config::get("biomart")$organism),
  filtered_counts = filter_genes(md = clean_md, count_df = counts),
  cqn_counts = cqn(filtered_counts, biomart_results)
)
