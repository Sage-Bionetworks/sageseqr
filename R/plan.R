plan <- drake::drake_plan(
  import_metadata = get_data(metadata_id,
                             metadata_version),
  import_counts = get_data(counts_id,
                    counts_version),
  counts = tibble::column_to_rownames(import_counts, var = gene_id_input),
  clean_md = clean_covariates(md = import_metadata, factors = factor_input,
                              continuous = continuous_input),
  #covar_correlation = CovariateAnalysis::getAssociationStatistics(clean_md, PVAL = 0.05),
  biomart_results = get_biomart(count_df = counts,
                                gene_id = id_input_for_biomart,
                                synid = biomart_id,
                                version = biomart_version,
                                filters = filters_input,
                                host = host_input,
                                organism = organism_input),
  filtered_counts = filter_genes(md = clean_md, count_df = counts),
  cqn_counts = cqn(filtered_counts, biomart_results)
)
