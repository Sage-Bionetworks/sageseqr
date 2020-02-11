plan <- drake_plan(
  metadata = get_data(config::get("metadata")$synID,
                config::get("counts")$synID),
  counts = get_data(config::get("counts")$synID,
                    config::get("counts")$version),
  clean_md = clean_covariates(md = metadata, factors = config::get("factors"),
                              continuous = config::get("continuous")),
  covar_correlation = CovariateAnalysis::getAssociationStatistics(clean_md, PVAL = 0.05),
  geneids = convert_geneids(count_matrix = counts),
  biomart_results = biomart_obj(geneids$ensembl_gene_id, host = "sep2019.archive.ensembl.org", organism = "hsa"),
  filtered_counts = filter_genes(md = clean_md, count_matrix = counts)
)

config <- drake::drake_config(plan)
drake::vis_drake_graph(config, targets_only = TRUE)


