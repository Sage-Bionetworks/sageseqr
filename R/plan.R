plan <- drake_plan(
  import_metadata = get_data(config::get("metadata")$synID,
                config::get("metadata")$version),
  import_counts = get_data(config::get("counts")$synID,
                    config::get("counts")$version),
  counts = tibble::column_to_rownames(import_counts, var = config::get("gene id")),
  clean_md = clean_covariates(md = import_metadata, factors = config::get("factors"),
                              continuous = config::get("continuous")),
  #covar_correlation = CovariateAnalysis::getAssociationStatistics(clean_md, PVAL = 0.05),
  geneids = convert_geneids(count_matrix = counts),
  biomart_results = get_biomart(geneids$ensembl_gene_id, host = "uswest.ensembl.org", organism = "hsa"),# Ensembl Release 99 (January 2020)
  filtered_counts = filter_genes(md = clean_md, count_matrix = as.data.frame(counts))
)

config <- drake::drake_config(plan)
drake::vis_drake_graph(config, targets_only = TRUE)


