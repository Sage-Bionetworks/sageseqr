plan <- drake_plan(
  import_metadata = get_data(config::get("metadata")$synID,
                config::get("metadata")$version),
  import_counts = get_data(config::get("counts")$synID,
                    config::get("counts")$version),
  counts = tibble::column_to_rownames(import_counts, var = config::get("gene id")),
  clean_md = clean_covariates(md = import_metadata, factors = config::get("factors"),
                              continuous = config::get("continuous")),
  #covar_correlation = CovariateAnalysis::getAssociationStatistics(clean_md, PVAL = 0.05),
  geneids = convert_geneids(count_df = counts),
  biomart_results = get_biomart(geneids$ensembl_gene_id, host = "ensembl.org", organism = "hsa"),# Ensembl Release 99 (January 2020)
  filtered_counts = filter_genes(md = clean_md, count_df = counts),
  cqn_counts = cqn(filtered_counts, biomart_results)
)

make_the_plan <- drake::drake_plan(plan)
drake::vis_drake_graph(make_the_plan, targets_only = TRUE)

