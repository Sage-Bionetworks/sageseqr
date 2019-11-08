Sys.setenv(R_CONFIG_ACTIVE = "default")

processing_plan <- drake_plan(
  raw_data = targets("config.yml"),
  clean_md = clean_covariates(metadata, factors = c("batch", "sex", "race", "spanish", "cogdx", "diagnosis","apoe4"),
                              continuous = c("rincontinuous", "age_death", "pmi", "educ", "pct_pf_reads_aligned", "pct_coding_bases", "pct_intergenic_bases", "pct_intronic_bases", "pct_ribosomal_bases"),
                              sample_variable = c("sampleid")),
  covar_correlation = CovariateAnalysis::getAssociationStatistics(clean_md, PVAL = 0.05),
  bioMart_query = get_biomart(count_matrix, organism = "hsa", host = "dec2016.archive.ensembl.org")

)
config <- drake_config(processing_plan)
drake::vis_drake_graph(config)
