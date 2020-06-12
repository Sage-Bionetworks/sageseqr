#' Execute the drake RNA-seq plan
#'
#' This function wraps the \code{"drake::plan()"}.
#'
#' @param metadata_id Synapse ID to clean metadata file with sample identifiers in a
#' column and variables of interest as column names. There cannot be any missing values.
#' @param metadata_version Optionally, include Synapse file version number. If omitted,
#' current version will be downloaded.
#' @param counts_id Synapse ID to counts data frame with identifiers to the metadata as
#' column names and gene ids in a column.
#' @param counts_version Optionally, include Synapse file version number.
#' @param gene_id_input Column name of the gene ids in the counts_id file.
#' @param factor_input Vector of factor variables. Variables must be present
#' in the metadata as column names.
#' @param continuous_input Vector of continuous variables. Variables must be
#' present in the metadata as column names.
#' @param biomart_id Synapse ID to biomart object.
#' @param biomart_version Optionally, include Synapse file version number.
#' @inheritParams get_biomart
#' @export
rnaseq_plan <- function(metadata_id, metadata_version, counts_id,
                        counts_version, gene_id_input, factor_input,
                        continuous_input, gene_id,
                        biomart_id, biomart_version, host, filters,
                        organism){
  drake::drake_plan(
    import_metadata = get_data(!!metadata_id,
                               !!metadata_version),
    import_counts = get_data(!!counts_id,
                             !!counts_version),
    counts = tibble::column_to_rownames(import_counts,
                                        var = !!gene_id_input),
    clean_md = clean_covariates(md = import_metadata,
                                factors = !!factor_input,
                                continuous = !!continuous_input),
    biomart_results = get_biomart(count_df = counts,
                                gene_id = !!gene_id,
                                synid = !!biomart_id,
                                version = !!biomart_version,
                                filters = !!filters,
                                host = !!host,
                                organism = !!organism),
    filtered_counts = filter_genes(md = clean_md,
                                   count_df = counts),
    cqn_counts = cqn(filtered_counts,
                     biomart_results)
)
}
