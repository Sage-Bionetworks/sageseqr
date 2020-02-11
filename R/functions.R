#'Detect file type and download data
#'
#'This function takes synIds and version number to download any
#'rectangular file type from Synapse.
#'
#'@param synid A character vector with a Synapse Id.
#'@param version Optional. A numeric vector with the Synapse file
#'version number.
#'@export
#'@return A tibble.
#'@examples
#'\dontrun{
#'file <- get_data(synID = "syn1234", version = 7)
#'
#'}
get_data <- function(synid, version = NULL) {
  df <- as_tibble(data.table::fread(synapser::synGet(synid,
                                                     version = as.numeric(version)
                                                     )$path))
  df
}
#'Coerce objects to type factors
#'
#'@inheritParams md
#'@param factors A vector of factor variables
coerce_factors <- function(md, factors) {
  md[, factors] <- lapply(md[, factors, drop = FALSE], factor)
  md
}
#'Coerce objects to type numeric
#'
#'@inheritParams md
#'@param continuous A vector of continuous variables
coerce_continuous <- function(md, continuous) {
  subset <- md[,continuous]
  test_coercion <- lapply(subset, function(x) class(type.convert(x)))
  if (all(test_coercion %in% c("integer", "numeric"))) {
    md[, continuous] <- lapply(md[, continuous, drop = FALSE], function(x) as.numeric(x))
    md
  } else {
    mismatched <- continuous[which(!test_coercion %in% c("integer", "numeric"))]
    stop(glue::glue("Variable {mismatched} can not be coerced to numeric."))
  }
}
#'Create covariate matrix from tidy metadata data frame.
#'
#'This function takes a tidy format. Coerces vectors to correct type.
#'
#'@inheritParams md
#'@inheritParams factors
#'@inheritParams continuous
#'
#'@export
#'@return A data frame with coerced variables.
#'@examples
clean_covariates <- function(md, factors, continuous) {
  if (missing(factors) | missing(continuous)) {
    stop("Factor and continuous variables are required.")
  } else if (length(intersect(factors, continuous)) != 0) {
    stop("Variables are present in both the continuous and factor arguments. Variables must be designated
         numeric or factor, not both.")
  } else if (!(factors %in% colnames(md) | !(continuous %in% colnames(md)))) {
    stop("Variables provided are not present in the metadata.")
  } else {
    md <- coerce_factors(md, factors)
    md <- coerce_continuous(md, continuous)
    md
  }
}
#'
#'Explore metadata by variable.
#'
#'This function produces boxplots from the variables provided.
#'
#'@inheritParams md
#'@param vars A vector of variables to visualize
#'
#'@export
#'@return
#'@examples
#' TO DO: test this
boxplot_vars <- function(md, include_vars, x_var) {
  md %>%
    dplyr::select(., c(include_vars, x_var)) %>%
    gather(key, value, -x_var) %>%
    ggplot(aes(x = x_var, y = value)) +
    geom_boxplot() +
    theme(legend.position = "top") +
    facet_grid(key ~ !! x_var, scales = "free")
}
#'Get available Ensembl dataset
#'
#'Helper function to search relative Ensembl datasets by partial organism names.
#'
#'@param organism A character vector of the organism name. This argument takes partial strings. For example,"hsa" will match "hsapiens_gene_ensembl".
#'@param host An optional character vector specifying the release version. This specification is highly recommended for a reproducible workflow. (see \code{"biomaRt::listEnsemblArchives()"})
biomart_obj <- function(organism, host) {
  message("Connecting to BioMart ...")
  ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = host)
  ds <- listDatasets(ensembl)[, "dataset"]
  ds <- grep(paste0("^", organism), ds, value = TRUE)
  if (length(ds) == 0) {
    stop(paste("Mart not found for:", organism))
  } else if (length(ds) > 1) {
    message("Found several marts")
    sapply(ds, function(d) message(paste(which(ds == d), d, sep = ": ")))
    n <- readline(paste0("Choose mart (1-", length(ds), ") : "))
    ds <- ds[as.integer(n)]
  }
  ensembl <- biomaRt::useDataset(ds, mart = ensembl)
  ensembl
}
#'Get Ensembl biomaRt object

#'Get GC content, gene Ids, gene symbols, gene biotypes, gene lengths
#'and other metadata from Ensembl BioMart.
#'
#'@param gene_ids Ensembl gene Ids. Transcript Ids must be converted to
#'gene Ids.
#'@param host An optional character vector specifying the release version.
#'This specification is highly recommended for a reproducible workflow.
#'(see \code{"biomaRt::listEnsemblArchives()"})
#'@param organism A character vector of the organism name. This argument
#'takes partial strings. For example,"hsa" will match "hsapiens_gene_ensembl".
get_biomart <- function(gene_ids, host, organism) {
    id_type <- "ensembl_gene_id"
    ensembl <- biomart_obj(organism, host)
    message(paste0("Downloading sequence",
                   ifelse(length(gene_ids) > 1, "s", ""), " ..."))
    if (length(gene_ids) > 100)
      message("This may take a few minutes ...")
    attrs <- c(id_type, "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end")
    coords <- biomaRt::getBM(filters = id_type, attributes = attrs, values = gene_ids, mart = ensembl)
    gene_ids <- unique(coords[, id_type])
    coords <- sapply(gene_ids, function(i) {
      i.coords <- coords[coords[, 1] == i, 3:5]
      g <- GenomicRanges::GRanges(i.coords[, 1], IRanges::IRanges(i.coords[, 2], i.coords[, 3]))
      g
    })
    length <- plyr::ldply(coords[gene_ids], function(x) sum(IRanges::width(x)), .id = "ensembl_gene_id") %>%
      dplyr::rename(gene_length = V1)

    gc_content <- getBM(filters = id_type,
                     attributes = c(id_type, "hgnc_symbol", "percentage_gene_gc_content",
                                    "gene_biotype", "chromosome_name"), values = gene_ids,
                     mart = ensembl)

    biomart_results <- dplyr::full_join(gc_content, length)

    #Downstream analysis requires the removal of duplicate ensembl gene Ids,
    #removed HGNC symbol will be printed in the console
    removed_hgnc_symbol <- biomart_results %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::filter(row_number(hgnc_symbol) != 1)
    vec <- glue::glue_collapse(glue::glue("{vec}"), sep = ", ", last = " and ")
    message(glue::glue("The following HGNC symbols were removed due to duplicate gene Ids:{vec}"))

    biomart_results <- biomart_results %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::filter(row_number(hgnc_symbol) == 1)

    biomart_results
}
#'Filter genes
#'
#'Remove genes that have less than 1 counts per million (cpm) in at least 50%
#'of samples per condition. If a biomaRt object is provided, gene lengths and
#'gene GC content is required and genes with missing values are also removed.
#'
#'@param md A data frame with sample identifiers in a column.
#'@param count_matrix A matrix with sample identifiers as column names and gene
#'Ids as rownames.
filter_genes <- function(md, count_matrix) {
  genes_to_analyze <- md %>%
    plyr::dlply(.(diagnosis), .fun = function(md, counts){
      processed_counts <- CovariateAnalysis::getGeneFilteredGeneExprMatrix(counts,
                                                                          MIN_GENE_CPM = 1,
                                                                          MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM = 0.5)
      processed_counts$filteredExprMatrix$gen
    }, count_matrix)
  genes_to_analyze <- unlist(genes_to_analyze) %>%
    unique()
  processed_counts <- CovariateAnalysis::getGeneFilteredGeneExprMatrix(counts[genes_to_analyze, ],
                                                                       MIN_GENE_CPM = 0,
                                                                       MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM = 0)
  #convert transcript Ids to gene Ids with convert_geneids()
  processed_counts$filteredExprMatrix$genes <- convert_geneids(processed_counts$filteredExprMatrix$counts)
  rownames(processed_counts$filteredExprMatrix$counts) <- convert_geneids(filtered_counts$filteredExprMatrix$counts)$ensembl_gene_id
  processed_counts
}
#'Get gene Ids
#'
#'@inheritParams count_matrix
#'
convert_geneids <- function(count_matrix) {
  if (all(grepl("\\.", rownames(count_matrix)))) {
    geneids <- tibble(ids = rownames(count_matrix)) %>%
      tidyr::separate(ids, c("ensembl_gene_id", "position"), sep = "\\.")
    geneids
  } else {
    geneids <- rownames(count_matrix)
  }
}
#' Conditional Quantile Normalization (CQN)
#'
#' CQN library normalization method is applied.
cqn <- function(filtered_counts, biomart_results) {
  normalized_counts <- cqn::cqn(filtered_counts$filteredExprMatrix$counts,
                                x = biomart_results[biomart_results$ensembl_gene_id %in%
                                                      filtered_counts$filteredExprMatrix$genes$ensembl_gene_id,
                                                    "percentage_gene_gc_content"],
                                lengths = biomart_results[biomart_results$ensembl_gene_id %in%
                                                            filtered_counts$filteredExprMatrix$genes$ensembl_gene_id,
                                                          "gene_length"],
                                lengthMethod = "smooth",
                                verbose = FALSE
                                )

}
#'
