#'
#'Detect file type and download data
#'
#'This function takes synIds and version number to download any rectangular file type from Synapse.
#'
#'@param synID A character vector with a Synapse Id.
#'@param version Optional. A numeric vector with the Synapse file version number.
#'@export
#'@return A tibble.
#'@examples
#'\dontrun{
#'file <- get_data(synID = "syn1234", version = 7)
#'
#'}
get_data <- function(synID, version = NULL){
  df <- as_tibble(data.table::fread(synapser::synGet(synID,
                                                     version = as.numeric(version)
                                                     )$path))
  df
}
#'
#'Create covariate matrix from tidy metadata data frame.
#'
#'This function takes a tidy format. Coerces objects to correct type.
#'
#'@param md
#'@param factors A vector of factor variables
#'@param continuous A vector of continuous variables
#' Remove primary variable for now - A vector of primary variables
#'@param sample_variable A character vector with the column name corresponding to sample identifiers
#'
#'@export
#'@return A data frame with coerced variables.
#'@examples
clean_covariates <- function(md, factors, continuous, sample_variable){
  md[, factors] <- lapply(md[, factors], factor)
  md[, continuous] <- lapply(md[, continuous], as.numeric)
  md
}
#'
#'Explore metadata by variable.
#'
#'This function produces boxplots from the variables provided.
#'
#'@param md
#'@param vars A vector of variables to visualize
#'
#'@export
#'@return
#'@examples
#' TO DO: test this
boxplot_vars <- function(md, include_vars, x_var){
  md %>%
    dplyr::select(., c(include_vars, x_var)) %>%
    gather(key, value, -x_var) %>%
    ggplot(aes(x = x_var, y = value)) +
    geom_boxplot() +
    theme(legend.position = "top") +
    facet_grid(key ~ !! x_var, scales = "free")
}
#'Compute biomaRt object
#'
#'Get GC content, gene IDs, gene symbols, gene biotypes and other metadata from ENSEMBL BioMart.
#'
#'@param count_matrix
#'@param host An optional character vector specifying the release version. This specification is highly recommended for a reproducible workflow!
#'@param organism A character vector of the organism name.
#'
get_biomart <- function(count_matrix, host = NULL, organism){
    id_type <- "ensembl_gene_id"
    id <- rownames(count_matrix)
    message("Connecting to BioMart ...")
    ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = host)
    ds <- listDatasets(ensembl)[, "dataset"]
    ds <- grep(paste0("^", org), ds, value = TRUE)
    if (length(ds) == 0){
      stop(paste("Mart not found for:", org))
    } else if (length(ds) > 1) {
      message("Found several marts")
      sapply(ds, function(d) message(paste(which(ds == d), d, sep = ": ")))
      n <- readline(paste0("Choose mart (1-", length(ds), ") : "))
      ds <- ds[as.integer(n)]
    }
    ensembl <- useDataset(ds, mart = ensembl)
    message(paste0("Downloading sequence",
                   ifelse(length(id) > 1, "s", ""), " ..."))
    if (length(id) > 100)
      message("This may take a few minutes ...")
    attrs <- c(id_type, "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end")
    coords <- getBM(filters = id_type, attributes = attrs, values = id, mart = ensembl)
    id <- unique(coords[, id_type])
    coords <- sapply(id, function(i) {
      i.coords <- coords[coords[, 1] == i, 3:5]
      g <- IRanges::GRanges(i.coords[, 1], IRanges::IRanges(i.coords[, 2], i.coords[, 3]))
      return(g)
    })
    coords <- lapply(coords[id], reduce)
    length <- plyr::ldply(coords, function(x) sum(IRanges::width(x)), .id = "ensembl_gene_id") %>%
      dplyr::rename(gne.length = V1)

    gc_content <- getBM(filters = id_type,
                     attributes = c(id_type, "hgnc_symbol", "percentage_gene_gc_content",
                                    "gene_biotype", "chromosome_name"), values = id,
                     mart = ensembl)

    res <- dplyr::full_join(gc_content, length)
    return(res)

}
#' Filter genes
#'
#' Remove genes that have less than 1 counts per million (cpm) in at least 50% of samples per condition. If a biomaRt object
#' is provided, gene lengths and gene GC content is required and genes with missing values are also removed.
#'
#'
#'@param md A data frame with sample identifiers in a column.
#'@param count_matrix A matrix with sample identifiers as columnnames and gene IDs as rownames.
#'@param condition_variables A vector that corresponds to metadata variables.
#'@param biomaRt_parameters An optional data frame to require genes have biomaRt computed GC content and gene lengths.
#'
#'
filter_genes <- function(md, count_matrix, condition_variables, biomaRt_parameters = NULL) {
  if (!is.null(biomaRt_parameters)){
    exists <- biomaRt_parameters #feed a list
    count_matrix <- count_matrix[]
  } else {
    count_matrix
  }
  genes_to_analyze <- md %>%
    plyr::dlply(.(diagnosis), .fun = function(md, counts){
      processed_counts = CovariateAnalysis::getGeneFilteredGeneExprMatrix(counts,
                                                                          MIN_GENE_CPM = 1,
                                                                          MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM = 0.5)
      processed_counts$filteredExprMatrix$genes
    }, count_matrix)


  foo <- unlist(genesToAnalyze) %>%
    unique() %>%
    intersect(GENE.GC$ensembl_gene_id[!is.na(GENE.GC$percentage_gc_content)]) %>%
    intersect(GENE.LEN$ensembl_gene_id[!is.na(GENE.LEN$gene.length)]) %>%
    setdiff(c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"))
  processed = getGeneFilteredGeneExprMatrix(COUNT[genesToAnalyze, ],
                                                   MIN_GENE_CPM=0,
                                                   MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0)
}

