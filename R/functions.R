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
#'@param md  A data frame with sample identifiers in a column and relevant experimental covariates.
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
#'This function takes a tidy format. Coerces vectors to correct type. Only include
#'covariates that have 2 or more levels.
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
  } else if (any(!(factors %in% colnames(md))) | any(!(all(continuous %in% colnames(md))))) {
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
  if (is.null(config::get("biomart")$synID)) {
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
  } else {
    biomart_results <- get_data(config::get("biomart")$synID,
             config::get("biomart")$version) %>%
      tibble::column_to_rownames(var = config::get("biomart")$`gene id`)
    required_variables <- c("gene_length", "percentage_gene_gc_content")
    if (!(required_variables %in% colnames(biomart_results))) {
      vars <- glue::glue_collapse(setdiff(required_variables, colnames(biomart_results)),
                                  sep = ", ",
                                  last = " and ")
      message(glue::glue("Warning: {vars} missing from biomart object.
                         This information is required for Conditional
                         Quantile Normalization"))
    }
    return(biomart_results)
  }
}
#'Filter genes
#'
#'Remove genes that have less than 1 counts per million (cpm) in at least 50%
#'of samples per condition. If a biomaRt object is provided, gene lengths and
#'gene GC content is required and genes with missing values are also removed.
#'
#'@inheritParams md
#'@param count_df A counts data frame with sample identifiers as rownames.
#'
filter_genes <- function(md, count_df) {
  genes_to_analyze <- md %>%
    plyr::dlply(.(diagnosis), .fun = function(md, counts){
      processed_counts <- CovariateAnalysis::getGeneFilteredGeneExprMatrix(counts,
                                                                          MIN_GENE_CPM = 1,
                                                                          MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM = 0.5)
      processed_counts$filteredExprMatrix$gen
    }, count_df)
  genes_to_analyze <- unlist(genes_to_analyze) %>%
    unique()
  processed_counts <- CovariateAnalysis::getGeneFilteredGeneExprMatrix(count_df[genes_to_analyze, ],
                                                                       MIN_GENE_CPM = 0,
                                                                       MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM = 0)
  # Convert transcript Ids to gene Ids in counts and gene list with convert_geneids()
  processed_counts$filteredExprMatrix$genes <- convert_geneids(processed_counts$filteredExprMatrix$counts)
  rownames(processed_counts$filteredExprMatrix$counts) <- convert_geneids(processed_counts$filteredExprMatrix$counts)$ensembl_gene_id

  processed_counts
}
#'Get gene Ids
#'
#'@inheritParams count_df
#'
convert_geneids <- function(count_df) {
  if (all(grepl("\\.", rownames(count_df)))) {
    geneids <- tibble(ids = rownames(count_df)) %>%
      tidyr::separate(ids, c("ensembl_gene_id", "position"), sep = "\\.")
    geneids
  } else {
    geneids <- rownames(count_df)
  }
}
#' Conditional Quantile Normalization (CQN)
#'
#' Normalize counts by CQN. By providing a biomart object, the systematic effect of GC content
#' is removed and gene length (in bp) variation is accounted for. Genes with missing GC content
#' or gene lengths will be removed from the counts matrix.
#'
#' @param filtered_counts
#' @param biomart_results
#'
cqn <- function(filtered_counts, biomart_results) {

  required_variables <- c("gene_length", "percentage_gene_gc_content")

  if (!all(required_variables %in% colnames(biomart_results))) {
    message(glue::glue("Error:{setdiff(required_variables,
                       colnames(biomart_results))} missing from
                       biomart object. This information is required
                       for Conditional Quantile Normalization"))
  } else {

    counts <- filtered_counts$filteredExprMatrix$counts

    genes_to_analyze <- intersect(rownames(counts), rownames(biomart_results))

    to_normalize <- subset(counts, rownames(counts) %in% genes_to_analyze)
    gc_length <- subset(biomart_results, rownames(biomart_results) %in% genes_to_analyze)

    normalized_counts <- cqn::cqn(to_normalize,
                                  x = gc_length[, "percentage_gene_gc_content"],
                                  lengths = gc_length[, "gene_length"],
                                  lengthMethod = "smooth",
                                  verbose = FALSE
    )

    normalized_counts$E <- normalized_counts$y + normalized_counts$offset

    return(normalized_counts)

  }


}
