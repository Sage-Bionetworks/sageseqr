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
  df <- tibble::as_tibble(data.table::fread(synapser::synGet(synid,
                                                     version = as.numeric(version)
                                                     )$path))
  df
}
#'Coerce objects to type factors.
#'
#'@param md  A data frame with sample identifiers in a column and relevant experimental covariates.
#'@param factors A vector of factor variables.
coerce_factors <- function(md, factors) {
  md[, factors] <- lapply(md[, factors, drop = FALSE], factor)
  md
}
#'Coerce objects to type numeric.
#'
#'@inheritParams coerce_factors
#'@param continuous A vector of continuous variables.
#'@importFrom magrittr %>%
coerce_continuous <- function(md, continuous) {
  subset <- md[,continuous]
  test_coercion <- lapply(subset, function(x) class(utils::type.convert(x)))
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
#'covariates that have 2 or more levels. Sample identifiers are stored as rownames.
#'
#'@inheritParams coerce_factors
#'@inheritParams coerce_continuous
#'@param sample_identifier The name of the column with the sample identifiers that map to the gene counts data frame.
#'
#'@export
#'@return A data frame with coerced variables.
#'@examples
#'data <- tibble::tribble(
#'  ~individualID, ~diagnosis, ~RIN,
#'  "ind5436", "control", 7.7,
#'  "ind234", "disease", 7.1
#'  )
#' clean_covariates(data, factors = c("individualID", "diagnosis"),
#' continuous = c("RIN"),
#' sample_identifier = c("individualID"))
clean_covariates <- function(md, factors, continuous, sample_identifier) {
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
    md <- tibble::column_to_rownames(md, var = sample_identifier)
    md
  }
}
#'Explore metadata by variable.
#'
#'This function produces boxplots from the variables provided.
#'
#'@inheritParams coerce_factors
#'@param include_vars A vector of variables to visualize
#'@param x_var Variable to plot on the x-axis.
#'@importFrom rlang .data
#'
#'@export
#'@return A boxplot with mutiple groups defined by the include_vars argument.
boxplot_vars <- function(md, include_vars, x_var) {
  df <- dplyr::select(md, !!include_vars, !!x_var) %>%
    tidyr::pivot_longer(-!!x_var, names_to = "key", values_to = "value")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]],
                                        y = .data$value,
                                        group = .data[[x_var]])) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[x_var]])) +
    ggplot2::facet_wrap(key ~ ., scales = "free")
  p
}
#'Get available Ensembl dataset
#'
#'Helper function to search relative Ensembl datasets by partial organism names.
#'
#'@param organism A character vector of the organism name. This argument takes partial strings. For example,"hsa" will match "hsapiens_gene_ensembl".
#'@param host An optional character vector specifying the release version. This specification is highly recommended for a reproducible workflow. (see \code{"biomaRt::listEnsemblArchives()"})
#'@export
biomart_obj <- function(organism, host) {
  message("Connecting to BioMart ...")
  ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = host)
  ds <- biomaRt::listDatasets(ensembl)[, "dataset"]
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
#' Subset counts data frame
#'
#' The biomaRt query requires only gene identifiers to be passed as input. Some alignment
#' output combines additional metadata with the counts file. This function will remove
#' extraneous rows and is especially meant to address output from the STAR aligner.
#' @inheritParams get_biomart
parse_counts <- function(count_df){
  if (any(grepl("N_", row.names(count_df)))) {
    ind <- which(grepl("N_", row.names(count_df)))
    df_parsed <- count_df[-ind,]
    df_parsed
  } else {
    count_df
  }
}
#'Get Ensembl biomaRt object
#'
#'Get GC content, gene Ids, gene symbols, gene biotypes, gene lengths
#'and other metadata from Ensembl BioMart. Object returned contains gene Ids
#'as rownames.
#'
#'@param count_df A counts data frame with sample identifiers as rownames.
#'@inheritParams get_data
#'@param gene_id Column name of gene Ids
#'@param filters A character vector listing biomaRt query filters.
#'(For a list of filters see \code{"biomaRt::listFilters()"})
#'@param host An optional character vector specifying the release version.
#'This specification is highly recommended for a reproducible workflow.
#'(see \code{"biomaRt::listEnsemblArchives()"})
#'@param organism A character vector of the organism name. This argument
#'takes partial strings. For example,"hsa" will match "hsapiens_gene_ensembl".
#'@importFrom rlang .data
#'@export
get_biomart <- function(count_df, gene_id, synid, version, host, filters, organism) {
  if (is.null(config::get("biomart")$synID)) {
    # Get available datset from Ensembl
    ensembl <- biomart_obj(organism, host)

    # Check for extraneous rows
    count_df <- parse_counts(count_df)

    # Parse gene IDs to use in query
    gene_ids <- convert_geneids(count_df)

    message(paste0("Downloading sequence",
                   ifelse(length(gene_ids) > 1, "s", ""), " ..."))

    if (length(gene_ids) > 100)
      message("This may take a few minutes ...")

    attrs <- c(filters, "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end")
    coords <- biomaRt::getBM(filters = filters,
                             attributes = attrs,
                             values = gene_ids,
                             mart = ensembl,
                             useCache = FALSE)
    gene_ids <- unique(coords[, filters])

    coords <- sapply(gene_ids, function(i) {
      i.coords <- coords[coords[, 1] == i, 3:5]
      g <- GenomicRanges::GRanges(i.coords[, 1], IRanges::IRanges(i.coords[, 2], i.coords[, 3]))
      g
    })

    length <- plyr::ldply(coords[gene_ids], function(x) sum(IRanges::width(x)), .id = "ensembl_gene_id") %>%
      dplyr::rename(gene_length = .data$V1)

    gc_content <- biomaRt::getBM(filters = filters,
                                 attributes = c(filters, "hgnc_symbol", "percentage_gene_gc_content",
                                                "gene_biotype", "chromosome_name"),
                                 values = gene_ids,
                                 mart = ensembl,
                                 useCache = FALSE)

    biomart_results <- dplyr::full_join(gc_content, length)

    # Duplicate Ensembl Ids are collapsed into a single entry
    biomart_results <- collapse_duplicate_hgnc_symbol(biomart_results)

    # Biomart IDs as rownames
    biomart_results <- tibble::column_to_rownames(biomart_results, var = gene_id)

    biomart_results
  } else {
    # Download biomart object from syndID specified in config.yml
    biomart_results <- get_data(synid, version)

    # Biomart IDs as rownames
    biomart_results <- tibble::column_to_rownames(biomart_results, var = gene_id)

    # Gene metadata required for count CQN
    required_variables <- c("gene_length", "percentage_gene_gc_content")

    if (!all(required_variables %in% colnames(biomart_results))) {
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
#'Duplicate HGNC
#'
#'Count normalization requires Ensembl Ids to be unique. In rare cases, there are more
#'than one HGNC symbol per gene Id. This function collapses the duplicate entries into
#'a single entry by appending the HGNC symbols in a comma separated list.
#'@param biomart_results Output of \code{"sageseqr::get_biomart()"}.
#'@importFrom rlang .data
#'
#'@export
collapse_duplicate_hgnc_symbol <- function(biomart_results){
  biomart_results %>%
    dplyr::group_by(.data$ensembl_gene_id) %>%
    dplyr::mutate(hgnc_symbol = paste(.data$hgnc_symbol, collapse = ", ")) %>%
    unique()
}
#' Filter genes
#'
#' Filter genes with low expression. This function is more permissive by setting conditions
#' that corresponds to metadata variables. The gene matrix is split by condition and the
#' counts per million (CPM) for a given condition is computed by \code{"sageseqr::simple_filter()"}.
#'
#' @inheritParams coerce_factors
#' @inheritParams get_biomart
#' @inheritParams simple_filter
#' @param conditions Conditions to bin gene counts that correspond to variables in `md`.
#' @param clean_metadata A data frame with sample identifiers as rownames and variables as
#' factors or numeric as determined by \code{"sageseqr::clean_covariates()"}.
#' @importFrom magrittr %>%
#' @export
filter_genes <- function(clean_metadata, count_df, conditions,
                         cpm_threshold, conditions_threshold) {
  if (!any(conditions %in% colnames(clean_metadata))) {
    stop("Conditions are missing from the metadata.")
  }

  # Check for extraneous rows
  count_df <- parse_counts(count_df)

  split_data <- split(clean_metadata,
                      f = as.list(clean_metadata[, conditions, drop = F]),
                      drop = T)

  map_genes <- purrr::map(split_data, function(x) {
    simple_filter(count_df[, rownames(x), drop = F],
                  cpm_threshold,
                  conditions_threshold)
    })

  genes_to_analyze <- unlist(map_genes) %>%
    unique() %>%
    sort()

  processed_counts <- count_df[genes_to_analyze,]

  # Convert transcript Ids to gene Ids in counts with convert_geneids()
  rownames(processed_counts) <- convert_geneids(processed_counts)

  processed_counts
}
#' Filter genes with low expression
#'
#' `simple_filter` converts a gene counts matrix into counts per million (CPM) and identifies
#' the genes that meet the minimum CPM threshold in a percentage of samples. The minimum CPM
#' threshold and percent threshold is user defined. The function returns a list of genes.
#'
#' @inheritParams get_biomart
#' @param cpm_threshold The minimum number of CPM allowed.
#' @param conditions_threshold Percentage of samples that should contain the minimum CPM.
#' @examples
#'\dontrun{
#'gene_list <- simple_filter(count_df = counts, cpm_threshold = 1, condition_threshold = 0.5)
#'}
simple_filter <- function(count_df, cpm_threshold, conditions_threshold) {
  cpm  <- edgeR::cpm(count_df)
  cpm[is.nan(cpm)] <- 0
  fraction <- rowMeans(cpm >= cpm_threshold)
  keep <- fraction >= conditions_threshold
  genes <- keep[keep]
  genes <- names(genes)
  genes
}
#'Get gene Ids
#'
#'@inheritParams get_biomart
#'@importFrom rlang .data
#'
#'@export
convert_geneids <- function(count_df) {
  if (any(grepl("\\.", rownames(count_df)))) {
    geneids <- tibble::tibble(ids = rownames(count_df)) %>%
      tidyr::separate(.data$ids, c("ensembl_gene_id", "position"), sep = "\\.")
    geneids$ensembl_gene_id
  } else {
    rownames(count_df)
  }
}
#' Conditional Quantile Normalization (CQN)
#'
#' Normalize counts by CQN. By providing a biomart object, the systematic effect of GC content
#' is removed and gene length (in bp) variation is accounted for. Genes with missing GC content
#' or gene lengths will be removed from the counts matrix.
#'
#'@param filtered_counts A counts data frame with genes removed that have low expression.
#'@inheritParams collapse_duplicate_hgnc_symbol
#'@export
cqn <- function(filtered_counts, biomart_results) {

  required_variables <- c("gene_length", "percentage_gene_gc_content")

  if (!all(required_variables %in% colnames(biomart_results))) {
    message(glue::glue("Error:{setdiff(required_variables,
                       colnames(biomart_results))} missing from
                       biomart object. This information is required
                       for Conditional Quantile Normalization"))
  } else {

    genes_to_analyze <- intersect(rownames(filtered_counts), rownames(biomart_results))

    to_normalize <- subset(filtered_counts, rownames(filtered_counts) %in% genes_to_analyze)
    gc_length <- subset(biomart_results, rownames(biomart_results) %in% genes_to_analyze)

    normalized_counts <- suppressWarnings(cqn::cqn(to_normalize,
                                  x = gc_length[, "percentage_gene_gc_content"],
                                  lengths = gc_length[, "gene_length"],
                                  lengthMethod = "smooth",
                                  verbose = FALSE
    )
    )

    normalized_counts$E <- normalized_counts$y + normalized_counts$offset

    return(normalized_counts)
  }
}
#'@importFrom quantreg rq
#'@export
quantreg::rq
#'@importFrom mclust Mclust
#'@export
mclust::Mclust
#'@importFrom mclust mclustBIC
#'@export
mclust::mclustBIC
