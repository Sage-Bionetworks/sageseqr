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
#' Create complete URL to Synapse entity
#'
#' This function creates the url reference to entities in Synapse.
#' @inheritParams get_data
#' @export
complete_path <- function(synid, version = NULL) {
  if (is.null(version)) {
    glue::glue("https://www.synapse.org/#!Synapse:{synid}")
  } else {
    glue::glue("https://www.synapse.org/#!Synapse:{synid}.{version}")
  }
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
  sagethemes::import_lato()
  df <- dplyr::select(md, !!include_vars, !!x_var) %>%
    tidyr::pivot_longer(-!!x_var, names_to = "key", values_to = "value")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]],
                                        y = .data$value,
                                        group = .data[[x_var]])) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[x_var]])) +
    ggplot2::facet_wrap(key ~ ., scales = "free") +
    sagethemes::scale_fill_sage_d() +
    sagethemes::theme_sage()
  p
}
#' Explore metadata by gene expression on the sex chromosomes.
#'
#' This function plots expression of X and Y marker genes, XIST and UTY respectively, and
#' colors each sample by the sex- or gender-specific labeling from the metadata. This is a
#' handy check to determine if samples were swapped or mislabeled.
#'
#' @inheritParams collapse_duplicate_hgnc_symbol
#' @inheritParams filter_genes
#' @param sex_var Column name of the sex or gender-specific metadata.
#' @export
plot_sexcheck <- function(clean_metadata, count_df, biomart_results, sex_var) {
  md <- tibble::rownames_to_column(clean_metadata, var = "sampleId") %>%
    dplyr::select(.data$sampleId, !!sex_var)
  counts <- tibble::rownames_to_column(count_df, var = "geneId")
  results <- dplyr::select(biomart_results, .data$hgnc_symbol, .data$chromosome_name)
  results <- dplyr::filter(results, .data$hgnc_symbol %in% c("XIST", "UTY"))
  results <- tibble::rownames_to_column(results, var = "geneId")
  results <- dplyr::left_join(results, counts)
  results <- tidyr::pivot_longer(results,
                                 -dplyr::all_of(c("geneId", "chromosome_name", "hgnc_symbol")),
                                 names_to = "sampleId",
                                 values_to = "counts(log)") %>%
    dplyr::mutate(`counts(log)` = log(.data$`counts(log)`),
                  `counts(log)` = ifelse(.data$`counts(log)` == -Inf, 0, .data$`counts(log)`))
  results <- dplyr::left_join(results, md, "sampleId")
  results <- tidyr::spread(results, key = .data$hgnc_symbol, value = .data$`counts(log)`) %>%
    dplyr::mutate(UTY = ifelse(is.na(.data$UTY), 0, .data$UTY),
                  XIST = ifelse(is.na(.data$XIST), 0, .data$XIST))
  p <- ggplot2::ggplot(results, ggplot2::aes(x = .data$XIST, y = .data$UTY)) +
    ggplot2::geom_point(ggplot2::aes(color = .data[[sex_var]])) +
    sagethemes::scale_color_sage_d() +
    sagethemes::theme_sage()
  p <- list(plot = p,
            sex_specific_counts = results)
  p
}
#' Conditionally wrap plot_sexcheck for drake
#'
#' Work around to expose plot_sexcheck to testing and export but also leverage
#' drakes function for skipping targets conditionally (see \code{"drake::cancel_if()"}).
#' @inheritParams plot_sexcheck
#' @export
conditional_plot_sexcheck <- function(clean_metadata, count_df, biomart_results, sex_var) {
  drake::cancel_if(is.null(sex_var))
  plot_sexcheck(clean_metadata, count_df, biomart_results, sex_var)
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
  if (class(conditions) == "list") {
    conditions <- unique(conditions[[1]])
  } else {
    conditions <- unique(conditions)
  }

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
#' Formula Builder
#'
#' If a linear mixed model is used, all categorical variables must be
#' modeled as a random effect. This function identifies the variable class
#' to scale numeric variables and model categorical variables as a random
#' effect. Additionally, variables are scaled to account for
#' multiple variables that might have an order of magnitude difference.
#'
#' @param model_variables Optional. Vector of variables to include in the linear (mixed) model.
#' If not supplied, the model will include all variables in \code{md}.
#' @param primary_variable Vector of variables that will be collapsed into a single
#' fixed effect interaction term.
#' @inheritParams coerce_factors
#' @export
build_formula <- function(md, model_variables = NULL, primary_variable) {

  if (!(all(purrr::map_lgl(md, function(x) inherits(x, c("numeric", "factor")))))) {
    stop("Use sageseqr::clean_covariates() to coerce variables into factor and numeric types.")
  }
  # Update metadata to reflect variable subset
  if (!is.null(model_variables)) {
    vars <- c(model_variables, primary_variable)
    md <- dplyr::select(md, dplyr::all_of(vars))
  }
  # Variables of factor or numeric class are required
    col_type <- dplyr::select(md, -primary_variable) %>%
      dplyr::summarise_all(class) %>%
      tidyr::pivot_longer(tidyr::everything(), names_to = "variable", values_to = "class")

  # Categorical or factor variables are modeled as a random effect by (1|variable)
  # Numeric variables are scaled to account for when the spread of data values differs
  # by an order of magnitude
  formula <- purrr::map(1:length(col_type$class), function(i){
    switch(col_type$class[i],
           "factor" = glue::glue("(1|", col_type$variable[i], ")"),
           "numeric" = glue::glue("scale(", col_type$variable[i], ")")
    )
  })

  # List with multiple values can cause glue_collapse to fail. This step is conservative
  # as it unlists all lists at this step.
  if (class(primary_variable) == "list") {
    primary_variable <- primary_variable[[1]]
  }

  # A new categorical is created to model an interaction between two variables
  interaction_term <- glue::glue_collapse({primary_variable}, "_")

  object <- list(metadata = md %>%
                   tidyr::unite(!!interaction_term, dplyr::all_of(primary_variable), sep = "_"),
                formula = formula(
                  glue::glue("~ {interaction_term}+",
                             glue::glue_collapse(formula, sep = "+")
                             )
                  ),
                formula_non_intercept = formula(
                  glue::glue("~ 0+{interaction_term}+",
                             glue::glue_collapse(formula, sep = "+")
                             )
                  ),
                primary_variable = interaction_term,
                variables = unlist(formula)
  )

  # Resolve dropped class type of factor without losing samples as rownames
  object$metadata[[interaction_term]] <- as.factor(object$metadata[[interaction_term]])

  object
}
#' Differential Expression with Dream
#'
#' Differential expression testing is performed with \code{"variancePartition::dream()"} to increase power
#' and decrease false positives for repeated measure designs. The `primary_variable` is modeled as a fixed
#' effect. If you wish to test an interaction term, `primary_variable` can take multiple variable names.
#'
#' @param cqn_counts A counts data frame normalized by CQN.
#' @inheritParams cqn
#' @inheritParams coerce_factors
#' @inheritParams build_formula
#' @inheritParams cqn
#' @export
differential_expression <- function(filtered_counts, cqn_counts, md, model_variables = NULL,
                                    primary_variable, biomart_results) {
  metadata_input <- build_formula(md, model_variables, primary_variable)
  gene_expression <- edgeR::DGEList(filtered_counts)
  gene_expression <- edgeR::calcNormFactors(gene_expression)
  voom_gene_expression <- variancePartition::voomWithDreamWeights(counts = gene_expression,
                                                                  formula = metadata_input$formula,
                                                                  data = metadata_input$metadata)
  voom_gene_expression$E <- cqn_counts

  de_conditions <- levels(metadata_input$metadata[[metadata_input$primary_variable]])

  conditions_for_contrast <- purrr::map_chr(de_conditions,
                                        function(x) glue::glue({metadata_input$primary_variable}, x)
                                        )

  setup_coefficients <- gtools::combinations(n = length(conditions_for_contrast),
                                             r = 2,
                                             v = conditions_for_contrast
                                             )

  contrasts <- lapply(seq_len(nrow(setup_coefficients)),
                      function(i) variancePartition::getContrast(exprObj = voom_gene_expression,
                                                                 formula = metadata_input$formula_non_intercept,
                                                                 data = metadata_input$metadata,
                                                                 coefficient = as.vector(setup_coefficients[i,])
                                                                 )
                      )

  contrasts <- as.data.frame(contrasts)

  df <- as.data.frame(setup_coefficients)

  names(contrasts) <- stringr::str_glue_data(df, "{V1} - {V2}")

  fit_contrasts = variancePartition::dream(exprObj = voom_gene_expression,
                                           formula = metadata_input$formula_non_intercept,
                                           data = metadata_input$metadata,
                                           L = contrasts
                                           )

  de <- lapply(names(contrasts), function(i, fit){
    genes <- limma::topTable(fit, coef = i, number = Inf, sort.by = "logFC")
    genes <- tibble::rownames_to_column(genes, var = "ensembl_gene_id")
  }, fit_contrasts)

  names(de) <- names(contrasts)

  de <- data.table::rbindlist(de, idcol = "Comparison") %>%
    dplyr::mutate(Comparison = gsub(metadata_input$primary_variable, "", .data$Comparison),
                  Direction = .data$logFC/abs(.data$logFC),
                  Direction = factor(.data$Direction, c(-1,1), c("-1" = "down", "1" = "up")),
                  Direction = ifelse(.data$`adj.P.Val` > 0.5 | abs(.data$logFC) < log2(1.2),
                                     "none",
                                     .data$Direction),
                  Direction = as.character(.data$Direction)
    )

  # Join metadata from biomart to differential expression results
  biomart_results <- tibble::rownames_to_column(biomart_results, var = "ensembl_gene_id")

  de <- dplyr::left_join(de, biomart_results, by = c("ensembl_gene_id"))

  object <- list(voom_object = voom_gene_expression,
                 contrasts_to_plot = contrasts,
                 fits = fit_contrasts,
                 differential_expression = de,
                 primary_variable = metadata_input$primary_variable,
                 formula = metadata_input$formula_non_intercept)

  object
}
#' Wrapper for Differential Expression Analysis
#'
#' This function will pass multiple conditions to test to \code{"sagseqr::differential_expression()"}.
#'
#' @param conditions A list of conditions to test as `primary_variable`
#' in \code{"sagseqr::differential_expression()"}.
#' @inheritParams differential_expression
#' @export
wrap_de <- function(conditions, filtered_counts, cqn_counts, md, model_variables, biomart_results) {
  purrr::map(conditions,
             function(x) differential_expression(filtered_counts, cqn_counts, md, model_variables,
                                                 primary_variable = x, biomart_results))

}
#' Stepwise Regression
#'
#' This function performs multivariate forward stepwise regression evaluated by multivariate Bayesian Information
#' Critera (BIC) by wrapping \code{"mvIC::mvForwardStepwise()"}.
#'
#' @inheritParams differential_expression
#' @inheritParams build_formula
#' @export
stepwise_regression <- function(md, model_variables = NULL, primary_variable, cqn_counts) {
  metadata_input <- build_formula(md, model_variables, primary_variable)
  model <- mvIC::mvForwardStepwise(exprObj = cqn_counts$E,
                                   baseFormula = metadata_input$formula,
                                   data = metadata_input$metadata,
                                   variables = array(metadata_input$variables)
  )
  model
}
