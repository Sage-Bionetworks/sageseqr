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
  # login to Synapse
  synapser::synLogin()

  df <- tibble::as_tibble(
    data.table::fread(
      synapser::synGet(
        synid,
        version = as.numeric(version)
        )$path,
      header = TRUE
      )
    )
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
#'Get Gene Length and GC: biomart
#'
#'Get gene length and GC content of `gene` from data frame (`df`) containing
#'gene IDs in a column labeled the value of `column_id`. `df` must have a column
#'labeled sequence where the sequence info for exons/or transcripts are and
#'feature start and stop coordinates in columns labled the value of `start` and
#'`stop`. This will return a vector of `gene`, length, and GC content as a
#'percentage.
#'
#'@param gene a gene name corresponding to a value in `df[,column_id]`
#'@param column_id the column name containg gene feature IDs
#'@param df the dataframe of exon/transcript annotations
#'@param start column name of feature start coordinate
#'@param end column name of feature end coordinate
#'@export
biomart_stats <- function(gene, column_id, df, start, end) {
  df <- df[ df[,column_id] == gene,]
  # Split the sequence into a list object
  seq <- strsplit(df$sequence, '')

  # Name the vector of bases
  for (i in 1:length(seq)) {
    names(seq[[i]]) <- df[i,start]:df[i,end]
  }

  # Unlist and remove overlapping bases
  bases <- unlist(seq)
  bases <- bases[!duplicated(names(bases))]

  # gene length:
  length <- length(bases)

  # Calculate base stats and GC content
  base_freq <- BSgenome::alphabetFrequency(Biostrings::DNAString(
    paste0(bases,collapse = '')
  ))

  gc <- sum(base_freq[c('G','C')])/sum(base_freq[c('A','G','C','T')])*100

  return(c('ID'= gene,
           'gene_length'= length,
           'percentage_gene_gc_content'=gc)
  )
}
#'Get Gene Length and GC: GTF/Fasta
#'
#'Get gene GC content of `gene` from data frame (`df`) containing
#'gene IDs in a column labeled the value of `column_id`. `df` must have a column
#'labeled sequence where the sequence info for exons/or transcripts are and
#'feature start and stop coordinates in columns labled the value of `start` and
#'`stop`. This will return a vector of `gene`, length, and GC content as a
#'percentage.
#'
#'@param i A character of a gene id.
#'@param data a dataframe of exons transcrtipts or genes, imported from a GTF file.
#'one line per-feature, data$V1 corresponds to the featur chormosome, data$V4 to
#' start coordinate, data$V5 to stop coordinate, and data$V7 to the sequence strand.
#'@param genome a list of DNAStringSet objects named by their contig (chromosome)
#'name.
#'@export
gtf_stats <- function( i, data, genome){
  gene_model <- data[data$gene_id == i,]

  gr <- GenomicRanges::reduce(
    GenomicRanges::GRanges(
      paste0(
        names(genome),':',
        gene_model$V4, '-',gene_model$V5, ':',
        gene_model$V7
      )))

  # calculate the length
  len <- sum(GenomicRanges::width(gr))

  ## Get Exon Sequence
  seq <- BSgenome::getSeq(genome, gr)

  ## Get Base Frequency
  alphafreq <- BSgenome::alphabetFrequency(seq)
  totalfreq <- colSums(alphafreq)

  ## get G/C content as a fraction
  gc <- (sum(totalfreq[c('G','C')]) /
           sum(totalfreq[c('A','G','C','T')])) * 100
  return(c(i, gc, len))
}
#'Get Ensembl biomaRt object
#'
#'Get GC content, gene Ids, gene symbols, gene biotypes, gene lengths
#'and other metadata from Ensembl BioMart. Object returned contains gene Ids
#'as rownames.
#'
#'@param count_df A counts data frame with sample identifiers as column names
#'and gene Ids are rownames.
#'@inheritParams get_data
#'@param filters A character vector listing biomaRt query filters.
#'(For a list of filters see \code{"biomaRt::listFilters()"})
#'@param host An optional character vector specifying the release version.
#'This specification is highly recommended for a reproducible workflow.
#'(see \code{"biomaRt::listEnsemblArchives()"})
#'@param organism A character vector of the organism name. This argument
#'takes partial strings. For example,"hsa" will match "hsapiens_gene_ensembl".
#'@param custom Defaults to NULL. If TRUE, the GC and Gene Length parameters will
#'only consider exotic regions.
#'@param gtfID Defaults to NULL. A character vector with a Synapse ID
#'corresponding to a gtf formatted gene annotation file.
#'@param gtfVersion Optional. A numeric vector with the gene GTF Synapse file
#' version number.
#'@param fastaID Defaults to NULL. A character vector with a Synapse ID
#'corresponding to a FASTA formatted genome annotation file.
#'@param fastaVersion Optional. A numeric vector with the genome FASTA Synapse
#'file version number.
#'@param cores An integer of cores to specify in the parallel backend (eg. 4).
#'@param isexon Defaults to FALSE. If TRUE, the GC and Gene Length parameters will
#'only consider exotic regions.
#'@importFrom rlang .data
#'@export
get_biomart <- function(count_df, synid, version, host, filters, organism,
                        custom, gtfID, gtfVersion = NULL, fastaID,
                        fastaVersion = NULL, cores = NULL, isexon = FALSE) {
  if (is.null(synid)) {
    if(is.null(cores)){
      cores = parallel::detectCores()-1
    }

    # Check for extraneous rows
    count_df <- parse_counts(count_df)

    # Parse gene IDs to use in query
    gene_ids <- convert_geneids(count_df)

    # Create Biomart Object
    if (is.null(custom) | isFALSE(custom)) {
      # use biomaRt to construct the biomart_object

      # Get available datset from Ensembl
      ensembl <- biomart_obj(organism, host)

      message(
        paste0(
          "Downloading sequence",
          ifelse(length(gene_ids) > 1, "s", ""),
          " ..."
        )
      )

      if (length(gene_ids) > 100) {
        message("This may take a few minutes ...")
      }
      if (isTRUE(isexon)) {
        # calculate length and GC only for exonic regions
        # set feature, start, and stop equal to exon values
        feature <- "ensembl_exon_id"
        start <- "exon_chrom_start"
        end <- "exon_chrom_end"
        typ <- "gene_exon"
      } else {
        # calculate length and GC only for entire transcript
        # set feature, start, and stop equal to transcript values
        feature <- "ensembl_transcript_id"
        start <- "transcript_start"
        end <- "transcript_end"
        typ <- "transcript_exon_intron"
      }

      # pull feature level coordinates by either exon or transcript coordinates
      #"hgnc_symbol", "gene_biotype",
      attrs <- c(filters, feature, "chromosome_name", start, end)
      coords <- biomaRt::getBM(
        filters = filters,
        attributes = attrs,
        values = gene_ids,
        mart = ensembl,
        useCache = FALSE
      )

      #Pull the feature Sequences from Biomart
      seqs <- biomaRt::getSequence(
        id = coords[,feature],
        type = feature,
        seqType = typ,
        mart = ensembl,
        useCache = FALSE
      )

      # Insert the sequences into coords
      row.names(seqs) <- seqs[,feature]
      coords$sequence <- NA
      coords$sequence <- seqs[ coords[,feature], ][,typ]

      coords[,start] <- as.numeric(coords[,start])
      coords[,end] <- as.numeric(coords[,end])

      # Biomart gene column name
      gene_value <- names(coords)[1]

      # Calculate GC content and Gene Length
      cl <- snow::makeCluster(cores, outfile="log.log")
      len_gc <- data.frame()

      for (chr in names(table(coords$chromosome_name))) {

        df <- coords[coords$chromosome_name == chr,]
        genes <- as.list(gene_ids[gene_ids %in%
                            coords[coords$chromosome_name == chr,]$ensembl_gene_id
        ])

        message(paste0("Starting Chromosome: ", chr))
        foo <- as.data.frame(do.call(rbind, parallel::parLapply(
          cl = cl,
          genes,
          biomart_stats,
          column_id = gene_value,
          df = df,
          start = start,
          end = end
        )))
        len_gc <- as.data.frame(rbind(len_gc,foo))
      }

      colnames(len_gc)[1] <- gene_value
      snow::stopCluster(cl)
      rm(cl)

      # Pull gene info
      attrs <- c(filters, "hgnc_symbol", "gene_biotype", "chromosome_name")
      gene_info <- biomaRt::getBM(
        filters = filters,
        attributes = attrs,
        values = gene_ids,
        mart = ensembl,
        useCache = FALSE
      )

      # Join Data:
      biomart_results <- gene_info %>%
        dplyr::full_join(len_gc, by = gene_value)
      biomart_results <- biomart_results[, c(gene_value, 'hgnc_symbol',
                                             'percentage_gene_gc_content', 'gene_biotype',
                                             'chromosome_name', 'gene_length'
      )
      ]
      biomart_results$percentage_gene_gc_content <- as.numeric(biomart_results$percentage_gene_gc_content)
      biomart_results$gene_length <- as.numeric(biomart_results$gene_length)
    } else {
      # use custom specified GTF and FASTA from synapse

      # Load GTF and FASTA from synapse
      synapser::synLogin()
      gtf <- synapser::synGet(gtfID, version = gtfVersion)
      genome <- synapser::synGet(fastaID, version = fastaVersion)

      biom <- GenomicTools.fileHandler::importGTF(file=gtf$path,
                                                  level="gene",
                                                  features=c("gene_id",
                                                             "gene_name",
                                                             "gene_type"))
      biom <- biom[ , c("V1", "gene_id", "gene_name", "gene_type")]

      # set feature level
      if (isTRUE(isexon)) {
        # calculate length and GC only for exonic regions
        lvl <- "exon"
      }else{
        # calculate length and GC only for entire transcript
        lvl <- "transcript"
      }
      # pull features from gtf file
      feature_gtf <- GenomicTools.fileHandler::importGTF(
        file=gtf$path,
        level=lvl
      )

      #Calc GC content
      dna <- Biostrings::readDNAStringSet(genome$path)

      # assemble chromosomes chromosome
      chroms <- list()
      for(contig in names(table(biom$V1))){
        chroms[[contig]] <- dna[grep(contig,names(dna),value=T)][1]
      }

      # Translation for fasta chrom name alterations
      region_names <- rep(NA, length(names(chroms)))
      names(region_names) <- names(chroms)
      for(nam in names(region_names)){
        region_names[nam] <- names(chroms[[nam]])
      }

      # Gene named vector with of chromosome designations
      #gene_chr <- biom$V1
      #names(gene_chr) <- biom$gene_id
      gene_chr <- biom$gene_id
      #genome <- chroms[gene_chr[i]]

      cl <- snow::makeCluster(cores, outfile = "log.log")
      stats <- matrix(NA,0,3)
      for (element in names(table(biom$V1))) {
        calcs <- do.call(rbind, parallel::parLapply(
          cl = cl,
          as.list(biom[biom$V1 == element, ]$gene_id),
          gtf_stats,
          data = as.data.frame(feature_gtf)[feature_gtf$V1 == element,],
          genome = chroms[[element]]
        ))
        stats <- rbind(stats,calcs)
      }
      snow::stopCluster(cl)
      rm(cl)

      stats <- as.data.frame(stats)
      colnames(stats) <- c(filters, 'percentage_gene_gc_content',
                           'gene_length')
      colnames(biom) <- c( 'chromosome_name', filters, 'hgnc_symbol',
                           'gene_biotype')
      biomart_results <- biom %>%
        dplyr::full_join(stats, by = 'ensembl_gene_id')

      #Remove PAR_Y
      biomart_results <- biomart_results[
        !(grepl('_PAR_Y', biomart_results$ensembl_gene_id)),
      ]

      # Clean chr from chromosomes
      biomart_results$chromosome_name <- gsub(
        'chr', '', biomart_results$chromosome_name
      )

      # Clean ENSGs
      biomart_results <- as.data.frame(biomart_results)
      biomart_results[,filters] <- do.call(
        rbind,
        strsplit(biomart_results[,filters], '[.]'))[,1]
      biomart_results <- biomart_results[ ,c(filters, 'hgnc_symbol',
                                             'percentage_gene_gc_content',
                                             'gene_biotype',
                                             'chromosome_name', 'gene_length')]
      biomart_results <- as.data.frame(biomart_results)
      biomart_results$percentage_gene_gc_content <- as.numeric(biomart_results$percentage_gene_gc_content)
      biomart_results$gene_length <- as.integer(biomart_results$gene_length)
    }
    # Finish cleaning bioMart Object
    # Duplicate Ensembl Ids are collapsed into a single entry
    biomart_results <- collapse_duplicate_hgnc_symbol(biomart_results)

    # Biomart IDs as rownames
    biomart_results <- tibble::column_to_rownames(
      biomart_results, var = filters
    )
  }else{

    # Download biomart object from syndID specified in config.yml
    biomart_results <- get_data(synid, version)

    # Biomart IDs as rownames
    biomart_results <- tibble::column_to_rownames(
      biomart_results, var = filters
    )

    # Gene metadata required for count CQN
    required_variables <- c("gene_length", "percentage_gene_gc_content")
    if (!all(required_variables %in% colnames(biomart_results))) {
      vars <- glue::glue_collapse(
        setdiff(required_variables, colnames(biomart_results)),
        sep = ", ",
        last = " and "
      )
      message(glue::glue("Warning: {vars} missing from biomart object.
                         This information is required for Conditional
                         Quantile Normalization"))
    }
  }
  return(biomart_results)
}
#'Check for consistent identifiers
#'
#'Subsequent steps assume unique identifiers are present within both the counts
#'matrix and the metadata. This function checks whether samples in the counts
#'matrix match samples in the metadata, and vice versa.versa.
#'@inheritParams get_biomart
#'@param md A data frame with sample identifiers as row names and relevant experimental covariates.
#'@export
check_mismatch <- function(md, count_df) {
  # Isolate the SampleID's from the metadata and counts matrix
  md_sampleid <- rownames(md)
  cnt_sampleid <- colnames(count_df)

  # Ensures each sample in the metadata has a corresponding column in the counts matrix and vice versa
  if ((length(setdiff(md_sampleid, cnt_sampleid)) > 0) || (length(setdiff(cnt_sampleid, md_sampleid)) > 0)){
    stop("All sample identifiers in the counts matrix and metadata file must match")
  }
}

#'Duplicate HGNC
#'
#'Count normalization requires Ensembl Ids to be unique. In rare cases, there are more
#'than one HGNC symbol per gene Id. This function collapses the duplicate entries into
#'a single entry by appending the HGNC symbols in a comma separated list.
#'@param biomart_results Output of \code{"sageseqr::get_biomart()"}. Gene Ids are
#'stored as rownames.
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
#' @param conditions Optional. Conditions to bin gene counts that correspond to
#' variables in `md`.
#' @param clean_metadata A data frame with sample identifiers as rownames and variables as
#' factors or numeric as determined by \code{"sageseqr::clean_covariates()"}.
#' @importFrom magrittr %>%
#' @export
filter_genes <- function(clean_metadata, count_df,
                         cpm_threshold, conditions_threshold,
                         conditions = NULL) {

  if (class(conditions) == "list") {
    conditions <- unique(conditions[[1]])
  } else {
    conditions <- unique(conditions)
  }

  # Check for extraneous rows
  count_df <- parse_counts(count_df)

  if (!is.null(conditions)) {
    if (!any(conditions %in% colnames(clean_metadata))) {
      stop("Conditions are missing from the metadata.")
    }

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
  } else {
    genes_to_analyze <- simple_filter(
      count_df,
      cpm_threshold,
      conditions_threshold
    )
  }

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
#' Dropped Genes Finder
#'
#' Finds gene IDs in the filtered counts which were removed in cqn normalization
#' most likely as a function if no GC or length estimates being present in the
#' user specified biomart object.
#' @param filtered_counts The target containing counts after low gene
#' expression has been removed. Defaults to target name constrained by
#'  \code{"targets::tar_make()"}.
#' @param cqn_counts A counts data frame normalized by CQN.
#'
#' @export
dropped_genes <- function (filtered_counts, cqn_counts) {
  missing <- row.names(filtered_counts)[
    !(row.names(filtered_counts) %in% row.names(cqn_counts))
  ]
  if (length(missing) == 0) {
    dropped <- NULL
  }else{
    dropped <- missing
  }
  return(dropped)
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
#' @param exclude_variables Vector of variables to exclude from testing.
#' @inheritParams coerce_factors
#' @export
build_formula <- function(md, primary_variable, model_variables = NULL,
                          exclude_variables = NULL) {
  if (!(all(purrr::map_lgl(md, function(x) inherits(x, c("numeric", "factor")))))) {
    stop("Use sageseqr::clean_covariates() to coerce variables into factor and numeric types.")
  }

  if (!is.null(model_variables)) {
    md <- dplyr::select(md, dplyr::all_of(c(model_variables, primary_variable)))
  }

  if (!is.null(exclude_variables)) {
    if (exclude_variables %in% colnames(md)) {
      stop("exclude_variables and model_variables are the same.")
    }
    md <- dplyr::select(md, -dplyr::all_of(exclude_variables))
  }

  # Update metadata to reflect variable subset

  # Variables of factor or numeric class are required
  col_type <- dplyr::select(md, -dplyr::all_of(primary_variable)) %>%
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
                formula_base_model = formula(
                  glue::glue("~ {interaction_term}")
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
#' @param p_value_threshold Numeric. P-values are adjusted by Benjamini and
#' Hochberg (BH) false discovery rate (FDR). Significant genes are those with an
#' adjusted p-value greater than this threshold.
#' @param fold_change_threshold Numeric. Significant genes are those with a
#' fold-change greater than this threshold.
#' @param cores An integer of cores to specify in the parallel backend (eg. 4).
#' @inheritParams cqn
#' @inheritParams coerce_factors
#' @inheritParams build_formula
#' @inheritParams cqn
#' @export
#' @return A named list with \code{"variancePartition::voomWithDreamWeights()"}
#'  normalized counts, contrasts from \code{"variancePartition::getContrasts()"},
#'  linear mixed model fits from \code{"variancePartition::dream()"}, differential
#'  expression results from \code{"limma::topTable()"} and gene feature-specific
#'  metadata, the response variable, the model formula fit to compute differential
#'  expression results.
#'  The list names are: \code{"list(voom_object, contrasts_to_plot, fits, differential_expression,
#'  primary_variable, formula)"}
differential_expression <- function(filtered_counts, cqn_counts, md,
                                    primary_variable, biomart_results,
                                    p_value_threshold, fold_change_threshold,
                                    model_variables = NULL,
                                    exclude_variables = NULL,
                                    cores = NULL) {
  # force order of samples in metadata to match order of samples in counts.
  # Required by variancePartition
  if(is.null(cores)){
    cores = parallel::detectCores()-1
  }
  md <- md[match(colnames(filtered_counts),rownames(md)),]

  metadata_input <- build_formula(md, primary_variable, model_variables)
  gene_expression <- edgeR::DGEList(filtered_counts)
  gene_expression <- edgeR::calcNormFactors(gene_expression)
  voom_gene_expression <- variancePartition::voomWithDreamWeights(counts = gene_expression,
                                                                  formula = metadata_input$formula,
                                                                  data = metadata_input$metadata,
                                                                  BPPARAM = BiocParallel::SnowParam(cores)
                                                                 )
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
                                           L = contrasts,
                                           BPPARAM = BiocParallel::SnowParam(cores)
                                           )

  de <- lapply(names(contrasts), function(i, fit){
    genes <- limma::topTable(fit, coef = i, number = Inf, sort.by = "logFC")
    genes <- tibble::rownames_to_column(genes, var = "ensembl_gene_id")
  }, fit_contrasts)

  names(de) <- names(contrasts)

  de <- data.table::rbindlist(de, idcol = "Comparison") %>%
    dplyr::mutate(Comparison = gsub(metadata_input$primary_variable, "", .data$Comparison),
                  Direction = .data$logFC/abs(.data$logFC),
                  Direction = ifelse(.data$Direction == -1,"down", .data$Direction),
                  Direction = ifelse(.data$Direction == 1, "up", .data$Direction),
                  Direction = ifelse(.data$`adj.P.Val` > p_value_threshold | abs(.data$logFC) < log2(fold_change_threshold),
                                     "none",
                                     .data$Direction)
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
#' @param conditions A named list of conditions to test as `primary_variable`
#' in \code{"sagseqr::differential_expression()"}.
#' @param dropped a vector of gene names to drop from filtered counts, as they
#' were not cqn normalized
#' @inheritParams differential_expression
#' @inheritParams build_formula
#' @export
wrap_de <- function(conditions, filtered_counts, cqn_counts, md, dropped,
                    biomart_results, p_value_threshold, fold_change_threshold,
                    model_variables = names(md), cores = NULL) {

  if (!is.null(dropped)) {
    filtered_counts <- filtered_counts[
      !(row.names(filtered_counts) %in% dropped),
    ]
  }
  purrr::map(
    conditions,
    function(x) differential_expression(
      filtered_counts,
      cqn_counts,
      md,
      primary_variable = x,
      biomart_results,
      p_value_threshold,
      fold_change_threshold,
      model_variables,
      cores = cores
      )
    )
}
#' Stepwise Regression
#'
#' This function performs multivariate forward stepwise regression evaluated by multivariate Bayesian Information
#' Critera (BIC) by wrapping \code{"mvIC::mvForwardStepwise()"}.
#'
#' @inheritParams differential_expression
#' @inheritParams build_formula
#' @param skip Defaults to NULL. If TRUE, this step will be skipped in the
#' targets plan.
#' @return Table with BIC criteria for exclusion or inclusion of variables in
#' the model, linear (mixed) model formula and vector of variables to include.
#' @export
stepwise_regression <- function(md, primary_variable, cqn_counts,
                                model_variables = names(md),
                                skip = NULL) {
  # skip stepwise generation if skip = TRUE
  if(isTRUE(skip)) {
    return("Skipping stepwise regression model generation...")
  } else {
  metadata_input <- build_formula(md, primary_variable, model_variables)
  model <- mvIC::mvForwardStepwise(exprObj = cqn_counts$E,
                                   baseFormula = metadata_input$formula_base_model,
                                   data = metadata_input$metadata,
                                   variables = array(metadata_input$variables)
  )

  to_visualize <- model$trace %>%
    dplyr::select(.data$iter, .data$variable, .data$delta,
                  .data$score, .data$isBest, .data$isAdded, .data$m) %>%
    dplyr::rename(iteration = .data$iter,
                  `best (tested against baseline)` = .data$isBest,
                  `added to model` = .data$isAdded,
                  `(effective) number of parameters estimated` = .data$m
                  )

  # return vector of variables that map to the metadata columns. This requires
  # some cleaning of extraneous characters. Do not include primary variable(s)

  to_include <- model$trace$variable[model$trace$isAdded == "yes" & model$trace$iter != 0]
  to_include <- gsub("\\)|\\(1\\||scale\\(", "", to_include)

  output <- list(
    to_visualize = to_visualize,
    formula = model$formula,
    variables_in_model = to_include
  )

  return(output)
  }
}
#' Summarize Biotypes
#'
#' Computes the fraction of genes of a particular biotype. The number of genes
#' must be above 100 to be summarized.
#'
#' @inheritParams collapse_duplicate_hgnc_symbol
#' @inheritParams cqn
#' @export
summarize_biotypes <- function(filtered_counts, biomart_results) {
  biomart_results[rownames(filtered_counts),] %>%
    dplyr::group_by(.data$gene_biotype) %>%
    dplyr::summarise(fraction = dplyr::n()) %>%
    dplyr::filter(.data$fraction > 100) %>%
    dplyr::mutate(fraction = .data$fraction/dim(filtered_counts)[1])
}
#' Prepare output
#'
#' Store data in temporary files to prepare for Synapse upload.
#' @param target The object to be stored.
#' @param data_name An identifier to embed in the file name.
#' @param rowname Optional. If applicable, the name of the variable to store
#' rownames.
#' @export
prepare_results <- function(target, data_name, rowname = NULL) {
  if (!is.null(rowname)) {
    target <- tibble::rownames_to_column(target, rowname)
  }

  # the file name will contain the primary name of the target
  tmp <- fs::file_temp(data_name, ext = ".tsv")

  utils::write.table(
    target,
    file = tmp,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
    )
  tmp
}
#' Store output to Synapse
#'
#' Store transformed data, markdown and html report to a Folder in Synapse.
#'
#' @param parent_id A Synapse Id that corresponds to a project or
#' folder to store output.
#' @param cqn_counts The target containing Conditional Quantile Normalized
#'  (CQN) counts. Defaults to target name constrained by
#'  \code{"targets::tar_make"}.
#' @param clean_md The target containing the metadata data frame.
#' Defaults to target name constrained by \code{"targets::tar_make"}.
#' @param filtered_counts The target containing counts after low gene
#' expression has been removed. Defaults to target name constrained by
#'  \code{"targets::tar_make()"}.
#' @param biomart_results The target containing gene annotations from
#' biomart. Defaults to target name constrained by
#' \code{"targets::tar_make"}.
#' @param de_results The target containing differential expression gene
#' lists. Defaults to target name constrained in the _targets.R file.
#' @param residualized_counts The target containing counts adjusted for batch
#' effects. Defaults to target name constrained in the _targets.R file.
#' @param report The target containing the rendered html document. Defaults
#' to target name constrained in the _targets.R file.
#' @param rownames A list of variables to store rownames ordered by `metadata`,
#' `filtered_counts`, `biomart_results`, `cqn_counts`. If not applicable,
#' set as NULL.
#' @param syn_names A list of human-readable names for the Synapse entities
#' ordered
#' by `metadata`, `filtered_counts`, `biomart_results`, `cqn_counts`.
#' @param inputs A character vector of Synapse Ids to create provenance between
#' output files and input files.
#' @param activity_provenance A phrase to describe the data transformation for
#' provenance.
#' @param data_names A list of identifiers to embed in the file name ordered
#' by `clean_md`, `filtered_counts`, `biomart_results`, `cqn_counts`,
#' `de_results`.
#' @param config_file Optional. Path to configuration file.
#' @param report_name Name of output markdown file.
#' @export
store_results <- function(clean_md = clean_md,
                          filtered_counts = filtered_counts,
                          biomart_results = biomart_results,
                          cqn_counts = cqn_counts$E,
                          de_results = de,
                          residualized_counts = residualized_counts,
                          report = report,
                          syn_names, data_names,
                          parent_id, inputs, activity_provenance,
                          rownames = NULL, config_file = NULL,
                          report_name = NULL) {

  # include sageseqr package version in Synapse provenance
  ver <- utils::packageVersion("sageseqr")
  description <- glue::glue(
    "analyzed with sageseqr {ver}"
  )

  # nest targets targets in a list. Every time a new target is to-be stored, it
  # must be added as an argument to this function and then added to this list.
  targets <- list(
    clean_md,
    filtered_counts,
    biomart_results,
    as.data.frame(cqn_counts)
  )

  # parse differential expression gene list
  de_results <- purrr::map(de_results, function(x) x$differential_expression)

  # parse residualized counts
  parse_residual_matrix <- purrr::map(residualized_counts, function(x) x$output)

  # append differential expression data frames already nested in list
  targets <- append(targets, parse_residual_matrix)
  targets <- append(targets, de_results)

  # append null rownames for differential expression object
  rownames <- append(rownames, rep(list(NULL), length(residualized_counts) + length(de_results)))

  mash <- list(
    target = targets,
    rowname = rownames,
    data_name = data_names
  )

  file_location <- purrr::pmap(
    mash,
    prepare_results
    )

  mash <- list(
    parent = parent_id,
    syn_names = syn_names,
    paths = file_location
    )

  file_to_upload <- purrr::pmap(
    mash,
    function(paths, parent, syn_names) synapser::File(
      path = paths,
      parent = parent,
      name = syn_names
      )
  )

  # nest input Synapse Ids to apply the same provenance to all files
  inputs <- list(inputs)

  mash <- list(
    files = file_to_upload,
    inputs = inputs,
    activity_provenance = activity_provenance,
    activity_description = description
    )

  # login to Synapse
  synapser::synLogin()

  for_provenance <- purrr::pmap(
    mash,
    function(
      files, inputs, activity_provenance, activity_description
      ) synapser::synStore(
        obj = files,
        used = inputs,
        activityName = activity_provenance,
        activityDescription = activity_description,
        forceVersion = FALSE
      )
  )

  if (!is.null(config_file)) {
    file <- synapser::File(
      path = config_file,
      parent = parent_id,
      name = "Configuration file"
    )

    config_provenance <- synapser::synStore(
      obj = file,
      activityName = activity_provenance
      )

    for_provenance <- append(for_provenance, config_provenance)
  }

  if (!is.null(report_name)) {
    path <- glue::glue("{getwd()}/{report_name}.html")

    used_ids <- unlist(
      purrr::map(
        for_provenance,
        ~.x$get("id")
        )
    )

    file <- synapser::File(
      path = path,
      parent = parent_id
    )

    markdown_provenance <- synapser::synStore(
      obj = file,
      used = used_ids,
      activityName = activity_provenance,
      activityDescription = description
    )
  }
  message(glue::glue("Files uploaded to {parent_id}"))
}
#' Provenance helper
#'
#' Collapse Synapse Ids and version.
#' @param metadata_id Synapse ID to clean metadata file with sample identifiers
#' in a column and variables of interest as column names. There cannot be any
#' missing values.
#' @param counts_id Synapse ID to counts data frame with identifiers to the
#' metadata as column names and gene ids in a column.
#' @param metadata_version Optionally, include Synapse file version number. If
#' omitted, current version will be downloaded.
#' @param counts_version Optionally, include Synapse file version number.
#' @param biomart_id Synapse ID to biomart object.
#' @param biomart_version Optionally, include Synapse file version number.
#' @export
provenance_helper <- function(metadata_id,  counts_id, metadata_version = NULL,
                              counts_version = NULL, biomart_id = NULL,
                              biomart_version = NULL) {

  if (is.null(metadata_version)) {
    ids <- metadata_id
  } else {
    ids <- glue::glue("{metadata_id}.{metadata_version}")
  }

  if (is.null(counts_version)) {
    ids <- c(ids, counts_id)
  } else {
    ids <- c(ids, glue::glue("{counts_id}.{counts_version}"))
  }

  if (is.null(biomart_version) & is.null(biomart_id)) {
    ids
  } else if (is.null(biomart_version) & !is.null(biomart_id)) {
    ids <- c(ids, biomart_id)
  } else {
    ids <- c(ids, glue::glue("{biomart_id}.{biomart_version}"))
  }
  ids
}
#' Initialize differential expression analysis workflow
#'
#' The `sageseqr` package provides a `targets` workflow to string together the
#' data processing and computational steps of RNA-seq differential expression
#' analysis. This funciton copies a markdown document and _targets.R file to your
#' working directory. The _targets.R file in your working directory is required
#' for the workflow to run.
#'
#' @export
start_de <- function() {
  # copy sageseqr-report.Rmd markdown to working directory
  if (!file.exists("sageseqr-report.Rmd")) {
    fs::file_copy(system.file("sageseqr-report.Rmd", package = "sageseqr"),
                  new_path = getwd())
  }

  # copy _targets.R file to working directory
  if (!file.exists("_targets.R")) {
    fs::file_copy(system.file("_targets.R", package = "sageseqr"),
                  new_path = getwd())
  }
}
#' Compute residualized counts matrix
#'
#' Residuals of the best fit linear regression model are computed for each
#'  observation. Batch effects are adjusted for in the returned counts matrix
#'  while preserving the effect of the predictor variable.
#'
#' Counts are normalized prior to linear modeling to compute residuals. A
#' precision weight is assigned to each gene feature to estimate the
#' mean-variance relationship. Counts normalized by conditional quantile
#'  normalization (CQN) are used in place of log2 normalized counts.
#'
#' @inheritParams cqn
#' @inheritParams differential_expression
#' @inheritParams filter_genes
#' @param dropped a vector of gene names to drop from filtered counts, as they
#' were not cqn normalized
#' @export
compute_residuals <- function(clean_metadata, filtered_counts, dropped,
                              cqn_counts = cqn_counts$E, primary_variable,
                              model_variables = NULL, cores = NULL)  {

  if (!is.null(dropped)) {
    filtered_counts <- filtered_counts[
      !(row.names(filtered_counts) %in% dropped),
    ]
  }

  # set the number of cores if not specified in the config
  if(is.null(cores)){
    cores = parallel::detectCores()-1
  }
  # force order of samples in metadata to match order of samples in counts.
  # Required by variancePartition
  clean_metadata <- clean_metadata[match(colnames(filtered_counts),rownames(clean_metadata)),]

  metadata_input <- build_formula(clean_metadata, primary_variable, model_variables)
  # Estimate voom weights with DREAM
  gene_expression <- edgeR::DGEList(filtered_counts)
  gene_expression <- edgeR::calcNormFactors(gene_expression)
  voom <- variancePartition::voomWithDreamWeights(
    counts = gene_expression,
    formula = metadata_input$formula,
    data = metadata_input$metadata,
    BPPARAM = BiocParallel::SnowParam(cores)
  )

  # fit linear model using weights and best model
  voom$E <- cqn_counts
  adjusted_fit <- variancePartition::dream(
    exprObj = voom,
    formula = metadata_input$formula,
    data = metadata_input$metadata,
    computeResiduals = TRUE,
    BPPARAM = BiocParallel::SnowParam(cores)
  )

  # compute residual matrix
  residual_gene_expression <- stats::residuals(adjusted_fit)

  # calculate weighted residuals and add back signal from predictor
  variables_to_add_back <- grep(
    metadata_input$primary_variable,
    colnames(adjusted_fit$design),
    value = TRUE
    )
  output <- residual_gene_expression +
    adjusted_fit$coefficients[
      ,variables_to_add_back
      ] %*% t(
        adjusted_fit$design[
          ,variables_to_add_back]
        )

  # save gene features
  output <- tibble::rownames_to_column(as.data.frame(output), var = "feature")

  return(
    list(
      output = output,
      signal = variables_to_add_back,
      adjusted_fit = adjusted_fit,
      formula = metadata_input$formula
      )
  )
}
