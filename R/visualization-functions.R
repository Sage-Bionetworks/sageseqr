#' Compute principal component analysis (PCA) and plot correlations
#'
#' Identify principal components (PCs) of normalized gene counts and correlate
#' these PCs to interesting covariates. This function wraps
#' `correlate_and_plot()` to visualize, with a heatmap, the relationship between
#' PCs and covariates that meet a false discovery rate (FDR) threshold and
#' return a list of significant covariates.
#'
#' @param percent_p_value_cutoff The p-value threshold in percent.
#' @param scaled Defaults to TRUE. Variables scaled to have unit
#' variance before the analysis takes place.
#' @param normalized_counts A counts data frame normalized by CQN, TMM, or
#' another preferred method, with genes as rownames.
#' @param plot_covariates_vs_pca Defaults to TRUE. If false, no plot is
#' returned.
#' @inheritParams filter_genes
#' @inheritParams correlate_and_plot
#' @return A list.
#' \itemize{
#'   \item significant_covariates - A vector of covariates where correlation
#'   p-value meets the FDR threshold.
#'   \item pc_results - A customized heatmap of significant covariates and PCs
#'   correlated.
#'   \item effects_significant_vars - A vector correlation values.
#' }
#' @export
run_pca_and_plot_correlations <- function(normalized_counts, clean_metadata,
                                      scaled = TRUE,
                                      percent_p_value_cutoff = 1.0,
                                      correlation_type = "pearson",
                                      plot_covariates_vs_pca = TRUE,
                                      maximum_fdr = 0.1) {
  md <- clean_metadata

  pca_results <- run_pca(normalized_counts = normalized_counts,
                   scaled = scaled,
                   percent_p_value_cutoff = percent_p_value_cutoff)

  sample_pc_values <- pca_results$sample_pc_values
  pve <- pca_results$pve

  npca <- ncol(sample_pc_values)

  # append percent variance explained to the PC matrix
  colnames(sample_pc_values) <- paste(
    colnames(sample_pc_values), " (", sprintf("%.2f", pve[1:npca]), "%)", sep=""
    )

  # subset metadata to remove variables or covariates that have missingness
  complete <- md[, which(
    apply(md,
          2,
          function(dat) all(!is.na(dat))
          ))]
  exclude_missing_data <- setdiff(colnames(md), colnames(complete))

  loop <- FALSE
  if (plot_covariates_vs_pca) {
    loop <- unique(c(loop, TRUE))
    }

  for (plot_all in loop) {
    p <- correlate_and_plot(
      sample_pc_values,
      md,
      correlation_type,
      weights = pve[1:dim(sample_pc_values)[2]],
      plot_all,
      exclude_missing_data,
      maximum_fdr = 0.1
      )
  }
  return(
    list(
      significant_covariates = p$significant_covariates,
      pc_results = p$plot,
      effects_significant_vars = p$effects_significant_vars
      )
    )
}
#'
#' Principal Component Analysis (PCA) of samples
#'
#' The counts matrix is transposed to compute principal components of each
#' sample. The data is centered, scaled and rotated by default. The principal
#' components (PCs) where the proportion of variance explained (PVE) meets
#' the \code{"percent_p_value_cutoff"} are returned. The default percent cutoff
#' is 1%.
#'
#' @inheritParams run_pca_and_plot_correlations
#' @export
#' @return A list with significant PCs rotated and PVE.
#' \itemize{
#'   \item sample_pc_values - A matrix of PCs with samples as rows.
#'   \item pve - A numeric vector of PVE greater than
#'   \code{"percent_p_value_cutoff"}.
#' }
run_pca <- function(normalized_counts, scaled = TRUE,
                    percent_p_value_cutoff = 1.0) {

  # estimate variance in data by PC, variables are zero centered and rotated
  # variables are returned:
  pca_results <- stats::prcomp(
    t(normalized_counts), center = TRUE, scale.= scaled, retx = TRUE
    )

  # examine how much variance is explained by PCs and consider those with
  # proportion of variance explained (PVE) >= (percent_p_value_cutoff %):
  pc_variance <- pca_results$sdev^2
  pve <- 100 * (pc_variance / sum(pc_variance))
  npca <- max(1, length(which(pve >= percent_p_value_cutoff)))

  sample_pc_values <- pca_results$x[, 1:npca, drop = FALSE]

  list(sample_pc_values = sample_pc_values,
       pve = pve)
}
#' Apply false discovery rate (FDR) threshold to correlation matrix.
#'
#' Adjusts a given a set of p-values using FDR and returns a vector of logical
#' values indicating whether the value is less than or equal to the threshold
#' FDR.
#'
#' @inheritParams correlate_and_plot
#' @param correlation_values A data frame with p-values in column `pvalue`.
#' @param mask_indices A logical vector to indicate whether a covariate has
#' missing values and should be skipped when computing FDR.
#' @examples
#'\dontrun{
#' correlation_values <- data.frame(
#'   compare = "PC1 (23.16%)",
#'   covariates = "batch",
#'   pvalue = 0.56519736,
#'   r = "0.05729811"
#'   )
#'}
#'@export
threshold_fdr <- function(correlation_values, mask_indices = NULL,
                          maximum_fdr = 0.1) {
  if (is.null(mask_indices)) {mask_indices = 1:nrow(correlation_values)}

  fdr = rep(1.0, nrow(correlation_values))
  fdr[mask_indices] = stats::p.adjust(
    correlation_values$pvalue[mask_indices],
    method = "fdr"
    )

  return (fdr <= maximum_fdr)
}
#' Calculate correlation between covariates and principal components (PCs)
#'
#' Identify covariates that significantly, as defined by a false discovery rate
#' (FDR) threshold, correlate with PCs computed from gene counts. This function
#' wraps `plot_pcs_with_covariates()` to return a heatmap to visualize this
#' relationship.
#'
#' @param principal_components Center, scaled and rotated principal components
#' matrix.
#' @inheritParams filter_genes
#' @param correlation_type Allowed values are "pearson", "spearman" and
#' "kendall". See \code{"psych::corr.test(method)"}.
#' @param weights A numeric vector of the proportion of variance explained
#' (PVE). Defaults to NULL.
#' @param to_plot Logical indicating whether correlation values should be
#' plotted. Defaults to FALSE.
#' @param exclude_missing_data A vector with column names to exclude from
#' `clean_metadata`.
#' @param maximum_fdr Maximum allowable false discovery rate (FDR). Defaults to
#' 0.1.
#' @return A list.
#' \itemize{
#'   \item plot - A customized heatmap of significant covariates and PCs
#'   correlated.
#'   \item significant_covariates - A vector of covariates where correlation
#'   p-value meets the FDR threshold.
#'   \item effects_significant_vars - A vector correlation values.
#' }
#' @export
correlate_and_plot <- function(principal_components, clean_metadata,
                                   correlation_type, weights = NULL,
                                   to_plot = FALSE,
                                   exclude_missing_data = NULL,
                                   maximum_fdr = 0.1) {
  # Identify factors and continuous variables
  md <- clean_metadata
  factors <- names(md)[sapply(md, is.factor)]
  continuous <- setdiff(names(md), factors)

  # Convert factor variables to numeric vector
  md[, factors] <- lapply(
    md[, factors],
    function(x) {
      x <- as.numeric(unclass(x))
    }
  )

  # Compute intraclass correlation (ICC)
  # between principal_components and factor covariates
  if (length(factors) > 0) {
    factor_cor <- apply(
      expand.grid(colnames(principal_components), factors),
      1,
      association_statistics_for_both,
      cbind(
        principal_components,
        md[rownames(principal_components), factors]
        ),
      p_value_cutoff = maximum_fdr
      )
    factor_estimate <- matrix(
      factor_cor['Estimate',],
      nrow = length(colnames(principal_components)),
      ncol = length(factors))
    factor_pval <- matrix(
      factor_cor['Pval',],
      nrow = length(colnames(principal_components)),
      ncol = length(factors))

    rownames(factor_estimate) <- colnames(principal_components)
    colnames(factor_estimate) <- factors

    rownames(factor_pval) <- colnames(principal_components)
    colnames(factor_pval) <- factors
  } else {
    factor_estimate <- NULL
    factor_pval <- NULL
  }

  # Calculate correlation between principal_components and continuous covariates
  if (length(continuous) > 0) {
    continuous_cor <- psych::corr.test(principal_components,
                          md[, continuous],
                          use = "pairwise.complete.obs",
                          method = correlation_type,
                          adjust = "none")
    continuous_estimate <- continuous_cor$r
    continuous_pval <- continuous_cor$p

    rownames(continuous_estimate) <- colnames(principal_components)
    colnames(continuous_estimate) <- continuous

    rownames(continuous_pval) <- colnames(principal_components)
    colnames(continuous_pval) <- continuous
  } else {
    continuous_estimate <- NULL
    continuous_pval <- NULL
  }

  estimate <- cbind(factor_estimate, continuous_estimate)
  pval <- cbind(factor_pval, continuous_pval)

  effects_significant_vars <- estimate
  # assign estimates (correlation and ICC to PCs) with a p-value greater than
  # the threshold FDR to 0
  effects_significant_vars[pval > maximum_fdr] <- 0
  effects_significant_vars <- colSums(
    abs(effects_significant_vars)*replicate(
      dim(effects_significant_vars)[2],
      weights/sum(weights)
      )
    )
  effects_significant_vars <- effects_significant_vars[
    order(
      abs(effects_significant_vars),
      decreasing = TRUE
      )
    ]

  cor_mat = reshape2::melt(pval, varnames = c("compare", "covariates"))
  colnames(cor_mat)[colnames(cor_mat) == "value"] = "pvalue"

  cor_mat$compare <- factor(cor_mat$compare, levels = rownames(pval))
  cor_mat$covariates <- factor(cor_mat$covariates, levels = colnames(pval))

  cor_mat$r <- reshape2::melt(estimate)$value

  calcFDRrows <- rep(TRUE, nrow(cor_mat))
  mark_missing <- NULL

  # if exclude_missing_data is not null, identify which rows to skip when
  # FDR computing
  if (!is.null(exclude_missing_data)) {
    calcFDRrows <- !(cor_mat$covariates %in% exclude_missing_data)
    mark_missing <- intersect(colnames(md), exclude_missing_data)
  }

  # determine entries that pass the threshold of "significance"
  significant_correlations <- threshold_fdr(
    cor_mat,
    mask_indices = calcFDRrows,
    maximum_fdr = 0.1
    )
  sorted <- sort(
    unique(cor_mat$covariates[significant_correlations])
    )

  #psg = possibly significant correlations
  psg <- threshold_fdr(cor_mat)
  # Mark only those incomplete covariates that would be significant
  # in the context of all covariates:
  psg <- psg & !calcFDRrows

  plot_rows <- 1:nrow(cor_mat)
  if (!to_plot) {
    # Plot all correlations for:
    # 1) Covariates with at least one significant correlation
    # 2) Excluded covariates
    plot_rows <- (cor_mat$covariates %in% sorted) | !calcFDRrows
  }
  plot <- stats::na.omit(cor_mat[plot_rows, ])

  for (mark in c("significant_correlations", "psg")) {
    use_mark <- get(mark)[plot_rows]
    if (length(which(use_mark)) > 0) {
      plot[, mark] <- use_mark[
        setdiff(
          1:length(use_mark),
          as.numeric(attr(plot, "na.action"))
          )
        ]
    }
  }
  # plot if to_plot = TRUE
  if (!plyr::empty(plot)) {
    plot <- plot_pcs_with_covariates(
      plot,
      glue::glue("*FDR <= {maximum_fdr}")
      )
  } else {
    plot <- NULL
  }

  return(
    list(
      plot = plot,
      significant_covariates = as.character(sorted),
      effects_significant_vars = effects_significant_vars
      )
    )
}
#' Explore relationship between principal components (PCs) and covariates
#'
#'
#' @param correlation_values A data frame that includes the p-value of the
#' correlation of a covariate to a PC as `pvalue`. Additionally, the data frame
#' should include a logical column `significant_correlations` to indicate
#' whether the FDR threshold is met for that variable.
#' @param note_threshold A customized string to note method used to adjust
#' p-values for multiple comparisons. Defaults to NULL.
#' @examples
#'\dontrun{
#' correlation_values <- data.frame(
#'   compare = "PC1 (23.16%)",
#'   covariates = "batch",
#'   pvalue = 0.56519736,
#'   r = "0.05729811",
#'   significant_correlations = FALSE
#'   )
#' note_threshold <- "FDR <= 0.1"
#'}
#'@export
plot_pcs_with_covariates <- function(correlation_values,
                               note_threshold = NULL) {

  # Reverse the Y-axis
  correlation_values$compare <- factor(correlation_values$compare,
                                       levels = rev(
                                         levels(
                                           correlation_values$compare
                                           )
                                         )
                                       )

  # Set threshold from -1 to 1
  plot <- correlation_values
  b <- c(-1,0,1)

  # Add * notation for p-values that meet the FDR
  plot$significant_correlations[plot$significant_correlations == TRUE] <- "*"
  plot$significant_correlations[plot$significant_correlations == FALSE] <- ""

  # Plot correlation
  p <- ggplot2::ggplot(
    plot,
    ggplot2::aes(
      .data$covariates,
      .data$compare)
    ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data$r),
      colour = "white"
      )
  p <- p + ggplot2::geom_text(
    ggplot2::aes(
      label = .data$significant_correlations
      )
    )
  p <- p + ggplot2::scale_fill_gradientn(
    limits = c(-1,1),
    colours = rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
    "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
    "#4393C3", "#2166AC", "#053061")),
    breaks = b,
    labels = format(b)) +
    ggplot2::labs(
      x = "",
      y = "",
      caption = note_threshold) +
    sagethemes::theme_sage() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 12)
    )


  return(p)
}
#' Visualize relationship between variables
#'
#' Output a correlation matrix to visualize the relationship between
#' continuous and factor variables. For factor variables, the Cramer's
#' V estimate and the P-value is computed from the Pearson chi-square
#' score between two variables.
#'
#' @param p_value_cutoff Set the p-value threshold.
#' @inheritParams filter_genes
#' @export
get_association_statistics <- function(clean_metadata, p_value_cutoff = 0.05) {

  # Identify factors and continuous variables
  md <- clean_metadata
  factors <- names(md)[sapply(md, is.factor)]
  continuous <- setdiff(names(md), factors)

  # Convert factor variables to numeric vector
  md[, factors] <- lapply(
    md[, factors],
    function(x) {
      x <- as.numeric(unclass(x))
      }
    )

  # Find association between factor variables pair-wise when there are more than 1
  # factor variable
  if (length(factors) != 0 & length(factors) != 1) {
    factor_cor <- apply(
      expand.grid(factors, factors),
      1,
      association_statistic_for_factors,
      md[, factors]
      )

    factor_estimate <- matrix(
      factor_cor['Estimate',],
      nrow = length(factors),
      ncol = length(factors)
      )
    factor_pval <- matrix(
      factor_cor['Pval',],
      nrow = length(factors),
      ncol = length(factors)
      )

    colnames(factor_estimate) <- factors
    rownames(factor_estimate) <- factors

    colnames(factor_pval) <- factors
    rownames(factor_pval) <- factors

  } else if (length(factors) == 1) {
    factor_estimate <- as.data.frame(1)
    colnames(factor_estimate) <- factors
    rownames(factor_estimate) <- factors

    factor_pval <- as.data.frame(0)
    colnames(factor_pval) <- factors
    rownames(factor_pval) <- factors

  } else {
    factor_estimate <- NULL
    factor_pval <- NULL
  }

  # Find correlation between continuous variables
  if (length(continuous) > 1) {
    tmp <- apply(
      md[, continuous, drop = F],
      2,
      as.numeric
      )
    continuous_cor = psych::corr.test(tmp, use = 'pairwise.complete.obs')
    continuous_estimate = continuous_cor$r
    continuous_pval = continuous_cor$p

  } else if (length(continuous) == 1) {
    continuous_estimate <- as.data.frame(1)
    colnames(continuous_estimate) <- continuous
    rownames(continuous_estimate) <- continuous

    continuous_pval <- as.data.frame(0)
    colnames(continuous_pval) <- continuous
    rownames(continuous_pval) <- continuous

  } else {
    continuous_estimate <- NULL
    continuous_pval <- NULL
  }

  # Find interclass correlation between factor and continuous variables
  if (length(factors) > 0 & length(continuous) > 0) {
    cor <- apply(
      expand.grid(factors, continuous),
      1,
      association_statistics_for_both,
      md[, c(factors, continuous)]
      )

    estimate <- matrix(
      cor['Estimate',],
      nrow = length(factors),
      ncol = length(continuous)
      )
    pval <- matrix(
      cor['Pval',],
      nrow = length(factors),
      ncol = length(continuous)
      )

    colnames(estimate) <- continuous
    rownames(estimate) <- factors

    colnames(pval) <- continuous
    rownames(pval) <- factors

  } else {
    estimate <- NULL
    pval<- NULL
  }

  # Combine all estimates that are significant
  if (length(factors) != 0 & length(continuous) == 0) {
    cor_estimate = factor_estimate
    cor_pval = factor_pval

  } else if (length(factors) == 0 & length(continuous) != 0) {
    cor_estimate = continuous_estimate
    cor_pval = continuous_pval

  } else if (length(factors) != 0 & length(continuous) != 0) {
    cor_estimate = rbind(cbind(factor_estimate, estimate),
                         cbind(t(estimate),continuous_estimate)
                         )

    cor_pval = rbind(cbind(factor_pval, pval),
                     cbind(t(pval), continuous_pval)
                     )

  } else {
    cor_estimate = NULL
    cor_pval = NULL
  }

  # plot heatmap
  p <- cor_estimate
  p[cor_pval > p_value_cutoff] <- 0

  return(
    list(estimate = cor_estimate, pval = cor_pval, plot = p)
    )
}
#' Compute the association statistic for factor variables
#'
#' Returns the Cramer's V estimate and the P-value from the Pearson
#' chi-square score between two variables.
#'
#' @inheritParams filter_genes
#' @inheritParams clean_covariates
#' @param na_action Defaults to removing rows with missing values.
#' @export
#'
association_statistic_for_factors <- function(factors, clean_metadata, na_action = "remove") {
  md <- clean_metadata
  if (na_action == "remove")
    md <- stats::na.omit(md[, factors])

  fac1 <- as.factor(md[,1])
  fac2 <- as.factor(md[,2])

  stats <- vcd::assocstats(x = stats::xtabs(~ fac1 + fac2))

  return(
    c(
      Estimate = stats$cramer,
      Pval = stats$chisq_tests['Pearson','P(> X^2)']
      )
    )
}
#' Find intraclass correlation (ICC) between factor and continuous covariates
#'
#' Fit a one-way ANOVA fixed effects model.
#' @param variables Variables to for ICC computation.Defaults to the column names
#' of \code{"md"}.
#' @inheritParams filter_genes
#' @inheritParams association_statistic_for_factors
#' @inheritParams get_association_statistics
#' @export
association_statistics_for_both <- function(variables = names(clean_metadata), clean_metadata,
                                            p_value_cutoff = 0.05, na_action = "remove") {
  md <- clean_metadata
  if (na_action == "remove")
    md <- stats::na.omit(md[, variables])

  stats <- psych::ICC(md[, variables], alpha = p_value_cutoff)

  pval <- summary(stats::aov(md[, variables[1]] ~ md[, variables[2]]))[[1]][["Pr(>F)"]][1]

  return(
    c(
      Estimate = stats$results['Single_raters_absolute','ICC'],
      Pval = pval
      )
    )
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
  sex_specific <- count_df[
    grepl(
      paste0(
        rownames(
          biomart_results[
            biomart_results$hgnc_symbol %in% c("UTY", "XIST"),
            ]),
        collapse = "|"
        ),
      rownames(count_df)
      ),
    ]
  rownames(sex_specific) <- convert_geneids(sex_specific)
  sex_specific <- tibble::rownames_to_column(sex_specific, var = "geneId")
  results <- dplyr::select(
    biomart_results,
    .data$hgnc_symbol,
    .data$chromosome_name
    )
  results <- dplyr::filter(results, .data$hgnc_symbol %in% c("XIST", "UTY"))
  results <- tibble::rownames_to_column(results, var = "geneId")
  results <- dplyr::left_join(results, sex_specific)
  results <- tidyr::pivot_longer(
    results,
    -dplyr::all_of(c("geneId", "chromosome_name", "hgnc_symbol")),
    names_to = "sampleId",
    values_to = "counts(log)"
    ) %>%
    dplyr::mutate(
      `counts(log)` = log(.data$`counts(log)`),
      `counts(log)` = ifelse(
        .data$`counts(log)` == -Inf, 0, .data$`counts(log)`
        )
      )
  results <- tidyr::pivot_wider(
    dplyr::select(results, -chromosome_name, -geneId),
    names_from = "hgnc_symbol",
    values_from = "counts(log)"
    )
  results <- dplyr::left_join(results, md, "sampleId")

  p <- ggplot2::ggplot(results, ggplot2::aes(x = .data$XIST, y = .data$UTY)) +
    ggplot2::geom_point(ggplot2::aes(color = .data[[sex_var]])) +
    sagethemes::scale_color_sage_d() +
    sagethemes::theme_sage() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))

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
#' Explore samples that are outliers
#'
#' Samples that are z standard deviations (SDs) from the mean of principal
#' components (PC) 1 and 2 are identified as outliers.
#'
#' @param z Allowable number of standard deviations (SDs) from the mean.
#' Defaults to 4.
#' @inheritParams cqn
#' @inheritParams filter_genes
#' @param color Discrete variable in `clean_metadata` differentiated
#' by color.
#' @param shape Discrete variable in `clean_metadata` differentiated
#' by shape.
#' @param size Continuous variable in `clean_metadata` differentiated
#' by size.
#' @export
identify_outliers <- function(filtered_counts, clean_metadata,
                              color, shape, size, z = 4) {
  PC <- stats::prcomp(limma::voom(
    filtered_counts),
    scale. = TRUE,
    center = TRUE)

  # Plot first 2 PCs
  data <- data.frame(SampleID = rownames(PC$rotation),
                         PC1 = PC$rotation[,1],
                         PC2 = PC$rotation[,2])

  # Percentage from each PC
  eigen <- PC$sdev^2
  pc1 <- eigen[1]/sum(eigen)
  pc2 <- eigen[2]/sum(eigen)

  # Identify outliers - samples 4SDs from the mean
  outliers <- as.character(
    data$SampleID[
      c(
        which(data$PC1 < mean(data$PC1) - z*stats::sd(data$PC1)),
        which(data$PC1 > mean(data$PC1) + z*stats::sd(data$PC1))
        ),
      drop = TRUE
      ]
    )
  outliers <- c(outliers,
                as.character(
                  data$SampleID[
                    c(
                      which(data$PC2 < mean(data$PC2) - z*stats::sd(data$PC2)),
                      which(data$PC2 > mean(data$PC2) + z*stats::sd(data$PC2))
                      ),
                    drop = TRUE
                    ]
                  )
                )

  plotdata <- dplyr::left_join(
    data,
    tibble::rownames_to_column(clean_metadata, "SampleID")
    ) %>%
    dplyr::mutate(label = .data$SampleID) %>%
    dplyr::mutate(label = ifelse((.data$label %in% outliers), .data$label, NA))

  p <- ggplot2::ggplot(plotdata, ggplot2::aes(x = .data$PC1, y = .data$PC2))

  p <- p + ggplot2::geom_point(ggplot2::aes(
    color = .data[[color]],
    size = .data[[size]],
    shape = .data[[shape]]
    )
  )

  p <- p + sagethemes::scale_color_sage_d() +
    sagethemes::theme_sage() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      family = "Lato",
      size = 4,
      hjust = 0
      )

  return(
    list(
      plot = p,
      outliers = outliers
      )
  )
}
