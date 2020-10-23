# Function to run principal component analysis and plot correlations
runPCAandPlotCorrelations <- function(genesBySamples, samplesByCovariates, dataName, isKeyPlot=FALSE,
                                      SCALE_DATA_FOR_PCA = TRUE, MIN_PVE_PCT_PC = 1.0, CORRELATION_TYPE = "pearson",
                                      ALSO_PLOT_ALL_COVARS_VS_PCA = TRUE, MAX_NUM_LEVELS_PER_COVAR = 50) {

  title = paste(ifelse(SCALE_DATA_FOR_PCA, "S", "Un-s"), "caled ", dataName, " ", " data in PCA; PVE >= ", MIN_PVE_PCT_PC, "%; ", CORRELATION_TYPE, " correlations ", sep="")
  writeLines(paste("\nRunning PCA and calculating correlations for:\n", title, sep=""))

  pcaRes <- runPCA(genesBySamples=genesBySamples,
                   SCALE_DATA_FOR_PCA=SCALE_DATA_FOR_PCA,
                   MIN_PVE_PCT_PC=MIN_PVE_PCT_PC)

  samplePCvals <- pcaRes$samplePCvals
  pve <- pcaRes$pve

  npca <- ncol(samplePCvals)

  colnames(samplePCvals) = paste(colnames(samplePCvals), " (", sprintf("%.2f", pve[1:npca]), "%)", sep="")

  # Find covariates without any missing data
  samplesByFullCovariates = samplesByCovariates[, which(apply(samplesByCovariates, 2,
                                                              function(dat) all(!is.na(dat))))]
  EXCLUDE_VARS_FROM_FDR = setdiff(colnames(samplesByCovariates), colnames(samplesByFullCovariates))

  add_PC_res = list()
  significantCovars = c()

  LOOP_PLOT_ALL_COVARS = FALSE
  if (ALSO_PLOT_ALL_COVARS_VS_PCA) { LOOP_PLOT_ALL_COVARS = unique(c(LOOP_PLOT_ALL_COVARS, TRUE)) }

  for (PLOT_ALL_COVARS in LOOP_PLOT_ALL_COVARS) {
    corrRes = calcCompleteCorAndPlot(samplePCvals,
                                     samplesByCovariates,
                                     CORRELATION_TYPE,
                                     title,
                                     WEIGHTS = pve[1:dim(samplePCvals)[2]],
                                     PLOT_ALL_COVARS,
                                     EXCLUDE_VARS_FROM_FDR)
    add_PC_res[[length(add_PC_res)+1]] = list(plotData=corrRes$plot, isKeyPlot=(isKeyPlot && !PLOT_ALL_COVARS))
    if (!PLOT_ALL_COVARS) {
      significantCovars = corrRes$significantCovars
      Effects.significantCovars = corrRes$Effects.significantCovars
    }
  }

  return(list(significantCovars=significantCovars, PC_res=add_PC_res, Effects.significantCovars = Effects.significantCovars))
}
#'
#'# Function to run principal component analysis
runPCA <- function(genesBySamples, SCALE_DATA_FOR_PCA = TRUE, MIN_PVE_PCT_PC = 1.0) {

  # estimate variance in data by PC:
  pca.res <- prcomp(t(genesBySamples), center=TRUE, scale.=SCALE_DATA_FOR_PCA, retx=TRUE)

  # examine how much variance is explained by PCs, and consider those with PVE >= (MIN_PVE_PCT_PC %):
  pc.var <- pca.res$sdev^2
  pve <- 100 * (pc.var / sum(pc.var))
  npca <- max(1,length(which(pve >= MIN_PVE_PCT_PC)))

  samplePCvals <- pca.res$x[, 1:npca, drop=FALSE]

  list(samplePCvals=samplePCvals, pve=pve)
}
#'# Function to calculate correlation and plot
calcCompleteCorAndPlot <- function(COMPARE_data, COVAR_data, correlationType, title,
                                   WEIGHTS = NULL, PLOT_ALL_COVARS=FALSE, EXCLUDE_VARS_FROM_FDR=NULL, MAX_FDR = 0.1) {

  # require(plyr)

  # Get factor and continuous covariates
  FactorCovariates <- colnames(COVAR_data)[sapply(COVAR_data,is.factor)]
  ContCovariates <- setdiff(colnames(COVAR_data),FactorCovariates)

  # Convert factor covariates to numeric vector
  COVAR_data[,FactorCovariates] <- apply(COVAR_data[,FactorCovariates],2,
                                         function(x){x <- unclass(x)})

  # Calculate correlation between compare_data and factor covariates
  if (length(FactorCovariates) > 0){
    comb <- expand.grid(colnames(COMPARE_data),FactorCovariates)
    factCont_cor <- apply(comb,1,
                          getFactorContAssociationStatistics,
                          cbind(COMPARE_data,COVAR_data[rownames(COMPARE_data),FactorCovariates]),
                          alpha=MAX_FDR)
    factCont_cor_vals <- matrix(factCont_cor['Estimate',],
                                nrow = length(colnames(COMPARE_data)),
                                ncol = length(FactorCovariates))
    factCont_cor_p <- matrix(factCont_cor['Pval',],
                             nrow = length(colnames(COMPARE_data)),
                             ncol = length(FactorCovariates))

    rownames(factCont_cor_vals) <- colnames(COMPARE_data)
    colnames(factCont_cor_vals) <- FactorCovariates

    rownames(factCont_cor_p) <- colnames(COMPARE_data)
    colnames(factCont_cor_p) <- FactorCovariates
  } else {
    factCont_cor_vals <- NULL
    factCont_cor_p <- NULL
  }

  # Calculate correlation between compare_data and factor covariates
  if (length(ContCovariates) > 0){
    cont_cor <- corr.test(COMPARE_data,
                          COVAR_data[,ContCovariates],
                          use='pairwise.complete.obs',
                          method=correlationType,
                          adjust="none")
    cont_cor_vals <- cont_cor$r
    cont_cor_p <- cont_cor$p

    rownames(cont_cor_vals) <- colnames(COMPARE_data)
    colnames(cont_cor_vals) <- ContCovariates

    rownames(cont_cor_p) <- colnames(COMPARE_data)
    colnames(cont_cor_p) <- ContCovariates
  } else {
    cont_cor_vals <- NULL
    cont_cor_p <- NULL
  }

  all_cor_vals = cbind(factCont_cor_vals,cont_cor_vals)
  all_cor_p = cbind(factCont_cor_p,cont_cor_p)

  Effects.significantCovars = all_cor_vals
  Effects.significantCovars[all_cor_p>MAX_FDR] = 0
  Effects.significantCovars = colSums(abs(Effects.significantCovars)*replicate(dim(Effects.significantCovars)[2],WEIGHTS/sum(WEIGHTS)))
  Effects.significantCovars = Effects.significantCovars[order(abs(Effects.significantCovars),decreasing=T)]

  cor_mat = melt(all_cor_p, varnames=c("COMPARE", "COVAR"))
  colnames(cor_mat)[colnames(cor_mat) == "value"] = "pvalue"

  cor_mat$COMPARE = factor(cor_mat$COMPARE, levels=rownames(all_cor_p))
  cor_mat$COVAR = factor(cor_mat$COVAR, levels=colnames(all_cor_p))

  cor_mat$r = melt(all_cor_vals)$value

  calcFDRrows = rep(TRUE, nrow(cor_mat))
  markColumnsAsMissing = NULL
  if (!is.null(EXCLUDE_VARS_FROM_FDR)) {
    calcFDRrows = !(cor_mat$COVAR %in% EXCLUDE_VARS_FROM_FDR)
    markColumnsAsMissing = intersect(colnames(COVAR_data), EXCLUDE_VARS_FROM_FDR)
  }

  # Entries that pass the threshold of "significance":
  markSignificantCorrelations = corMatFDRthreshFunc(cor_mat, indicesMask=calcFDRrows, MAX_FDR = 0.1)
  significantCorrelatedCovars = sort(unique(cor_mat$COVAR[markSignificantCorrelations]))

  markPotentialSignificantCorrelations = corMatFDRthreshFunc(cor_mat)
  # Specially mark only those incomplete covariates that would be significant in the context of all covariates:
  markPotentialSignificantCorrelations = markPotentialSignificantCorrelations & !calcFDRrows

  plotRows = 1:nrow(cor_mat)
  if (!PLOT_ALL_COVARS) {
    # Plot all correlations for:
    # 1) Covariates with at least one significant correlation
    # 2) Excluded covariates
    plotRows = (cor_mat$COVAR %in% significantCorrelatedCovars) | !calcFDRrows
  }
  plotCor = na.omit(cor_mat[plotRows, ])

  for (markCor in c("markSignificantCorrelations", "markPotentialSignificantCorrelations")) {
    useMarkCor = get(markCor)[plotRows]
    if (length(which(useMarkCor)) > 0) {
      plotCor[, markCor] = useMarkCor[ setdiff(1:length(useMarkCor), as.numeric(attr(plotCor, "na.action"))) ]
    }
  }

  if (!plyr::empty(plotCor)){
    plot = plotCorWithCompare(plotCor, title, paste("FDR <= ", MAX_FDR, sep=""), markColumnsAsMissing)
  } else{
    plot = NULL
  }

  return(list(plot=plot, significantCovars=as.character(significantCorrelatedCovars), Effects.significantCovars = Effects.significantCovars))
}
#' Visualize relationship between variables
#'
#' Output a correlation matrix to visualize the relationship between
#' continuous and factor variables. For factor variables, the Cramer's
#' V estimate and the P-value is computed from the Pearson chi-square
#' score between two variables.
#'
#' @param p_value_cutoff Set the p-value threshold.
#' @inheritParams coerce_factors
#' @export
get_association_statistics <- function(md, p_value_cutoff = 0.05) {

  # Identify factors and continuous variables
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
    list(estimate = cor_estimate, pval = cor_pval, plot = p))
}
#' Compute the association statistic for factor variables
#'
#' Returns the Cramer's V estimate and the P-value from the Pearson
#' chi-square score between two variables.
#'
#' @inheritParams coerce_factors
#' @param na_action Defaults to removing rows with missing values.
#' @export
#'
association_statistic_for_factors <- function(factors, md, na_action = "remove") {
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
#' @inheritParams coerce_factors
#' @inheritParams association_statistic_for_factors
#' @inheritParams get_association_statistics
#' @export
association_statistics_for_both <- function(variables = names(md), md,
                                            p_value_cutoff = 0.05, na_action = "remove") {
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
