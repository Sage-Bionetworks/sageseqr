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
  tmp <- cor_estimate
  tmp[cor_pval > p_value_cutoff] <- 0
  col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                 "#4393C3", "#2166AC", "#053061")))
  p <- corrplot::corrplot(tmp, col = col2(200), tl.col = "#000000")

  return(
    list(estimate = cor_estimate, pval = cor_pval, plot = p))
}
#' Compute the association statistic for factor variables
#'
#' Returns the Cramer's V estimate and the P-value from the Pearson
#' chi-square score between two variables.
#'
#' @inheritParams coerce_factors
#' @export
#'
association_statistic_for_factors <- function(factors, md, na.action = "remove") {
  if (na.action == "remove")
    md <- na.omit(md[, factors])

  fac1 = as.factor(md[,1])
  fac2 = as.factor(md[,2])

  stats = vcd::assocstats(stats::xtabs(~ fac1 + fac2))

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
#'
#' @inheritParams coerce_factors
#' @export
association_statistics_for_both <- function(variables = names(md), md,
                                            na.action = "remove", alpha = 0.05) {
  if (na.action == "remove")
    md = na.omit(md[, variables])

  stats = psych::ICC(md[, variables], alpha = alpha)

  pval = summary(aov(md[, variables[1]] ~md[, variables[2]]))[[1]][["Pr(>F)"]][1]

  return(
    c(
      Estimate = stats$results['Single_raters_absolute','ICC'],
      Pval = pval
      )
    )
}
