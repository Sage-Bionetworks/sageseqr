% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{compute_residuals}
\alias{compute_residuals}
\title{Compute residualized counts matrix}
\usage{
compute_residuals(
  clean_metadata,
  filtered_counts,
  dropped,
  cqn_counts = cqn_counts$E,
  primary_variable,
  random_effect = NULL,
  model_variables = NULL,
  is_num = NULL,
  num_var = NULL,
  cores = NULL
)
}
\arguments{
\item{clean_metadata}{A data frame with sample identifiers as rownames and variables as
factors or numeric as determined by \code{"sageseqr::clean_covariates()"}.}

\item{filtered_counts}{A counts data frame with genes removed that have low expression.}

\item{dropped}{a vector of gene names to drop from filtered counts, as they
were not cqn normalized}

\item{cqn_counts}{A counts data frame normalized by CQN.}

\item{primary_variable}{Vector of variables that will be collapsed into a single
fixed effect interaction term.}

\item{random_effect}{A vector of variables to consider as random effects instead
of fixed effects.}

\item{model_variables}{Optional. Vector of variables to include in the linear (mixed) model.
If not supplied, the model will include all variables in \code{md}.}

\item{is_num}{Is there a numerical covariate to use as an interaction with the primary variable(s). default= NULL}

\item{num_var}{A numerical metadata column to use in an inaction with the primary variable(s). default= NULL}

\item{cores}{An integer of cores to specify in the parallel backend (eg. 4).}
}
\description{
Residuals of the best fit linear regression model are computed for each
 observation. Batch effects are adjusted for in the returned counts matrix
 while preserving the effect of the predictor variable.
}
\details{
Counts are normalized prior to linear modeling to compute residuals. A
precision weight is assigned to each gene feature to estimate the
mean-variance relationship. Counts normalized by conditional quantile
 normalization (CQN) are used in place of log2 normalized counts.
}
