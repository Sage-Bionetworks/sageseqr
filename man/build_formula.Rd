% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{build_formula}
\alias{build_formula}
\title{Formula Builder}
\usage{
build_formula(
  md,
  primary_variable,
  model_variables = NULL,
  exclude_variables = NULL
)
}
\arguments{
\item{md}{A data frame with sample identifiers in a column and relevant experimental covariates.}

\item{primary_variable}{Vector of variables that will be collapsed into a single
fixed effect interaction term.}

\item{model_variables}{Optional. Vector of variables to include in the linear (mixed) model.
If not supplied, the model will include all variables in \code{md}.}

\item{exclude_variables}{Vector of variables to exclude from testing.}
}
\description{
If a linear mixed model is used, all categorical variables must be
modeled as a random effect. This function identifies the variable class
to scale numeric variables and model categorical variables as a random
effect. Additionally, variables are scaled to account for
multiple variables that might have an order of magnitude difference.
}
