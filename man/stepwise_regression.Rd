% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{stepwise_regression}
\alias{stepwise_regression}
\title{Stepwise Regression}
\usage{
stepwise_regression(
  md,
  primary_variable,
  cqn_counts,
  model_variables = names(md),
  skip = NULL,
  random_effect = NULL,
  add_model = NULL
)
}
\arguments{
\item{md}{A data frame with sample identifiers in a column and relevant experimental covariates.}

\item{primary_variable}{Vector of variables that will be collapsed into a single
fixed effect interaction term.}

\item{cqn_counts}{A counts data frame normalized by CQN.}

\item{model_variables}{Optional. Vector of variables to include in the linear (mixed) model.
If not supplied, the model will include all variables in \code{md}.}

\item{skip}{Defaults to NULL. If TRUE, this step will be skipped in the
targets plan.}

\item{random_effect}{A vector of variables to consider as random effects instead
of fixed effects.}

\item{add_model}{Optional. User Speciffied variables to add to the null model 
apriori to model generation. (Default = NULL)}
}
\value{
Table with BIC criteria for exclusion or inclusion of variables in
the model, linear (mixed) model formula and vector of variables to include.
}
\description{
This function performs multivariate forward stepwise regression evaluated by multivariate Bayesian Information
Critera (BIC) by wrapping \code{"mvIC::mvForwardStepwise()"}.
}
