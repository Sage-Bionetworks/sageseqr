<!-- badges: start -->
  [![R build status](https://github.com/kelshmo/sageseqr/workflows/R-CMD-check/badge.svg)](https://github.com/kelshmo/sageseqr/actions)
<!-- badges: end -->
# Installation 
`remotes::install_github("Sage-Bionetworks/sageseqr")`

# RNA-seq normalization workflow in R

The `sageseqr` package integrates the [`drake` R package](https://github.com/ropensci/drake/), the [`config` package for R](https://cran.r-project.org/web/packages/config/vignettes/introduction.html), and [Synapse](https://www.synapse.org/). `drake` tracks dependency relationships in the workflow and only updates data when it has changed. A `config` file allows inputs and parameters to be explicitly defined in one location. Synapse is a data repository that allows sensitive data to be [stored and shared responsibly](https://docs.synapse.org/articles/article_index.html#governance). 

The workflow takes RNA-seq gene counts and sample metadata as inputs, normalizes counts by [CQN](https://bioconductor.org/packages/release/bioc/html/cqn.html), removes outliers based on a user-defined threshold, empirically selects meaningful covariates and returns differential expression analysis results.

# Access to Data 

Anyone can create a [Synapse account](https://docs.synapse.org/articles/getting_started.html) and access public data in a variety of disciplines: [Alzheimer's Disease Knowledge portal](https://adknowledgeportal.synapse.org/), [CommonMind Consoritum](https://www.synapse.org/#!Synapse:syn2759792/wiki/69613).   
