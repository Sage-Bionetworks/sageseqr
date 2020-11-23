<!-- badges: start -->
  [![R build status](https://github.com/kelshmo/sageseqr/workflows/R-CMD-check/badge.svg)](https://github.com/kelshmo/sageseqr/actions)
<!-- badges: end -->
# Installation 
`remotes::install_github("Sage-Bionetworks/sageseqr")`

# RNA-seq normalization workflow in R

The `sageseqr` package integrates the [`drake` R package](https://github.com/ropensci/drake/), the [`config` package for R](https://cran.r-project.org/web/packages/config/vignettes/introduction.html), and [Synapse](https://www.synapse.org/). `drake` tracks dependency relationships in the workflow and only updates data when it has changed. A `config` file allows inputs and parameters to be explicitly defined in one location. Synapse is a data repository that allows sensitive data to be [stored and shared responsibly](https://docs.synapse.org/articles/article_index.html#governance). 

The workflow takes RNA-seq gene counts and sample metadata as inputs, normalizes counts by conditional quantile normalization [(CQN)](https://bioconductor.org/packages/release/bioc/html/cqn.html), removes outliers based on a user-defined threshold, empirically selects meaningful covariates and returns differential expression analysis results. The data is also visualized in several ways to help you understand meaningful trends. The visualizations include a heatmap identifying highly correlated covariates, a sample-specific x and y marker gene check, boxplots visualizing the distribution of continuous variables and a principal component analysis (PCA) to visualize sample distribution.

# The Targets

The series of steps that make up the workflow are called targets. The target objects are stored in a cache and can either be read or loaded into your environment with the `drake` functions `readd` or `loadd`. Source code for each target can be visualized by setting `show_source = TRUE` with `loadd` and `readd`. The targets are called by the `sageseqr` `rnaseq_plan` function and are:

Raw data: 
- `import_metadata`
- `import_counts`
- `biomart_results` - the complete list of genes with biomaRt annotations.

Exploratory data visualizations:
- `gene_coexpression` - the distribution of correlated gene counts.
- `boxplots` - the distribution of continuous variables.
- `sex_plot` - the distribution of samples by x and y marker genes.
- `correlation_plot` - the correlation of covariates.
- `significant_covariates_plot` - the correlation of covariates to gene 
                                  expression.
- `outliers` - the clustering of samples by PCA.

Transformed or normalized data:
- `clean_md` - metadata with factor and numeric types.
- `filtered_counts` - counts matrix with low gene expression removed.
- `biotypes` - gene proportions summarized by biotype.
- `cqn_counts` - CQN normalized counts. 
- `model` - Model selected by multivariate forward stepwise regression 
            (evaluated by Bayesian Information Criteria (BIC)).

# Access to Data 

Anyone can create a [Synapse account](https://docs.synapse.org/articles/getting_started.html) and access public data in a variety of disciplines: [Alzheimer's Disease Knowledge portal](https://adknowledgeportal.synapse.org/), [CommonMind Consoritum](https://www.synapse.org/#!Synapse:syn2759792/wiki/69613).   
