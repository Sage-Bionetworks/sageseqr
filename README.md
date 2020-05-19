<!-- badges: start -->
  [![R build status](https://github.com/kelshmo/sageseqr/workflows/R-CMD-check/badge.svg)](https://github.com/kelshmo/sageseqr/actions)
<!-- badges: end -->
# sage-seqr
RNASeq normalization pipeline with R.

This repository provides R markdowns for performing covariate adjustments and differential expression analysis of the RNA-seq data from AMP-AD and a Docker environment, thus all package dependencies are loaded. 

The results is an environment ready for analysis. 

# Install Docker 

1. [Find your supported platform and download Docker.](https://docs.docker.com/v17.12/install/#supported-platforms)

# Run Container

This Docker container is built on [`rocker/tidyverse`](https://hub.docker.com/r/rocker/tidyverse/) producing a debian stable work environment. 

1. To get container for the first time: `docker pull kelshmo/covarr`
2. Suggested configuration for running container: 
```
docker run -d -p 8787:8787 -e PASSWORD=<insert password> kelshmo/covarr #replace <password> with your password
```
- `-d` flag allows you the retain use of the terminal while the Docker image is running 
- `-p` specifies port of choice
- `-e` designates a password to login to RStudio

If you are having trouble with your Docker, this [documentation](https://ropenscilabs.github.io/r-docker-tutorial/02-Launching-Docker.html) (or [this documentation](https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image)) may be useful.

3. Now, RStudio has launched invisibly. Open a browser and enter your ip address followed by :8787. You should be greeted by the RStudio welcome screen.

username: rstudio
password: <password>
  
4. You can `docker stop <container name>` to pause your work and return to it later with `docker start <container name>`.

* To find your container name `docker ps -a`

# Clone the Repo

`git clone https://github.com/kelshmo/sageseqr.git` to get the `sageseqr` tool in your Docker repository.

# Access to Data 

Data access requires a [Synapse account](https://docs.synapse.org/articles/getting_started.html) and access to the [AMP-AD Knowledge portal](https://www.synapse.org/#!Synapse:syn2580853/wiki/409854).   

You may also be interested in exploring novel Alzheimer's disease targets in the [Agora portal](https://agora.ampadportal.org/genes).

# Run your own data

1. Update the default parameters to the Synapse ID where your data is stored and variables you want to test. The full list of configurable options are: 
```
counts:
  synID: Synapse ID to counts data frame with identifiers to the metadata as column names and gene ids in a column.
  version: Optionally, include Synapse file version number (e.g. 3).
  gene id: Column name that corresponds to the gene ids (e.g. feature).
metadata:
  synID: Synapse ID to cleaned metadata file with sample identifiers in a column and variables of interest as column names.
  version: Optionally, include Synapse file version number (e.g. 3).
  sample id: Column name that corresponds to the sample ids (e.g. donorid).
biomart:
  synID: This input is optional. If left blank, Ensembl will be queried with the gene ids provided in the counts. Otherwise, you may provide the Synapse ID to gene metadata from Ensembl. This must include gene length and GC content in order to implement Conditional Quantile Normalization. 
  version: Optionally, include Synapse file version number (e.g. 3).
  gene id: Column name that corresponds to the gene ids (e.g. ensembl_gene_id)
factors: List of factor variables in brackets. Variables must be present in the metadata as column names (e.g. [ "donorid", "source"]).
continuous: List of continuous variables in brackets. Variables must be present in the metadata as column names (e.g[ "rin", "rin2"]).
```
2. `initialize.R` sets the config file and loads the package dependencies, functions and drake plan. Set `Sys.setenv(R_CONFIG_ACTIVE = "default")` in `initialize.R`. 
