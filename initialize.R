Sys.setenv(R_CONFIG_ACTIVE = "default")

# TO DO: how to setup cache for end user
# make_cache

source("R/functions.R")
source("R/plan.R")

# Login to Synapse. Make a Synapse account and use synaper to login: https://r-docs.synapse.org/articles/manageSynapseCredentials.html
synapser::synLogin()

# Run the analysis
plan <- rnaseq_plan(metadata_id = config::get("metadata")$synID,
                    metadata_version = config::get("metadata")$version,
                    counts_id = config::get("counts")$synID,
                    counts_version = config::get("counts")$version,
                    gene_id_input = config::get("counts")$`gene id`,
                    factor_input = config::get("factors"),
                    continuous_input = config::get("continuous"),
                    gene_id = config::get("biomart")$`gene id`,
                    biomart_id = config::get("biomart")$synID,
                    biomart_version = config::get("biomart")$version,
                    filters = config::get("biomart")$filters,
                    host = config::get("biomart")$host,
                    organism = config::get("biomart")$organism
)

drake::make(plan)

# Visualize the results
drake::vis_drake_graph(plan,
                       targets_only = TRUE)
