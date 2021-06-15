Sys.setenv(R_CONFIG_ACTIVE = "ihab")

# TO DO: how to setup cache for end user
# make_cache

source("R/packages.R")
source("R/functions.R")

# Login to Synapse. Make a Synapse account and use synaper to login: https://r-docs.synapse.org/articles/manageSynapseCredentials.html
synapser::synLogin()

source("R/plan.R")

# Run the analysis
drake::make(plan)


# Visualize the results
drake::vis_drake_graph( drake::drake_config(plan), targets_only = TRUE)