Sys.setenv(R_CONFIG_ACTIVE = "default")

# TO DO: how to setup cache for end user
# make_cache

source("R/functions.R")
source("R/plan.R")

# Login to Synapse. Make a Synapse account and use synaper to login: https://r-docs.synapse.org/articles/manageSynapseCredentials.html
synapser::synLogin()

# Check for changes to the config file
config <- config::get()
source("check_targets.R")

# Run the analysis
drake::make(execute_plan(config = config))

# Visualize the results
drake::vis_drake_graph(execute_plan(config = config),
                       targets_only = TRUE)
