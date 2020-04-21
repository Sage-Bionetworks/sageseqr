Sys.setenv(R_CONFIG_ACTIVE = "mayo")

# TO DO: how to setup cache for end user
# make_cache

source("R/packages.R")
source("R/functions.R")

synLogin()

source("R/plan.R")

# Run the analysis
drake::make(plan)

# Visualize the results
drake::vis_drake_graph(plan, targets_only = TRUE)
