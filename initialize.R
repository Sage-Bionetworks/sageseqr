Sys.setenv(R_CONFIG_ACTIVE = "default")

# TO DO: how to setup cache for end user
# make_cache

source("R/packages.R")
source("R/functions.R")

synLogin()

source("R/plan.R")

drake::make(plan)
