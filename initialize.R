Sys.setenv(R_CONFIG_ACTIVE = "mayo")

# TO DO: how to setup cache for end user
# make_cache

synLogin()

source("R/packages.R")
source("R/functions.R")

drake::make(plan)
