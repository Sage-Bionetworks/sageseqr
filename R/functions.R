#'
#'Detect file type and download data
#'
#'This function takes synIds and version number to download a series of rectangular files into a nested tibble. The file contents are
#'stored nested in a variable filecontents.
#'
#'@param inputs rectangular data where id = synId and version = file version number are required variables. Tidyverse argument required.
#'@export
#'@return A tibble.
#'@examples
#'\dontrun{
#'
#'inputs <- tibble(object = c("count", metadata"), id = c("syn5678", "syn2345"), version = c("1", "8"))
#'
#'input_data <- get_data(inputs)
#'}
get_data <- function(inputs){
  inputs %>%
    mutate(thefile = purrr::map2(id, 
                                 version, 
                                 function(x,y) as_tibble(data.table::fread(synapser::synGet(x, version = y)$path))))
}
#'
#'Create covariate matrix from tidy metadata data frame. 
#'
#'This function takes a tidy format.
clean_covariates <- function(df, factors, continuous){
  
}