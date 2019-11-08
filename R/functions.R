#'
#'Detect file type and download data
#'
#'This function takes synIds and version number to download any rectangular file type from Synapse.
#'
#'@param synID A character vector with a Synapse Id.
#'@param version Optional. A numeric vector with the Synapse file version number.
#'@export
#'@return A tibble.
#'@examples
#'\dontrun{
#'file <- get_data(synID = "syn1234", version = 7)
#'
#'}
get_data <- function(synID, version = NULL){
  df <- as_tibble(data.table::fread(synapser::synGet(synID, version = as.numeric(version))$path))
  df
}
#'
#'Create covariate matrix from tidy metadata data frame.
#'
#'This function takes a tidy format. Coerces objects to correct type.
#'
#'@param md
#'@param factors A vector of factor variables
#'@param continuous A vector of continuous variables
#' Remove primary variable for now - A vector of primary variables
#'
#'@export
#'@return A data frame with coerced variables.
#'@examples
clean_covariates <- function(md, factors, continuous){
  if(!missing(factors) & !missing(continuous)) {
    md[,factors] <- lapply(md[,factors], factor)
    md[,continuous] <- lapply(md[,continuous], as.numeric)
    md
  } else {
    message("Factor and continuous variables are required.")
  }

}
#'
#'Explore metadata by variable.
#'
#'This function produces boxplots from the variables provided.
#'
#'@param md
#'@param vars A vector of variables to visualize
#'
#'@export
#'@return
#'@examples
#' TO DO: test this
boxplot_vars <- function(md, include_vars, x_var){
  md %>%
    dplyr::select(., c(include_vars, x_var)) %>%
    gather(key, value, -x_var) %>%
    ggplot(aes(x = x_var, y = value)) +
    geom_boxplot() +
    theme(legend.position = 'top') +
    facet_grid(key ~ !! x_var, scales = "free")
}
#' Filter genes
#'
#'@param md A data frame with sample identifiers in a column
#'@param count_matrix
#'@param sample_variable
#'
#'

