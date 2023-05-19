#' @title Arguments for EEM
#' @description
#' Extract arguments necessary to run EEM from interaction matrix
#'
#' @param interaction_matrix interaction signs matrix, can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second
#' @param bounds_growth_rate vector of 2 elements containing lower and upper bounds for growth rates
#' @param upper_interaction_strength coefficient representing interaction strength. Default:1
#' @param model model representing species interactions, default "GLV" (Generalized Lokta Voltera). options include "Baker", "Adams" and "customized"
#' @return A list of arguments defining the problem.
#' @examples
#' library(EEMtoolbox)
#' args_function(dingo_matrix, c(-5,5)) #automatically loads an example of interaction matrix as dingo_matrix
#' @return list:
#' model: considered model
#' n_species: number of species
#' skip_parameters: vector of 0-1 indicating which parameters to skip
#' num_params: number of non-zero parameters
#' lower: vector of length num_params representing lower bound of each parameter
#' upper: vector of length num_params representing upper bound of each parameter
#' @export
args_function <- function(interaction_matrix,
                          bounds_growth_rate,
                          upper_interaction_strength=1,
                          model="GLV"){
  args <- list()

  args$model <- model

  if (class(interaction_matrix)[1]=="matrix"){
    args$n_species <- ncol(interaction_matrix) #number of species in ecosystem network
  } else {#
    args$n_species <- ncol(interaction_matrix[[1]]) #number of species in ecosystem network
  }

  non_zero_params <- EEMtoolbox::get_nonzero_parameters(interaction_matrix)
  #define global arguments
  args$skip_parameters <- non_zero_params$skip_parameters

  #number of parameters excluding 0 parameters
  args$num_params <- args$n_species+args$n_species^2 - sum( args$skip_parameters)

  # define uniform prior bounds
  #lower lim
  args$lower <- c(rep(1,args$n_species)*bounds_growth_rate[1],
                  non_zero_params$lower_interaction_bound*upper_interaction_strength)

  args$upper <- c(rep(1,args$n_species)*bounds_growth_rate[2],
                  non_zero_params$upper_interaction_bound*upper_interaction_strength) #upper lim


  return(args)
}
