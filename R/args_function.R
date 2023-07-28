#' @title Arguments for EEM
#' @description
#' Extract arguments necessary to run EEM from interaction matrix
#' @param interaction_matrix interaction signs matrix. If model is GLV or Gompertz it can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second.     #if model is Baker, the interaction_matrix has to be a list of two lists, the first list contains matrices defining lower and upper bounds of alphas, the second list contains matrices defining lower and upper bounds of betas
#' @param bounds_growth_rate vector of 2 elements containing lower and upper bounds for growth rates
#' @param upper_interaction_strength coefficient representing interaction strength. Default:1
#' @param model model representing species interactions, default "GLV" (Generalized Lotka Volterra). options include "Baker", "Adams" and "customized"
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
  #define global arguments
  n_species <- n_species_function(interaction_matrix, model)#number of species
  args <- EEMtoolbox::get_nonzero_parameters(interaction_matrix, n_species, model)
  args$model <- model#considered model
  args$n_species <- n_species#number of species

  # define uniform prior bounds
  #lower lim
  args$lower <- c(rep(1,args$n_species)*bounds_growth_rate[1],
                  args$lower_interaction_bound*upper_interaction_strength)

  #upper lim
  args$upper <- c(rep(1,args$n_species)*bounds_growth_rate[2],
                  args$upper_interaction_bound*upper_interaction_strength)
  return(args)
}
