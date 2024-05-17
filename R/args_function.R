#' @title Arguments for EEM
#' @description
#' Extract arguments necessary to run EEM from interaction matrix
#' @param interaction_matrix interaction signs matrix. If model is GLV or Gompertz it can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second.     #if model is Baker, the interaction_matrix has to be a list of two lists, the first list contains matrices defining lower and upper bounds of alphas, the second list contains matrices defining lower and upper bounds of betas
#' @param upper_bounds_growth_rate upper bound of growth rates. Input can be one number (same upper bound for all species) or a vector of growth rates upper bounds for each species. Default 5
#' @param lower_bounds_growth_rate lower bound of growth rates. Input can be one number (same lower bound for all species) or a vector of growth rates lower bounds for each species. Default 0
#' @param upper_interaction_strength coefficient representing interaction strength. Default:1
#' @param model model representing species interactions, default "GLV" (Generalized Lotka Volterra). options include "Baker", "Adams" and "customized"
#' @return A list of arguments defining the problem.
#' @examples
#' library(EEMtoolbox)
#' args_function(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' @return list:
#' model: considered model
#' n_species: number of species
#' skip_parameters: vector of 0-1 indicating which parameters to skip
#' num_params: number of non-zero parameters
#' lower: vector of length num_params representing lower bound of each parameter
#' upper: vector of length num_params representing upper bound of each parameter
#' @export
args_function <- function(interaction_matrix,
                          upper_bounds_growth_rate=5,
                          lower_bounds_growth_rate=0,
                          upper_interaction_strength=1,
                          model="GLV"){
  #define global arguments
  n_species <- EEMtoolbox::n_species_function(interaction_matrix, model)#number of species
  args <- EEMtoolbox::get_nonzero_parameters(interaction_matrix, n_species, model)
  args$model <- model#considered model
  args$n_species <- n_species#number of species

  # define uniform prior bounds
  #lower lim
  if (length(lower_bounds_growth_rate)==1){
    args$lower <- c(lower_bounds_growth_rate*rep(1,n_species),
                    args$lower_interaction_bound)
  } else {
    args$lower <- c(lower_bounds_growth_rate,
                    args$lower_interaction_bound)
  }


  #upper lim
  if (length(upper_bounds_growth_rate)==1){
    args$upper <- c(upper_bounds_growth_rate*rep(1,n_species),
                    args$upper_interaction_bound)
  } else {
    args$upper <- c(upper_bounds_growth_rate,
                    args$upper_interaction_bound)
  }


  return(args)
}
