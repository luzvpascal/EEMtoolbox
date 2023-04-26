#' @title Arguments for EEM
#' @description
#' Extract arguments necessary to run EEM from interaction matrix
#'
#' @param interaction_matrix interaction signs matrix, can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second
#' @param bounds_growth_rate vector of 2 elements containing lower and upper bounds for growth rates
#' @param upper_interaction_strength vector of 2 elements containing lower and upper bounds for growth rates
#' @param model model representing species interactions, default "GLV" (Generalized Lokta Voltera). options include "Baker", "Adams" and "customized"
#' @return A list of arguments defining the problem.
#' @export
args_function <- function(interaction_matrix,
                          bounds_growth_rate,
                          upper_interaction_strength=1,
                          model="GLV"){
  args <- list()

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

  if (model=="Baker"){

  }

  return(args)
}
