#' @title Arguments for EEM
#' @description
#' Extract number of species from interaction matrix
#' @param interaction_matrix interaction signs matrix, can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second
#' @param model model representing species interactions, default "GLV" (Generalized Lotka Volterra). options include "Baker", "Adams" and "customized"
#' @return A list of arguments defining the problem.
#' @examples
#' n_species_function(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' @return n_species: number of species
#' @export

n_species_function <- function(interaction_matrix,model="GLV"){

  if (model=="GLV"|model=="Gompertz"){
    #if model is GLV or Gompertz, the interaction_matrix can be input as a single
    #matrix of interactions or as a list of matrices defining lower and upper bounds
    #for interaction terms lower first and upper second
    if (class(interaction_matrix)[1]=="matrix"){
      return(ncol(interaction_matrix)) #number of species in ecosystem network
    } else {#
      return(ncol(interaction_matrix[[1]])) #number of species in ecosystem network
    }
  }
  else if (model == "Baker"){
    #if model is Baker, the interaction_matrix can be:
    # a list of two matrices (alphas and betas), the default lower bound is 0 for alphas,
    # and the default upper bound is 0 for betas
    # OR a list of two lists
    # the first list contains matrices defining lower and upper bounds of alphas
    # the second list contains matrices defining lower and upper bounds of betas
    if (class(interaction_matrix[[1]])[1]=="matrix"){
      return(ncol(interaction_matrix[[1]])) #number of species in ecosystem network
    } else {#
      return(ncol(interaction_matrix[[1]][[1]]))
    }
  }
}
