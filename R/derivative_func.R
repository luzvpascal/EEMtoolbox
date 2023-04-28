#' @title Derivatives of models
#' @description
#' Derivatives functions to solve ODE systems
#' @param interaction_matrix interaction signs matrix, can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second
#' @param bounds_growth_rate vector of growth rates
#' @param current_abundance vector of current species abundance
#' @param model model representing species interactions. Default "GLV" (Generalized Lokta Voltera). options include "Baker", "Gompertz" and "customized"
#' @examples
#' library(EEMtoolbox)
#' EEM(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' @return list: part_vals: ensemble of parameters, marginal distributions
#' @export
derivative_func <- function(interaction_matrix,
                       growth_rate,
                       current_abundance,
                       model="GLV"){

  if (model=="GLV"){
    return(current_abundance*growth_rate + (interaction_matrix%*%current_abundance)*current_abundance)
  }
  if (model=="Baker"){
    A <- matrix(pmax(0, interaction_matrix), ncol=ncol(interaction_matrix))
    B <- matrix(pmin(0, interaction_matrix), ncol=ncol(interaction_matrix))

    P <- diag(A)
    M <- A-P
    return(current_abundance*growth_rate*(1-exp(-M%*%current_abundance-P))
           +(B%*%current_abundance)*current_abundance)
  }
  if (model=="Gompertz"){
    return(current_abundance*growth_rate +
             (interaction_matrix%*%log(current_abundance))*current_abundance)
  }

}
