#' @title Derivatives of models
#' @description
#' Derivatives functions to solve ODE systems
#' @param interaction_matrix_value interaction matrix parametrized. If model is GLV or Gompertz, input a single matrix. If model is Baker, input a list of 2 matrices: the first element is the alphas and the second the betas.
#' @param growth_rate vector of growth rates
#' @param current_abundance vector of current species abundance
#' @param model model representing species interactions. Default "GLV" (Generalized Lokta Voltera). Options include "Baker", "Gompertz" and "customized"
#' @examples
#' library(EEMtoolbox)
#' output <- EEM(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' interaction_matrix_value <- output[[1]]$interaction_matrix
#' growth_rate <- output[[1]]$growthrates
#' current_abundance <- rep(10, 8) # 8 species, initial abundance is 10 for all
#' derivative_func(interaction_matrix_value, growth_rate, current_abundance)
#' @return list: part_vals: ensemble of parameters, marginal distributions
#' @export
derivative_func <- function(interaction_matrix_value,
                       growth_rate,
                       current_abundance,
                       model="GLV"){

  if (model=="GLV"){
    return(current_abundance*growth_rate + (interaction_matrix_value%*%current_abundance)*current_abundance)
  }
  if (model=="Baker"){
    A <- interaction_matrix_value[[1]]
    B <- interaction_matrix_value[[2]]

    P <- diag(A)
    M <- A-diag(P)
    return(current_abundance*growth_rate*(1-exp(-M%*%current_abundance-P))
           +(B%*%current_abundance)*current_abundance)
  }
  if (model=="Gompertz"){
    return(current_abundance*growth_rate +
             (interaction_matrix_value%*%log(current_abundance))*current_abundance)
  }
  if (model == "customized"){
    print("please provide a function to compute derivatives")
    return(0)
  }

}
