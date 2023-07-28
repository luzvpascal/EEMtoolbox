#' @title Derivatives of models
#' @description
#' Derivatives functions to solve ODE systems
#' @param Time a float of current studied time
#' @param State vector of current species abundance
#' @param Pars a list like object containing: model: string of model representing species interactions. Options include "GLV" (Generalized Lotka Volterra), "Baker", "Gompertz". interaction_matrix_value interaction matrix parametrized. If model is GLV or Gompertz, input a single matrix. If model is Baker, input a list of 2 matrices: the first element is the alphas and the second the betas. growth_rate: vector of growth rates
#' @examples
#' library(EEMtoolbox)
#' output <- EEM(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' current_abundance <- rep(10, 8) # 8 species, initial abundance is 10 for all
#' Pars = list(model="GLV",
#'   interaction_matrix_value=output[[1]]$interaction_matrix,
#'   growth_rate <- output[[1]]$growthrates)
#' derivative_func(Time=0, State=current_abundance, Pars=Pars)
#' @return list: part_vals: ensemble of parameters, marginal distributions
#' @export
derivative_func <- function(Time,
                            State,
                            Pars){

  with(as.list(c(State, Pars)), {
    if (model=="GLV"){
      dState <- State*growth_rate + (interaction_matrix_value%*%State)*State
      return(list(dState))
    }
    if (model=="Baker"){
        A <- interaction_matrix_value[[1]]
        B <- interaction_matrix_value[[2]]

        P <- diag(A)
        M <- A-diag(P)
        dState <- State*growth_rate*(1-exp(-M%*%State-P))+(B%*%State)*State
        return(list(dState))
    }
    if (model=="Gompertz"){
        dState <- State*growth_rate +
          (interaction_matrix_value%*%log(State))*State
        return(list(dState))
    }
  })
}
