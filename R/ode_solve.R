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
#'
ode_solve <- function(parameters, t_window, IC){
  t = seq(from = t_window[1],to = t_window[2], by = 0.01)
  y = matrix(0,length(t),length(IC))
  y[1,] = IC
  for (i in 2:length(t)){
    y[i,] = pmax(y[i-1,] + (t[2]-t[1])*fun(t[i-1],y[i-1,]), 0)
  }
  return(list(t,y))
}
