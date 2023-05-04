#' @title Solve ODE
#' @description
#' Solves ODEs
#' @param interaction_matrix parametrized interaction matrix
#' @param growth_rate parametrized growth rate vector
#' @param t_window time window to solve ODE
#' @param initial_condition vector of initial species abundances
#' @param model model representing species interactions. Default "GLV" (Generalized Lokta Voltera). options include "Baker", "Gompertz" and "customized"
#' @examples
#' library(EEMtoolbox)
#' EEM(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' @return list: part_vals: ensemble of parameters, marginal distributions
#' @export

ode_solve <- function(interaction_matrix, growth_rate, t_window, initial_condition, model = "GLV"){
  time_steps <- seq(from = t_window[1],to = t_window[2],by = 0.01)#vector of time steps
  y <-  matrix(0,length(t),length(initial_condition)) #matrix of abudnances

  y[1,] = initial_condition #setting initial condition

  for (i in 2:length(time_steps)){
    y[i,] = pmax(y[i-1,] + (time_steps[2]-time_steps[1])*EEMtoolbox::derivative_func(interaction_matrix,
                                                                   growth_rate,
                                                                   y[i-1,],
                                                                   model),
                 0)
  }
  return(list(time_steps=time_steps,
              y=y))
}
