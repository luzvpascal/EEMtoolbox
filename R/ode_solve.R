#' @title Solve ODE
#' @description
#' Solves ODEs
#' @param initial_condition vector of initial species abundances
#' @param interaction_matrix_value interaction matrix parametrized. If model is GLV or Gompertz, input a single matrix. If model is Baker, input a list of 2 matrices: the first element is the alphas and the second the betas.
#' @param growth_rate parametrized growth rate vector.
#' @param t_window time window to solve ODE
#' @param time_step_len length of each time step, default = 0.01
#' @param model model representing species interactions. Default "GLV" (Generalized Lokta Voltera). options include "Baker", "Gompertz" and "customized"
#' @param derivative derivative function. Default \link[EEMtoolbox]{derivative_func}
#' @examples
#' library(EEMtoolbox)
#' output <- EEM(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' interaction_matrix_value <- output[[1]]$interaction_matrix
#' growth_rate <- output[[1]]$growthrates
#' initial_condition <- rep(10, 8) # 8 species, initial abundance is 10 for all
#' ode_solve(initial_condition, interaction_matrix_value, growth_rate, t_window=c(0,1))
#' @return list: time_steps: vector of time steps of the simulation
#' abundances: matrix of abundances, each row corresponds to the values of abundances at a given time step
#' @export

ode_solve <- function(initial_condition,
                      interaction_matrix_value,
                      growth_rate,
                      t_window,
                      time_step_len=0.01,
                      model = "GLV",
                      derivative=EEMtoolbox::derivative_func){
  time_steps <- seq(from = t_window[1],to = t_window[2],by = time_step_len)#vector of time steps
  y <-  matrix(0,nrow=length(time_steps),ncol=length(initial_condition)) #matrix of abudnances

  y[1,] = initial_condition #setting initial condition

  for (i in 2:length(time_steps)){
    y[i,] = pmax(y[i-1,] + time_step_len*derivative(interaction_matrix_value,
                                                    growth_rate,
                                                    y[i-1,],
                                                    model),
                 0)
  }
  return(list(time_steps=time_steps,
              abundances=y))
}
