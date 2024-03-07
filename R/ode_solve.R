#' @title Solve ODE
#' @description
#' Solves ODEs
#' @param initial_condition vector of initial species abundances
#' @param interaction_matrix_value interaction matrix parametrized. If model is GLV or Gompertz, input a single matrix. If model is Baker, input a list of 2 matrices: the first element is the alphas and the second the betas.
#' @param growth_rate parametrized growth rate vector.
#' @param t_window time window to solve ODE
#' @param time_step_len length of each time step, default = 0.01
#' @param model model representing species interactions. Default "GLV" (Generalized Lotka Volterra). options include "Baker" and "Gompertz".
#' @param derivative derivative function. Default \link[EEMtoolbox]{derivative_func}
#' @examples
#' library(EEMtoolbox)
#' output <- EEM(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' interaction_matrix_value <- output[[1]]$interaction_matrix
#' growth_rate <- output[[1]]$growthrates
#' initial_condition <- rep(10, 8) # 8 species, initial abundance is 10 for all
#' ode_solve(initial_condition, interaction_matrix_value, growth_rate, t_window=c(0,1))
#' @return matrix: first column (time) vector of time steps of the simulation, the following columns correspond to the values of abundances at a given time step
#' @export

ode_solve <- function(initial_condition,
                      interaction_matrix_value,
                      growth_rate,
                      t_window,
                      time_step_len=0.01,
                      model = "GLV",
                      derivative=EEMtoolbox::derivative_func){

  time_steps <- seq(from = t_window[1],to = t_window[2],by = time_step_len)#vector of time steps
  pars  <- list(interaction_matrix_value = interaction_matrix_value,
                growth_rate  = growth_rate,
                model=model)

  out <- deSolve::ode(y=initial_condition,
                      times=time_steps,
                      func=derivative,
                      parms=pars)

  return(out)
}
