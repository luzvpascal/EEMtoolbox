#' @title Solve ODE
#' @description
#' Solves ODEs
#' @param parameters list like object of ensemble of parameters (outputs of  \link[EEMtoolbox]{EEM})
#' @param initial_condition vector of initial species abundances
#' @param t_window time window to solve ODE
#' @param time_step_len length of each time step, default = 0.01
#' @param model model representing species interactions. Default "GLV" (Generalized Lotka Volterra). options include "Baker" and "Gompertz".
#' @param derivative derivative function. Default \link[EEMtoolbox]{derivative_func}
#' @examples
#' library(EEMtoolbox)
#' output <- EEM(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' projections(output, t_window=c(0,1))
#' @return matrix: first column (time) vector of time steps of the simulation, the following columns correspond to the values of abundances at a given time step
#' @export

projections <- function(parameters,
                        initial_condition,
                        t_window,
                        time_step_len=0.01,
                        model = "GLV",
                        derivative=EEMtoolbox::derivative_func){

  ode_solve_it <- function(pars,
                           initial_condition,
                           t_window,
                           time_step_len,
                           model,
                           derivative){
    return(EEMtoolbox::ode_solve(interaction_matrix_value = pars$interaction_matrix,
                                 growth_rate=pars$growthrates,
                                 initial_condition,
                                 t_window,
                                 time_step_len,
                                 model,
                                 derivative))
  }

  abundance <- lapply(outputs,ode_solve_it, model=model,
                                            initial_condition=initial_condition,
                                            t_window=t_window,
                                            time_step_len=time_step_len,
                                            derivative)

  return(abundance)
}
