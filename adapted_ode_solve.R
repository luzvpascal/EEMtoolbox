adapted_ode_solve <- function(initial_condition,
                              interaction_matrix_value,
                              growth_rate,
                              t_window,
                              time_step_len = 0.01,
                              model = "GLV",
                              derivative = EEMtoolbox::derivative_func,
                              recruitment_pars,
                              recruitment_event,
                              recruitment_times){

  time_steps <- seq(from = t_window[1],to = t_window[2],by = time_step_len) #vector of time steps

  pars  <- list(interaction_matrix_value = interaction_matrix_value,
                growth_rate  = growth_rate,
                model = model,
                recruitment_pars = recruitment_pars) #added, parameters related to recruitement

  out <- deSolve::ode(y = initial_condition,
                      times = time_steps,
                      func = derivative,
                      parms = pars,
                      events = list(func = recruitment_event,
                                    time = recruitment_times))

  return(out)
}
