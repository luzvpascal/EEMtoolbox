adapted_ode_solve <- function(initial_condition,
                              interaction_matrix_value,
                              growth_rate,
                              t_window,
                              time_step_len,
                              model = "GLV",
                              derivative = EEMtoolbox::derivative_func,
                              recruitment_pars,
                              recruitment_event,
                              recruitment_times,
                              extinction_threshold){

  time_steps <- round(seq(from = t_window[1],
                          to = t_window[2],
                          by = time_step_len),
                      nchar(strsplit(
                        as.character(time_step_len), "\\.")[[1]][2])) #vector of time steps -> to make sure it is actually the correct sequence, there is a bug with the seq() of r

  pars  <- c(list(interaction_matrix_value = interaction_matrix_value,
                  growth_rate  = growth_rate,
                  model = model),
             recruitment_pars, #added, parameters related to recruitment
             list(extinction_threshold = extinction_threshold))

  out <- deSolve::ode(y = initial_condition,
                      times = time_steps,
                      func = derivative,
                      parms = pars,
                      events = list(func = recruitment_event,
                                    time = recruitment_times)) # added

  return(out)
}
