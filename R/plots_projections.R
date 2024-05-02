#' @title Solve ODE and plot solutions
#' @description
#' Solves ODEs and plot solutions
#' @param parameters list like object of ensemble of parameters (outputs of  \link[EEMtoolbox]{EEM})
#' @param initial_condition vector of initial species abundances
#' @param t_window time window to solve ODE
#' @param time_step_len length of each time step, default = 0.01
#' @param model model representing species interactions. Default "GLV" (Generalized Lotka Volterra). options include "Baker" and "Gompertz".
#' @param derivative derivative function. Default \link[EEMtoolbox]{derivative_func}
#' @param species_names vector of strings for names of species. If NA plots only display species index number, . Default NA.
#' @examples
#' library(EEMtoolbox)
#' output <- EEM(matrix(c(-1,-1,1,-1),ncol=2)) #automatically loads an example of interaction matrix as dingo_matrix
#' plots_projections(output,  c(1,1), t_window=c(0,1))
#' @return ggplot of abundances per species
#' @export

plots_projections <- function(parameters,
                        initial_condition,
                        t_window,
                        time_step_len=0.01,
                        model = "GLV",
                        derivative=EEMtoolbox::derivative_func,
                        species_names=NA){

  ode_solve_it <- function(pars,
                           initial_condition,
                           t_window,
                           time_step_len,
                           model,
                           derivative){
    return(as.data.frame(EEMtoolbox::ode_solve(interaction_matrix_value = pars$interaction_matrix,
                                 growth_rate=pars$growthrates,
                                 initial_condition,
                                 t_window,
                                 time_step_len,
                                 model,
                                 derivative)))
  }

  abundance <- lapply(parameters,ode_solve_it, model=model,
                                            initial_condition=initial_condition,
                                            t_window=t_window,
                                            time_step_len=time_step_len,
                                            derivative)

  abundance <- dplyr::bind_rows(abundance)
  if (!is.na(species_names[1])){
    names(abundance) <- c("time", species_names)
  }

  abundance <- tidyr::pivot_longer(abundance,
                        !time, names_to = c("species"), values_to = "pop")

  p <- ggplot2::ggplot(abundance, ggplot2::aes(x = time,
                                               y = pop,
                                               color = species,
                                               fill = species)) +
    ggplot2::stat_summary(geom = "line", fun = mean,linewidth = 1.5) +
    ggplot2::stat_summary(geom = "ribbon", fun.data = function(x) {
      quantiles <- quantile(x, c(0.025, 0.975))
      data.frame(ymin = quantiles[1], ymax = quantiles[2])
    }, alpha = 0.2) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Species"),
           color = ggplot2::guide_legend(title = "Species")) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Years") +
    ggplot2::ylab("Abundance") +
    ggplot2::ylim(c(0, quantile(abundance$pop, 0.975))) +
    ggplot2::facet_wrap(~species)

  return(p)
}
