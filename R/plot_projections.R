#' @title Solve ODE and plot solutions
#' @description
#' Solves ODEs and plot solutions
#' @param parameters list like object of ensemble of parameters (outputs of  \link[EEMtoolbox]{EEM})
#' @param initial_condition vector of initial species abundances. If parameter scaled is TRUE, the parameter initial_condition should be scaled to steady state
#' @param t_window time window to solve ODE
#' @param time_step_len length of each time step, default = 0.01
#' @param model model representing species interactions. Default "GLV" (Generalized Lotka Volterra). options include "Bimler-Baker" and "Gompertz".
#' @param derivative derivative function. Default \link[EEMtoolbox]{derivative_func}
#' @param scaled Boolean indicating if projections should be scaled to steady state. If true, the parameter initial_condition should be scaled too. Default FALSE
#' @param species_names vector of strings for names of species. If NA plots only display species index number, . Default NA.
#' @param average Boolean indicating if projections should be averaged, or indivially ploted
#' @examples
#' library(EEMtoolbox)
#' output <- EEM(matrix(c(-1,-1,1,-1),ncol=2)) #automatically loads an example of interaction matrix as dingo_matrix
#' plot_projections(output,  c(1,1), t_window=c(0,1))
#' @return ggplot of abundances per species
#' @export
plot_projections <- function(parameters,
                        initial_condition,
                        t_window,
                        time_step_len=0.01,
                        model = "GLV",
                        derivative=EEMtoolbox::derivative_func,
                        scaled=FALSE,
                        species_names=NA,
                        average=TRUE){

  abundance <- calculate_projections(parameters,
                                      initial_condition,
                                      t_window,
                                      time_step_len,
                                      model,
                                      derivative,
                                      scaled,
                                      species_names)
  # Add either the ribbon or individual trajectories based on the average variable
  if (average) {
    abundance <- dplyr::group_by(abundance, time, species)
    abundance <- dplyr::summarise(abundance,
                                  median_pop = median(pop),
                                  upper = quantile(pop, 0.975),
                                  lower = quantile(pop, 0.025)
                                  )

    p <- ggplot2::ggplot(abundance) +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Species"),
                      color = ggplot2::guide_legend(title = "Species")) +
      ggplot2::xlab("Time") +
      ggplot2::ylab("Abundance") +
      ggplot2::facet_wrap(~ species) +
      ggplot2::geom_ribbon(ggplot2::aes(x = time,
                                        ymin = lower,
                                        ymax = upper,
                                        color = species,
                                        fill = species),
                           alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(x = time,
                                      y = median_pop,
                                      color = species)
                         ,linewidth = 1.5)
  } else {
    p <- ggplot2::ggplot(abundance) +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Species"),
                      color = ggplot2::guide_legend(title = "Species")) +
      ggplot2::xlab("Time") +
      ggplot2::ylab("Abundance") +
      ggplot2::facet_wrap(~ species) +
      ggplot2::geom_line(ggplot2::aes(group=sim), alpha = 0.4)+
      ggplot2::stat_summary(geom = "line", fun = mean, col="black", linewidth = 1)
  }
  return(p)
}
