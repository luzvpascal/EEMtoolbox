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

  ode_solve_it <- function(pars,
                           initial_condition,
                           t_window,
                           time_step_len,
                           model,
                           derivative,
                           scaled){
    #if scaled: find steady state, initial condition is scaled to steady state
    if (scaled) {
      if (model %in% c("GLV", "Gompertz")){
        steady_state <- solve(pars$interaction_matrix,-pars$growthrates)
        if (model=="Gompertz"){
          steady_state <- exp(steady_state)
        }
      } else if (model =="Bimler-Baker"){
        r <- pars$growthrates
        A <- pars$interaction_matrix_alphas
        B <- pars$interaction_matrix_betas

        R <- r
        P <- diag(A)
        M <- A-diag(P)

        fn <- function(N) {
          output <- R*(1-exp(-M%*%N-P))+B%*%N #change to positive
          return(output)
        }
        sol <- nleqslv::nleqslv(rep(100,length(r)), fn) #we give a large positive warmstart
        steady_state <- sol$x
      }
      initial <- initial_condition*steady_state #initial condition is scaled to steady state
    } else {
      initial <- initial_condition #non-scaled to steady state
    }
    ## solve ODE
    if (model %in% c("GLV", "Gompertz")){
      projections <- as.data.frame(EEMtoolbox::ode_solve(interaction_matrix_value = pars$interaction_matrix,
                                                         growth_rate=pars$growthrates,
                                                         initial,
                                                         t_window,
                                                         time_step_len,
                                                         model,
                                                         derivative))
    } else if (model =="Bimler-Baker") {
      projections <-as.data.frame(EEMtoolbox::ode_solve(interaction_matrix_value =
                                            list(pars$interaction_matrix_alphas,
                                                 pars$interaction_matrix_betas),
                                          growth_rate=pars$growthrates,
                                          initial,
                                          t_window,
                                          time_step_len,
                                          model,
                                          derivative))
    }

    if (scaled){#if scaled to steady state, scale abundances
      projections[,-1] <- projections[,-1]/matrix(steady_state,
                                                  ncol=length(steady_state),
                                                  nrow=nrow(projections),
                                                  byrow=TRUE)
    }
    return(projections)
  }

  abundance <- lapply(parameters,ode_solve_it, model=model,
                                            initial_condition=initial_condition,
                                            t_window=t_window,
                                            time_step_len=time_step_len,
                                            derivative,
                                            scaled)

  abundance <- dplyr::bind_rows(abundance, .id="sim")
  if (!is.na(species_names[1])){
    names(abundance) <- c("sim", "time", species_names)
  }

  ## remove projections that could not be solved ####
  remove_indexes <-  dplyr::group_by(abundance, sim)
  remove_indexes <-  dplyr::summarise(remove_indexes, max_time=max(time))
  remove_indexes <-  dplyr::filter(remove_indexes, max_time < t_window[2])
  remove_indexes <- remove_indexes$sim

  if (length(remove_indexes)>0){
    print("The ODE could not be solved for parameter sets (index):")
    print(remove_indexes)
    print("These parameter sets will be removed from the abundance predictions")

    #remove parameter sets#
    # abundance <- dplyr::filter(abundance, !(sim %in% remove_indexes))
    abundance <- abundance[which(!(abundance$sim %in% remove_indexes)),]
  }

  #pivot for plotting
  abundance <- tidyr::pivot_longer(abundance,
                                   !c(time,sim), names_to = c("species"), values_to = "pop")

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
