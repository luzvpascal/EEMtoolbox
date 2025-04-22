#' @title Calculate projections by solving ODEs
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
#' @param mode "recruitment" or "removal". Default "recruitment"
#' @param init_intervention_amount number of individuals in initial interventions, default = 0
#' @param init_intervention_interval initial intervention interval, default = 0
#' @param sustain_intervention_amount number of individuals in subsequent interventions, default = 0
#' @param sustain_intervention_interval subsequent intervention interval, default = NA
#' @param sustain_intervention_threshold abundance threshold where subsequent interventions stop, default = NA
#' @param introduced_species_index indicates which species in the index is the introduced one, default = 1
#' @examples
#' library(EEMtoolbox)
#' output <- EEM(matrix(c(-1,-1,1,-1),ncol=2)) #automatically loads an example of interaction matrix as dingo_matrix
#' calculate_projections(output,  c(1,1), t_window=c(0,1))
#' @return dataset of species abundances over time
#' @export
adapted_calculate_projections <-
  function(parameters,
           initial_condition,
           t_window,
           time_step_len=0.01,
           model = "GLV",
           derivative=EEMtoolbox::derivative_func,
           scaled=FALSE,
           species_names=NA,
           # Recruitment parameters:
           mode = "recruitment", # "recruitment" or "removal"
           init_intervention_amount = 0,
           init_intervention_timepoints = NA,
           sustain_intervention_amount = 0,
           sustain_intervention_timepoints = NA,
           sustain_intervention_threshold = NA,
           introduced_species_index = 1,
           multiplier = 1) {

    #Create a loading bar object to see where you're at
  pb <- txtProgressBar(min = 0,
                       max = length(parameters),
                       style = 3)

  # We'll pass the recruitment schedule via the parameters list.
  # Define the intervention times over the time window:
  recruitment_times <-
    sort(unique(c(init_intervention_timepoints, sustain_intervention_timepoints)))

  recruitment_event <- function(time, y, pars) {
    # Check if current time is an initial intervention time:
    # If so, add the initial intervention amount to the introduced species:
    if (time %in% pars$init_intervention_timepoints) {
      y[pars$introduced_species_index] <-
        #in this case, y is the current vector of species abundance at a given time
        y[pars$introduced_species_index] + pars$init_intervention_amount
    }
if (mode == "removal") {
  if (time %in% pars$sustain_intervention_timepoints &&
      y[pars$introduced_species_index] > pars$sustain_intervention_threshold) {
    y[pars$introduced_species_index] <-
      y[pars$introduced_species_index] + pars$sustain_intervention_amount
  }
} else if (mode == "recruitment") {
  # Check if current time is a sustaining intervention time AND abundance is below threshold:
  # If so, add the sustaining intervention amount to the introduced species:
  if (time %in% pars$sustain_intervention_timepoints &&
      y[pars$introduced_species_index] < pars$sustain_intervention_threshold) {
    y[pars$introduced_species_index] <-
      y[pars$introduced_species_index] + pars$sustain_intervention_amount
  }
}
    return(y)
  }

  #parameters we will need for the ode_solve_it function
  recruitment_times <- sort(unique(c(init_intervention_timepoints,
                                     sustain_intervention_timepoints))) # not sure if we need this!

  # Create a recruitment parameter list to pass to the ODE solver. This way, we
  # don't need to write them one after the other but rather keep them all
  # together.
  recruitment_pars <-
    list(init_intervention_timepoints = init_intervention_timepoints,
         sustain_intervention_timepoints = sustain_intervention_timepoints,
         init_intervention_amount = init_intervention_amount,
         sustain_intervention_amount = sustain_intervention_amount,
         sustain_intervention_threshold = sustain_intervention_threshold,
         introduced_species_index = introduced_species_index)

  ode_solve_it <- function(pars,
                           initial_condition,
                           t_window,
                           time_step_len,
                           model,
                           derivative,
                           scaled,
                           recruitment_event, #added -> function defined earlier
                           recruitment_times, #added -> we defined earlier. Do we need this?
                           recruitment_pars) { #added -> all the parameters related to recruitement, packed together.
    #if scaled: find steady state, initial condition is scaled to steady state
    if (scaled) { #if scaled == TRUE
      if (model %in% c("GLV", "Gompertz")) {
        steady_state <- solve(pars$interaction_matrix,-pars$growthrates)
        if (model == "Gompertz") {
          steady_state <- exp(steady_state)
        }
      } else if (model == "Bimler-Baker") {
        r <- pars$growthrates
        A <- pars$interaction_matrix_alphas
        B <- pars$interaction_matrix_betas

        R <- r
        P <- diag(A)
        M <- A - diag(P)

        fn <- function(N) {
          output <- R*(1 - exp(-M %*% N - P)) + B %*% N #change to positive
          return(output)
        }
        sol <- nleqslv::nleqslv(rep(100,length(r)), fn) #we give a large positive warmstart
        steady_state <- sol$x
      }
      initial <- initial_condition * steady_state #initial condition is scaled to steady state -> interpret initial as a multiplier
    } else {
      initial <- initial_condition #non-scaled to steady state
    }

    ## solve ODE
    if (model %in% c("GLV", "Gompertz")) {
      projections <-
        as.data.frame(adapted_ode_solve( #use adapted ode_solve function
          interaction_matrix_value = pars$interaction_matrix,
          growth_rate = pars$growthrates,
          initial,
          t_window,
          time_step_len,
          model,
          derivative,
          recruitment_pars = recruitment_pars,
          recruitment_event = recruitment_event,
          recruitment_times = recruitment_times))
    } else if (model == "Bimler-Baker") {
      projections <-
        as.data.frame(adapted_ode_solve( #use adapted ode_solve function
          interaction_matrix_value = list(pars$interaction_matrix_alphas,
                                          pars$interaction_matrix_betas),
          growth_rate = pars$growthrates,
          initial,
          t_window,
          time_step_len,
          model,
          derivative,
          recruitment_pars = recruitment_pars,
          recruitment_event = recruitment_event,
          recruitment_times = recruitment_times))
    }

    if (scaled) { #if scaled to steady state, scale abundances
      projections[,-1] <- projections[,-1]/matrix(steady_state,
                                                  ncol = length(steady_state),
                                                  nrow = nrow(projections),
                                                  byrow = TRUE)
    }
    return(projections)
  }

  abundance <- lapply(seq_len(length(parameters)),
                      function(i) {
                        setTxtProgressBar(pb, i)
                        ode_solve_it(parameters[[i]],
                                     model = model,
                                     initial_condition = initial_condition,
                                     t_window = t_window,
                                     time_step_len = time_step_len,
                                     derivative,
                                     scaled,
                                     recruitment_event = recruitment_event, # added
                                     recruitment_times = recruitment_times, # added
                                     recruitment_pars = recruitment_pars) # added
                      })

  abundance <- dplyr::bind_rows(abundance, .id = "sim")

  if (!is.na(species_names[1])) {
    names(abundance) <- c("sim", "time", species_names)
  }

  ## remove projections that could not be solved ####
  remove_indexes <-  dplyr::group_by(abundance, sim)
  remove_indexes <-  dplyr::summarise(remove_indexes, max_time = max(time))
  remove_indexes <-  dplyr::filter(remove_indexes, max_time < t_window[2])
  remove_indexes <- remove_indexes$sim

  if (length(remove_indexes) > 0) {
    print("Following parameter sets could not be solved and were removed:")
    print(remove_indexes)

    #remove parameter sets#
    # abundance <- dplyr::filter(abundance, !(sim %in% remove_indexes)) #this line was already in green
    abundance <- abundance[which(!(abundance$sim %in% remove_indexes)),]
  }

  #pivot for plotting
  abundance <- tidyr::pivot_longer(abundance,
                                   !c(time,sim), #select all columns except time and sim
                                   names_to = c("species"),
                                   values_to = "pop")
  if (scaled == FALSE) {
    abundance$pop <- abundance$pop*multiplier #multiply by the multiplier
    for (i in 1:nrow(abundance)) {
      if (abundance$pop[i] < 0) {
        abundance$pop[i] <- 0
        #here add something that since this point on the abundance has to stay 0 for this sim unless individuals are added OR do this earlier so that the GLV directly calculates from 0
      }
      }
    }

  return(abundance)
}
