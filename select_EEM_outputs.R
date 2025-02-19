select_EEM_outputs <- function(ensemble,
                               target_lower,
                               target_upper,
                               sim_args) {
  # ensemble: a list of parameter sets from EEM.
  # target_equilibrium: the vector of desired species abundances (e.g., c(10000, 100, ...)).
  # tolerance: a threshold on the difference (relative error) between the normalized computed equilibrium and target.
  # sim_args: the simulation arguments, created by (e.g.) args_function, needed for reconstruction.

  # Prepare an empty list to collect selected parameter sets
  selected <- list()

  for (i in seq_along(ensemble)) {
    param_set <- ensemble[[i]]

    # Reconstruct the parameters: growth rates and interaction matrix.
    reconstructed <- EEMtoolbox::reconstruct_matrix_growthrates(param_set, sim_args)
    r <- reconstructed$growthrates$growthrates
    A <- reconstructed$growthrates$interaction_matrix

    # Compute the equilibrium: A %*% N + r = 0  =>  N = solve(A, -r)
    current_eq <- tryCatch(solve(A, -r),
                           error = function(e) rep(NA, sim_args$n_species))

    if (any(is.na(current_eq))) {
      cat("Parameter set", i, "produced NA equilibrium.\n")
      next
    } else {
      cat("Parameter set", i, "equilibrium:", current_eq, "\n")
    }

    # Check for NA values (if the system is unsolvable, skip this set)
    if (any(is.na(current_eq))) {
      next
    }

    # If the equilibrium is within the bounds, select this parameter set.
    if (length(which(current_eq > target_lower)) == sim_args$n_species &&
        length(which(current_eq < target_upper)) == sim_args$n_species) {
      selected[[length(selected) + 1]] <- param_set
    }
  }

  return(selected)
}
