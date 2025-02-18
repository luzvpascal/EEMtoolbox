adapted_discrepancy_continuous_sum <-
function(data, 
                                               target_equilibrium = NULL, 
                                               tolerance = 0.1) {
  n_species <- length(data) / 2
  equilibrium_points <- data[seq_len(n_species)]
  # computed equilibrium abundances
  stability_eigenvalues <- data[(n_species + 1):length(data)]
  # stability eigenvalues
  
  # Feasibility penalty: sum of negative equilibrium parts.
  measured_infeasibility <- abs(sum(pmin(0, equilibrium_points)))
  # Stability penalty: sum of positive eigenvalues.
  measured_instability <- sum(pmax(0, stability_eigenvalues))
  #until here, this is the same function as the original one
  
  
  # Initialize target penalty. is 0 if no target equilibrium is provided
  target_penalty <- 0
  
  # If a target equilibrium is provided, compute a penalty for deviation.
  if (!is.null(target_equilibrium)) {
    #if target equilibrium is not equal to the number of species
    if (length(target_equilibrium) != n_species) {
      stop("Length of target_equilibrium must equal the number of species.")
    }
    
   # For each species, compute the relative difference:
    # relative_diff[i] = |computed[i] - target[i]| / target[i]
    relative_diff <- abs(equilibrium_points - target_equilibrium) / target_equilibrium
    
    # For each species, if the relative difference exceeds tol_percent, then penalty is the excess.
    penalties <- pmax(0, relative_diff - tolerance)
    
    # Sum over species to get the overall target penalty.
    target_penalty <- sum(penalties)
  }
  
  # Total discrepancy is the sum of feasibility, instability, and target penalties.
  summ <- measured_infeasibility + measured_instability + target_penalty
  return(summ)
}
