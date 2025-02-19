#' @title Adapted disc_func for EEM
#' @description
#' Discrepancy measure: summarises how infeasible or instable this system is + adds penalty for wrong equilibria
#' @param data vector of the equilibrium points followed by the real eigenvalues of the jacobian matrix, as returned by \link[EEMtoolbox]{summarise_ecosystem_features}
#' @param target_lower vector of the lower bound of the target equilibrium
#' @param target_upper vector of the upper bound of the target equilibrium
#' @return summary statistic (discrepancy measure).
#' @export
adapted_discrepancy_continuous_sum <- function(data,
                                               target_lower = NULL,
                                               target_upper = NULL) {
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

  # If target lower and upper bounds are provided, calculate per-species penalties.
  if (!is.null(target_lower) && !is.null(target_upper)) {
    if (length(target_lower) != n_species || length(target_upper) != n_species)
    {
      stop("target_lower and target_upper must be vectors of
           length equal to the number of species.")
    }

    # For each species, if the computed equilibrium is below target_lower, add (target_lower - computed).
    lower_penalty <- pmax(0, target_lower - equilibrium_points)
    # For each species, if the computed equilibrium is above target_upper, add (computed - target_upper).
    upper_penalty <- pmax(0, equilibrium_points - target_upper)

    # Sum the penalties across all species.
    target_penalty <- sum(lower_penalty + upper_penalty)
  }

  # Total discrepancy is the sum of feasibility, instability, and target penalties.
  summ <- measured_infeasibility + measured_instability + target_penalty
  return(summ)
}
