#' @title Discrepancy function
#' @description
#' Discrepancy measure: summarises how infeasible or instable this system is
#' @param data vector of the equilibrium points followed by the real eigenvalues of the jacobian matrix, as returned by \link[EEMtoolbox]{summarise_ecosystem_features}
#' @return summary statistic (discrepancy measure).
#' @export

discrepancy_continuous_sum <- function(data){

  #extract information
  n_species <- length(data)/2
  equilibrium_points <- data[seq(n_species)]
  stability_eigenvalues <- data[seq(n_species+1,length(data))]

  #summarise in terms of feasibility and stability
  measured_infeasibility <- abs(sum(pmin(0,equilibrium_points)))
  measured_instability <- sum(pmax(0,stability_eigenvalues))
  summ <- measured_infeasibility+measured_instability
  return(summ)
}
