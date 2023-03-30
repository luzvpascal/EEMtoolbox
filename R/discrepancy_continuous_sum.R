discrepancy_continuous_sum <- function(data){
##
# Calculation of summary statistics
#the data passed in is the equilibrium points followed by the real
#eigenvalues of the jacobian matrix.
#summarises how infeasible or instable this system is

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
