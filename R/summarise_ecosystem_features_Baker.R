summarise_ecosystem_features_GLV <- function(parameters,sim_args){
  ## to simulate this ecosystem, we simply calculate the equilibrium
  # abundances and the stability. for the Generalized Lokta Volterra model
  
  ##extract information
  n_species <- sim_args$n_species
  r <- parameters[seq(n_species)]
  
  #To optimize
  #reconstruct A matrix from array of parameters
  A_nonzero_values <- parameters[seq(n_species+1, length(parameters))]
  A_values <- rep(0, n_species^2)
  A_values[as.logical(1-sim_args$skip_parameters)] <- A_nonzero_values
  A <- matrix(A_values,ncol=n_species,nrow=n_species)
  
  # FEASIBILITY CHECK
  #find equilibrium abundances for feasibility
  equilibrium_points <- solve(A,-r)
  
  # STABILITY CHECK
  #calculate Jacobian for stability ###CHECK THIS PROCESS
  jacobian <- A*matrix(equilibrium_points, ncol=n_species, nrow=n_species)
  diag(jacobian) <- diag(jacobian) + A%*%(equilibrium_points) + r
  #check stability
  stability_eigenvalues <- Re(eigen(jacobian)$values)
  
  #return features
  summarised_features <- c(equilibrium_points,
                              stability_eigenvalues)
  return(summarised_features)
}
