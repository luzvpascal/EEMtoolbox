summarise_ecosystem_features_Adams <- function(parameters,sim_args){
  ## to simulate this ecosystem, we simply calculate the equilibrium
  # abundances and the stability. for the Generalized Lokta Volterra model
  
  ##extract information
  n_species <- sim_args$n_species
  
  #Ai = ln(ri)
  r <- parameters[seq(n_species)]
  A <- log(r)
  
  #reconstruct B matrix from array of parameters
  B_nonzero_values <- parameters[seq(n_species+1, length(parameters))]
  B_values <- rep(0, n_species^2)
  B_values[as.logical(1-sim_args$skip_parameters)] <- B_nonzero_values
  B <- matrix(B_values,ncol=n_species,nrow=n_species)
  
  # FEASIBILITY CHECK
  #find equilibrium abundances for feasibility
  equilibrium_points_X <- solve(diag(n_species)-B,A)
  equilibrium_points <- exp(equilibrium_points_X)
  
  # STABILITY CHECK
  #calculate Jacobian for stability
  Ni_matrix <- matrix(equilibrium_points, nrow = n_species, ncol=n_species)
  Nj_matrix <- t(Ni_matrix)
  jacobian <- B*Ni_matrix/Nj_matrix - diag(n_species)
  #check stability
  stability_eigenvalues <- Re(eigen(jacobian)$values)
  
  #return features
  summarised_features <- c(equilibrium_points,
                           stability_eigenvalues)
  return(summarised_features)
}
