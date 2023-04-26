#' @title Summary of ecosystem features for Gompertz model
#' @description
#' Tests the feasibility and stability of a vector of sampled parameters for the Gompertz model
#'
#' @param parameters a vector of sampled parameters
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @return vector of values: first half are the steady states (indicating feasibility) and second half the eigen values of Jacobian (indicating stability)
#' @export

summarise_ecosystem_features_Gompertz <- function(parameters,sim_args){

  ##extract information
  n_species <- sim_args$n_species

  #reconstruct interaction matrix and vector of growth rates from array of parameters
  reconstruct <- EEMtoolbox::reconstruct_matrix_growthrates(parameters,sim_args)
  r <- reconstruct$growthrates
  B <- reconstruct$interaction_matrix

  # FEASIBILITY CHECK
  #find equilibrium abundances for feasibility
  equilibrium_points_X <- solve(-B,r)
  equilibrium_points <- exp(equilibrium_points_X)

  # STABILITY CHECK
  #calculate Jacobian for stability
  # Ni_matrix <- matrix(equilibrium_points, nrow = n_species, ncol=n_species)
  # Nj_matrix <- t(Ni_matrix)
  jacobian <- B #*Ni_matrix/Nj_matrix
  #check stability
  stability_eigenvalues <- Re(eigen(jacobian)$values)

  #return features
  summarised_features <- c(equilibrium_points,
                           stability_eigenvalues)
  return(summarised_features)
}
