#' @title Summary of ecosystem features for Baker model
#' @description
#' Tests the feasibility and stability of a vector of sampled parameters for the Baker model
#'
#' @param parameters a vector of sampled parameters
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @return vector of values: first half are the steady states (indicating feasibility) and second half the eigen values of Jacobian (indicating stability)
#' @export

summarise_ecosystem_features_Baker <- function(parameters,sim_args){

  ##extract information
  n_species <- sim_args$n_species

  #reconstruct interaction matrix and vector of growth rates from array of parameters
  reconstruct <- EEMtoolbox::reconstruct_matrix_growthrates(parameters,sim_args)
  r <- reconstruct$growthrates
  A <- reconstruct$interaction_matrix

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
