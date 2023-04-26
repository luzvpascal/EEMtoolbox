#' @title Reconstruct interaction matrix and growth rates vector from set of parameters
#'
#' @param parameters a vector of sampled parameters
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @return A list of 3 elements growthrates: a vector of arguments defining the problem. alphas_matrix=matrix of alpha elements and interaction_matrix: reconstructed interaction matrix
#' @export

reconstruct_matrix_growthrates_Baker <-function(parameters,sim_args){
  ##extract information
  n_species <- sim_args$n_species
  r <- parameters[seq(n_species)]

  #reconstruct Matrix of alphas from array of parameters
  A_nonzero_values <- parameters[seq(n_species+1, n_species+n_species**2)]
  A_values <- rep(0, n_species^2)
  A_values[as.logical(1-sim_args$skip_parameters_alphas)] <- A_nonzero_values
  A <- matrix(A_values,ncol=n_species,nrow=n_species)

  B_nonzero_values <- parameters[seq(n_species+n_species**2+1, length(parameters))]
  B_values <- rep(0, n_species^2)
  B_values[as.logical(1-sim_args$skip_parameters)] <- B_nonzero_values
  B <- matrix(B_values,ncol=n_species,nrow=n_species)

  return(list(growthrates=r, alphas_matrix=A, interaction_matrix=B))
}
