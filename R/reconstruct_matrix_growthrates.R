#' @title Reconstruct interaction matrix and growth rates vector from set of parameters
#'
#' @param parameters a vector of sampled parameters
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @return A list of 2 elements growthrates: a vector of arguments defining the problem. interaction_matrix: reconstructed interaction matrix
#' @export

reconstruct_matrix_growthrates <-function(parameters,sim_args){
  ##extract information
  n_species <- sim_args$n_species
  r <- parameters[seq(n_species)]

  #reconstruct interaction matrix from array of parameters
  A_nonzero_values <- parameters[seq(n_species+1, length(parameters))]
  A_values <- rep(0, n_species^2)
  A_values[as.logical(1-sim_args$skip_parameters)] <- A_nonzero_values
  A <- matrix(A_values,ncol=n_species,nrow=n_species)

  return(list(growthrates=r, interaction_matrix=A))
}
