#' @title Reconstruct interaction matrix and growth rates vector from set of parameters
#'
#' @param parameters a vector of sampled parameters
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @return A list of 2 elements
#' growthrates: a vector of arguments defining the problem,
#' interaction_matrix: reconstructed interaction matrix.
#' If model is Bimler-Baker,
#' interaction_matrix_alphas: reconstructed interaction matrix of alphas
#' interaction_matrix_betas: reconstructed interaction matrix of betas
#' @export

reconstruct_matrix_growthrates <-function(parameters,sim_args){
  ##extract information
  n_species <- sim_args$n_species
  r <- parameters[seq(n_species)]

  #reconstruct interaction matrix from array of parameters
  if (sim_args$model=="GLV"|sim_args$model=="Gompertz"){

    A_nonzero_values <- parameters[seq(n_species+1, length(parameters))]
    A_values <- rep(0, n_species^2)
    A_values[as.logical(sim_args$keep_parameters)] <- A_nonzero_values
    A <- matrix(A_values,ncol=n_species,nrow=n_species)

    return(list(growthrates=r, interaction_matrix=A))

  } else if (sim_args$model=="Bimler-Baker"){

    A_nonzero_values <- parameters[n_species+seq(1, sim_args$num_params_alphas)]
    A_values <- rep(0, n_species^2)
    A_values[as.logical(sim_args$keep_parameters_alphas)] <- A_nonzero_values
    A <- matrix(A_values,ncol=n_species,nrow=n_species)

    B_nonzero_values <- parameters[n_species+sim_args$num_params_alphas+
                                     seq(1, sim_args$num_params_betas)]
    B_values <- rep(0, n_species^2)
    B_values[as.logical(sim_args$keep_parameters_betas)] <- B_nonzero_values
    B <- matrix(B_values,ncol=n_species,nrow=n_species)

    return(list(growthrates=r, interaction_matrix_alphas=A,
                interaction_matrix_betas=B))

  }

}
