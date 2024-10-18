#' @title Summary of ecosystem features for the Bimler-Baker model
#' @description
#' Tests the feasibility and stability of a vector of sampled parameters for the Bimler-Baker model
#'
#' @param parameters a vector of sampled parameters
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @return vector of values: first half are the steady states (indicating feasibility) and second half the eigen values of Jacobian (indicating stability)
#' @export
#' @import nleqslv

summarise_ecosystem_features_Baker <- function(parameters,sim_args){

  ##extract information
  n_species <- sim_args$n_species

  #reconstruct interaction matrix and vector of growth rates from array of parameters
  # reconstruct <- EEMtoolbox::reconstruct_matrix_growthrates_Baker(parameters,sim_args)
  reconstruct <- EEMtoolbox::reconstruct_matrix_growthrates(parameters,sim_args)
  r <- reconstruct$growthrates
  A <- reconstruct$interaction_matrix_alphas
  B <- reconstruct$interaction_matrix_betas

  # FEASIBILITY CHECK
  #find equilibrium abundances for feasibility
  R <- r
  P <- diag(A)
  M <- A-diag(P)

  fn <- function(N) {
    output <- R*(1-exp(-M%*%N-P))+B%*%N #change to positive
    return(output)
  }

  sol <- nleqslv::nleqslv(rep(100,n_species), fn) #we give a large positive warmstart

  equilibrium_points <- sol$x
  # STABILITY CHECK
  #calculate Jacobian for stability ###CHECK THIS PROCESS
  diag_elements <- diag(B)*equilibrium_points#change to positive

  R_matrix <- matrix(R, ncol = length(R), nrow = length(R))
  N_matrix <- matrix(equilibrium_points,
                     ncol = length(equilibrium_points),
                     nrow = length(equilibrium_points))
  exp_matrix <- matrix(exp(-M%*%equilibrium_points-P),
                       ncol = length(R), nrow = length(R))

  jacobian <- R_matrix*N_matrix*A*exp_matrix + B*N_matrix#change to positive

  diag(jacobian) <- diag_elements
  #check stability
  stability_eigenvalues <- Re(eigen(jacobian)$values)

  #return features
  summarised_features <- c(equilibrium_points,
                           stability_eigenvalues)
  return(summarised_features)
}
