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
  # reconstruct <- EEMtoolbox::reconstruct_matrix_growthrates_Baker(parameters,sim_args)
  reconstruct <- EEMtoolbox::reconstruct_matrix_growthrates(parameters,sim_args)
  r <- reconstruct$growthrates
  mat <- reconstruct$interaction_matrix

  A <- pmax(0, mat)
  B <- pmin(0, mat)

  # FEASIBILITY CHECK
  #find equilibrium abundances for feasibility
  R <- r
  P <- diag(A)
  M <- A-P

  fn <- function(N) {
    output <- R*(1-exp(-M%*%N-P))-B%*%N
    return(output)
  }

  sol <- nleqslv::nleqslv(rep(0,n_species), fn)

  # STABILITY CHECK  TODO
  #calculate Jacobian for stability ###CHECK THIS PROCESS
  diag_elements <- -diag(B)*N_sol

  R_matrix <- matrix(R, ncol = length(R), nrow = length(R))
  N_matrix <- matrix(N_sol, ncol = length(N_sol), nrow = length(N_sol))

  jabobian <- R_matrix*N_matrix*A*exp(-M%*%N_sol-P)-B*N_matrix*t(N_matrix)

  diag(jabobian) <- diag_elements
  #check stability
  stability_eigenvalues <- Re(eigen(jacobian)$values)

  #return features
  summarised_features <- c(equilibrium_points,
                           stability_eigenvalues)
  return(summarised_features)
}
