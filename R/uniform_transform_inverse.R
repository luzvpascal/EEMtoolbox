#' @title inverse transformation function
#' @description
#' This function uses a reverts the transformation for uniform distributions
#' @param theta vector of parameters to transformed back
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @return vector of transformed parameters
#' @export
#'
uniform_transform_inverse <- function(theta_trans,prior_args){
  # This function inverses the transform for uniform distributions.
  lower_mat <- matrix(prior_args$lower, nrow=nrow(theta_trans), ncol=ncol(theta_trans),
                      byrow=TRUE)
  upper_mat <- matrix(prior_args$upper, nrow=nrow(theta_trans), ncol=ncol(theta_trans),
                      byrow=TRUE)

  theta <- (lower_mat+upper_mat*exp(theta_trans))/(1+exp(theta_trans))
  return(theta)
}
