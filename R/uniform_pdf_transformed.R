#' @title density of transformed uniform distribution
#' @description
#' A short description...
#'
#' @param theta_trans vector of parameters to transform
#' @return vector of transformed parameters
#' @export

uniform_pdf_transformed <- function(theta_trans){
  # This function calculates the density of the transformed uniform
  # distribution
  theta <- prod(abs(exp(theta_trans)/((exp(theta_trans)+1)^2)))
  return(theta)
}
