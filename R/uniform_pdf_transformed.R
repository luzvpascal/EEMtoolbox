uniform_pdf_transformed <- function(theta_trans){
  # This function calculates the density of the transformed uniform
  # distribution
  theta <- prod(abs(exp(theta_trans)/((exp(theta_trans)+1)^2)))
  return(theta)
}
