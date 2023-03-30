uniform_transform <- function(theta,prior_args){
  # This function uses a transform for uniform distributions, so that all
  # values are within the uniform bounds.
  lower_mat <- matrix(prior_args$lower, nrow=nrow(theta), ncol=ncol(theta),
                      byrow=TRUE)
  upper_mat <- matrix(prior_args$upper, nrow=nrow(theta), ncol=ncol(theta),
                      byrow=TRUE)
  return(log((theta-lower_mat)/(upper_mat-theta)))

}