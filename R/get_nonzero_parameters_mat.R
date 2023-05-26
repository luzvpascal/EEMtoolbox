#' @title List of non-zero parameters
#' @description
#' Array of indicating if a parameter in the interaction matrix is skipped (yes=1) or not (no=0), as well as the interaction sign (upper and lower bounds)
#' @param interaction_matrix interaction signs matrix. If model is GLV or Gompertz it can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second.     #if model is Baker, the interaction_matrix has to be a list of two lists, the first list contains matrices defining lower and upper bounds of alphas, the second list contains matrices defining lower and upper bounds of betas
#' @return A list
#' keep_parameters : parameter skipped (yes=1) or not (no=0)
#' lower_interaction_bound: lower bounds of each parameter sign
#' upper_interaction_bound: upper bounds of each parameter sign
#' @export
get_nonzero_parameters_mat <- function(interaction_matrix){
  # This function returns an array of whether a parameter in the interaction
  # matrix is skipped (yes=1) or not (no=0), as well as the interaction sign
  # where a term is not skipped (interaction_terms_nonzero).

  if (class(interaction_matrix)[1]=="matrix"){
    #determine which interaction terms are 0
    n_species <- ncol(interaction_matrix)

    #get all of the interaction terms as an array
    interaction_terms <- matrix(c(interaction_matrix),nrow=1,ncol=n_species^2)#check dimensions

    keep_parameters <- c(interaction_terms!=0)

    #extract the non-zero terms
    interaction_terms_nonzero <- interaction_terms[which((keep_parameters))]
    # define uniform bounds on interaction terms
    lower_interaction_bound <- pmin(interaction_terms_nonzero, 0)
    upper_interaction_bound <- pmax(interaction_terms_nonzero, 0)


  } else {#list of upper and lower bound
    n_species <- ncol(interaction_matrix[[1]])

    #get all of the interaction terms as an array
    lower_interaction_terms <- c(interaction_matrix[[1]])
    lower_keep_parameters <- c(lower_interaction_terms!=0)

    upper_interaction_terms <- c(interaction_matrix[[2]])
    upper_keep_parameters <- c(upper_interaction_terms!=0)

    #extract the non-zero terms
    keep_parameters <- lower_keep_parameters|upper_keep_parameters

    # define uniform bounds on interaction terms
    lower_interaction_bound <- lower_interaction_terms[as.logical(keep_parameters)]
    upper_interaction_bound <- upper_interaction_terms[as.logical(keep_parameters)]
  }
  return(list(keep_parameters=as.numeric(keep_parameters),
              lower_interaction_bound=lower_interaction_bound,
              upper_interaction_bound=upper_interaction_bound))

}
