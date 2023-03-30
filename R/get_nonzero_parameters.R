get_nonzero_parameters <- function(interaction_matrix){
  # This function returns an array of whether a parameter in the interaction
  # matrix is skipped (yes=1) or not (no=0), as well as the interaction sign
  # where a term is not skipped (interaction_terms_nonzero). 
  
  if (class(interaction_matrix)[1]=="matrix"){
    #determine which interaction terms are 0
    n_species <- ncol(interaction_matrix)
    
    #get all of the interaction terms as an array
    interaction_terms <- matrix(c(interaction_matrix),nrow=1,ncol=n_species^2)#check dimensions
    
    skip_parameters <- c(interaction_terms==0)
    
    #extract the non-zero terms
    interaction_terms_nonzero <- interaction_terms[which(!(skip_parameters))]
    # define uniform bounds on interaction terms
    lower_interaction_bound <- pmin(interaction_terms_nonzero, 0)
    upper_interaction_bound <- pmax(interaction_terms_nonzero, 0)
    
        
  } else {#list of upper and lower bound
    n_species <- ncol(interaction_matrix[[1]])
    
    #get all of the interaction terms as an array
    lower_interaction_terms <- c(interaction_matrix[[1]])
    lower_skip_parameters <- c(lower_interaction_terms==0)
    
    upper_interaction_terms <- c(interaction_matrix[[2]])
    upper_skip_parameters <- c(upper_interaction_terms==0)
    
    #extract the non-zero terms
    skip_parameters <- lower_skip_parameters*upper_skip_parameters
    
    # define uniform bounds on interaction terms
    lower_interaction_bound <- lower_interaction_terms[!as.logical(skip_parameters)]
    upper_interaction_bound <- upper_interaction_terms[!as.logical(skip_parameters)]
  }
  return(list(skip_parameters=as.numeric(skip_parameters),
              lower_interaction_bound=lower_interaction_bound,
              upper_interaction_bound=upper_interaction_bound))

}
