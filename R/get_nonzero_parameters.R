#' @title List of non-zero parameters
#' @description
#' Array of indicating if a parameter in the interaction matrix is skipped (yes=1) or not (no=0), as well as the interaction sign (upper and lower bounds)
#' @param interaction_matrix interaction signs matrix. If model is GLV or Gompertz it can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second.     #if model is Baker, the interaction_matrix has to be a list of two lists, the first list contains matrices defining lower and upper bounds of alphas, the second list contains matrices defining lower and upper bounds of betas
#' @param model model representing species interactions, default "GLV" (Generalized Lokta Voltera). options include "Baker", "Adams" and "customized"
#' @return A list
#' keep_parameters : parameter kept (yes=1) or not (no=0)
#' lower_interaction_bound: lower bounds of each parameter sign
#' upper_interaction_bound: upper bounds of each parameter sign
#' num_params: total number of unknown parameters
#'
#' Baker model only: keep_parameters_alphas : parameter kept (yes=1) or not (no=0)
#' Baker model only: keep_parameters_betas : parameter kept (yes=1) or not (no=0)

#' Baker model only: num_params_alphas: total number of unknown alpha parameters
#' Baker model only: num_params_betas: total number of unknown beta parameters

#' Baker model only: lower_interaction_bound_alphas: lower bounds of each alpha parameter sign
#' Baker model only: upper_interaction_bound_alphas: upper bounds of each alpha parameter sign

#' Baker model only: lower_interaction_bound_betas: lower bounds of each beta parameter sign
#' Baker model only: upper_interaction_bound_betas: upper bounds of each beta parameter sign

#' @export
get_nonzero_parameters <- function(interaction_matrix,n_species,model="GLV"){
  # This function returns an array of whether a parameter in the interaction
  # matrix is skipped (yes=1) or not (no=0), as well as the interaction sign
  # where a term is not skipped (interaction_terms_nonzero).

  if (model=="GLV"|model=="Gompertz"){
    list_nonzero_parameters <- EEMtoolbox::get_nonzero_parameters_mat(interaction_matrix)

    return(list(keep_parameters=as.numeric(list_nonzero_parameters$keep_parameters),
                lower_interaction_bound=list_nonzero_parameters$lower_interaction_bound,
                upper_interaction_bound=list_nonzero_parameters$upper_interaction_bound,
                num_params= n_species+ sum(as.numeric(list_nonzero_parameters$keep_parameters))
                )
           )
  } else if (model=="Baker"){
    list_nonzero_alphas <- EEMtoolbox::get_nonzero_parameters_mat(interaction_matrix[[1]])
    list_nonzero_betas <- EEMtoolbox::get_nonzero_parameters_mat(interaction_matrix[[2]])

    return(list(
      keep_parameters_alphas=as.numeric(list_nonzero_alphas$keep_parameters),
      keep_parameters_betas=as.numeric(list_nonzero_betas$keep_parameters),

      lower_interaction_bound_alphas=list_nonzero_alphas$lower_interaction_bound,
      upper_interaction_bound_alphas=list_nonzero_alphas$upper_interaction_bound,

      lower_interaction_bound_betas=list_nonzero_betas$lower_interaction_bound,
      upper_interaction_bound_betas=list_nonzero_betas$upper_interaction_bound,

      lower_interaction_bound=c(list_nonzero_alphas$lower_interaction_bound,
                                list_nonzero_betas$lower_interaction_bound),

      upper_interaction_bound=c(list_nonzero_alphas$upper_interaction_bound,
                                list_nonzero_betas$upper_interaction_bound),

      num_params= n_species+
        sum(as.numeric(list_nonzero_alphas$keep_parameters)) +
        sum(as.numeric(list_nonzero_betas$keep_parameters)),

      num_params_alphas = sum(as.numeric(list_nonzero_alphas$keep_parameters)),
      num_params_betas = sum(as.numeric(list_nonzero_betas$keep_parameters))
    )
    )
  }

}
