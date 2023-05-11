#' @title Summary of ecosystem features
#' @description
#' Tests the feasibility and stability of a vector of sampled parameters
#' @param parameters a vector of sampled parameters
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @return vector of values: first half are the steady states (indicating feasibility) and second half the eigen values of Jacobian (indicating stability)
#' @export

summarise_ecosystem_features <- function(parameters,sim_args){
  ## to simulate this ecosystem, we simply calculate the equilibrium
  # abundances and the stability. for the Generalized Lokta Volterra model

  if (sim_args$model=="GLV"){
    return(EEMtoolbox::summarise_ecosystem_features_GLV(parameters,sim_args))
  } else if (sim_args$model=="Baker"){
    return(EEMtoolbox::summarise_ecosystem_features_Baker(parameters,sim_args))
  } else if (sim_args$model=="Gompertz"){
    return(EEMtoolbox::summarise_ecosystem_features_Gompertz(parameters,sim_args))
  }
}
