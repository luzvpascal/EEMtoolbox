#' @title Uniform sampling of parameters
#' @description
#' Uniform sampling of parameters between lower and upper bounds defined in sim_args
#'
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @return a vector of sampled parameters
#' @export
#' @import stats

uniform_sampler <- function(sim_args){stats::runif(sim_args$num_params, sim_args$lower,sim_args$upper)} #uniform distribution
