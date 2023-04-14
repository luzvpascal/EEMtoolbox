#' @title Generation of model ensembles
#' @description
#' Generation of model ensembles based on generalized Lotka Volterra, and the other two model structures, generating algorithms include Approximate Bayesian Computation methods and standard ensemble ecosystem modelling (Baker et al., 2017)
#' @param interaction_matrix interaction signs matrix, can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second
#' @param bounds_growth_rate vector of 2 elements containing lower and upper bounds for growth rates. Default c(-5,5)
#' @param n_ensemble Number of desired ensemble members. Default to 10
#' @param model model representing species interactions. Default "GLV" (Generalized Lokta Voltera). options include "Baker", "Gompertz" and "customized"
#' @param algorithm algorithm used for sampling. Default "SMC-ABC" (Vollert et al., 2023) options include "standard EEM"
#' @param summ_func function calculating equilibrium points and real parts of the Jacobians eigenvalues to summarise ecosystem features. Default =summarise_ecosystem_features_GLV. Options include summarise_ecosystem_features_Baker (automatically chosen if model="Baker") and summarise_ecosystem_features_Gompertz, (automatically chosen if model="Gompertz"). Needs to be defined if model="customized" chosen.
#' @param disc_func summary statistic (discrepancy measure). Default discrepancy_continuous_sum
#' @param sampler sampling function that generates random vectors from the joint prior distribution. Default EEMtoolbox::sampler function (uniform)
#' @param trans_f transform of prior parameter space to ensure unbounded support for MCMC sampling. Default EEMtoolbox::uniform_transform
#' @param trans_finv inverse of trans_f function. Default EEMtoolbox::uniform_transform_inverse
#' @param pdf joint probability density function. Default EEMtoolbox::uniform_pdf_transformed
#' @param n_particles number of particles in the sample. Default 10000
#' @param mcmc_trials number of MCMC steps to try before selecting appropriate number. Default 10
#' @param dist_final target discrepancy threshold. Default 0. If zero, p_acc_min is used to determine stopping criteria.
#' @param a tuning parameter for adaptive selection of discrepancy threshold sequence. Defalut 0.6
#' @param c tuning parameter for choosing the number of MCMC iterations in move step. Default 0.01
#' @param p_acc_min minimum acceptable acceptance rate in the MCMC interations before exit. Default 0.0001
#' @param output_prior logical. If set to TRUE, algorithm returns prior distributions of parameters ensemble of parameters. Default FALSE
#' @param output_args logical. If set to TRUE, algorithm returns output from EEMtoolbox::args_function for this problem
#' @param output_matrix logical. If set to TRUE, algorithm returns interaction matrix and growthrates
#' @examples
#' library(EEMtoolbox)
#' EEM(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' @return list: part_vals: ensemble of parameters, marginal distributions
#' @export
EEM <- function(interaction_matrix,
                bounds_growth_rate=c(-5,5),
                n_ensemble=10,
                model="GLV",
                algorithm="SMC-ABC",
                summ_func=EEMtoolbox::summarise_ecosystem_features_GLV,
                disc_func=EEMtoolbox::discrepancy_continuous_sum,
                sampler=EEMtoolbox::uniform_sampler,
                trans_f=EEMtoolbox::uniform_transform,
                trans_finv=EEMtoolbox::uniform_transform_inverse,
                pdf=EEMtoolbox::uniform_pdf_transformed,
                n_particles=10000,
                mcmc_trials=10,
                dist_final=0,
                a=0.6,
                c=0.01,
                p_acc_min=0.0001,
                output_prior=FALSE,
                output_args=FALSE,
                output_matrix=FALSE
                ){
  # TESTS if inputs are correct ###########
  #interaction_matrix tests
  stopifnot(class(interaction_matrix)[1]=="matrix"|class(interaction_matrix)=="list")

  if (class(interaction_matrix)[1]=="matrix"){
    stopifnot(nrow(interaction_matrix)==ncol(interaction_matrix))
  } else {
    stopifnot(nrow(interaction_matrix[[1]])==ncol(interaction_matrix[[1]]))
    stopifnot(nrow(interaction_matrix[[2]])==ncol(interaction_matrix[[2]]))
    stopifnot(nrow(interaction_matrix[[1]])==nrow(interaction_matrix[[2]]))
  }

  print("need to test values of interaction_matrix: positive and negative?")
  #n_ensemble tests
  stopifnot(is.numeric(n_ensemble),
            (n_ensemble)>0)
  #model tests
  stopifnot(((model=="GLV")|(model=="Baker")|(model=="Gompertz")|(model=="customized")))
  #algorithm
  stopifnot(((algorithm=="SMC-ABC")|(algorithm=="standard EEM")))
  #summ_func
  stopifnot(class(summ_func)=="function")
  #disc_func
  stopifnot(class(disc_func)=="function")
  #sampler
  stopifnot(class(sampler)=="function")
  #trans_f
  stopifnot(class(trans_f)=="function")
  #trans_finv
  stopifnot(class(trans_finv)=="function")
  #pdf
  stopifnot(class(pdf)=="function")
  #n_particles
  stopifnot(is.numeric(n_particles),
            (n_particles)>n_ensemble)
  #mcmc_trials
  stopifnot(is.numeric(mcmc_trials),
            (n_particles)>mcmc_trials)
  #dist_final
  stopifnot(is.numeric(dist_final),
            (dist_final)>=0)
  #a
  stopifnot(is.numeric(a),
            (a)>=0)
  #c
  stopifnot(is.numeric(c),
            (c)>=0)
  #p_acc_min
  stopifnot(is.numeric(p_acc_min),
            (p_acc_min)>=0)
  #output_prior
  stopifnot(class(output_prior)=="logical")
  #output_args
  stopifnot(class(output_args)=="logical")

  # SETTING values to parameters ####
  if (model == "Baker"){
    summ_func <- EEMtoolbox::summarise_ecosystem_features_Baker
  } else if (model == "Gompertz"){
    summ_func <- EEMtoolbox::summarise_ecosystem_features_Gompertz
    bounds_growth_rate <- c(1.1,1.3)
    print("r values bounded between 1 and 2?")
  }

  # Defining special arguments ####
  sim_args <- EEMtoolbox::args_function(interaction_matrix,
                        bounds_growth_rate,
                        model=model)

  ## RUNNING search algorithms####
  if (algorithm == "SMC-ABC"){
    print('Begin SMC-ABC search method')
    outputs <- EEMtoolbox::EEM_SMC_method(sim_args,
                                          summ_func,
                                          disc_func,
                                          sampler,
                                          trans_f,
                                          trans_finv,
                                          pdf,
                                          n_particles,
                                          mcmc_trials,
                                          dist_final,
                                          a,
                                          c,
                                          p_acc_min)
  } else if ((algorithm=="standard EEM")){
    print('Begin standard search method')
    outputs <- EEMtoolbox::EEM_standard_method(sim_args,
                                               summ_func,
                                               disc_func,
                                               sampler,
                                               n_particles)
  }
  if (output_matrix){
    output_function <- apply(outputs$part_vals, 1, EEMtoolbox::reconstruct_matrix_growthrates)
    return(output_function)
  } else {
    output_function <- list()
    output_function$part_vals <- outputs$part_vals
    if (output_prior){
      output_function$prior_sample <- outputs$prior_sample
    }
    if (output_args){
      output_function$args <- args
    }
    return(output_function)
  }

}
