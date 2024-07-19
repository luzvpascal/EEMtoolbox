#' @title Generation of model ensembles using standard EEM sampling method
#' @description
#' Generation of model ensembles using standard EEM sampling method
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @param summ_func function calculating equilibrium points and real parts of the Jacobians eigenvalues to summarise ecosystem features.
#' @param disc_func summary statistic (discrepancy measure).
#' @param sampler sampling function that generates random vectors from the joint prior distribution.
#' @param trans_f transform of prior parameter space to ensure unbounded support for MCMC sampling.
#' @param n_particles number of particles in the sample.
#' @param n_ensemble Number of desired ensemble members. Default to 5000
#' @param n_cores Number of cores desired to be used for sampling. Default set to 1 core (sequential sampling).
#' @examples
#' library(EEMtoolbox)
#'
#' EEM(dingo_matrix) #automatically loads an example of interaction matrix as dingo_matrix
#' @return list: sims=number of simulations
#' part_vals=parameter values
#' part_s=discrepancy value
#' prior_sample=prior distribution
#' @export
#' @import parallel
#' @import doParallel
#' @import foreach
EEM_standard_method <- function(sim_args,
                                summ_func,
                                disc_func,
                                sampler,
                                trans_f,
                                n_particles,
                                n_ensemble=5000,
                                n_cores=1L){
  # initial prior rejection algorithm: EEM_standard_method
  start <- Sys.time()
  outputs <- EEMtoolbox::EEM_standard_method_iteration(sim_args,
                                                       summ_func,
                                                       disc_func,
                                                       sampler,
                                                       trans_f,
                                                       n_particles,
                                                       n_cores)
  end <- Sys.time()

  if (sum(outputs$part_s==0)>=n_ensemble){
    idx <- which(outputs$part_s==0)[seq(n_ensemble)]
    return(list(sims=outputs$sims,
                part_s = outputs$part_s[idx],
                part_vals = outputs$part_vals[idx,],
                part_sim = outputs$part_sim[idx,],
                prior_sample=outputs$prior_sample))
  } else {
    acceptance_rate <- sum(outputs$part_s==0)/n_particles
    time.taken <- round(end - start,2)
    units_time <- units(end - start)
    estimated_iterations <- n_ensemble/(acceptance_rate*n_particles)

    print(paste('Estimated acceptance rate:', acceptance_rate))
    if (acceptance_rate == 0){
      print(paste('Time to generate ', n_ensemble,
                  ' ensembles might be very long with the standard search method'))
      print("We recommend using the SMC-ABC algorithm to sample ensembles.")
      print("See documentation of the function EEM for the input 'algorithm'")
    } else {
      print(paste('Estimated time to generate ', n_ensemble,
                  ' ensembles: ', round(estimated_iterations*time.taken,2),
                  units_time, sep=""))
    }

    idx <- which(outputs$part_s==0)
    sims <- outputs$sims
    part_s <- outputs$part_s[idx]
    part_vals <- outputs$part_vals[idx,]
    part_sim <- outputs$part_sim[idx,]
    prior_sample <- outputs$prior_sample

    while(nrow(part_vals)<=n_ensemble){
      #run standard iteration until generating at least n_ensemble ensembles
      outputs <- EEMtoolbox::EEM_standard_method_iteration(sim_args,
                                                           summ_func,
                                                           disc_func,
                                                           sampler,
                                                           trans_f,
                                                           n_particles,
                                                           n_cores)
      idx <- which(outputs$part_s==0)
      sims <- sims+outputs$sims
      part_s <- c(part_s, outputs$part_s[idx])
      part_vals <- rbind(part_vals,outputs$part_vals[idx,])
      part_sim <- rbind(part_sim,outputs$part_sim[idx,])
      prior_sample <- rbind(prior_sample,outputs$prior_sample)
    }
    return(list(sims=sims,
                part_vals=part_vals,
                part_sim=part_sim,
                part_s=part_s,
                prior_sample=prior_sample))
  }
}
