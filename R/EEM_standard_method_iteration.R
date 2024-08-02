#' @title Generation of model ensembles using standard EEM sampling method
#' @description
#' Generation of model ensembles using standard EEM sampling method
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @param summ_func function calculating equilibrium points and real parts of the Jacobians eigenvalues to summarise ecosystem features.
#' @param disc_func summary statistic (discrepancy measure).
#' @param sampler sampling function that generates random vectors from the joint prior distribution.
#' @param trans_f transform of prior parameter space to ensure unbounded support for MCMC sampling.
#' @param n_ensemble Number of desired ensemble members. Default to 5000
#' @param n_cores Number of cores desired to be used for sampling. Default set to 1 core (sequential sampling).
#' @return list: sims=number of simulations
#' part_vals=parameter values
#' part_s=discrepancy value
#' prior_sample=prior distribution
#' @export
#' @import parallel
#' @import doParallel
#' @import foreach
EEM_standard_method_iteration <- function(sim_args,
                                summ_func,
                                disc_func,
                                sampler,
                                trans_f,
                                n_ensemble=5000,
                                n_cores=1L){
  # sample prior
  available_cores <- n_cores_function(n_cores)
  cl <- parallel::makeCluster(available_cores[1])
  doParallel::registerDoParallel(cl)
  part_vals <- foreach::foreach(i = 1:n_ensemble, .combine="rbind") %dopar% {
    sampler(sim_args)
  }
  parallel::stopCluster(cl)#stop cluster
  rm(cl)
  part_vals <- unname(part_vals)

  # simulate model
  cl <- parallel::makeCluster(available_cores[1])
  doParallel::registerDoParallel(cl)
  part_sim <- foreach::foreach(i = 1:n_ensemble, .combine = "c") %dopar% {
    summ_func(part_vals[i,], sim_args)
  }

  parallel::stopCluster(cl) #stop cluster
  rm(cl)
  part_sim <- matrix(unname(part_sim), nrow=n_ensemble, byrow = TRUE)

  #simulation
  # evaluate the discrepancy metric
  part_s <- sapply(seq(n_ensemble),
                   function(x, part_sim) disc_func(part_sim[x,]),
                   part_sim = part_sim) #summary

  #save prior sample
  prior_sample <- part_vals
  #track the number of simulations
  sims <- n_ensemble

  # sort the particles
  ix <- order(part_s)
  part_s <- part_s[ix]
  part_vals <- part_vals[ix,]
  part_sim <- part_sim[ix,]

  return(list(sims=sims,
              part_vals=part_vals,
              part_sim=part_sim,
              part_s=part_s,
              prior_sample=prior_sample))

}
