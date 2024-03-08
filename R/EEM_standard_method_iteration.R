#' @title Generation of model ensembles using standard EEM sampling method
#' @description
#' Generation of model ensembles using standard EEM sampling method
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @param summ_func function calculating equilibrium points and real parts of the Jacobians eigenvalues to summarise ecosystem features.
#' @param disc_func summary statistic (discrepancy measure).
#' @param sampler sampling function that generates random vectors from the joint prior distribution.
#' @param trans_f transform of prior parameter space to ensure unbounded support for MCMC sampling.
#' @param n_particles number of particles in the sample.
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
                                n_particles){
  # sample prior
  part_vals <- t(sapply(seq(n_particles),
                        function(x, sim_args) sampler(sim_args), sim_args=sim_args))
  # simulate model
  cores <- parallel::detectCores()
  if (cores[1] >= 2){
    cl <- parallel::makeCluster(cores[1]-1) #not to overload your computer
  } else {
    cl <- parallel::makeCluster(1) #not to overload your computer
  }
  doParallel::registerDoParallel(cl)
  part_sim <- foreach::foreach(i = 1:n_particles) %dopar% {
    summ_func(part_vals[i,], sim_args)
  }
  #stop cluster
  parallel::stopCluster(cl)
  rm(cl)
  part_sim <- matrix(unlist(part_sim), nrow=n_particles, byrow = TRUE)
  #simulation
  # evaluate the discrepancy metric
  part_s <- sapply(seq(n_particles),
                   function(x, part_sim) disc_func(part_sim[x,]),
                   part_sim = part_sim) #summary

  #save prior sample
  prior_sample <- part_vals
  #track the number of simulations
  sims <- n_particles

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
