#' @title Generation of model ensembles Sequential Monte Carlo - Approximate Bayesian Computation.
#' @description
#' Generation of model ensembles using Adaptive sequential Monte Carlo for approximate Bayesian computation
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @param summ_func function calculating equilibrium points and real parts of the Jacobians eigenvalues to summarise ecosystem features.
#' @param disc_func summary statistic (discrepancy measure).
#' @param sampler sampling function that generates random vectors from the joint prior distribution.
#' @param trans_f transform of prior parameter space to ensure unbounded support for MCMC sampling.
#' @param trans_finv inverse of trans_f function.
#' @param pdf joint probability density function.
#' @param mcmc_trials number of MCMC steps to try before selecting appropriate number.
#' @param dist_final target discrepancy threshold. Default 0. If zero, p_acc_min is used to determine stopping criteria.
#' @param a tuning parameter for adaptive selection of discrepancy threshold sequence.
#' @param c tuning parameter for choosing the number of MCMC iterations in move step.
#' @param p_acc_min minimum acceptable acceptance rate in the MCMC interations before exit.
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

EEM_SMC_method <- function(sim_args,
                           summ_func,
                           disc_func,
                           sampler,
                           trans_f,
                           trans_finv,
                           pdf,
                           mcmc_trials,
                           dist_final,
                           a,
                           c,
                           p_acc_min,
                           n_ensemble=5000,
                           n_cores=1L){

  # initial prior rejection algorithm: EEM_standard_method
  outputs <- EEMtoolbox::EEM_standard_method_iteration(sim_args,
                                                       summ_func,
                                                       disc_func,
                                                       sampler,
                                                       trans_f,
                                                       n_ensemble,
                                                       n_cores)

  n_sets_correct <- sum(outputs$part_s==0)
  if (n_sets_correct>=n_ensemble){
    print("Exiting SMC ABC search as number of desired ensemble members attained")
    print("modify n_ensemble if necessary")
    return(outputs)
  }
  print(paste("Number of parameter sets found so far:", n_sets_correct, "/", n_ensemble))

  part_vals <- outputs$part_vals # sample prior
  part_sim <- outputs$part_sim # simulate model
  part_s <- outputs$part_s # evaluate the discrepancy metric
  prior_sample <- outputs$prior_sample #prior

  # values for adaptive steps
  num_drop <- floor(n_ensemble*a) #Number of particles dropped each iteration
  num_keep <- n_ensemble-num_drop #Number of particles kept each iteration

  # determine next disprepacy threshold
  dist_max <- part_s[n_ensemble]
  dist_next <- part_s[num_keep]
  print(paste('Current maximum discrepancy', round(dist_max, digits = 2),'now trying for ',
              round(dist_next, digits = 2), 'want to get to ',dist_final))

  # iterate toward target discrepancy
  #track the number of simulations
  sims <- n_ensemble
  # transform the parameters
  part_vals <- trans_f(part_vals,sim_args) #part vals is transformed

  while (dist_max > dist_final){
    # compute the covariance matrix (of particles that remain) required
    # for the Independent MH move step
    cov_matrix <- (2.38^2)*cov(part_vals[seq(num_keep),])/(dim(part_vals)[2])

    ###########
    # resample#
    ###########

    resample <- sample(seq(num_keep), n_ensemble-num_keep, TRUE) #duplicate good ones
    part_vals[seq(num_keep+1,n_ensemble),] <- part_vals[resample,]
    part_s[seq(num_keep+1,n_ensemble)] <- part_s[resample]
    part_sim[seq(num_keep+1,n_ensemble),] <- part_sim[resample,]


    ############
    # move step#
    ############
    available_cores <- n_cores_function(n_cores)
    cl <- parallel::makeCluster(available_cores[1])
    doParallel::registerDoParallel(cl)
    # print("mcmc1")

    #parallel for loop
    mcmc_outcome <- foreach::foreach(i=(num_keep+1):n_ensemble,
                     .packages =c('EEMtoolbox')) %dopar% {
      EEMtoolbox::MCMC(i,
                       sim_args,
                       mcmc_trials,
                       dist_next,
                       part_vals,
                       part_s,
                       part_sim,
                       cov_matrix,
                       summ_func,
                       disc_func,
                       trans_finv,
                       pdf)
    }
    parallel::stopCluster(cl)#stop cluster
    rm(cl)

    part_vals[seq(num_keep+1,n_ensemble),] <- matrix(unlist(lapply(mcmc_outcome, `[[`, 1)),
                             nrow=length(mcmc_outcome), byrow = TRUE)
    part_s[seq(num_keep+1,n_ensemble)] <- (unlist(lapply(mcmc_outcome, `[[`, 2)))
    part_sim[seq(num_keep+1,n_ensemble),] <- matrix(unlist(lapply(mcmc_outcome, `[[`, 3)),
                            nrow=length(mcmc_outcome), byrow = TRUE)
    i_acc <- (unlist(lapply(mcmc_outcome, `[[`, 4)))
    sims_mcmc <- (unlist(lapply(mcmc_outcome, `[[`, 5)))

    ################
    # end move step#
    ################

    #################################################
    # determine number of MCMC iterations to perfrom#
    #################################################
    acc_rate <- sum(i_acc)/(mcmc_trials*(n_ensemble-num_keep))
    mcmc_iters <-   floor(log(c)/log(1-acc_rate)+1)

    ############
    # move step#
    ############
    available_cores <- n_cores_function(n_cores)
    cl <- parallel::makeCluster(available_cores[1])
    doParallel::registerDoParallel(cl)
    # print("mcmc2")
    #run for loop parallel
    mcmc_outcome2 <- foreach::foreach(i=(num_keep+1):n_ensemble,
                  .packages =c('EEMtoolbox')) %dopar% {
      EEMtoolbox::MCMC(i,
                       sim_args,
                       (mcmc_iters-mcmc_trials),
                       dist_next,
                       part_vals,
                       part_s,
                       part_sim,
                       cov_matrix,
                       summ_func,
                       disc_func,
                       trans_finv,
                       pdf)
    }
    parallel::stopCluster(cl) #stop cluster
    rm(cl)

    part_vals[seq(num_keep+1,n_ensemble),] <- matrix(unlist(lapply(mcmc_outcome, `[[`, 1)),
                                                      nrow=length(mcmc_outcome), byrow = TRUE)
    part_s[seq(num_keep+1,n_ensemble)] <- (unlist(lapply(mcmc_outcome, `[[`, 2)))
    part_sim[seq(num_keep+1,n_ensemble),] <- matrix(unlist(lapply(mcmc_outcome, `[[`, 3)),
                                                     nrow=length(mcmc_outcome), byrow = TRUE)
    i_acc <- (unlist(lapply(mcmc_outcome, `[[`, 4)))
    sims_mcmc <- (unlist(lapply(mcmc_outcome, `[[`, 5)))
    ################
    # end move step#
    ################

    num_mcmc_iters <- max(0, mcmc_iters - mcmc_trials) + mcmc_trials
    p_acc <- sum(i_acc)/(num_mcmc_iters*(n_ensemble-num_keep))

    sims <- sims + sum(sims_mcmc)
    mcmc_trials <- ceiling(mcmc_iters/2)

    # compute the next distance and maximum distance
    # sort the particles
    ix <- order(part_s)
    part_s <- part_s[ix]
    part_vals <- part_vals[ix,]
    part_sim <- part_sim[ix,]

    # if most of the particles are under the final target then don't
    # drop them
    if (sum((part_s > dist_final)) < num_drop){
      num_drop <- sum((part_s > dist_final))
      num_keep <- n_ensemble-num_drop
    }

    #calculate distances
    dist_max <- part_s[n_ensemble]
    dist_next <- part_s[num_keep]

    # check to see if we reach desired tolerance at next iteration
    if (dist_next < dist_final){
      dist_next <- dist_final
    }

    n_sets_correct <- sum(part_s==0)
    print(paste("Number of parameter sets found so far:", min(n_sets_correct,n_ensemble)
                "/", n_ensemble))
    print(paste('The next distance is', round(dist_next, digits = 2),
                ' and the maximum distance is ', round(dist_max, digits = 2), ' and the number to drop is ', num_drop))

    print(paste("Number of MCMC runs so far:", sims))
    #if we are not accepting enough particles - give up!
    # print(paste("p_acc", p_acc))
    if (p_acc < p_acc_min){
        print('Getting out as MCMC acceptance rate is below acceptable threshold')
        break
    }
  }

  #transform back
  for (i in seq(n_ensemble)){
      part_vals[i,] <- trans_finv(matrix(part_vals[i,], nrow=1),sim_args)
  }

  return(list(sims=sims,
              part_vals=part_vals,
              part_s=part_s,
              prior_sample=prior_sample))
}
