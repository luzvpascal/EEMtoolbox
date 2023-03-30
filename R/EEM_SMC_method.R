EEM_SMC_method <- function(summ_func,
                           args,
                           disc_func,
                           sampler,
                           trans_f,
                           trans_finv,
                           pdf,
                           n_particles,
                           dist_final,
                           a,
                           c,
                           p_acc_min,
                           mcmc_trials){
  ## Adaptive sequential Monte Carlo for approximate Bayesian computation.####
  #  
  # Parameters: 
  #  sim_func  - user defined function that simulates the model of interest given
  #              a vector of parameter values.
  #	 args      - any additional arguments required
  #  dist_func - distance function/discrepancy metric for comparing simulated and 
  #               true data.
  #  sampler    - a sampling function that generates random vectors from
  #                            the joint prior distribution
  #  pdf        - the joint probability density function
  #  trans_f    - transform of prior parameter space to ensure 
  #                            unbounded support for MCMC sampling.
  #  trans_finv - inverse of transform.
  #  n_particles- number of particles for SMC sampler.
  #  dist_final - target discrepancy threshold. If zero, then p_acc_min is used to
  #               determine stopping criteria.
  #  a          - tuning parameter for adaptive selection of discrepancy threshold 
  #               sequence. 
  #  c          - tuning parameter for choosing the number of MCMC iterations in 
  #               move step.
  #  p_acc_min  - minimum acceptable acceptance rate in the MCMC interations. If
  #               zero the dist_final is used to determine stopping criteria. 
  #
  # Returns:
  #    part_val  - parameter values for each particle.
  #    part_sim  - summary statistics for each particle.
  #    part_c    - discrepancy metric value for each particle.
  #    sim       - total number of model simulations performed.
  #    dist_t    - smallest discrepacy threshold reached.
  #    p_acc_min - smallest MCMC acceptance rate reached.
  #
  ## Authors:
  #     Christopher Drovandi (c.drovandi@qut.edu.au)
  #           School of Mathematical Sciences
  #           Science Faculty 
  #           Queensland University of Technology
  #
  #     David J. Warne (david.warne@qut.edu.au)
  #           School of Mathematical Sciences
  #           Science Faculty
  #           Queensland University of Technology
  ###########################################################################
  
  # values for adaptive steps
  num_drop <- floor(n_particles*a) #Number of particles dropped each iteration
  num_keep <- n_particles-num_drop #Number of particles kept each iteration
  
  # initial prior rejection algorithm 
  # sample prior
  part_vals <- t(sapply(seq(n_particles),
         function(x, args) sampler(args), args=args))
  # simulate model
  cores=detectCores()
  cl <- makeCluster(cores[1]-2) #not to overload your computer
  registerDoParallel(cl)
  part_sim <- foreach::foreach(i = 1:n_particles) %dopar% {
    summ_func(part_vals[i,], args)
  }
  #stop cluster
  stopCluster(cl)
  part_sim <- matrix(unlist(part_sim), nrow=n_particles, byrow = TRUE)
  #simulation
  # evaluate the discrepancy metric
  part_s <- sapply(seq(n_particles),
                   function(x, part_sim) disc_func(part_sim[x,]),
                   part_sim = part_sim) #summary
  
  print(paste('Estimated acceptance rate <-',sum(part_s==0)/n_particles))
  
  #save prior sample
  prior_sample <- part_vals
  #track the number of simulations
  sims <- n_particles
  # transform the parameters
  part_vals <- trans_f(part_vals,args)        #part vals is transformed
  
  # sort the particles
  ix <- order(part_s)
  part_s <- part_s[ix]
  part_vals <- part_vals[ix,]
  part_sim <- part_sim[ix,]
  
  # determine next disprepacy threshold
  dist_max <- part_s[n_particles]
  dist_next <- part_s[num_keep]
  print(paste('Current maximum discrepancy',dist_max,'now trying for ', dist_next, 'want to get to ',dist_final))
  
  # interate toward target discrepancy
  while (dist_max > dist_final){
    # compute the covariance matrix (of particles that remain) required 
    # for the Independent MH move step
    cov_matrix <- (2.38^2)*cov(part_vals[seq(num_keep),])/(dim(part_vals)[2])
    
    ###########
    # resample#
    ###########
    
    resample <- sample(seq(num_keep), n_particles-num_keep, TRUE) #duplicate good ones
    part_vals[seq(num_keep+1,n_particles),] <- part_vals[resample,]
    part_s[seq(num_keep+1,n_particles)] <- part_s[resample]
    part_sim[seq(num_keep+1,n_particles),] <- part_sim[resample,]
    
    
    ############
    # move step#
    ############
    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1]-2) #not to overload your computer
    registerDoParallel(cl)
    print("mcmc1")
    mcmc_outcome <- foreach::foreach(i=(num_keep+1):n_particles, 
                     .export=c('MCMC')) %dopar% {
      MCMC(i, mcmc_trials, part_vals, part_s, part_sim, 
           cov_matrix, trans_finv, pdf, summ_func, disc_func, num_keep, dist_next, args)
    }
    #stop cluster
    stopCluster(cl)
    part_vals[seq(num_keep+1,n_particles),] <- matrix(unlist(lapply(mcmc_outcome, `[[`, 1)), 
                             nrow=length(mcmc_outcome), byrow = TRUE)
    part_s[seq(num_keep+1,n_particles)] <- (unlist(lapply(mcmc_outcome, `[[`, 2)))
    part_sim[seq(num_keep+1,n_particles),] <- matrix(unlist(lapply(mcmc_outcome, `[[`, 3)), 
                            nrow=length(mcmc_outcome), byrow = TRUE)
    i_acc <- (unlist(lapply(mcmc_outcome, `[[`, 4)))
    sims_mcmc <- (unlist(lapply(mcmc_outcome, `[[`, 5)))
    
    ################
    # end move step#
    ################
    
    #################################################
    # determine number of MCMC iterations to perfrom#
    #################################################
    acc_rate <- sum(i_acc)/(mcmc_trials*(n_particles-num_keep))
    mcmc_iters <-   floor(log(c)/log(1-acc_rate)+1)
    
    ############
    # move step#
    ############
    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1]-2) #not to overload your computer
    registerDoParallel(cl)
    print("mcmc2")
    mcmc_outcome2 <- foreach::foreach(i=(num_keep+1):n_particles, 
                  .export=c('MCMC')) %dopar% {
                    #, 'uniform_transform_inverse','uniform_pdf_transformed', 'summ_func','disc_func'
      MCMC(i, (mcmc_iters-mcmc_trials), part_vals, part_s, part_sim, 
           cov_matrix, trans_finv, pdf, summ_func, disc_func, num_keep, dist_next, args)
    }
    #stop cluster
    stopCluster(cl)
    part_vals[seq(num_keep+1,n_particles),] <- matrix(unlist(lapply(mcmc_outcome, `[[`, 1)), 
                                                      nrow=length(mcmc_outcome), byrow = TRUE)
    part_s[seq(num_keep+1,n_particles)] <- (unlist(lapply(mcmc_outcome, `[[`, 2)))
    part_sim[seq(num_keep+1,n_particles),] <- matrix(unlist(lapply(mcmc_outcome, `[[`, 3)), 
                                                     nrow=length(mcmc_outcome), byrow = TRUE)
    i_acc <- (unlist(lapply(mcmc_outcome, `[[`, 4)))
    sims_mcmc <- (unlist(lapply(mcmc_outcome, `[[`, 5)))
    ################
    # end move step#
    ################
    
    num_mcmc_iters <- max(0, mcmc_iters - mcmc_trials) + mcmc_trials
    p_acc <- sum(i_acc)/(num_mcmc_iters*(n_particles-num_keep))
    
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
      num_keep <- n_particles-num_drop
    }
   
    #calculate distances
    dist_max <- part_s[n_particles]
    dist_next <- part_s[num_keep]
    
    # check to see if we reach desired tolerance at next iteration
    if (dist_next < dist_final){
      dist_next <- dist_final
    }
    
    print(paste('The next distance is', dist_next, ' and the maximum distance is ', dist_max, ' and the number to drop is ', num_drop))
    
    #if we are not accepting enough particles - give up!
    print(paste("p_acc", p_acc))
    if (p_acc < p_acc_min){
        print('Getting out as MCMC acceptance rate is below acceptable threshold')
        break
    }
  }
  
  #transform back
  for (i in seq(n_particles)){
      part_vals[i,] <- trans_finv(matrix(part_vals[i,], nrow=1),args)
  }
  
  return(list(sims=sims, 
              part_vals=part_vals, 
              part_s=part_s,
              prior_sample=prior_sample))
}
