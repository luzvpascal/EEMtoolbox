MCMC <- function(i, mcmc_trials, part_vals, part_s, part_sim, 
                 cov_matrix, trans_finv, pdf, summ_func, disc_func, num_keep, dist_next, args){
  # documentation
  sims_mcmc <-  0
  i_acc <- 0
  for (r in seq(mcmc_trials)){
    # Gaussian random walk
    part_vals_prop <- MASS::mvrnorm(n=1,mu=part_vals[i,],Sigma=cov_matrix)
    
    # Transform back to calculate prior probs and discrepancy
    prop <- trans_finv(matrix(part_vals_prop, nrow=1), args)
    
    # Calculate prior probabilities
    prior_curr <- pdf(matrix(part_vals[i,], nrow=1))
    prior_prop <- pdf(matrix(part_vals_prop, nrow=1))
    
    # early rejection (assumes symmetric proposals)
    if(!((is.nan(prior_prop/prior_curr))| (runif(1) > prior_prop/prior_curr))){
      break
    }
    
    #find proposal discrepancy
    part_sim_prop <- summ_func(prop, args)
    dist_prop <- disc_func(part_sim_prop)
    
    sims_mcmc=sims_mcmc+1
    # ABC part of the acceptance probability
    
    #Accept a particle if it is within the target distance.
    if (dist_prop <= dist_next) {
      # then the metropolis-hastings ratio is satisfied
      part_vals[i,] <- part_vals_prop 
      part_s[i] <- dist_prop 
      part_sim[i,] <- part_sim_prop
      i_acc <- i_acc + 1
    }
  }
  return(list(part_vals=part_vals[i,],
              part_s=part_s[i],
              part_sim=part_sim[i,],
              i_acc=i_acc,
              sims_mcmc=sims_mcmc))
}
