EEM_standard_method <- function(summ_func,
                           args,
                           disc_func,
                           sampler,
                           n_particles,
                           n_ensemble){
  ## Standard EEM.####
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

  # sort the particles
  ix <- order(part_s)
  part_s <- part_s[ix]
  part_vals <- part_vals[ix,]
  part_sim <- part_sim[ix,]
  
  return(list(sims=sims,
              part_vals=part_vals,
              part_s = part_s,
              prior_sample=prior_sample))
  
}
