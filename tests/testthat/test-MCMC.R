test_that("MCMC works", {
  # sample prior
  n_particles <- 10
  model <- "GLV"
  interaction_matrix <-
  sim_args <- EEMtoolbox::args_function(interaction_matrix,
                                        bounds_growth_rate,
                                        model=model)

  part_vals <- t(sapply(seq(n_particles),
                        function(x, sim_args) sampler(sim_args), sim_args=sim_args))
  # simulate model
  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores[1]-2) #not to overload your computer
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

  EEMtoolbox::MCMC()
})
