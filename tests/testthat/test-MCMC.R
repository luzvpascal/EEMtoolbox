test_that("MCMC works", {
  # sample prior
  n_particles <- 10
  num_keep <- 5
  mcmc_trials <- 100

  model <- "GLV"
  interaction_matrix <- matrix(c(-1, -1, -1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0,
                                  1, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, 0, -1, 0, 0,
                                  0, 0, 1, 1, -1, 1, 1, 0, 0, 0, 0, 1, 1, -1, 0, 1, 0,
                                  0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, 1, 0, -1), nrow=8)

  sim_args <- EEMtoolbox::args_function(interaction_matrix,
                                        upper_bounds_growth_rate = 5,
                                        lower_bounds_growth_rate = 0,
                                        model=model)

  part_vals <- t(sapply(seq(n_particles),
                        function(x, sim_args) EEMtoolbox::uniform_sampler(sim_args), sim_args=sim_args))
  # simulate model
  part_sim <- list()
  for (i in seq(n_particles)){
    part_sim[[i]] <- EEMtoolbox::summarise_ecosystem_features(part_vals[i,], sim_args)
  }
  part_sim <- matrix(unlist(part_sim), nrow=n_particles, byrow = TRUE)

  #simulation
  # evaluate the discrepancy metric
  part_s <- sapply(seq(n_particles),
                   function(x, part_sim) EEMtoolbox::discrepancy_continuous_sum(part_sim[x,]),
                   part_sim = part_sim) #summary

  # sort the particles
  ix <- order(part_s)
  part_s <- part_s[ix]
  part_vals <- part_vals[ix,]
  part_sim <- part_sim[ix,]

  ##MCMC
  cov_matrix <- (2.38^2)*cov(part_vals[seq(num_keep),])/(dim(part_vals)[2])

  output <- EEMtoolbox::MCMC(num_keep+1,
                   sim_args,
                   mcmc_trials,
                   dist_next=part_s[num_keep],
                   part_vals,
                   part_s,
                   part_sim,
                   cov_matrix,
                   summ_func=EEMtoolbox::summarise_ecosystem_features,
                   disc_func=EEMtoolbox::discrepancy_continuous_sum,
                   trans_finv = EEMtoolbox::uniform_transform_inverse,
                   pdf=EEMtoolbox::uniform_pdf_transformed)

  #' part_s: discrepancy measure for particle i
  #' part_sim: summary of ecosystem features for particle i
  #' i_acc: number of times particle movement accepted
  #' sims_mcmc: number of successful walks
  expect_true(output$part_s <= part_s[num_keep+1])
})
