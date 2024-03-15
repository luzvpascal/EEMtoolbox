test_that("Input interaction matrix or as list works", {
  dingo_matrix <- matrix(c(-1, -1, -1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0,
                           1, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, 0, -1, 0, 0,
                           0, 0, 1, 1, -1, 1, 1, 0, 0, 0, 0, 1, 1, -1, 0, 1, 0,
                           0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, 1, 0, -1), nrow=8)


  #expect same arguments for different types of inputs
  mat_list <- list(pmin(dingo_matrix,0), pmax(dingo_matrix,0))

  args_sim_mat <- EEMtoolbox::args_function(dingo_matrix,
                                            upper_bounds_growth_rate = 5,
                                            lower_bounds_growth_rate = -5,
                                            upper_interaction_strength = 1,
                                            model="GLV")
  args_sim_list <- EEMtoolbox::args_function(mat_list,
                                             upper_bounds_growth_rate = 5,
                                             lower_bounds_growth_rate = -5,
                                             upper_interaction_strength = 1,
                                             model="GLV")

  testthat::expect_equal(args_sim_mat$model,
                         args_sim_list$model)

  testthat::expect_equal(args_sim_mat$n_species,
                         args_sim_list$n_species)

  testthat::expect_equal(args_sim_mat$skip_parameters,
                         args_sim_list$skip_parameters)

  testthat::expect_equal(args_sim_mat$lower,
                         args_sim_list$lower)

  testthat::expect_equal(args_sim_mat$upper,
                         args_sim_list$upper)

  ##

  args_sim_mat_2 <- EEMtoolbox::args_function(dingo_matrix,
                                            upper_bounds_growth_rate = rep(5, nrow(dingo_matrix)),
                                            lower_bounds_growth_rate = rep(0, nrow(dingo_matrix)),
                                            upper_interaction_strength = 1,
                                            model="GLV")
  args_sim_list_2 <- EEMtoolbox::args_function(mat_list,
                                             upper_bounds_growth_rate = rep(5, nrow(dingo_matrix)),
                                             lower_bounds_growth_rate = rep(0, nrow(dingo_matrix)),
                                             upper_interaction_strength = 1,
                                             model="GLV")

  testthat::expect_equal(args_sim_mat_2$model,
                         args_sim_list_2$model)

  testthat::expect_equal(args_sim_mat_2$n_species,
                         args_sim_list_2$n_species)

  testthat::expect_equal(args_sim_mat_2$skip_parameters,
                         args_sim_list_2$skip_parameters)

  testthat::expect_equal(args_sim_mat_2$lower,
                         args_sim_list_2$lower)

  testthat::expect_equal(args_sim_mat_2$upper,
                         args_sim_list_2$upper)
})
