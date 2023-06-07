test_that("get_nonzero_parameters works", {

  ## GLV ##
  interaction_matrix <- matrix(c(-1,1,1,-1),ncol=2)
  non_zero <- EEMtoolbox::get_nonzero_parameters(interaction_matrix,
                                                 n_species = ncol(interaction_matrix),
                                                 model="GLV")

  interaction_matrix_list <- list(matrix(pmin(interaction_matrix,0), ncol=2),
                                  matrix(pmax(interaction_matrix,0), ncol=2))
  non_zero_list <- EEMtoolbox::get_nonzero_parameters(interaction_matrix_list,
                                                      n_species = ncol(interaction_matrix),
                                                      model="GLV")
  expect_equal(non_zero$keep_parameters, non_zero_list$keep_parameters)
  expect_equal(non_zero$lower_interaction_bound, non_zero_list$lower_interaction_bound)
  expect_equal(non_zero$upper_interaction_bound, non_zero_list$upper_interaction_bound)

  ## Gompertz ##
  interaction_matrix <- matrix(c(-1,1,1,-1),ncol=2)
  non_zero <- EEMtoolbox::get_nonzero_parameters(interaction_matrix,
                                                 n_species = ncol(interaction_matrix),
                                                 model="Gompertz")

  interaction_matrix_list <- list(matrix(pmin(interaction_matrix,0), ncol=2),
                                  matrix(pmax(interaction_matrix,0), ncol=2))
  non_zero_list <- EEMtoolbox::get_nonzero_parameters(interaction_matrix_list,
                                                      n_species = ncol(interaction_matrix),
                                                      model="Gompertz")
  expect_equal(non_zero$keep_parameters, non_zero_list$keep_parameters)
  expect_equal(non_zero$lower_interaction_bound, non_zero_list$lower_interaction_bound)
  expect_equal(non_zero$upper_interaction_bound, non_zero_list$upper_interaction_bound)

  ## Baker ##
  interaction_matrix_list_alphas <- list(matrix(pmin(interaction_matrix,0), ncol=2),
                                  matrix(pmax(interaction_matrix,0), ncol=2))
  interaction_matrix_list_betas <- interaction_matrix_list_alphas

  non_zero_list <- EEMtoolbox::get_nonzero_parameters(
    list(interaction_matrix_list_alphas, interaction_matrix_list_betas),
    n_species = ncol(interaction_matrix),
    model="Baker")

  expect_equal(non_zero_list$keep_parameters_alphas,
               non_zero_list$keep_parameters_betas)

  expect_equal(non_zero_list$lower_interaction_bound_alphas,
               non_zero_list$lower_interaction_bound_betas)

  expect_equal(non_zero_list$upper_interaction_bound_alphas,
               non_zero_list$upper_interaction_bound_betas)

})
