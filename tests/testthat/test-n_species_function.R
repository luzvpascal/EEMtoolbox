test_that("n_species_function works", {

  interaction_matrix <- matrix(c(-1,1,1,-1),ncol=2)

  ## GLV
  n_species_GLV <- EEMtoolbox::n_species_function(interaction_matrix,
                                      model="GLV")
  expect_equal(n_species_GLV, 2)

  ## Gompertz
  n_species_gompertz <- EEMtoolbox::n_species_function(interaction_matrix,
                                      model="Gompertz")
  expect_equal(n_species_gompertz, 2)

  ## Baker
  interaction_matrix_list_alphas <- list(matrix(pmin(interaction_matrix,0), ncol=2),
                                         matrix(pmax(interaction_matrix,0), ncol=2))
  interaction_matrix_list_betas <- interaction_matrix_list_alphas

  interaction_matrix_list <- list(interaction_matrix_list_alphas,
                                  interaction_matrix_list_betas)

  n_species_baker <- EEMtoolbox::n_species_function(interaction_matrix_list,
                                           model="Bimler-Baker")
  expect_equal(n_species_baker, 2)
  n_species_baker <- EEMtoolbox::n_species_function(interaction_matrix_list_alphas,
                                           model="Bimler-Baker")
  expect_equal(n_species_baker, 2)

})
