test_that("summarise ecosystem features works", {
  interaction_matrix <- matrix(c(-1,1,1,-1),ncol = 2)

  ## GLV
  model_test <- "GLV"
  outputs <- EEMtoolbox::EEM(interaction_matrix,
                             model = model_test,
                             n_ensemble=10,
                             output_matrix = FALSE,
                             output_args = TRUE,
                             output_discrepancy = TRUE,
                             output_prior = TRUE)

  summary_values <- EEMtoolbox::summarise_ecosystem_features(
    outputs$part_vals[1,],sim_args = outputs$sim_args)

  summary_values2 <- EEMtoolbox::summarise_ecosystem_features_GLV(
    outputs$part_vals[1,],sim_args = outputs$sim_args)

  expect_true(unique(summary_values==summary_values2))

  ## Gompertz
  model_test <- "Gompertz"
  outputs <- EEMtoolbox::EEM(interaction_matrix,
                             model = model_test,
                             n_ensemble=10,
                             output_matrix = FALSE,
                             output_args = TRUE,
                             output_discrepancy = TRUE,
                             output_prior = TRUE)

  summary_values <- EEMtoolbox::summarise_ecosystem_features(
    outputs$part_vals[1,],sim_args = outputs$sim_args)

  summary_values2 <- EEMtoolbox::summarise_ecosystem_features_Gompertz(
    outputs$part_vals[1,],sim_args = outputs$sim_args)

  expect_true(unique(summary_values==summary_values2))

  ## Baker
  model_test = "Baker"
  interaction_matrix_list_alphas <- list(matrix(pmin(interaction_matrix,0), ncol=2),
                                         matrix(pmax(interaction_matrix,0), ncol=2))
  interaction_matrix_list_betas <- interaction_matrix_list_alphas

  outputs <- EEMtoolbox::EEM(interaction_matrix = list(interaction_matrix_list_alphas,
                                                       interaction_matrix_list_betas),
                             model = model_test,
                             n_ensemble=10,
                             output_matrix = FALSE,
                             output_args = TRUE,
                             output_discrepancy = TRUE,
                             output_prior = TRUE)

  summary_values <- EEMtoolbox::summarise_ecosystem_features(
    outputs$part_vals[1,],sim_args = outputs$sim_args)

  summary_values2 <- EEMtoolbox::summarise_ecosystem_features_Baker(
    outputs$part_vals[1,],sim_args = outputs$sim_args)

  expect_true(unique(summary_values==summary_values2))

})
