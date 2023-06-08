test_that("ode_solve works", {

  #GLV and Gompertz
  interaction_matrix <- matrix(c(-1,1,1,-1),ncol = 2)

  for (model_test in c("GLV", "Gompertz")){
    outputs <- EEMtoolbox::EEM(interaction_matrix,
                               model = model_test,
                               n_ensemble=1,
                               output_matrix = FALSE,
                               output_args = TRUE,
                               output_discrepancy = TRUE,
                               output_prior = TRUE)

    initcond <- EEMtoolbox::summarise_ecosystem_features(
      parameters = outputs$part_vals[1,],
      sim_args = outputs$sim_args)[seq(2)]

    test_values <- EEMtoolbox::reconstruct_matrix_growthrates(
      outputs$part_vals[1,],sim_args = outputs$sim_args)

    output_pred <- EEMtoolbox::ode_solve(
      interaction_matrix=test_values$interaction_matrix,
      growth_rate=test_values$growthrates,
      t_window = c(0,10),
      model = model_test,
      initial_condition =initcond-runif(2,min=-0.5,max=0.5)*initcond
    )

    expect_true(max(abs((initcond-output_pred$abundances[nrow(output_pred$abundances),])/initcond))<= 0.1)
  }


  #Baker
  model_test = "Baker"
  interaction_matrix_list_alphas <- list(matrix(pmin(interaction_matrix,0), ncol=2),
                                         matrix(pmax(interaction_matrix,0), ncol=2))
  interaction_matrix_list_betas <- interaction_matrix_list_alphas

  outputs <- EEMtoolbox::EEM(interaction_matrix = list(interaction_matrix_list_alphas,
                                                       interaction_matrix_list_betas),
                             model = model_test,
                             n_ensemble=1,
                             output_matrix = FALSE,
                             output_args = TRUE,
                             output_discrepancy = TRUE,
                             output_prior = TRUE)

  initcond <- EEMtoolbox::summarise_ecosystem_features(
    parameters = outputs$part_vals[1,],
    sim_args = outputs$sim_args)[seq(2)]

  test_values <- EEMtoolbox::reconstruct_matrix_growthrates(
    outputs$part_vals[1,],sim_args = outputs$sim_args)

  output_pred <- EEMtoolbox::ode_solve(
    interaction_matrix=list(test_values$interaction_matrix_alphas,
                            test_values$interaction_matrix_betas),
    growth_rate=test_values$growthrates,
    t_window = c(0,10),
    model = model_test,
    initial_condition =initcond-runif(2,min=-0.5,max=0.5)*initcond
  )

  expect_true(max(abs((initcond-output_pred$abundances[nrow(output_pred$abundances),])/initcond))<= 0.1)
})
