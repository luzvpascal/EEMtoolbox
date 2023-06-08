test_that("ode_solve works", {

  #GLV
  matrix_test <- matrix(c(-1,1,1,-1),ncol = 2)

  for (model_test in c("GLV", "Gompertz")){
    outputs <- EEMtoolbox::EEM(matrix_test,
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


})
