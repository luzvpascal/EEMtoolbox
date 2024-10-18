test_that("multiplication works", {

  interaction_matrix <- matrix(c(-1,1,1,-1),ncol = 2)
  for (model_test in c("GLV", "Gompertz")){
    outputs <- EEMtoolbox::EEM(interaction_matrix,
                               model = model_test,
                               output_matrix = FALSE,
                               output_args = TRUE,
                               output_discrepancy = TRUE,
                               output_prior = TRUE)

    test_values <- EEMtoolbox::reconstruct_matrix_growthrates(
      outputs$part_vals[1,],sim_args = outputs$sim_args)

    expect_true(class(test_values)=="list")
    expect_true(class(test_values$interaction_matrix)[1]=="matrix")
    expect_true(class(test_values$growthrates)=="numeric")
  }

  ## Baker
  model_test = "Bimler-Baker"
  interaction_matrix_list_alphas <- list(matrix(0, ncol=2,nrow=2),
                                         matrix(pmax(interaction_matrix,0), ncol=2))
  interaction_matrix_list_betas <-list(matrix(pmin(interaction_matrix,0), ncol=2),
                                       matrix(0, ncol=2,nrow=2))
  outputs <- EEMtoolbox::EEM(interaction_matrix = list(interaction_matrix_list_alphas,
                                                       interaction_matrix_list_betas),
                             model = model_test,
                             algorithm = "EEM-SMC",
                             output_matrix = FALSE,
                             output_args = TRUE,
                             output_discrepancy = TRUE,
                             output_prior = TRUE)


  test_values <- EEMtoolbox::reconstruct_matrix_growthrates(
    outputs$part_vals[1,],sim_args = outputs$sim_args)


  expect_true(class(test_values)=="list")
  expect_true(class(test_values$interaction_matrix_alphas)[1]=="matrix")
  expect_true(class(test_values$interaction_matrix_betas)[1]=="matrix")
  expect_true(class(test_values$growthrates)=="numeric")

})
