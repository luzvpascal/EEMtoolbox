test_that("EEM tests", {
  interaction_matrix <- matrix(c(-1,1,1,-1),ncol = 2)

  ## GLV
  model_test <- "GLV"
  #return parameter values SMC ABC
  outputs <- EEMtoolbox::EEM(interaction_matrix,
                             model = model_test,
                             algorithm = "EEM-SMC",
                             n_ensemble=5000,
                             output_matrix = FALSE,
                             output_args = TRUE,
                             output_discrepancy = TRUE,
                             output_prior = TRUE)

  expect_true(class(outputs$part_vals)[1]=="matrix")
  expect_true(class(outputs$prior_sample)[1]=="matrix")
  expect_true(class(outputs$part_s)[1]=="numeric")
  expect_true(class(outputs$sim_args)[1]=="list")

  interaction_matrix <- matrix(c(-1,1,1,-1),ncol = 2)

  ## GLV
  model_test <- "GLV"
  #return parameter values SMC ABC
  outputs <- EEMtoolbox::EEM(EEMtoolbox::dingo_matrix,
                             model = model_test,
                             algorithm = "EEM-SMC",
                             n_ensemble=5000,
                             output_matrix = FALSE,
                             output_args = TRUE,
                             output_discrepancy = TRUE,
                             output_prior = TRUE)

  expect_true(class(outputs$part_vals)[1]=="matrix")
  expect_true(class(outputs$prior_sample)[1]=="matrix")
  expect_true(class(outputs$part_s)[1]=="numeric")
  expect_true(class(outputs$sim_args)[1]=="list")

  #return parameter values standard EEM
  outputs <- EEMtoolbox::EEM(interaction_matrix,
                             model = model_test,
                             algorithm = "standard-EEM",
                             output_matrix = FALSE,
                             output_args = TRUE,
                             output_discrepancy = TRUE,
                             output_prior = TRUE)

  expect_true(class(outputs$part_vals)[1]=="matrix")
  expect_true(class(outputs$prior_sample)[1]=="matrix")
  expect_true(class(outputs$part_s)[1]=="numeric")
  expect_true(class(outputs$sim_args)[1]=="list")


  ## return matrix SMC ABC
  outputs <- EEMtoolbox::EEM(interaction_matrix,
                             model = model_test,
                             algorithm = "EEM-SMC",
                             output_matrix = TRUE)

  expect_true(class(outputs)[1]=="list")
  expect_true(class(outputs[[1]])[1]=="list")
  expect_true(class(outputs[[1]]$interaction_matrix)[1]=="matrix")
  expect_true(class(outputs[[1]]$growthrates)[1]=="numeric")

  ## return matrix standard EEM
  outputs <- EEMtoolbox::EEM(interaction_matrix,
                             model = model_test,
                             algorithm = "standard-EEM",
                             output_matrix = TRUE)

  expect_true(class(outputs)[1]=="list")
  expect_true(class(outputs[[1]])[1]=="list")
  expect_true(class(outputs[[1]]$interaction_matrix)[1]=="matrix")
  expect_true(class(outputs[[1]]$growthrates)[1]=="numeric")

})
