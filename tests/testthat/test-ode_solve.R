test_that("ode_solve works", {

  #GLV
  model_test <- "GLV"
  initcond <- c(67.82167, 11.63940) #artificial test
  test_values <- list()
  test_values$interaction_matrix <- matrix(c(-0.06967885, 0.1433530, 0.13648650, -0.7074722), ncol = 2, byrow = TRUE)
  test_values$growthrates <-c( 3.057193, -1.022190)

  output_pred <- EEMtoolbox::ode_solve(
    interaction_matrix=test_values$interaction_matrix,
    growth_rate=test_values$growthrates,
    t_window = c(0,10),
    model = model_test,
    initial_condition =initcond-runif(2,min=-0.5,max=0.5)*initcond
  )
  expect_true(max(abs((initcond-output_pred$abundances[nrow(output_pred$abundances),])/initcond))<= 0.1)

  #Gompertz
  model_test <- "Gompertz"
  initcond <- c( 0.187454942, 0.003951098) #artificial test
  test_values <- list()
  test_values$interaction_matrix <- matrix(c( -0.7038178, 0.02047591, 0.6381929, -0.48880103), ncol = 2, byrow = TRUE)
  test_values$growthrates <-c( -1.065035, -1.636435)

  output_pred <- EEMtoolbox::ode_solve(
    interaction_matrix=test_values$interaction_matrix,
    growth_rate=test_values$growthrates,
    t_window = c(0,10),
    model = model_test,
    initial_condition =initcond-runif(2,min=-0.5,max=0.5)*initcond
  )
  expect_true(max(abs((initcond-output_pred$abundances[nrow(output_pred$abundances),])/initcond))<= 0.1)


  #Baker
  model_test = "Baker"
  initcond <- c(17.06103, 12.35617)

  test_values <- list()
  test_values$interaction_matrix_alphas <- matrix(c(-0.8185450, 0.1201095,0.9095789, -0.1474479), ncol = 2, byrow = TRUE)
  test_values$interaction_matrix_betas <- matrix(c(-0.3942834, 0.6435232,0.2984477,-0.7150782), ncol = 2, byrow = TRUE)
  test_values$growthrates <-c(-2.519711, 3.743805)

  output_pred <- EEMtoolbox::ode_solve(
    interaction_matrix=list(test_values$interaction_matrix_alphas,
                            test_values$interaction_matrix_betas),
    growth_rate=test_values$growthrates,
    t_window = c(0,10),
    model = model_test,
    initial_condition =initcond-runif(2,min=-0.5,max=0.5)*initcond
  )

  expect_true(max(abs((initcond-output_pred$abundances[nrow(output_pred$abundances),])/initcond))<= 0.2)
})
