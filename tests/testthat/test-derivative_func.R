test_that("derivative_func works", {
  interaction_matrix <- matrix(c(-1,1,1,-1),ncol=2)
  growth_rate <- c(1,1)
  current_abundance <- c(1,1)

  ## GLV
  glv_deriv <- EEMtoolbox::derivative_func(interaction_matrix,
                               growth_rate,
                               current_abundance,
                               model="GLV")
  expect_equal(c(glv_deriv), c(1,1))


  ## Gompertz
  gompertz_deriv <- EEMtoolbox::derivative_func(interaction_matrix,
                                             growth_rate,
                                             current_abundance,
                                             model="Gompertz")
  expect_equal(c(gompertz_deriv), c(1,1))

  ## Baker
  diag(interaction_matrix) <- 0
  baker_deriv <- EEMtoolbox::derivative_func(list(interaction_matrix,
                                                  interaction_matrix),
                                             c(0,0),
                                             current_abundance,
                                             model="Baker")
  expect_equal(c(baker_deriv), c(1,1))

  ## customized
  customized_deriv <- EEMtoolbox::derivative_func(interaction_matrix,
                                             c(0,0),
                                             current_abundance,
                                             model="customized")
  expect_equal(customized_deriv, 0)

})
