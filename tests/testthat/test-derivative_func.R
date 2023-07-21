test_that("derivative_func works", {
  interaction_matrix <- matrix(c(-1,1,1,-1),ncol=2)
  growth_rate <- c(1,1)
  current_abundance <- c(1,1)

  ## GLV
  Pars <- list(interaction_matrix_value=interaction_matrix,
               growth_rate=growth_rate,
               model="GLV")
  glv_deriv <- EEMtoolbox::derivative_func(Time=0,
                               State=current_abundance,
                               Pars=Pars)
  expect_equal(unlist(glv_deriv), c(1,1))


  ## Gompertz
  Pars <- list(interaction_matrix_value=interaction_matrix,
               growth_rate=growth_rate,
               model="Gompertz")
  gompertz_deriv <- EEMtoolbox::derivative_func(Time=0,
                                                State=current_abundance,
                                                Pars=Pars)
  expect_equal(c(unlist(gompertz_deriv)), c(1,1))

  ## Baker
  # diag(interaction_matrix) <- 0
  Pars <- list(interaction_matrix_value=list(interaction_matrix,
                                             interaction_matrix),
               growth_rate=growth_rate,
               model="Baker")
  baker_deriv <- EEMtoolbox::derivative_func(Time=0,
                                             State=current_abundance,
                                             Pars=Pars)
  expect_equal(c(unlist(baker_deriv)), c(0,0))

})
