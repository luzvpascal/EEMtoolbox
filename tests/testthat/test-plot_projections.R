test_that("plot_projections works", {
  #non-scaled average
  output <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol=2))
  plot <- EEMtoolbox::plot_projections(output,  c(1,1), t_window=c(0,1))


  expect_true(class(plot)[1]=="gg")
  expect_true(class(plot)[2]=="ggplot")

  #scaled average
  plot <- EEMtoolbox::plot_projections(output,  c(1,1), t_window=c(0,1),scaled=TRUE)

  expect_true(class(plot)[1]=="gg")
  expect_true(class(plot)[2]=="ggplot")

  #non-scaled average
  output <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol=2))
  plot <- EEMtoolbox::plot_projections(output,  c(1,1), t_window=c(0,1),average=FALSE)


  expect_true(class(plot)[1]=="gg")
  expect_true(class(plot)[2]=="ggplot")

  #scaled average
  plot <- EEMtoolbox::plot_projections(output,  c(1,1), t_window=c(0,1),
                                       average=FALSE,scaled=TRUE)

  expect_true(class(plot)[1]=="gg")
  expect_true(class(plot)[2]=="ggplot")

  # Bimler-baker
  model_test = "Bimler-Baker"
  interaction_matrix_alphas <- matrix(c(1, 1,1, 1), ncol = 2, byrow = TRUE)
  interaction_matrix_betas <- matrix(c(-1, -1,-1, -1), ncol = 2, byrow = TRUE)
  output <- EEMtoolbox::EEM(list(interaction_matrix_alphas,interaction_matrix_betas),
                            model=model_test)
  plot <- EEMtoolbox::plot_projections(output,  c(1,1), t_window=c(0,1), model=model_test)
  plot2 <- EEMtoolbox::plot_projections(output,  c(1,1), t_window=c(0,1), model=model_test,scaled=TRUE)

  expect_true(class(plot)[1]=="gg")
  expect_true(class(plot)[2]=="ggplot")
  expect_true(class(plot2)[1]=="gg")
  expect_true(class(plot2)[2]=="ggplot")

})
