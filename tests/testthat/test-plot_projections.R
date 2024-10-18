test_that("plot_projections works", {
  output <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol=2))
  plot <- EEMtoolbox::plot_projections(output,  c(1,1), t_window=c(0,1))


  expect_true(class(plot)[1]=="gg")
  expect_true(class(plot)[2]=="ggplot")

  plot <- EEMtoolbox::plot_projections(output,  c(1,1), t_window=c(0,1),scaled=TRUE)


  expect_true(class(plot)[1]=="gg")
  expect_true(class(plot)[2]=="ggplot")


})
