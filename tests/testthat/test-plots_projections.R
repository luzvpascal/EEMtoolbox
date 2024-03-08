test_that("plots_projections works", {
  output <- EEMtoolbox::EEM(matrix(c(-1,-1,1,-1),ncol=2))
  plot <- EEMtoolbox::plots_projections(output,  c(1,1), t_window=c(0,1))


  expect_true(class(plot)[1]=="gg")
  expect_true(class(plot)[2]=="ggplot")

})
