test_that("plot_posterior_density works", {

  dingo_matrix <- matrix(c(-1, -1, -1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0,
                           1, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, 0, -1, 0, 0,
                           0, 0, 1, 1, -1, 1, 1, 0, 0, 0, 0, 1, 1, -1, 0, 1, 0,
                           0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, 1, 0, -1), nrow=8)

  output <- EEMtoolbox::EEM(dingo_matrix,
                output_prior=TRUE,
                output_discrepancy=TRUE,
                output_matrix=FALSE)
  ix <- which(output$part_s==0) #indexes of interest
  prior_sample <- output$prior_sample
  posterior_sample <- output$part_vals[ix,]
  param_names <- seq(ncol(prior_sample))
  plot <- EEMtoolbox::plot_posterior_density(prior_sample,posterior_sample,param_names)

  expect_true(class(plot)[1]=="gg")
  expect_true(class(plot)[2]=="ggplot")

})
