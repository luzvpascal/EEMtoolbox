test_that("plot_posterior_density works", {

  interaction_matrix <- matrix(c(-1,1,1,-1),ncol = 2)

  output <- EEMtoolbox::EEM(interaction_matrix,
                            output_prior=TRUE,
                            output_discrepancy=TRUE,
                            output_matrix=FALSE)
  prior_sample <- output$prior_sample
  posterior_sample <- output$part_vals
  param_names <- seq(ncol(prior_sample))
  plot <- EEMtoolbox::plot_posterior_density(prior_sample,posterior_sample,param_names)

  expect_true(class(plot)[1]=="gg")
  expect_true(class(plot)[2]=="ggplot")

})
