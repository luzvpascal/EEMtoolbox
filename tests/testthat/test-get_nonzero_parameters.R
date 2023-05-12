test_that("get_nonzero_parameters works", {

  interaction_matrix <- matrix(c(-1,1,1,-1),ncol=2)
  non_zero <- EEMtoolbox::get_nonzero_parameters(interaction_matrix )

  interaction_matrix_list <- list(matrix(pmin(interaction_matrix,0), ncol=2),
                                  matrix(pmax(interaction_matrix,0), ncol=2))
  non_zero_list <- EEMtoolbox::get_nonzero_parameters(interaction_matrix_list)
  expect_equal(non_zero$skip_parameters, non_zero_list$skip_parameters)
  expect_equal(non_zero$lower_interaction_bound, non_zero_list$lower_interaction_bound)
  expect_equal(non_zero$upper_interaction_bound, non_zero_list$upper_interaction_bound)
})
