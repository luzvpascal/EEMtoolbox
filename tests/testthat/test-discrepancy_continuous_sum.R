test_that("test discrepancy function", {

  test_vect <- rep(0,8)
  expect_equal(EEMtoolbox::discrepancy_continuous_sum(test_vect), 0)

  test_vect2 <- c(rep(1,4), rep(-1,4)) #test stability and feasibility
  expect_equal(EEMtoolbox::discrepancy_continuous_sum(test_vect2), 0)

  test_vect3 <- c(rep(-1,4), rep(-1,4)) #infeasible
  expect_equal(EEMtoolbox::discrepancy_continuous_sum(test_vect3), 4)

  test_vect4 <- c(rep(1,4), rep(1,4)) #unstable
  expect_equal(EEMtoolbox::discrepancy_continuous_sum(test_vect4), 4)

  test_vect5 <- c(rep(-1,4), rep(1,4)) #infeasible and unstable
  expect_equal(EEMtoolbox::discrepancy_continuous_sum(test_vect5), 8)

})
