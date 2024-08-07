test_that("n_cores_function works", {
  #sequential
  ncores <- EEMtoolbox::n_cores_function()
  expect_true(ncores==1L)

  #sequential
  ncores <- EEMtoolbox::n_cores_function(2L)
  expect_true(ncores==2L)
})
