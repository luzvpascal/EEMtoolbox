test_that("n_cores_function works", {
  #sequential
  ncores <- EEMtoolbox::n_cores_function()
  cl <- parallel::makeCluster(ncores[1])
  doParallel::registerDoParallel(cl)
  start <- Sys.time()
  sum_test <- foreach::foreach(i = 1:10000) %dopar% {
    prod(i,i)
  }
  end <- Sys.time()
  #stop cluster
  parallel::stopCluster(cl)
  rm(cl)

  time_sequential <- end - start

  #parallel
  ncores_parallel <- EEMtoolbox::n_cores_function(2L)
  cl_parallel <- parallel::makeCluster(ncores_parallel[1])
  doParallel::registerDoParallel(cl_parallel)
  start_parallel <- Sys.time()
  sum_test <- foreach::foreach(i = 1:10000) %dopar% {
    prod(i,i)
  }
  end_parallel <- Sys.time()
  #stop cluster
  parallel::stopCluster(cl_parallel)
  rm(cl_parallel)

  time_parallel <- end_parallel - start_parallel

  expect_true(time_parallel<=time_sequential)
})
