#' @title Number of cores
#' @description
#' Gets number of cores available for parallel tasks
#' @param n_cores Number of cores desired for parallel sampling. Default set to 1 core (sequential sampling).
#' @export
n_cores_function <- function(n_cores=1L){
  #setup parallel backend to use many processors
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    max_cores <- 2L
  } else {
    # use all cores in devtools::test()
    max_cores <- parallelly::availableCores(omit = 1) #do not overload your computer
  }
  return(min(max_cores, n_cores))
}
