#' @title Arguments for EEM
#' @description
#' Extract arguments necessary to run EEM from interaction matrix
#'
#' @param prior_sample interaction signs matrix, can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second
#' @param posterior_sample vector of 2 elements containing lower and upper bounds for growth rates
#' @param sim_args vector of 2 elements containing lower and upper bounds for growth rates
#' @param param_names model representing species interactions, default "GLV" (Generalized Lokta Voltera). options include "Baker", "Adams" and "customized"
#' @return A list of arguments defining the problem.
#' @export
plots_posterior_density <- function(prior_sample, posterior_sample, sim_args,param_names){
  #figure layout settings
  n_params <- ncol(posterior_sample)

  #locally extract prior bounds
  prior_lowers <- sim_args$lower
  prior_uppers <- sim_args$upper

  #for each parameter
  data_densities <- data.frame(xx=numeric(),
                               yy=numeric(),
                               step=character(),
                               param=character())
  for (i in seq(n_params)){
    dens_prior <- stats::density(prior_sample[,i],bw = "nrd")
    dens_prior <- data.frame(xx=dens_prior$x, yy=dens_prior$y)
    # dens_prior <- reflection_function(dens, prior_lowers[i], prior_uppers[i])
    dens_prior$step <- "prior"

    dens_posterior <- stats::density(posterior_sample[,i],bw = "nrd")
    # dens_posterior <- reflection_function(dens_post, prior_lowers[i], prior_uppers[i])
    dens_posterior <- data.frame(xx=dens_posterior$x, yy=dens_posterior$y)
    dens_posterior$step <- "posterior"
    data_dens_plot_i <- rbind(dens_prior,dens_posterior)
    data_dens_plot_i$param <- param_names[i]
    data_densities <- rbind(data_densities, data_dens_plot_i)
  }

  write.csv(data_densities, "data_densities.csv", row.names = FALSE)
  g <- ggplot2::ggplot(data_densities, ggplot2::aes(x=xx,
                             y=yy,
                             group = interaction(step),
                             col = step))+
    ggplot2::geom_line()+
    ggplot2::facet_wrap(~param,scales = "free")+
    ggplot2::theme(text = ggplot2::element_text(size = 5))
  return(g)
}
