#' @title Arguments for EEM
#' @description
#' Extract arguments necessary to run EEM from interaction matrix
#'
#' @param prior_sample matrix of prior sample of parameters.
#' @param posterior_sample matrix of posterior sample of parameters.
#' @param param_names vector of names for each parameter.
#' @examples
#' library(EEMtoolbox)
#' output <- EEM(dingo_matrix,  #automatically loads an example of interaction matrix as dingo_matrix
#'               output_prior=TRUE,
#'               output_discrepancy=TRUE,
#'               output_matrix=FALSE)
#' ix <- which(output$part_s==0) #indexes of interest
#' prior_sample <- output$prior_sample
#' posterior_sample <- output$part_vals[ix,]
#' param_names <- seq(ncol(prior_sample))
#' plots_posterior_density(prior_sample,posterior_sample,param_names)
#' @return ggplot2 figure
#'
#' @export
plots_posterior_density <- function(prior_sample,
                                    posterior_sample,
                                    param_names){
  #figure layout settings
  n_params <- ncol(posterior_sample)

  #for each parameter
  data_densities <- data.frame(xx=numeric(),
                               yy=numeric(),
                               step=character(),
                               param=character())
  for (i in seq(n_params)){
    dens_prior <- stats::density(prior_sample[,i])
    dens_prior <- data.frame(xx=dens_prior$x, yy=dens_prior$y)
    dens_prior$step <- "prior"

    dens_posterior <- stats::density(posterior_sample[,i],bw = "nrd")
    dens_posterior <- data.frame(xx=dens_posterior$x, yy=dens_posterior$y)
    dens_posterior$step <- "posterior"
    data_dens_plot_i <- rbind(dens_prior,dens_posterior)
    data_dens_plot_i$param <- param_names[i]
    data_densities <- rbind(data_densities, data_dens_plot_i)
  }

  g <- ggplot2::ggplot(data_densities, ggplot2::aes(x=xx,
                             y=yy,
                             group = interaction(step),
                             col = step))+
    ggplot2::geom_line()+
    ggplot2::facet_wrap(~param,scales = "free")+
    ggplot2::theme(text = ggplot2::element_text(size = 10))+
    ggplot2::labs(x="Value",y="Density")
  return(g)
}
