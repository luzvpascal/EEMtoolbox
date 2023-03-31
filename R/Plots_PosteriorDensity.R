#' @export
Plots_PosteriorDensity <- function(prior_sample, posterior_sample, prior_args,param_names){
  #figure layout settings
  n_params <- ncol(posterior_sample)

  #locally extract prior bounds
  prior_lowers <- prior_args$lower
  prior_uppers <- prior_args$upper

  #check if we have zero-parameters included in sample
  # if (n_params != length(prior_args$lower)){#we have zero parameters
  #   interaction_sample_nonzero <- matrix(0,nrow(posterior_sample),length(prior_args$lower))
  #   count =0
  #   for (i in seq(length(prior_args$skip_parameters))){
  #     if (prior_args$skip_parameters[i]==0){#if values are non-zero
  #       count <- count+1
  #       interaction_sample_nonzero[,count] <- posterior_sample[,prior_args$n_species+i]
  #     }
  #   }
  #   posterior_sample <- c(posterior_sample[,seq(prior_args$n_species)], interaction_sample_nonzero)
  #   n_params <- ncol(posterior_sample)
  # }

  #for each parameter
  data_densities <- data.frame(xx=numeric(),
                               xx=numeric(),
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
  g <- ggplot2::ggplot(data_densities, aes(x=xx,
                             y=yy,
                             group = interaction(step),
                             col = step))+
    ggplot2::geom_line()+
    ggplot2::facet_wrap(~param,scales = "free")+
    ggplot2::theme(text = element_text(size = 5))
  return(g)
}
