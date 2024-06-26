% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_posterior_density.R
\name{plot_posterior_density}
\alias{plot_posterior_density}
\title{plot_posterior_density}
\usage{
plot_posterior_density(prior_sample, posterior_sample, param_names)
}
\arguments{
\item{prior_sample}{matrix of prior sample of parameters.}

\item{posterior_sample}{matrix of posterior sample of parameters.}

\item{param_names}{vector of names for each parameter.}
}
\value{
ggplot2 figure
}
\description{
Extract arguments necessary to run EEM from interaction matrix
}
\examples{
library(EEMtoolbox)
output <- EEM(dingo_matrix,  #automatically loads an example of interaction matrix as dingo_matrix
              output_prior=TRUE,
              output_discrepancy=TRUE,
              output_matrix=FALSE)
prior_sample <- output$prior_sample
posterior_sample <- output$part_vals
param_names <- seq(ncol(prior_sample))
plot_posterior_density(prior_sample,posterior_sample,param_names)
}
