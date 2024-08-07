% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EEM_standard_method_iteration.R
\name{EEM_standard_method_iteration}
\alias{EEM_standard_method_iteration}
\title{Generation of model ensembles using standard EEM sampling method}
\usage{
EEM_standard_method_iteration(
  sim_args,
  summ_func,
  disc_func,
  sampler,
  trans_f,
  n_ensemble = 5000,
  n_cores = 1L
)
}
\arguments{
\item{sim_args}{a list of arguments as returned by \link[EEMtoolbox]{args_function}}

\item{summ_func}{function calculating equilibrium points and real parts of the Jacobians eigenvalues to summarise ecosystem features.}

\item{disc_func}{summary statistic (discrepancy measure).}

\item{sampler}{sampling function that generates random vectors from the joint prior distribution.}

\item{trans_f}{transform of prior parameter space to ensure unbounded support for MCMC sampling.}

\item{n_ensemble}{Number of desired ensemble members. Default to 5000}

\item{n_cores}{Number of cores desired to be used for sampling. Default set to 1 core (sequential sampling).}
}
\value{
list: sims=number of simulations
part_vals=parameter values
part_s=discrepancy value
prior_sample=prior distribution
}
\description{
Generation of model ensembles using standard EEM sampling method
}
