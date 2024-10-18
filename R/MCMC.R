#' @title Monte Carlo Markov Chain
#' @description
#' Runs MCMC
#' @param i line number of particle to be moved
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @param mcmc_trials number of MCMC steps to try before selecting appropriate number.
#' @param dist_next next objective distance
#' @param part_vals matrix of current values of parameters of each particle (each particle represented by a row)
#' @param part_s vector of discrepancy measures of each particle (each particle represented by a row)
#' @param part_sim matrix of current simulation values: equilibriums and eigenvalues of jacobian (each particle represented by a row)
#' @param cov_matrix covariance matrix
#' @param summ_func function calculating equilibrium points and real parts of the Jacobians eigenvalues to summarise ecosystem features. Default =summarise_ecosystem_features_GLV. Options include summarise_ecosystem_features_Baker (automatically chosen if model="Bimler-Baker") and summarise_ecosystem_features_Gompertz, (automatically chosen if model="Gompertz"). Needs to be defined if model="customized" chosen.
#' @param disc_func summary statistic (discrepancy measure).
#' @param trans_finv inverse of trans_f function.
#' @param pdf joint probability density function.
#' @return list
#' part_vals: updated value of parameter values for particle i
#' part_s: discrepancy measure for particle i
#' part_sim: summary of ecosystem features for particle i
#' i_acc: number of times particle movement accepted
#' sims_mcmc: number of successful walks
#' @export
#' @import MASS

MCMC <- function(i,
                 sim_args,
                 mcmc_trials,
                 dist_next,
                 part_vals,
                 part_s,
                 part_sim,
                 cov_matrix,
                 summ_func,
                 disc_func,
                 trans_finv,
                 pdf){
  # documentation
  sims_mcmc <-  0
  i_acc <- 0
  for (r in seq(mcmc_trials)){
    # Gaussian random walk
    part_vals_prop <- MASS::mvrnorm(n=1,mu=part_vals[i,],Sigma=cov_matrix)

    # Transform back to calculate prior probs and discrepancy
    prop <- trans_finv(matrix(part_vals_prop, nrow=1), sim_args)

    # Calculate prior probabilities
    prior_curr <- pdf(matrix(part_vals[i,], nrow=1))
    prior_prop <- pdf(matrix(part_vals_prop, nrow=1))

    # early rejection (assumes symmetric proposals)
    if(((is.nan(prior_prop/prior_curr))| (runif(1) > prior_prop/prior_curr))){
      next
    }

    #find proposal discrepancy
    part_sim_prop <- summ_func(prop, sim_args)
    dist_prop <- disc_func(part_sim_prop)

    sims_mcmc=sims_mcmc+1
    # ABC part of the acceptance probability

    #Accept a particle if it is within the target distance.
    if (dist_prop <= dist_next) {
      # then the metropolis-hastings ratio is satisfied
      part_vals[i,] <- part_vals_prop
      part_s[i] <- dist_prop
      part_sim[i,] <- part_sim_prop
      i_acc <- i_acc + 1
    }
  }
  return(list(part_vals=part_vals[i,],
              part_s=part_s[i],
              part_sim=part_sim[i,],
              i_acc=i_acc,
              sims_mcmc=sims_mcmc))
}
