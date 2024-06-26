<!-- badges: start -->
[![codecov](https://codecov.io/gh/luzvpascal/EEMtoolbox/branch/main/graph/badge.svg?token=MH9JLF9HEQ)](https://codecov.io/gh/luzvpascal/EEMtoolbox)
[![Build Status](https://app.travis-ci.com/luzvpascal/EEMtoolbox.svg?branch=main)](https://app.travis-ci.com/luzvpascal/EEMtoolbox)
[![R-CMD-check](https://github.com/luzvpascal/EEMtoolbox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/luzvpascal/EEMtoolbox/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# EEMtoolbox
EEMtoolbox is an R-package that efficiently generates an ensemble of plausible quantitative models that describe an ecosystem from a species interaction network.

EEMtoolbox supports three different models that represent species interactions: Generalized Lokta Voltera, Baker model and Gompertz model. Following Baker et al., (2017), the generated models must be feasible (coexistence: positive equilibrium abundances) and stable (negative eigenvalues of Jacobian). Customized models can also be provided by the user ([click here](#customizing-input-model)).

Our package generates ensemble members in two possible ways: standard EEM (Baker et al., 2017) and EEM-SMC (Vollert et al. in preparation). Both methods can generate representative and equivalent ensembles.The standard EEM method samples the parameter space until the desired number of ensemble members is generated, which has proven to be efficient for small networks. The EEM-SMC method takes advantage of Approximate Bayesian Computation methods (Drovandi and Pettitt 2011), which can speed up the generation of ensemble members specially for large networks.

## Table of contents
* [Installation and usage](#installation-and-usage)
* [Three ecosystem interactions models](#three-ecosystem-interactions-models)

# Installation and usage
## Installation

Users need to make sure that their version of R is at least 4.3.1. We recommend running the following code in R to preinstall all the necessary packages:
``` r
packages_to_install <- c("deSolve", "doParallel", "doSNOW", "foreach", "ggplot2", "MASS", "nleqslv", "parallel", "stats")

# Install packages if not already installed
install_if_not_installed <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

# Apply the function to install packages
invisible(lapply(packages_to_install, install_if_not_installed))
```

To install EEMtoolbox, run the following line
``` r
devtools::install_github("luzvpascal/EEMtoolbox", host = "https://api.github.com")
```

## Running EEM
The main function of EEMtoolbox is `EEM`. This function inputs an `interaction_matrix` and outputs ensemble members (for the Generalized Lokta Voltera model by default).
```r
library(EEMtoolbox)
outputs <- EEM(dingo_matrix) #dingo_matrix is included in the package
```
**Other arguments of `EEM` function**

`interaction_matrix`: interaction signs matrix, can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second

`bounds_growth_rate`: vector of 2 elements containing lower and upper bounds for growth rates. Default c(-5,5).

`n_ensemble`: Number of desired ensemble members. Default to 10.

`model`: model representing species interactions. Default "GLV" (Generalized Lokta Voltera). Options include "Baker", "Gompertz" and "customized".

`algorithm`: algorithm used for sampling. Default "standard EEM".

## Predicting species abundances
```r
library(tidyverse)
system <- outputs[[1]] #select first ensemble returned by EEM.
output_pred <- EEMtoolbox::ode_solve(interaction_matrix=system$interaction_matrix,
                                     growth_rate=test_values$growthrates,
                                     t_window = c(0,10),
                                     model = model_test,
                                     initial_condition =initcond-runif(8,min=-0.5,max=0.5)*initcond
                                    )
```

## Ploting predictions
```r
abundance <- as.data.frame(output_pred$y)
names(abundance) <- species_list
abundance$time <- seq(nrow(abundance))

abundance <- abundance %>%
  pivot_longer(!time, names_to = "species", values_to = "pop")

p <- ggplot(abundance, aes(x=time, y=pop, color=species, fill = species)) +
  stat_summary(geom = "line", fun = mean) +
  guides(fill = guide_legend(title = "Species Group"), color = guide_legend(title="Species Group")) +
  theme_bw() +
  xlab("Years") +
  ylab("Abundance") +
  scale_x_continuous(breaks = seq(1,nrow(output_pred$y), length.out=11)-1,
                     labels = (seq(1,nrow(output_pred$y), length.out=11)-1)/100)
p
```

# Three ecosystem interactions models
Ecosystem dynamics can be modelled using ordinary differential equations (ODE), which seek to predict species abundances over time. Here we present the three types of models represented as ODEs that are supported by our package. For each model, we provide the feasibility and stability requirements.

## Generalized Lotka Volterra model
The Generalized Lokta Voltera equations are

$$
\frac{dn_i}{dt} = \left[ r_i + \sum_{i=1}^N \alpha_{i,j} n_j(t) \right] n_i(t)
$$

where:
- $n_i(t)$ is the abundance of the ith ecosystem node at time $t$
- $r_i$ is the intrinsic growth rate of the i-th species
- $N$ is the total number of species
- $\alpha_{i,j}$ is the per-capita interaction strength representing the effect of the j-th species the i-th.

These equations can be rewritten as:

$$
\frac{d\mathbf{n}}{dt} = \left[\mathbf{r} + \mathbf{A} \mathbf{n}\right]\circ \mathbf{n}
$$

where:
- $\mathbf{n}$ is the vector of abundances
- $\mathbf{r}$ is the vector of intrinsic growth rates
- $\mathbf{A}$ is the $N \times N$ interaction matrix of per-capita interaction strengths, where $A_{i,j}=\alpha_{i,j}$.
- $\circ$ is the Hadamard (element wise) product.

### Feasibility
The solution to $\frac{d\mathbf{n}^{\*}}{dt} = \left[\mathbf{r} + \mathbf{A} \mathbf{n}^{\*}\right]\circ \mathbf{n}^{\*} = 0$ is

$$
\mathbf{n}^{\*} = - \mathbf{A}^{-1} \mathbf{r}.
$$

The feasibility condition is verified if $\mathbf{n}^{\*} > 0$.

### Stability
To verify the (Lyapunov) stability of the solution $\mathbf{n}^{\*}$, we need to calculate the eigen values of the Jacobian matrix $J$ at equilibrium $\mathbf{n}^{\*}$. The Jacobian matrix $J$ is defined as:

$$
J_{i,j} = \left(\frac{\partial f_i}{\partial n_j}\right)(n=n^{\*}), \quad \text{where  }f_i(n_i) = \left[ r_i + \sum_{i=1}^N \alpha_{i,j} n_j(t) \right] n_i(t)
$$

The Jacobian is thus defined as:

$$
J_{i,j} = \alpha_{i,j} n_i^{\*}(t), \quad \forall i \neq j
$$

$$
J_{i,i} = r_i + \sum_{i=1}^N \alpha_{i,j} n_j^{\*}(t) + \alpha_{i,i} n_i^{\*}(t)
$$

The system is considered stable if the real part of each eigenvalue ($\lambda_i$) is negative, i.e. $Re(\lambda_i) \leq 0 $.

## Baker model

### Feasibility

### Stability

## Gompertz model
The ODEs of the Gompertz model are

$$
\frac{dn_i}{dt} = \left[ r_i + \sum_{i=1}^N \alpha_{i,j} log(n_j(t)) \right] n_i(t)
$$

where:
- $n_i(t)$ is the abundance of the ith ecosystem node at time $t$
- $r_i$ is the intrinsic growth rate of the i-th species
- $N$ is the total number of species
- $\alpha_{i,j}$ is the per-capita interaction strength representing the effect of the j-th species the i-th.

### Feasibility

The solution of the Gompertz model ODEs is:

$$
\mathbf{n}^{\*} = e^{- \mathbf{A}^{-1} \mathbf{r}}.
$$

This solution $\mathbf{n}^{\*}$ is always positive, therefore the ecosystem is always feasible in the Gompertz model.

### Stability

The Jacobian matrix is given by:

$$
J_{i,j} = \alpha_{i,j} \frac{n_i^{\*}}{n_j^{\*}}.
$$

The system is considered stable if the real part of each eigenvalue ($\lambda_i$) is negative, i.e. $Re(\lambda_i) \leq 0 $.

