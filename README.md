<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/luzvpascal/EEMtoolbox/branch/main/graph/badge.svg)](https://app.codecov.io/gh/luzvpascal/EEMtoolbox?branch=main)
  <!-- badges: end -->

# EEMtoolbox
EEMtoolbox is an R-package that efficiently integrates species interaction networks in order to propose an ensemble of plausible quantitative models that describe an ecosystem.

## Installation
To install EEMtoolbox, we recommend running the following line
``` r
devtools::install_github("luzvpascal/EEMtoolbox", host = "https://api.github.com")
library(smsPOMDP)
```

## Running EEM
The main function of EEMtoolbox is `EEM`.
```r
EEM(interaction_matrix,
    bounds_growth_rate=c(-5,5),
    n_ensemble=10,
    model="GLV",
    algorithm="SMC-ABC",
    summ_func=EEMtoolbox::summarise_ecosystem_features,
    disc_func=EEMtoolbox::discrepancy_continuous_sum,
    sampler=EEMtoolbox::uniform_sampler,
    trans_f=EEMtoolbox::uniform_transform,
    trans_finv=EEMtoolbox::uniform_transform_inverse,
    pdf=EEMtoolbox::uniform_pdf_transformed,
    n_particles=10000,
    mcmc_trials=10,
    dist_final=0,
    a=0.6,
    c=0.01,
    p_acc_min=0.0001,
    output_prior=FALSE,
    output_args=FALSE,
    output_discrepancy=FALSE,
    output_matrix=TRUE)
```

`interaction_matrix`: interaction signs matrix, can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second

`bounds_growth_rate`: vector of 2 elements containing lower and upper bounds for growth rates. Default c(-5,5)

`n_ensemble`: Number of desired ensemble members. Default to 10

`model`: model representing species interactions. Default "GLV" (Generalized Lokta Voltera). options include "Baker", "Gompertz" and "customized"

`algorithm`: algorithm used for sampling. Default "SMC-ABC" (Vollert et al., 2023) options include "standard EEM"
