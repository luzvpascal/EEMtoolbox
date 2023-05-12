<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/luzvpascal/EEMtoolbox/branch/main/graph/badge.svg)](https://app.codecov.io/gh/luzvpascal/EEMtoolbox?branch=main)
  <!-- badges: end -->

# EEMtoolbox
EEMtoolbox is an R-package that efficiently generates an ensemble of plausible quantitative models that describe an ecosystem from a species interaction network. 

EEMtoolbox supports three different models that represent species interactions: Generalized Lokta Voltera, Baker model and Gompertz model. Following the requirements from (Baker et al., 2017), the generated models verify ecosytem coexistence (feasibility: steady states positive) and stability (eigen values of Jacobian negative). Customized models can also be provided by the user [(click here)](#customizing-input-model)

Our package includes the implementation of recent work from our project team (Vollert et al. in preparation), which has demonstrated that Approximate Bayesian Computation methods (Drovandi and Pettitt 2011) can speed up the generation

## Installation
To install EEMtoolbox, run the following line
``` r
devtools::install_github("luzvpascal/EEMtoolbox", host = "https://api.github.com")
```

## Running EEM
The main function of EEMtoolbox is `EEM`. This function generates ensemble members
```r
library(EEMtoolbox)
EEM(dingo_matrix) #dingo_matrix is included in the package
```

`interaction_matrix`: interaction signs matrix, can be input as a single matrix of interactions or as a list of matrices defining lower and upper bounds for interaction terms lower first and upper second

`bounds_growth_rate`: vector of 2 elements containing lower and upper bounds for growth rates. Default c(-5,5)

`n_ensemble`: Number of desired ensemble members. Default to 10

`model`: model representing species interactions. Default "GLV" (Generalized Lokta Voltera). options include "Baker", "Gompertz" and "customized"

`algorithm`: algorithm used for sampling. Default "SMC-ABC" (Vollert et al., 2023) options include "standard EEM"

## Predicting species abundances 

## Customizing input model 

## Customizing search algorithm
