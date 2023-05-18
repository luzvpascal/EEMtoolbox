<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/luzvpascal/EEMtoolbox/branch/main/graph/badge.svg)](https://app.codecov.io/gh/luzvpascal/EEMtoolbox?branch=main)
  <!-- badges: end -->

# EEMtoolbox
EEMtoolbox is an R-package that efficiently generates an ensemble of plausible quantitative models that describe an ecosystem from a species interaction network. 

EEMtoolbox supports three different models that represent species interactions: Generalized Lokta Voltera, Baker model and Gompertz model. Following the requirements from (Baker et al., 2017), the generated models verify ecosytem coexistence (feasibility: steady states positive) and stability (eigen values of Jacobian negative). Customized models can also be provided by the user ([click here](#customizing-input-model)).

Our package generates ensemble members in two possible ways: standard EEM (Baker et al., 2017) and EEM-SMC (Vollert et al. in preparation). The standard EEM method uniformly samples the parameter space until the desired number of ensemble members is generated, which has proven to be efficient for small networks. The EEM-SMC method takes advantage of Approximate Bayesian Computation methods (Drovandi and Pettitt 2011), which can speed up the generation of ensemble members specially for large networks.

## Installation
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

`algorithm`: algorithm used for sampling. Default "SMC-ABC" (Vollert et al., 2023) options include "standard EEM".

## Predicting species abundances 
```r
library(tidyverse)

index <- 1
initcond <- summarise_ecosystem_features(parameters = outputs$part_vals[index,],
                                         sim_args = outputs$sim_args)
discrepancy_continuous_sum(initcond)
initcond <- initcond[seq(8)]
test_values <- EEMtoolbox::reconstruct_matrix_growthrates(outputs$part_vals[index,],sim_args = outputs$sim_args)

output_pred <- EEMtoolbox::ode_solve(interaction_matrix=test_values$interaction_matrix,
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


## Customizing input model 

## Customizing search algorithm
