<!-- badges: start -->
[![codecov](https://codecov.io/gh/luzvpascal/EEMtoolbox/branch/main/graph/badge.svg?token=MH9JLF9HEQ)](https://codecov.io/gh/luzvpascal/EEMtoolbox)
[![Build Status](https://app.travis-ci.com/luzvpascal/EEMtoolbox.svg?branch=main)](https://app.travis-ci.com/luzvpascal/EEMtoolbox)
[![R-CMD-check](https://github.com/luzvpascal/EEMtoolbox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/luzvpascal/EEMtoolbox/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/621079578.svg)](https://doi.org/10.5281/zenodo.14032063)
<!-- badges: end -->

# EEMtoolbox
EEMtoolbox is an R-package that efficiently generates an ensemble of plausible quantitative models that describe an ecosystem from a species interaction network.

EEMtoolbox supports three different models that represent species interactions: Generalized Lokta Voltera, Baker model and Gompertz model. Following Baker et al., (2017), the generated models must be feasible (coexistence: positive equilibrium abundances) and stable (negative eigenvalues of Jacobian). Customized models can also be provided by the user.

Our package generates ensemble members in two possible ways: standard EEM (Baker et al., 2017) and EEM-SMC (Vollert et al. in preparation). Both methods can generate representative and equivalent ensembles.The standard EEM method samples the parameter space until the desired number of ensemble members is generated, which has proven to be efficient for small networks. The EEM-SMC method takes advantage of Approximate Bayesian Computation methods (Drovandi and Pettitt 2011), which can speed up the generation of ensemble members specially for large networks.

## Installation

Users need to make sure that their version of R is at least 4.3.1. We recommend running the following code in R to preinstall all the necessary packages:
``` r
packages_to_install <- c("deSolve", "doParallel", "doSNOW","dplyr", "foreach", "ggplot2", "magrittr", "MASS", "nleqslv", "parallel", "parallelly","stats","tidyr")

# Install packages if not already installed
install_if_not_installed <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

# Apply the function to install packages
invisible(lapply(packages_to_install, install_if_not_installed))
```

To install EEMtoolbox, download the repository as a zip file.
Rename the file as EEMtoolbox.zip, then run the following code
``` r
devtools::install_local("path_to_file/EEMtoolbox.zip", repos = NULL, type = "win.binary")
```
