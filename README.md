dfSEDI
================

# dfSEDI

The **dfSEDI** package implements semiparametric efficient data
integration methods for dual-frame sampling as described in Morikawa &
Kim (202x).

The package provides the following main estimators:

- `Eff` : semiparametric efficient estimator (DML2, K-fold)
- `Eff_S` : sub-efficient estimator (based on Remark 6, DML2, K-fold)
- `Eff_P` : parametric efficient estimator (working model)

Internally, the package also implements basic estimators such as

- `df_estimate_P` : probability-only estimator
- `df_estimate_NP` : non-probability-only estimator (Chang & Kott type)
- `df_estimate_NP_P` : combined NP ∪ P estimator (Chang & Kott type)

For each estimator, sandwich-based variance, standard error, and 95%
confidence intervals are computed from influence-function-based
contributions.

## Installation

You can install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("KMorikawaISU/dfSEDI")
```

## Example: one simulated dual-frame dataset

This section shows a minimal example using **one simulated dataset** (no
Monte Carlo repetition). We reproduce a simple dual-frame scenario
similar to that in the paper.

The data-generating mechanism is provided as an example script in

`inst/examples/dualframe_simulation.R`

which defines the function `generate_dualframe_population()`.

``` r
library(dfSEDI)

## Load the example script shipped with the package
## This script defines `generate_dualframe_population()`
example_file <- system.file("examples", "dualframe_simulation.R",
                            package = "dfSEDI")
source(example_file)

set.seed(1)

## 1. Generate one artificial dual-frame population
N <- 10000
dat <- generate_dualframe_population(N = N)

## The data frame `dat` must contain at least:
##   - X   : matrix (or AsIs-matrix) of covariates (n x p)
##   - y   : outcome variable
##   - d_np: indicator of inclusion in the non-probability sample
##   - d_p : indicator of inclusion in the probability sample
##   - pi_p: design inclusion probability for the P-sample
##   - pi_np: true inclusion probability for the NP-sample (for simulations)

## 2. Semiparametric efficient estimator Eff (DML2, K-fold)
fit_eff <- Eff(
  dat         = dat,
  K           = 2,    # DML2 with 2-fold cross-fitting
  phi_start   = NULL, # default: intercept = logit(mean(d_np)), other components = 0
  max_restart = 10,
  progress    = TRUE  # show simple progress in the console
)

fit_eff$theta  # point estimate
fit_eff$se     # sandwich standard error
fit_eff$ci     # 95% confidence interval

## 3. Sub-efficient estimator Eff_S (Remark 6)
fit_effS <- Eff_S(
  dat      = dat,
  K        = 2,
  progress = TRUE
)

fit_effS$theta
fit_effS$se
fit_effS$ci

## 4. Parametric efficient estimator Eff_P (working model)
##    Here we let the function choose default starting values.
fit_effP <- Eff_P(
  dat       = dat,
  phi_start = NULL,  # default starting values
  eta4_star = 0,
  max_iter  = 20,
  progress  = TRUE
)

fit_effP$theta
fit_effP$se
fit_effP$ci
```

Each of these objects has the following structure:

- `theta` : point estimate of the target parameter (e.g., population
  mean)
- `var` : sandwich variance estimate
- `se` : sandwich standard error (sqrt of `var`)
- `ci` : 95% sandwich-based confidence interval (numeric vector of
  length 2)
- `phi` : estimated parameter vector for the non-probability inclusion
  model (for `Eff` and `Eff_P`)
- `info` : list of meta-information (e.g., K, starting values, estimator
  type)

This interface allows users to call each estimator separately and obtain
its point estimate, sandwich standard error, and confidence interval on
a single dual-frame dataset, without running any Monte Carlo
experiments.
