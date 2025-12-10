dfSEDI
================

# dfSEDI

The **dfSEDI** package implements semiparametric efficient data
integration methods for dual-frame sampling as described in Morikawa &
Kim (202x).

The package provides the following estimators:

- `Eff` : semiparametric efficient estimator (DML2, K-fold)
- `Eff_S` : sub-efficient estimator (based on Remark 6)
- `Eff_P` : parametric efficient estimator (working model)

Internally, the package also implements basic estimators such as

- `P` : probability-only estimator
- `NP` : non-probability-only estimator
- `NP_P`: combined NP ∪ P estimator

For each estimator, sandwich-based variance, standard error, and 95%
confidence intervals are computed from influence-function-based
contributions.

## Installation

You can install the development version from GitHub:

    # install.packages("devtools")
    devtools::install_github("KMorikawaISU/dfSEDI")

Replace `"yourname"` by your GitHub user name.

## Example: one simulated dataset

This section shows a minimal example using **one simulated dataset** (no
repetition). We reproduce the scenario \[O2\] × \[NP1\] × P used in the
paper.

We first generate an artificial dual-frame population and then compute
the three main estimators: `Eff`, `Eff_S`, and `Eff_P`.

    library(dfSEDI)

    set.seed(1)

    # 1. Generate one artificial dual-frame population
    N <- 10000
    dat <- generate_dualframe_population(N = N)

    # dat must contain at least:
    #   x, y, d_np, d_p, pi_p, pi_np

    # 2. Semiparametric efficient estimator Eff (DML2, K-fold)
    fit_eff <- Eff(
      dat       = dat,
      K         = 2,                      # DML2 with 2 folds
      phi_start = c(-2.15, -0.5, -0.75),  # initial values for pi_NP
      max_restart = 10
    )

    fit_eff$theta  # point estimate
    fit_eff$se     # sandwich standard error
    fit_eff$ci     # 95% confidence interval

    # 3. Sub-efficient estimator Eff_S (Remark 6)
    fit_effS <- Eff_S(
      dat = dat,
      K   = 2
    )

    fit_effS$theta
    fit_effS$se
    fit_effS$ci

    # 4. Parametric efficient estimator Eff_P (working model)
    fit_effP <- Eff_P(
      dat       = dat,
      phi_start = c(-2.15, -0.5, -0.75),
      eta4_star = 0
    )

    fit_effP$theta
    fit_effP$se
    fit_effP$ci

Each of these objects has the following structure:

- `theta` : point estimate
- `var` : sandwich variance estimate
- `se` : sandwich standard error (sqrt of `var`)
- `ci` : 95% sandwich-based confidence interval (numeric vector of
  length 2)
- `phi` : estimated parameter vector for the non-probability inclusion
  model (Eff and Eff_P only)
- `info` : list of meta-information (e.g., K, starting values, estimator
  type)

This interface allows users to call each estimator separately and obtain
its point estimate, sandwich standard error, and confidence interval
without running any Monte Carlo experiments.
