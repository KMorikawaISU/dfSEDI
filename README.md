dfSEDI
================

# dfSEDI

dfSEDI implements semiparametric efficient data integration methods for
dual-frame sampling (Morikawa & Kim, 202x).

Main user-facing estimators:

- Eff : semiparametric efficient estimator (supports DML1 / DML2,
  K-fold)
- Eff_S : sub-efficient estimator (Remark 6-type, K-fold)
- Eff_P : parametric efficient estimator (working model)

For comparison (mainly for simulations), the package also includes:

- df_estimate_P : probability-only estimator (HT-type)
- df_estimate_NP : non-probability-only estimator (Chang & Kott-type)
- df_estimate_NP_P : NP union P estimator (Chang & Kott-type)

All estimators return a point estimate for theta = E(Y), and
sandwich-type SE/CI computed from influence-function / pseudo-outcome
contributions.

------------------------------------------------------------------------

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("KMorikawaISU/dfSEDI")
```

Load the package:

``` r
library(dfSEDI)
```

------------------------------------------------------------------------

## Data format

All core functions assume a data frame dat that contains at least:

- X : covariates used in the NP propensity and outcome models (either a
  matrix stored as a single column, e.g. X = I(matrix(…)), or a data
  frame that can be converted to a numeric matrix) Alternatively, a
  numeric vector x can be provided (1D covariate case).
- y : outcome variable
- d_np : indicator for inclusion in the non-probability sample (0/1)
- d_p : indicator for inclusion in the probability sample (0/1)
- pi_p : design inclusion probability for the probability sample
- pi_np : (optional; simulations only) true NP inclusion probability

Internally, dfSEDI extracts the covariate matrix via a helper: - if
dat$X exists, it is used (after coercion to numeric matrix)
- otherwise, if dat$x exists, it is treated as a 1-column matrix

------------------------------------------------------------------------

## Example: one simulated dataset (no Monte Carlo)

The package ships an example script in inst/examples:

- inst/examples/dualframe_simulation.R

It defines: - generate_dualframe_population(N)

Example:

``` r
library(dfSEDI)

# Load the bundled example script (defines generate_dualframe_population)
example_file <- system.file("examples", "dualframe_simulation.R", package = "dfSEDI")
source(example_file)

set.seed(18)
N <- 10000
dat <- generate_dualframe_population(N = N)

str(dat)
```

------------------------------------------------------------------------

## Eff: DML1 (Recommend)

DML1 conceptually computes fold-specific estimates and aggregates them:

- split data into folds S_1, …, S_K
- for each fold k:
  - estimate nuisances (including pi_p) on training data (-k)
  - compute fold-specific estimates on test fold k
- aggregate fold-specific estimates across k

DML1 example:

``` r
fit_eff_dml1 <- Eff(
  dat         = dat,
  K           = 3,
  type        = 1,      # DML1
  phi_start   = NULL,
  max_restart = 10,
  progress    = TRUE
)

fit_eff_dml1$phi
fit_eff_dml1$theta
fit_eff_dml1$se
fit_eff_dml1$ci
fit_eff_dml1$info
```

------------------------------------------------------------------------

## Eff: DML2 (takes much more time than DML1)

DML2 conceptually estimates nuisances on the full data and computes a
single estimate.

- type = 1 (or “DML1”) : DML1
- type = 2 (or “DML2”) : DML2

DML2 example:

``` r
fit_eff_dml2 <- Eff(
  dat         = dat,
  K           = 2,
  type        = 2,      # DML2
  phi_start   = NULL,
  max_restart = 10,
  progress    = TRUE
)

fit_eff_dml2$phi
fit_eff_dml2$theta
fit_eff_dml2$se
fit_eff_dml2$ci
fit_eff_dml2$info
```

Note: the current version focuses on inference for theta. Standard
errors for phi are not provided (phi_se may be NULL).

------------------------------------------------------------------------

## Eff_S (sub-efficient)

``` r
fit_effS <- Eff_S(
  dat      = dat,
  K        = 2,
  progress = TRUE
)

fit_effS$theta
fit_effS$se
fit_effS$ci
fit_effS$info
```

------------------------------------------------------------------------

## Eff_P (parametric working model)

``` r
fit_effP <- Eff_P(
  dat       = dat,
  phi_start = NULL,
  eta4_star = 0,
  max_iter  = 20,
  progress  = TRUE
)

fit_effP$phi
fit_effP$theta
fit_effP$se
fit_effP$ci
fit_effP$info
```

------------------------------------------------------------------------

## Returned object structure

Eff / Eff_S / Eff_P return a list with (at least):

- theta : point estimate of theta = E(Y)
- var : sandwich variance estimate for theta
- se : sandwich standard error for theta
- ci : 95% confidence interval for theta
- phi : estimated NP propensity parameter vector (Eff and Eff_P)
- info : meta information (K, type, convergence, etc.)

phi standard errors are not currently included (phi_se may be NULL).
