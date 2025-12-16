dfSEDI
================

# dfSEDI

dfSEDI implements semiparametric efficient data integration methods for
dual-frame sampling (Morikawa & Kim, 202x).

Main user-facing estimators:

- `Eff` : semiparametric efficient estimator (supports DML1 / DML2,
  K-fold)
- `Eff_S` : sub-efficient estimator (Remark 6-type, K-fold)
- `Eff_P` : parametric efficient estimator (working model)

For comparison (mainly for simulations), the package also includes:

- `df_estimate_P` : probability-only estimator (HT-type)
- `df_estimate_NP` : non-probability-only estimator (Chang & Kott-type)
- `df_estimate_NP_P` : NP union P estimator (Chang & Kott-type)

All estimators return a point estimate for $\theta = E(Y)$, and
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

All core functions assume a data frame `dat` that contains at least:

- `X` : covariates used in the *parametric* NP inclusion model and also
  as inputs to nuisance models  
  - can be a **numeric matrix** stored as a single column,
  e.g. `X = I(matrix(...))`, or  
  - a **data.frame** (recommended if you have categorical variables),
  which will be converted internally  
  - alternatively, a numeric vector `x` can be provided (1D covariate
  case)
- `y` : outcome variable
- `d_np` : indicator for inclusion in the non-probability sample (0/1)
- `d_p` : indicator for inclusion in the probability sample (0/1)
- `pi_p` : design inclusion probability for the probability sample
- `pi_np` : (optional; simulations only) true NP inclusion probability

### Mixed continuous/discrete covariates (automatic)

dfSEDI estimates nuisance functions (e.g., $\mu(X)=E[Y\mid X]$, $\pi_P$,
$\eta_4^*(X)$, $h_4^*(X)$) using kernel ridge regression via
`kernlab::gausspr`.

If `X` contains **categorical/discrete** columns, dfSEDI automatically
uses a **product kernel**:

- RBF (Gaussian) kernel on continuous columns
- Delta (Kronecker) kernel on categorical columns (only same-category
  pairs interact)

This is equivalent to fitting separate KRR models *within each category
combination*, and it often yields **large speedups** because the kernel
matrix becomes effectively block-diagonal.

#### How categorical columns are detected

- If `dat$X` is a **data.frame**, columns of type `factor`, `character`,
  or `logical` are treated as categorical.
- Numeric columns with “few” unique values may also be treated as
  categorical by a heuristic.

**Recommended practice:** supply `X` as a `data.frame` and explicitly
mark categorical variables as `factor`.

> Note (parametric part): In the NP inclusion model for $\phi$, `X`
> enters linearly as numeric columns.  
> If you have multi-level categorical variables and want a proper
> parametric specification, dummy-encode them yourself (e.g., with
> `model.matrix`) before passing them as `X`.

------------------------------------------------------------------------

## Example: one simulated dataset (no Monte Carlo)

The package ships an example script in `inst/examples`:

- `inst/examples/dualframe_simulation.R`

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

If you want `X` to contain an explicit categorical column (recommended
in real data), you can do:

``` r
dat$X <- data.frame(
  x1 = dat$X[, 1],
  x2 = factor(dat$X[, 2])
)
```

------------------------------------------------------------------------

## Eff: DML1

DML1 conceptually computes fold-specific estimates and aggregates them:

- split data into folds $S_1, \ldots, S_K$
- for each fold $k$:
  - estimate nuisances (including $\pi_p$) on training data (all folds
    except $k$)
  - compute fold-specific estimates on test fold $k$
- aggregate fold-specific estimates across $k$

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

## Eff: DML2

- `type = 1` (or `"DML1"`) : DML1
- `type = 2` (or `"DML2"`) : DML2 (default)

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

Eff returns joint sandwich-based standard errors for both $\theta$ and
$\phi$ (`phi_se`, `phi_ci`).

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

`Eff_S` is based on a sub-efficient influence function that only
involves the outcome regression $E(Y \mid X)$.

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

`Eff_P` assumes parametric working models for both the sampling
mechanism and the outcome regression. It returns point estimates and
standard errors for $\theta$; standard errors for $\phi$ are not
currently implemented in `Eff_P`.

------------------------------------------------------------------------

## Returned object structure

`Eff`, `Eff_S`, and `Eff_P` return a list with (at least):

- `theta` : point estimate of $\theta = E(Y)$
- `var` : sandwich variance estimate for `theta`
- `se` : sandwich standard error for `theta`
- `ci` : 95% confidence interval for `theta`

For the propensity parameter $\phi$:

- `Eff` (DML1 / DML2) returns
  - `phi` : estimate of the NP propensity parameter
  - `phi_var` : joint sandwich variance matrix for `phi`
  - `phi_se` : standard errors for each component of `phi`
  - `phi_ci` : 95% confidence intervals for each component of `phi`
- `Eff_P` returns `phi` but currently leaves `phi_var`, `phi_se`, and
  `phi_ci` as `NULL`.

The `info` element stores meta information (K, type, convergence code,
etc.) that can be useful for diagnostics.
