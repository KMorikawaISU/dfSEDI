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

If you plan to use mixed continuous/discrete covariates via `x_cont` /
`x_disc` (Li & Racine-type regression), you also need the
[`np`](https://cran.r-project.org/package=np) package:

``` r
install.packages("np")
```

Load the package:

``` r
library(dfSEDI)
```

------------------------------------------------------------------------

## Data format

All core functions assume a data frame `dat` that contains at least:

- `X` : covariates used in the parametric NP propensity model and some
  nuisance models  
  (either a matrix stored as a single column, e.g. `X = I(matrix(...))`,
  or a data frame that can be converted to a numeric matrix).  
  Alternatively, a numeric vector `x` can be provided (1D covariate
  case).  
  If neither `X` nor `x` is present, `df_get_X()` will construct `X`
  from `x_cont` / `x_disc` (see below) by coercing them to numeric.
- `y` : outcome variable
- `d_np` : indicator for inclusion in the non-probability sample (0/1)
- `d_p` : indicator for inclusion in the probability sample (0/1)
- `pi_p` : design inclusion probability for the probability sample
- `pi_np` : (optional; simulations only) true NP inclusion probability

Optionally, you can split the covariates used in the **nonparametric**
nuisance models into continuous and discrete parts:

- `x_cont` : continuous covariates for Li & Racine–type regression  
  (vector, matrix, or data frame; converted to numeric internally)
- `x_disc` : discrete / categorical covariates  
  (vector, matrix, or data frame; columns are converted to `factor`)

If `x_disc` is present, dfSEDI uses Li & Racine (2007) mixed-type kernel
regression for the nuisance functions via the
[`np`](https://cran.r-project.org/package=np) package.  
If `x_disc` is *not* present, dfSEDI falls back to Gaussian kernels via
`kernlab::gausspr` (continuous-only case), which matches older versions
of the package.

Internally, `df_get_X(dat)` constructs the parametric design matrix as
follows:

- if `dat$X` exists, it is used (after coercion to a numeric matrix);
- otherwise, if `dat$x` exists, it is treated as a 1-column matrix;
- otherwise, if `x_cont` and/or `x_disc` are present, they are combined
  into a numeric matrix (discrete covariates are coerced to numeric
  codes).

For nonparametric regression,
`df_get_np_design(dat, include_y = FALSE/TRUE)` builds a `data.frame`
with continuous variables as numeric and discrete variables as factors,
which is then passed to `np::npregbw()` / `np::npreg()`.

------------------------------------------------------------------------

## Example: one simulated dataset (no Monte Carlo)

The package ships an example script in `inst/examples`:

- `inst/examples/dualframe_simulation.R`

It defines:

- `generate_dualframe_population(N)`

where the covariates satisfy:

- `x1` : continuous covariate (stored in `x_cont`)
- `x2` : binary covariate (stored in `x_disc` and also included in `X`)

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

## Eff: DML1 (recommended)

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

## Eff: DML2 (takes more time than DML1)

DML2 conceptually estimates nuisances on the full data and computes a
single estimate.

- `type = 1` (or `"DML1"`) : DML1
- `type = 2` (or `"DML2"`) : DML2

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
