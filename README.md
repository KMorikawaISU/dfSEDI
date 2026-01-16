dfSEDI
================

# dfSEDI

dfSEDI implements semiparametric efficient data integration methods for
dual-frame sampling (Morikawa & Kim, 202x).

Main user-facing estimators:

- `Eff` : semiparametric efficient estimator (**default: DML1**,
  supports DML1 / DML2, K-fold)
- `Eff_S` : sub-efficient estimator (Remark 6-type, K-fold)
- `Eff_P` : parametric efficient estimator (working model)

For comparison (mainly for simulations), the package also includes:

- `df_estimate_P` : probability-only estimator (HT-type)
- `df_estimate_NP` : non-probability-only estimator (Chang & Kott-type;
  **requires `base_fun`**)
- `df_estimate_NP_P` : NP union P estimator (Chang & Kott-type;
  **requires `base_fun`**)

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
  - a **data.frame** (recommended if you have categorical variables)  
  - alternatively, a numeric vector `x` can be provided (1D covariate
  case)
- `y` : outcome variable (continuous or binary)
- `d_np` : indicator for inclusion in the non-probability sample (0/1)
- `d_p` : indicator for inclusion in the probability sample (0/1)
- `pi_p` : design inclusion probability for the probability sample
- `pi_np` : (optional; simulations only) true NP inclusion probability

------------------------------------------------------------------------

## `pi_x_cols` (use only a subset of X in the NP response model)

By default, dfSEDI uses **all** columns of `X` in the parametric
non-probability inclusion (response-mechanism) model:

$$
\pi_{np}(X,Y) = \text{expit}\left(\phi^\top (1, X, Y)\right).
$$

In some applications, you may want to estimate the response mechanism
using **only a subset of** `X` (e.g., a small set of key design
variables), while still using the **full** `X` for nuisance regressions
such as $\kappa(X)$ and $\eta(X)$. You can do this by passing
`pi_x_cols` to `Eff()` (and also to `Eff_P()`, `df_estimate_NP()`,
`df_estimate_NP_P()` if needed).

- `pi_x_cols = NULL` (default): use **all** columns of `X` (same
  behavior as earlier versions).
- `pi_x_cols` can be:
  - an **integer** vector of column indices (1-based),
  - a **logical** vector of length `ncol(df_get_X(dat))`,
  - a **character** vector:
    - exact column names of the internal design matrix `df_get_X(dat)`,
      **or**
    - (when `dat$X` is a `data.frame`) original term/variable names like
      `"METRO"`; dfSEDI will include **all dummy columns** created for
      that term.

Important:

- `pi_x_cols` affects **only** the construction of $(1, X_{\pi}, Y)$
  inside $\pi_{np}$ and the score terms involving $\phi$.
- It does **not** restrict the covariates used in the later nuisance
  regressions; those still use the full `X` by default.

Example:

``` r
# Suppose dat$X is a data.frame that includes METRO, REGION9, and race variables:
fit <- Eff(dat, K = 2, type = 2, pi_x_cols = c("METRO", "REGION9", "q10_1"))

# You can inspect what columns were actually used for the pi_np model:
fit$info$pi_x_cols_used
```

### Important: handling the $(d_{np}, d_p) = (0,0)$ pattern (`x_info`)

In the paper’s notation, there are four response patterns
$(\delta_{NP}, \delta_P)\in\{0,1\}^2$. In many *practical*
data-integration tasks you only have units that appear in *either* the
probability sample or the non-probability sample, i.e. you only observe
the union sample $\delta_{NP}\cup\delta_P = 1$. In that common case,
there are **no rows** with $(d_{np}, d_p)=(0,0)$.

The paper also notes that the theory covers the case where $X$ is
unavailable when $(\delta_{NP},\delta_P)=(0,0)$; in that case, it
suffices to set the augmentation functions $\eta(X)$ and $\kappa(X)$ to
zero.

dfSEDI exposes this choice via the flag `x_info`:

- `x_info = FALSE` (**recommended for most real-data use cases**):  
  Use the simplified formulation with $\eta(X)=\kappa(X)=0$.  
  This mode is intended for the common situation where you only have the
  union sample, or you do not have covariates $X$ for $(0,0)$ units. It
  also avoids computations involving the $(0,0)$ pattern and can be
  substantially faster.

- `x_info = TRUE`:  
  Use the full formulation as written in the main text of the paper.  
  This is appropriate if you explicitly include $(0,0)$ units in `dat`
  **and** you have $X$ for those units (e.g., simulation studies, or
  settings with a population register providing $X$).

> Note on population size `N` when `x_info = FALSE`:
>
> If you pass a **union-sample-only** dataset to
> `Eff(..., x_info = FALSE)` (i.e., `dat` has no $(d_{np}, d_p)=(0,0)$
> rows), dfSEDI can still target the full-frame mean $\theta = E(Y)$
> **if you provide the population size** via `N = ...`. Internally, the
> estimator is computed on a population-total scale and then rescaled by
> `nrow(dat) / N`.
>
> If you omit `N`, `Eff(..., x_info = FALSE)` returns the *un-rescaled*
> average over the observed union sample. In that case, you can recover
> the population-mean scale by multiplying `theta`, `se`, and `ci` by
> `nrow(dat_union) / N`.

> Practical tip: If your input data include only sampled units (union
> sample), you can create it as:
>
> ``` r
> dat_union <- subset(dat, d_np == 1 | d_p == 1)
> ```

------------------------------------------------------------------------

## `prob_only` in `Eff()` (estimating the augmentation terms)

When `x_info = TRUE`, `Eff()` estimates the optimal augmentation
functions $\eta(X)$ and $\kappa(X)$ via cross-fitted nonparametric
regression (DML).

In the paper these nuisance functions admit two equivalent
representations (see Eq. (9)/(11) and Eq. (10)/(12) in the current
draft):

- **Probability-sample-only** representation: conditional expectations
  are taken over the probability sample ($\delta_P = 1$), involving
  $O_{NP\cup P} / \pi_P$.
- **Combined-sample** representation: conditional expectations are taken
  over the union sample ($\delta_{NP}\cup\delta_P = 1$), involving
  $O_{NP\cup P} / \pi_{NP\cup P}$.

dfSEDI exposes this choice via `prob_only`:

- `prob_only = FALSE` (default): use the **combined-sample** version
  (10)/(12).
- `prob_only = TRUE`: use the **probability-sample-only** version
  (9)/(11).

Notes:

- `prob_only` only affects the estimation of $\eta(X)$ and $\kappa(X)$.
  The target parameter and the main estimating equations are unchanged.
- If `x_info = FALSE`, the augmentation terms are fixed to zero and
  `prob_only` is ignored.
- Trade-off: `prob_only = TRUE` uses fewer training points (only the
  probability sample), so it can be noisier when $n_P$ is small;
  `prob_only = FALSE` uses more training points (union sample), but
  depends more on the combined-sample weights.

------------------------------------------------------------------------

## `nonpara_method` (nuisance regression back-end)

dfSEDI estimates several nuisance functions by regression (e.g.,
$\mu(X)=E[Y\mid X]$, $\bar\pi_P(L)$ imputation, and the DML nuisance
terms $\eta(X)$ / $\kappa(X)$). You can now choose the regression
back-end via `nonpara_method`:

- `nonpara_method = "KRR"` (**default**): kernel ridge regression with
  an RBF/Gaussian kernel applied to **all** columns of the design
  matrix.
- `nonpara_method = "mixed KRR"`: mixed-type KRR with a product kernel:
  - RBF kernel on continuous columns
  - delta (Kronecker) kernel on categorical columns (fits separate
    models within each categorical cell)
- `nonpara_method = "logistic KRR"`: kernel logistic regression (used to
  estimate $P(Y=1\mid X)$ when `y` is binary).
- `nonpara_method = "glmnet_linear"`: elastic-net regression via
  `glmnet` with Gaussian family.
- `nonpara_method = "glmnet_logistic"`: elastic-net logistic regression
  via `glmnet` with binomial family (used for $P(Y=1\mid X)$ when `y` is
  binary).
- `nonpara_method = "RF_cont"`: random forest **regression** (continuous
  outcome), using `ranger` if available (fast) and falling back to
  `randomForest` if not.
- `nonpara_method = "RF_binom"`: random forest **classification** for
  $P(Y=1\mid X)$, using `ranger` if available and falling back to
  `randomForest` if not.

Important notes:

- `glmnet_*` requires the `glmnet` package
  (`install.packages("glmnet")`).
- `RF_*` requires either `ranger` (**recommended**) or `randomForest`
  (`install.packages("ranger")`).
- For **continuous** regression tasks, the binomial-only methods
  (`"logistic KRR"`, `"glmnet_logistic"`, `"RF_binom"`) automatically
  fall back to their continuous counterparts with a warning.
- For **binary** `y` in `Eff()`, selecting `"logistic KRR"`,
  `"glmnet_logistic"`, or `"RF_binom"` uses the **binary-mixture**
  strategy:
  1)  estimate $p(X)=P(Y=1\mid X)$, then  
  2)  compute the conditional expectations by mixing the $Y=0$ and $Y=1$
      expressions.
- For **binary** `y` in `Eff()`, selecting `"KRR"`, `"mixed KRR"`,
  `"glmnet_linear"`, or `"RF_cont"` uses the **direct** strategy (treat
  `y` as numeric 0/1) and directly regresses the required
  numerator/denominator pseudo-outcomes.

Example usage:

``` r
# Efficient estimator (DML): choose the nuisance regression back-end
fit_eff_krr     <- Eff(dat, K = 2, x_info = TRUE, prob_only = FALSE, nonpara_method = "KRR")
fit_eff_mixed   <- Eff(dat, K = 2, x_info = TRUE, prob_only = FALSE, nonpara_method = "mixed KRR")
fit_eff_glmnet  <- Eff(dat, K = 2, x_info = TRUE, prob_only = FALSE, nonpara_method = "glmnet_linear")
fit_eff_rf      <- Eff(dat, K = 2, x_info = TRUE, prob_only = FALSE, nonpara_method = "RF_cont")

# If y is binary and you want to use a probability model P(Y=1|X):
fit_eff_logit   <- Eff(dat, K = 2, x_info = TRUE, prob_only = FALSE, nonpara_method = "logistic KRR")
fit_eff_glmnetL <- Eff(dat, K = 2, x_info = TRUE, prob_only = FALSE, nonpara_method = "glmnet_logistic")
fit_eff_rfL     <- Eff(dat, K = 2, x_info = TRUE, prob_only = FALSE, nonpara_method = "RF_binom")
```

## Mixed continuous/discrete covariates (automatic)

> Note: The mixed-type (RBF × delta) product-kernel behavior below is
> used when you set `nonpara_method = "mixed KRR"`. If you set
> `nonpara_method = "KRR"` (plain RBF on all columns) or `"glmnet_*"` /
> `"RF_*"`, dfSEDI will not use the delta kernel.

dfSEDI estimates nuisance functions (e.g., $\mu(X)=E[Y\mid X]$, $\pi_P$,
$\kappa(X)$, $\eta(X)$) using kernel methods.

- If `X` is continuous-only, dfSEDI uses Gaussian/RBF kernels (via
  `kernlab::gausspr`).
- If `X` contains **categorical/discrete** columns, dfSEDI automatically
  uses a **product kernel**:
  - RBF (Gaussian) kernel on continuous columns
  - Delta (Kronecker) kernel on categorical columns (only same-category
    pairs interact)

This is equivalent to fitting separate kernel models *within each
category combination*, and it often yields **large speedups** because
the kernel matrix becomes effectively block-diagonal.

#### How categorical columns are detected

- If `dat$X` is a **data.frame**, columns of type `factor`, `character`,
  or `logical` are treated as categorical.
- Numeric columns with “few” unique values may also be treated as
  categorical by a heuristic.

**Recommended practice:** supply `X` as a `data.frame` and explicitly
mark categorical variables as `factor`.

> Note (parametric part for $\phi$): If `X` is a `data.frame` with
> factor columns, dfSEDI internally uses dummy encoding (`model.matrix`)
> for the NP inclusion model so multi-level categorical variables enter
> as dummies automatically.

#### Small-cell warning (categorical cells)

With a delta kernel, very small categorical cells can make nuisance
estimates unstable. dfSEDI can warn (or ask interactively) when cell
sizes are too small. You can control this behavior via options, e.g.:

``` r
options(dfSEDI.min_cell_size = 30)                 # threshold
options(dfSEDI.small_cell_action = "ask")          # "ask" / "warn" / "stop" / "none"
options(dfSEDI.allow_small_cells = FALSE)          # bypass flag
```

------------------------------------------------------------------------

## Chang & Kott-type NP estimators: specifying `base_fun`

`df_estimate_NP()` and `df_estimate_NP_P()` solve a Chang & Kott-type
estimating equation for the NP propensity parameter $\phi$ using a
user-supplied basis function:

- `base_fun` must be a function of the form
  `base_fun <- function(X) { ... }`
- `X` is the $n \times p$ design matrix extracted from `dat` (via
  `df_get_X(dat)`)
- `base_fun(X)` must return a **numeric matrix** with shape:

$$
B = \texttt{base\\_fun}(X) \in \mathbb{R}^{n \times (p+2)}.
$$

dfSEDI performs strict checks for safety:

- `nrow(base_fun(X))` must equal `nrow(X)`
- `ncol(base_fun(X))` must equal `ncol(X) + 2`
- all entries must be finite numeric values (no `NA`, `Inf`, `-Inf`)

If any check fails, the estimator stops with an error (rather than
returning silent nonsense).

### Example `base_fun` choices

Quadratic augmentation (example):

``` r
base_fun <- function(X) {
  cbind(1, X, X[, 1]^2)
}
```

> Tip: If you use a closure like this, make sure `base_fun` is defined
> for the same `dat` you pass into `df_estimate_NP()` /
> `df_estimate_NP_P()`.

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
dat <- generate_dualframe_population(N = N, Scenario = 1)

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

### Example A (typical): union sample only (`x_info = FALSE`)

If you want to mimic the most common real-data situation where you only
observe the union sample:

``` r
dat_union <- subset(dat, d_np == 1 | d_p == 1)

# Eff uses DML1 by default (faster).
fit_eff_union <- Eff(
  dat         = dat_union,
  K           = 2,
  x_info      = FALSE,
  N           = N,      # population size (frame size)
  progress    = TRUE
)

fit_eff_union$theta
fit_eff_union$se
fit_eff_union$ci
fit_eff_union$info
```

### Example B (simulation / full info): include $(0,0)$ rows (`x_info = TRUE`)

For simulation studies where you keep the full population-like data
structure:

``` r
# Eff uses DML1 by default (faster).
fit_eff_full <- Eff(
  dat         = dat,
  K           = 2,
  x_info      = TRUE,
  progress    = TRUE
)

fit_eff_full$theta
fit_eff_full$se
fit_eff_full$ci
fit_eff_full$info
```

------------------------------------------------------------------------

## Monte Carlo simulation (run_mc)

The example script `inst/examples/dualframe_simulation.R` also defines:

- `run_mc(B, N, Scenario, K, seed_start, show_progress, progress_each_fit, pi_p_offset, ...)`
- `summarize_mc(res, theta_true = 0)`

`run_mc()` returns a long-format data frame (one row per replication ×
estimator). With the current example script, the `estimator` column
takes values:

- `P`
- `NP`
- `NP_P`
- `Eff1` (Eff with `type = 1`; DML1)
- `Eff2` (Eff with `type = 2`; DML2)
- `Eff1_union` (Eff on `dat_union` with `type = 1`; DML1)
- `Eff2_union` (Eff on `dat_union` with `type = 2`; DML2)
- `Eff_S`
- `Eff_P`

The typical columns are:

- `Scenario`, `rep`, `estimator`
- `theta`, `se`, `ci_l`, `ci_u`
- `n_np`, `n_p`, `n_union`
- `error`

### Serial MC (built-in progress bar)

``` r
library(dfSEDI)

example_file <- system.file("examples", "dualframe_simulation.R", package = "dfSEDI")
source(example_file)

res <- run_mc(
  B = 100,
  N = 10000,
  Scenario = 1,
  K = 2,
  seed_start = 1,
  show_progress = TRUE,
  progress_each_fit = FALSE
)

summarize_mc(res, theta_true = 0)
```

### Parallel MC with a progress bar (Windows)

Windows does not support fork-based parallelism (`mclapply`). Use a
PSOCK cluster with `pbapply`:

Needed packages (Suggests): `pbapply` (and base R `parallel`)

``` r
# install.packages("pbapply")  # if needed

library(dfSEDI)
library(parallel)
library(pbapply)

example_file <- system.file("examples", "dualframe_simulation.R", package = "dfSEDI")
source(example_file)

B <- 200
N <- 10000
Scenario <- 1
K <- 2

seeds <- 1 + seq_len(B) - 1

n_workers <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_workers)
on.exit(parallel::stopCluster(cl), add = TRUE)

pbapply::pboptions(type = "txt")

parallel::clusterEvalQ(cl, library(dfSEDI))

# Export the functions defined by the example script
parallel::clusterExport(
  cl,
  varlist = c(
    "normalize_scenario",
    "generate_dualframe_population",
    "make_base_fun",
    "safe_fit",
    "extract_row",
    "fit_all_estimators_once",
    "N",
    "Scenario",
    "K"
  ),
  envir = environment()
)

one_rep <- function(seed) {
  set.seed(seed)
  dat <- generate_dualframe_population(N = N, Scenario = Scenario)
  fits <- fit_all_estimators_once(dat = dat, Scenario = Scenario, K = K, progress_each = FALSE)
  rbind(
    extract_row(fits$P,          "P",          seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP,         "NP",         seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP_P,       "NP_P",       seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_S,      "Eff_S",      seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_P,      "Eff_P",      seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff1_union, "Eff1_union", seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff2_union, "Eff2_union", seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff1,       "Eff1",       seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff2,       "Eff2",       seed, Scenario, fits$n_np, fits$n_p, fits$n_union)
  )
}

out_list <- pbapply::pblapply(seeds, one_rep, cl = cl)
res_par <- do.call(rbind, out_list)

summarize_mc(res_par, theta_true = 0)
```

### Parallel MC with a progress bar (macOS / Linux)

We recommend using the same PSOCK + `pbapply` approach as in the Windows
section. This is more robust than fork-based parallelism and avoids the
common `scheduled cores ... did not deliver results` issue.

Needed packages (Suggests): `pbapply` (and base R `parallel`)

``` r
# install.packages("pbapply")  # if needed

library(dfSEDI)
library(parallel)
library(pbapply)

example_file <- system.file("examples", "dualframe_simulation.R", package = "dfSEDI")
source(example_file)

B <- 200
N <- 10000
Scenario <- 1
K <- 2

seeds <- 1 + seq_len(B) - 1

n_workers <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_workers)
on.exit(parallel::stopCluster(cl), add = TRUE)

pbapply::pboptions(type = "txt")

parallel::clusterEvalQ(cl, library(dfSEDI))

parallel::clusterExport(
  cl,
  varlist = c(
    "normalize_scenario",
    "generate_dualframe_population",
    "make_base_fun",
    "safe_fit",
    "extract_row",
    "fit_all_estimators_once",
    "N",
    "Scenario",
    "K"
  ),
  envir = environment()
)

one_rep <- function(seed) {
  set.seed(seed)
  dat <- generate_dualframe_population(N = N, Scenario = Scenario)
  fits <- fit_all_estimators_once(dat = dat, Scenario = Scenario, K = K, progress_each = FALSE)
  rbind(
    extract_row(fits$P,          "P",          seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP,         "NP",         seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP_P,       "NP_P",       seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_S,      "Eff_S",      seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_P,      "Eff_P",      seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff1_union, "Eff1_union", seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff2_union, "Eff2_union", seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff1,       "Eff1",       seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff2,       "Eff2",       seed, Scenario, fits$n_np, fits$n_p, fits$n_union)
  )
}

out_list <- pbapply::pblapply(seeds, one_rep, cl = cl)
res_par <- do.call(rbind, out_list)

summarize_mc(res_par, theta_true = 0)
```

Notes: - In the parallel examples above, we run one replication per
seed, and the progress bar is shown in the master process. - If you want
to list optional packages in DESCRIPTION, add: - `pbapply` (progress bar
for PSOCK clusters; works on Windows/macOS/Linux) - `glmnet` (needed
when you use `nonpara_method = "glmnet_*"`) - `ranger` (recommended;
needed when you use `nonpara_method = "RF_*"`) - (optional fallback)
`randomForest` can also be used if `ranger` is not installed
