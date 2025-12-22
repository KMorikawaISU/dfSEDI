############################################################
## dualframe_simulation.R
##
## Simulation runner for dfSEDI dual-frame estimators
##
## Observation pattern (fixed; matches Table 1 / your intended setup):
##   - y is observed on the union sample only:
##       (d_np, d_p) = (0,0)  ->  y = NA   (and y_true is kept)
##       otherwise            ->  y observed
##   - pi_p is observed only when d_p == 1:
##       d_p == 0             ->  pi_p = NA (and pi_p_true is kept)
##       d_p == 1             ->  pi_p observed
##
## Notes:
##   - The core code (dualframe_core.R) replaces y=NA on (0,0) units with 0 internally so
##     score contributions stay finite and (0,0) units are not dropped from means.
##   - Eff / Eff_P / Eff_S only *impute* pi_p when pi_p is NA (they do not overwrite if pi_p
##     is already filled with finite values).
############################################################

suppressPackageStartupMessages({
  if (requireNamespace("dfSEDI", quietly = TRUE)) {
    library(dfSEDI)
  }
})

# If you are NOT using the dfSEDI package and instead source dualframe_core.R,
# ensure Eff / Eff_S / Eff_P exist.
if (!exists("Eff", mode = "function")) {
  warning("Eff() not found. Did you source dualframe_core.R (or load dfSEDI)?", call. = FALSE)
}

safe_fit <- function(expr) {
  tryCatch(expr, error = function(e) {
    list(error = conditionMessage(e))
  })
}

get_scenario_params <- function(Scenario) {
  if (Scenario == 1) {
    return(list(mu_y = function(x, z) 2 * x + z,
                pi_np = function(x, z) plogis(-1 + x + 0.5 * z),
                pi_p  = function(x, z) plogis(-1 + 1.5 * x - 0.5 * z)))
  } else if (Scenario == 2) {
    return(list(mu_y = function(x, z) 2 * x + z,
                pi_np = function(x, z) plogis(-1 + x + 0.5 * z),
                pi_p  = function(x, z) plogis(-1 + 1.5 * x - 0.5 * z)))
  } else if (Scenario == 3) {
    return(list(mu_y = function(x, z) 2 * x + z,
                pi_np = function(x, z) plogis(-1 + x^2 + 0.5 * z),
                pi_p  = function(x, z) plogis(-1 + 1.5 * x - 0.5 * z)))
  } else if (Scenario == 4) {
    return(list(mu_y = function(x, z) 2 * x + z,
                pi_np = function(x, z) plogis(-1 + x^2 + 0.5 * z),
                pi_p  = function(x, z) plogis(-1 + 1.5 * x^2 - 0.5 * z)))
  } else {
    stop("Unknown Scenario: ", Scenario)
  }
}

# Data generator --------------------------------------------------------------
generate_dualframe_population <- function(N, Scenario, pi_p_offset = 0) {

  params <- get_scenario_params(Scenario)

  x <- runif(N, -1, 1)
  z <- rbinom(N, 1, 0.5)
  eps <- rnorm(N, 0, 1)

  mu_y <- params$mu_y(x, z)
  y_true <- mu_y + eps
  y_obs  <- y_true

  # Sampling probabilities
  pi_np <- params$pi_np(x, z)
  d_np <- rbinom(N, 1, pi_np)

  pi_p_true <- params$pi_p(x, z)
  pi_p_true <- pmin(pmax(pi_p_offset + pi_p_true, 0), 1)
  d_p <- rbinom(N, 1, pi_p_true)

  # --- Observation pattern (fixed) ---
  # 1) y is missing for (0,0)
  y_obs[d_np == 0 & d_p == 0] <- NA

  # 2) pi_p observed only when d_p==1
  pi_p_obs <- pi_p_true
  pi_p_obs[d_p == 0] <- NA

  X <- cbind(x = x, z = z)

  data.frame(
    X = I(X),
    y = y_obs,
    y_true = y_true,
    d_np = d_np,
    d_p  = d_p,
    pi_p = pi_p_obs,
    pi_p_true = pi_p_true,
    pi_np = pi_np
  )
}

# One-rep fit --------------------------------------------------------------
fit_all_estimators_once <- function(dat, K = 2, Eff_type = 2, x_info = TRUE) {

  # Full frame size (population size)
  N_pop <- nrow(dat)

  # Union-only data (drop (0,0)) for estimators that assume no X for (0,0)
  d_np <- as.numeric(dat$d_np)
  d_p  <- as.numeric(dat$d_p)

  dat_union <- dat[d_np == 1 | d_p == 1, , drop = FALSE]

  # base_fun required by NP / NP_P / Chang-Kott estimators (uses X only)
  base_fun <- function(Xmat) {
    Xmat <- as.matrix(Xmat)
    cbind(1, Xmat, Xmat[, 1]^2)
  }

  # 1) P-only mean estimator
  est_P <- safe_fit(df_estimate_P(dat))

  # 2) NP-only mean estimator (Chang-Kott with p.use=FALSE)
  est_NP <- safe_fit(df_estimate_NP(dat, base_fun = base_fun))

  # 3) NP+P mean estimator (Chang-Kott with p.use=TRUE)
  est_NP_P <- safe_fit(df_estimate_NP_P(dat, base_fun = base_fun))

  # Efficient estimator Eff on full data
  est_Eff <- safe_fit(Eff(dat, K = K, type = Eff_type, x_info = x_info))

  # Efficient estimator Eff on union-only data (x_info=FALSE) with scaling by N
  est_Eff_union <- safe_fit(Eff(dat_union, K = K, type = Eff_type, x_info = FALSE, N = N_pop))

  # Sub-efficient estimator Eff_S
  est_Eff_S <- safe_fit(Eff_S(dat, K = K, x_info = x_info))

  # Parametric efficient estimator Eff_P
  est_Eff_P <- safe_fit(Eff_P(dat, x_info = x_info))

  # Collect results into a single-row data.frame
  extract_theta <- function(obj) if (is.list(obj) && !("error" %in% names(obj))) obj$theta else NA_real_
  extract_se    <- function(obj) if (is.list(obj) && !("error" %in% names(obj))) obj$se    else NA_real_

  out <- data.frame(
    theta_P       = extract_theta(est_P),
    se_P          = extract_se(est_P),

    theta_NP      = extract_theta(est_NP),
    se_NP         = extract_se(est_NP),

    theta_NP_P    = extract_theta(est_NP_P),
    se_NP_P       = extract_se(est_NP_P),

    theta_Eff     = extract_theta(est_Eff),
    se_Eff        = extract_se(est_Eff),

    theta_EffU    = extract_theta(est_Eff_union),
    se_EffU       = extract_se(est_Eff_union),

    theta_Eff_S   = extract_theta(est_Eff_S),
    se_Eff_S      = extract_se(est_Eff_S),

    theta_Eff_P   = extract_theta(est_Eff_P),
    se_Eff_P      = extract_se(est_Eff_P)
  )

  out
}

# One Monte Carlo replicate (long format)
mc_one_rep_long <- function(seed, N, Scenario, K = 2, Eff_type = 2, x_info = TRUE, pi_p_offset = 0,
                            show_progress = FALSE) {

  set.seed(seed)

  dat <- generate_dualframe_population(N = N, Scenario = Scenario, pi_p_offset = pi_p_offset)

  if (isTRUE(show_progress)) {
    cat(sprintf("Rep seed=%d: fitting estimators...\n", seed))
    flush.console()
  }

  res <- fit_all_estimators_once(dat, K = K, Eff_type = Eff_type, x_info = x_info)

  res$seed <- seed
  res
}

# Run Monte Carlo ----------------------------------------------------------
run_mc <- function(B = 100,
                   N = 10000,
                   Scenario = 1,
                   K = 2,
                   Eff_type = 2,
                   x_info = TRUE,
                   seed_start = 1,
                   show_progress = FALSE,
                   progress_each_fit = FALSE,
                   pi_p_offset = 0,
                   parallel = c("none", "multicore", "snow"),
                   n_cores = NULL) {

  parallel <- match.arg(parallel)

  seeds <- seed_start + seq_len(B) - 1L

  if (parallel == "none") {
    res_list <- vector("list", B)
    for (b in seq_len(B)) {
      res_list[[b]] <- mc_one_rep_long(
        seed = seeds[b],
        N = N,
        Scenario = Scenario,
        K = K,
        Eff_type = Eff_type,
        x_info = x_info,
        pi_p_offset = pi_p_offset,
        show_progress = progress_each_fit
      )
      if (show_progress) {
        cat(sprintf("Completed %d/%d\n", b, B))
        flush.console()
      }
    }
  } else if (parallel == "multicore") {
    if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
    res_list <- parallel::mclapply(
      seeds,
      function(s) mc_one_rep_long(
        seed = s,
        N = N,
        Scenario = Scenario,
        K = K,
        Eff_type = Eff_type,
        x_info = x_info,
        pi_p_offset = pi_p_offset,
        show_progress = progress_each_fit
      ),
      mc.cores = n_cores
    )
  } else if (parallel == "snow") {
    if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c(
        "safe_fit",
        "get_scenario_params",
        "generate_dualframe_population",
        "fit_all_estimators_once",
        "mc_one_rep_long",
        "df_estimate_P",
        "df_estimate_NP",
        "df_estimate_NP_P",
        "Eff",
        "Eff_S",
        "Eff_P"
      ),
      envir = environment()
    )

    res_list <- parallel::parLapply(
      cl,
      seeds,
      function(s) mc_one_rep_long(
        seed = s,
        N = N,
        Scenario = Scenario,
        K = K,
        Eff_type = Eff_type,
        x_info = x_info,
        pi_p_offset = pi_p_offset,
        show_progress = progress_each_fit
      )
    )
  }

  do.call(rbind, res_list)
}

# Summaries ---------------------------------------------------------------
summarize_mc <- function(res, theta_true = NULL) {
  if (!is.data.frame(res)) stop("summarize_mc(): res must be a data.frame.")

  cols <- setdiff(names(res), "seed")

  out <- data.frame(
    estimator = cols,
    mean = NA_real_,
    sd = NA_real_,
    rmse = NA_real_,
    bias = NA_real_
  )

  for (j in seq_along(cols)) {
    nm <- cols[j]
    v <- res[[nm]]
    v <- as.numeric(v)
    v <- v[is.finite(v)]

    out$mean[j] <- mean(v)
    out$sd[j]   <- stats::sd(v)

    if (!is.null(theta_true) && is.finite(theta_true)) {
      out$bias[j] <- out$mean[j] - theta_true
      out$rmse[j] <- sqrt(mean((v - theta_true)^2))
    }
  }

  out
}
