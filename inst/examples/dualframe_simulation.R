############################################################
## inst/examples/dualframe_simulation.R
##
## Dual-frame simulation + MC runner for dfSEDI
##
## - Scenarios S1/S2/S3 (paper-style DGP)
## - Fit: P / NP / NP_P / Eff / Eff_union_dat / Eff_S / Eff_P
## - MC: run_mc(B=..., ...) returns long-format results
## - Summary: summarize_mc(res, theta_true=...)
## - Parallel: run_mc(..., parallel="multicore" or "snow", n_cores=...)
##
## Notes:
## - This script defines functions only; it does NOT run MC by default.
## - For Windows, use parallel="snow" (multicore is not supported).
##
## IMPORTANT (Eff default in dfSEDI):
## - Recent versions of dfSEDI use DML2 as the default for Eff().
## - To make simulations reproducible across versions, this script
##   *explicitly passes* `type = Eff_type` to Eff().
## - The default Eff_type in this script is set to 2 (DML2).
############################################################

# dfSEDI is the main dependency
library(dfSEDI)

############################################################
## 0) Utilities
############################################################

#------------------------------------------------------------
# Normalize scenario input: 1/2/3 or "S1"/"S2"/"S3"
#------------------------------------------------------------
normalize_scenario <- function(Scenario) {
  if (is.character(Scenario)) {
    s <- toupper(Scenario)
    s <- sub("^S", "", s)
    Scenario <- suppressWarnings(as.integer(s))
  }
  Scenario <- as.integer(Scenario)
  if (!(Scenario %in% 1:3)) {
    stop("Scenario must be 1, 2, 3 or 'S1', 'S2', 'S3'.")
  }
  Scenario
}

#------------------------------------------------------------
# Safe wrapper: return NA outputs instead of stopping the whole loop
#------------------------------------------------------------
safe_fit <- function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      list(theta = NA_real_, se = NA_real_, ci = c(NA_real_, NA_real_), error = e$message)
    }
  )
}

#------------------------------------------------------------
# Convert a single fit object into one row (long format)
#------------------------------------------------------------
extract_row <- function(fit, estimator, rep, Scenario, n_np, n_p, n_union) {
  ci <- fit$ci
  if (is.null(ci) || length(ci) != 2) ci <- c(NA_real_, NA_real_)

  data.frame(
    Scenario  = paste0("S", normalize_scenario(Scenario)),
    rep       = rep,
    estimator = estimator,
    theta     = as.numeric(fit$theta),
    se        = as.numeric(fit$se),
    ci_l      = as.numeric(ci[1]),
    ci_u      = as.numeric(ci[2]),
    n_np      = n_np,
    n_p       = n_p,
    n_union   = n_union,
    error     = if (!is.null(fit$error)) fit$error else NA_character_,
    stringsAsFactors = FALSE
  )
}

############################################################
## 1) Data generating process (paper-style scenarios)
############################################################
#------------------------------------------------------------
# DGP summary
#
# S1: Outcome O1 and NP mechanism NP1, available L=(X,Y)
#   X ~ N(0,1)
#   Y = -exp(-2) + cos(2X) + 0.5X + eps, eps~N(0,0.5^2)
#   logit P(d_np=1 | x,y) = -2.15 - 0.5 x - 0.75 y
#
# S2: Outcome O2 and NP mechanism NP1, available L=(X,Y)
#   Y = 0.8X + eps
#
# S3: Outcome O3 and NP mechanism NP2, available L=(X,Z,Y)
#   X ~ N(0,1), Z ~ Bern(0.5)
#   Y = 0.2 + 0.8X - 0.4Z + eps
#   logit P(d_np=1 | x,z,y) = -1.5 - 0.3 cos(2y) - 0.1 (y-1)^2
#
# Probability sample (all scenarios):
#   logit P(d_p=1 | x) = -3 - 0.25 (x-2)^2
#   pi_p = pi_p_offset + logistic(lin_p)
#------------------------------------------------------------
generate_dualframe_population <- function(N, Scenario = 1, pi_p_offset = 0.005) {
  sc <- normalize_scenario(Scenario)

  # Covariates
  x <- rnorm(N, 0, 1)          # X ~ N(0,1)
  z <- rbinom(N, 1, 0.5)       # Z ~ Bern(0.5) (generated for all; used in S3)

  # Outcome error: eps ~ N(0, (1/2)^2) => sd = 0.5
  eps <- rnorm(N, mean = 0, sd = 0.5)

  # Outcome model
  if (sc == 1) {
    mu_y <- -exp(-2) + cos(2 * x) + 0.5 * x
  } else if (sc == 2) {
    mu_y <- 0.8 * x
  } else {
    mu_y <- 0.2 + 0.8 * x - 0.4 * z
  }
  y <- mu_y + eps

  # Non-probability sampling mechanism
  if (sc %in% c(1, 2)) {
    lin_np <- -2.15 - 0.5 * x - 0.75 * y
  } else {
    lin_np <- -1.5 - 0.3 * cos(2 * y) - 0.1 * (y - 1)^2
  }
  pi_np <- plogis(lin_np)
  d_np  <- rbinom(N, 1, pi_np)

  # Probability sampling design
  lin_p <- -3.0 - 0.25 * (x - 2)^2
  pi_p  <- pi_p_offset + plogis(lin_p)
  pi_p  <- pmin(pmax(pi_p, 0), 1)
  d_p   <- rbinom(N, 1, pi_p)

  # Match "available data" dimension by scenario
  # S1/S2: X is n x 1 (x only)
  # S3:    X is n x 2 (x and z)
  X <- if (sc %in% c(1, 2)) {
    matrix(x, ncol = 1)
  } else {
    cbind(x, z)
  }

  data.frame(
    X     = I(X),   # store as a matrix-column (AsIs class is OK; dfSEDI strips internally)
    y     = y,
    d_np  = d_np,
    d_p   = d_p,
    pi_p  = pi_p,
    pi_np = pi_np
  )
}

############################################################
## 2) Default basis function for Chang & Kott-type NP estimators
############################################################
# - S1/S2: g(x) = (1, x, x^2)^T  -> returns n x (p+2) with p=1 => 3 columns
# - S3:    g(x) = (1, x, x^2, z)^T -> p=2 => 4 columns
make_base_fun <- function(Scenario) {
  sc <- normalize_scenario(Scenario)

  if (sc %in% c(1, 2)) {
    base_fun <- function(X) {
      X <- as.matrix(X)
      cbind(1, X[, 1], X[, 1]^2)
    }
  } else {
    base_fun <- function(X) {
      X <- as.matrix(X)
      cbind(1, X[, 1], X[, 1]^2, X[, 2])
    }
  }
  base_fun
}

############################################################
## 3) Fit all requested estimators once (single dataset)
############################################################
#   P, NP, NP_P, Eff, Eff(union_dat), Eff_S, Eff_P
#
# Arguments:
# - K         : #folds used inside Eff/Eff_S
# - Eff_type  : 1 (DML1) or 2 (DML2). Default = 2 (matches current Eff default).
# - x_info    : TRUE (full formula), FALSE (union-only formulation)
#
fit_all_estimators_once <- function(dat,
                                    Scenario,
                                    K = 2,
                                    Eff_type = 2,
                                    x_info = TRUE,
                                    progress_each = FALSE) {

  base_fun <- make_base_fun(Scenario)

  n_np    <- sum(dat$d_np == 1)
  n_p     <- sum(dat$d_p  == 1)
  n_union <- sum(dat$d_np == 1 | dat$d_p == 1)

  # (1) Probability-only (HT-type)
  fit_P <- safe_fit(df_estimate_P(dat))

  # (2) Non-probability-only and (3) NP union P (Chang & Kott-type)
  fit_NP   <- safe_fit(df_estimate_NP(dat, base_fun = base_fun))
  fit_NP_P <- safe_fit(df_estimate_NP_P(dat, base_fun = base_fun))

  # (4) Efficient estimator (full data allowed): x_info is user-controlled here
  fit_Eff <- safe_fit(Eff(
    dat         = dat,
    K           = K,
    type        = Eff_type,        # 1=DML1, 2=DML2
    phi_start   = NULL,
    max_restart = 10,
    progress    = progress_each,
    x_info      = x_info
  ))

  # (5) Efficient estimator using union sample only:
  #     if x_info=FALSE and you want population mean E(Y), pass N.
  dat_union <- dat[dat$d_np == 1 | dat$d_p == 1, , drop = FALSE]

  eff_union_args <- list(
    dat         = dat_union,
    K           = K,
    type        = Eff_type,
    phi_start   = NULL,
    max_restart = 10,
    progress    = progress_each,
    x_info      = FALSE
  )
  if ("N" %in% names(formals(Eff))) {
    eff_union_args$N <- nrow(dat)  # true frame size in simulations
  }
  fit_Eff_union <- safe_fit(do.call(Eff, eff_union_args))

  # (6) Sub-efficient
  fit_EffS <- safe_fit(Eff_S(
    dat      = dat,
    K        = K,
    progress = progress_each,
    x_info   = x_info
  ))

  # (7) Parametric efficient (working model)
  fit_EffP <- safe_fit(Eff_P(
    dat       = dat,
    phi_start = NULL,
    eta4_star = 0,
    max_iter  = 20,
    progress  = progress_each,
    x_info    = x_info
  ))

  list(
    P         = fit_P,
    NP        = fit_NP,
    NP_P      = fit_NP_P,
    Eff       = fit_Eff,
    Eff_union = fit_Eff_union,
    Eff_S     = fit_EffS,
    Eff_P     = fit_EffP,
    n_np      = n_np,
    n_p       = n_p,
    n_union   = n_union
  )
}

############################################################
## 4) One MC replication: returns long-format rows
############################################################
mc_one_rep_long <- function(seed,
                            rep_id,
                            N,
                            Scenario,
                            K,
                            Eff_type,
                            x_info,
                            pi_p_offset,
                            progress_each_fit = FALSE) {

  set.seed(seed)

  dat <- generate_dualframe_population(
    N = N,
    Scenario = Scenario,
    pi_p_offset = pi_p_offset
  )

  fits <- fit_all_estimators_once(
    dat = dat,
    Scenario = Scenario,
    K = K,
    Eff_type = Eff_type,
    x_info = x_info,
    progress_each = progress_each_fit
  )

  rbind(
    extract_row(fits$P,         "P",             rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP,        "NP",            rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP_P,      "NP_P",          rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff,       "Eff",           rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_union, "Eff_union_dat", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_S,     "Eff_S",         rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_P,     "Eff_P",         rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union)
  )
}

############################################################
## 5) MC runner (serial / parallel)
############################################################
run_mc <- function(B,
                   N = 10000,
                   Scenario = 1,
                   K = 2,
                   Eff_type = 2,              # 1=DML1, 2=DML2 (default)
                   x_info = TRUE,
                   seed_start = 1,
                   show_progress = TRUE,
                   progress_each_fit = FALSE,
                   pi_p_offset = 0.005,
                   parallel = c("none", "multicore", "snow"),
                   n_cores = NULL) {

  parallel <- match.arg(parallel)

  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) stop("B must be a positive integer.")

  sc <- normalize_scenario(Scenario)

  seeds <- seed_start + seq_len(B) - 1L
  rep_ids <- seq_len(B)

  # worker function (must be self-contained)
  worker <- function(i) {
    mc_one_rep_long(
      seed = seeds[i],
      rep_id = rep_ids[i],
      N = N,
      Scenario = sc,
      K = K,
      Eff_type = Eff_type,
      x_info = x_info,
      pi_p_offset = pi_p_offset,
      progress_each_fit = progress_each_fit
    )
  }

  # ---- Serial ----
  if (parallel == "none") {
    out_list <- vector("list", B)

    if (show_progress) {
      pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      on.exit(close(pb), add = TRUE)
    }

    for (b in seq_len(B)) {
      out_list[[b]] <- worker(b)
      if (show_progress) utils::setTxtProgressBar(pb, b)
    }

    return(do.call(rbind, out_list))
  }

  # ---- Parallel (base 'parallel' package; recommended in R) ----
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for parallel != 'none'.", call. = FALSE)
  }

  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores(logical = TRUE)
    if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
    # leave 1 core free by default
    n_cores <- max(1L, as.integer(n_cores) - 1L)
  }
  n_cores <- as.integer(n_cores)
  if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L

  if (parallel == "multicore") {
    # NOTE: not supported on Windows
    out_list <- parallel::mclapply(seq_len(B), worker, mc.cores = n_cores)
    return(do.call(rbind, out_list))
  }

  # snow cluster (works on Windows/macOS/Linux)
  cl <- parallel::makeCluster(n_cores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  # Make sure workers can see dfSEDI + functions defined in this script
  parallel::clusterEvalQ(cl, library(dfSEDI))

  # Export everything the worker uses (functions + objects)
  parallel::clusterExport(
    cl,
    varlist = c(
      "normalize_scenario",
      "generate_dualframe_population",
      "make_base_fun",
      "safe_fit",
      "extract_row",
      "fit_all_estimators_once",
      "mc_one_rep_long",
      "seeds",
      "rep_ids",
      "N",
      "sc",
      "K",
      "Eff_type",
      "x_info",
      "pi_p_offset",
      "progress_each_fit"
    ),
    envir = environment()
  )

  out_list <- parallel::parLapply(cl, seq_len(B), function(i) {
    mc_one_rep_long(
      seed = seeds[i],
      rep_id = rep_ids[i],
      N = N,
      Scenario = sc,
      K = K,
      Eff_type = Eff_type,
      x_info = x_info,
      pi_p_offset = pi_p_offset,
      progress_each_fit = progress_each_fit
    )
  })

  do.call(rbind, out_list)
}

############################################################
## 6) MC summary
############################################################
# Summarize MC results (mean, bias, RMSE, coverage)
# theta_true defaults to 0 (true mean is 0 in the bundled scenarios)
summarize_mc <- function(res, theta_true = 0) {
  split_list <- split(res, res$estimator)

  sum_df <- do.call(rbind, lapply(split_list, function(d) {
    ok <- is.finite(d$theta) & is.finite(d$se) & is.finite(d$ci_l) & is.finite(d$ci_u)
    dd <- d[ok, , drop = FALSE]

    if (nrow(dd) == 0) {
      return(data.frame(
        estimator  = d$estimator[1],
        n_total    = nrow(d),
        n_ok       = 0,
        mean_theta = NA_real_,
        bias       = NA_real_,
        sd         = NA_real_,
        rmse       = NA_real_,
        mean_se    = NA_real_,
        coverage   = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      estimator  = dd$estimator[1],
      n_total    = nrow(d),
      n_ok       = nrow(dd),
      mean_theta = mean(dd$theta),
      bias       = mean(dd$theta) - theta_true,
      sd         = stats::sd(dd$theta),
      rmse       = sqrt(mean((dd$theta - theta_true)^2)),
      mean_se    = mean(dd$se),
      coverage   = mean(dd$ci_l <= theta_true & theta_true <= dd$ci_u),
      stringsAsFactors = FALSE
    )
  }))

  rownames(sum_df) <- NULL
  sum_df[order(sum_df$estimator), ]
}

############################################################
## 7) Single dataset quick test (no MC)
############################################################
main <- function(N = 10000,
                 K = 2,
                 Scenario = 1,
                 Eff_type = 2,
                 x_info = TRUE,
                 pi_p_offset = 0.005) {

  dat <- generate_dualframe_population(
    N = N,
    Scenario = Scenario,
    pi_p_offset = pi_p_offset
  )

  fits <- fit_all_estimators_once(
    dat = dat,
    Scenario = Scenario,
    K = K,
    Eff_type = Eff_type,
    x_info = x_info,
    progress_each = TRUE
  )

  cat("\n=== Eff ===\n")
  print(fits$Eff$phi)
  print(fits$Eff$theta)
  print(fits$Eff$se)
  print(fits$Eff$ci)

  invisible(list(dat = dat, fits = fits))
}

############################################################
## 8) Usage examples (commented out)
############################################################

# ---- Example: single run ----
# source(system.file("examples", "dualframe_simulation.R", package="dfSEDI"))
# out <- main(N=10000, K=2, Scenario=1, Eff_type=2, x_info=TRUE)
#
# ---- Example: MC (serial; default Eff_type=2 i.e. DML2) ----
# res <- run_mc(B=100, N=10000, Scenario=1, K=2, Eff_type=2, x_info=TRUE, parallel="none")
# summarize_mc(res, theta_true=0)
#
# ---- Example: MC (parallel, macOS/Linux) ----
# res <- run_mc(B=200, N=10000, Scenario=1, K=2, Eff_type=2, x_info=TRUE,
#               parallel="multicore", n_cores=4)
# summarize_mc(res)
#
# ---- Example: MC (parallel, Windows/macOS/Linux) ----
# res <- run_mc(B=200, N=10000, Scenario=1, K=2, Eff_type=2, x_info=TRUE,
#               parallel="snow", n_cores=4)
# summarize_mc(res)
