############################################################
## inst/examples/dualframe_simulation.R
############################################################

library(dfSEDI)

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
# Data generating process (paper-style scenarios)
#
# S1: Outcome O1 and NP mechanism NP1, available L=(X,Y)
#   X ~ N(0,1)
#   Y = -exp(-2) + cos(2X) + 0.5X + eps
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
    X     = I(X),   # store as a matrix-column
    y     = y,
    d_np  = d_np,
    d_p   = d_p,
    pi_p  = pi_p,
    pi_np = pi_np
  )
}

#------------------------------------------------------------
# Default basis function for Chang & Kott-type NP estimators
# - S1/S2: g(x) = (1, x, x^2)^T  -> returns n x (p+2) with p=1 => 3 columns
# - S3:    g(x) = (1, x, x^2, z)^T -> p=2 => 4 columns
#------------------------------------------------------------
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

#------------------------------------------------------------
# Fit all requested estimators once
#   P, NP, NP_P, Eff, Eff(union_dat), Eff_S, Eff_P
#------------------------------------------------------------
fit_all_estimators_once <- function(dat, Scenario, K = 2, progress_each = FALSE) {
  base_fun <- make_base_fun(Scenario)

  n_np    <- sum(dat$d_np == 1)
  n_p     <- sum(dat$d_p  == 1)
  n_union <- sum(dat$d_np == 1 | dat$d_p == 1)

  # (1) Probability-only (HT-type)
  fit_P <- safe_fit(df_estimate_P(dat))

  # (2) Non-probability-only and (3) NP union P (Chang & Kott-type)
  fit_NP   <- safe_fit(df_estimate_NP(dat, base_fun = base_fun))
  fit_NP_P <- safe_fit(df_estimate_NP_P(dat, base_fun = base_fun))

  # (4) Efficient (full data with (0,0) rows allowed): x_info = TRUE
  fit_Eff <- safe_fit(Eff(
    dat         = dat,
    K           = K,
    phi_start   = NULL,
    max_restart = 10,
    progress    = progress_each,
    x_info      = TRUE
  ))

  # (5) Efficient using union sample only: x_info = FALSE
  dat_union <- dat[dat$d_np == 1 | dat$d_p == 1, , drop = FALSE]
  fit_Eff_union <- safe_fit(Eff(
    dat         = dat_union,
    K           = K,
    phi_start   = NULL,
    max_restart = 10,
    progress    = progress_each,
    x_info      = FALSE
  ))

  # (6) Sub-efficient
  fit_EffS <- safe_fit(Eff_S(
    dat      = dat,
    K        = K,
    progress = progress_each,
    x_info   = TRUE
  ))

  # (7) Parametric efficient (working model)
  fit_EffP <- safe_fit(Eff_P(
    dat       = dat,
    phi_start = NULL,
    eta4_star = 0,
    max_iter  = 20,
    progress  = progress_each,
    x_info    = TRUE
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
## One Monte Carlo replication (returns fits + theta vector)
############################################################
simulate_one_replication <- function(seed,
                                     N        = 10000,
                                     K        = 2,
                                     Scenario = 1,
                                     progress = FALSE) {

  set.seed(seed)
  dat <- generate_dualframe_population(N = N, Scenario = Scenario)

  fits <- fit_all_estimators_once(dat = dat, Scenario = Scenario, K = K, progress_each = progress)

  theta_vec <- c(
    P             = fits$P$theta,
    NP            = fits$NP$theta,
    NP_P          = fits$NP_P$theta,
    Eff           = fits$Eff$theta,
    Eff_union_dat = fits$Eff_union$theta,
    Eff_S         = fits$Eff_S$theta,
    Eff_P         = fits$Eff_P$theta
  )

  list(
    theta = theta_vec,
    fits  = fits,
    dat   = dat
  )
}

############################################################
## Monte Carlo runner (long-format output)
############################################################
run_mc <- function(B,
                   N = 10000,
                   Scenario = 1,
                   K = 2,
                   seed_start = 1,
                   show_progress = TRUE,
                   progress_each_fit = FALSE,
                   pi_p_offset = 0.005) {

  out_list <- vector("list", B)

  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (b in seq_len(B)) {
    seed <- seed_start + b - 1
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
      progress_each = progress_each_fit
    )

    out_list[[b]] <- rbind(
      extract_row(fits$P,         "P",             b, Scenario, fits$n_np, fits$n_p, fits$n_union),
      extract_row(fits$NP,        "NP",            b, Scenario, fits$n_np, fits$n_p, fits$n_union),
      extract_row(fits$NP_P,      "NP_P",          b, Scenario, fits$n_np, fits$n_p, fits$n_union),
      extract_row(fits$Eff,       "Eff",           b, Scenario, fits$n_np, fits$n_p, fits$n_union),
      extract_row(fits$Eff_union, "Eff_union_dat", b, Scenario, fits$n_np, fits$n_p, fits$n_union),
      extract_row(fits$Eff_S,     "Eff_S",         b, Scenario, fits$n_np, fits$n_p, fits$n_union),
      extract_row(fits$Eff_P,     "Eff_P",         b, Scenario, fits$n_np, fits$n_p, fits$n_union)
    )

    if (show_progress) setTxtProgressBar(pb, b)
  }

  do.call(rbind, out_list)
}

#------------------------------------------------------------
# Optional: summarize MC results (mean, bias, RMSE, coverage)
# theta_true defaults to 0 (true mean is 0 in the bundled scenarios)
#------------------------------------------------------------
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
      sd         = sd(dd$theta),
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
## main(): single dataset test, no Monte Carlo
############################################################
main <- function(N = 10000, K = 2, Scenario = 1, x_info = TRUE) {

  dat <- generate_dualframe_population(N = N, Scenario = Scenario)

  sc <- normalize_scenario(Scenario)

  # Optional starting value for phi (useful for debugging only)
  # For S1/S2, the true NP mechanism is linear in (1, x, y).
  # For S3, the mechanism is nonlinear in y, so there is no "true phi" under the working model.
  phi_start <- if (sc %in% c(1, 2)) c(-2.15, -0.5, -0.75) else NULL

  # 1) Efficient estimator Eff
  fit_eff <- Eff(
    dat         = dat,
    K           = K,
    phi_start   = phi_start,
    max_restart = 10,
    progress    = TRUE,
    x_info      = x_info
  )

  cat("\n=== Eff ===\n")
  print(fit_eff$phi)
  print(fit_eff$theta)
  print(fit_eff$se)
  print(fit_eff$ci)

  # 2) Sub-efficient estimator Eff_S
  fit_effS <- Eff_S(
    dat      = dat,
    K        = K,
    progress = TRUE,
    x_info   = x_info
  )

  cat("\n=== Eff_S ===\n")
  print(fit_effS$theta)
  print(fit_effS$se)
  print(fit_effS$ci)

  # 3) Parametric efficient estimator Eff_P
  fit_effP <- Eff_P(
    dat       = dat,
    phi_start = phi_start,
    eta4_star = 0,
    max_iter  = 20,
    progress  = TRUE,
    x_info    = x_info
  )

  cat("\n=== Eff_P ===\n")
  print(fit_effP$phi)
  print(fit_effP$theta)
  print(fit_effP$se)
  print(fit_effP$ci)

  invisible(list(
    dat     = dat,
    Eff     = fit_eff,
    Eff_S   = fit_effS,
    Eff_P   = fit_effP
  ))
}
