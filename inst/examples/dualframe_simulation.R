############################################################
## dualframe_simulation.R
##
## Dual-frame simulation + MC runner for dfSEDI
##
## - Scenarios S1/S2/S3 (paper-style DGP)
## - Fit: P / NP / NP_P / Eff1 / Eff2 / Eff1_union / Eff2_union / Eff_S / Eff_P
## - MC: run_mc(B=..., ...) returns long-format results
## - Summary: summarize_mc(res, theta_true=...)
## - Parallel:
##     * macOS/Linux: parallel="multicore" (fast; no progress bar in base R)
##     * Windows/macOS/Linux: parallel="snow" (PSOCK; progress bar + load balancing)
##
## Notes
## - Eff_type defaults to 2 (DML2) to match the current Eff() default in dfSEDI.
## - The snow implementation uses internal parallel:::sendCall/recvOneResult to enable
##   dynamic scheduling + a progress bar. This avoids the common "idle CPU" issue.
##
## NEW (this revision):
## - MC output now also includes phi estimates and (when available) their SE/CI.
##   Columns added (fixed width up to 4 coefficients):
##     phi_1..phi_4, phi_se_1..phi_se_4, phi_ci_l_1..phi_ci_l_4, phi_ci_u_1..phi_ci_u_4
##   Notes:
##   - Eff returns phi + phi_se + phi_ci (lower/upper).
##   - NP / NP_P return phi but (currently) do not return phi_se/phi_ci -> filled with NA.
##   - Estimators without phi -> all phi fields are NA.
############################################################

library(dfSEDI)

############################################################
## 0) Utilities
############################################################

normalize_scenario <- function(Scenario) {
  if (is.character(Scenario)) {
    s <- toupper(Scenario)
    s <- sub("^S", "", s)
    Scenario <- suppressWarnings(as.integer(s))
  }
  Scenario <- as.integer(Scenario)
  if (!(Scenario %in% 1:3)) stop("Scenario must be 1, 2, 3 or 'S1', 'S2', 'S3'.")
  Scenario
}

safe_fit <- function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      list(
        theta = NA_real_,
        se    = NA_real_,
        ci    = c(NA_real_, NA_real_),
        phi   = NA_real_,
        phi_se = NA_real_,
        phi_ci = matrix(NA_real_, nrow = 2, ncol = 1,
                        dimnames = list(c("lower","upper"), NULL)),
        error = e$message
      )
    }
  )
}

extract_row <- function(fit, estimator, rep_id, Scenario, n_np, n_p, n_union) {
  ci <- fit$ci
  if (is.null(ci) || length(ci) != 2) ci <- c(NA_real_, NA_real_)

  scenario_lab <- tryCatch(
    paste0("S", normalize_scenario(Scenario)),
    error = function(e) NA_character_
  )

  theta_val <- if (!is.null(fit$theta) && length(fit$theta) == 1) as.numeric(fit$theta) else NA_real_
  se_val    <- if (!is.null(fit$se)   && length(fit$se)   == 1) as.numeric(fit$se)   else NA_real_

  out <- data.frame(
    Scenario  = scenario_lab,
    rep       = as.integer(rep_id),
    estimator = as.character(estimator),
    theta     = theta_val,
    se        = se_val,
    ci_l      = as.numeric(ci[1]),
    ci_u      = as.numeric(ci[2]),
    n_np      = as.integer(n_np),
    n_p       = as.integer(n_p),
    n_union   = as.integer(n_union),
    error     = if (!is.null(fit$error) && length(fit$error) > 0) as.character(fit$error)[1] else NA_character_,
    stringsAsFactors = FALSE
  )

  p_phi_max <- 4L

  phi    <- if (!is.null(fit$phi))    as.numeric(fit$phi)    else numeric(0)
  phi_se <- if (!is.null(fit$phi_se)) as.numeric(fit$phi_se) else numeric(0)

  phi    <- c(phi,    rep(NA_real_, p_phi_max))[seq_len(p_phi_max)]
  phi_se <- c(phi_se, rep(NA_real_, p_phi_max))[seq_len(p_phi_max)]

  phi_ci_l <- rep(NA_real_, p_phi_max)
  phi_ci_u <- rep(NA_real_, p_phi_max)
  if (!is.null(fit$phi_ci) && is.matrix(fit$phi_ci) && nrow(fit$phi_ci) == 2) {
    phi_ci_l <- c(as.numeric(fit$phi_ci[1, ]), rep(NA_real_, p_phi_max))[seq_len(p_phi_max)]
    phi_ci_u <- c(as.numeric(fit$phi_ci[2, ]), rep(NA_real_, p_phi_max))[seq_len(p_phi_max)]
  }

  for (j in seq_len(p_phi_max)) {
    out[[paste0("phi_",      j)]] <- phi[j]
    out[[paste0("phi_se_",   j)]] <- phi_se[j]
    out[[paste0("phi_ci_l_", j)]] <- phi_ci_l[j]
    out[[paste0("phi_ci_u_", j)]] <- phi_ci_u[j]
  }

  out
}


############################################################
## 1) Data generating process (paper-style scenarios)
############################################################

generate_dualframe_population <- function(N, Scenario = 1, pi_p_offset = 0.005) {
  sc <- normalize_scenario(Scenario)

  x <- rnorm(N, 0, 1)
  z <- rbinom(N, 1, 0.5)
  eps <- rnorm(N, mean = 0, sd = 0.5)

  if (sc == 1) {
    mu_y <- -exp(-2) + cos(2 * x) + 0.5 * x
  } else if (sc == 2) {
    mu_y <- 0.8 * x
  } else {
    mu_y <- 0.2 + 0.8 * x - 0.4 * z
  }
  y <- mu_y + eps

  if (sc %in% c(1, 2)) {
    lin_np <- -2.15 - 0.5 * x - 0.75 * y
  } else {
    lin_np <- -1.5 - 0.3 * cos(2 * y) - 0.1 * (y - 1)^2
  }
  pi_np <- plogis(lin_np)
  d_np  <- rbinom(N, 1, pi_np)

  lin_p <- -3.0 - 0.25 * (x - 2)^2
  pi_p  <- pi_p_offset + plogis(lin_p)
  pi_p  <- pmin(pmax(pi_p, 0), 1)
  d_p   <- rbinom(N, 1, pi_p)

  # Requested: pi_p is observed only for d_p==1 units.
  pi_p[d_p == 0] <- NA_real_

  X <- if (sc %in% c(1, 2)) {
    matrix(x, ncol = 1)
  } else {
    cbind(x, z)
  }

  data.frame(
    X     = I(X),   # matrix-column; AsIs is OK (dfSEDI strips internally)
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
## 3) Fit all estimators once (single dataset)
############################################################

fit_all_estimators_once <- function(dat, Scenario, K = 2, progress_each = FALSE) {
  sc <- normalize_scenario(Scenario)
  base_fun <- make_base_fun(sc)

  n_np    <- sum(dat$d_np == 1)
  n_p     <- sum(dat$d_p == 1)
  n_union <- sum(dat$d_np == 1 | dat$d_p == 1)

  # Initialization for pi_np = expit((1, X, y)' phi)
  if (sc %in% c(1, 2)) {
    phi_start_true <- c(-2.15, -0.5, -0.75)
  } else {
    # Scenario 3 uses a nonlinear pi_np DGP; use a simple finite start.
    phi_start_true <- c(-1.5, 0, 0, 0)
  }

  phi_start_NP    <- phi_start_true
  phi_start_NP_P  <- phi_start_true
  phi_start_Eff1  <- phi_start_true
  phi_start_Eff2  <- phi_start_true

  fit_P   <- safe_fit(df_estimate_P(dat))
  fit_NP  <- safe_fit(df_estimate_NP(dat, base_fun = base_fun, phi_start = phi_start_NP))
  fit_NP_P <- safe_fit(df_estimate_NP_P(dat, base_fun = base_fun, phi_start = phi_start_NP_P))

  fit_Eff1 <- safe_fit(Eff(dat = dat, K = K, x_info = TRUE, type = 1, phi_start = phi_start_Eff1, progress = progress_each))
  fit_Eff2 <- safe_fit(Eff(dat = dat, K = K, x_info = TRUE, type = 2, phi_start = phi_start_Eff2, progress = progress_each))

  dat_union <- subset(dat, d_np == 1 | d_p == 1)

  eff_union_args <- list(
    dat       = dat_union,
    K         = K,
    x_info    = FALSE,
    progress  = progress_each,
    type      = 1,
    phi_start = phi_start_Eff1
  )
  if ("N" %in% names(formals(Eff))) {
    eff_union_args$N <- nrow(dat)
  }
  fit_Eff1_union <- safe_fit(do.call(Eff, eff_union_args))

  eff_union_args$type      <- 2
  eff_union_args$phi_start <- phi_start_Eff2
  fit_Eff2_union <- safe_fit(do.call(Eff, eff_union_args))

  fit_Eff_S <- safe_fit(Eff_S(dat = dat, K = K, x_info = TRUE, progress = progress_each))
  fit_Eff_P <- safe_fit(Eff_P(dat = dat, x_info = TRUE, phi_start = phi_start_true, progress = progress_each))

  list(
    P          = fit_P,
    NP         = fit_NP,
    NP_P       = fit_NP_P,
    Eff1       = fit_Eff1,
    Eff2       = fit_Eff2,
    Eff1_union = fit_Eff1_union,
    Eff2_union = fit_Eff2_union,
    Eff_S      = fit_Eff_S,
    Eff_P      = fit_Eff_P,
    n_np       = n_np,
    n_p        = n_p,
    n_union    = n_union
  )
}


############################################################
## 4) One MC replication -> long-format rows
############################################################

mc_one_rep_long <- function(
    seed,
    rep_id,
    N,
    Scenario,
    K = 2,
    Eff_type = 2,
    x_info = TRUE,
    pi_p_offset = 0.005,
    progress_each_fit = FALSE
) {
  # Eff_type is kept for backward compatibility; we return both Eff1 and Eff2.
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
    x_info = x_info,
    progress_each = progress_each_fit
  )

  rbind(
    extract_row(fits$P,         "P",          rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP,        "NP",         rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP_P,      "NP_P",       rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_S,     "Eff_S",      rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_P,     "Eff_P",      rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff1_union,"Eff1_union", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff2_union,"Eff2_union", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff1,      "Eff1",       rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff2,      "Eff2",       rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union)
  )
}


############################################################
## 5) Snow helper: dynamic scheduling + progress bar
############################################################

df_par_lapply_lb_progress <- function(cl, X, fun, show_progress = TRUE) {
  if (!requireNamespace("parallel", quietly = TRUE)) stop("Package 'parallel' is required.")

  # Use internal functions explicitly; avoids "not exported" errors from parallel::sendCall
  sendCall      <- parallel:::sendCall
  recvOneResult <- parallel:::recvOneResult

  X <- as.list(X)
  n <- length(X)
  if (n == 0L) return(list())

  nworkers <- length(cl)
  nworkers <- min(nworkers, n)

  pb <- NULL
  if (isTRUE(show_progress)) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  next_i <- 1L
  for (w in seq_len(nworkers)) {
    sendCall(cl[[w]], fun, list(X[[next_i]]), tag = next_i)
    next_i <- next_i + 1L
  }

  out <- vector("list", n)
  done <- 0L

  while (done < n) {
    res <- recvOneResult(cl)
    idx <- res$tag
    out[[idx]] <- res$value
    done <- done + 1L
    if (!is.null(pb)) utils::setTxtProgressBar(pb, done)

    if (next_i <= n) {
      sendCall(cl[[res$node]], fun, list(X[[next_i]]), tag = next_i)
      next_i <- next_i + 1L
    }
  }

  out
}

############################################################
## 6) MC runner (serial / parallel)
############################################################
# Output:
# - returns a data.frame in long format
# - rows = 9 * B
# - columns (core): Scenario, rep, estimator, theta, se, ci_l, ci_u, n_np, n_p, n_union, error
# - plus phi fields: phi_1..phi_4, phi_se_1..phi_se_4, phi_ci_l_1..phi_ci_l_4, phi_ci_u_1..phi_ci_u_4
run_mc <- function(B,
                   N = 10000,
                   Scenario = 1,
                   K = 2,
                   Eff_type = 2,              # 2=DML2 (default), 1=DML1
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

  worker_i <- function(i) {
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

  if (parallel == "none") {
    out_list <- vector("list", B)
    pb <- NULL
    if (isTRUE(show_progress)) {
      pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      on.exit(try(close(pb), silent = TRUE), add = TRUE)
    }
    for (b in seq_len(B)) {
      out_list[[b]] <- worker_i(b)
      if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
    }
    return(do.call(rbind, out_list))
  }

  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for parallel != 'none'.", call. = FALSE)
  }

  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores(logical = TRUE)
    if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
    n_cores <- max(1L, as.integer(n_cores) - 1L)
  }
  n_cores <- as.integer(n_cores)
  if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
  n_cores <- min(n_cores, B)

  if (parallel == "multicore") {
    if (isTRUE(show_progress)) {
      message("parallel='multicore': show_progress is ignored (no progress bar in base R for mclapply).")
    }
    out_list <- parallel::mclapply(seq_len(B), worker_i, mc.cores = n_cores, mc.preschedule = FALSE)
    return(do.call(rbind, out_list))
  }

  cl <- parallel::makeCluster(n_cores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

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
    envir = environment(run_mc)
  )

  idxs <- seq_len(B)
  out_list <- df_par_lapply_lb_progress(
    cl = cl,
    X = idxs,
    fun = worker_i,
    show_progress = show_progress
  )

  do.call(rbind, out_list)
}

############################################################
## 7) MC summary
############################################################

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
