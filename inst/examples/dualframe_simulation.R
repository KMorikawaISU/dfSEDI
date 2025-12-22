############################################################
## dualframe_simulation.R
##
## Dual-frame simulation + MC runner for dfSEDI
##
## - Scenarios S1/S2/S3 (paper-style DGP)
## - Fit: P / NP / NP_P / Eff / Eff_union_dat / Eff_S / Eff_P
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

extract_row <- function(fit, estimator, rep, Scenario, n_np, n_p, n_union) {
  # Robust row extractor: always returns a 1-row data.frame even if `fit` is NULL/partial
  scalar_num <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_real_)
    suppressWarnings(as.numeric(x)[1])
  }
  scalar_int <- function(x) {
    v <- scalar_num(x)
    if (!is.finite(v)) return(NA_integer_)
    as.integer(round(v))
  }
  scalar_chr <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_character_)
    as.character(x)[1]
  }

  if (is.null(fit) || !is.list(fit)) {
    fit <- list(
      theta = NA_real_, se = NA_real_, ci = c(NA_real_, NA_real_),
      phi = NULL, phi_se = NULL, phi_ci = NULL,
      error = "fit is NULL or not a list"
    )
  }

  ci <- fit$ci
  if (is.null(ci) || length(ci) < 2L) ci <- c(NA_real_, NA_real_)

  # Always output up to 6 phi components (unused ones are NA)
  p_phi_max <- 6L

  phi <- fit$phi
  phi_vec <- rep(NA_real_, p_phi_max)
  if (!is.null(phi) && length(phi) > 0L) {
    p_here <- min(length(phi), p_phi_max)
    phi_vec[seq_len(p_here)] <- as.numeric(phi)[seq_len(p_here)]
  }

  phi_se <- fit$phi_se
  phi_se_vec <- rep(NA_real_, p_phi_max)
  if (!is.null(phi_se) && length(phi_se) > 0L) {
    p_here <- min(length(phi_se), p_phi_max)
    phi_se_vec[seq_len(p_here)] <- as.numeric(phi_se)[seq_len(p_here)]
  }

  phi_ci <- fit$phi_ci
  phi_ci_l <- rep(NA_real_, p_phi_max)
  phi_ci_u <- rep(NA_real_, p_phi_max)
  if (!is.null(phi_ci) && is.matrix(phi_ci) && nrow(phi_ci) >= 2L && ncol(phi_ci) >= 1L) {
    p_here <- min(ncol(phi_ci), p_phi_max)
    phi_ci_l[seq_len(p_here)] <- as.numeric(phi_ci[1, seq_len(p_here)])
    phi_ci_u[seq_len(p_here)] <- as.numeric(phi_ci[2, seq_len(p_here)])
  }

  out <- data.frame(
    Scenario  = paste0("S", normalize_scenario(Scenario)),
    rep       = scalar_int(rep),
    estimator = as.character(estimator),
    theta     = scalar_num(fit$theta),
    se        = scalar_num(fit$se),
    ci_l      = scalar_num(ci[1]),
    ci_u      = scalar_num(ci[2]),
    n_np      = scalar_int(n_np),
    n_p       = scalar_int(n_p),
    n_union   = scalar_int(n_union),

    phi_1 = phi_vec[1], phi_2 = phi_vec[2], phi_3 = phi_vec[3],
    phi_4 = phi_vec[4], phi_5 = phi_vec[5], phi_6 = phi_vec[6],

    phi_se_1 = phi_se_vec[1], phi_se_2 = phi_se_vec[2], phi_se_3 = phi_se_vec[3],
    phi_se_4 = phi_se_vec[4], phi_se_5 = phi_se_vec[5], phi_se_6 = phi_se_vec[6],

    phi_ci_l_1 = phi_ci_l[1], phi_ci_l_2 = phi_ci_l[2], phi_ci_l_3 = phi_ci_l[3],
    phi_ci_l_4 = phi_ci_l[4], phi_ci_l_5 = phi_ci_l[5], phi_ci_l_6 = phi_ci_l[6],

    phi_ci_u_1 = phi_ci_u[1], phi_ci_u_2 = phi_ci_u[2], phi_ci_u_3 = phi_ci_u[3],
    phi_ci_u_4 = phi_ci_u[4], phi_ci_u_5 = phi_ci_u[5], phi_ci_u_6 = phi_ci_u[6],

    error     = scalar_chr(fit$error),
    stringsAsFactors = FALSE
  )
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

fit_all_estimators_once <- function(dat,
                                    Scenario = 1,
                                    K = 2,
                                    phi_start_true = TRUE,
                                    Eff_type = 2,
                                    progress_each = FALSE) {

  # --- choose phi_start (examples: provide "true" / near-true starts) ---
  if (isTRUE(phi_start_true)) {
    Scenario <- normalize_scenario(Scenario)

    if (Scenario == 1) {
      phi_start_NP   <- c(-2.15, -0.5, 0)              # for Chang-Kott base_fun (X only)
      phi_start_NP_P <- phi_start_NP
      phi_start_Eff  <- c(-2.15, -0.5, -0.75)          # true for Scenario 1 (L = (1, X, y))
    } else if (Scenario == 2) {
      phi_start_NP   <- c(-2.15, -0.5, 0)
      phi_start_NP_P <- phi_start_NP
      phi_start_Eff  <- c(-2.15, -0.5, -0.75)          # true for Scenario 2 (same NP model)
    } else if (Scenario == 3) {
      phi_start_NP   <- c(-1.5, 0, 0, 0)               # (not "true" in Scenario 3; just stable)
      phi_start_NP_P <- phi_start_NP
      phi_start_Eff  <- c(-1.5, 0, 0, 0)               # (not "true" in Scenario 3; just stable)
    }
  } else {
    phi_start_NP   <- NULL
    phi_start_NP_P <- NULL
    phi_start_Eff  <- NULL
  }

  # --- sample sizes ---
  n_np    <- sum(dat$d_np == 1, na.rm = TRUE)
  n_p     <- sum(dat$d_p == 1,  na.rm = TRUE)
  n_union <- sum((dat$d_np == 1 | dat$d_p == 1), na.rm = TRUE)

  # --- basis function for Chang & Kott-type estimators ---
  base_fun <- make_base_fun(Scenario)

  # --- P only ---
  fit_P <- safe_fit(
    df_estimate_P(dat = dat)
  )

  # --- NP only (Chang & Kott type) ---
  fit_NP <- safe_fit(
    df_estimate_NP(
      dat      = dat,
      base_fun = base_fun,
      phi_start = phi_start_NP
    )
  )

  # --- NP_P (Chang & Kott type) ---
  # In the example DGP we set pi_p = NA when d_p == 0. For NP_P we need finite pi_p to compute the union probability.
  # Here we impute pi_p by a simple *parametric* linear regression of 1/pi_p on L=(X, y) using P-sample units.
  dat_np_p <- dat

  if (any(!is.finite(dat_np_p$pi_p))) {

    get_X_mat <- function(dat_local) {
      if ("X" %in% names(dat_local)) {
        Xraw <- dat_local$X
      } else if ("x" %in% names(dat_local)) {
        Xraw <- dat_local$x
      } else {
        stop("dat must contain column 'X' (matrix/data.frame) or 'x' (numeric vector).")
      }

      # strip AsIs if present
      if (inherits(Xraw, "AsIs")) class(Xraw) <- setdiff(class(Xraw), "AsIs")

      if (is.data.frame(Xraw)) {
        mm <- stats::model.matrix(~ ., data = Xraw)
        if ("(Intercept)" %in% colnames(mm)) {
          mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
        }
        Xmat <- mm
      } else if (is.matrix(Xraw)) {
        Xmat <- Xraw
      } else {
        Xmat <- matrix(as.numeric(Xraw), ncol = 1)
      }

      Xmat <- as.matrix(Xmat)
      storage.mode(Xmat) <- "numeric"
      Xmat
    }

    X_all <- get_X_mat(dat_np_p)
    L_all <- cbind(X_all, y = as.numeric(dat_np_p$y))
    colnames(L_all) <- c(paste0("x", seq_len(ncol(L_all) - 1L)), "y")

    idx_obs <- which(dat_np_p$d_p == 1 & is.finite(dat_np_p$pi_p))
    idx_mis <- which(!is.finite(dat_np_p$pi_p))

    pi_obs <- dat_np_p$pi_p[idx_obs]
    mean_pi <- mean(pi_obs, na.rm = TRUE)
    if (!is.finite(mean_pi) || mean_pi <= 0 || mean_pi >= 1) mean_pi <- 1e-3

    if (length(idx_mis) > 0L) {
      if (length(idx_obs) >= 5L) {
        df_obs <- as.data.frame(L_all[idx_obs, , drop = FALSE])
        df_obs$inv_pi <- 1 / pmax(pi_obs, 1e-8)

        lm_fit <- try(stats::lm(inv_pi ~ ., data = df_obs), silent = TRUE)
        if (!inherits(lm_fit, "try-error")) {
          df_mis <- as.data.frame(L_all[idx_mis, , drop = FALSE])
          inv_hat <- try(stats::predict(lm_fit, newdata = df_mis), silent = TRUE)

          if (!inherits(inv_hat, "try-error")) {
            pi_hat <- 1 / inv_hat
            bad <- !is.finite(pi_hat) | pi_hat <= 0 | pi_hat >= 1
            pi_hat[bad] <- mean_pi
            pi_hat <- pmin(pmax(pi_hat, 1e-8), 1 - 1e-8)
            dat_np_p$pi_p[idx_mis] <- pi_hat
          } else {
            dat_np_p$pi_p[idx_mis] <- mean_pi
          }
        } else {
          dat_np_p$pi_p[idx_mis] <- mean_pi
        }
      } else {
        # Not enough P-sample units to fit a regression; fallback (user said 0 is OK).
        dat_np_p$pi_p[idx_mis] <- mean_pi
      }
    }

    # final safety
    dat_np_p$pi_p[!is.finite(dat_np_p$pi_p)] <- mean_pi
    dat_np_p$pi_p <- pmin(pmax(as.numeric(dat_np_p$pi_p), 1e-8), 1 - 1e-8)
  }

  fit_NP_P <- safe_fit(
    df_estimate_NP_P(
      dat      = dat_np_p,
      base_fun = base_fun,
      phi_start = phi_start_NP_P
    )
  )

  # --- Efficient (Eff): output both type=1 and type=2 ---
  fit_Eff_type1 <- safe_fit(
    Eff(
      dat      = dat,
      K        = K,
      type     = 1,
      phi_start = phi_start_Eff,
      x_info   = TRUE,
      progress = progress_each
    )
  )

  fit_Eff_type2 <- safe_fit(
    Eff(
      dat      = dat,
      K        = K,
      type     = 2,
      phi_start = phi_start_Eff,
      x_info   = TRUE,
      progress = progress_each
    )
  )

  # --- Efficient on union-only data (x_info=FALSE) ---
  dat_union <- subset(dat, d_np == 1 | d_p == 1)

  eff_union_args <- list(
    dat      = dat_union,
    K        = K,
    phi_start = phi_start_Eff,
    x_info   = FALSE,
    progress = progress_each
  )
  if ("N" %in% names(formals(Eff))) eff_union_args$N <- nrow(dat)

  eff_union_args$type <- 1
  fit_Eff_union_type1 <- safe_fit(do.call(Eff, eff_union_args))

  eff_union_args$type <- 2
  fit_Eff_union_type2 <- safe_fit(do.call(Eff, eff_union_args))

  # --- Sub-efficient + parametric efficient ---
  fit_EffS <- safe_fit(
    Eff_S(
      dat      = dat,
      K        = K,
      x_info   = TRUE,
      progress = progress_each
    )
  )

  fit_EffP <- safe_fit(
    Eff_P(
      dat    = dat,
      x_info = TRUE
    )
  )

  # Backward-compatible aliases (so older example code using fits$Eff / fits$Eff_union still works)
  fit_Eff_default <- if (isTRUE(Eff_type == 1)) fit_Eff_type1 else fit_Eff_type2
  fit_Eff_union_default <- if (isTRUE(Eff_type == 1)) fit_Eff_union_type1 else fit_Eff_union_type2

  list(
    P         = fit_P,
    NP        = fit_NP,
    NP_P      = fit_NP_P,

    Eff       = fit_Eff_default,
    Eff_union = fit_Eff_union_default,

    Eff_type1 = fit_Eff_type1,
    Eff_type2 = fit_Eff_type2,

    Eff_union_type1 = fit_Eff_union_type1,
    Eff_union_type2 = fit_Eff_union_type2,

    Eff_S     = fit_EffS,
    Eff_P     = fit_EffP,

    n_np      = n_np,
    n_p       = n_p,
    n_union   = n_union
  )
}



############################################################
## 4) One MC replication -> long-format rows
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
    extract_row(fits$P, "P", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP, "NP", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP_P, "NP_P", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_type1, "Eff_type1", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_type2, "Eff_type2", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_union_dat_type1, "Eff_union_dat_type1", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_union_dat_type2, "Eff_union_dat_type2", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_S, "Eff_S", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_P, "Eff_P", rep_id, Scenario, fits$n_np, fits$n_p, fits$n_union)
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
run_mc <- function(
    B = 100,
    N = 10000,
    Scenario = 1,
    K = 2,
    seed_start = 1,
    show_progress = TRUE,
    progress_each_fit = FALSE,
    pi_p_offset = 0.005,
    parallel = c("none", "multicore", "snow"),
    n_cores = max(1L, parallel::detectCores() - 1L),
    ...
) {

  parallel <- match.arg(parallel)
  sc <- normalize_scenario(Scenario)

  seeds   <- seed_start + seq_len(B) - 1L
  rep_ids <- seq_len(B)

  estimators_all <- c(
    "P", "NP", "NP_P",
    "Eff_type1", "Eff_type2",
    "Eff_union_dat_type1", "Eff_union_dat_type2",
    "Eff_S", "Eff_P"
  )

  make_error_fit <- function(msg) {
    list(
      theta = NA_real_,
      se    = NA_real_,
      ci    = c(NA_real_, NA_real_),
      phi   = NULL,
      phi_se = NULL,
      phi_ci = NULL,
      error = msg
    )
  }

  make_error_rows <- function(rep_id, sc, msg) {
    do.call(rbind, lapply(estimators_all, function(est) {
      extract_row(make_error_fit(msg), est, rep_id, sc, NA_real_, NA_real_, NA_real_)
    }))
  }

  worker_i <- function(i) {
    tryCatch({
      mc_one_rep_long(
        rep_id        = rep_ids[i],
        seed          = seeds[i],
        N             = N,
        Scenario      = sc,
        K             = K,
        progress_each = progress_each_fit,
        pi_p_offset   = pi_p_offset,
        ...
      )
    }, error = function(e) {
      make_error_rows(rep_ids[i], sc, conditionMessage(e))
    })
  }

  idxs <- seq_len(B)

  if (parallel == "none") {
    if (isTRUE(show_progress)) {
      pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      on.exit(close(pb), add = TRUE)
    }

    out_list <- vector("list", B)
    for (i in idxs) {
      out_list[[i]] <- worker_i(i)
      if (isTRUE(show_progress)) utils::setTxtProgressBar(pb, i)
    }

  } else if (parallel == "multicore") {

    message("parallel='multicore': show_progress is ignored (no progress bar in base R for mclapply).")

    out_list <- parallel::mclapply(
      idxs,
      worker_i,
      mc.cores = n_cores,
      mc.preschedule = FALSE
    )

  } else if (parallel == "snow") {

    n_workers <- max(1L, as.integer(n_cores))
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, library(dfSEDI))

    # Export everything needed by `worker_i()` and its callees
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
        "progress_each_fit",
        "pi_p_offset",
        "estimators_all",
        "make_error_fit",
        "make_error_rows",
        "worker_i"
      ),
      envir = environment()
    )

    if (isTRUE(show_progress)) {
      out_list <- df_par_lapply_lb_progress(cl, idxs, worker_i)
    } else {
      out_list <- parallel::parLapplyLB(cl, idxs, worker_i)
    }

  }

  res <- do.call(rbind, out_list)
  rownames(res) <- NULL
  res
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
