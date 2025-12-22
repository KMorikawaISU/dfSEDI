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

extract_row <- function(fit, estimator, rep, Scenario, n_np, n_p, n_union, p_phi_max = 4L) {
  ci <- fit$ci
  if (is.null(ci) || length(ci) != 2) ci <- c(NA_real_, NA_real_)

  # phi (fixed-width output up to p_phi_max)
  phi     <- rep(NA_real_, p_phi_max)
  phi_se  <- rep(NA_real_, p_phi_max)
  phi_ci_l <- rep(NA_real_, p_phi_max)
  phi_ci_u <- rep(NA_real_, p_phi_max)

  if (!is.null(fit$phi)) {
    phi0 <- as.numeric(fit$phi)
    k <- min(length(phi0), p_phi_max)
    if (k > 0) phi[1:k] <- phi0[1:k]
  }

  if (!is.null(fit$phi_se)) {
    se0 <- as.numeric(fit$phi_se)
    k <- min(length(se0), p_phi_max)
    if (k > 0) phi_se[1:k] <- se0[1:k]
  }

  if (!is.null(fit$phi_ci) && is.matrix(fit$phi_ci) && nrow(fit$phi_ci) == 2) {
    lo0 <- as.numeric(fit$phi_ci[1, ])
    up0 <- as.numeric(fit$phi_ci[2, ])
    k <- min(length(lo0), length(up0), p_phi_max)
    if (k > 0) {
      phi_ci_l[1:k] <- lo0[1:k]
      phi_ci_u[1:k] <- up0[1:k]
    }
  }

  out <- data.frame(
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

  # Append phi columns (fixed width)
  for (j in seq_len(p_phi_max)) out[[paste0("phi_", j)]] <- phi[j]
  for (j in seq_len(p_phi_max)) out[[paste0("phi_se_", j)]] <- phi_se[j]
  for (j in seq_len(p_phi_max)) out[[paste0("phi_ci_l_", j)]] <- phi_ci_l[j]
  for (j in seq_len(p_phi_max)) out[[paste0("phi_ci_u_", j)]] <- phi_ci_u[j]

  out
}


# Create a "failed replication" long-format data.frame (all NA rows),
# used as a robust fallback in parallel execution when a worker crashes.
make_na_fit <- function(msg = NA_character_, p_phi_max = 4L) {
  list(
    theta  = NA_real_,
    se     = NA_real_,
    ci     = c(NA_real_, NA_real_),
    phi    = rep(NA_real_, p_phi_max),
    phi_se = rep(NA_real_, p_phi_max),
    phi_ci = matrix(NA_real_, nrow = 2, ncol = p_phi_max,
                    dimnames = list(c("lower","upper"), NULL)),
    error  = msg
  )
}

make_na_rep_long <- function(rep_id,
                             Scenario,
                             n_np = NA_integer_,
                             n_p = NA_integer_,
                             n_union = NA_integer_,
                             p_phi_max = 4L,
                             err = NA_character_) {

  ests <- c(
    "P",
    "NP",
    "NP_P",
    "Eff_type1",
    "Eff_type2",
    "Eff_union_dat_type1",
    "Eff_union_dat_type2",
    "Eff_S",
    "Eff_P"
  )

  fits <- lapply(ests, function(e) make_na_fit(msg = err, p_phi_max = p_phi_max))
  names(fits) <- ests

  do.call(
    rbind,
    lapply(ests, function(e) {
      extract_row(fits[[e]], e, rep_id, Scenario, n_np, n_p, n_union, p_phi_max = p_phi_max)
    })
  )
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
                                    Eff_type = 2,
                                    x_info = TRUE,
                                    progress_each = FALSE) {

  Scenario <- normalize_scenario(Scenario)

  # Scenario-specific basis function (for Chang & Kott-style NP/NP_P)
  base_fun <- make_base_fun(Scenario)

  # "True" phi as starting values (used in the examples)
  # - For NP/NP_P: phi corresponds to base_fun(X) (dimension p_x + 2)
  # - For Eff: phi corresponds to L = (1, X, y) (dimension 1 + p_x + 1)
  X <- as.matrix(dat$X)
  p_x <- ncol(X)
  p_phi_np <- ncol(base_fun(X))

  phi_start_NP <- rep(0, p_phi_np)
  phi_start_Eff <- rep(0, 1 + p_x + 1)

  if (Scenario %in% c(1, 2)) {
    # DGP uses lin_np = -2.15 - 0.5*x - 0.75*y (Scenario 1/2)
    # NP/NP_P base_fun does not include y, so we use the x-coefficient part as a start.
    phi_start_NP[1] <- -2.15
    if (p_phi_np >= 2) phi_start_NP[2] <- -0.5

    # Eff uses (1, X, y)
    phi_start_Eff <- c(-2.15, -0.5, -0.75)
  } else if (Scenario == 3) {
    # DGP uses lin_np = -1.5 - 0.3*cos(2*y) - 0.1*(y-1)^2 (Scenario 3)
    # As a simple start, use only the intercept.
    phi_start_NP[1] <- -1.5
    phi_start_Eff[1] <- -1.5
  }

  # Safe wrapper
  safe_fit <- function(expr) {
    out <- try(expr, silent = TRUE)
    if (inherits(out, "try-error")) {
      return(list(phi = NA, theta = NA_real_, var = NA_real_, se = NA_real_, ci = c(NA_real_, NA_real_)))
    }
    out
  }

  # Union data (only units observed in either sample)
  dat_union <- dat[as.numeric(dat$d_p) == 1 | as.numeric(dat$d_np) == 1, , drop = FALSE]

  # --- Estimators ---
  fit_P <- safe_fit(df_estimate_P(dat))

  fit_NP <- safe_fit(df_estimate_NP(dat, base_fun = base_fun, phi_start = phi_start_NP))

  # NP_P needs pi_p for union-sample units; in this simulation, pi_p is NA when d_p == 0,
  # so we impute pi_p for (d_np == 1, d_p == 0) before calling df_estimate_NP_P.
  dat_np_p <- dat

  # NP_P needs pi_p for union-sample units; in realistic data pi_p is missing when d_p == 0.
  # We impute pi_p for (d_np, d_p) = (1, 0) units using a simple *parametric* propensity model for d_p:
  #   - if Y is continuous : linear regression  (LPM) for d_p ~ (X, y)
  #   - if Y is binary     : logistic regression for d_p ~ (X, y)
  # We do this with cross-fitting over K folds (for stability), and then fill pi_p in the target cells.
  # If that fails, we fall back to impute_pi_p_crossfit() when available; otherwise mean imputation.
  mis_pi_p <- which(dat_np_p$d_np == 1 & dat_np_p$d_p == 0 & !is.finite(as.numeric(dat_np_p$pi_p)))
  if (length(mis_pi_p) > 0) {
    imputed <- FALSE

    # --- 1) Parametric cross-fit imputation ---
    imputed <- tryCatch({
      yv <- as.numeric(dat_np_p$y)
      y_uniq <- sort(unique(yv[is.finite(yv)]))
      is_bin_y <- (length(y_uniq) <= 2) && all(y_uniq %in% c(0, 1))

      X <- dat_np_p$X
      X_df <- if (is.matrix(X)) as.data.frame(X) else if (is.data.frame(X)) X else data.frame(x = X)
      if (is.null(colnames(X_df))) colnames(X_df) <- paste0("x", seq_len(ncol(X_df)))
      df_mod <- data.frame(d_p = as.numeric(dat_np_p$d_p), y = yv, X_df)

      n <- nrow(df_mod)
      if (!is.finite(K) || K < 2) {
        folds <- rep(1L, n)
      } else {
        folds <- sample(rep(seq_len(K), length.out = n))
      }

      pred <- rep(NA_real_, n)
      for (k in sort(unique(folds))) {
        te <- which(folds == k)
        tr <- which(folds != k)
        # If K < 2, this is just the full-sample fit (tr == all, te == all).
        fit_k <- tryCatch({
          if (is_bin_y) {
            suppressWarnings(stats::glm(d_p ~ ., data = df_mod[tr, , drop = FALSE], family = stats::binomial()))
          } else {
            stats::lm(d_p ~ ., data = df_mod[tr, , drop = FALSE])
          }
        }, error = function(e) NULL)

        if (!is.null(fit_k)) {
          if (is_bin_y) {
            pred[te] <- as.numeric(stats::predict(fit_k, newdata = df_mod[te, , drop = FALSE], type = "response"))
          } else {
            pred[te] <- as.numeric(stats::predict(fit_k, newdata = df_mod[te, , drop = FALSE]))
          }
        }
      }

      # If some folds failed, fall back to a global fit for those rows.
      bad <- which(!is.finite(pred))
      if (length(bad) > 0) {
        fit_g <- tryCatch({
          if (is_bin_y) {
            suppressWarnings(stats::glm(d_p ~ ., data = df_mod, family = stats::binomial()))
          } else {
            stats::lm(d_p ~ ., data = df_mod)
          }
        }, error = function(e) NULL)
        if (!is.null(fit_g)) {
          if (is_bin_y) {
            pred[bad] <- as.numeric(stats::predict(fit_g, newdata = df_mod[bad, , drop = FALSE], type = "response"))
          } else {
            pred[bad] <- as.numeric(stats::predict(fit_g, newdata = df_mod[bad, , drop = FALSE]))
          }
        }
      }

      pred <- pmin(pmax(pred, 1e-6), 1 - 1e-6)
      dat_np_p$pi_p[mis_pi_p] <- pred[mis_pi_p]
      TRUE
    }, error = function(e) FALSE)

    # --- 2) Fallback: kernel cross-fit imputation (if available) ---
    if (!imputed || any(!is.finite(as.numeric(dat_np_p$pi_p[mis_pi_p])))) {
      if (exists("impute_pi_p_crossfit", mode = "function")) {
        dat_np_p <- impute_pi_p_crossfit(dat_np_p, K = K, sigma = NULL)
        imputed <- TRUE
      }
    }

    # --- 3) Fallback: mean imputation / zero ---
    if (any(!is.finite(as.numeric(dat_np_p$pi_p[mis_pi_p])))) {
      mu_pi <- mean(as.numeric(dat_np_p$pi_p[is.finite(as.numeric(dat_np_p$pi_p))]), na.rm = TRUE)
      if (!is.finite(mu_pi)) mu_pi <- 0
      dat_np_p$pi_p[mis_pi_p] <- mu_pi
    }
  }

  fit_NP_P <- safe_fit(df_estimate_NP_P(dat_np_p, base_fun = base_fun, phi_start = phi_start_NP_P))

  # Eff: output BOTH types
  fit_Eff_type1 <- safe_fit(Eff(dat, K = K, type = 1, x_info = x_info, phi_start = phi_start_Eff, progress = progress_each))
  fit_Eff_type2 <- safe_fit(Eff(dat, K = K, type = 2, x_info = x_info, phi_start = phi_start_Eff, progress = progress_each))

  # Eff on union-only data (x_info = FALSE)
  fit_Eff_union_dat_type1 <- safe_fit(Eff(dat_union, K = K, type = 1, x_info = FALSE, phi_start = phi_start_Eff, progress = progress_each))
  fit_Eff_union_dat_type2 <- safe_fit(Eff(dat_union, K = K, type = 2, x_info = FALSE, phi_start = phi_start_Eff, progress = progress_each))

  # Oracles / semi-oracles
  fit_Eff_S <- safe_fit(efficient_estimator_semioracle(dat))
  fit_Eff_P <- safe_fit(efficient_estimator_oracle_p(dat))

  list(
    P = fit_P,
    NP = fit_NP,
    NP_P = fit_NP_P,
    Eff_type1 = fit_Eff_type1,
    Eff_type2 = fit_Eff_type2,
    Eff_union_dat_type1 = fit_Eff_union_dat_type1,
    Eff_union_dat_type2 = fit_Eff_union_dat_type2,
    Eff_S = fit_Eff_S,
    Eff_P = fit_Eff_P
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
    tryCatch(
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
      ),
      error = function(e) {
        make_na_rep_long(
          rep_id   = rep_ids[i],
          Scenario = sc,
          p_phi_max = 4L,
          err = paste0("MC replication failed: ", conditionMessage(e))
        )
      }
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

    out_list <- parallel::mclapply(
      X = seq_len(B),
      FUN = worker_i,
      mc.cores = n_cores,
      mc.preschedule = FALSE
    )

    # Robust fallback: if a forked worker crashes, mclapply can return NULL elements
    # and emit "jobs did not deliver results". Replace such elements by NA rows
    # so downstream code (summarize_mc) still works.
    if (length(out_list) < B) {
      if (length(out_list) == 0L) {
        out_list <- vector("list", B)
      }
      if ((length(out_list) + 1L) <= B) {
        for (i in (length(out_list) + 1L):B) {
          out_list[[i]] <- make_na_rep_long(
            rep_id   = rep_ids[i],
            Scenario = sc,
            p_phi_max = 4L,
            err = "multicore job did not deliver results"
          )
        }
      }
    }

    for (i in seq_len(B)) {
      if (is.null(out_list[[i]]) || !is.data.frame(out_list[[i]]) || nrow(out_list[[i]]) == 0L) {
        # Try re-running this replication serially on the master process.
        # This often succeeds when the failure was specific to fork-based parallelism.
        out_i <- tryCatch(worker_i(i), error = function(e) NULL)

        if (!is.null(out_i) && is.data.frame(out_i) && nrow(out_i) > 0L) {
          out_list[[i]] <- out_i
        } else {
          out_list[[i]] <- make_na_rep_long(
            rep_id   = rep_ids[i],
            Scenario = sc,
            p_phi_max = 4L,
            err = "multicore job did not deliver results"
          )
        }
      }
    }

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
      "make_na_fit",
      "make_na_rep_long",
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

  idxs <- seq_len(B)
  out_list <- df_par_lapply_lb_progress(
    cl = cl,
    X = idxs,
    fun = worker_i,
    show_progress = show_progress
  )

  # Robust fallback: if a worker dies, fill its result with NA rows
  if (length(out_list) < B) {
    if (length(out_list) == 0L) {
      out_list <- vector("list", B)
    }
    if ((length(out_list) + 1L) <= B) {
      for (i in (length(out_list) + 1L):B) {
        out_list[[i]] <- make_na_rep_long(
          rep_id   = rep_ids[i],
          Scenario = sc,
          p_phi_max = 4L,
          err = "snow worker did not deliver results"
        )
      }
    }
  }

  for (i in seq_len(B)) {
    if (is.null(out_list[[i]]) || !is.data.frame(out_list[[i]]) || nrow(out_list[[i]]) == 0L) {
      # Try re-running serially on the master process (PSOCK workers can still die).
      out_i <- tryCatch(worker_i(i), error = function(e) NULL)

      if (!is.null(out_i) && is.data.frame(out_i) && nrow(out_i) > 0L) {
        out_list[[i]] <- out_i
      } else {
        out_list[[i]] <- make_na_rep_long(
          rep_id   = rep_ids[i],
          Scenario = sc,
          p_phi_max = 4L,
          err = "snow worker did not deliver results"
        )
      }
    }
  }

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
