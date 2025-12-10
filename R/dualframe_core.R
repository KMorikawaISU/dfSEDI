############################################################
## dualframe_core.R
##
## Core estimation functions for dual-frame data integration:
##   - Eff     : semiparametric efficient estimator (DML2, K-fold)
##   - Eff_S   : sub-efficient estimator (Remark 6)
##   - Eff_P   : parametric efficient estimator (working model)
##
## Internal helpers:
##   - P, NP, NP+P estimators (df_estimate_P, df_estimate_NP, df_estimate_NP_P)
##   - DML2 cross-fitting machinery
##   - Sandwich variance / SE / 95% CI from per-observation contributions
##
## IMPORTANT:
##   - x can be multi-dimensional.
##   - User specifies which columns of dat are covariates via `x_cols`.
##   - All design matrices use X = dat[ , x_cols], y = dat$y.
############################################################

############################################################
## Helper: extract X as a numeric matrix
############################################################

#' Extract covariate matrix X from dat
#'
#' @param dat    data.frame with covariates and y, d_np, d_p, pi_p, pi_np
#' @param x_cols character or integer vector indicating covariate columns
#' @return numeric matrix with nrow(dat) rows
df_get_X <- function(dat, x_cols) {
  if (is.null(x_cols)) {
    if ("x" %in% names(dat)) {
      X <- dat[, "x", drop = FALSE]
    } else {
      stop("x_cols is NULL and no column named 'x' found in dat.")
    }
  } else if (is.character(x_cols)) {
    missing <- setdiff(x_cols, names(dat))
    if (length(missing) > 0L) {
      stop("The following x_cols are not in dat: ",
           paste(missing, collapse = ", "))
    }
    X <- dat[, x_cols, drop = FALSE]
  } else {
    ## assume numeric column indices
    X <- dat[, x_cols, drop = FALSE]
  }
  as.matrix(X)
}

############################################################
## Sandwich variance + 95% CI from contributions
############################################################

#' Compute sandwich variance and 95% CI from per-observation contributions
#'
#' @param contrib Numeric vector of length n (influence / pseudo-outcome).
#' @param level   Confidence level (default 0.95).
#' @return list(theta, var, se, ci)
df_sandwich_from_contrib <- function(contrib, level = 0.95) {
  n <- length(contrib)
  if (n <= 1L || all(!is.finite(contrib))) {
    return(list(
      theta = NA_real_,
      var   = NA_real_,
      se    = NA_real_,
      ci    = c(NA_real_, NA_real_)
    ))
  }

  theta_hat <- mean(contrib)
  s2        <- stats::var(contrib)         # sample variance (n-1 denominator)
  var_hat   <- s2 / n                      # Var(mean(contrib))
  se_hat    <- sqrt(var_hat)
  alpha     <- 1 - level
  z         <- stats::qnorm(1 - alpha / 2)
  ci        <- c(theta_hat - z * se_hat,
                 theta_hat + z * se_hat)

  list(theta = theta_hat,
       var   = var_hat,
       se    = se_hat,
       ci    = ci)
}

############################################################
## Basis functions (for Chang & Kott estimating eq.)
## We use g2 = g3 = g4 = (1, X, y) for general dimension.
############################################################

df_g2 <- function(dat, x_cols) {
  X <- df_get_X(dat, x_cols)
  cbind(1, X, dat$y)
}
df_g3 <- df_g2
df_g4 <- df_g2

############################################################
## Chang & Kott type estimating equation for pi_NP(phi)
############################################################

pi_np.est_simple <- function(dat, h4, x_cols, p.use = TRUE) {
  function(phi) {
    X    <- df_get_X(dat, x_cols)
    y    <- dat$y
    l    <- cbind(1, X, y)      # logistic design for pi_NP
    pi_p <- dat$pi_p
    pi_np <- 1 / (1 + exp(-l %*% phi))

    d_p  <- dat$d_p
    d_np <- dat$d_np

    if (p.use) {
      d_set4 <- matrix(
        1 - (d_p + d_np - d_p * d_np) /
          (pi_p + pi_np - pi_p * pi_np),
        nrow(dat), 1
      )
    } else {
      d_set4 <- matrix(1 - d_np / pi_np, nrow(dat), 1)
    }

    H4     <- h4(dat, x_cols)              # (n x p_phi)
    est_eq <- as.vector(crossprod(H4, d_set4))
    est_eq
  }
}

############################################################
## Efficient score helper functions (Section 4)
############################################################

h4_prob_denom_function <- function(pi_np, pi_p, phi) {
  pi_np_p <- pi_p + pi_np - pi_p * pi_np
  (1 - pi_np_p) / (pi_np_p)^2
}

h4_prob_numer_function <- function(pi_np, pi_p, phi, y) {
  pi_np_p <- pi_p + pi_np - pi_p * pi_np
  - y * (1 - pi_np_p) / (pi_np_p)^2
}

eta4_prob_numer_function <- function(pi_np, pi_p, phi, l) {
  pi_np_p <- pi_p + pi_np - pi_p * pi_np
  - l * pi_np * (1 - pi_np_p) / (pi_np_p)^2
}

############################################################
## Kernel regression: E(Y|X), pi_P, eta4*, h4*
############################################################

# Nonparametric regression mu(x) = E[Y | X = x] using P-only data
regression_expectation_kernlab <- function(dat, newdata, x_cols, sigma = NULL) {
  X_all <- df_get_X(dat, x_cols)
  idx   <- dat$d_p == 1 & dat$d_np == 0
  X_obs <- X_all[idx, , drop = FALSE]
  y_obs <- dat$y[idx]

  if (length(y_obs) < 2L) {
    return(rep(mean(y_obs), nrow(newdata)))
  }

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(X_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)
  reg_model  <- kernlab::gausspr(x = X_obs,
                                 y = y_obs,
                                 kernel = rbf_kernel)

  X_new   <- df_get_X(newdata, x_cols)
  reg_pred <- kernlab::predict(reg_model, X_new)
  as.numeric(reg_pred)
}

# Nonparametric regression for pi_P (inverse probability regression)
pi_p_estimation_kernlab <- function(dat, new_l, x_cols, sigma = NULL) {
  X_all <- df_get_X(dat, x_cols)
  idx   <- dat$d_p == 1
  X_obs <- X_all[idx, , drop = FALSE]
  y_obs <- dat$y[idx]
  l_obs <- cbind(X_obs, y_obs)
  pi_p_obs <- dat$pi_p[idx]

  if (length(pi_p_obs) < 2L) {
    return(rep(mean(pi_p_obs), nrow(new_l)))
  }

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(l_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)
  reg_model  <- kernlab::gausspr(x = l_obs,
                                 y = 1 / pi_p_obs,
                                 kernel = rbf_kernel)
  reg_pred <- kernlab::predict(reg_model, as.matrix(new_l))

  as.numeric(1 / reg_pred)
}

# Conditional expectation for eta4*(L;phi) given X (for efficient phi)
estimate_conditional_expectation_kernlab_phi <- function(dat, phi, newdata,
                                                         x_cols, sigma = NULL) {
  X_all <- df_get_X(dat, x_cols)
  idx   <- dat$d_p == 1 | dat$d_np == 1
  X_obs <- X_all[idx, , drop = FALSE]
  y_obs <- dat$y[idx]
  l_obs <- cbind(1, X_obs, y_obs)
  pi_p_obs <- dat$pi_p[idx]
  pi_np_obs <- 1 / (1 + exp(-l_obs %*% phi))

  Denom_vals <- h4_prob_denom_function(pi_np_obs, pi_p_obs, phi)

  n_obs <- nrow(l_obs)
  p_dim <- ncol(l_obs)
  Numer_mat <- matrix(NA_real_, nrow = n_obs, ncol = p_dim)
  for (i in seq_len(n_obs)) {
    Numer_mat[i, ] <- eta4_prob_numer_function(
      pi_np_obs[i], pi_p_obs[i], phi = phi, l = l_obs[i, ]
    )
  }

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(X_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)

  eta4_denom_model <- kernlab::gausspr(x = X_obs,
                                       y = Denom_vals,
                                       kernel = rbf_kernel)
  eta4_numer_model <- lapply(
    seq_len(p_dim),
    function(j) kernlab::gausspr(
      x = X_obs,
      y = Numer_mat[, j],
      kernel = rbf_kernel
    )
  )

  X_new <- df_get_X(newdata, x_cols)
  eta4_denom_pred <- kernlab::predict(eta4_denom_model, X_new)
  eta4_numer_pred <- sapply(
    eta4_numer_model,
    function(m) kernlab::predict(m, X_new)
  )

  sweep(eta4_numer_pred, 1, eta4_denom_pred, FUN = "/")
}

# Conditional expectation for h4*(X;phi) given X (for efficient theta)
estimate_conditional_expectation_kernlab_theta <- function(dat, phi, newdata,
                                                           x_cols, sigma = NULL) {
  X_all <- df_get_X(dat, x_cols)
  idx   <- dat$d_p == 1 | dat$d_np == 1
  X_obs <- X_all[idx, , drop = FALSE]
  y_obs <- dat$y[idx]
  l_obs <- cbind(1, X_obs, y_obs)
  pi_p_obs <- dat$pi_p[idx]
  pi_np_obs <- 1 / (1 + exp(-l_obs %*% phi))

  Denom_vals <- h4_prob_denom_function(pi_np_obs, pi_p_obs, phi)
  Numer_vals <- h4_prob_numer_function(pi_np_obs, pi_p_obs, phi, y_obs)

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(X_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)

  h4_denom_model <- kernlab::gausspr(x = X_obs,
                                     y = Denom_vals,
                                     kernel = rbf_kernel)
  h4_numer_model <- kernlab::gausspr(x = X_obs,
                                     y = Numer_vals,
                                     kernel = rbf_kernel)

  X_new <- df_get_X(newdata, x_cols)
  h4_denom_pred <- kernlab::predict(h4_denom_model, X_new)
  h4_numer_pred <- kernlab::predict(h4_numer_model, X_new)

  as.numeric(h4_numer_pred / h4_denom_pred)
}

############################################################
## Efficient phi estimating equation (Theorem 2)
############################################################

estimating_equation_optimal_phi <- function(dat1, dat2, x_cols,
                                            eta4_star = NULL) {
  X1   <- df_get_X(dat1, x_cols)
  y1   <- dat1$y
  d_p1 <- dat1$d_p
  d_np1 <- dat1$d_np
  pi_p1 <- dat1$pi_p
  l1    <- cbind(1, X1, y1)

  function(phi) {
    if (is.null(eta4_star)) {
      eta4_star_local <- estimate_conditional_expectation_kernlab_phi(
        dat2, phi, dat1, x_cols = x_cols
      )
    } else {
      eta4_star_local <- eta4_star
    }

    pi_np1   <- as.numeric(1 / (1 + exp(-l1 %*% phi)))
    pi_np_p1 <- pi_np1 + pi_p1 - pi_np1 * pi_p1

    term1 <- l1 * pi_np1 * (1 - pi_np1) / pi_np_p1 *
      (d_np1 / pi_np1 * (pi_p1 - d_p1) -
         d_p1 / (1 - pi_np1) * (1 - d_np1 / pi_np1))

    term2_coef <- 1 -
      (d_p1 * d_np1) / (pi_p1 * pi_np1) -
      1 / pi_np_p1 *
      (d_np1 * (1 - d_p1 / pi_p1) + d_p1 * (1 - d_np1 / pi_np1))

    est_eq <- term1 + term2_coef * eta4_star_local

    apply(est_eq, 2, mean)
  }
}

############################################################
## Efficient theta contributions (given phi, h4*)
############################################################

# Per-observation contributions for efficient theta estimator
efficient_theta_contrib <- function(dat1, dat2, phi, x_cols, h4_star = NULL) {
  X1   <- df_get_X(dat1, x_cols)
  y1   <- dat1$y
  d_p1 <- dat1$d_p
  d_np1 <- dat1$d_np
  pi_p1 <- dat1$pi_p
  l1    <- cbind(1, X1, y1)

  pi_np1   <- as.numeric(1 / (1 + exp(-l1 %*% phi)))
  pi_np_p1 <- pi_np1 + pi_p1 - pi_np1 * pi_p1

  if (is.null(h4_star)) {
    h4_star_local <- estimate_conditional_expectation_kernlab_theta(
      dat2, phi, dat1, x_cols = x_cols
    )
  } else {
    h4_star_local <- h4_star
  }

  term1 <- y1 * (d_p1 / pi_p1 +
                   (1 - d_p1 / pi_p1) / pi_np_p1 *
                   (d_np1 - d_p1 * (d_np1 - pi_np1)))

  term2_coef <- 1 -
    (d_p1 * d_np1) / (pi_p1 * pi_np1) -
    1 / pi_np_p1 *
    (d_np1 * (1 - d_p1 / pi_p1) + d_p1 * (1 - d_np1 / pi_np1))

  contrib <- term1 - term2_coef * h4_star_local
  contrib
}

# Simple wrapper: returns only mean, kept for internal use
optimal_theta <- function(dat1, dat2, phi, x_cols, h4_star = NULL) {
  contrib <- efficient_theta_contrib(dat1, dat2, phi, x_cols, h4_star)
  mean(contrib)
}

############################################################
## Sub-efficient theta contributions (Remark 6)
############################################################

subefficient_contrib <- function(dat, mu_hat) {
  d_np <- dat$d_np
  d_p  <- dat$d_p
  y    <- dat$y
  pi_p <- dat$pi_p

  pseudo <- d_np * y +
    (1 - d_np) * (d_p / pi_p * y + (1 - d_p / pi_p) * mu_hat)

  pseudo
}

subefficient_score_theta <- function(dat, mu_hat) {
  mean(subefficient_contrib(dat, mu_hat))
}

############################################################
## DML2 helpers: K-fold split & pi_P cross-fitting
############################################################

make_folds <- function(n, K) {
  fold_id <- sample(rep(1:K, length.out = n))
  split(seq_len(n), fold_id)
}

impute_pi_p_crossfit <- function(dat, folds, x_cols,
                                 sigma = NULL, progress = FALSE) {
  dat_cf <- dat

  if (progress) {
    cat("Cross-fitting pi_P ...\n")
    pb <- utils::txtProgressBar(min = 0, max = length(folds),
                                style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (k in seq_along(folds)) {
    idx_k     <- folds[[k]]
    dat_fold  <- dat_cf[idx_k, ]
    dat_train <- dat_cf[-idx_k, ]

    X_fold <- df_get_X(dat_fold, x_cols)
    l_fold <- cbind(X_fold, dat_fold$y)
    d_p_fold  <- dat_fold$d_p
    d_np_fold <- dat_fold$d_np

    pi_p_mis_local <- which(d_np_fold == 1 & d_p_fold == 0)
    if (length(pi_p_mis_local) > 0) {
      tilde_pi_p <- pi_p_estimation_kernlab(
        dat_train,
        l_fold[pi_p_mis_local, , drop = FALSE],
        x_cols = x_cols,
        sigma  = sigma
      )
      dat_cf$pi_p[idx_k[pi_p_mis_local]] <- tilde_pi_p
    }

    if (progress) {
      utils::setTxtProgressBar(pb, k)
    }
  }
  dat_cf
}

############################################################
## Internal P / NP / NP+P estimators
############################################################

# P-only estimator (Horvitz–Thompson type)
df_estimate_P <- function(dat) {
  contrib <- dat$d_p * dat$y / dat$pi_p
  df_sandwich_from_contrib(contrib)
}

# NP-only estimator (Chang & Kott, p.use = FALSE)
df_estimate_NP <- function(dat,
                           x_cols,
                           phi_start = NULL,
                           max_iter  = 20) {

  if (is.null(phi_start)) {
    p_x  <- ncol(df_get_X(dat, x_cols))
    phi_start <- c(-log(1 / mean(dat$d_np) - 1),
                   rep(0, p_x),
                   0)
  }

  phi_est <- list(termcd = 99)
  k <- 0
  while (phi_est$termcd > 2 && k < max_iter) {
    init <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
    phi_est <- nleqslv::nleqslv(
      init,
      pi_np.est_simple(dat, df_g4, x_cols, p.use = FALSE)
    )
    k <- k + 1
  }

  if (phi_est$termcd > 2) {
    phi_hat   <- rep(NA_real_, length(phi_start))
    contrib   <- rep(NA_real_, nrow(dat))
    theta_res <- df_sandwich_from_contrib(contrib)
  } else {
    phi_hat <- phi_est$x
    X       <- df_get_X(dat, x_cols)
    l       <- cbind(1, X, dat$y)
    pi_np   <- as.numeric(1 / (1 + exp(-l %*% phi_hat)))
    contrib <- dat$d_np * dat$y / pi_np
    theta_res <- df_sandwich_from_contrib(contrib)
  }

  list(phi = c(phi_hat),
       theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci)
}

# NP∪P estimator (Chang & Kott, p.use = TRUE)
df_estimate_NP_P <- function(dat,
                             x_cols,
                             phi_start = NULL,
                             max_iter  = 20) {

  if (is.null(phi_start)) {
    p_x  <- ncol(df_get_X(dat, x_cols))
    phi_start <- c(-log(1 / mean(dat$d_np) - 1),
                   rep(0, p_x),
                   0)
  }

  phi_est <- list(termcd = 99)
  k <- 0
  while (phi_est$termcd > 2 && k < max_iter) {
    init <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
    phi_est <- nleqslv::nleqslv(
      init,
      pi_np.est_simple(dat, df_g4, x_cols, p.use = TRUE)
    )
    k <- k + 1
  }

  if (phi_est$termcd > 2) {
    phi_hat   <- rep(NA_real_, length(phi_start))
    contrib   <- rep(NA_real_, nrow(dat))
    theta_res <- df_sandwich_from_contrib(contrib)
  } else {
    phi_hat <- phi_est$x
    X       <- df_get_X(dat, x_cols)
    l       <- cbind(1, X, dat$y)
    pi_np   <- as.numeric(1 / (1 + exp(-l %*% phi_hat)))
    denom   <- pi_np + dat$pi_p - pi_np * dat$pi_p
    contrib <- dat$y * (dat$d_np + dat$d_p - dat$d_np * dat$d_p) / denom
    theta_res <- df_sandwich_from_contrib(contrib)
  }

  list(phi = c(phi_hat),
       theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci)
}

############################################################
## Efficient estimator Eff (DML2, K-fold)
############################################################

efficient_estimator_dml2 <- function(dat,
                                     x_cols,
                                     phi_start   = NULL,
                                     K           = 2,
                                     max_restart = 10,
                                     progress    = FALSE) {

  N_total <- nrow(dat)
  if (is.null(phi_start)) {
    p_x  <- ncol(df_get_X(dat, x_cols))
    phi_start <- c(-log(1 / mean(dat$d_np) - 1),
                   rep(0, p_x),
                   0)
  }

  folds  <- make_folds(N_total, K)
  dat_cf <- impute_pi_p_crossfit(dat, folds, x_cols,
                                 sigma = NULL,
                                 progress = progress)

  obj_phi <- function(phi) {
    eq_agg <- rep(0, length(phi))
    for (k in seq_along(folds)) {
      idx_k     <- folds[[k]]
      dat_test  <- dat_cf[idx_k, ]
      dat_train <- dat_cf[-idx_k, ]

      ee_fun_k <- estimating_equation_optimal_phi(
        dat_test, dat_train, x_cols = x_cols, eta4_star = NULL
      )
      eq_k   <- ee_fun_k(phi)
      eq_agg <- eq_agg + length(idx_k) / N_total * eq_k
    }
    sum(eq_agg^2)
  }

  attempt <- 0
  res     <- list(convergence = 1)

  if (progress) {
    cat("Solving for phi (Eff, random restarts) ...\n")
    pb_phi <- utils::txtProgressBar(min = 0, max = max_restart,
                                    style = 3)
  }

  while (attempt < max_restart && res$convergence != 0) {
    attempt <- attempt + 1
    init    <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
    res_try <- try(stats::optim(par = init, fn = obj_phi), silent = TRUE)
    if (!inherits(res_try, "try-error")) {
      res <- res_try
    }
    if (progress) {
      utils::setTxtProgressBar(pb_phi, attempt)
    }
  }

  if (progress) {
    close(pb_phi)
  }

  if (res$convergence != 0) {
    phi_hat   <- rep(NA_real_, length(phi_start))
    contrib   <- rep(NA_real_, N_total)
    theta_res <- df_sandwich_from_contrib(contrib)
    return(list(
      phi   = phi_hat,
      theta = theta_res$theta,
      var   = theta_res$var,
      se    = theta_res$se,
      ci    = theta_res$ci
    ))
  }

  phi_hat <- res$par

  if (progress) {
    cat("\nCross-fitting h4*(X) for Eff ...\n")
    pb_h4 <- utils::txtProgressBar(min = 0, max = length(folds),
                                   style = 3)
  }

  h4_star_all <- numeric(N_total)
  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(N_total), idx_test)
    dat_test  <- dat_cf[idx_test, ]
    dat_train <- dat_cf[idx_train, ]

    h4_star_k <- estimate_conditional_expectation_kernlab_theta(
      dat_train, phi_hat, dat_test,
      x_cols = x_cols, sigma = NULL
    )
    h4_star_all[idx_test] <- h4_star_k

    if (progress) {
      utils::setTxtProgressBar(pb_h4, k)
    }
  }

  if (progress) {
    close(pb_h4)
  }

  contrib   <- efficient_theta_contrib(dat_cf, dat_cf, phi_hat,
                                       x_cols = x_cols,
                                       h4_star = h4_star_all)
  theta_res <- df_sandwich_from_contrib(contrib)

  list(phi   = phi_hat,
       theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci)
}

############################################################
## Sub-efficient estimator Eff_S (DML2, K-fold)
############################################################

subefficient_estimator_dml2 <- function(dat,
                                        x_cols,
                                        K        = 2,
                                        progress = FALSE) {

  n     <- nrow(dat)
  folds <- make_folds(n, K)

  X_all <- df_get_X(dat, x_cols)
  idx   <- dat$d_p == 1 & dat$d_np == 0
  X_obs <- X_all[idx, , drop = FALSE]

  if (nrow(X_obs) >= 2L) {
    dist_matrix <- as.matrix(stats::dist(X_obs))
    sigma_mu    <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  } else {
    sigma_mu <- NULL
  }

  if (progress) {
    cat("Cross-fitting mu(x) = E[Y|X] for Eff_S ...\n")
    pb_mu <- utils::txtProgressBar(min = 0, max = length(folds),
                                   style = 3)
  }

  mu_hat_all <- numeric(n)

  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(n), idx_test)

    dat_train <- dat[idx_train, ]
    dat_test  <- dat[idx_test,  ]

    mu_hat_k <- regression_expectation_kernlab(
      dat_train, dat_test, x_cols = x_cols, sigma = sigma_mu
    )
    mu_hat_all[idx_test] <- mu_hat_k

    if (progress) {
      utils::setTxtProgressBar(pb_mu, k)
    }
  }

  if (progress) {
    close(pb_mu)
  }

  contrib   <- subefficient_contrib(dat, mu_hat_all)
  theta_res <- df_sandwich_from_contrib(contrib)

  list(theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci)
}

############################################################
## Efficient parametric estimator Eff_P (working model)
############################################################

efficient_parametric_estimator <- function(dat,
                                           x_cols,
                                           phi_start = NULL,
                                           eta4_star = 0,
                                           max_iter  = 20,
                                           progress  = FALSE) {

  if (is.null(phi_start)) {
    p_x  <- ncol(df_get_X(dat, x_cols))
    phi_start <- c(-log(1 / mean(dat$d_np) - 1),
                   rep(0, p_x),
                   0)
  }

  phi_est_para_fit <- list(termcd = 99)
  k3 <- 0

  if (progress) {
    cat("Solving for phi (Eff_P) ...\n")
    pb_phi <- utils::txtProgressBar(min = 0, max = max_iter,
                                    style = 3)
  }

  while (phi_est_para_fit$termcd > 2 && k3 < max_iter) {
    phi_est_para_fit <- nleqslv::nleqslv(
      phi_start,
      estimating_equation_optimal_phi(dat, dat,
                                      x_cols = x_cols,
                                      eta4_star = eta4_star)
    )
    k3 <- k3 + 1
    if (progress) {
      utils::setTxtProgressBar(pb_phi, k3)
    }
  }

  if (progress) {
    close(pb_phi)
  }

  if (phi_est_para_fit$termcd > 2) {
    phi_est_para <- rep(NA_real_, length(phi_start))
    contrib      <- rep(NA_real_, nrow(dat))
    theta_res    <- df_sandwich_from_contrib(contrib)
  } else {
    phi_est_para <- phi_est_para_fit$x

    ## Working model for E[Y|X]: simple linear regression using all x_cols
    subset_data <- dat[dat$d_np == 1 | dat$d_p == 1, ]
    X_sub       <- df_get_X(subset_data, x_cols)
    df_sub      <- data.frame(y = subset_data$y,
                              X_sub)
    colnames(df_sub)[-1] <- paste0("x", seq_len(ncol(X_sub)))
    formula_str <- paste("y ~", paste(colnames(df_sub)[-1], collapse = " + "))
    lm_fit <- stats::lm(stats::as.formula(formula_str), data = df_sub)

    X_all  <- df_get_X(dat, x_cols)
    df_all <- data.frame(X_all)
    colnames(df_all) <- paste0("x", seq_len(ncol(X_all)))
    predicted_values <- stats::predict(lm_fit, newdata = df_all)

    contrib   <- efficient_theta_contrib(dat, dat, phi_est_para,
                                         x_cols = x_cols,
                                         h4_star = predicted_values)
    theta_res <- df_sandwich_from_contrib(contrib)
  }

  list(phi   = phi_est_para,
       theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci)
}

############################################################
## Public user-facing wrappers: Eff, Eff_S, Eff_P
############################################################

#' Semiparametric efficient estimator Eff (DML2, K-fold)
#'
#' @param dat A data.frame containing at least:
#'   covariates (specified by x_cols), y, d_np, d_p, pi_p, pi_np.
#' @param x_cols Character or integer vector indicating covariate columns.
#'   Default is "x" for backward compatibility.
#' @param K Number of folds for DML2 cross-fitting (default 2).
#' @param phi_start Initial values for the logistic model of pi_NP.
#'   If NULL (default), a suitable length is constructed automatically
#'   based on the dimension of x_cols.
#' @param max_restart Maximum number of random restarts in the optimization
#'   of the phi estimating equation.
#' @param progress Logical; if TRUE, show console progress bars.
#'
#' @return A list with components theta, var, se, ci, phi, info.
#' @export
Eff <- function(dat,
                x_cols     = "x",
                K          = 2,
                phi_start  = NULL,
                max_restart = 10,
                progress   = interactive()) {

  res <- efficient_estimator_dml2(
    dat         = dat,
    x_cols      = x_cols,
    phi_start   = phi_start,
    K           = K,
    max_restart = max_restart,
    progress    = progress
  )

  res$info <- list(
    type        = "Eff",
    x_cols      = x_cols,
    K           = K,
    phi_start   = phi_start,
    max_restart = max_restart,
    progress    = progress
  )

  res
}

#' Sub-efficient estimator Eff_S (Remark 6, DML2, K-fold)
#'
#' @param dat A data.frame containing at least:
#'   covariates, y, d_np, d_p, pi_p, pi_np.
#' @param x_cols Character or integer vector indicating covariate columns.
#' @param K Number of folds for cross-fitting (default 2).
#' @param progress Logical; if TRUE, show console progress bar.
#'
#' @return A list with components theta, var, se, ci, info.
#' @export
Eff_S <- function(dat,
                  x_cols   = "x",
                  K        = 2,
                  progress = interactive()) {

  res <- subefficient_estimator_dml2(dat, x_cols = x_cols,
                                     K = K, progress = progress)

  res$info <- list(
    type     = "Eff_S",
    x_cols   = x_cols,
    K        = K,
    progress = progress
  )

  res
}

#' Parametric efficient estimator Eff_P (working model)
#'
#' @param dat A data.frame containing at least:
#'   covariates, y, d_np, d_p, pi_p, pi_np.
#' @param x_cols Character or integer vector indicating covariate columns.
#' @param phi_start Initial values for the logistic model of pi_NP.
#'   If NULL, the length is set automatically based on x_cols.
#' @param eta4_star Constant (or vector) used as eta4* in the estimating
#'   equation for phi; typically 0 for the working model.
#' @param max_iter Maximum number of attempts (outer iterations) in solving
#'   the phi estimating equation.
#' @param progress Logical; if TRUE, show console progress bar.
#'
#' @return A list with components theta, var, se, ci, phi, info.
#' @export
Eff_P <- function(dat,
                  x_cols   = "x",
                  phi_start = NULL,
                  eta4_star = 0,
                  max_iter  = 20,
                  progress  = interactive()) {

  res <- efficient_parametric_estimator(
    dat       = dat,
    x_cols    = x_cols,
    phi_start = phi_start,
    eta4_star = eta4_star,
    max_iter  = max_iter,
    progress  = progress
  )

  res$info <- list(
    type      = "Eff_P",
    x_cols    = x_cols,
    phi_start = phi_start,
    eta4_star = eta4_star,
    max_iter  = max_iter,
    progress  = progress
  )

  res
}
