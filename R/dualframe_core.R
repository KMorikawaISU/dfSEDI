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
############################################################

############################################################
## Sandwich variance + 95% CI from contributions
############################################################

#' Compute sandwich variance and 95% CI from per-observation contributions
#'
#' @param contrib Numeric vector of length n (influence / pseudo-outcome).
#' @param level Confidence level (default 0.95).
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
## Basic basis functions (for Chang & Kott estimating eq.)
############################################################

df_g2 <- function(dat) {
  cbind(1, dat$x, dat$y)
}
df_g3 <- df_g2
df_g4 <- function(dat) {
  cbind(1, dat$x, dat$x^2)
}

############################################################
## Chang & Kott type estimating equation for pi_NP(phi)
############################################################

pi_np.est_simple <- function(dat, h2, h3, h4, p.use = TRUE) {
  function(phi) {
    x    <- dat$x
    y    <- dat$y
    l    <- cbind(1, x, y)
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

    est_eq <- apply(t(h4(dat)) %*% d_set4, 1, sum)
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
regression_expectation_kernlab <- function(dat, new_x, sigma = NULL) {
  x_obs <- dat$x[dat$d_p == 1 & dat$d_np == 0]
  y_obs <- dat$y[dat$d_p == 1 & dat$d_np == 0]

  if (length(x_obs) < 2L) {
    # Fallback: too few observations, use simple mean
    return(rep(mean(y_obs), length(new_x)))
  }

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(x_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)
  reg_model  <- kernlab::gausspr(x = as.matrix(x_obs),
                                 y = y_obs,
                                 kernel = rbf_kernel)

  as.numeric(predict(reg_model, as.matrix(new_x)))
}

# Nonparametric regression for pi_P (inverse probability regression)
pi_p_estimation_kernlab <- function(dat, new_l, sigma = NULL) {
  x_obs    <- dat$x[dat$d_p == 1]
  y_obs    <- dat$y[dat$d_p == 1]
  l_obs    <- cbind(x_obs, y_obs)
  pi_p_obs <- dat$pi_p[dat$d_p == 1]

  if (length(x_obs) < 2L) {
    # Fallback: use mean pi_P
    return(rep(mean(pi_p_obs), nrow(new_l)))
  }

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(l_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)
  reg_model  <- kernlab::gausspr(x = as.matrix(l_obs),
                                 y = 1 / pi_p_obs,
                                 kernel = rbf_kernel)
  reg_pred <- predict(reg_model, as.matrix(new_l))

  as.numeric(1 / reg_pred)
}

# Conditional expectation for eta4*(L;phi) given X (for efficient phi)
estimate_conditional_expectation_kernlab_phi <- function(dat, phi, new_x, sigma = NULL) {
  x_obs    <- dat$x[dat$d_p == 1 | dat$d_np == 1]
  y_obs    <- dat$y[dat$d_p == 1 | dat$d_np == 1]
  l_obs    <- cbind(1, x_obs, y_obs)
  pi_p_obs <- dat$pi_p[dat$d_p == 1 | dat$d_np == 1]
  pi_np_obs <- 1 / (1 + exp(-l_obs %*% phi))

  Denom_vals <- h4_prob_denom_function(pi_np_obs, pi_p_obs, phi)

  n_obs <- length(x_obs)
  p_dim <- ncol(l_obs)
  Numer_mat <- matrix(NA_real_, nrow = n_obs, ncol = p_dim)
  for (i in seq_len(n_obs)) {
    Numer_mat[i, ] <- eta4_prob_numer_function(
      pi_np_obs[i], pi_p_obs[i], phi = phi, l = l_obs[i, ]
    )
  }

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(x_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)

  eta4_denom_model <- kernlab::gausspr(x = as.matrix(x_obs),
                                       y = Denom_vals,
                                       kernel = rbf_kernel)
  eta4_numer_model <- lapply(
    seq_len(p_dim),
    function(j) kernlab::gausspr(
      x = as.matrix(x_obs),
      y = Numer_mat[, j],
      kernel = rbf_kernel
    )
  )

  eta4_denom_pred <- predict(eta4_denom_model, as.matrix(new_x))
  eta4_numer_pred <- sapply(
    eta4_numer_model,
    function(m) predict(m, as.matrix(new_x))
  )

  sweep(eta4_numer_pred, 1, eta4_denom_pred, FUN = "/")
}

# Conditional expectation for h4*(X;phi) given X (for efficient theta)
estimate_conditional_expectation_kernlab_theta <- function(dat, phi, new_x, sigma = NULL) {
  x_obs    <- dat$x[dat$d_p == 1 | dat$d_np == 1]
  y_obs    <- dat$y[dat$d_p == 1 | dat$d_np == 1]
  l_obs    <- cbind(1, x_obs, y_obs)
  pi_p_obs <- dat$pi_p[dat$d_p == 1 | dat$d_np == 1]
  pi_np_obs <- 1 / (1 + exp(-l_obs %*% phi))

  Denom_vals <- h4_prob_denom_function(pi_np_obs, pi_p_obs, phi)
  Numer_vals <- h4_prob_numer_function(pi_np_obs, pi_p_obs, phi, y_obs)

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(x_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)

  h4_denom_model <- kernlab::gausspr(x = as.matrix(x_obs),
                                     y = Denom_vals,
                                     kernel = rbf_kernel)
  h4_numer_model <- kernlab::gausspr(x = as.matrix(x_obs),
                                     y = Numer_vals,
                                     kernel = rbf_kernel)

  h4_denom_pred <- predict(h4_denom_model, as.matrix(new_x))
  h4_numer_pred <- predict(h4_numer_model, as.matrix(new_x))

  as.numeric(h4_numer_pred / h4_denom_pred)
}

############################################################
## Efficient phi estimating equation (Theorem 2)
############################################################

estimating_equation_optimal_phi <- function(dat1, dat2, eta4_star = NULL) {
  x1    <- dat1$x
  y1    <- dat1$y
  d_p1  <- dat1$d_p
  d_np1 <- dat1$d_np
  pi_p1 <- dat1$pi_p
  l1    <- cbind(1, x1, y1)

  function(phi) {
    if (is.null(eta4_star)) {
      eta4_star_local <- estimate_conditional_expectation_kernlab_phi(
        dat2, phi, x1
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
efficient_theta_contrib <- function(dat1, dat2, phi, h4_star = NULL) {
  x1    <- dat1$x
  y1    <- dat1$y
  d_p1  <- dat1$d_p
  d_np1 <- dat1$d_np
  pi_p1 <- dat1$pi_p
  l1    <- cbind(1, x1, y1)

  pi_np1   <- as.numeric(1 / (1 + exp(-l1 %*% phi)))
  pi_np_p1 <- pi_np1 + pi_p1 - pi_np1 * pi_p1

  if (is.null(h4_star)) {
    h4_star_local <- estimate_conditional_expectation_kernlab_theta(
      dat2, phi, x1
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
optimal_theta <- function(dat1, dat2, phi, h4_star = NULL) {
  contrib <- efficient_theta_contrib(dat1, dat2, phi, h4_star)
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

impute_pi_p_crossfit <- function(dat, folds, sigma = NULL) {
  dat_cf <- dat
  for (k in seq_along(folds)) {
    idx_k     <- folds[[k]]
    dat_fold  <- dat_cf[idx_k, ]
    dat_train <- dat_cf[-idx_k, ]

    l_fold    <- cbind(dat_fold$x, dat_fold$y)
    d_p_fold  <- dat_fold$d_p
    d_np_fold <- dat_fold$d_np

    pi_p_mis_local <- which(d_np_fold == 1 & d_p_fold == 0)
    if (length(pi_p_mis_local) > 0) {
      tilde_pi_p <- pi_p_estimation_kernlab(
        dat_train,
        l_fold[pi_p_mis_local, , drop = FALSE],
        sigma = sigma
      )
      dat_cf$pi_p[idx_k[pi_p_mis_local]] <- tilde_pi_p
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
                           phi_start = NULL,
                           max_iter = 20) {

  if (is.null(phi_start)) {
    phi_start <- c(-log(1 / mean(dat$d_np) - 1), 0, 0)
  }

  phi_est <- list(termcd = 99)
  k <- 0
  while (phi_est$termcd > 2 && k < max_iter) {
    init <- phi_start + stats::runif(3, -0.1, 0.1)
    phi_est <- nleqslv::nleqslv(
      init,
      pi_np.est_simple(dat, df_g2, df_g3, df_g4, p.use = FALSE)
    )
    k <- k + 1
  }

  if (phi_est$termcd > 2) {
    phi_hat   <- rep(NA_real_, 3)
    contrib   <- rep(NA_real_, nrow(dat))
    theta_res <- df_sandwich_from_contrib(contrib)
  } else {
    phi_hat <- phi_est$x
    l       <- cbind(1, dat$x, dat$y)
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
                             phi_start = NULL,
                             max_iter = 20) {

  if (is.null(phi_start)) {
    phi_start <- c(-log(1 / mean(dat$d_np) - 1), 0, 0)
  }

  phi_est <- list(termcd = 99)
  k <- 0
  while (phi_est$termcd > 2 && k < max_iter) {
    init <- phi_start + stats::runif(3, -0.1, 0.1)
    phi_est <- nleqslv::nleqslv(
      init,
      pi_np.est_simple(dat, df_g2, df_g3, df_g4, p.use = TRUE)
    )
    k <- k + 1
  }

  if (phi_est$termcd > 2) {
    phi_hat   <- rep(NA_real_, 3)
    contrib   <- rep(NA_real_, nrow(dat))
    theta_res <- df_sandwich_from_contrib(contrib)
  } else {
    phi_hat <- phi_est$x
    l       <- cbind(1, dat$x, dat$y)
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
                                     phi_start   = c(-2.15, -0.5, -0.75),
                                     K           = 2,
                                     max_restart = 10) {

  N_total <- nrow(dat)

  # K-fold & pi_P cross-fitting
  folds  <- make_folds(N_total, K)
  dat_cf <- impute_pi_p_crossfit(dat, folds)

  # DML2 objective for phi
  obj_phi <- function(phi) {
    eq_agg <- rep(0, length(phi))
    for (k in seq_along(folds)) {
      idx_k     <- folds[[k]]
      dat_test  <- dat_cf[idx_k, ]
      dat_train <- dat_cf[-idx_k, ]

      ee_fun_k <- estimating_equation_optimal_phi(
        dat_test, dat_train, eta4_star = NULL
      )
      eq_k   <- ee_fun_k(phi)
      eq_agg <- eq_agg + length(idx_k) / N_total * eq_k
    }
    sum(eq_agg^2)
  }

  # Optimize phi with random restarts
  attempt <- 0
  res     <- list(convergence = 1)
  while (attempt < max_restart && res$convergence != 0) {
    attempt <- attempt + 1
    init    <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
    res_try <- try(stats::optim(par = init, fn = obj_phi), silent = TRUE)
    if (!inherits(res_try, "try-error")) {
      res <- res_try
    }
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

  # Cross-fitted h4*(X) for all observations
  h4_star_all <- numeric(N_total)
  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(N_total), idx_test)
    dat_test  <- dat_cf[idx_test, ]
    dat_train <- dat_cf[idx_train, ]

    h4_star_k <- estimate_conditional_expectation_kernlab_theta(
      dat_train, phi_hat, dat_test$x
    )
    h4_star_all[idx_test] <- h4_star_k
  }

  # Efficient contributions for all observations
  contrib   <- efficient_theta_contrib(dat_cf, dat_cf, phi_hat,
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

subefficient_estimator_dml2 <- function(dat, K = 2) {

  n     <- nrow(dat)
  folds <- make_folds(n, K)

  x_obs <- dat$x[dat$d_p == 1 & dat$d_np == 0]
  if (length(x_obs) >= 2L) {
    dist_matrix <- as.matrix(stats::dist(x_obs))
    sigma_mu    <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  } else {
    sigma_mu <- NULL
  }

  mu_hat_all <- numeric(n)

  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(n), idx_test)

    dat_train <- dat[idx_train, ]
    dat_test  <- dat[idx_test,  ]

    mu_hat_k <- regression_expectation_kernlab(
      dat_train, dat_test$x, sigma = sigma_mu
    )
    mu_hat_all[idx_test] <- mu_hat_k
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
                                           phi_start = c(-2.15, -0.5, -0.75),
                                           eta4_star = 0) {

  phi_est_para <- rep(NA_real_, 3)

  phi_est_para_fit <- list(termcd = 99)
  k3 <- 0
  while (phi_est_para_fit$termcd > 2 && k3 < 20) {
    phi_est_para_fit <- nleqslv::nleqslv(
      phi_start,
      estimating_equation_optimal_phi(dat, dat, eta4_star = eta4_star)
    )
    k3 <- k3 + 1
  }

  if (phi_est_para_fit$termcd > 2) {
    phi_est_para <- rep(NA_real_, 3)
    contrib      <- rep(NA_real_, nrow(dat))
    theta_res    <- df_sandwich_from_contrib(contrib)
  } else {
    phi_est_para <- phi_est_para_fit$x

    subset_data <- dat[dat$d_np == 1 | dat$d_p == 1, ]
    lm_fit <- stats::lm(y ~ x, data = subset_data)
    predicted_values <- stats::predict(lm_fit, newdata = data.frame(x = dat$x))

    contrib   <- efficient_theta_contrib(dat, dat, phi_est_para,
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
#' Computes the semiparametric efficient estimator (Eff) using DML2 with
#' K-fold cross-fitting. It returns the point estimate, sandwich variance,
#' standard error and 95% confidence interval, together with the estimated
#' phi parameter.
#'
#' @param dat A data.frame containing at least:
#'   x, y, d_np, d_p, pi_p, pi_np.
#' @param K Number of folds for DML2 cross-fitting (default 2).
#' @param phi_start Initial values for the logistic model of pi_NP.
#' @param max_restart Maximum number of random restarts in the optimization
#'   of the phi estimating equation.
#'
#' @return A list with components:
#'   theta, var, se, ci, phi, info.
#' @export
Eff <- function(dat,
                K           = 2,
                phi_start   = c(-2.15, -0.5, -0.75),
                max_restart = 10) {

  res <- efficient_estimator_dml2(
    dat         = dat,
    phi_start   = phi_start,
    K           = K,
    max_restart = max_restart
  )

  res$info <- list(
    type        = "Eff",
    K           = K,
    phi_start   = phi_start,
    max_restart = max_restart
  )

  res
}

#' Sub-efficient estimator Eff_S (Remark 6, DML2, K-fold)
#'
#' Computes the sub-efficient estimator (Eff_S) based on Remark 6
#' using DML2 with K-fold cross-fitting for the regression
#' mu(x) = E[Y | X = x].
#'
#' @param dat A data.frame containing at least:
#'   x, y, d_np, d_p, pi_p, pi_np.
#' @param K Number of folds for cross-fitting (default 2).
#'
#' @return A list with components:
#'   theta, var, se, ci, info.
#' @export
Eff_S <- function(dat,
                  K = 2) {

  res <- subefficient_estimator_dml2(dat, K = K)

  res$info <- list(
    type = "Eff_S",
    K    = K
  )

  res
}

#' Parametric efficient estimator Eff_P (working model)
#'
#' Computes the parametric efficient estimator (Eff_P) under a working
#' parametric model for the efficient score, using the same estimating
#' equation for phi and a linear regression model for E[Y | X].
#'
#' @param dat A data.frame containing at least:
#'   x, y, d_np, d_p, pi_p, pi_np.
#' @param phi_start Initial values for the logistic model of pi_NP.
#' @param eta4_star Constant (or vector) used as eta4* in the estimating
#'   equation for phi; typically 0 for the working model.
#'
#' @return A list with components:
#'   theta, var, se, ci, phi, info.
#' @export
Eff_P <- function(dat,
                  phi_start = c(-2.15, -0.5, -0.75),
                  eta4_star = 0) {

  res <- efficient_parametric_estimator(
    dat       = dat,
    phi_start = phi_start,
    eta4_star = eta4_star
  )

  res$info <- list(
    type      = "Eff_P",
    phi_start = phi_start,
    eta4_star = eta4_star
  )

  res
}
