############################################################
## dualframe_core.R  (multi-dimensional X, DML1 & DML2)
##
## Core estimation functions for semiparametric efficient
## dual-frame data integration with general matrix-valued X.
##
## Public user-facing functions:
##   - Eff   : semiparametric efficient estimator
##             (DML2 by default, DML1 optional)
##   - Eff_S : sub-efficient estimator (Remark 6, DML2)
##   - Eff_P : parametric efficient estimator (working model)
##
## Basic estimators (for simulation / comparison):
##   - df_estimate_P    : probability-only estimator
##   - df_estimate_NP   : NP-only estimator (Chang & Kott type)
##   - df_estimate_NP_P : NP ∪ P estimator (Chang & Kott type)
##
## Assumptions on dat:
##   - either:
##       * dat$X is a numeric matrix (n x p) of covariates
##     or
##       * dat$x is a numeric vector and we treat it as 1D X
##   - dat$y, dat$d_np, dat$d_p, dat$pi_p, dat$pi_np exist and are numeric
##
## External packages (must be available in the session or declared
## in the package NAMESPACE / DESCRIPTION):
##   - nleqslv
##   - kernlab
############################################################

############################################################
## 0. Utility: extract X and sandwich variance
############################################################

# Extract covariate matrix X from dat
# Priority:
#   1) dat$X (matrix / AsIs-matrix / data.frame)
#   2) dat$x (numeric vector -> 1D matrix)
df_get_X <- function(dat) {
  if ("X" %in% names(dat)) {
    X <- dat$X
    if (is.data.frame(X)) {
      X <- as.matrix(X)
    }
    if (!is.matrix(X)) {
      stop("df_get_X(): `X` must be a matrix or convertible to matrix.")
    }
    X <- apply(X, 2, as.numeric)  # ensure numeric
    return(X)
  }
  if ("x" %in% names(dat)) {
    return(matrix(as.numeric(dat$x), ncol = 1))
  }
  stop("df_get_X(): cannot extract X. Provide either a matrix column `X` or a numeric column `x`.")
}

# Sandwich variance + 95% CI from per-observation contributions
#   contrib_i ≈ ψ(Z_i; θ, η)
#   θ_hat = mean(contrib)
#   Var(θ_hat) ≈ Var(contrib)/n
df_sandwich_from_contrib <- function(contrib, level = 0.95) {
  contrib <- as.numeric(contrib)
  contrib <- contrib[is.finite(contrib)]
  n <- length(contrib)

  if (n <= 1L) {
    return(list(
      theta = NA_real_,
      var   = NA_real_,
      se    = NA_real_,
      ci    = c(NA_real_, NA_real_)
    ))
  }

  theta_hat <- mean(contrib)
  s2        <- stats::var(contrib)
  var_hat   <- s2 / n
  se_hat    <- sqrt(var_hat)

  alpha <- 1 - level
  z     <- stats::qnorm(1 - alpha / 2)
  ci    <- c(theta_hat - z * se_hat,
             theta_hat + z * se_hat)

  list(theta = theta_hat,
       var   = var_hat,
       se    = se_hat,
       ci    = ci)
}

############################################################
## 1. Basis functions for NP logistic model (multi-X)
############################################################

g2 <- function(dat) {
  X <- df_get_X(dat)
  cbind(1, X, as.numeric(dat$y))
}
g3 <- g2
g4 <- function(dat) {
  X <- df_get_X(dat)
  cbind(1, X, as.numeric(dat$y))
}

############################################################
## 2. Chang & Kott type estimating eq. for pi_NP(phi)
############################################################

pi_np.est_simple <- function(dat, h2, h3, h4, p.use = TRUE) {
  function(phi) {
    phi <- as.numeric(phi)
    X    <- df_get_X(dat)
    y    <- as.numeric(dat$y)
    l    <- cbind(1, X, y)
    pi_p <- as.numeric(dat$pi_p)

    eta   <- as.numeric(l %*% phi)
    pi_np <- 1 / (1 + exp(-eta))

    d_p  <- as.numeric(dat$d_p)
    d_np <- as.numeric(dat$d_np)

    if (p.use) {
      denom  <- pi_p + pi_np - pi_p * pi_np
      d_set4 <- matrix(1 - (d_p + d_np - d_p * d_np) / denom,
                       nrow = nrow(dat), ncol = 1)
    } else {
      d_set4 <- matrix(1 - d_np / pi_np,
                       nrow = nrow(dat), ncol = 1)
    }

    H4     <- as.matrix(h4(dat))   # n x dim(phi)
    est_eq <- t(H4) %*% d_set4     # dim(phi) x 1
    as.numeric(est_eq)
  }
}

############################################################
## 3. Efficient score helper functions
############################################################

h4_prob_denom_function <- function(pi_np, pi_p, phi) {
  pi_np_p <- pi_p + pi_np - pi_p * pi_np
  (1 - pi_np_p) / (pi_np_p ^ 2)
}

h4_prob_numer_function <- function(pi_np, pi_p, phi, y) {
  pi_np_p <- pi_p + pi_np - pi_p * pi_np
  - y * (1 - pi_np_p) / (pi_np_p ^ 2)
}

eta4_prob_numer_function <- function(pi_np, pi_p, phi, l) {
  pi_np_p <- pi_p + pi_np - pi_p * pi_np
  - as.numeric(pi_np) * (1 - pi_np_p) / (pi_np_p ^ 2) * as.numeric(l)
}

############################################################
## 4. Kernel regression: E(Y|X), pi_P, eta4*, h4* (multi-X)
############################################################

# mu(x) = E[Y|X=x] using P-only data
regression_expectation_kernlab <- function(dat, new_X, sigma = NULL) {
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 & as.numeric(dat$d_np) == 0)
  X_obs <- X_all[idx, , drop = FALSE]
  y_obs <- as.numeric(dat$y[idx])

  if (length(y_obs) < 2L) {
    return(rep(mean(y_obs), nrow(new_X)))
  }

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(X_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)
  reg_model  <- kernlab::gausspr(
    x      = as.matrix(X_obs),
    y      = y_obs,
    kernel = rbf_kernel
  )

  as.numeric(predict(reg_model, as.matrix(new_X)))
}

# Nonparametric regression for pi_P (inverse probability regression)
#   Train on P-sample (d_p == 1) and predict 1/pi_p for new_L
pi_p_estimation_kernlab <- function(dat, new_L, sigma = NULL) {
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1)

  X_obs    <- X_all[idx, , drop = FALSE]
  y_obs    <- as.numeric(dat$y[idx])
  L_obs    <- cbind(X_obs, y_obs)
  pi_p_obs <- as.numeric(dat$pi_p[idx])

  if (length(pi_p_obs) < 2L) {
    return(rep(mean(pi_p_obs), nrow(new_L)))
  }

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(L_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)
  reg_model  <- kernlab::gausspr(
    x      = as.matrix(L_obs),
    y      = 1 / pi_p_obs,
    kernel = rbf_kernel
  )
  reg_pred <- predict(reg_model, as.matrix(new_L))

  as.numeric(1 / reg_pred)
}

# eta4*(L;phi) given X
estimate_conditional_expectation_kernlab_phi <- function(dat, phi, new_X,
                                                         sigma = NULL) {
  phi   <- as.numeric(phi)
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 | as.numeric(dat$d_np) == 1)

  X_obs    <- X_all[idx, , drop = FALSE]
  y_obs    <- as.numeric(dat$y[idx])
  pi_p_obs <- as.numeric(dat$pi_p[idx])

  l_obs  <- cbind(1, X_obs, y_obs)
  l_obs  <- as.matrix(l_obs)
  eta    <- as.numeric(l_obs %*% phi)
  pi_np  <- 1 / (1 + exp(-eta))

  Denom_vals <- h4_prob_denom_function(pi_np, pi_p_obs, phi)

  p_dim     <- length(phi)
  n_obs     <- nrow(l_obs)
  Numer_mat <- matrix(NA_real_, nrow = n_obs, ncol = p_dim)
  for (i in seq_len(n_obs)) {
    Numer_mat[i, ] <- eta4_prob_numer_function(
      pi_np[i], pi_p_obs[i], phi = phi, l = l_obs[i, ]
    )
  }

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(X_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)

  eta4_denom_model <- kernlab::gausspr(
    x      = as.matrix(X_obs),
    y      = as.numeric(Denom_vals),
    kernel = rbf_kernel
  )
  eta4_numer_model <- lapply(
    seq_len(p_dim),
    function(j) kernlab::gausspr(
      x      = as.matrix(X_obs),
      y      = as.numeric(Numer_mat[, j]),
      kernel = rbf_kernel
    )
  )

  eta4_denom_pred <- predict(eta4_denom_model, as.matrix(new_X))
  eta4_numer_pred <- sapply(
    eta4_numer_model,
    function(m) predict(m, as.matrix(new_X))
  )

  sweep(eta4_numer_pred, 1, eta4_denom_pred, "/")
}

# h4*(X;phi) given X
estimate_conditional_expectation_kernlab_theta <- function(dat, phi, new_X,
                                                           sigma = NULL) {
  phi   <- as.numeric(phi)
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 | as.numeric(dat$d_np) == 1)

  X_obs    <- X_all[idx, , drop = FALSE]
  y_obs    <- as.numeric(dat$y[idx])
  pi_p_obs <- as.numeric(dat$pi_p[idx])

  l_obs  <- cbind(1, X_obs, y_obs)
  l_obs  <- as.matrix(l_obs)
  eta    <- as.numeric(l_obs %*% phi)
  pi_np  <- 1 / (1 + exp(-eta))

  Denom_vals <- h4_prob_denom_function(pi_np, pi_p_obs, phi)
  Numer_vals <- h4_prob_numer_function(pi_np, pi_p_obs, phi, y_obs)

  if (is.null(sigma)) {
    dist_matrix <- as.matrix(stats::dist(X_obs))
    sigma <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  }

  rbf_kernel <- kernlab::rbfdot(sigma = sigma)

  h4_denom_model <- kernlab::gausspr(
    x      = as.matrix(X_obs),
    y      = as.numeric(Denom_vals),
    kernel = rbf_kernel
  )
  h4_numer_model <- kernlab::gausspr(
    x      = as.matrix(X_obs),
    y      = as.numeric(Numer_vals),
    kernel = rbf_kernel
  )

  h4_denom_pred <- predict(h4_denom_model, as.matrix(new_X))
  h4_numer_pred <- predict(h4_numer_model, as.matrix(new_X))

  as.numeric(h4_numer_pred / h4_denom_pred)
}

############################################################
## 5. Efficient phi estimating equation (multi-X)
############################################################

# Returns a function of phi:
#   ee(phi) = mean_i s_phi(z_i; phi, eta4_star^{(-k)})
# where dat1 plays the role of "test fold" and dat2 is "training"
estimating_equation_optimal_phi <- function(dat1,
                                            dat2,
                                            eta4_star = NULL) {
  X1    <- df_get_X(dat1)
  y1    <- as.numeric(dat1$y)
  d_p1  <- as.numeric(dat1$d_p)
  d_np1 <- as.numeric(dat1$d_np)
  pi_p1 <- as.numeric(dat1$pi_p)
  l1    <- as.matrix(cbind(1, X1, y1))

  function(phi) {
    phi <- as.numeric(phi)
    if (is.null(eta4_star)) {
      eta4_star_local <- estimate_conditional_expectation_kernlab_phi(
        dat2, phi, new_X = X1, sigma = NULL
      )
    } else {
      eta4_star_local <- eta4_star
    }

    eta      <- as.numeric(l1 %*% phi)
    pi_np1   <- 1 / (1 + exp(-eta))
    pi_np_p1 <- pi_np1 + pi_p1 - pi_np1 * pi_p1

    term1 <- l1 * pi_np1 * (1 - pi_np1) / pi_np_p1 *
      (d_np1 / pi_np1 * (pi_p1 - d_p1) -
         d_p1 / (1 - pi_np1) * (1 - d_np1 / pi_np1))

    term2_coef <- 1 -
      (d_p1 * d_np1) / (pi_p1 * pi_np1) -
      1 / pi_np_p1 *
      (d_np1 * (1 - d_p1 / pi_p1) +
         d_p1 * (1 - d_np1 / pi_np1))

    est_eq <- term1 + term2_coef * eta4_star_local

    colMeans(est_eq)
  }
}

############################################################
## 6. Efficient theta contributions (multi-X)
############################################################

efficient_theta_contrib <- function(dat1,
                                    dat2,  # kept for interface; not used if h4_star given
                                    phi,
                                    h4_star = NULL) {
  phi   <- as.numeric(phi)
  X1    <- df_get_X(dat1)
  y1    <- as.numeric(dat1$y)
  d_p1  <- as.numeric(dat1$d_p)
  d_np1 <- as.numeric(dat1$d_np)
  pi_p1 <- as.numeric(dat1$pi_p)
  l1    <- as.matrix(cbind(1, X1, y1))

  eta      <- as.numeric(l1 %*% phi)
  pi_np1   <- 1 / (1 + exp(-eta))
  pi_np_p1 <- pi_np1 + pi_p1 - pi_np1 * pi_p1

  if (is.null(h4_star)) {
    # Non-cross-fitted version (rarely used now)
    h4_star_local <- estimate_conditional_expectation_kernlab_theta(
      dat2, phi, new_X = X1, sigma = NULL
    )
  } else {
    h4_star_local <- as.numeric(h4_star)
  }

  term1 <- y1 * (d_p1 / pi_p1 +
                   (1 - d_p1 / pi_p1) / pi_np_p1 *
                   (d_np1 - d_p1 * (d_np1 - pi_np1)))

  term2_coef <- 1 -
    (d_p1 * d_np1) / (pi_p1 * pi_np1) -
    1 / pi_np_p1 *
    (d_np1 * (1 - d_p1 / pi_p1) +
       d_p1 * (1 - d_np1 / pi_np1))

  as.numeric(term1 - term2_coef * h4_star_local)
}

optimal_theta <- function(dat1, dat2, phi, h4_star = NULL) {
  contrib <- efficient_theta_contrib(dat1, dat2, phi, h4_star)
  mean(contrib)
}

############################################################
## 7. Sub-efficient contributions (Eff_S)
############################################################

subefficient_contrib <- function(dat, mu_hat) {
  d_np <- as.numeric(dat$d_np)
  d_p  <- as.numeric(dat$d_p)
  y    <- as.numeric(dat$y)
  pi_p <- as.numeric(dat$pi_p)
  mu   <- as.numeric(mu_hat)

  d_np * y +
    (1 - d_np) * (d_p / pi_p * y + (1 - d_p / pi_p) * mu)
}

############################################################
## 8. K-fold split & pi_P cross-fitting (multi-X)
############################################################

make_folds <- function(n, K) {
  fold_id <- sample(rep(1:K, length.out = n))
  split(seq_len(n), fold_id)
}

# Cross-fitting of pi_p:
#   For each fold k:
#     - training set T_k = complement of S_k
#     - for NP-only units in S_k, predict pi_p from T_k via kernel regression
#   This is used both in DML1 and DML2 so that pi_p is treated
#   as a nuisance estimated on the same training sets as eta (h4*, eta4*).
impute_pi_p_crossfit <- function(dat,
                                 folds,
                                 sigma    = NULL,
                                 progress = FALSE) {
  dat_cf <- dat

  if (progress) {
    cat("Step 1/3: cross-fitting pi_P ...\n")
    pb <- utils::txtProgressBar(min = 0, max = length(folds), style = 3)
  }

  for (k in seq_along(folds)) {
    idx_k     <- folds[[k]]
    dat_fold  <- dat_cf[idx_k, ]
    dat_train <- dat_cf[-idx_k, ]

    X_f   <- df_get_X(dat_fold)
    y_f   <- as.numeric(dat_fold$y)
    d_p_f  <- as.numeric(dat_fold$d_p)
    d_np_f <- as.numeric(dat_fold$d_np)

    L_f  <- cbind(X_f, y_f)
    L_f  <- as.matrix(L_f)

    # NP-only units in this fold: need pi_p_hat from training T_k
    pi_p_mis <- which(d_np_f == 1 & d_p_f == 0)
    if (length(pi_p_mis) > 0) {
      tilde_pi <- pi_p_estimation_kernlab(
        dat_train,
        new_L = L_f[pi_p_mis, , drop = FALSE],
        sigma = sigma
      )
      dat_cf$pi_p[idx_k[pi_p_mis]] <- tilde_pi
    }

    if (progress) {
      utils::setTxtProgressBar(pb, k)
    }
  }

  if (progress) {
    close(pb)
    cat("\n")
  }

  dat_cf
}

############################################################
## 9. Basic estimators: P, NP, NP+P
############################################################

df_estimate_P <- function(dat) {
  contrib   <- as.numeric(dat$d_p) * as.numeric(dat$y) / as.numeric(dat$pi_p)
  theta_res <- df_sandwich_from_contrib(contrib)
  list(theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci)
}

df_estimate_NP <- function(dat,
                           phi_start = NULL,
                           max_iter  = 20) {
  X <- df_get_X(dat)
  p <- ncol(X)

  if (is.null(phi_start)) {
    phi_start <- c(-log(1 / mean(dat$d_np) - 1),
                   rep(0, p),
                   0)
  }

  phi_est <- list(termcd = 99)
  k <- 0
  while (phi_est$termcd > 2 && k < max_iter) {
    init <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
    phi_est <- nleqslv::nleqslv(
      init,
      pi_np.est_simple(dat, g2, g3, g4, p.use = FALSE)
    )
    k <- k + 1
  }

  if (phi_est$termcd > 2) {
    phi_hat   <- rep(NA_real_, length(phi_start))
    contrib   <- rep(NA_real_, nrow(dat))
    theta_res <- df_sandwich_from_contrib(contrib)
  } else {
    phi_hat <- as.numeric(phi_est$x)
    l       <- as.matrix(cbind(1, X, as.numeric(dat$y)))
    eta     <- as.numeric(l %*% phi_hat)
    pi_np   <- 1 / (1 + exp(-eta))

    contrib   <- as.numeric(dat$d_np) * as.numeric(dat$y) / pi_np
    theta_res <- df_sandwich_from_contrib(contrib)
  }

  list(phi   = phi_hat,
       theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci)
}

df_estimate_NP_P <- function(dat,
                             phi_start = NULL,
                             max_iter  = 20) {
  X <- df_get_X(dat)
  p <- ncol(X)

  if (is.null(phi_start)) {
    phi_start <- c(-log(1 / mean(dat$d_np) - 1),
                   rep(0, p),
                   0)
  }

  phi_est <- list(termcd = 99)
  k <- 0
  while (phi_est$termcd > 2 && k < max_iter) {
    init <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
    phi_est <- nleqslv::nleqslv(
      init,
      pi_np.est_simple(dat, g2, g3, g4, p.use = TRUE)
    )
    k <- k + 1
  }

  if (phi_est$termcd > 2) {
    phi_hat   <- rep(NA_real_, length(phi_start))
    contrib   <- rep(NA_real_, nrow(dat))
    theta_res <- df_sandwich_from_contrib(contrib)
  } else {
    phi_hat <- as.numeric(phi_est$x)
    l       <- as.matrix(cbind(1, X, as.numeric(dat$y)))
    eta     <- as.numeric(l %*% phi_hat)
    pi_np   <- 1 / (1 + exp(-eta))

    pi_p   <- as.numeric(dat$pi_p)
    d_np   <- as.numeric(dat$d_np)
    d_p    <- as.numeric(dat$d_p)
    denom  <- pi_np + pi_p - pi_np * pi_p
    contrib <- as.numeric(dat$y) *
      (d_np + d_p - d_np * d_p) / denom
    theta_res <- df_sandwich_from_contrib(contrib)
  }

  list(phi   = phi_hat,
       theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci)
}

############################################################
## 10. Efficient estimator Eff (DML2, multi-X)
############################################################

efficient_estimator_dml2 <- function(dat,
                                     phi_start   = NULL,
                                     K           = 2,
                                     max_restart = 10,
                                     progress    = FALSE) {
  N_total <- nrow(dat)
  X_all   <- df_get_X(dat)
  p_x     <- ncol(X_all)

  if (is.null(phi_start)) {
    phi_start <- c(
      -log(1 / mean(dat$d_np) - 1),
      rep(0, p_x),
      0
    )
  }

  # Folds for DML2
  folds  <- make_folds(N_total, K)

  # Cross-fit pi_p on the same folds as eta
  dat_cf <- impute_pi_p_crossfit(dat, folds,
                                 sigma    = NULL,
                                 progress = progress)

  # Objective for phi (DML2): sum_k w_k * ee_k(phi)
  obj_phi <- function(phi) {
    phi <- as.numeric(phi)
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
    sum(eq_agg ^ 2)
  }

  if (progress) {
    cat("Step 2/3: solving phi (Eff, DML2, Nelder-Mead) ...\n")
  }

  attempt <- 0
  res     <- list(convergence = 1)
  while (attempt < max_restart && res$convergence != 0) {
    attempt <- attempt + 1
    if (progress) {
      cat(sprintf("  attempt %d / %d ...\n", attempt, max_restart))
      flush.console()
    }

    init <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
    res_try <- try(
      stats::optim(par = init,
                   fn  = obj_phi,
                   method = "Nelder-Mead"),
      silent = TRUE
    )
    if (!inherits(res_try, "try-error")) {
      res <- res_try
      if (res$convergence == 0) break
    }
  }

  if (res$convergence != 0 && progress) {
    cat("  NOTE: optim did not fully converge; using last iterate as phi_hat.\n")
  }

  phi_hat <- as.numeric(res$par)

  # Cross-fitting h4*(X) for theta (Eff)
  if (progress) {
    cat("Step 3/3: cross-fitting h4*(X) for Eff ...\n")
    pb <- utils::txtProgressBar(min = 0, max = length(folds), style = 3)
  }

  h4_star_all <- numeric(N_total)
  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(N_total), idx_test)

    dat_test  <- dat_cf[idx_test,  ]
    dat_train <- dat_cf[idx_train, ]
    X_test    <- df_get_X(dat_test)

    h4_star_k <- estimate_conditional_expectation_kernlab_theta(
      dat_train, phi_hat, new_X = X_test, sigma = NULL
    )
    h4_star_all[idx_test] <- h4_star_k

    if (progress) {
      utils::setTxtProgressBar(pb, k)
    }
  }

  if (progress) {
    close(pb)
    cat("\n")
  }

  contrib   <- efficient_theta_contrib(dat_cf, dat_cf, phi_hat,
                                       h4_star = h4_star_all)
  theta_res <- df_sandwich_from_contrib(contrib)

  list(phi   = phi_hat,
       theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci,
       info  = list(
         type        = "Eff",
         dml_type    = "DML2",
         K           = K,
         phi_start   = phi_start,
         max_restart = max_restart,
         progress    = progress
       ))
}

############################################################
## 11. Efficient estimator Eff (DML1, multi-X)
############################################################

efficient_estimator_dml1 <- function(dat,
                                     phi_start   = NULL,
                                     K           = 2,
                                     max_restart = 10,
                                     progress    = FALSE) {
  N_total <- nrow(dat)
  X_all   <- df_get_X(dat)
  p_x     <- ncol(X_all)

  if (is.null(phi_start)) {
    phi_start <- c(
      -log(1 / mean(dat$d_np) - 1),
      rep(0, p_x),
      0
    )
  }

  # Folds for DML1
  folds <- make_folds(N_total, K)

  # Cross-fit pi_p on the same folds as eta (DML1 version)
  dat_cf <- impute_pi_p_crossfit(dat, folds,
                                 sigma    = NULL,
                                 progress = progress)

  phi_mat      <- matrix(NA_real_, nrow = K, ncol = length(phi_start))
  theta_k      <- rep(NA_real_, K)
  contrib_all  <- rep(NA_real_, N_total)

  if (progress) {
    cat("DML1: per-fold estimation for phi and theta ...\n")
  }

  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(N_total), idx_test)

    dat_test  <- dat_cf[idx_test,  ]
    dat_train <- dat_cf[idx_train, ]

    # 1) Solve phi^{(k)} using only test fold S_k with nuisances from T_k
    ee_fun_k <- estimating_equation_optimal_phi(
      dat_test, dat_train, eta4_star = NULL
    )

    attempt <- 0
    phi_sol <- list(termcd = 99)
    while (attempt < max_restart && phi_sol$termcd > 2) {
      attempt  <- attempt + 1
      phi_init <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
      phi_try  <- try(
        nleqslv::nleqslv(phi_init, ee_fun_k),
        silent = TRUE
      )
      if (!inherits(phi_try, "try-error")) {
        phi_sol <- phi_try
      }
    }

    if (phi_sol$termcd > 2) {
      if (progress) {
        cat(sprintf("  fold %d: phi did not converge.\n", k))
      }
      next
    }

    phi_k <- as.numeric(phi_sol$x)
    phi_mat[k, ] <- phi_k

    # 2) Compute theta^{(k)} on S_k using nuisances from T_k
    X_test <- df_get_X(dat_test)
    h4_star_k <- estimate_conditional_expectation_kernlab_theta(
      dat_train, phi_k, new_X = X_test, sigma = NULL
    )

    contrib_k <- efficient_theta_contrib(
      dat_test, dat_train, phi_k, h4_star = h4_star_k
    )
    theta_k[k] <- mean(contrib_k)

    # store per-observation contributions for sandwich
    contrib_all[idx_test] <- contrib_k

    if (progress) {
      cat(sprintf("  fold %d/%d done.\n", k, K))
      flush.console()
    }
  }

  # Aggregate phi over successful folds
  valid_phi <- apply(!is.na(phi_mat), 1, all)
  if (!any(valid_phi)) {
    phi_hat <- rep(NA_real_, length(phi_start))
  } else {
    phi_hat <- colMeans(phi_mat[valid_phi, , drop = FALSE])
  }

  # Aggregate theta as mean of all contributions (consistent with fold-wise mean)
  theta_res <- df_sandwich_from_contrib(contrib_all)

  list(phi   = phi_hat,
       theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci,
       info  = list(
         type        = "Eff",
         dml_type    = "DML1",
         K           = K,
         phi_start   = phi_start,
         max_restart = max_restart,
         progress    = progress
       ))
}

############################################################
## 12. Sub-efficient estimator Eff_S (multi-X, DML2-style)
############################################################

subefficient_estimator_dml2 <- function(dat,
                                        K        = 2,
                                        progress = FALSE) {
  n     <- nrow(dat)
  folds <- make_folds(n, K)

  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 & as.numeric(dat$d_np) == 0)
  if (length(idx) >= 2L) {
    dist_matrix <- as.matrix(stats::dist(X_all[idx, , drop = FALSE]))
    sigma_mu    <- 1 / stats::median(dist_matrix[upper.tri(dist_matrix)])
  } else {
    sigma_mu <- NULL
  }

  if (progress) {
    cat("Cross-fitting mu(X) for Eff_S ...\n")
    pb <- utils::txtProgressBar(min = 0, max = length(folds), style = 3)
  }

  mu_hat_all <- numeric(n)
  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(n), idx_test)

    dat_train <- dat[idx_train, ]
    dat_test  <- dat[idx_test,  ]
    X_test    <- df_get_X(dat_test)

    mu_hat_k <- regression_expectation_kernlab(
      dat_train, new_X = X_test, sigma = sigma_mu
    )
    mu_hat_all[idx_test] <- mu_hat_k

    if (progress) {
      utils::setTxtProgressBar(pb, k)
    }
  }

  if (progress) {
    close(pb)
    cat("\n")
  }

  contrib   <- subefficient_contrib(dat, mu_hat_all)
  theta_res <- df_sandwich_from_contrib(contrib)

  list(theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci,
       info  = list(
         type     = "Eff_S",
         K        = K,
         progress = progress
       ))
}

############################################################
## 13. Parametric efficient estimator Eff_P (multi-X)
############################################################

efficient_parametric_estimator <- function(dat,
                                           phi_start = NULL,
                                           eta4_star = 0,
                                           max_iter  = 20,
                                           progress  = FALSE) {
  X_all <- df_get_X(dat)
  p_x   <- ncol(X_all)
  n     <- nrow(dat)

  if (is.null(phi_start)) {
    phi_start <- c(
      -log(1 / mean(dat$d_np) - 1),
      rep(0, p_x),
      0
    )
  }

  phi_est <- list(termcd = 99)
  k <- 0
  while (phi_est$termcd > 2 && k < max_iter) {
    if (progress) {
      cat(sprintf("Solving phi for Eff_P, attempt %d / %d ...\n",
                  k + 1, max_iter))
    }
    phi_est <- nleqslv::nleqslv(
      phi_start,
      estimating_equation_optimal_phi(dat, dat, eta4_star = eta4_star)
    )
    k <- k + 1
  }

  if (phi_est$termcd > 2) {
    if (progress) {
      cat("  NOTE: phi did not fully converge in Eff_P.\n")
    }
    phi_hat   <- rep(NA_real_, length(phi_start))
    contrib   <- rep(NA_real_, n)
    theta_res <- df_sandwich_from_contrib(contrib)
  } else {
    phi_hat <- as.numeric(phi_est$x)

    subset_idx  <- which(as.numeric(dat$d_np) == 1 | as.numeric(dat$d_p) == 1)
    subset_data <- dat[subset_idx, ]
    X_sub       <- df_get_X(subset_data)
    df_sub      <- data.frame(y = as.numeric(subset_data$y), X_sub)
    colnames(df_sub)[-1] <- paste0("x", seq_len(ncol(X_sub)))

    lm_fit <- stats::lm(y ~ ., data = df_sub)

    X_all_full <- df_get_X(dat)
    df_full    <- data.frame(X_all_full)
    colnames(df_full) <- paste0("x", seq_len(ncol(X_all_full)))
    predicted_values <- stats::predict(lm_fit, newdata = df_full)

    contrib   <- efficient_theta_contrib(
      dat, dat, phi_hat,
      h4_star = predicted_values
    )
    theta_res <- df_sandwich_from_contrib(contrib)
  }

  list(phi   = phi_hat,
       theta = theta_res$theta,
       var   = theta_res$var,
       se    = theta_res$se,
       ci    = theta_res$ci,
       info  = list(
         type      = "Eff_P",
         phi_start = phi_start,
         eta4_star = eta4_star,
         max_iter  = max_iter,
         progress  = progress
       ))
}

############################################################
## 14. Public wrappers: Eff, Eff_S, Eff_P
############################################################
## Eff:
##   type = 1  or dml_type = 1 / "DML1" -> DML1
##   type = 2  or dml_type = 2 / "DML2" -> DML2 (default)
############################################################

Eff <- function(dat,
                K           = 2,
                phi_start   = NULL,
                max_restart = 10,
                type        = NULL,   # <- added for backward compatibility
                dml_type    = 2,
                progress    = interactive()) {

  # If `type` is provided, use it as alias for dml_type
  if (!is.null(type)) {
    dml_type <- type
  }

  # allow numeric 1/2 or character "DML1"/"DML2"
  if (is.numeric(dml_type)) {
    dml_type <- if (dml_type == 1) "DML1" else "DML2"
  }
  dml_type <- toupper(as.character(dml_type))
  if (!dml_type %in% c("DML1", "DML2")) {
    stop("Eff(): dml_type must be 1, 2, 'DML1', or 'DML2'.")
  }

  if (dml_type == "DML2") {
    efficient_estimator_dml2(
      dat         = dat,
      phi_start   = phi_start,
      K           = K,
      max_restart = max_restart,
      progress    = progress
    )
  } else {
    efficient_estimator_dml1(
      dat         = dat,
      phi_start   = phi_start,
      K           = K,
      max_restart = max_restart,
      progress    = progress
    )
  }
}

Eff_S <- function(dat,
                  K        = 2,
                  progress = interactive()) {
  subefficient_estimator_dml2(
    dat      = dat,
    K        = K,
    progress = progress
  )
}

Eff_P <- function(dat,
                  phi_start = NULL,
                  eta4_star = 0,
                  max_iter  = 20,
                  progress  = interactive()) {
  efficient_parametric_estimator(
    dat       = dat,
    phi_start = phi_start,
    eta4_star = eta4_star,
    max_iter  = max_iter,
    progress  = progress
  )
}
