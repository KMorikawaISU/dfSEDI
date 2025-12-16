############################################################
## dualframe_core.R  (multi-dimensional X, DML1 & DML2, JOINT sandwich)
##
## Core estimation functions for semiparametric efficient
## dual-frame data integration with general matrix-valued X.
##
## Public user-facing functions:
##   - Eff   : semiparametric efficient estimator (DML2 default, DML1 optional)
##   - Eff_S : sub-efficient estimator (Remark 6, DML2-style)
##   - Eff_P : parametric efficient estimator (working model)
##
## Basic estimators (for simulation / comparison):
##   - df_estimate_P    : probability-only estimator
##   - df_estimate_NP   : NP-only estimator (Chang & Kott type)
##   - df_estimate_NP_P : NP ∪ P estimator (Chang & Kott type)
##
## Key update in this version:
##   - (phi, theta) are handled JOINTLY in the sandwich variance
##     for BOTH DML1 and DML2 (phi_var/phi_se/phi_ci are returned).
##
## Assumptions on dat:
##   - either:
##       * dat$X is a matrix / AsIs-matrix / data.frame of covariates
##         (data.frame may include factor/character/logical; they are encoded to numeric codes)
##     or
##       * dat$x is a numeric vector and we treat it as 1D X
##   - dat$y, dat$d_np, dat$d_p, dat$pi_p exist and are numeric
##   - dat$pi_p may contain NA for NP-only units; it will be imputed via kernel regression
##
## External packages used via namespaces:
##   - nleqslv::nleqslv
##   - kernlab::gausspr, kernlab::rbfdot
##   - stats::*, utils::*
############################################################

############################################################
## 0. Utilities: extract X, robust kernlab prediction, mixed-type KRR, numerics
############################################################

# Extract covariate matrix X from dat.
# Priority:
#   1) dat$X (matrix / AsIs-matrix / data.frame)
#   2) dat$x (numeric vector -> 1D matrix)
#
# NOTE:
#   - If dat$X is a data.frame and contains factor/character/logical columns,
#     we encode them into integer codes (1..K) so downstream computations remain numeric.
#   - We attach an attribute "dfSEDI_categorical_cols" for columns that were
#     originally factor/character/logical in a data.frame (best-effort).
df_get_X <- function(dat) {
  if ("X" %in% names(dat)) {
    X <- dat$X
    # If X is an AsIs matrix-column created by I(X), it is typically still a matrix.
    # We do not need special handling, but keep this comment for clarity.

    # If X is a data.frame, convert column-by-column robustly
    if (is.data.frame(X)) {
      n <- nrow(X)
      p <- ncol(X)
      X_num <- matrix(NA_real_, nrow = n, ncol = p)
      colnames(X_num) <- colnames(X)

      cat_cols <- integer(0)

      for (j in seq_len(p)) {
        col <- X[[j]]

        if (is.factor(col)) {
          X_num[, j] <- as.numeric(col)   # factor -> integer codes as numeric
          cat_cols <- c(cat_cols, j)
        } else if (is.character(col)) {
          X_num[, j] <- as.numeric(factor(col))
          cat_cols <- c(cat_cols, j)
        } else if (is.logical(col)) {
          X_num[, j] <- as.numeric(col)
          cat_cols <- c(cat_cols, j)
        } else {
          # numeric/integer/etc.
          X_num[, j] <- suppressWarnings(as.numeric(col))
        }
      }

      attr(X_num, "dfSEDI_categorical_cols") <- unique(cat_cols)
      return(X_num)
    }

    # If X is a matrix (possibly non-numeric), coerce column-by-column
    if (is.matrix(X)) {
      n <- nrow(X)
      p <- ncol(X)
      X_num <- matrix(NA_real_, nrow = n, ncol = p)
      dimnames(X_num) <- dimnames(X)

      # If the matrix is already numeric/integer, this is fast and safe.
      if (is.numeric(X) || is.integer(X)) {
        storage.mode(X) <- "numeric"
        X_num <- X
        attr(X_num, "dfSEDI_categorical_cols") <- integer(0)
        return(X_num)
      }

      # Otherwise, attempt robust conversion.
      cat_cols <- integer(0)
      for (j in seq_len(p)) {
        col <- X[, j]
        if (is.factor(col)) {
          X_num[, j] <- as.numeric(col)
          cat_cols <- c(cat_cols, j)
        } else if (is.character(col)) {
          tmp <- suppressWarnings(as.numeric(col))
          # If non-numeric strings, encode as factor codes
          if (all(is.na(tmp) & !is.na(col))) {
            X_num[, j] <- as.numeric(factor(col))
            cat_cols <- c(cat_cols, j)
          } else {
            X_num[, j] <- tmp
          }
        } else if (is.logical(col)) {
          X_num[, j] <- as.numeric(col)
          cat_cols <- c(cat_cols, j)
        } else {
          X_num[, j] <- suppressWarnings(as.numeric(col))
        }
      }

      attr(X_num, "dfSEDI_categorical_cols") <- unique(cat_cols)
      return(X_num)
    }

    stop("df_get_X(): `X` must be a matrix or data.frame (or convertible).")
  }

  if ("x" %in% names(dat)) {
    X_num <- matrix(as.numeric(dat$x), ncol = 1)
    attr(X_num, "dfSEDI_categorical_cols") <- integer(0)
    return(X_num)
  }

  stop("df_get_X(): cannot extract X. Provide either a matrix/data.frame column `X` or a numeric column `x`.")
}

# Robust prediction wrapper for kernlab::gausspr objects.
# This avoids common failures where stats::predict() cannot dispatch on gausspr objects.
df_kernlab_predict <- function(model, newdata) {
  newdata <- as.matrix(newdata)

  # Try kernlab's predict (S4) first via namespace lookup
  pred <- tryCatch({
    if (requireNamespace("kernlab", quietly = TRUE) &&
        exists("predict", envir = asNamespace("kernlab"), inherits = FALSE)) {
      get("predict", envir = asNamespace("kernlab"))(model, newdata)
    } else {
      stats::predict(model, newdata)
    }
  }, error = function(e1) {
    # Fallback: base/stats predict
    tryCatch({
      stats::predict(model, newdata)
    }, error = function(e2) {
      stop(
        "dfSEDI: prediction failed for a kernlab::gausspr model.\n",
        "This usually happens when the predict method is not available in the current session.\n",
        "Please try running `library(kernlab)` once and re-run.\n",
        "Original error: ", conditionMessage(e2)
      )
    })
  })

  as.numeric(pred)
}

############################################################
## Mixed-type KRR helpers (continuous + categorical auto handling)
############################################################

# Heuristic detection of "categorical-like" columns among numeric columns:
# - unique levels <= max_levels
# - unique fraction <= max_unique_frac
df_detect_categorical_cols_heuristic <- function(X, max_levels = 20L, max_unique_frac = 0.2) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (p == 0L) return(integer(0))
  if (n <= 1L) return(integer(0))

  out <- logical(p)

  for (j in seq_len(p)) {
    x <- X[, j]
    x <- x[!is.na(x)]
    if (length(x) == 0L) {
      out[j] <- TRUE
      next
    }
    u <- length(unique(x))
    frac <- u / max(1, length(x))
    out[j] <- (u <= max_levels) && (frac <= max_unique_frac)
  }

  which(out)
}

# Infer categorical columns:
#  1) columns tagged by df_get_X() as originally factor/character/logical in data.frame
#  2) plus heuristic detection for numeric discrete columns (e.g., 0/1, small integers)
df_infer_categorical_cols <- function(X,
                                      max_levels = 20L,
                                      max_unique_frac = 0.2) {
  X <- as.matrix(X)
  cat_attr <- attr(X, "dfSEDI_categorical_cols")
  if (is.null(cat_attr)) cat_attr <- integer(0)
  cat_attr <- as.integer(cat_attr)

  cat_heur <- df_detect_categorical_cols_heuristic(
    X, max_levels = max_levels, max_unique_frac = max_unique_frac
  )

  sort(unique(c(cat_attr, cat_heur)))
}

# Estimate RBF kernel parameter sigma (kernlab's rbfdot uses exp(-sigma * ||x-x'||^2)).
# Uses median pairwise distance on (optionally) a random subsample to avoid O(n^2) blow-up.
df_estimate_rbf_sigma <- function(X_cont, max_n = 2000L) {
  X_cont <- as.matrix(X_cont)
  n <- nrow(X_cont)

  if (n < 2L) return(NULL)

  max_n <- as.integer(max_n)
  if (is.na(max_n) || max_n <= 1L) max_n <- n
  if (n > max_n) {
    idx <- sample.int(n, max_n)
    X_use <- X_cont[idx, , drop = FALSE]
  } else {
    X_use <- X_cont
  }

  d <- stats::dist(X_use)
  dvec <- as.numeric(d)
  dvec <- dvec[is.finite(dvec) & dvec > 0]

  if (length(dvec) == 0L) return(1)

  med <- stats::median(dvec)
  if (!is.finite(med) || med <= 0) return(1)

  1 / med
}

# Create a group key for categorical columns (row-wise).
df_make_cat_key <- function(X_cat) {
  X_cat <- as.matrix(X_cat)
  n <- nrow(X_cat)
  if (n == 0L) return(character(0))
  if (ncol(X_cat) == 0L) return(rep("ALL", n))

  # Convert to data.frame for do.call(paste, ...)
  dfc <- as.data.frame(X_cat, stringsAsFactors = FALSE)
  # Use a separator unlikely to appear in numeric print
  do.call(paste, c(dfc, sep = "\r"))
}

# Mixed-type KRR prediction:
# - Identify categorical cols (auto by default)
# - Split by categorical key (delta kernel)
# - Within each key, run RBF KRR on continuous columns using kernlab::gausspr
# - If a test key is unseen, fallback to "global continuous-only" model if possible, else global mean
df_krr_mixed_predict <- function(X_train,
                                 y_train,
                                 X_test,
                                 sigma = NULL,
                                 cat_cols = NULL,
                                 cat_max_levels = 20L,
                                 cat_max_unique_frac = 0.2,
                                 sigma_max_n = 2000L) {
  X_train <- as.matrix(X_train)
  X_test  <- as.matrix(X_test)
  y_train <- as.numeric(y_train)

  ntr <- nrow(X_train)
  nte <- nrow(X_test)

  if (nte == 0L) return(numeric(0))
  if (ntr == 0L) return(rep(NA_real_, nte))

  p <- ncol(X_train)
  if (p != ncol(X_test)) stop("df_krr_mixed_predict(): X_train and X_test must have same number of columns.")

  if (is.null(cat_cols)) {
    cat_cols <- df_infer_categorical_cols(
      X_train,
      max_levels = cat_max_levels,
      max_unique_frac = cat_max_unique_frac
    )
  } else {
    cat_cols <- as.integer(cat_cols)
    cat_cols <- cat_cols[cat_cols >= 1L & cat_cols <= p]
    cat_cols <- sort(unique(cat_cols))
  }

  cont_cols <- setdiff(seq_len(p), cat_cols)

  # If no categorical cols detected, fall back to the original all-continuous behavior.
  if (length(cat_cols) == 0L) {
    if (length(y_train) < 2L) return(rep(mean(y_train), nte))

    if (is.null(sigma)) {
      sigma <- df_estimate_rbf_sigma(X_train, max_n = sigma_max_n)
      if (is.null(sigma)) sigma <- 1
    }

    rbf_kernel <- kernlab::rbfdot(sigma = sigma)
    reg_model  <- kernlab::gausspr(
      x      = as.matrix(X_train),
      y      = y_train,
      kernel = rbf_kernel
    )
    return(df_kernlab_predict(reg_model, as.matrix(X_test)))
  }

  # Mixed: split by categorical key
  key_train <- df_make_cat_key(X_train[, cat_cols, drop = FALSE])
  key_test  <- df_make_cat_key(X_test[,  cat_cols, drop = FALSE])

  # Prepare global fallback model (continuous-only, ignoring categories)
  global_mean <- mean(y_train)
  global_model <- NULL

  if (length(cont_cols) > 0L && length(y_train) >= 2L) {
    if (is.null(sigma)) {
      sigma <- df_estimate_rbf_sigma(X_train[, cont_cols, drop = FALSE], max_n = sigma_max_n)
      if (is.null(sigma)) sigma <- 1
    }
    rbf_kernel <- kernlab::rbfdot(sigma = sigma)

    global_model <- tryCatch(
      kernlab::gausspr(
        x = as.matrix(X_train[, cont_cols, drop = FALSE]),
        y = y_train,
        kernel = rbf_kernel
      ),
      error = function(e) NULL
    )
  } else {
    rbf_kernel <- NULL
  }

  # Fit per-key models
  uniq_keys <- unique(key_train)
  models <- vector("list", length(uniq_keys))
  names(models) <- uniq_keys
  means  <- setNames(rep(NA_real_, length(uniq_keys)), uniq_keys)

  for (k in uniq_keys) {
    idx <- which(key_train == k)
    y_g <- y_train[idx]
    means[k] <- mean(y_g)

    if (length(cont_cols) == 0L || length(y_g) < 2L) {
      models[[k]] <- NULL
    } else {
      X_g <- X_train[idx, cont_cols, drop = FALSE]
      # If sigma not set (e.g. cont cols exist but global sigma couldn't be computed), compute locally
      sigma_g <- sigma
      if (is.null(sigma_g)) {
        sigma_g <- df_estimate_rbf_sigma(X_g, max_n = sigma_max_n)
        if (is.null(sigma_g)) sigma_g <- 1
      }
      kern_g <- kernlab::rbfdot(sigma = sigma_g)

      models[[k]] <- tryCatch(
        kernlab::gausspr(
          x = as.matrix(X_g),
          y = y_g,
          kernel = kern_g
        ),
        error = function(e) NULL
      )
    }
  }

  # Predict
  pred <- rep(global_mean, nte)

  for (k in unique(key_test)) {
    idx <- which(key_test == k)

    if (k %in% uniq_keys) {
      # Seen key: use per-key model/mean
      if (length(cont_cols) == 0L) {
        pred[idx] <- means[[k]]
      } else if (is.null(models[[k]])) {
        pred[idx] <- means[[k]]
      } else {
        X_t <- X_test[idx, cont_cols, drop = FALSE]
        pred[idx] <- df_kernlab_predict(models[[k]], as.matrix(X_t))
      }
    } else {
      # Unseen key: fallback
      if (!is.null(global_model) && length(cont_cols) > 0L) {
        X_t <- X_test[idx, cont_cols, drop = FALSE]
        pred[idx] <- df_kernlab_predict(global_model, as.matrix(X_t))
      } else {
        pred[idx] <- global_mean
      }
    }
  }

  as.numeric(pred)
}

# Convenience: compute (cat_cols, sigma) once for repeated calls on the same X_train
df_krr_mixed_prepare <- function(X_train,
                                 sigma = NULL,
                                 cat_cols = NULL,
                                 cat_max_levels = 20L,
                                 cat_max_unique_frac = 0.2,
                                 sigma_max_n = 2000L) {
  X_train <- as.matrix(X_train)
  p <- ncol(X_train)

  if (is.null(cat_cols)) {
    cat_cols <- df_infer_categorical_cols(
      X_train,
      max_levels = cat_max_levels,
      max_unique_frac = cat_max_unique_frac
    )
  } else {
    cat_cols <- as.integer(cat_cols)
    cat_cols <- cat_cols[cat_cols >= 1L & cat_cols <= p]
    cat_cols <- sort(unique(cat_cols))
  }

  cont_cols <- setdiff(seq_len(p), cat_cols)

  sigma_use <- sigma
  if (is.null(sigma_use) && length(cont_cols) > 0L) {
    sigma_use <- df_estimate_rbf_sigma(X_train[, cont_cols, drop = FALSE], max_n = sigma_max_n)
    if (is.null(sigma_use)) sigma_use <- 1
  }

  list(cat_cols = cat_cols, cont_cols = cont_cols, sigma = sigma_use)
}

############################################################
## Numeric Jacobian & sandwich utilities
############################################################

# Numeric Jacobian (central differences) for a vector-valued function.
# Returns a matrix with dim: length(fun(x0)) x length(x0).
df_numeric_jacobian <- function(fun, x0, eps = NULL) {
  x0 <- as.numeric(x0)
  f0 <- as.numeric(fun(x0))
  p  <- length(x0)
  m  <- length(f0)

  if (is.null(eps)) {
    # Scale-aware step size
    eps <- sqrt(.Machine$double.eps) * pmax(1, abs(x0))
  } else {
    eps <- rep_len(eps, p)
  }

  J <- matrix(NA_real_, nrow = m, ncol = p)
  for (j in seq_len(p)) {
    x_plus  <- x0; x_plus[j]  <- x_plus[j]  + eps[j]
    x_minus <- x0; x_minus[j] <- x_minus[j] - eps[j]
    f_plus  <- as.numeric(fun(x_plus))
    f_minus <- as.numeric(fun(x_minus))
    J[, j]  <- (f_plus - f_minus) / (2 * eps[j])
  }
  J
}

# Compute joint sandwich variance for an M-estimator solving:
#   mean_i s_i(beta) = 0
# where s_i is an n x q matrix of stacked scores.
# A is q x q Jacobian of mean score at beta_hat.
df_joint_sandwich <- function(score_mat, A_hat) {
  score_mat <- as.matrix(score_mat)
  n <- nrow(score_mat)
  q <- ncol(score_mat)

  if (n <= q) {
    return(list(var = matrix(NA_real_, q, q)))
  }

  sbar <- colMeans(score_mat)
  S    <- sweep(score_mat, 2, sbar, "-")

  B_hat <- crossprod(S) / n

  V <- tryCatch({
    A_inv <- solve(A_hat)
    A_inv %*% B_hat %*% t(A_inv) / n
  }, error = function(e) {
    matrix(NA_real_, q, q)
  })

  list(var = V, A = A_hat, B = B_hat)
}

############################################################
## 1. Basis functions for NP logistic model (multi-X)
############################################################

# l = (1, X, y) for NP inclusion model
g2 <- function(dat) {
  X <- df_get_X(dat)
  cbind(1, X, as.numeric(dat$y))
}
g3 <- g2
g4 <- g2

############################################################
## 2. Chang & Kott type estimating eq. for pi_NP(phi)
############################################################

pi_np.est_simple <- function(dat, h2, h3, h4, p.use = TRUE) {
  function(phi) {
    phi <- as.numeric(phi)
    X    <- df_get_X(dat)
    y    <- as.numeric(dat$y)
    l    <- as.matrix(cbind(1, X, y))
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
##    (UPDATED: auto mixed-type KRR for categorical columns in X)
############################################################

# mu(x) = E[Y|X=x] using (d_p==1 & d_np==0) data
regression_expectation_kernlab <- function(dat, new_X, sigma = NULL) {
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 & as.numeric(dat$d_np) == 0)
  X_obs <- X_all[idx, , drop = FALSE]
  y_obs <- as.numeric(dat$y[idx])

  if (length(y_obs) < 1L) {
    return(rep(NA_real_, nrow(new_X)))
  }
  if (length(y_obs) < 2L) {
    return(rep(mean(y_obs), nrow(new_X)))
  }

  cat_cols <- df_infer_categorical_cols(X_obs)

  prep <- df_krr_mixed_prepare(X_obs, sigma = sigma, cat_cols = cat_cols)
  df_krr_mixed_predict(
    X_train = X_obs,
    y_train = y_obs,
    X_test  = as.matrix(new_X),
    sigma   = prep$sigma,
    cat_cols = prep$cat_cols
  )
}

# Estimate pi_P(L) via regression of 1/pi_P on L=(X,y) using d_p==1 units
pi_p_estimation_kernlab <- function(dat, new_L, sigma = NULL) {
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 & is.finite(as.numeric(dat$pi_p)))

  X_obs    <- X_all[idx, , drop = FALSE]
  y_obs    <- as.numeric(dat$y[idx])
  L_obs    <- cbind(X_obs, y_obs)
  pi_p_obs <- as.numeric(dat$pi_p[idx])

  if (length(pi_p_obs) < 1L) {
    return(rep(NA_real_, nrow(new_L)))
  }
  if (length(pi_p_obs) < 2L) {
    return(rep(mean(pi_p_obs), nrow(new_L)))
  }

  # categorical columns are the same as in X, last column (y) is continuous
  cat_cols_X <- df_infer_categorical_cols(X_obs)
  cat_cols_L <- cat_cols_X

  prep <- df_krr_mixed_prepare(L_obs, sigma = sigma, cat_cols = cat_cols_L)
  reg_pred <- df_krr_mixed_predict(
    X_train = L_obs,
    y_train = 1 / pi_p_obs,
    X_test  = as.matrix(new_L),
    sigma   = prep$sigma,
    cat_cols = prep$cat_cols
  )

  as.numeric(1 / reg_pred)
}

# eta4*(X;phi): E[eta4_numer(L;phi) | X] / E[h4_denom(L;phi) | X]
estimate_conditional_expectation_kernlab_phi <- function(dat, phi, new_X, sigma = NULL) {
  phi   <- as.numeric(phi)
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 | as.numeric(dat$d_np) == 1)

  X_obs    <- X_all[idx, , drop = FALSE]
  y_obs    <- as.numeric(dat$y[idx])
  pi_p_obs <- as.numeric(dat$pi_p[idx])

  l_obs  <- as.matrix(cbind(1, X_obs, y_obs))
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

  if (nrow(X_obs) < 2L) {
    denom_pred <- rep(mean(Denom_vals), nrow(new_X))
    numer_pred <- matrix(rep(colMeans(Numer_mat), each = nrow(new_X)),
                         nrow = nrow(new_X), ncol = p_dim, byrow = TRUE)
    return(sweep(numer_pred, 1, denom_pred, "/"))
  }

  cat_cols <- df_infer_categorical_cols(X_obs)
  prep <- df_krr_mixed_prepare(X_obs, sigma = sigma, cat_cols = cat_cols)

  denom_pred <- df_krr_mixed_predict(
    X_train = X_obs,
    y_train = as.numeric(Denom_vals),
    X_test  = as.matrix(new_X),
    sigma   = prep$sigma,
    cat_cols = prep$cat_cols
  )

  numer_pred <- sapply(
    seq_len(p_dim),
    function(j) df_krr_mixed_predict(
      X_train = X_obs,
      y_train = as.numeric(Numer_mat[, j]),
      X_test  = as.matrix(new_X),
      sigma   = prep$sigma,
      cat_cols = prep$cat_cols
    )
  )

  sweep(numer_pred, 1, denom_pred, "/")
}

# h4*(X;phi): E[h4_numer(L;phi) | X] / E[h4_denom(L;phi) | X]
estimate_conditional_expectation_kernlab_theta <- function(dat, phi, new_X, sigma = NULL) {
  phi   <- as.numeric(phi)
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 | as.numeric(dat$d_np) == 1)

  X_obs    <- X_all[idx, , drop = FALSE]
  y_obs    <- as.numeric(dat$y[idx])
  pi_p_obs <- as.numeric(dat$pi_p[idx])

  l_obs  <- as.matrix(cbind(1, X_obs, y_obs))
  eta    <- as.numeric(l_obs %*% phi)
  pi_np  <- 1 / (1 + exp(-eta))

  Denom_vals <- h4_prob_denom_function(pi_np, pi_p_obs, phi)
  Numer_vals <- h4_prob_numer_function(pi_np, pi_p_obs, phi, y_obs)

  if (nrow(X_obs) < 2L) {
    denom_pred <- rep(mean(Denom_vals), nrow(new_X))
    numer_pred <- rep(mean(Numer_vals), nrow(new_X))
    return(as.numeric(numer_pred / denom_pred))
  }

  cat_cols <- df_infer_categorical_cols(X_obs)
  prep <- df_krr_mixed_prepare(X_obs, sigma = sigma, cat_cols = cat_cols)

  denom_pred <- df_krr_mixed_predict(
    X_train = X_obs,
    y_train = as.numeric(Denom_vals),
    X_test  = as.matrix(new_X),
    sigma   = prep$sigma,
    cat_cols = prep$cat_cols
  )
  numer_pred <- df_krr_mixed_predict(
    X_train = X_obs,
    y_train = as.numeric(Numer_vals),
    X_test  = as.matrix(new_X),
    sigma   = prep$sigma,
    cat_cols = prep$cat_cols
  )

  as.numeric(numer_pred / denom_pred)
}

############################################################
## 5. Efficient scores: per-observation contributions for phi and theta
############################################################

# Per-observation score contributions for phi (n x p_phi)
# This is the "est_eq" matrix inside estimating_equation_optimal_phi().
df_score_phi_contrib <- function(dat1, phi, eta4_star_local) {
  phi <- as.numeric(phi)

  X1    <- df_get_X(dat1)
  y1    <- as.numeric(dat1$y)
  d_p1  <- as.numeric(dat1$d_p)
  d_np1 <- as.numeric(dat1$d_np)
  pi_p1 <- as.numeric(dat1$pi_p)

  l1 <- as.matrix(cbind(1, X1, y1))

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

  eta4_star_local <- as.matrix(eta4_star_local)
  if (nrow(eta4_star_local) != nrow(dat1)) {
    stop("df_score_phi_contrib(): eta4_star_local must have nrow(dat1).")
  }

  term1 + term2_coef * eta4_star_local
}

# Per-observation contribution for theta (n-vector), given phi and h4_star
efficient_theta_contrib <- function(dat1, phi, h4_star_local) {
  phi <- as.numeric(phi)

  X1    <- df_get_X(dat1)
  y1    <- as.numeric(dat1$y)
  d_p1  <- as.numeric(dat1$d_p)
  d_np1 <- as.numeric(dat1$d_np)
  pi_p1 <- as.numeric(dat1$pi_p)

  l1 <- as.matrix(cbind(1, X1, y1))

  eta      <- as.numeric(l1 %*% phi)
  pi_np1   <- 1 / (1 + exp(-eta))
  pi_np_p1 <- pi_np1 + pi_p1 - pi_np1 * pi_p1

  h4_star_local <- as.numeric(h4_star_local)

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

############################################################
## 6. DML helpers: folds and pi_P cross-fitting (used for both DML1/DML2)
############################################################

make_folds <- function(n, K) {
  fold_id <- sample(rep(1:K, length.out = n))
  split(seq_len(n), fold_id)
}

# Cross-fit imputation of pi_p for units with (d_np==1 & d_p==0) where pi_p may be NA.
# Note: We impute only for the fold's held-out units (standard cross-fitting).
impute_pi_p_crossfit <- function(dat, folds, sigma = NULL, progress = FALSE) {
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

    L_f  <- as.matrix(cbind(X_f, y_f))

    pi_p_mis <- which(d_np_f == 1 & d_p_f == 0 & !is.finite(as.numeric(dat_fold$pi_p)))
    if (length(pi_p_mis) > 0) {
      tilde_pi <- pi_p_estimation_kernlab(
        dat_train,
        new_L = L_f[pi_p_mis, , drop = FALSE],
        sigma = sigma
      )
      dat_cf$pi_p[idx_k[pi_p_mis]] <- tilde_pi
    }

    if (progress) utils::setTxtProgressBar(pb, k)
  }

  if (progress) {
    close(pb)
    cat("\n")
  }

  dat_cf
}

############################################################
## 7. Basic estimators: P, NP, NP+P (for comparison)
############################################################

df_sandwich_from_contrib <- function(contrib, level = 0.95) {
  contrib <- as.numeric(contrib)
  contrib <- contrib[is.finite(contrib)]
  n <- length(contrib)

  if (n <= 1L) {
    return(list(theta = NA_real_, var = NA_real_, se = NA_real_, ci = c(NA_real_, NA_real_)))
  }

  theta_hat <- mean(contrib)
  s2        <- stats::var(contrib)
  var_hat   <- s2 / n
  se_hat    <- sqrt(var_hat)

  z  <- stats::qnorm(1 - (1 - level) / 2)
  ci <- c(theta_hat - z * se_hat, theta_hat + z * se_hat)

  list(theta = theta_hat, var = var_hat, se = se_hat, ci = ci)
}

df_estimate_P <- function(dat) {
  contrib   <- as.numeric(dat$d_p) * as.numeric(dat$y) / as.numeric(dat$pi_p)
  theta_res <- df_sandwich_from_contrib(contrib)
  list(theta = theta_res$theta, var = theta_res$var, se = theta_res$se, ci = theta_res$ci)
}

df_estimate_NP <- function(dat, phi_start = NULL, max_iter = 20) {
  X <- df_get_X(dat)
  p <- ncol(X)

  if (is.null(phi_start)) {
    phi_start <- c(-log(1 / mean(dat$d_np) - 1), rep(0, p), 0)
  }

  phi_est <- list(termcd = 99)
  it <- 0
  while (phi_est$termcd > 2 && it < max_iter) {
    init <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
    phi_est <- nleqslv::nleqslv(init, pi_np.est_simple(dat, g2, g3, g4, p.use = FALSE))
    it <- it + 1
  }

  if (phi_est$termcd > 2) {
    phi_hat   <- rep(NA_real_, length(phi_start))
    theta_res <- df_sandwich_from_contrib(rep(NA_real_, nrow(dat)))
  } else {
    phi_hat <- as.numeric(phi_est$x)
    l       <- as.matrix(cbind(1, X, as.numeric(dat$y)))
    pi_np   <- 1 / (1 + exp(-as.numeric(l %*% phi_hat)))
    contrib <- as.numeric(dat$d_np) * as.numeric(dat$y) / pi_np
    theta_res <- df_sandwich_from_contrib(contrib)
  }

  list(phi = phi_hat, theta = theta_res$theta, var = theta_res$var, se = theta_res$se, ci = theta_res$ci)
}

df_estimate_NP_P <- function(dat, phi_start = NULL, max_iter = 20) {
  X <- df_get_X(dat)
  p <- ncol(X)

  if (is.null(phi_start)) {
    phi_start <- c(-log(1 / mean(dat$d_np) - 1), rep(0, p), 0)
  }

  phi_est <- list(termcd = 99)
  it <- 0
  while (phi_est$termcd > 2 && it < max_iter) {
    init <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
    phi_est <- nleqslv::nleqslv(init, pi_np.est_simple(dat, g2, g3, g4, p.use = TRUE))
    it <- it + 1
  }

  if (phi_est$termcd > 2) {
    phi_hat   <- rep(NA_real_, length(phi_start))
    theta_res <- df_sandwich_from_contrib(rep(NA_real_, nrow(dat)))
  } else {
    phi_hat <- as.numeric(phi_est$x)
    l       <- as.matrix(cbind(1, X, as.numeric(dat$y)))
    pi_np   <- 1 / (1 + exp(-as.numeric(l %*% phi_hat)))

    pi_p  <- as.numeric(dat$pi_p)
    d_np  <- as.numeric(dat$d_np)
    d_p   <- as.numeric(dat$d_p)
    denom <- pi_np + pi_p - pi_np * pi_p

    contrib <- as.numeric(dat$y) * (d_np + d_p - d_np * d_p) / denom
    theta_res <- df_sandwich_from_contrib(contrib)
  }

  list(phi = phi_hat, theta = theta_res$theta, var = theta_res$var, se = theta_res$se, ci = theta_res$ci)
}

############################################################
## 8. Efficient estimator Eff (DML2): JOINT (phi, theta) sandwich
############################################################

efficient_estimator_dml2 <- function(dat,
                                     phi_start   = NULL,
                                     K           = 2,
                                     max_restart = 10,
                                     progress    = FALSE) {
  n <- nrow(dat)
  X <- df_get_X(dat)
  p_x <- ncol(X)
  p_phi <- 1 + p_x + 1  # (1, X, y)

  if (is.null(phi_start)) {
    phi_start <- c(-log(1 / mean(dat$d_np) - 1), rep(0, p_x), 0)
  }

  folds <- make_folds(n, K)

  # Step 1: cross-fit pi_p (impute only if pi_p has NA for NP-only)
  dat_cf <- impute_pi_p_crossfit(dat, folds, sigma = NULL, progress = progress)

  # Objective for phi: squared norm of aggregated mean score
  obj_phi <- function(phi) {
    phi <- as.numeric(phi)
    eq_agg <- rep(0, length(phi))

    for (k in seq_along(folds)) {
      idx_test  <- folds[[k]]
      idx_train <- setdiff(seq_len(n), idx_test)

      dat_test  <- dat_cf[idx_test,  ]
      dat_train <- dat_cf[idx_train, ]

      X_test <- df_get_X(dat_test)
      eta4_k <- estimate_conditional_expectation_kernlab_phi(dat_train, phi, new_X = X_test)

      sphi_k <- df_score_phi_contrib(dat_test, phi, eta4_k)
      eq_k   <- colMeans(sphi_k)

      w_k <- length(idx_test) / n
      eq_agg <- eq_agg + w_k * eq_k
    }

    sum(eq_agg^2)
  }

  if (progress) {
    cat("Step 2/3: solving phi (Eff, DML2, Nelder-Mead) ...\n")
    flush.console()
  }

  attempt <- 0
  res     <- list(convergence = 1, par = phi_start)

  while (attempt < max_restart && res$convergence != 0) {
    attempt <- attempt + 1
    init <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)

    res_try <- try(
      stats::optim(par = init, fn = obj_phi, method = "Nelder-Mead"),
      silent = TRUE
    )
    if (!inherits(res_try, "try-error")) {
      res <- res_try
      if (progress) {
        cat(sprintf("  attempt %d/%d done (convergence=%d)\n", attempt, max_restart, res$convergence))
        flush.console()
      }
      if (res$convergence == 0) break
    }
  }

  phi_hat <- as.numeric(res$par)

  # Step 3: cross-fit eta4*(X) and h4*(X) at phi_hat, then compute theta_hat
  if (progress) {
    cat("Step 3/3: cross-fitting nuisances for joint scores ...\n")
    pb <- utils::txtProgressBar(min = 0, max = length(folds), style = 3)
  }

  eta4_all <- matrix(NA_real_, nrow = n, ncol = p_phi)
  h4_all   <- rep(NA_real_, n)

  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(n), idx_test)

    dat_test  <- dat_cf[idx_test,  ]
    dat_train <- dat_cf[idx_train, ]

    X_test <- df_get_X(dat_test)

    eta4_k <- estimate_conditional_expectation_kernlab_phi(dat_train, phi_hat, new_X = X_test)
    h4_k   <- estimate_conditional_expectation_kernlab_theta(dat_train, phi_hat, new_X = X_test)

    eta4_all[idx_test, ] <- eta4_k
    h4_all[idx_test]     <- h4_k

    if (progress) utils::setTxtProgressBar(pb, k)
  }

  if (progress) {
    close(pb)
    cat("\n")
  }

  contrib_theta <- efficient_theta_contrib(dat_cf, phi_hat, h4_all)
  theta_hat     <- mean(contrib_theta)

  # Build stacked scores s_i(beta_hat) with nuisance fixed at phi_hat
  s_phi   <- df_score_phi_contrib(dat_cf, phi_hat, eta4_all)
  s_theta <- as.numeric(contrib_theta - theta_hat)

  score_mat <- cbind(s_phi, s_theta)

  # Build Jacobian A_hat (holding eta4_all and h4_all fixed)
  g_phi_mean <- function(phi) colMeans(df_score_phi_contrib(dat_cf, phi, eta4_all))
  A11 <- df_numeric_jacobian(g_phi_mean, phi_hat)

  g_theta_mean <- function(phi) {
    mean(efficient_theta_contrib(dat_cf, phi, h4_all) - theta_hat)
  }
  A21 <- df_numeric_jacobian(g_theta_mean, phi_hat)  # 1 x p_phi

  A_hat <- rbind(
    cbind(A11, rep(0, p_phi)),
    cbind(A21, -1)
  )

  js <- df_joint_sandwich(score_mat, A_hat)
  V  <- js$var

  phi_var <- V[1:p_phi, 1:p_phi, drop = FALSE]
  th_var  <- V[p_phi + 1, p_phi + 1]

  phi_se <- sqrt(diag(phi_var))
  th_se  <- sqrt(th_var)

  z <- stats::qnorm(0.975)
  phi_ci <- rbind(phi_hat - z * phi_se, phi_hat + z * phi_se)
  rownames(phi_ci) <- c("lower", "upper")
  theta_ci <- c(theta_hat - z * th_se, theta_hat + z * th_se)

  list(
    phi    = phi_hat,
    phi_var= phi_var,
    phi_se = phi_se,
    phi_ci = phi_ci,
    theta  = theta_hat,
    var    = th_var,
    se     = th_se,
    ci     = theta_ci,
    info   = list(
      type        = "Eff",
      dml_type    = "DML2",
      K           = K,
      phi_start   = phi_start,
      max_restart = max_restart,
      progress    = progress,
      convergence = res$convergence
    )
  )
}

############################################################
## 9. Efficient estimator Eff (DML1): per-fold (phi_k, theta_k) + JOINT sandwich combine
############################################################

efficient_estimator_dml1 <- function(dat,
                                     phi_start   = NULL,
                                     K           = 2,
                                     max_restart = 10,
                                     progress    = FALSE) {
  n <- nrow(dat)
  X <- df_get_X(dat)
  p_x <- ncol(X)
  p_phi <- 1 + p_x + 1

  if (is.null(phi_start)) {
    phi_start <- c(-log(1 / mean(dat$d_np) - 1), rep(0, p_x), 0)
  }

  folds <- make_folds(n, K)

  # Cross-fit pi_p once (impute for each held-out fold if needed)
  dat_cf <- impute_pi_p_crossfit(dat, folds, sigma = NULL, progress = FALSE)

  phi_k_mat   <- matrix(NA_real_, nrow = K, ncol = p_phi)
  theta_k     <- rep(NA_real_, K)
  var_k_list  <- vector("list", K)
  w_k         <- rep(NA_real_, K)  # weights = n_k / n

  if (progress) {
    cat("DML1: per-fold estimation for (phi, theta) and joint sandwich ...\n")
    flush.console()
  }

  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(n), idx_test)
    w_k[k]    <- length(idx_test) / n

    dat_test  <- dat_cf[idx_test,  ]
    dat_train <- dat_cf[idx_train, ]

    # Solve phi_k on the test fold score with nuisance trained on training fold
    ee_fun_k <- function(phi) {
      X_test <- df_get_X(dat_test)
      eta4_k <- estimate_conditional_expectation_kernlab_phi(dat_train, phi, new_X = X_test)
      colMeans(df_score_phi_contrib(dat_test, phi, eta4_k))
    }

    attempt <- 0
    sol <- list(termcd = 99)
    while (attempt < max_restart && sol$termcd > 2) {
      attempt <- attempt + 1
      init <- phi_start + stats::runif(length(phi_start), -0.1, 0.1)
      sol_try <- try(nleqslv::nleqslv(init, ee_fun_k), silent = TRUE)
      if (!inherits(sol_try, "try-error")) sol <- sol_try
    }

    if (sol$termcd > 2) {
      if (progress) cat(sprintf("  fold %d/%d: phi did not converge.\n", k, K))
      next
    }

    phi_k <- as.numeric(sol$x)
    phi_k_mat[k, ] <- phi_k

    # Nuisance on test fold at phi_k (trained on training fold)
    X_test <- df_get_X(dat_test)
    eta4_k <- estimate_conditional_expectation_kernlab_phi(dat_train, phi_k, new_X = X_test)
    h4_k   <- estimate_conditional_expectation_kernlab_theta(dat_train, phi_k, new_X = X_test)

    # Scores on test fold
    s_phi_k   <- df_score_phi_contrib(dat_test, phi_k, eta4_k)
    contrib_k <- efficient_theta_contrib(dat_test, phi_k, h4_k)
    theta_k[k] <- mean(contrib_k)
    s_theta_k <- as.numeric(contrib_k - theta_k[k])

    score_k <- cbind(s_phi_k, s_theta_k)

    # Jacobian on fold k (holding eta4_k and h4_k fixed)
    g_phi_mean_k <- function(phi) colMeans(df_score_phi_contrib(dat_test, phi, eta4_k))
    A11_k <- df_numeric_jacobian(g_phi_mean_k, phi_k)

    g_theta_mean_k <- function(phi) {
      mean(efficient_theta_contrib(dat_test, phi, h4_k) - theta_k[k])
    }
    A21_k <- df_numeric_jacobian(g_theta_mean_k, phi_k)

    A_k <- rbind(
      cbind(A11_k, rep(0, p_phi)),
      cbind(A21_k, -1)
    )

    js_k <- df_joint_sandwich(score_k, A_k)
    var_k_list[[k]] <- js_k$var

    if (progress) {
      cat(sprintf("  fold %d/%d done.\n", k, K))
      flush.console()
    }
  }

  valid <- which(is.finite(theta_k) & apply(is.finite(phi_k_mat), 1, all))
  if (length(valid) == 0) {
    return(list(
      phi = rep(NA_real_, p_phi),
      phi_var = matrix(NA_real_, p_phi, p_phi),
      phi_se = rep(NA_real_, p_phi),
      phi_ci = matrix(NA_real_, nrow = 2, ncol = p_phi,
                      dimnames = list(c("lower","upper"), NULL)),
      theta = NA_real_,
      var = NA_real_,
      se = NA_real_,
      ci = c(NA_real_, NA_real_),
      info = list(type = "Eff", dml_type = "DML1", K = K, progress = progress)
    ))
  }

  # Weighted average of fold estimates (weights by fold size)
  w_use <- w_k[valid]
  w_use <- w_use / sum(w_use)

  phi_hat   <- as.numeric(t(w_use) %*% phi_k_mat[valid, , drop = FALSE])
  theta_hat <- as.numeric(sum(w_use * theta_k[valid]))

  # Combine variances of fold estimators: Var(sum w_k beta_k) ~ sum w_k^2 Var(beta_k)
  V_hat <- matrix(0, nrow = p_phi + 1, ncol = p_phi + 1)
  for (j in seq_along(valid)) {
    kk <- valid[j]
    V_hat <- V_hat + (w_use[j]^2) * var_k_list[[kk]]
  }

  phi_var <- V_hat[1:p_phi, 1:p_phi, drop = FALSE]
  th_var  <- V_hat[p_phi + 1, p_phi + 1]

  phi_se <- sqrt(diag(phi_var))
  th_se  <- sqrt(th_var)

  z <- stats::qnorm(0.975)
  phi_ci <- rbind(phi_hat - z * phi_se, phi_hat + z * phi_se)
  rownames(phi_ci) <- c("lower", "upper")
  theta_ci <- c(theta_hat - z * th_se, theta_hat + z * th_se)

  list(
    phi    = phi_hat,
    phi_var= phi_var,
    phi_se = phi_se,
    phi_ci = phi_ci,
    theta  = theta_hat,
    var    = th_var,
    se     = th_se,
    ci     = theta_ci,
    info   = list(
      type        = "Eff",
      dml_type    = "DML1",
      K           = K,
      phi_start   = phi_start,
      max_restart = max_restart,
      progress    = progress
    )
  )
}

############################################################
## 10. Sub-efficient estimator Eff_S (DML2-style): theta only
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

subefficient_estimator_dml2 <- function(dat, K = 2, progress = FALSE) {
  n <- nrow(dat)
  folds <- make_folds(n, K)

  X_all <- df_get_X(dat)
  idx_mu <- which(as.numeric(dat$d_p) == 1 & as.numeric(dat$d_np) == 0)

  # UPDATED: estimate sigma using continuous columns only
  cat_cols <- df_infer_categorical_cols(X_all)
  cont_cols <- setdiff(seq_len(ncol(X_all)), cat_cols)

  if (length(idx_mu) >= 2L && length(cont_cols) > 0L) {
    sigma_mu <- df_estimate_rbf_sigma(X_all[idx_mu, cont_cols, drop = FALSE], max_n = 2000L)
  } else {
    sigma_mu <- NULL
  }

  if (progress) {
    cat("Cross-fitting mu(X) for Eff_S ...\n")
    pb <- utils::txtProgressBar(min = 0, max = length(folds), style = 3)
  }

  mu_hat_all <- rep(NA_real_, n)

  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(n), idx_test)

    dat_train <- dat[idx_train, ]
    dat_test  <- dat[idx_test,  ]
    X_test    <- df_get_X(dat_test)

    mu_k <- regression_expectation_kernlab(dat_train, new_X = X_test, sigma = sigma_mu)
    mu_hat_all[idx_test] <- mu_k

    if (progress) utils::setTxtProgressBar(pb, k)
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
       info  = list(type = "Eff_S", K = K, progress = progress))
}

############################################################
## 11. Parametric efficient estimator Eff_P: theta + phi (phi via nleqslv), joint variance NOT implemented
##     (theta variance is computed from contributions as before)
############################################################

efficient_parametric_estimator <- function(dat,
                                           phi_start = NULL,
                                           eta4_star = 0,
                                           max_iter  = 20,
                                           progress  = FALSE) {
  X_all <- df_get_X(dat)
  p_x <- ncol(X_all)
  p_phi <- 1 + p_x + 1
  n <- nrow(dat)

  if (is.null(phi_start)) {
    phi_start <- c(-log(1 / mean(dat$d_np) - 1), rep(0, p_x), 0)
  }

  # Solve phi with estimating equation using eta4_star fixed (working model)
  ee_para <- function(phi) {
    X1    <- df_get_X(dat)
    y1    <- as.numeric(dat$y)
    d_p1  <- as.numeric(dat$d_p)
    d_np1 <- as.numeric(dat$d_np)
    pi_p1 <- as.numeric(dat$pi_p)
    l1    <- as.matrix(cbind(1, X1, y1))

    phi <- as.numeric(phi)
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

    est_eq <- term1 + term2_coef * 0
    colMeans(est_eq)
  }

  sol <- list(termcd = 99)
  it <- 0
  while (sol$termcd > 2 && it < max_iter) {
    it <- it + 1
    sol_try <- try(nleqslv::nleqslv(phi_start, ee_para), silent = TRUE)
    if (!inherits(sol_try, "try-error")) sol <- sol_try
    if (progress) {
      cat(sprintf("Solving phi for Eff_P: iter %d/%d (termcd=%s)\n", it, max_iter, sol$termcd))
      flush.console()
    }
  }

  if (sol$termcd > 2) {
    phi_hat <- rep(NA_real_, p_phi)
    theta_res <- df_sandwich_from_contrib(rep(NA_real_, n))
  } else {
    phi_hat <- as.numeric(sol$x)

    # Working outcome regression: linear model y ~ X on (d_np==1 | d_p==1)
    subset_idx <- which(as.numeric(dat$d_np) == 1 | as.numeric(dat$d_p) == 1)
    dat_sub <- dat[subset_idx, ]
    X_sub <- df_get_X(dat_sub)
    df_sub <- data.frame(y = as.numeric(dat_sub$y), X_sub)
    colnames(df_sub)[-1] <- paste0("x", seq_len(ncol(X_sub)))

    lm_fit <- stats::lm(y ~ ., data = df_sub)

    X_full <- df_get_X(dat)
    df_full <- data.frame(X_full)
    colnames(df_full) <- paste0("x", seq_len(ncol(X_full)))
    mu_hat <- stats::predict(lm_fit, newdata = df_full)

    contrib <- efficient_theta_contrib(dat, phi_hat, h4_star_local = mu_hat)
    theta_res <- df_sandwich_from_contrib(contrib)
  }

  list(
    phi    = phi_hat,
    phi_var= NULL,  # joint variance for Eff_P not implemented in this file
    phi_se = NULL,
    phi_ci = NULL,
    theta  = theta_res$theta,
    var    = theta_res$var,
    se     = theta_res$se,
    ci     = theta_res$ci,
    info   = list(type = "Eff_P", phi_start = phi_start, eta4_star = eta4_star, max_iter = max_iter, progress = progress)
  )
}

############################################################
## 12. Public wrappers: Eff, Eff_S, Eff_P
############################################################

# Eff: choose DML1 or DML2 via dml_type (or alias `type`)
Eff <- function(dat,
                K           = 2,
                phi_start   = NULL,
                max_restart = 10,
                type        = NULL,
                dml_type    = 2,
                progress    = interactive()) {

  if (!is.null(type)) dml_type <- type

  if (is.numeric(dml_type)) dml_type <- if (dml_type == 1) "DML1" else "DML2"
  dml_type <- toupper(as.character(dml_type))
  if (!dml_type %in% c("DML1", "DML2")) stop("Eff(): dml_type must be 1, 2, 'DML1', or 'DML2'.")

  if (dml_type == "DML2") {
    efficient_estimator_dml2(dat, phi_start = phi_start, K = K, max_restart = max_restart, progress = progress)
  } else {
    efficient_estimator_dml1(dat, phi_start = phi_start, K = K, max_restart = max_restart, progress = progress)
  }
}

Eff_S <- function(dat, K = 2, progress = interactive()) {
  subefficient_estimator_dml2(dat, K = K, progress = progress)
}

Eff_P <- function(dat,
                  phi_start = NULL,
                  eta4_star = 0,
                  max_iter  = 20,
                  progress  = interactive()) {
  efficient_parametric_estimator(dat, phi_start = phi_start, eta4_star = eta4_star, max_iter = max_iter, progress = progress)
}
