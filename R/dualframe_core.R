############################################################
## dualframe_core.R
##
## Mixed continuous/discrete covariates:
##  - Nonparametric nuisances use delta x RBF (mixed-type KRR)
##  - Small categorical cell warning / ask-to-continue (English)
##
## Parametric part for phi:
##  - If X is a data.frame with categorical columns, we use model.matrix(~ .)
##    (drop intercept) so categorical vars enter as L-1 dummies
##
## Binary Y support (requested):
##  - Eff     : unchanged (still uses Gaussian / mixed kernel regression for eta4*, h4*)
##  - Eff_S   : if y is binary and logit=TRUE (or logit=NULL auto-detect), mu(X)=P(Y=1|X)
##              is estimated by kernel logistic regression (RBF + ridge, IRLS)
##  - Eff_P   : if y is binary and logit=TRUE (or logit=NULL auto-detect), mu(X) is estimated
##              by parametric logistic regression (glm binomial)
##
## NEW: x_info flag (speed / (0,0)-absent setting)
##  - x_info=TRUE  (default): keep the original behavior (estimate h4*, eta4* by KRR for Eff)
##  - x_info=FALSE: assume no (d_np,d_p)=(0,0) units are available (or X unavailable there).
##                  We then set h4*=eta4*=0 and skip their estimation entirely (faster).
##                  Any (0,0) rows accidentally included are dropped inside wrappers.
##
## NEW (population size + Hajek option; output-preserving by default):
##  - df_estimate_P(), df_estimate_NP(), df_estimate_NP_P() now accept:
##      N     : population size (frame size). Required if dat does NOT include (d_p,d_np)=(0,0) rows
##              and hajek=FALSE.
##      hajek : if TRUE, compute the Hajek (ratio) version of the estimator (does not require N).
##  - Default behavior is unchanged when full frame data are provided (N is inferred as nrow(dat)).
##
## IMPORTANT (Eff + x_info=FALSE):
##  - When x_info=FALSE and the data do not include (d_p,d_np)=(0,0) rows (e.g., a union-only dataset),
##    Eff() can still be interpreted as targeting the full-frame mean E(Y) ONLY if the population size N
##    is supplied. In that case, theta (and its SE/CI) is scaled by n_obs / N, where n_obs = nrow(dat)
##    after dropping (0,0) rows.
##
## Stability updates:
##  - RBF sigma estimation uses deterministic subsampling (reproducible)
##  - Probability clipping for pi_p, pi_np, pi_np_p to avoid Inf/NaN
##  - Denominator floor for ratio nuisances eta4*/h4* to avoid division by ~0
##  - DML1: safe fallback and fold skipping if estimating equations are non-finite
##  - Sandwich variance: drop non-finite score rows; robust inversion with ridge;
##    theta variance fallback to contrib-based sandwich when needed
##
## IMPORTANT (Sigma rule):
##  - RBF sigma rule (legacy / old version):
##      sigma = 1 / median(||x_i - x_j||)
##    where the median is taken over pairwise distances of a deterministic subsample (for speed).
##
## Output-preserving speed updates (this revision):
##  - Cache df_krr_mixed_prepare() per fold for eta4*/h4*
##  - Vectorize Numer_mat construction in eta4*(X;phi)
##  - Faster mixed prediction: pre-split indices; fit global model only if needed
##  - df_get_X fast-path for numeric matrices
##  - Internal fast score evaluators for DML objectives/jacobians (avoid repeated df_get_X)
##
## Hotfix (Scenario S1/S2 gausspr failure):
##  - If dat is created by data.frame(X = I(matrix(...))), then dat$X has class "AsIs".
##    kernlab::gausspr is an S4 generic and can error for x="AsIs" even if the object is matrix-like.
##  - We strip only the "AsIs" class inside df_apply_x_info() and df_get_X().
############################################################

############################################################
## 0. Utilities
############################################################

.dfSEDI_state <- new.env(parent = emptyenv())
.dfSEDI_state$allow_small_cells   <- FALSE
.dfSEDI_state$seen_smallcell_keys <- character(0)

df_ask_yes_no <- function(question, default = FALSE) {
  if (requireNamespace("utils", quietly = TRUE) &&
      exists("askYesNo", envir = asNamespace("utils"), inherits = FALSE)) {
    ans <- get("askYesNo", envir = asNamespace("utils"))(question, default = default)
    if (is.na(ans)) return(isTRUE(default))
    return(isTRUE(ans))
  }

  prompt <- if (isTRUE(default)) " [Y/n] " else " [y/N] "
  ans <- readline(paste0(question, prompt))
  ans <- trimws(tolower(ans))
  if (ans == "") return(isTRUE(default))
  ans %in% c("y", "yes")
}

# Clip probabilities into (eps, 1-eps)
df_clip_prob <- function(p, eps = NULL) {
  if (is.null(eps)) eps <- getOption("dfSEDI.prob_clip", 1e-8)
  eps <- as.numeric(eps)[1]
  if (!is.finite(eps) || eps <= 0 || eps >= 0.5) eps <- 1e-8
  p <- as.numeric(p)
  p <- pmin(pmax(p, eps), 1 - eps)
  p
}

# Floor a numeric vector at a positive value
df_floor_pos <- function(x, floor = NULL) {
  if (is.null(floor)) floor <- getOption("dfSEDI.denom_floor", 1e-10)
  floor <- as.numeric(floor)[1]
  if (!is.finite(floor) || floor <= 0) floor <- 1e-10
  x <- as.numeric(x)
  x[!is.finite(x)] <- floor
  x <- pmax(x, floor)
  x
}


# Multi-start optimization helper
# - Generate n_candidates initial values: center + Unif[-width, width]
# - Evaluate objective at all candidates
# - Try optim() from the best candidates in order (up to n_try)
# - Stop at the first successful optimization (convergence == 0)
df_multistart_optim <- function(obj_fun,
                                center,
                                n_candidates = 100L,
                                n_try = 10L,
                                width = 0.5,
                                method = "Nelder-Mead",
                                control = NULL) {
  center <- as.numeric(center)
  tol_obj <- getOption("dfSEDI.phi_obj_tol", 1e-4)
  if (!is.finite(tol_obj) || tol_obj <= 0) tol_obj <- 1e-4

  if (!all(is.finite(center))) stop("df_multistart_optim(): 'center' must be finite.")

  p <- length(center)
  n_candidates <- as.integer(n_candidates)[1]
  if (!is.finite(n_candidates) || n_candidates < 1) n_candidates <- 100L
  n_try <- as.integer(n_try)[1]
  if (!is.finite(n_try) || n_try < 1) n_try <- 10L
  n_try <- min(n_try, n_candidates)

  width <- as.numeric(width)[1]
  if (!is.finite(width) || width < 0) width <- 0.5
  if (is.null(control)) control <- list()

  # Candidate initial values: center + Unif[-width, width]
  init_mat <- matrix(stats::runif(n_candidates * p, min = -width, max = width),
                     nrow = n_candidates, ncol = p)
  init_mat <- sweep(init_mat, 2, center, "+")

  # Evaluate objective
  vals <- rep(Inf, n_candidates)
  for (i in seq_len(n_candidates)) {
    v <- tryCatch(obj_fun(init_mat[i, ]), error = function(e) Inf)
    if (!is.finite(v)) v <- Inf
    vals[i] <- v
  }

  ord <- order(vals)

  # Try optim() from the best starts (up to n_try)
  for (j in seq_len(n_try)) {
    idx <- ord[j]
    start <- init_mat[idx, ]
    res <- tryCatch(
      stats::optim(par = start, fn = obj_fun, method = method, control = control),
      error = function(e) NULL
    )
    if (!is.null(res) &&
        is.list(res) &&
        is.finite(res$value) &&
        all(is.finite(res$par)) &&
        (isTRUE(res$convergence == 0) || res$value <= tol_obj)) {
      return(list(ok = TRUE, res = res, start = start, start_rank = j, init_vals = vals, init_order = ord))
    }
  }

  list(ok = FALSE, res = NULL, init_vals = vals, init_order = ord)
}

# Prepare pi_p when it is observed only for d_p==1 units.
# - For union-sample units (d_p==1 or d_np==1), pi_p must be finite (impute before if needed).
# - For non-union units (d_p==0 & d_np==0), pi_p can be NA; we fill with a harmless value so
#   algebraic expressions remain finite (these rows have zero influence on terms that need pi_p).
df_prepare_pi_p <- function(pi_p, d_p, d_np, context = "dfSEDI") {
  pi_p <- as.numeric(pi_p)
  d_p  <- as.numeric(d_p)
  d_np <- as.numeric(d_np)

  union <- (d_p == 1 | d_np == 1)

  if (any(union & !is.finite(pi_p))) {
    stop(
      context, ": `pi_p` is missing (NA/Inf) for some union-sample units (d_p==1 or d_np==1).\n",
      "Impute pi_p for non-P units in the union sample first (e.g., impute_pi_p_crossfit() or a model-based imputation).",
      call. = FALSE
    )
  }

  fill_val <- mean(pi_p[is.finite(pi_p)], na.rm = TRUE)
  if (!is.finite(fill_val)) fill_val <- 0.5

  pi_p[!union & !is.finite(pi_p)] <- fill_val
  df_clip_prob(pi_p)
}


# Parametric (linear) imputation of pi_p for NP-only union units (d_np==1, d_p==0).
# We mirror the *target* of the nonparametric imputer used in Eff:
#   regress 1/pi_p on L = (X, y) using P-sample units (d_p==1), then predict for missing.
# Cross-fitting is intentionally NOT used here (parametric model).
df_impute_pi_p_lm <- function(dat,
                              context = "dfSEDI",
                              min_train = 5L) {
  if (!is.data.frame(dat)) return(dat)
  if (!("pi_p" %in% names(dat))) return(dat)
  if (!("d_p" %in% names(dat))) return(dat)
  if (!("d_np" %in% names(dat))) return(dat)

  d_p  <- as.numeric(dat$d_p)
  d_np <- as.numeric(dat$d_np)
  pi_p <- as.numeric(dat$pi_p)

  # Only need pi_p for union-sample units; in typical "pi_p observed only when d_p==1"
  # setups, the missing cases are (d_np, d_p) = (1, 0).
  idx_mis <- which(d_p == 0 & d_np == 1 & !is.finite(pi_p))
  if (length(idx_mis) == 0L) return(dat)

  idx_tr <- which(d_p == 1 & is.finite(pi_p))
  if (length(idx_tr) < min_train) {
    # Fallback: mean-impute (best-effort).
    fill <- mean(pi_p[idx_tr], na.rm = TRUE)
    if (!is.finite(fill)) fill <- 0.5
    pi_p[idx_mis] <- fill
    dat$pi_p <- pi_p
    return(dat)
  }

  # Design matrix for regression inputs
  X <- df_get_X(dat)
  if (!("y" %in% names(dat))) {
    y <- rep(0, nrow(X))
  } else {
    y <- as.numeric(dat$y)
  }

  # L = (X, y)
  L <- cbind(X, y)

  # Train: 1/pi_p ~ L  (linear model)
  train_df <- as.data.frame(L[idx_tr, , drop = FALSE])
  colnames(train_df) <- c(paste0("x", seq_len(ncol(X))), "y")
  train_df$resp <- 1 / df_clip_prob(pi_p[idx_tr])

  fit <- tryCatch(
    stats::lm(resp ~ ., data = train_df),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    fill <- mean(pi_p[idx_tr], na.rm = TRUE)
    if (!is.finite(fill)) fill <- 0.5
    pi_p[idx_mis] <- fill
    dat$pi_p <- pi_p
    return(dat)
  }

  new_df <- as.data.frame(L[idx_mis, , drop = FALSE])
  colnames(new_df) <- c(paste0("x", seq_len(ncol(X))), "y")

  pred_resp <- tryCatch(
    stats::predict(fit, newdata = new_df),
    error = function(e) rep(NA_real_, nrow(new_df))
  )
  pred_resp <- as.numeric(pred_resp)

  # Guard: predicted 1/pi_p must be positive & finite
  resp_bar <- mean(train_df$resp, na.rm = TRUE)
  if (!is.finite(resp_bar) || resp_bar <= 0) resp_bar <- 1 / df_clip_prob(mean(pi_p[idx_tr], na.rm = TRUE))

  pred_resp[!is.finite(pred_resp) | pred_resp <= 0] <- resp_bar

  pi_hat <- 1 / pred_resp
  pi_hat <- df_clip_prob(pi_hat)

  pi_p[idx_mis] <- pi_hat
  dat$pi_p <- pi_p
  dat
}

# Apply x_info setting:
#  - x_info=TRUE : keep dat as is
#  - x_info=FALSE: drop any (d_p,d_np)=(0,0) rows (if present)
df_apply_x_info <- function(dat, x_info = TRUE) {
  if (!is.data.frame(dat)) return(dat)

  # Hotfix: strip only "AsIs" (from data.frame(X = I(matrix(...)))).
  # kernlab::gausspr is an S4 generic and can fail for x="AsIs".
  if ("X" %in% names(dat) && inherits(dat$X, "AsIs")) {
    x0 <- dat$X
    class(x0) <- setdiff(class(x0), "AsIs")
    dat$X <- x0
  }

  if (isTRUE(x_info)) return(dat)
  if (!all(c("d_p", "d_np") %in% names(dat))) return(dat)

  d_p  <- as.numeric(dat$d_p)
  d_np <- as.numeric(dat$d_np)
  keep <- !(d_p == 0 & d_np == 0)

  if (any(!keep, na.rm = TRUE)) {
    # silent by default (speed path); set option to warn if desired
    if (isTRUE(getOption("dfSEDI.warn_drop_00", FALSE))) {
      warning(sprintf("dfSEDI: x_info=FALSE dropped %d rows with (d_p,d_np)=(0,0).",
                      sum(!keep, na.rm = TRUE)), call. = FALSE)
    }
    dat <- dat[keep, , drop = FALSE]
  }

  # (After subsetting, ensure AsIs didn't persist)
  if ("X" %in% names(dat) && inherits(dat$X, "AsIs")) {
    x0 <- dat$X
    class(x0) <- setdiff(class(x0), "AsIs")
    dat$X <- x0
  }

  dat
}

# Binary Y helpers ------------------------------------------------------------

df_is_binary_y <- function(y) {
  y0 <- y
  y0 <- y0[!is.na(y0)]
  if (length(y0) == 0L) return(FALSE)

  if (is.logical(y0)) return(TRUE)

  if (is.factor(y0)) {
    return(length(levels(y0)) == 2L)
  }

  if (is.character(y0)) {
    # try numeric coercion; if fails, treat as factor-ish
    tmp <- suppressWarnings(as.numeric(y0))
    if (all(is.na(tmp))) {
      return(length(unique(y0)) == 2L)
    } else {
      u <- sort(unique(tmp[is.finite(tmp)]))
      return(length(u) == 2L)
    }
  }

  if (is.numeric(y0) || is.integer(y0)) {
    u <- sort(unique(as.numeric(y0[is.finite(y0)])))
    return(length(u) == 2L)
  }

  FALSE
}

# Convert binary y to 0/1 numeric.
# - factor: second level is treated as success (1) by default
# - numeric: larger unique value is treated as success (1)
df_as_binary01 <- function(y, success_level = NULL) {
  if (is.logical(y)) {
    y01 <- as.numeric(y)
    return(list(y01 = y01, mapping = "logical(TRUE)=1,FALSE=0"))
  }

  if (is.factor(y)) {
    if (length(levels(y)) != 2L) stop("df_as_binary01(): factor y must have exactly 2 levels.")
    lev <- levels(y)
    if (is.null(success_level)) success_level <- lev[2]
    y01 <- as.numeric(y == success_level)
    return(list(y01 = y01, mapping = paste0("factor(success='", success_level, "')")))
  }

  if (is.character(y)) {
    tmp <- suppressWarnings(as.numeric(y))
    if (all(is.na(tmp))) {
      f <- factor(y)
      return(df_as_binary01(f, success_level = success_level))
    } else {
      return(df_as_binary01(tmp, success_level = success_level))
    }
  }

  if (is.numeric(y) || is.integer(y)) {
    ynum <- as.numeric(y)
    u <- sort(unique(ynum[is.finite(ynum)]))
    if (length(u) != 2L) stop("df_as_binary01(): numeric y must have exactly 2 unique finite values.")
    success_val <- if (is.null(success_level)) u[2] else as.numeric(success_level)
    y01 <- as.numeric(ynum == success_val)
    return(list(y01 = y01, mapping = paste0("numeric(success=", success_val, ", other=", u[1], ")")))
  }

  stop("df_as_binary01(): unsupported y type.")
}

df_sigmoid <- function(z) {
  z <- as.numeric(z)
  out <- rep(NA_real_, length(z))
  idx1 <- which(z >= 0)
  idx0 <- which(z < 0)
  out[idx1] <- 1 / (1 + exp(-z[idx1]))
  ez <- exp(z[idx0])
  out[idx0] <- ez / (1 + ez)
  out
}

df_infer_logit_flag <- function(logit, y) {
  if (is.null(logit)) {
    return(df_is_binary_y(y))
  }
  isTRUE(logit)
}

# Extract parametric design matrix X
df_get_X <- function(dat) {
  if ("X" %in% names(dat)) {
    X <- dat$X

    # Hotfix: strip only "AsIs" (from data.frame(X = I(matrix(...)))).
    # kernlab::gausspr is an S4 generic and can fail for x="AsIs".
    if (inherits(X, "AsIs")) {
      class(X) <- setdiff(class(X), "AsIs")
    }

    if (is.data.frame(X)) {
      X_df <- X

      for (nm in names(X_df)) {
        col <- X_df[[nm]]

        if (is.character(col)) {
          tmp <- suppressWarnings(as.numeric(col))
          if (all(is.finite(tmp) | is.na(col))) {
            X_df[[nm]] <- tmp
          } else {
            X_df[[nm]] <- factor(col)
          }
        } else if (is.logical(col)) {
          X_df[[nm]] <- factor(col)
        } else if (is.factor(col)) {
          lev <- levels(col)
          lev_num <- suppressWarnings(as.numeric(lev))
          if (length(lev_num) == length(lev) && all(is.finite(lev_num))) {
            ord <- order(lev_num)
            X_df[[nm]] <- factor(col, levels = lev[ord])
          }
        }
      }

      mm <- stats::model.matrix(~ ., data = X_df)
      assign_vec <- attr(mm, "assign")

      if ("(Intercept)" %in% colnames(mm)) {
        keep <- colnames(mm) != "(Intercept)"
        mm <- mm[, keep, drop = FALSE]
        assign_vec <- assign_vec[keep]
      }

      mm <- as.matrix(mm)
      storage.mode(mm) <- "numeric"

      term_labels <- attr(stats::terms(~ ., data = X_df), "term.labels")
      fac_names <- names(X_df)[sapply(X_df, is.factor)]
      is_cat_term <- term_labels %in% fac_names
      cat_term_idx <- which(is_cat_term)
      cat_cols_mm <- which(assign_vec %in% cat_term_idx)

      attr(mm, "dfSEDI_categorical_cols") <- as.integer(cat_cols_mm)
      return(mm)
    }

    if (is.data.frame(X)) X <- as.matrix(X)
    if (!is.matrix(X)) stop("df_get_X(): `X` must be a matrix or convertible to matrix.")

    # speed: if already numeric matrix, return as-is (just enforce numeric storage + attr)
    if (is.matrix(X) && (is.double(X) || is.integer(X))) {
      X_num <- X
      storage.mode(X_num) <- "numeric"

      # ensure AsIs is not present (S4 dispatch safety)
      if (inherits(X_num, "AsIs")) {
        class(X_num) <- setdiff(class(X_num), "AsIs")
      }

      if (is.null(attr(X_num, "dfSEDI_categorical_cols"))) {
        attr(X_num, "dfSEDI_categorical_cols") <- integer(0)
      }
      return(X_num)
    }

    X_num <- apply(X, 2, function(col) suppressWarnings(as.numeric(col)))
    X_num <- as.matrix(X_num)
    storage.mode(X_num) <- "numeric"
    if (inherits(X_num, "AsIs")) {
      class(X_num) <- setdiff(class(X_num), "AsIs")
    }
    attr(X_num, "dfSEDI_categorical_cols") <- integer(0)
    return(X_num)
  }

  if ("x" %in% names(dat)) {
    X_num <- matrix(as.numeric(dat$x), ncol = 1)
    attr(X_num, "dfSEDI_categorical_cols") <- integer(0)
    return(X_num)
  }

  stop("df_get_X(): cannot extract X. Provide either `X` or `x`.")
}

# Build l = [1, X, y] as numeric matrix (used in fast score paths)
df_make_l_matrix <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  l <- as.matrix(cbind(1, X, y))
  storage.mode(l) <- "numeric"
  l
}

# Evaluate and validate base_fun(X) for pi_np estimation (strict checks).
# Expected convention in this file: X is n x p and base_fun(X) returns n x (p+2).
df_eval_base_fun <- function(base_fun,
                             X,
                             expected_ncol = NULL,
                             context = "base_fun") {
  if (!is.function(base_fun)) {
    stop(context, ": `base_fun` must be a function.", call. = FALSE)
  }

  X <- as.matrix(X)
  n <- nrow(X)

  b_raw <- tryCatch(
    base_fun(X),
    error = function(e) {
      stop(context, ": base_fun(X) failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (is.data.frame(b_raw)) b_raw <- as.matrix(b_raw)

  # Allow a length-n vector (treated as a single-column basis)
  if (is.null(dim(b_raw))) {
    if (length(b_raw) == n) {
      b_raw <- matrix(b_raw, nrow = n, ncol = 1)
    } else {
      stop(
        context, ": base_fun(X) returned a vector of length ", length(b_raw),
        " but expected length ", n, " (nrow(X)).",
        call. = FALSE
      )
    }
  }

  if (!is.matrix(b_raw)) {
    stop(context, ": base_fun(X) must return a matrix (or a data.frame coercible to matrix).",
         call. = FALSE)
  }

  if (nrow(b_raw) != n) {
    stop(
      context, ": base_fun(X) must have the same number of rows as X.\n",
      "nrow(X) = ", n, ", but nrow(base_fun(X)) = ", nrow(b_raw), ".",
      call. = FALSE
    )
  }

  if (ncol(b_raw) < 1L) {
    stop(context, ": base_fun(X) must have at least 1 column.", call. = FALSE)
  }

  if (!is.null(expected_ncol)) {
    expected_ncol <- as.integer(expected_ncol)[1]
    if (is.finite(expected_ncol) && expected_ncol >= 1L) {
      if (ncol(b_raw) != expected_ncol) {
        stop(
          context, ": unexpected number of columns in base_fun(X).\n",
          "Expected ", expected_ncol, " columns, but got ", ncol(b_raw), ".",
          call. = FALSE
        )
      }
    }
  }

  b <- as.matrix(b_raw)
  storage.mode(b) <- "numeric"

  if (any(!is.finite(b))) {
    idx <- which(!is.finite(b), arr.ind = TRUE)
    show <- idx[seq_len(min(5L, nrow(idx))), , drop = FALSE]
    show_txt <- paste(apply(show, 1, function(rc) paste0("(", rc[1], ",", rc[2], ")")), collapse = ", ")
    stop(
      context, ": base_fun(X) must return finite numeric values (no NA/Inf).\n",
      "Example non-finite positions (row,col): ", show_txt, ".",
      call. = FALSE
    )
  }

  b
}

df_kernlab_predict <- function(model, newdata) {
  newdata <- as.matrix(newdata)

  pred <- tryCatch({
    if (requireNamespace("kernlab", quietly = TRUE) &&
        exists("predict", envir = asNamespace("kernlab"), inherits = FALSE)) {
      get("predict", envir = asNamespace("kernlab"))(model, newdata)
    } else {
      stats::predict(model, newdata)
    }
  }, error = function(e1) {
    tryCatch({
      stats::predict(model, newdata)
    }, error = function(e2) {
      stop(
        "dfSEDI: prediction failed for a kernlab::gausspr model.\n",
        "Please try running `library(kernlab)` once and re-run.\n",
        "Original error: ", conditionMessage(e2)
      )
    })
  })

  as.numeric(pred)
}

############################################################
## Mixed-type KRR helpers
############################################################

df_detect_categorical_cols_heuristic <- function(X,
                                                 max_levels = 20L,
                                                 max_unique_frac = 0.2,
                                                 always_categorical_levels = 3L) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  if (p == 0L || n <= 1L) return(integer(0))

  out <- logical(p)
  for (j in seq_len(p)) {
    x <- X[, j]
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      out[j] <- FALSE
      next
    }
    u <- length(unique(x))
    frac <- u / max(1, length(x))
    out[j] <- (u <= always_categorical_levels) || ((u <= max_levels) && (frac <= max_unique_frac))
  }

  which(out)
}

df_infer_categorical_cols <- function(X,
                                      max_levels = 20L,
                                      max_unique_frac = 0.2,
                                      always_categorical_levels = 3L) {
  X <- as.matrix(X)
  cat_attr <- attr(X, "dfSEDI_categorical_cols")
  if (is.null(cat_attr)) cat_attr <- integer(0)
  cat_attr <- as.integer(cat_attr)

  cat_heur <- df_detect_categorical_cols_heuristic(
    X,
    max_levels = max_levels,
    max_unique_frac = max_unique_frac,
    always_categorical_levels = always_categorical_levels
  )

  sort(unique(c(cat_attr, cat_heur)))
}

# Deterministic subsampling for sigma estimation (reproducible).
# Legacy sigma rule (old version): sigma = 1 / median(||x_i-x_j||)
df_estimate_rbf_sigma <- function(X_cont, max_n = 2000L) {
  X_cont <- as.matrix(X_cont)
  n <- nrow(X_cont)
  if (n < 2L) return(NULL)

  max_n <- as.integer(max_n)
  if (is.na(max_n) || max_n < 2L) max_n <- n

  if (n > max_n) {
    idx <- unique(pmax(1L, pmin(n, as.integer(round(seq.int(1, n, length.out = max_n))))))
    X_use <- X_cont[idx, , drop = FALSE]
  } else {
    X_use <- X_cont
  }

  d <- stats::dist(X_use)
  dv <- as.numeric(d)
  dv <- dv[is.finite(dv) & dv > 0]
  if (length(dv) == 0L) return(1)

  med <- stats::median(dv)
  if (!is.finite(med) || med <= 0) return(1)

  1 / med
}

df_make_cat_key <- function(X_cat) {
  X_cat <- as.matrix(X_cat)
  n <- nrow(X_cat)
  if (n == 0L) return(character(0))
  if (ncol(X_cat) == 0L) return(rep("ALL", n))
  dfc <- as.data.frame(X_cat, stringsAsFactors = FALSE)
  do.call(paste, c(dfc, sep = "\r"))
}

df_check_small_categorical_cells <- function(X_train,
                                             cat_cols,
                                             context = "mixed-type kernel regression") {
  X_train <- as.matrix(X_train)
  p <- ncol(X_train)
  cat_cols <- sort(unique(as.integer(cat_cols)))
  cat_cols <- cat_cols[cat_cols >= 1L & cat_cols <= p]
  if (length(cat_cols) == 0L) return(invisible(TRUE))

  if (isTRUE(getOption("dfSEDI.allow_small_cells", FALSE))) return(invisible(TRUE))
  if (isTRUE(.dfSEDI_state$allow_small_cells)) return(invisible(TRUE))

  min_cell <- getOption("dfSEDI.min_cell_size", 30L)
  min_cell <- as.integer(min_cell)
  if (!is.finite(min_cell) || min_cell <= 0L) return(invisible(TRUE))

  action <- getOption("dfSEDI.small_cell_action",
                      if (interactive()) "ask" else "warn")
  action <- tolower(as.character(action))[1]
  if (!action %in% c("ask", "warn", "stop", "none")) {
    action <- if (interactive()) "ask" else "warn"
  }
  if (action == "none") return(invisible(TRUE))

  key <- df_make_cat_key(X_train[, cat_cols, drop = FALSE])
  tab <- table(key)
  if (length(tab) == 0L) return(invisible(TRUE))

  min_n <- min(tab)
  n_small <- sum(tab < min_cell)
  if (min_n >= min_cell) return(invisible(TRUE))

  id <- paste0(context, "|n=", nrow(X_train),
               "|cats=", paste(cat_cols, collapse = ","),
               "|thr=", min_cell,
               "|min=", as.integer(min_n),
               "|nsmall=", as.integer(n_small),
               "|ncells=", length(tab))
  if (id %in% .dfSEDI_state$seen_smallcell_keys) return(invisible(TRUE))
  .dfSEDI_state$seen_smallcell_keys <- unique(c(.dfSEDI_state$seen_smallcell_keys, id))

  q10 <- as.numeric(stats::quantile(as.numeric(tab), probs = 0.10, names = FALSE, type = 7))
  q25 <- as.numeric(stats::quantile(as.numeric(tab), probs = 0.25, names = FALSE, type = 7))
  q50 <- as.numeric(stats::quantile(as.numeric(tab), probs = 0.50, names = FALSE, type = 7))

  msg <- paste0(
    "dfSEDI warning: very small categorical cell sizes detected in ", context, ".\n",
    "With a delta kernel, mixed-type KRR effectively fits separate models within each categorical cell.\n",
    "Small cells can make nuisance estimates unstable.\n",
    "Training sample summary: n = ", nrow(X_train),
    ", #cells = ", length(tab),
    ", min cell size = ", as.integer(min_n),
    ", cells <", min_cell, " = ", as.integer(n_small), ".\n",
    "Cell-size quantiles: 10%=", round(q10, 2),
    ", 25%=", round(q25, 2),
    ", 50%=", round(q50, 2), ".\n",
    "Categorical columns (indices in the design matrix): ", paste(cat_cols, collapse = ","), ".\n",
    "You can control this via: options(dfSEDI.min_cell_size=...), options(dfSEDI.small_cell_action='ask'/'warn'/'stop'/'none'),\n",
    "or bypass via: options(dfSEDI.allow_small_cells=TRUE)."
  )

  if (action == "warn") {
    warning(msg, call. = FALSE)
    return(invisible(TRUE))
  }
  if (action == "stop") {
    stop(msg, call. = FALSE)
  }

  warning(msg, call. = FALSE)
  if (!interactive()) {
    stop(
      paste0(
        "Non-interactive session detected. To proceed, set options(dfSEDI.small_cell_action='warn') ",
        "or options(dfSEDI.allow_small_cells=TRUE)."
      ),
      call. = FALSE
    )
  }

  ok <- df_ask_yes_no("Continue anyway?", default = FALSE)
  if (!ok) stop("Aborted by user due to small categorical cell sizes.", call. = FALSE)

  .dfSEDI_state$allow_small_cells <- TRUE
  invisible(TRUE)
}

df_krr_mixed_prepare <- function(X_train,
                                 sigma = NULL,
                                 cat_cols = NULL,
                                 cat_max_levels = 20L,
                                 cat_max_unique_frac = 0.2,
                                 always_categorical_levels = 3L,
                                 sigma_max_n = 2000L,
                                 context = "mixed-type kernel regression") {
  X_train <- as.matrix(X_train)
  p <- ncol(X_train)

  if (is.null(cat_cols)) {
    cat_cols <- df_infer_categorical_cols(
      X_train,
      max_levels = cat_max_levels,
      max_unique_frac = cat_max_unique_frac,
      always_categorical_levels = always_categorical_levels
    )
  } else {
    cat_cols <- sort(unique(as.integer(cat_cols)))
    cat_cols <- cat_cols[cat_cols >= 1L & cat_cols <= p]
  }

  cont_cols <- setdiff(seq_len(p), cat_cols)

  if (length(cat_cols) > 0L) {
    df_check_small_categorical_cells(X_train, cat_cols, context = context)
  }

  sigma_use <- sigma
  if (is.null(sigma_use) && length(cont_cols) > 0L) {
    sigma_use <- df_estimate_rbf_sigma(X_train[, cont_cols, drop = FALSE], max_n = sigma_max_n)
    if (is.null(sigma_use)) sigma_use <- 1
  }

  list(cat_cols = cat_cols, cont_cols = cont_cols, sigma = sigma_use)
}

df_krr_mixed_predict <- function(X_train,
                                 y_train,
                                 X_test,
                                 sigma,
                                 cat_cols,
                                 cont_cols,
                                 sigma_max_n = 2000L) {
  X_train <- as.matrix(X_train)
  X_test  <- as.matrix(X_test)
  y_train <- as.numeric(y_train)

  ntr <- nrow(X_train)
  nte <- nrow(X_test)
  if (nte == 0L) return(numeric(0))
  if (ntr == 0L) return(rep(NA_real_, nte))

  if (length(cat_cols) == 0L) {
    if (length(y_train) < 2L) return(rep(mean(y_train), nte))
    if (is.null(sigma)) {
      sigma <- df_estimate_rbf_sigma(X_train, max_n = sigma_max_n)
      if (is.null(sigma)) sigma <- 1
    }
    rbf_kernel <- kernlab::rbfdot(sigma = sigma)
    reg_model  <- kernlab::gausspr(x = X_train, y = y_train, kernel = rbf_kernel)
    return(df_kernlab_predict(reg_model, X_test))
  }

  key_train <- df_make_cat_key(X_train[, cat_cols, drop = FALSE])
  key_test  <- df_make_cat_key(X_test[,  cat_cols, drop = FALSE])

  split_train <- split(seq_len(ntr), key_train)
  split_test  <- split(seq_len(nte), key_test)
  uniq_keys <- names(split_train)

  global_mean <- mean(y_train)

  # Fit global model only if needed (i.e., test contains keys not in train)
  need_global <- length(setdiff(names(split_test), uniq_keys)) > 0L

  global_model <- NULL
  if (isTRUE(need_global) && length(cont_cols) > 0L && length(y_train) >= 2L) {
    sig_g <- sigma
    if (is.null(sig_g)) {
      sig_g <- df_estimate_rbf_sigma(X_train[, cont_cols, drop = FALSE], max_n = sigma_max_n)
      if (is.null(sig_g)) sig_g <- 1
    }
    ker_g <- kernlab::rbfdot(sigma = sig_g)
    global_model <- tryCatch(
      kernlab::gausspr(x = X_train[, cont_cols, drop = FALSE], y = y_train, kernel = ker_g),
      error = function(e) NULL
    )
  }

  means <- vapply(split_train, function(idx) mean(y_train[idx]), numeric(1))
  models <- setNames(vector("list", length(uniq_keys)), uniq_keys)

  for (k in uniq_keys) {
    idx <- split_train[[k]]
    y_g <- y_train[idx]

    if (length(cont_cols) == 0L || length(y_g) < 2L) {
      models[[k]] <- NULL
    } else {
      X_g <- X_train[idx, cont_cols, drop = FALSE]
      sig_k <- sigma
      if (is.null(sig_k)) {
        sig_k <- df_estimate_rbf_sigma(X_g, max_n = sigma_max_n)
        if (is.null(sig_k)) sig_k <- 1
      }
      ker_k <- kernlab::rbfdot(sigma = sig_k)
      models[[k]] <- tryCatch(
        kernlab::gausspr(x = X_g, y = y_g, kernel = ker_k),
        error = function(e) NULL
      )
    }
  }

  pred <- rep(global_mean, nte)

  for (k in names(split_test)) {
    idx <- split_test[[k]]

    if (k %in% uniq_keys) {
      if (length(cont_cols) == 0L || is.null(models[[k]])) {
        pred[idx] <- means[[k]]
      } else {
        X_t <- X_test[idx, cont_cols, drop = FALSE]
        pred[idx] <- df_kernlab_predict(models[[k]], X_t)
      }
    } else {
      if (!is.null(global_model) && length(cont_cols) > 0L) {
        X_t <- X_test[idx, cont_cols, drop = FALSE]
        pred[idx] <- df_kernlab_predict(global_model, X_t)
      } else {
        pred[idx] <- global_mean
      }
    }
  }

  as.numeric(pred)
}

############################################################
## Kernel logistic regression (binary mu) for Eff_S
############################################################

# RBF kernel matrix: exp(-sigma * ||x-x'||^2)
df_rbf_kernel <- function(X1, X2 = NULL, sigma = 1) {
  X1 <- as.matrix(X1)
  if (is.null(X2)) X2 <- X1
  X2 <- as.matrix(X2)
  sigma <- as.numeric(sigma)[1]
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1

  # squared distances: ||a-b||^2 = ||a||^2 + ||b||^2 - 2 a'b
  a2 <- rowSums(X1^2)
  b2 <- rowSums(X2^2)
  G  <- X1 %*% t(X2)
  D2 <- outer(a2, b2, "+") - 2 * G
  D2[D2 < 0] <- 0
  exp(-sigma * D2)
}

df_solve_ridge <- function(A, b, ridge_seq = NULL) {
  A <- as.matrix(A)
  b <- as.matrix(b)

  if (is.null(ridge_seq)) {
    ridge_seq <- getOption("dfSEDI.solve_ridge_seq", c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2))
  }
  ridge_seq <- as.numeric(ridge_seq)
  ridge_seq <- ridge_seq[is.finite(ridge_seq) & ridge_seq >= 0]
  if (length(ridge_seq) == 0L) ridge_seq <- c(0, 1e-8, 1e-6, 1e-4, 1e-2)

  for (lam in ridge_seq) {
    A_try <- A
    if (lam > 0) A_try <- A_try + diag(lam, nrow(A_try))
    sol <- tryCatch(solve(A_try, b), error = function(e) NULL)
    if (!is.null(sol) && all(is.finite(sol))) return(sol)
  }
  NULL
}

# Kernel logistic regression fit (IRLS) with intercept:
# minimize -loglik(y | f) + (lambda/2) * alpha' K alpha,  f = K alpha + b
df_klogit_fit <- function(X, y01, sigma = NULL, lambda = NULL, max_iter = NULL, tol = NULL) {
  X <- as.matrix(X)
  y01 <- as.numeric(y01)

  n <- nrow(X)
  if (n < 2L) return(list(ok = FALSE))

  if (is.null(sigma)) sigma <- df_estimate_rbf_sigma(X, max_n = min(2000L, n))
  if (is.null(sigma)) sigma <- 1

  if (is.null(lambda)) lambda <- getOption("dfSEDI.klr_lambda", 1e-2)
  lambda <- as.numeric(lambda)[1]
  if (!is.finite(lambda) || lambda <= 0) lambda <- 1e-2

  if (is.null(max_iter)) max_iter <- getOption("dfSEDI.klr_max_iter", 50L)
  max_iter <- as.integer(max_iter)
  if (!is.finite(max_iter) || max_iter < 1L) max_iter <- 50L

  if (is.null(tol)) tol <- getOption("dfSEDI.klr_tol", 1e-6)
  tol <- as.numeric(tol)[1]
  if (!is.finite(tol) || tol <= 0) tol <- 1e-6

  # If only one class present, return constant prob
  u <- sort(unique(y01[is.finite(y01)]))
  if (length(u) < 2L) {
    p0 <- df_clip_prob(mean(y01), eps = getOption("dfSEDI.prob_clip", 1e-8))
    return(list(ok = TRUE, type = "constant", p = p0, X = X, sigma = sigma))
  }

  # deterministic subsample if too large (optional)
  max_n <- getOption("dfSEDI.klr_max_n", 2000L)
  max_n <- as.integer(max_n)
  if (!is.finite(max_n) || max_n < 2L) max_n <- n
  if (n > max_n) {
    idx <- unique(pmax(1L, pmin(n, as.integer(round(seq.int(1, n, length.out = max_n))))))
    X <- X[idx, , drop = FALSE]
    y01 <- y01[idx]
    n <- nrow(X)
  }

  K <- df_rbf_kernel(X, X, sigma = sigma)
  jitter <- getOption("dfSEDI.klr_jitter", 1e-10)
  jitter <- as.numeric(jitter)[1]
  if (!is.finite(jitter) || jitter < 0) jitter <- 0
  if (jitter > 0) K <- K + diag(jitter, n)

  # init
  pbar <- mean(y01)
  pbar <- df_clip_prob(pbar)
  b <- as.numeric(stats::qlogis(pbar))
  alpha <- rep(0, n)

  eps_w <- getOption("dfSEDI.klr_w_floor", 1e-6)
  eps_w <- as.numeric(eps_w)[1]
  if (!is.finite(eps_w) || eps_w <= 0) eps_w <- 1e-6

  for (it in seq_len(max_iter)) {
    f <- as.numeric(K %*% alpha + b)
    p <- df_sigmoid(f)
    p <- df_clip_prob(p)

    w <- p * (1 - p)
    w <- pmax(w, eps_w)

    z <- f + (y01 - p) / w

    WK <- sweep(K, 1, w, "*")              # diag(w) K
    KWK <- K %*% WK                        # K diag(w) K
    Kw1 <- as.numeric(K %*% w)             # K diag(w) 1
    K_wz <- as.numeric(K %*% (w * z))      # K diag(w) z

    one_W_K <- colSums(WK)                 # 1' diag(w) K
    one_W_1 <- sum(w)
    one_W_z <- sum(w * z)

    A11 <- KWK + lambda * K
    A12 <- matrix(Kw1, nrow = n, ncol = 1)
    A21 <- matrix(one_W_K, nrow = 1, ncol = n)
    A22 <- matrix(one_W_1, nrow = 1, ncol = 1)

    A <- rbind(cbind(A11, A12),
               cbind(A21, A22))
    rhs <- c(K_wz, one_W_z)

    sol <- df_solve_ridge(A, rhs)
    if (is.null(sol)) {
      return(list(ok = FALSE))
    }

    alpha_new <- as.numeric(sol[1:n])
    b_new     <- as.numeric(sol[n + 1])

    diff <- max(abs(alpha_new - alpha), abs(b_new - b))
    alpha <- alpha_new
    b <- b_new
    if (is.finite(diff) && diff < tol) break
  }

  list(ok = TRUE, type = "klogit", alpha = alpha, b = b, X = X, sigma = sigma)
}

df_klogit_predict <- function(fit, X_new) {
  X_new <- as.matrix(X_new)
  if (is.null(fit) || !isTRUE(fit$ok)) return(rep(NA_real_, nrow(X_new)))

  if (isTRUE(fit$type == "constant")) {
    p0 <- df_clip_prob(fit$p)
    return(rep(p0, nrow(X_new)))
  }

  K_new <- df_rbf_kernel(fit$X, X_new, sigma = fit$sigma)
  f_new <- as.numeric(t(K_new) %*% fit$alpha + fit$b)
  p_new <- df_sigmoid(f_new)
  df_clip_prob(p_new)
}

# Mixed-type kernel logistic (delta on categorical, RBF on continuous).
df_klogit_mixed_predict <- function(X_train,
                                    y01_train,
                                    X_test,
                                    sigma,
                                    cat_cols,
                                    cont_cols,
                                    sigma_max_n = 2000L) {
  X_train <- as.matrix(X_train)
  X_test  <- as.matrix(X_test)
  y01_train <- as.numeric(y01_train)

  ntr <- nrow(X_train)
  nte <- nrow(X_test)
  if (nte == 0L) return(numeric(0))
  if (ntr == 0L) return(rep(NA_real_, nte))

  # global fallback (fit only if needed)
  global_mean <- df_clip_prob(mean(y01_train))

  key_train <- if (length(cat_cols) > 0L) df_make_cat_key(X_train[, cat_cols, drop = FALSE]) else rep("ALL", ntr)
  key_test  <- if (length(cat_cols) > 0L) df_make_cat_key(X_test[,  cat_cols, drop = FALSE]) else rep("ALL", nte)

  split_train <- split(seq_len(ntr), key_train)
  split_test  <- split(seq_len(nte), key_test)
  uniq_keys <- names(split_train)

  need_global <- length(setdiff(names(split_test), uniq_keys)) > 0L

  global_fit <- NULL
  if (isTRUE(need_global) && length(cont_cols) > 0L) {
    Xc <- X_train[, cont_cols, drop = FALSE]
    if (is.null(sigma)) sigma <- df_estimate_rbf_sigma(Xc, max_n = sigma_max_n)
    if (is.null(sigma)) sigma <- 1
    global_fit <- df_klogit_fit(Xc, y01_train, sigma = sigma)
    if (!isTRUE(global_fit$ok)) global_fit <- NULL
  }

  if (length(cat_cols) == 0L) {
    if (is.null(global_fit)) return(rep(global_mean, nte))
    return(df_klogit_predict(global_fit, X_test[, cont_cols, drop = FALSE]))
  }

  cell_fit <- setNames(vector("list", length(uniq_keys)), uniq_keys)
  cell_mean <- vapply(split_train, function(idx) df_clip_prob(mean(y01_train[idx])), numeric(1))

  for (k in uniq_keys) {
    idx <- split_train[[k]]
    yk <- y01_train[idx]

    if (length(cont_cols) == 0L) {
      cell_fit[[k]] <- NULL
      next
    }

    uk <- unique(yk[is.finite(yk)])
    if (length(uk) < 2L || length(yk) < 3L) {
      cell_fit[[k]] <- NULL
      next
    }

    Xk <- X_train[idx, cont_cols, drop = FALSE]
    sig_k <- sigma
    if (is.null(sig_k)) sig_k <- df_estimate_rbf_sigma(Xk, max_n = sigma_max_n)
    if (is.null(sig_k)) sig_k <- 1

    fitk <- df_klogit_fit(Xk, yk, sigma = sig_k)
    if (!isTRUE(fitk$ok)) fitk <- NULL
    cell_fit[[k]] <- fitk
  }

  pred <- rep(global_mean, nte)

  for (k in names(split_test)) {
    idx <- split_test[[k]]

    if (k %in% uniq_keys) {
      if (is.null(cell_fit[[k]]) || length(cont_cols) == 0L) {
        pred[idx] <- cell_mean[[k]]
      } else {
        Xt <- X_test[idx, cont_cols, drop = FALSE]
        pred[idx] <- df_klogit_predict(cell_fit[[k]], Xt)
      }
    } else {
      if (!is.null(global_fit) && length(cont_cols) > 0L) {
        Xt <- X_test[idx, cont_cols, drop = FALSE]
        pred[idx] <- df_klogit_predict(global_fit, Xt)
      } else {
        pred[idx] <- global_mean
      }
    }
  }

  df_clip_prob(pred)
}

# mu(x) = P(Y=1 | X=x) via mixed-type kernel logistic regression
regression_expectation_klogit <- function(dat, new_X, sigma = NULL, success_level = NULL) {
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 & as.numeric(dat$d_np) == 0)
  X_obs <- X_all[idx, , drop = FALSE]
  y_obs <- dat$y[idx]

  if (length(y_obs) < 1L) return(rep(NA_real_, nrow(new_X)))

  bin <- df_as_binary01(y_obs, success_level = success_level)
  y01 <- as.numeric(bin$y01)

  # if only one class, return constant
  u <- unique(y01[is.finite(y01)])
  if (length(u) < 2L) {
    p0 <- df_clip_prob(mean(y01))
    return(rep(p0, nrow(new_X)))
  }

  prep <- df_krr_mixed_prepare(
    X_train = X_obs,
    sigma   = sigma,
    context = "mu(X) = P(Y=1 | X) kernel logistic regression"
  )

  df_klogit_mixed_predict(
    X_train  = X_obs,
    y01_train= y01,
    X_test   = as.matrix(new_X),
    sigma    = prep$sigma,
    cat_cols = prep$cat_cols,
    cont_cols= prep$cont_cols
  )
}

############################################################
## Numeric Jacobian + sandwich
############################################################

df_numeric_jacobian <- function(fun, x0, eps = NULL) {
  x0 <- as.numeric(x0)
  f0 <- as.numeric(fun(x0))
  p  <- length(x0)
  m  <- length(f0)

  if (is.null(eps)) {
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

# Robust joint sandwich:
#  - drop rows with any non-finite score entries
#  - robust inversion (solve with ridge; qr.solve; MASS::ginv)
df_joint_sandwich <- function(score_mat, A_hat) {
  score_mat <- as.matrix(score_mat)
  A_hat <- as.matrix(A_hat)
  q <- ncol(score_mat)

  if (ncol(A_hat) != q || nrow(A_hat) != q) {
    return(list(var = matrix(NA_real_, q, q), A = A_hat, B = matrix(NA_real_, q, q), n_eff = 0L))
  }

  keep <- apply(score_mat, 1, function(r) all(is.finite(r)))
  score_use <- score_mat[keep, , drop = FALSE]
  n <- nrow(score_use)

  if (n <= q) {
    return(list(var = matrix(NA_real_, q, q), A = A_hat, B = matrix(NA_real_, q, q), n_eff = n))
  }

  if (any(!is.finite(A_hat))) {
    return(list(var = matrix(NA_real_, q, q), A = A_hat, B = matrix(NA_real_, q, q), n_eff = n))
  }

  sbar <- colMeans(score_use)
  S <- sweep(score_use, 2, sbar, "-")
  B_hat <- crossprod(S) / n

  A_inv <- NULL
  ridge_seq <- getOption("dfSEDI.sandwich_ridge_seq", c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2))
  ridge_seq <- as.numeric(ridge_seq)
  ridge_seq <- ridge_seq[is.finite(ridge_seq) & ridge_seq >= 0]
  if (length(ridge_seq) == 0L) ridge_seq <- c(0, 1e-8, 1e-6, 1e-4, 1e-2)

  for (lam in ridge_seq) {
    A_try <- A_hat
    if (lam > 0) A_try <- A_try + diag(lam, nrow(A_try))
    A_inv_try <- tryCatch(solve(A_try), error = function(e) NULL)
    if (!is.null(A_inv_try) && all(is.finite(A_inv_try))) {
      A_inv <- A_inv_try
      break
    }
  }

  if (is.null(A_inv)) {
    A_inv <- tryCatch(qr.solve(A_hat), error = function(e) NULL)
    if (!is.null(A_inv) && !all(is.finite(A_inv))) A_inv <- NULL
  }

  if (is.null(A_inv) && requireNamespace("MASS", quietly = TRUE)) {
    A_inv <- MASS::ginv(A_hat)
    if (!all(is.finite(A_inv))) A_inv <- NULL
  }

  V <- if (is.null(A_inv)) {
    matrix(NA_real_, q, q)
  } else {
    A_inv %*% B_hat %*% t(A_inv) / n
  }

  list(var = V, A = A_hat, B = B_hat, n_eff = n)
}

############################################################
## 2. Chang & Kott estimating eq.
############################################################

pi_np.est_simple <- function(dat, base_fun, p.use = TRUE) {

  if (!is.data.frame(dat)) {
    stop("pi_np.est_simple(): `dat` must be a data.frame.", call. = FALSE)
  }
  if (!is.function(base_fun)) {
    stop("pi_np.est_simple(): `base_fun` must be a function.", call. = FALSE)
  }

  # Required columns in dat
  need <- c("d_np")
  if (isTRUE(p.use)) need <- c(need, "d_p", "pi_p")
  if (!all(need %in% names(dat))) {
    miss <- setdiff(need, names(dat))
    stop("pi_np.est_simple(): missing required columns in `dat`: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  # Cache objects that do not depend on phi
  X <- df_get_X(dat)

  # Enforce the convention: X is n x p  => base_fun(X) must be n x (p+2)
  expected_k <- ncol(X) + 2L

  b_mat <- df_eval_base_fun(
    base_fun,
    X,
    expected_ncol = expected_k,
    context = "pi_np.est_simple()/base_fun"
  )

  n <- nrow(b_mat)

  d_np <- as.numeric(dat$d_np)
  if (length(d_np) != n) {
    stop("pi_np.est_simple(): length(d_np) must equal nrow(dat).", call. = FALSE)
  }
  if (any(!is.finite(d_np))) {
    stop("pi_np.est_simple(): `d_np` contains NA/Inf.", call. = FALSE)
  }
  if (any(!(d_np %in% c(0, 1)))) {
    stop("pi_np.est_simple(): `d_np` must be coded as 0/1 (or logical).", call. = FALSE)
  }

  d_p <- NULL
  pi_p <- NULL
  if (isTRUE(p.use)) {
    d_p <- as.numeric(dat$d_p)
    if (length(d_p) != n) {
      stop("pi_np.est_simple(): length(d_p) must equal nrow(dat).", call. = FALSE)
    }
    if (any(!is.finite(d_p))) {
      stop("pi_np.est_simple(): `d_p` contains NA/Inf.", call. = FALSE)
    }
    if (any(!(d_p %in% c(0, 1)))) {
      stop("pi_np.est_simple(): `d_p` must be coded as 0/1 (or logical).", call. = FALSE)
    }

    # pi_p can be NA for (d_p,d_np)=(0,0); fill those to keep algebra finite.
    # For union units (d_p==1 or d_np==1), pi_p must be finite (impute first if needed).
    pi_p <- df_prepare_pi_p(dat$pi_p, d_p = d_p, d_np = d_np, context = "pi_np.est_simple()")
    if (length(pi_p) != n) {
      stop("pi_np.est_simple(): length(pi_p) must equal nrow(dat).", call. = FALSE)
    }
  }

  function(phi) {
    phi <- as.numeric(phi)

    if (length(phi) != ncol(b_mat)) {
      stop(
        sprintf("pi_np.est_simple(): length(phi)=%d but ncol(base_fun(X))=%d.",
                length(phi), ncol(b_mat)),
        call. = FALSE
      )
    }
    if (any(!is.finite(phi))) {
      stop("pi_np.est_simple(): `phi` contains NA/Inf.", call. = FALSE)
    }

    eta   <- as.numeric(b_mat %*% phi)
    pi_np <- df_clip_prob(plogis(eta))

    if (isTRUE(p.use)) {
      denom  <- df_clip_prob(pi_p + pi_np - pi_p * pi_np)
      d_union <- d_p + d_np - d_p * d_np
      d_set4 <- 1 - d_union / denom
    } else {
      d_set4 <- 1 - d_np / pi_np
    }

    as.numeric(crossprod(b_mat, d_set4))
  }
}

############################################################
## 3. Helper functions
############################################################

h4_prob_denom_function <- function(pi_np, pi_p, phi) {
  pi_np <- df_clip_prob(pi_np)
  pi_p  <- df_clip_prob(pi_p)
  pi_np_p <- df_clip_prob(pi_p + pi_np - pi_p * pi_np)
  (1 - pi_np_p) / (pi_np_p ^ 2)
}

h4_prob_numer_function <- function(pi_np, pi_p, phi, y) {
  pi_np <- df_clip_prob(pi_np)
  pi_p  <- df_clip_prob(pi_p)
  pi_np_p <- df_clip_prob(pi_p + pi_np - pi_p * pi_np)
  - y * (1 - pi_np_p) / (pi_np_p ^ 2)
}

eta4_prob_numer_function <- function(pi_np, pi_p, phi, l) {
  pi_np <- df_clip_prob(pi_np)
  pi_p  <- df_clip_prob(pi_p)
  pi_np_p <- df_clip_prob(pi_p + pi_np - pi_p * pi_np)
  - as.numeric(pi_np) * (1 - pi_np_p) / (pi_np_p ^ 2) * as.numeric(l)
}

############################################################
## 4. Mixed-type KRR nuisance regressions (Eff)
############################################################

regression_expectation_kernlab <- function(dat, new_X, sigma = NULL) {
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 & as.numeric(dat$d_np) == 0)
  X_obs <- X_all[idx, , drop = FALSE]
  y_obs <- as.numeric(dat$y[idx])

  if (length(y_obs) < 1L) return(rep(NA_real_, nrow(new_X)))
  if (length(y_obs) < 2L) return(rep(mean(y_obs), nrow(new_X)))

  prep <- df_krr_mixed_prepare(
    X_train = X_obs,
    sigma   = sigma,
    context = "mu(X) = E[Y | X] regression"
  )

  df_krr_mixed_predict(
    X_train  = X_obs,
    y_train  = y_obs,
    X_test   = as.matrix(new_X),
    sigma    = prep$sigma,
    cat_cols = prep$cat_cols,
    cont_cols= prep$cont_cols
  )
}

pi_p_estimation_kernlab <- function(dat, new_L, sigma = NULL) {
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 & is.finite(as.numeric(dat$pi_p)))

  X_obs    <- X_all[idx, , drop = FALSE]
  y_obs    <- as.numeric(dat$y[idx])
  L_obs    <- cbind(X_obs, y_obs)
  pi_p_obs <- df_clip_prob(as.numeric(dat$pi_p[idx]))

  if (length(pi_p_obs) < 1L) return(rep(NA_real_, nrow(new_L)))
  if (length(pi_p_obs) < 2L) return(rep(mean(pi_p_obs), nrow(new_L)))

  cat_cols_X <- df_infer_categorical_cols(X_obs)
  cat_cols_L <- cat_cols_X  # do NOT treat y as categorical

  prep <- df_krr_mixed_prepare(
    X_train = L_obs,
    sigma   = sigma,
    cat_cols = cat_cols_L,
    context = "pi_P(L) regression / imputation"
  )

  reg_pred <- df_krr_mixed_predict(
    X_train  = L_obs,
    y_train  = 1 / pi_p_obs,
    X_test   = as.matrix(new_L),
    sigma    = prep$sigma,
    cat_cols = prep$cat_cols,
    cont_cols= prep$cont_cols
  )

  reg_pred <- as.numeric(reg_pred)
  pi_hat <- rep(NA_real_, length(reg_pred))
  bad <- (!is.finite(reg_pred)) | (reg_pred <= 0)

  pi_hat[bad] <- mean(pi_p_obs)
  pi_hat[!bad] <- 1 / reg_pred[!bad]

  df_clip_prob(pi_hat)
}

# --- NEW: core evaluators that accept pre-extracted (X_obs, l_obs, pi_p_obs, prep) ---
estimate_conditional_expectation_kernlab_phi_core <- function(phi,
                                                              new_X,
                                                              X_obs,
                                                              l_obs,
                                                              pi_p_obs,
                                                              prep) {
  phi <- as.numeric(phi)
  new_X <- as.matrix(new_X)
  X_obs <- as.matrix(X_obs)
  l_obs <- as.matrix(l_obs)
  pi_p_obs <- df_clip_prob(as.numeric(pi_p_obs))

  eta   <- as.numeric(l_obs %*% phi)
  pi_np <- df_clip_prob(1 / (1 + exp(-eta)))

  Denom_vals <- h4_prob_denom_function(pi_np, pi_p_obs, phi)

  # vectorized Numer_mat: Numer_mat[i,] = -pi_np[i] * Denom_vals[i] * l_obs[i,]
  coef <- -pi_np * as.numeric(Denom_vals)
  Numer_mat <- l_obs * coef

  p_dim <- length(phi)

  if (nrow(X_obs) < 2L) {
    denom_pred <- rep(mean(Denom_vals), nrow(new_X))
    denom_pred <- df_floor_pos(denom_pred)
    numer_pred <- matrix(rep(colMeans(Numer_mat), each = nrow(new_X)),
                         nrow = nrow(new_X), ncol = p_dim, byrow = TRUE)
    return(sweep(numer_pred, 1, denom_pred, "/"))
  }

  denom_pred <- df_krr_mixed_predict(
    X_train  = X_obs,
    y_train  = as.numeric(Denom_vals),
    X_test   = new_X,
    sigma    = prep$sigma,
    cat_cols = prep$cat_cols,
    cont_cols= prep$cont_cols
  )
  denom_pred <- df_floor_pos(denom_pred)

  numer_pred <- sapply(
    seq_len(p_dim),
    function(j) df_krr_mixed_predict(
      X_train  = X_obs,
      y_train  = as.numeric(Numer_mat[, j]),
      X_test   = new_X,
      sigma    = prep$sigma,
      cat_cols = prep$cat_cols,
      cont_cols= prep$cont_cols
    )
  )

  sweep(numer_pred, 1, denom_pred, "/")
}

estimate_conditional_expectation_kernlab_theta_core <- function(phi,
                                                                new_X,
                                                                X_obs,
                                                                y_obs,
                                                                l_obs,
                                                                pi_p_obs,
                                                                prep) {
  phi <- as.numeric(phi)
  new_X <- as.matrix(new_X)
  X_obs <- as.matrix(X_obs)
  y_obs <- as.numeric(y_obs)
  l_obs <- as.matrix(l_obs)
  pi_p_obs <- df_clip_prob(as.numeric(pi_p_obs))

  eta   <- as.numeric(l_obs %*% phi)
  pi_np <- df_clip_prob(1 / (1 + exp(-eta)))

  Denom_vals <- h4_prob_denom_function(pi_np, pi_p_obs, phi)
  # vectorized Numer_vals: -y * Denom_vals
  Numer_vals <- -y_obs * as.numeric(Denom_vals)

  if (nrow(X_obs) < 2L) {
    denom_pred <- rep(mean(Denom_vals), nrow(new_X))
    denom_pred <- df_floor_pos(denom_pred)
    numer_pred <- rep(mean(Numer_vals), nrow(new_X))
    return(as.numeric(numer_pred / denom_pred))
  }

  denom_pred <- df_krr_mixed_predict(
    X_train  = X_obs,
    y_train  = as.numeric(Denom_vals),
    X_test   = new_X,
    sigma    = prep$sigma,
    cat_cols = prep$cat_cols,
    cont_cols= prep$cont_cols
  )
  denom_pred <- df_floor_pos(denom_pred)

  numer_pred <- df_krr_mixed_predict(
    X_train  = X_obs,
    y_train  = as.numeric(Numer_vals),
    X_test   = new_X,
    sigma    = prep$sigma,
    cat_cols = prep$cat_cols,
    cont_cols= prep$cont_cols
  )

  as.numeric(numer_pred / denom_pred)
}

# Backward-compatible wrappers (same outputs as before)
estimate_conditional_expectation_kernlab_phi <- function(dat, phi, new_X, sigma = NULL, prep = NULL) {
  phi <- as.numeric(phi)
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 | as.numeric(dat$d_np) == 1)

  X_obs    <- X_all[idx, , drop = FALSE]
  y_obs    <- as.numeric(dat$y[idx])
  pi_p_obs <- df_clip_prob(as.numeric(dat$pi_p[idx]))

  l_obs  <- as.matrix(cbind(1, X_obs, y_obs))
  storage.mode(l_obs) <- "numeric"

  if (is.null(prep)) {
    prep <- df_krr_mixed_prepare(
      X_train = X_obs,
      sigma   = sigma,
      context = "eta4*(X;phi) nuisance regression"
    )
  } else {
    # keep cat/cont from prep; override sigma if explicitly provided
    if (!is.null(sigma)) prep$sigma <- sigma
  }

  estimate_conditional_expectation_kernlab_phi_core(
    phi     = phi,
    new_X   = as.matrix(new_X),
    X_obs   = X_obs,
    l_obs   = l_obs,
    pi_p_obs= pi_p_obs,
    prep    = prep
  )
}

estimate_conditional_expectation_kernlab_theta <- function(dat, phi, new_X, sigma = NULL, prep = NULL) {
  phi <- as.numeric(phi)
  X_all <- df_get_X(dat)
  idx   <- which(as.numeric(dat$d_p) == 1 | as.numeric(dat$d_np) == 1)

  X_obs    <- X_all[idx, , drop = FALSE]
  y_obs    <- as.numeric(dat$y[idx])
  pi_p_obs <- df_clip_prob(as.numeric(dat$pi_p[idx]))

  l_obs  <- as.matrix(cbind(1, X_obs, y_obs))
  storage.mode(l_obs) <- "numeric"

  if (is.null(prep)) {
    prep <- df_krr_mixed_prepare(
      X_train = X_obs,
      sigma   = sigma,
      context = "h4*(X;phi) nuisance regression"
    )
  } else {
    if (!is.null(sigma)) prep$sigma <- sigma
  }

  estimate_conditional_expectation_kernlab_theta_core(
    phi     = phi,
    new_X   = as.matrix(new_X),
    X_obs   = X_obs,
    y_obs   = y_obs,
    l_obs   = l_obs,
    pi_p_obs= pi_p_obs,
    prep    = prep
  )
}

############################################################
## 5. Scores
############################################################

# Original versions (kept for compatibility)
df_score_phi_contrib <- function(dat1, phi, eta4_star_local) {
  phi <- as.numeric(phi)

  X1    <- df_get_X(dat1)
  y1    <- as.numeric(dat1$y)
  d_p1  <- as.numeric(dat1$d_p)
  d_np1 <- as.numeric(dat1$d_np)

  # pi_p may be NA for (d_p,d_np)=(0,0); fill those. Union rows must be finite.
  pi_p1 <- df_prepare_pi_p(dat1$pi_p, d_p = d_p1, d_np = d_np1, context = "df_score_phi_contrib()")

  l1 <- as.matrix(cbind(1, X1, y1))

  eta      <- as.numeric(l1 %*% phi)
  pi_np1   <- df_clip_prob(1 / (1 + exp(-eta)))
  pi_np_p1 <- df_clip_prob(pi_np1 + pi_p1 - pi_np1 * pi_p1)

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

  out <- term1 + term2_coef * eta4_star_local
  out[!is.finite(out)] <- NA_real_
  out
}

efficient_theta_contrib <- function(dat1, phi, h4_star_local) {
  phi <- as.numeric(phi)

  X1    <- df_get_X(dat1)
  y1    <- as.numeric(dat1$y)
  d_p1  <- as.numeric(dat1$d_p)
  d_np1 <- as.numeric(dat1$d_np)

  # pi_p may be NA for (d_p,d_np)=(0,0); fill those. Union rows must be finite.
  pi_p1 <- df_prepare_pi_p(dat1$pi_p, d_p = d_p1, d_np = d_np1, context = "efficient_theta_contrib()")

  l1 <- as.matrix(cbind(1, X1, y1))

  eta      <- as.numeric(l1 %*% phi)
  pi_np1   <- df_clip_prob(1 / (1 + exp(-eta)))
  pi_np_p1 <- df_clip_prob(pi_np1 + pi_p1 - pi_np1 * pi_p1)

  h4_star_local <- as.numeric(h4_star_local)

  term1 <- y1 * (d_p1 / pi_p1 +
                   (1 - d_p1 / pi_p1) / pi_np_p1 *
                   (d_np1 - d_p1 * (d_np1 - pi_np1)))

  term2_coef <- 1 -
    (d_p1 * d_np1) / (pi_p1 * pi_np1) -
    1 / pi_np_p1 *
    (d_np1 * (1 - d_p1 / pi_p1) +
       d_p1 * (1 - d_np1 / pi_np1))

  out <- as.numeric(term1 - term2_coef * h4_star_local)
  out[!is.finite(out)] <- NA_real_
  out
}

# Fast versions (output-identical) used internally for DML objective/jacobians
df_score_phi_contrib_fast <- function(l1, d_p1, d_np1, pi_p1, phi, eta4_star_local) {
  phi <- as.numeric(phi)
  l1 <- as.matrix(l1)
  d_p1 <- as.numeric(d_p1)
  d_np1 <- as.numeric(d_np1)

  # pi_p may be NA for (d_p,d_np)=(0,0); fill those. Union rows must be finite.
  pi_p1 <- df_prepare_pi_p(pi_p1, d_p = d_p1, d_np = d_np1, context = "df_score_phi_contrib_fast()")

  eta      <- as.numeric(l1 %*% phi)
  pi_np1   <- df_clip_prob(1 / (1 + exp(-eta)))
  pi_np_p1 <- df_clip_prob(pi_np1 + pi_p1 - pi_np1 * pi_p1)

  term1 <- l1 * pi_np1 * (1 - pi_np1) / pi_np_p1 *
    (d_np1 / pi_np1 * (pi_p1 - d_p1) -
       d_p1 / (1 - pi_np1) * (1 - d_np1 / pi_np1))

  term2_coef <- 1 -
    (d_p1 * d_np1) / (pi_p1 * pi_np1) -
    1 / pi_np_p1 *
    (d_np1 * (1 - d_p1 / pi_p1) +
       d_p1 * (1 - d_np1 / pi_np1))

  eta4_star_local <- as.matrix(eta4_star_local)
  out <- term1 + term2_coef * eta4_star_local
  out[!is.finite(out)] <- NA_real_
  out
}

efficient_theta_contrib_fast <- function(l1, y1, d_p1, d_np1, pi_p1, phi, h4_star_local) {
  phi <- as.numeric(phi)
  l1 <- as.matrix(l1)
  y1 <- as.numeric(y1)
  d_p1 <- as.numeric(d_p1)
  d_np1 <- as.numeric(d_np1)

  # pi_p may be NA for (d_p,d_np)=(0,0); fill those. Union rows must be finite.
  pi_p1 <- df_prepare_pi_p(pi_p1, d_p = d_p1, d_np = d_np1, context = "efficient_theta_contrib_fast()")

  h4_star_local <- as.numeric(h4_star_local)

  eta      <- as.numeric(l1 %*% phi)
  pi_np1   <- df_clip_prob(1 / (1 + exp(-eta)))
  pi_np_p1 <- df_clip_prob(pi_np1 + pi_p1 - pi_np1 * pi_p1)

  term1 <- y1 * (d_p1 / pi_p1 +
                   (1 - d_p1 / pi_p1) / pi_np_p1 *
                   (d_np1 - d_p1 * (d_np1 - pi_np1)))

  term2_coef <- 1 -
    (d_p1 * d_np1) / (pi_p1 * pi_np1) -
    1 / pi_np_p1 *
    (d_np1 * (1 - d_p1 / pi_p1) +
       d_p1 * (1 - d_np1 / pi_np1))

  out <- as.numeric(term1 - term2_coef * h4_star_local)
  out[!is.finite(out)] <- NA_real_
  out
}

############################################################
## 6. DML helpers
############################################################

make_folds <- function(n, K) {
  fold_id <- sample(rep(1:K, length.out = n))
  split(seq_len(n), fold_id)
}

# Simple (non-crossfit) pi_p imputation for np-only union units.
# This is used by estimators that require pi_p on the union sample but are not cross-fitted.
impute_pi_p_simple <- function(dat, sigma = NULL) {
  if (!is.data.frame(dat)) stop("impute_pi_p_simple(): `dat` must be a data.frame.", call. = FALSE)
  if (!all(c("d_np", "d_p", "pi_p", "y") %in% names(dat))) {
    miss <- setdiff(c("d_np", "d_p", "pi_p", "y"), names(dat))
    stop("impute_pi_p_simple(): missing required columns in `dat`: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  d_np <- as.numeric(dat$d_np)
  d_p  <- as.numeric(dat$d_p)
  pi_p_raw <- as.numeric(dat$pi_p)

  idx_mis <- which(d_np == 1 & d_p == 0 & !is.finite(pi_p_raw))
  if (length(idx_mis) == 0L) return(dat)

  X_all <- df_get_X(dat)
  y_all <- as.numeric(dat$y)
  L_all <- as.matrix(cbind(X_all, y_all))

  tilde_pi <- pi_p_estimation_kernlab(dat, new_L = L_all[idx_mis, , drop = FALSE], sigma = sigma)
  dat$pi_p[idx_mis] <- tilde_pi
  dat
}

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
## 7. Basic estimators
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
  se_hat    <- sqrt(pmax(var_hat, 0))

  z  <- stats::qnorm(1 - (1 - level) / 2)
  ci <- c(theta_hat - z * se_hat, theta_hat + z * se_hat)

  list(theta = theta_hat, var = var_hat, se = se_hat, ci = ci, n_eff = n)
}

# --- NEW: population size + Hajek helpers (used only when N is supplied or hajek=TRUE) ---

df_parse_population_N <- function(N) {
  if (is.null(N)) return(NULL)
  if (length(N) == 0L) return(NULL)

  if (is.character(N)) {
    N0 <- gsub("[,[:space:]]", "", N[1])
    N <- suppressWarnings(as.numeric(N0))
  } else {
    N <- as.numeric(N)[1]
  }

  if (!is.finite(N) || N <= 0) {
    stop("Population size `N` must be a positive finite number.", call. = FALSE)
  }

  if (abs(N - round(N)) > 1e-8) {
    warning("`N` was not an integer; rounding to the nearest integer.", call. = FALSE)
  }

  as.numeric(round(N))
}

df_dat_has_frame_units <- function(dat) {
  if (!is.data.frame(dat)) return(FALSE)
  if (!all(c("d_p", "d_np") %in% names(dat))) return(FALSE)
  d_p  <- as.numeric(dat$d_p)
  d_np <- as.numeric(dat$d_np)
  any(d_p == 0 & d_np == 0, na.rm = TRUE)
}

df_resolve_population_N <- function(dat, N = NULL, require = FALSE, context = "estimator") {
  N0 <- df_parse_population_N(N)
  if (!is.null(N0)) {
    if (N0 < nrow(dat)) {
      stop(context, ": supplied `N` (", N0, ") is smaller than nrow(dat) (", nrow(dat), ").", call. = FALSE)
    }
    return(N0)
  }

  if (df_dat_has_frame_units(dat)) {
    return(as.numeric(nrow(dat)))
  }

  if (isTRUE(require)) {
    stop(
      context, ": full frame data were not provided (no (d_p,d_np)=(0,0) rows detected).\n",
      "Please supply the population size via `N=...` (e.g., N=10000), or set `hajek=TRUE`.",
      call. = FALSE
    )
  }

  NULL
}

# Sandwich like df_sandwich_from_contrib(), but assumes `n_missing` additional 0-contributions.
df_sandwich_from_contrib_with_zeros <- function(contrib, n_missing = 0L, level = 0.95) {
  contrib <- as.numeric(contrib)
  contrib <- contrib[is.finite(contrib)]
  n_obs <- length(contrib)

  n_missing <- as.integer(n_missing)
  if (!is.finite(n_missing) || n_missing < 0L) n_missing <- 0L

  n_eff <- n_obs + n_missing
  if (n_eff <= 1L) {
    return(list(theta = NA_real_, var = NA_real_, se = NA_real_, ci = c(NA_real_, NA_real_), n_eff = n_eff))
  }

  sum_c   <- sum(contrib)
  sumsq_c <- sum(contrib^2)

  theta_hat <- sum_c / n_eff

  # sample variance of vector length n_eff: [contrib, 0, ..., 0]
  s2 <- (sumsq_c - n_eff * theta_hat^2) / (n_eff - 1)
  s2 <- pmax(s2, 0)

  var_hat <- s2 / n_eff
  se_hat  <- sqrt(pmax(var_hat, 0))

  z  <- stats::qnorm(1 - (1 - level) / 2)
  ci <- c(theta_hat - z * se_hat, theta_hat + z * se_hat)

  list(theta = theta_hat, var = var_hat, se = se_hat, ci = ci, n_eff = n_eff)
}

# Hajek (ratio) variance proxy for weighted mean: sum(w*y)/sum(w)
df_hajek_mean_with_se <- function(y, w, level = 0.95) {
  y <- as.numeric(y)
  w <- as.numeric(w)

  keep <- is.finite(y) & is.finite(w) & (w > 0)
  y <- y[keep]
  w <- w[keep]

  if (length(y) <= 1L) {
    return(list(theta = NA_real_, var = NA_real_, se = NA_real_, ci = c(NA_real_, NA_real_), n_eff = length(y)))
  }

  den <- sum(w)
  if (!is.finite(den) || den <= 0) {
    return(list(theta = NA_real_, var = NA_real_, se = NA_real_, ci = c(NA_real_, NA_real_), n_eff = length(y)))
  }

  theta_hat <- sum(w * y) / den

  v <- sum((w^2) * (y - theta_hat)^2) / (den^2)
  v <- pmax(v, 0)
  se <- sqrt(v)

  z  <- stats::qnorm(1 - (1 - level) / 2)
  ci <- c(theta_hat - z * se, theta_hat + z * se)

  list(theta = theta_hat, var = v, se = se, ci = ci, n_eff = length(y))
}

df_estimate_P <- function(dat, N = NULL, hajek = FALSE) {
  d_p  <- as.numeric(dat$d_p)
  y    <- as.numeric(dat$y)
  pi_p_raw <- as.numeric(dat$pi_p)

  if (isTRUE(hajek)) {
    idx <- which(d_p == 1)
    pi_use <- df_clip_prob(pi_p_raw[idx])
    if (any(!is.finite(pi_use))) {
      stop("df_estimate_P(): `pi_p` is missing (NA/Inf) for some d_p==1 units.", call. = FALSE)
    }
    w <- 1 / pi_use
    hj <- df_hajek_mean_with_se(y = y[idx], w = w)
    return(list(theta = hj$theta, var = hj$var, se = hj$se, ci = hj$ci))
  }

  N_total <- df_resolve_population_N(dat, N = N, require = TRUE, context = "df_estimate_P()")
  n_missing <- as.integer(round(N_total - nrow(dat)))
  if (n_missing < 0L) {
    stop("df_estimate_P(): N < nrow(dat).", call. = FALSE)
  }

  contrib <- rep(0, length(y))
  idx <- which(d_p == 1)
  if (length(idx) > 0L) {
    pi_use <- df_clip_prob(pi_p_raw[idx])
    if (any(!is.finite(pi_use))) {
      stop("df_estimate_P(): `pi_p` is missing (NA/Inf) for some d_p==1 units.", call. = FALSE)
    }
    contrib[idx] <- y[idx] / pi_use
  }

  theta_res <- if (n_missing == 0L) {
    df_sandwich_from_contrib(contrib)
  } else {
    df_sandwich_from_contrib_with_zeros(contrib, n_missing = n_missing)
  }

  list(theta = theta_res$theta, var = theta_res$var, se = theta_res$se, ci = theta_res$ci)
}

df_estimate_NP <- function(dat,
                           base_fun,
                           phi_start = NULL,
                           max_iter = 20,
                           N = NULL,
                           hajek = FALSE) {

  # Validate base_fun
  if (!is.function(base_fun)) stop("df_estimate_NP(): base_fun must be a function.")

  X <- df_get_X(dat)
  expected_k <- ncol(X) + 2L

  # Evaluate base_fun once (also validates dimensions)
  b_mat <- df_eval_base_fun(base_fun, X, context = "df_estimate_NP")
  if (ncol(b_mat) != expected_k) {
    stop(sprintf("df_estimate_NP(): base_fun(X) must return n x (p+2) matrix. Expected %d columns, got %d.",
                 expected_k, ncol(b_mat)))
  }
  p_phi <- ncol(b_mat)

  d_np <- as.numeric(dat$d_np)
  y <- as.numeric(dat$y)
  l_mat <- df_make_l_matrix(X, y)

  # Center for initial values
  if (is.null(phi_start)) {
    p0 <- df_clip_prob(mean(d_np, na.rm = TRUE))
    phi_center <- c(stats::qlogis(p0), rep(0, p_phi - 1L))
  } else {
    phi_center <- as.numeric(phi_start)
    if (length(phi_center) != p_phi) stop("df_estimate_NP(): phi_start length mismatch.")
    if (any(!is.finite(phi_center))) stop("df_estimate_NP(): phi_start contains NA/Inf.")
  }

  # Estimating equation (p.use = FALSE)
  ee_phi <- function(phi) {
    eta <- as.numeric(l_mat %*% phi)
    pi_np <- df_clip_prob(stats::plogis(eta))
    d_set4 <- 1 - d_np / pi_np
    as.numeric(crossprod(b_mat, d_set4))
  }

  # Objective: squared norm of estimating equations
  big_penalty <- 1e100
  obj_phi <- function(phi) {
    eq <- tryCatch(ee_phi(phi), error = function(e) rep(NA_real_, p_phi))
    if (any(!is.finite(eq))) return(big_penalty)
    sum(eq^2)
  }

  # Generate 100 initial candidates around phi_center and try optim() up to 10 times
  n_try <- min(10L, as.integer(max_iter))
  ms <- df_multistart_optim(obj_fun = obj_phi,
                            center = phi_center,
                            n_candidates = 100L,
                            n_try = n_try,
                            width = 0.5,
                            method = "Nelder-Mead",
                            control = list(maxit = 500))

  if (!isTRUE(ms$ok)) {
    phi_hat <- rep(NA_real_, p_phi)
    theta_res <- df_sandwich_from_contrib(rep(NA_real_, nrow(dat)))
    return(list(phi = phi_hat,
                theta = theta_res$theta,
                var = theta_res$var,
                se = theta_res$se,
                ci = theta_res$ci))
  }

  phi_hat <- ms$res$par

  # Compute pi_np under phi_hat
  eta <- as.numeric(l_mat %*% phi_hat)
  pi_np <- df_clip_prob(stats::plogis(eta))

  # Theta
  if (isTRUE(hajek)) {
    idx <- which(d_np == 1)
    if (length(idx) == 0) {
      theta_res <- list(theta = NA_real_, var = NA_real_, se = NA_real_, ci = c(NA_real_, NA_real_))
    } else {
      w <- 1 / pi_np[idx]
      theta_res <- df_hajek_mean_with_se(y[idx], w)
    }
  } else {
    N_total <- df_resolve_population_N(dat, N = N, require = TRUE)
    n_missing <- N_total - nrow(dat)

    contrib <- rep(0, nrow(dat))
    idx <- which(d_np == 1)
    if (length(idx) > 0) contrib[idx] <- y[idx] / pi_np[idx]

    if (n_missing > 0) {
      theta_res <- df_sandwich_from_contrib_with_zeros(contrib, n_missing = n_missing)
    } else {
      theta_res <- df_sandwich_from_contrib(contrib)
    }
  }

  list(phi = phi_hat,
       theta = theta_res$theta,
       var = theta_res$var,
       se = theta_res$se,
       ci = theta_res$ci)
}

df_estimate_NP_P <- function(dat,
                             base_fun,
                             phi_start = NULL,
                             max_iter = 20,
                             N = NULL,
                             hajek = FALSE) {

  # Validate base_fun
  if (!is.function(base_fun)) stop("df_estimate_NP_P(): base_fun must be a function.")

  X <- df_get_X(dat)
  expected_k <- ncol(X) + 2L

  # Evaluate base_fun once (also validates dimensions)
  b_mat <- df_eval_base_fun(base_fun, X, context = "df_estimate_NP_P")
  if (ncol(b_mat) != expected_k) {
    stop(sprintf("df_estimate_NP_P(): base_fun(X) must return n x (p+2) matrix. Expected %d columns, got %d.",
                 expected_k, ncol(b_mat)))
  }
  p_phi <- ncol(b_mat)

  d_np <- as.numeric(dat$d_np)
  d_p <- as.numeric(dat$d_p)
  y <- as.numeric(dat$y)
  l_mat <- df_make_l_matrix(X, y)
  # If pi_p is missing for (d_np, d_p) = (1, 0), impute it via a parametric lm on 1/pi_p.
  dat <- df_impute_pi_p_lm(dat, context = "df_estimate_NP_P()")


  # pi_p can be missing for d_p == 0. We only require it for union-sample units.
  pi_p <- df_prepare_pi_p(dat)

  # Center for initial values
  if (is.null(phi_start)) {
    p0 <- df_clip_prob(mean(d_np, na.rm = TRUE))
    phi_center <- c(stats::qlogis(p0), rep(0, p_phi - 1L))
  } else {
    phi_center <- as.numeric(phi_start)
    if (length(phi_center) != p_phi) stop("df_estimate_NP_P(): phi_start length mismatch.")
    if (any(!is.finite(phi_center))) stop("df_estimate_NP_P(): phi_start contains NA/Inf.")
  }

  # Estimating equation (p.use = TRUE)
  ee_phi <- function(phi) {
    eta <- as.numeric(l_mat %*% phi)
    pi_np <- df_clip_prob(stats::plogis(eta))

    denom <- df_clip_prob(pi_p + pi_np - pi_p * pi_np)
    d_union <- d_p + d_np - d_p * d_np

    d_set4 <- 1 - d_union / denom
    as.numeric(crossprod(b_mat, d_set4))
  }

  # Objective: squared norm of estimating equations
  big_penalty <- 1e100
  obj_phi <- function(phi) {
    eq <- tryCatch(ee_phi(phi), error = function(e) rep(NA_real_, p_phi))
    if (any(!is.finite(eq))) return(big_penalty)
    sum(eq^2)
  }

  # Generate 100 initial candidates around phi_center and try optim() up to 10 times
  n_try <- min(10L, as.integer(max_iter))
  ms <- df_multistart_optim(obj_fun = obj_phi,
                            center = phi_center,
                            n_candidates = 100L,
                            n_try = n_try,
                            width = 0.5,
                            method = "Nelder-Mead",
                            control = list(maxit = 500))

  if (!isTRUE(ms$ok)) {
    phi_hat <- rep(NA_real_, p_phi)
    theta_res <- df_sandwich_from_contrib(rep(NA_real_, nrow(dat)))
    return(list(phi = phi_hat,
                theta = theta_res$theta,
                var = theta_res$var,
                se = theta_res$se,
                ci = theta_res$ci))
  }

  phi_hat <- ms$res$par

  # Compute pi_np under phi_hat
  eta <- as.numeric(l_mat %*% phi_hat)
  pi_np <- df_clip_prob(stats::plogis(eta))

  # Union denominator and indicator
  denom <- df_clip_prob(pi_p + pi_np - pi_p * pi_np)
  d_union <- d_p + d_np - d_p * d_np

  # Theta
  if (isTRUE(hajek)) {
    idx <- which(d_union == 1)
    if (length(idx) == 0) {
      theta_res <- list(theta = NA_real_, var = NA_real_, se = NA_real_, ci = c(NA_real_, NA_real_))
    } else {
      w <- 1 / denom[idx]
      theta_res <- df_hajek_mean_with_se(y[idx], w)
    }
  } else {
    N_total <- df_resolve_population_N(dat, N = N, require = TRUE)
    n_missing <- N_total - nrow(dat)

    contrib <- rep(0, nrow(dat))
    idx <- which(d_union == 1)
    if (length(idx) > 0) contrib[idx] <- y[idx] / denom[idx]

    if (n_missing > 0) {
      theta_res <- df_sandwich_from_contrib_with_zeros(contrib, n_missing = n_missing)
    } else {
      theta_res <- df_sandwich_from_contrib(contrib)
    }
  }

  list(phi = phi_hat,
       theta = theta_res$theta,
       var = theta_res$var,
       se = theta_res$se,
       ci = theta_res$ci)
}

############################################################
## 8. Efficient estimator Eff (DML2)
############################################################

efficient_estimator_dml2 <- function(dat,
                                     phi_start   = NULL,
                                     K           = 2,
                                     max_restart = 10,
                                     progress    = FALSE,
                                     x_info      = TRUE,
                                     N           = NULL) {

  had_frame_units <- df_dat_has_frame_units(dat)

  N_total <- NULL
  if (!isTRUE(x_info)) {
    N_total <- df_resolve_population_N(dat, N = N, require = FALSE, context = "Eff()/x_info=FALSE")
  }

  dat <- df_apply_x_info(dat, x_info)

  n <- nrow(dat)
  X <- df_get_X(dat)
  p_x <- ncol(X)
  p_phi <- 1 + p_x + 1

  theta_scale <- 1
  if (!isTRUE(x_info)) {
    if (is.null(N_total) && !isTRUE(had_frame_units)) {
      warning(
        "Eff(): x_info=FALSE but population size N was not provided and no (d_p,d_np)=(0,0) rows were detected.\n",
        "Returning the mean over the provided data (not the full-frame mean).\n",
        "To target the full-frame mean, call Eff(..., x_info=FALSE, N=<population size>).",
        call. = FALSE
      )
    }
    if (!is.null(N_total) && is.finite(N_total) && N_total > 0) {
      theta_scale <- n / N_total
    }
  }

  if (is.null(phi_start)) {
    p0 <- df_clip_prob(mean(as.numeric(dat$d_np), na.rm = TRUE))
    phi_start <- c(stats::qlogis(p0), rep(0, p_x), 0)
  }

  folds <- make_folds(n, K)
  dat_cf <- impute_pi_p_crossfit(dat, folds, sigma = NULL, progress = progress)

  # Build fold caches (output-preserving): precompute X_test/l_test and train prep for eta4*/h4*
  make_fold_cache <- function(idx_test, idx_train) {
    dat_test  <- dat_cf[idx_test,  , drop = FALSE]
    dat_train <- dat_cf[idx_train, , drop = FALSE]

    # test side (for fast scores)
    X_test <- df_get_X(dat_test)
    y_test <- as.numeric(dat_test$y)
    l_test <- df_make_l_matrix(X_test, y_test)

    test <- list(
      l    = l_test,
      y    = y_test,
      d_p  = as.numeric(dat_test$d_p),
      d_np = as.numeric(dat_test$d_np),
      pi_p = as.numeric(dat_test$pi_p)  # allow NA for (0,0); handled in df_prepare_pi_p()
    )

    if (!isTRUE(x_info)) {
      return(list(
        idx_test = idx_test,
        w = length(idx_test) / n,
        X_test = X_test,
        test = test,
        train = NULL
      ))
    }

    # train side for nuisances
    X_all_tr <- df_get_X(dat_train)
    idx_obs  <- which(as.numeric(dat_train$d_p) == 1 | as.numeric(dat_train$d_np) == 1)

    X_obs    <- X_all_tr[idx_obs, , drop = FALSE]
    y_obs    <- as.numeric(dat_train$y[idx_obs])
    pi_p_obs <- df_clip_prob(as.numeric(dat_train$pi_p[idx_obs]))
    l_obs    <- df_make_l_matrix(X_obs, y_obs)

    prep_eta4 <- df_krr_mixed_prepare(
      X_train = X_obs,
      sigma   = NULL,
      context = "eta4*(X;phi) nuisance regression"
    )
    prep_h4 <- df_krr_mixed_prepare(
      X_train = X_obs,
      sigma   = NULL,
      context = "h4*(X;phi) nuisance regression"
    )

    train <- list(
      X_obs     = X_obs,
      y_obs     = y_obs,
      pi_p_obs  = pi_p_obs,
      l_obs     = l_obs,
      prep_eta4 = prep_eta4,
      prep_h4   = prep_h4
    )

    list(
      idx_test = idx_test,
      w = length(idx_test) / n,
      X_test = X_test,
      test = test,
      train = train
    )
  }

  fold_cache <- vector("list", length(folds))
  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(n), idx_test)
    fold_cache[[k]] <- make_fold_cache(idx_test, idx_train)
  }

  big_penalty <- 1e100

  obj_phi <- function(phi) {
    phi <- as.numeric(phi)
    eq_agg <- rep(0, length(phi))

    for (k in seq_along(fold_cache)) {
      fc <- fold_cache[[k]]

      eta4_k <- if (isTRUE(x_info)) {
        estimate_conditional_expectation_kernlab_phi_core(
          phi      = phi,
          new_X    = fc$X_test,
          X_obs    = fc$train$X_obs,
          l_obs    = fc$train$l_obs,
          pi_p_obs = fc$train$pi_p_obs,
          prep     = fc$train$prep_eta4
        )
      } else {
        matrix(0, nrow = nrow(fc$test$l), ncol = p_phi)
      }

      sphi_k <- df_score_phi_contrib_fast(
        l1 = fc$test$l,
        d_p1 = fc$test$d_p,
        d_np1 = fc$test$d_np,
        pi_p1 = fc$test$pi_p,
        phi = phi,
        eta4_star_local = eta4_k
      )

      eq_k <- colMeans(sphi_k, na.rm = TRUE)
      if (any(!is.finite(eq_k))) return(big_penalty)

      eq_agg <- eq_agg + fc$w * eq_k
    }

    if (any(!is.finite(eq_agg))) return(big_penalty)
    sum(eq_agg^2)
  }

  if (progress) {
    cat("Step 2/3: solving phi (Eff, DML2, Nelder-Mead) ...\n")
    flush.console()
  }

  # --- Multi-start phi estimation (100 candidates -> try up to 10) ---
  tol_obj   <- getOption("dfSEDI.phi_obj_tol", 1e-4)
  maxit_phi <- as.integer(getOption("dfSEDI.phi_maxit", 2000L))
  if (!is.finite(tol_obj) || tol_obj <= 0) tol_obj <- 1e-4
  if (!is.finite(maxit_phi) || maxit_phi < 200L) maxit_phi <- 2000L

  n_cand <- 100L
  cand_mat <- matrix(rep(phi_start, each = n_cand), nrow = n_cand) +
    matrix(stats::runif(n_cand * p_phi, -0.5, 0.5), nrow = n_cand)
  obj_vals <- apply(cand_mat, 1L, obj_phi)
  ord <- order(obj_vals)

  max_try <- min(10L, as.integer(max_restart))
  if (!is.finite(max_try) || max_try < 1L) max_try <- 10L
  max_try <- min(max_try, 10L)

  res <- NULL
  success <- FALSE
  best <- NULL
  res  <- NULL
  for (attempt in seq_len(max_try)) {
    init <- as.numeric(cand_mat[ord[attempt], ])
    res_try <- tryCatch(
      stats::optim(
        par     = init,
        fn      = obj_phi,
        method  = "Nelder-Mead",
        control = list(maxit = maxit_phi)
      ),
      error = function(e) NULL
    )

    if (is.null(res_try) || !is.finite(res_try$value)) {
      if (progress) {
        cat(sprintf("  attempt %d/%d: failed (error/non-finite)\n", attempt, max_try))
      }
      next
    }

    if (progress) {
      cat(sprintf("  attempt %d/%d: convergence=%d, obj=%.6g\n",
                  attempt, max_try, res_try$convergence, res_try$value))
    }

    if (is.null(best) || res_try$value < best$value) {
      best <- res_try
    }

    if (isTRUE(res_try$convergence == 0) || res_try$value <= tol_obj) {
      res <- res_try
      success <- TRUE
      break
    }
  }

  if (!isTRUE(success) && !is.null(best) && is.finite(best$value) && best$value <= tol_obj) {
    res <- best
    success <- TRUE
  }

  if (!isTRUE(success) || is.null(res) || !is.finite(res$value) || res$value >= 0.99 * big_penalty) {
    if (progress) {
      cat("  phi optimization did not reach the tolerance; returning NA.\n")
    } else {
      warning("Eff(): phi optimization did not reach the tolerance; returning NA.")
    }

    return(list(
      phi      = rep(NA_real_, p_phi),
      phi_var  = matrix(NA_real_, nrow = p_phi, ncol = p_phi),
      phi_se   = rep(NA_real_, p_phi),
      phi_ci   = matrix(
        NA_real_,
        nrow = 2L,
        ncol = p_phi,
        dimnames = list(c("l", "u"), paste0("phi_", seq_len(p_phi)))
      ),
      theta    = NA_real_,
      var      = NA_real_,
      se       = NA_real_,
      ci       = c(NA_real_, NA_real_),
      info     = list(
        type            = "Eff",
        dml_type         = "DML2",
        K                = K,
        phi_start        = phi_start,
        max_restart      = max_restart,
        progress         = progress,
        x_info           = isTRUE(x_info),
        aug_terms        = if (isTRUE(x_info)) "estimated" else "fixed_zero",
        convergence      = if (!is.null(best)) best$convergence else NA_integer_,
        phi_obj          = if (!is.null(best)) best$value else NA_real_,
        theta_var_method = theta_var_method,
        sandwich_n_eff   = n,
        population_N     = N_total,
        n_obs            = n,
        theta_scale      = theta_scale
      )
    ))
  }

  phi_hat <- as.numeric(res$par)

  eta4_all <- matrix(NA_real_, nrow = n, ncol = p_phi)
  h4_all   <- rep(NA_real_, n)

  if (isTRUE(x_info)) {
    if (progress) {
      cat("Step 3/3: cross-fitting nuisances for joint scores ...\n")
      pb <- utils::txtProgressBar(min = 0, max = length(fold_cache), style = 3)
    }

    for (k in seq_along(fold_cache)) {
      fc <- fold_cache[[k]]
      idx_test <- fc$idx_test

      eta4_k <- estimate_conditional_expectation_kernlab_phi_core(
        phi      = phi_hat,
        new_X    = fc$X_test,
        X_obs    = fc$train$X_obs,
        l_obs    = fc$train$l_obs,
        pi_p_obs = fc$train$pi_p_obs,
        prep     = fc$train$prep_eta4
      )
      h4_k <- estimate_conditional_expectation_kernlab_theta_core(
        phi      = phi_hat,
        new_X    = fc$X_test,
        X_obs    = fc$train$X_obs,
        y_obs    = fc$train$y_obs,
        l_obs    = fc$train$l_obs,
        pi_p_obs = fc$train$pi_p_obs,
        prep     = fc$train$prep_h4
      )

      eta4_all[idx_test, ] <- eta4_k
      h4_all[idx_test]     <- h4_k

      if (progress) utils::setTxtProgressBar(pb, k)
    }

    if (progress) {
      close(pb)
      cat("\n")
    }
  } else {
    # x_info=FALSE: set augmentation terms to 0 and skip estimation (fast path)
    eta4_all[,] <- 0
    h4_all[] <- 0
    if (progress) {
      cat("Step 3/3: x_info=FALSE -> skip eta4*/h4* estimation (set to 0).\n")
      flush.console()
    }
  }

  # Precompute l for full sample for score/jacobian (fast path; output-identical)
  X_cf <- df_get_X(dat_cf)
  y_cf <- as.numeric(dat_cf$y)
  l_cf <- df_make_l_matrix(X_cf, y_cf)
  d_p_cf  <- as.numeric(dat_cf$d_p)
  d_np_cf <- as.numeric(dat_cf$d_np)
  pi_p_cf <- as.numeric(dat_cf$pi_p)  # allow NA for (0,0); handled in df_prepare_pi_p()

  contrib_theta <- efficient_theta_contrib_fast(l_cf, y_cf, d_p_cf, d_np_cf, pi_p_cf, phi_hat, h4_all)

  # Unscaled (n-based) mean
  theta_hat_n     <- mean(contrib_theta, na.rm = TRUE)
  theta_hat       <- theta_hat_n

  s_phi   <- df_score_phi_contrib_fast(l_cf, d_p_cf, d_np_cf, pi_p_cf, phi_hat, eta4_all)
  s_theta <- as.numeric(contrib_theta - theta_hat_n)

  score_mat <- cbind(s_phi, s_theta)

  g_phi_mean <- function(phi) colMeans(df_score_phi_contrib_fast(l_cf, d_p_cf, d_np_cf, pi_p_cf, phi, eta4_all), na.rm = TRUE)
  A11 <- df_numeric_jacobian(g_phi_mean, phi_hat)

  g_theta_mean <- function(phi) mean(efficient_theta_contrib_fast(l_cf, y_cf, d_p_cf, d_np_cf, pi_p_cf, phi, h4_all) - theta_hat_n, na.rm = TRUE)
  A21 <- df_numeric_jacobian(g_theta_mean, phi_hat)

  A_hat <- rbind(
    cbind(A11, rep(0, p_phi)),
    cbind(A21, -1)
  )

  js <- df_joint_sandwich(score_mat, A_hat)
  V  <- js$var

  phi_var <- V[1:p_phi, 1:p_phi, drop = FALSE]
  th_var_n  <- V[p_phi + 1, p_phi + 1]

  phi_se <- sqrt(pmax(diag(phi_var), 0))
  th_se_n  <- sqrt(pmax(th_var_n, 0))

  z <- stats::qnorm(0.975)
  phi_ci <- rbind(phi_hat - z * phi_se, phi_hat + z * phi_se)
  rownames(phi_ci) <- c("lower", "upper")

  theta_ci_n <- c(theta_hat_n - z * th_se_n, theta_hat_n + z * th_se_n)

  theta_var_method <- "joint_sandwich"

  # fallback (still n-based)
  if (!is.finite(th_se_n) || !is.finite(th_var_n)) {
    fb <- df_sandwich_from_contrib(contrib_theta)
    th_var_n <- fb$var
    th_se_n  <- fb$se
    theta_ci_n <- fb$ci
    theta_var_method <- "contrib_sandwich_fallback"
  }

  # Apply population scaling for theta when x_info=FALSE and N is available (or full frame was provided)
  theta_hat <- theta_hat_n * theta_scale
  th_var    <- th_var_n * (theta_scale^2)
  th_se     <- th_se_n * theta_scale
  theta_ci  <- theta_ci_n * theta_scale

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
      x_info      = isTRUE(x_info),
      aug_terms   = if (isTRUE(x_info)) "estimated" else "fixed_zero",
      convergence = res$convergence,
      theta_var_method = theta_var_method,
      sandwich_n_eff = js$n_eff,
      population_N = N_total,
      n_obs = n,
      theta_scale = theta_scale
    )
  )
}

############################################################
## 9. Efficient estimator Eff (DML1) with safe fallback
############################################################

efficient_estimator_dml1 <- function(dat,
                                     phi_start   = NULL,
                                     K           = 2,
                                     max_restart = 10,
                                     progress    = FALSE,
                                     x_info      = TRUE,
                                     N           = NULL) {

  had_frame_units <- df_dat_has_frame_units(dat)

  N_total <- NULL
  if (!isTRUE(x_info)) {
    N_total <- df_resolve_population_N(dat, N = N, require = FALSE, context = "Eff()/x_info=FALSE")
  }

  dat <- df_apply_x_info(dat, x_info)

  n <- nrow(dat)
  X <- df_get_X(dat)
  p_x <- ncol(X)
  p_phi <- 1 + p_x + 1

  theta_scale <- 1
  if (!isTRUE(x_info)) {
    if (is.null(N_total) && !isTRUE(had_frame_units)) {
      warning(
        "Eff(): x_info=FALSE but population size N was not provided and no (d_p,d_np)=(0,0) rows were detected.\n",
        "Returning the mean over the provided data (not the full-frame mean).\n",
        "To target the full-frame mean, call Eff(..., x_info=FALSE, N=<population size>).",
        call. = FALSE
      )
    }
    if (!is.null(N_total) && is.finite(N_total) && N_total > 0) {
      theta_scale <- n / N_total
    }
  }

  if (is.null(phi_start)) {
    p0 <- df_clip_prob(mean(as.numeric(dat$d_np), na.rm = TRUE))
    phi_start <- c(stats::qlogis(p0), rep(0, p_x), 0)
  }

  folds <- make_folds(n, K)
  dat_cf <- impute_pi_p_crossfit(dat, folds, sigma = NULL, progress = FALSE)

  phi_k_mat   <- matrix(NA_real_, nrow = K, ncol = p_phi)
  theta_k     <- rep(NA_real_, K)
  var_k_list  <- vector("list", K)
  w_k         <- rep(NA_real_, K)

  contrib_all <- rep(NA_real_, n)

  big_penalty <- 1e100

  if (progress) {
    cat("DML1: per-fold estimation for (phi, theta) and joint sandwich ...\n")
    flush.console()
  }

  for (k in seq_along(folds)) {
    idx_test  <- folds[[k]]
    idx_train <- setdiff(seq_len(n), idx_test)
    w_k[k]    <- length(idx_test) / n

    dat_test  <- dat_cf[idx_test,  , drop = FALSE]
    dat_train <- dat_cf[idx_train, , drop = FALSE]

    X_test <- df_get_X(dat_test)
    y_test <- as.numeric(dat_test$y)
    l_test <- df_make_l_matrix(X_test, y_test)

    d_p_test  <- as.numeric(dat_test$d_p)
    d_np_test <- as.numeric(dat_test$d_np)
    pi_p_test <- as.numeric(dat_test$pi_p)  # allow NA for (0,0); handled in df_prepare_pi_p()

    # train cache for eta4/h4 (only if x_info=TRUE)
    train_cache <- NULL
    if (isTRUE(x_info)) {
      X_all_tr <- df_get_X(dat_train)
      idx_obs  <- which(as.numeric(dat_train$d_p) == 1 | as.numeric(dat_train$d_np) == 1)

      X_obs    <- X_all_tr[idx_obs, , drop = FALSE]
      y_obs    <- as.numeric(dat_train$y[idx_obs])
      pi_p_obs <- df_clip_prob(as.numeric(dat_train$pi_p[idx_obs]))
      l_obs    <- df_make_l_matrix(X_obs, y_obs)

      prep_eta4 <- df_krr_mixed_prepare(
        X_train = X_obs,
        sigma   = NULL,
        context = "eta4*(X;phi) nuisance regression"
      )
      prep_h4 <- df_krr_mixed_prepare(
        X_train = X_obs,
        sigma   = NULL,
        context = "h4*(X;phi) nuisance regression"
      )

      train_cache <- list(
        X_obs     = X_obs,
        y_obs     = y_obs,
        pi_p_obs  = pi_p_obs,
        l_obs     = l_obs,
        prep_eta4 = prep_eta4,
        prep_h4   = prep_h4
      )
    }

    ee_fun_k <- function(phi) {
      eta4_k <- if (isTRUE(x_info)) {
        estimate_conditional_expectation_kernlab_phi_core(
          phi      = phi,
          new_X    = X_test,
          X_obs    = train_cache$X_obs,
          l_obs    = train_cache$l_obs,
          pi_p_obs = train_cache$pi_p_obs,
          prep     = train_cache$prep_eta4
        )
      } else {
        matrix(0, nrow = nrow(l_test), ncol = p_phi)
      }

      colMeans(df_score_phi_contrib_fast(l_test, d_p_test, d_np_test, pi_p_test, phi, eta4_k), na.rm = TRUE)
    }

    # --- phi estimation (optim only) ---
    # We solve for phi on this fold by minimizing || mean(score_phi) ||^2.
    # This avoids nleqslv (which can be slow/unstable in small samples).
    obj_k <- function(phi) {
      eq <- ee_fun_k(phi)
      if (any(!is.finite(eq))) return(big_penalty)
      sum(eq^2)
    }

    method_phi <- getOption("dfSEDI.phi_optim_method", "Nelder-Mead")
    method_phi <- as.character(method_phi)[1]
    if (!method_phi %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")) method_phi <- "Nelder-Mead"

    maxit_phi <- getOption("dfSEDI.phi_optim_maxit", 500L)
    maxit_phi <- as.integer(maxit_phi)[1]
    if (!is.finite(maxit_phi) || maxit_phi < 1L) maxit_phi <- 500L

    tol_obj <- getOption("dfSEDI.phi_obj_tol", 1e-4)
    tol_obj <- as.numeric(tol_obj)[1]
    if (!is.finite(tol_obj) || tol_obj <= 0) tol_obj <- 1e-4

    control_phi <- getOption("dfSEDI.phi_optim_control", list(maxit = maxit_phi))
    if (!is.list(control_phi)) control_phi <- list(maxit = maxit_phi)
    if (is.null(control_phi$maxit)) control_phi$maxit <- maxit_phi
    # --- Multi-start phi estimation (100 candidates -> try up to 10) ---
    n_cand <- 100L
    cand_mat <- matrix(rep(phi_start, each = n_cand), nrow = n_cand) +
      matrix(stats::runif(n_cand * p_phi, -0.5, 0.5), nrow = n_cand)
    obj_vals <- apply(cand_mat, 1L, obj_k)
    ord <- order(obj_vals)

    max_try <- min(10L, as.integer(max_restart))
    if (!is.finite(max_try) || max_try < 1L) max_try <- 10L
    max_try <- min(max_try, 10L)

    best <- NULL
    success <- FALSE
    for (a in seq_len(max_try)) {
      init <- cand_mat[ord[a], ]
      res_try <- try(stats::optim(par = init, fn = obj_k, method = method_phi, control = control_phi),
                     silent = TRUE)
      if (inherits(res_try, "try-error") || !is.finite(res_try$value)) {
        next
      }
      if (is.null(best) || res_try$value < best$value) best <- res_try
      if (progress) {
        cat(sprintf("  fold %d/%d: attempt %d/%d (conv=%d, obj=%.3g)\n",
                    k, K, a, max_try, res_try$convergence, res_try$value))
        flush.console()
      }
      if (isTRUE(res_try$convergence == 0) || (is.finite(res_try$value) && res_try$value <= tol_obj)) {
        success <- TRUE
        best <- res_try
        break
      }
    }

    if (!isTRUE(success) || is.null(best)) {
      if (progress) {
        cat(sprintf("  fold %d/%d: phi optimization failed after multiple starts.\n", k, K))
        flush.console()
      }
      next
    }


    if (!is.finite(best$value) || best$value >= 0.99 * big_penalty) {
      if (progress) {
        cat(sprintf("  fold %d/%d: phi failed (objective non-finite).\n", k, K))
        flush.console()
      } else {
        warning(sprintf("dfSEDI: fold %d/%d: phi failed (objective non-finite).", k, K),
                call. = FALSE)
      }
      next
    }

    phi_k <- as.numeric(best$par)

    if (progress) {
      cat(sprintf("  fold %d/%d: phi via optim (obj=%.3e).\n", k, K, best$value))
      flush.console()
    }

    phi_k_mat[k, ] <- phi_k

    if (isTRUE(x_info)) {
      eta4_k <- estimate_conditional_expectation_kernlab_phi_core(
        phi      = phi_k,
        new_X    = X_test,
        X_obs    = train_cache$X_obs,
        l_obs    = train_cache$l_obs,
        pi_p_obs = train_cache$pi_p_obs,
        prep     = train_cache$prep_eta4
      )
      h4_k <- estimate_conditional_expectation_kernlab_theta_core(
        phi      = phi_k,
        new_X    = X_test,
        X_obs    = train_cache$X_obs,
        y_obs    = train_cache$y_obs,
        l_obs    = train_cache$l_obs,
        pi_p_obs = train_cache$pi_p_obs,
        prep     = train_cache$prep_h4
      )
    } else {
      eta4_k <- matrix(0, nrow = nrow(l_test), ncol = p_phi)
      h4_k   <- rep(0, nrow(l_test))
    }

    s_phi_k   <- df_score_phi_contrib_fast(l_test, d_p_test, d_np_test, pi_p_test, phi_k, eta4_k)
    contrib_k <- efficient_theta_contrib_fast(l_test, y_test, d_p_test, d_np_test, pi_p_test, phi_k, h4_k)

    contrib_all[idx_test] <- contrib_k

    theta_k[k] <- mean(contrib_k, na.rm = TRUE)
    s_theta_k <- as.numeric(contrib_k - theta_k[k])

    score_k <- cbind(s_phi_k, s_theta_k)

    g_phi_mean_k <- function(phi) colMeans(df_score_phi_contrib_fast(l_test, d_p_test, d_np_test, pi_p_test, phi, eta4_k), na.rm = TRUE)
    A11_k <- df_numeric_jacobian(g_phi_mean_k, phi_k)

    g_theta_mean_k <- function(phi) mean(efficient_theta_contrib_fast(l_test, y_test, d_p_test, d_np_test, pi_p_test, phi, h4_k) - theta_k[k], na.rm = TRUE)
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
      info = list(type = "Eff", dml_type = "DML1", K = K, progress = progress,
                  x_info = isTRUE(x_info),
                  population_N = N_total,
                  n_obs = n,
                  theta_scale = theta_scale)
    ))
  }

  w_use <- w_k[valid]
  w_use <- w_use / sum(w_use)

  phi_hat   <- as.numeric(t(w_use) %*% phi_k_mat[valid, , drop = FALSE])

  theta_hat_n <- mean(contrib_all, na.rm = TRUE)

  valid_var <- valid[sapply(valid, function(kk) {
    Vkk <- var_k_list[[kk]]
    !is.null(Vkk) && all(is.finite(Vkk))
  })]

  V_hat <- matrix(NA_real_, nrow = p_phi + 1, ncol = p_phi + 1)
  if (length(valid_var) > 0) {
    w_var <- w_k[valid_var]
    w_var <- w_var / sum(w_var)
    V_hat <- matrix(0, nrow = p_phi + 1, ncol = p_phi + 1)
    for (j in seq_along(valid_var)) {
      kk <- valid_var[j]
      V_hat <- V_hat + (w_var[j]^2) * var_k_list[[kk]]
    }
  }

  phi_var <- V_hat[1:p_phi, 1:p_phi, drop = FALSE]
  th_var_n  <- V_hat[p_phi + 1, p_phi + 1]

  phi_se <- sqrt(pmax(diag(phi_var), 0))
  th_se_n  <- sqrt(pmax(th_var_n, 0))

  z <- stats::qnorm(0.975)
  phi_ci <- rbind(phi_hat - z * phi_se, phi_hat + z * phi_se)
  rownames(phi_ci) <- c("lower", "upper")
  theta_ci_n <- c(theta_hat_n - z * th_se_n, theta_hat_n + z * th_se_n)

  theta_var_method <- "fold_joint_sandwich"
  if (!is.finite(th_se_n) || !is.finite(th_var_n)) {
    fb <- df_sandwich_from_contrib(contrib_all)
    th_var_n <- fb$var
    th_se_n  <- fb$se
    theta_ci_n <- fb$ci
    theta_var_method <- "contrib_sandwich_fallback"
  }

  # Apply population scaling for theta
  theta_hat <- theta_hat_n * theta_scale
  th_var    <- th_var_n * (theta_scale^2)
  th_se     <- th_se_n * theta_scale
  theta_ci  <- theta_ci_n * theta_scale

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
      progress    = progress,
      x_info      = isTRUE(x_info),
      aug_terms   = if (isTRUE(x_info)) "estimated" else "fixed_zero",
      theta_var_method = theta_var_method,
      n_contrib_eff = sum(is.finite(contrib_all)),
      population_N = N_total,
      n_obs = n,
      theta_scale = theta_scale
    )
  )
}

############################################################
## 10. Sub-efficient estimator Eff_S (now supports binary via kernel logistic)
############################################################

subefficient_contrib <- function(dat, mu_hat) {
  d_np <- as.numeric(dat$d_np)
  d_p  <- as.numeric(dat$d_p)
  y    <- as.numeric(dat$y)
  mu   <- as.numeric(mu_hat)

  pi_p_raw <- as.numeric(dat$pi_p)

  w_p <- rep(0, length(y))  # this equals d_p / pi_p
  idx_p <- which(d_p == 1)
  if (length(idx_p) > 0L) {
    pi_use <- df_clip_prob(pi_p_raw[idx_p])
    if (any(!is.finite(pi_use))) {
      stop("Eff_S: `pi_p` is missing (NA/Inf) for some d_p==1 units.", call. = FALSE)
    }
    w_p[idx_p] <- 1 / pi_use
  }

  d_np * y +
    (1 - d_np) * (w_p * y + (1 - w_p) * mu)
}

subefficient_estimator_dml2 <- function(dat, K = 2, logit = NULL, progress = FALSE, x_info = TRUE) {
  dat <- df_apply_x_info(dat, x_info)

  n <- nrow(dat)
  folds <- make_folds(n, K)

  use_logit <- df_infer_logit_flag(logit, dat$y)
  if (isTRUE(use_logit) && !df_is_binary_y(dat$y)) {
    stop("Eff_S: logit=TRUE requires binary y (0/1 numeric, logical, or 2-level factor).", call. = FALSE)
  }

  X_all <- df_get_X(dat)
  idx_mu <- which(as.numeric(dat$d_p) == 1 & as.numeric(dat$d_np) == 0)

  sigma_mu <- NULL
  if (length(idx_mu) >= 2L) {
    cat_cols <- df_infer_categorical_cols(X_all)
    cont_cols <- setdiff(seq_len(ncol(X_all)), cat_cols)
    if (length(cont_cols) > 0L) {
      sigma_mu <- df_estimate_rbf_sigma(X_all[idx_mu, cont_cols, drop = FALSE], max_n = 2000L)
    }
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

    mu_k <- if (isTRUE(use_logit)) {
      regression_expectation_klogit(dat_train, new_X = X_test, sigma = sigma_mu)
    } else {
      regression_expectation_kernlab(dat_train, new_X = X_test, sigma = sigma_mu)
    }

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
       info  = list(type = "Eff_S", K = K, progress = progress,
                    x_info = isTRUE(x_info),
                    logit = use_logit,
                    mu_model = if (isTRUE(use_logit)) "kernel_logistic" else "kernel_ridge"))
}

############################################################
## 11. Parametric efficient estimator Eff_P (now supports binary via logistic glm)
############################################################

efficient_parametric_estimator <- function(dat,
                                           phi_start = NULL,
                                           eta4_star = 0,
                                           max_iter  = 20,
                                           logit     = NULL,
                                           progress  = FALSE,
                                           x_info    = TRUE) {
  dat <- df_apply_x_info(dat, x_info)

  # If pi_p is missing for np-only union units, impute it (single fit, no cross-fitting).
  if (any(as.numeric(dat$d_np) == 1 & as.numeric(dat$d_p) == 0 & !is.finite(as.numeric(dat$pi_p)))) {
    dat <- impute_pi_p_simple(dat, sigma = NULL)
  }

  X_all <- df_get_X(dat)
  p_x <- ncol(X_all)
  p_phi <- 1 + p_x + 1
  n <- nrow(dat)

  use_logit <- df_infer_logit_flag(logit, dat$y)
  if (isTRUE(use_logit) && !df_is_binary_y(dat$y)) {
    stop("Eff_P: logit=TRUE requires binary y (0/1 numeric, logical, or 2-level factor).", call. = FALSE)
  }

  if (is.null(phi_start)) {
    p0 <- df_clip_prob(mean(as.numeric(dat$d_np), na.rm = TRUE))
    phi_start <- c(stats::qlogis(p0), rep(0, p_x), 0)
  }

  # Prepare pi_p once (allows NA for (0,0) rows; union rows must be finite).
  d_p1  <- as.numeric(dat$d_p)
  d_np1 <- as.numeric(dat$d_np)
  pi_p1 <- df_prepare_pi_p(dat$pi_p, d_p = d_p1, d_np = d_np1, context = "Eff_P")

  ee_para <- function(phi) {
    X1    <- df_get_X(dat)
    y1    <- as.numeric(dat$y)
    l1    <- as.matrix(cbind(1, X1, y1))

    phi <- as.numeric(phi)
    eta      <- as.numeric(l1 %*% phi)
    pi_np1   <- df_clip_prob(1 / (1 + exp(-eta)))
    pi_np_p1 <- df_clip_prob(pi_np1 + pi_p1 - pi_np1 * pi_p1)

    term1 <- l1 * pi_np1 * (1 - pi_np1) / pi_np_p1 *
      (d_np1 / pi_np1 * (pi_p1 - d_p1) -
         d_p1 / (1 - pi_np1) * (1 - d_np1 / pi_np1))

    term2_coef <- 1 -
      (d_p1 * d_np1) / (pi_p1 * pi_np1) -
      1 / pi_np_p1 *
      (d_np1 * (1 - d_p1 / pi_p1) +
         d_p1 * (1 - d_np1 / pi_np1))

    est_eq <- term1 + term2_coef * 0
    colMeans(est_eq, na.rm = TRUE)
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
    mu_model <- NA_character_
  } else {
    phi_hat <- as.numeric(sol$x)

    if (!isTRUE(x_info)) {
      # x_info=FALSE: no need to fit working h4/mu model; set to 0 (fast path)
      mu_hat <- rep(0, n)
      mu_model <- "fixed_zero(x_info=FALSE)"
    } else {
      subset_idx <- which(as.numeric(dat$d_np) == 1 | as.numeric(dat$d_p) == 1)
      dat_sub <- dat[subset_idx, ]
      X_sub <- df_get_X(dat_sub)

      X_full <- df_get_X(dat)

      if (isTRUE(use_logit)) {
        y01_sub <- df_as_binary01(dat_sub$y)$y01
        df_sub <- data.frame(y = y01_sub, X_sub)
        colnames(df_sub)[-1] <- paste0("x", seq_len(ncol(X_sub)))

        glm_fit <- stats::glm(y ~ ., data = df_sub, family = stats::binomial())

        df_full <- data.frame(X_full)
        colnames(df_full) <- paste0("x", seq_len(ncol(X_full)))
        mu_hat <- as.numeric(stats::predict(glm_fit, newdata = df_full, type = "response"))
        mu_hat <- df_clip_prob(mu_hat)

        mu_model <- "parametric_logistic"
      } else {
        df_sub <- data.frame(y = as.numeric(dat_sub$y), X_sub)
        colnames(df_sub)[-1] <- paste0("x", seq_len(ncol(X_sub)))

        lm_fit <- stats::lm(y ~ ., data = df_sub)

        df_full <- data.frame(X_full)
        colnames(df_full) <- paste0("x", seq_len(ncol(X_full)))
        mu_hat <- as.numeric(stats::predict(lm_fit, newdata = df_full))

        mu_model <- "parametric_linear"
      }
    }

    contrib <- efficient_theta_contrib(dat, phi_hat, h4_star_local = mu_hat)
    theta_res <- df_sandwich_from_contrib(contrib)
  }

  list(
    phi    = phi_hat,
    phi_var= NULL,
    phi_se = NULL,
    phi_ci = NULL,
    theta  = theta_res$theta,
    var    = theta_res$var,
    se     = theta_res$se,
    ci     = theta_res$ci,
    info   = list(type = "Eff_P",
                  phi_start = phi_start,
                  eta4_star = eta4_star,
                  max_iter  = max_iter,
                  progress  = progress,
                  x_info    = isTRUE(x_info),
                  logit     = use_logit,
                  mu_model  = mu_model)
  )
}

############################################################
## 12. Public wrappers
############################################################

Eff <- function(dat,
                K           = 2,
                phi_start   = NULL,
                max_restart = 10,
                type        = 2, #default
                dml_type    = 1,
                progress    = interactive(),
                x_info      = TRUE,
                N           = NULL) {

  if (!is.null(type)) dml_type <- type

  if (is.numeric(dml_type)) dml_type <- if (dml_type == 1) "DML1" else "DML2"
  dml_type <- toupper(as.character(dml_type))
  if (!dml_type %in% c("DML1", "DML2")) stop("Eff(): dml_type must be 1, 2, 'DML1', or 'DML2'.")

  if (dml_type == "DML2") {
    efficient_estimator_dml2(dat, phi_start = phi_start, K = K, max_restart = max_restart,
                             progress = progress, x_info = x_info, N = N)
  } else {
    efficient_estimator_dml1(dat, phi_start = phi_start, K = K, max_restart = max_restart,
                             progress = progress, x_info = x_info, N = N)
  }
}

Eff_S <- function(dat, K = 2, logit = NULL, progress = interactive(), x_info = TRUE) {
  subefficient_estimator_dml2(dat, K = K, logit = logit, progress = progress, x_info = x_info)
}

Eff_P <- function(dat,
                  phi_start = NULL,
                  eta4_star = 0,
                  max_iter  = 20,
                  logit     = NULL,
                  progress  = interactive(),
                  x_info    = TRUE) {
  efficient_parametric_estimator(dat, phi_start = phi_start, eta4_star = eta4_star,
                                 max_iter = max_iter, logit = logit, progress = progress,
                                 x_info = x_info)
}
