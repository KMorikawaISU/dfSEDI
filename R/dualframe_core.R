# -----------------------------------------------------------
# Median heuristic for RBF kernel parameter in kernlab::rbfdot
#   k(x,x') = exp( -sigma * ||x-x'||^2 )
# Use: sigma = 1 / (2 * med^2), where med = median pairwise distance
# To avoid O(n^2) cost, we subsample up to max_n points.
# -----------------------------------------------------------
median_heuristic_rbf_sigma <- function(X, max_n = 1000L, seed = 1L) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n <= 1L) return(1)

  m <- min(n, as.integer(max_n))
  if (m < n) {
    set.seed(as.integer(seed))
    idx <- sample.int(n, size = m, replace = FALSE)
    Xs <- X[idx, , drop = FALSE]
  } else {
    Xs <- X
  }

  d <- as.numeric(stats::dist(Xs))
  d <- d[is.finite(d) & d > 0]
  if (length(d) == 0L) return(1)

  med <- stats::median(d)
  if (!is.finite(med) || med <= 0) return(1)

  1 / (2 * med^2)
}
