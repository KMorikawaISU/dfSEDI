############################################################
## dualframe_simulation.R
##
## Simple data generator for examples:
##   generate_dualframe_population()
############################################################

#' Generate a dual-frame population (simulation example)
#'
#' This function reproduces the scenario [O2] × [NP1] × P:
#'   - X ~ N(0, 1)
#'   - Y | X ~ N(X/2, 1/2^2)
#'   - Non-probability inclusion: logit pi_NP = -2.15 - 0.5 X - 0.75 Y
#'   - Probability inclusion:
#'       pi_P = 0.005 + 1 / (1 + exp(3 + (X - 2)^2 / 4))
#'
#' It is intended only for simulation / documentation examples.
#'
#' @param N Population size.
#' @return A data.frame with columns x, y, d_np, d_p, pi_p, pi_np.
#' @export
generate_dualframe_population <- function(N) {
  x <- stats::rnorm(N, 0, 1)
  y <- stats::rnorm(N, x / 2, 1 / 2)

  pi_np <- 1 / (1 + exp(2.15 + 0.5 * x + 0.75 * y))
  pi_p  <- 0.005 + 1 / (1 + exp(3 + (x - 2)^2 / 4))

  d_np <- stats::rbinom(N, 1, pi_np)
  d_p  <- stats::rbinom(N, 1, pi_p)

  data.frame(
    x     = x,
    y     = y,
    d_np  = d_np,
    d_p   = d_p,
    pi_p  = pi_p,
    pi_np = pi_np
  )
}
