############################################################
## inst/examples/dualframe_simulation.R
############################################################

library(dfSEDI)

generate_dualframe_population <- function(N) {
  # Continuous and discrete covariates
  x1 <- rnorm(N, 0, 1)          # continuous
  x2 <- rbinom(N, 1, 1 / 2)     # binary (discrete)

  # Full X for parametric components (logistic model, etc.)
  X  <- cbind(x1, x2)

  mu_y <- 0.5 * x1 - 0.4 * x2 + 0.2
  y    <- rnorm(N, mean = mu_y, sd = 0.5)

  lin_np <- -1.9 - 0.5 * x1 - 0.25 * x2 - 0.75 * y
  pi_np  <- 1 / (1 + exp(-lin_np))

  lin_p <- -3.0 - (x1 - 2)^2 / 4
  pi_p  <- 0.005 + 1 / (1 + exp(-lin_p))

  d_np <- rbinom(N, 1, pi_np)
  d_p  <- rbinom(N, 1, pi_p)

  data.frame(
    X      = I(X),   # X is a matrix-column (x1, x2) for parametric part
    y      = y,
    d_np   = d_np,
    d_p    = d_p,
    pi_p   = pi_p,
    pi_np  = pi_np
  )
}

############################################################
## One Monte Carlo replication
############################################################

simulate_one_replication <- function(seed,
                                     N        = 10000,
                                     K        = 2,
                                     progress = FALSE,
                                     x_info   = TRUE) {

  set.seed(seed)
  dat <- generate_dualframe_population(N)

  # If you want to emulate "no (0,0) units observed", you can set x_info=FALSE:
  # dat <- dat[dat$d_np == 1 | dat$d_p == 1, , drop = FALSE]
  # The estimators below also accept x_info and will drop (0,0) internally.
  # (But note: SM is then a sample mean, not a population mean.)

  # Simple mean (SM)
  res_SM <- df_sandwich_from_contrib(dat$y)

  # Basic estimators: P, NP, NP+P
  fit_P   <- df_estimate_P(dat)
  fit_NP  <- df_estimate_NP(dat)
  fit_NPP <- df_estimate_NP_P(dat)

  # Main estimators: Eff, Eff_S, Eff_P
  fit_Eff <- Eff(
    dat         = dat,
    K           = K,
    phi_start   = NULL,
    max_restart = 10,
    progress    = progress,
    x_info      = x_info
  )

  fit_EffS <- Eff_S(
    dat      = dat,
    K        = K,
    progress = progress,
    x_info   = x_info
  )

  fit_EffP <- Eff_P(
    dat       = dat,
    phi_start = NULL,
    eta4_star = 0,
    max_iter  = 20,
    progress  = progress,
    x_info    = x_info
  )

  # Collect theta (point estimates)
  theta_vec <- c(
    SM    = res_SM$theta,
    P     = fit_P$theta,
    NP    = fit_NP$theta,
    NP_P  = fit_NPP$theta,
    Eff_S = fit_EffS$theta,
    Eff_P = fit_EffP$theta,
    Eff   = fit_Eff$theta
  )

  list(theta = theta_vec,
       Eff   = fit_Eff,
       Eff_S = fit_EffS,
       Eff_P = fit_EffP)
}

############################################################
## main(): single dataset test, no Monte Carlo
############################################################

main <- function(N = 10000, K = 2, x_info = TRUE) {

  dat <- generate_dualframe_population(N)

  # For debugging: true phi for this DGP (for example)
  phi_true <- c(-2.15, -0.5, -0.75)

  # 1) Efficient estimator Eff
  fit_eff <- Eff(
    dat         = dat,
    K           = K,
    phi_start   = phi_true,
    max_restart = 10,
    progress    = TRUE,
    x_info      = x_info
  )

  cat("\n=== Eff ===\n")
  print(fit_eff$phi)
  print(fit_eff$theta)
  print(fit_eff$se)
  print(fit_eff$ci)

  # 2) Sub-efficient estimator Eff_S
  fit_effS <- Eff_S(
    dat      = dat,
    K        = K,
    progress = TRUE,
    x_info   = x_info
  )

  cat("\n=== Eff_S ===\n")
  print(fit_effS$theta)
  print(fit_effS$se)
  print(fit_effS$ci)

  # 3) Parametric efficient estimator Eff_P
  fit_effP <- Eff_P(
    dat       = dat,
    phi_start = phi_true,
    eta4_star = 0,
    max_iter  = 20,
    progress  = TRUE,
    x_info    = x_info
  )

  cat("\n=== Eff_P ===\n")
  print(fit_effP$phi)
  print(fit_effP$theta)
  print(fit_effP$se)
  print(fit_effP$ci)

  invisible(list(
    dat     = dat,
    Eff     = fit_eff,
    Eff_S   = fit_effS,
    Eff_P   = fit_effP
  ))
}
