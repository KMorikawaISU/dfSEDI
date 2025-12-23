
# mc_run.R
#
# Monte Carlo run (PSOCK + pbapply) and boxplots for theta and phi.
# This script assumes dfSEDI is installed and that the bundled example script
# 'inst/examples/dualframe_simulation.R' is available in the package.

# install.packages("devtools")
# devtools::install_github("KMorikawaISU/dfSEDI")

# install.packages("pbapply")  # if needed

# mc_run.R
#
# Monte Carlo run (PSOCK + pbapply) and boxplots for theta and phi.
# This script assumes dfSEDI is installed and that the bundled example script
# 'inst/examples/dualframe_simulation.R' is available in the package.

# install.packages("devtools")
# devtools::install_github("KMorikawaISU/dfSEDI")

# install.packages("pbapply")  # if needed

# mc_run.R
#
# Monte Carlo run (PSOCK + pbapply), save results, then load and plot.
# This script assumes dfSEDI is installed and that the bundled example script
# 'inst/examples/dualframe_simulation.R' is available in the package.
devtools::install_github("KMorikawaISU/dfSEDI")



library(dfSEDI)
library(parallel)
library(pbapply)

example_file <- system.file("examples", "dualframe_simulation.R", package = "dfSEDI")
source(example_file)

B <- 3
N <- 10000
Scenario <- 1
K <- 2

set.seed(18)
N <- 10000
dat <- generate_dualframe_population(N = N)

base_fun <- function(X) {
  cbind(1, X, X[, 1]^2)
}

seeds <- seed_start <- 1 + seq_len(B) - 1

n_workers <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_workers)
on.exit(parallel::stopCluster(cl), add = TRUE)

pbapply::pboptions(type = "txt")

parallel::clusterEvalQ(cl, library(dfSEDI))

# Export the functions defined by the example script
parallel::clusterExport(
  cl,
  varlist = c(
    "normalize_scenario",
    "generate_dualframe_population",
    "make_base_fun",
    "safe_fit",
    "extract_row",
    "fit_all_estimators_once",
    "N",
    "Scenario",
    "K"
  ),
  envir = environment()
)

one_rep <- function(seed) {
  set.seed(seed)
  dat <- generate_dualframe_population(N = N, Scenario = Scenario)
  fits <- fit_all_estimators_once(dat = dat, Scenario = Scenario, K = K, progress_each = FALSE)

  rbind(
    extract_row(fits$P,          "P",          seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP,         "NP",         seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP_P,       "NP_P",       seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_S,      "Eff_S",      seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_P,      "Eff_P",      seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff1_union, "Eff1_union", seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff2_union, "Eff2_union", seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff1,       "Eff1",       seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff2,       "Eff2",       seed, Scenario, fits$n_np, fits$n_p, fits$n_union)
  )
}

out_list <- pbapply::pblapply(seeds, one_rep, cl = cl)
res_par <- do.call(rbind, out_list)

summarize_mc(res_par, theta_true = 0)


# ------------------------
# Boxplots (theta, phi)
# ------------------------

library(dplyr)
library(ggplot2)

est_order <- c(
  "P", "NP", "NP_P",
  "Eff_S", "Eff_P",
  "Eff1_union", "Eff2_union",
  "Eff1", "Eff2"
)

res_plot <- res_par %>%
  mutate(
    Scenario = toupper(as.character(Scenario)),
    Scenario = ifelse(grepl("^S", Scenario), Scenario, paste0("S", Scenario)),
    estimator = factor(as.character(estimator), levels = est_order),
    Scenario  = factor(Scenario)
  )

# Theta boxplot
ggplot(res_plot, aes(x = estimator, y = theta)) +
  geom_boxplot(na.rm = TRUE, outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.25, na.rm = TRUE) +
  facet_wrap(~ Scenario, ncol = 1) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = expression(hat(theta)))


# Phi (y-coefficient): S1/S2 -> phi_3, S3 -> phi_4
res_phi <- res_plot %>%
  mutate(
    phi_target = case_when(
      Scenario %in% c("S1", "S2") ~ phi_3,
      Scenario %in% c("S3")       ~ phi_4,
      TRUE                        ~ NA_real_
    ),
    phi_target_name = case_when(
      Scenario %in% c("S1", "S2") ~ "phi_3",
      Scenario %in% c("S3")       ~ "phi_4",
      TRUE                        ~ NA_character_
    )
  ) %>%
  filter(estimator %in% c("NP", "NP_P", "Eff1_union", "Eff2_union", "Eff1", "Eff2"))

ggplot(res_phi, aes(x = estimator, y = phi_target)) +
  geom_boxplot(na.rm = TRUE, outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.25, na.rm = TRUE) +
  facet_wrap(~ Scenario, ncol = 1) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "phi (S1–S2: phi_3, S3: phi_4)")
