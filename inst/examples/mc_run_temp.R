devtools::install_github("KMorikawaISU/dfSEDI",force=TRUE)
# install.packages("devtools")
library(dplyr)
library(ggplot2)

#install.packages(c("ggplot2", "dplyr"))
#install.packages(c("ggthemes", "viridis"))

base_fun <- function(X) {
  cbind(1, X, X[, 1]^2)
}


# install.packages("pbapply")  # if needed

library(dfSEDI)
library(parallel)
library(pbapply)

example_file <- system.file("examples", "dualframe_simulation.R", package = "dfSEDI")
source(example_file)

B <- 30
N <- 10000
Scenario <- 1
K <- 2

seeds <- 1 + seq_len(B) - 1

n_workers <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_workers)
on.exit(parallel::stopCluster(cl), add = TRUE)

pbapply::pboptions(type = "txt")

parallel::clusterEvalQ(cl, library(dfSEDI))

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
    extract_row(fits$P,         "P",             seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP,        "NP",            seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$NP_P,      "NP_P",          seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff,       "Eff",           seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_union, "Eff_union_dat", seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_S,     "Eff_S",         seed, Scenario, fits$n_np, fits$n_p, fits$n_union),
    extract_row(fits$Eff_P,     "Eff_P",         seed, Scenario, fits$n_np, fits$n_p, fits$n_union)
  )
}

systime <- system.time(
  out_list <- pbapply::pblapply(seeds, one_rep, cl = cl)
)
res_par <- do.call(rbind, out_list)

summarize_mc(res_par, theta_true = 0)







#====theta plot====#


have_ggthemes <- requireNamespace("ggthemes", quietly = TRUE)
have_viridis  <- requireNamespace("viridis", quietly = TRUE)

df_theta <- res_par %>%
  filter(is.finite(theta)) %>%
  mutate(estimator = factor(estimator))

p_theta_jitter <- ggplot(df_theta, aes(x = estimator, y = theta, fill = estimator)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.85) +
  geom_jitter(aes(color = estimator), width = 0.15, alpha = 0.25, size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "theta") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

if (have_ggthemes) p_theta_jitter <- p_theta_jitter + ggthemes::theme_few()
if (have_viridis) {
  p_theta_jitter <- p_theta_jitter +
    scale_fill_viridis_d(option = "C") +
    scale_color_viridis_d(option = "C")
}

p_theta_jitter






