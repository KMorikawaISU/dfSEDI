
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

# ------------------------
# Settings
# ------------------------
B <- 500
N <- 10000
Scenario <- 1
K <- 2

seeds <- seed_start <- 1 + seq_len(B) - 1

# Output directory
sim_dir <- "mc_out"
dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

# Create a tag for filenames
tag <- sprintf(
  "S%s_N%s_B%s_K%s_%s",
  Scenario, N, B, K,
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

rds_file <- file.path(sim_dir, paste0("res_par_", tag, ".rds"))
csv_file <- file.path(sim_dir, paste0("res_par_", tag, ".csv"))
txt_file <- file.path(sim_dir, paste0("summary_", tag, ".txt"))

# ------------------------
# Load bundled example functions
# ------------------------
example_file <- system.file("examples", "dualframe_simulation.R", package = "dfSEDI")
source(example_file)

# ------------------------
# Parallel cluster
# ------------------------

#Note: When I set the number of workers to (#cores - 1L), my PC crashed.
n_workers <- max(1L, parallel::detectCores() - 2L)
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

# ------------------------
# One replicate
# ------------------------
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

# ------------------------
# Run MC
# ------------------------

library(dplyr)
library(purrr)

B <- 500  # 例：MCの回数

results <- purrr::map_dfr(seq_len(B), function(b){

  # 再現性が欲しければ（任意）
  # set.seed(123 + b)

  run_one_mc(
    mc_id        = b,
    frame        = frame,
    prob         = prob,
    nonprob      = nonprob,
    y_var        = y_var,
    x_cont_vars  = x_cont_vars,
    x_cat_vars   = x_cat_vars,
    np_groups    = np_groups,
    K_eff        = K_eff,
    Eff_dml_type = Eff_dml_type,
    K_effS       = K_effS,
    progress     = progress_each_fit,
    run_full_dat = TRUE,
    np_max_iter  = np_max_iter,
    use_hajek    = use_hajek
  )
})

# 保存
dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(sim_dir, "dfSEDI_results_withEff_mc.csv")
write.csv(results, out_csv, row.names = FALSE)


out_list <- pbapply::pblapply(seeds, one_rep, cl = cl)
res_par <- do.call(rbind, out_list)

# ------------------------
# Save results
# ------------------------
saveRDS(res_par, file = rds_file)
utils::write.csv(res_par, file = csv_file, row.names = FALSE)

# Also save text summary output (optional but handy)
zz <- file(txt_file, open = "wt")
sink(zz)
cat("MC settings\n")
cat(sprintf("B=%s, N=%s, Scenario=%s, K=%s\n\n", B, N, Scenario, K))
cat("seed_start:\n")
print(seed_start)
cat("\nseeds:\n")
print(seeds)
cat("\nsummary_mc:\n")
print(summarize_mc(res_par, theta_true = 0))
sink()
close(zz)

message("Saved RDS: ", rds_file)
message("Saved CSV: ", csv_file)
message("Saved TXT: ", txt_file)

