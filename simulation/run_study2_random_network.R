# =====================================================================
# Main simulation script
#
# Two settings:
#   Setting 1 — correctly specified: truth network = working network
#     CGa/a : W_int = raw G_A, W_dep = raw G_A  (broad interference model)
#     CGb/b : W_int = raw G_M, W_dep = raw G_M  (narrow dependence model)
#     CGc/c : W_int = raw G_M, W_dep = raw G_A  (separated model — reference)
#   Setting 2 — misspecified: truth always CGc, working model varies
#     CGc/a : estimate under G_A  (over-interference)
#     CGc/b : estimate under G_M  (under-dependence)
#
# Pipeline:
#   §6  — PL estimation, all scenarios (sequential) → pl_estimates.csv
#   §7  — Causal estimation with estimated params   (parallel)
#   §8  — Causal estimation with true params        (parallel)
#   §9  — Results table
#
# Saved outputs (results/):
#   networks.rds        — base_net + net_CGa/b/c
#   reps_CGa/b/c.rds    — generated replications per truth network
#   pl_estimates.csv    — PL estimates, all scenarios combined
#
# Usage: source("simulation/run_study2_random_network.R")   or step through section by section
# =====================================================================

source("R/utils.R");      source("simulation/sim_utils.R")
source("R/estimation.R"); source("R/causal.R")

# ── 1. Configuration ──────────────────────────────────────────────────
N       <- 200
S       <- 500
dep_deg <- 4
int_deg <- 1

dgp_burn <- 1500   # independent-chain burn-in; save sweep 1501
dgp_parallel <- TRUE
dgp_n_cores <- max(1L, parallel::detectCores() - 1L)

alpha  <- 0.5     # Bernoulli treatment allocation probability
K      <- 50      # inner-loop Gibbs iterations (K in manuscript)
m_star <- 20     # burn-in for causal Gibbs sampler (m* in Auto-G paper)

params  <- make_parameters()
out_dir <- file.path("results", sprintf("new_new_N%d_dep%g_int%g_beta3_%s_theta_%s_S%d",
                     N, dep_deg, int_deg,
                     gsub("\\.", "p", as.character(params$beta3)),
                     gsub("\\.", "p", as.character(params$theta)),
                     S))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_run_settings <- function(out_dir, params, config, base_net = NULL) {
  lines <- c(
    "# MultiNetCausalRcpp run settings",
    sprintf("timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("out_dir: %s", normalizePath(out_dir, winslash = "/", mustWork = FALSE)),
    "",
    "[configuration]",
    sprintf("N: %s", config$N),
    sprintf("S: %s", config$S),
    sprintf("dep_deg_input: %s", config$dep_deg),
    sprintf("int_deg_input: %s", config$int_deg),
    sprintf("dgp_burn: %s", config$dgp_burn),
    sprintf("dgp_parallel: %s", config$dgp_parallel),
    sprintf("dgp_n_cores: %s", config$dgp_n_cores),
    sprintf("alpha: %s", config$alpha),
    sprintf("K: %s", config$K),
    sprintf("m_star: %s", config$m_star),
    sprintf("base_seed_network: %s", config$base_seed_network),
    sprintf("seed_reps_CGa: %s", config$seed_reps_CGa),
    sprintf("seed_reps_CGb: %s", config$seed_reps_CGb),
    sprintf("seed_reps_CGc: %s", config$seed_reps_CGc),
    "",
    "[true_parameters]"
  )

  param_lines <- vapply(names(params), function(nm) {
    sprintf("%s: %s", nm, format(unname(params[[nm]]), scientific = FALSE, trim = TRUE))
  }, character(1L))
  lines <- c(lines, param_lines)

  if (!is.null(base_net)) {
    lines <- c(lines,
      "",
      "[realized_base_network]",
      sprintf("base_int_deg: %.6f", base_net$base_int_deg),
      sprintf("base_aux_deg: %.6f", base_net$base_aux_deg),
      sprintf("base_aux_only_deg: %.6f", base_net$base_aux_only_deg),
      sprintf("base_dep_deg: %.6f", base_net$base_dep_deg),
      sprintf("frac_isolated_int: %.6f", base_net$frac_isolated_int)
    )
  }

  writeLines(lines, file.path(out_dir, "run_settings.txt"))
}

run_config <- list(
  N = N, S = S, dep_deg = dep_deg, int_deg = int_deg,
  dgp_burn = dgp_burn, dgp_parallel = dgp_parallel, dgp_n_cores = dgp_n_cores,
  alpha = alpha, K = K, m_star = m_star,
  base_seed_network = 123,
  seed_reps_CGa = 1,
  seed_reps_CGb = 2,
  seed_reps_CGc = 3
)

# ── 2. Base network ───────────────────────────────────────────────────
base_net <- make_random_network(N, dep_deg, int_deg, seed = run_config$base_seed_network)
net_CGa  <- make_network_spec(base_net, "CGa")
net_CGb  <- make_network_spec(base_net, "CGb")
net_CGc  <- make_network_spec(base_net, "CGc")
save_run_settings(out_dir, params, run_config, base_net)
cat("Run settings saved to", file.path(out_dir, "run_settings.txt"), "\n")

# ── 3. Quick diagnostic plot ──────────────────────────────────────────
reps_check <- generate_replications(net_CGc, params, S = 1, seed = 123,
                                     burn = 50)
plot_networks(base_net, reps_check[[1]]$A, reps_check[[1]]$Y)

# ── 3b. Quick causal diagnostic ───────────────────────────────────────
# One rep, small K: verify pipeline before the full S=500 run.
cat("\nQuick causal diagnostic (CGc, 1 rep)...\n")
quick_causal_check(net_CGc, params)

# ── 4. Quick DGP timing + generate replications → save RDS ───────────
cat("\nCalibrating DGP generation time...\n")
quick_dgp_timing(net_CGc, params, S_total = S, burn = dgp_burn,
                 n_cores = dgp_n_cores, n_groups = 3)

cat("\nGenerating replications...\n")
reps_CGa <- generate_replications(net_CGa, params, S, seed = 1,
                                   burn = dgp_burn, parallel = dgp_parallel,
                                   n_cores = dgp_n_cores, label = "DGP CGa")
reps_CGb <- generate_replications(net_CGb, params, S, seed = 2,
                                   burn = dgp_burn, parallel = dgp_parallel,
                                   n_cores = dgp_n_cores, label = "DGP CGb")
reps_CGc <- generate_replications(net_CGc, params, S, seed = 3,
                                   burn = dgp_burn, parallel = dgp_parallel,
                                   n_cores = dgp_n_cores, label = "DGP CGc")

saveRDS(list(base_net = base_net,
             net_CGa  = net_CGa,
             net_CGb  = net_CGb,
             net_CGc  = net_CGc),
        file.path(out_dir, "networks.rds"))
saveRDS(reps_CGa, file.path(out_dir, "reps_CGa.rds"))
saveRDS(reps_CGb, file.path(out_dir, "reps_CGb.rds"))
saveRDS(reps_CGc, file.path(out_dir, "reps_CGc.rds"))
cat("Networks and replications saved to", out_dir, "\n")

# ── 5. PL estimation function ─────────────────────────────────────────
run_pl_estimation <- function(reps, work_net) {
  network_pl <- as_network(work_net)
  lapply(seq_along(reps), function(s) {
    rep <- reps[[s]]
    fit_pl(rep$Y, rep$A, rep$L, network_pl)
  })
}

# ── 6. PL estimation — all scenarios (sequential) ─────────────────────
cat("\n=== PL estimation ===\n")
pl_CGaa <- run_pl_estimation(reps_CGa, net_CGa)
pl_CGbb <- run_pl_estimation(reps_CGb, net_CGb)
pl_CGcc <- run_pl_estimation(reps_CGc, net_CGc)
pl_CGca <- run_pl_estimation(reps_CGc, net_CGa)
pl_CGcb <- run_pl_estimation(reps_CGc, net_CGb)

# Save: one CSV, rows = instances, cols = parameters
pl_to_df <- function(pl_list, scenario) {
  mat <- do.call(rbind, lapply(pl_list, params_flat_vector))
  data.frame(scenario = scenario, rep = seq_len(nrow(mat)), mat,
             check.names = FALSE)
}
pl_all <- rbind(
  pl_to_df(pl_CGaa, "CGaa"),
  pl_to_df(pl_CGbb, "CGbb"),
  pl_to_df(pl_CGcc, "CGcc"),
  pl_to_df(pl_CGca, "CGca"),
  pl_to_df(pl_CGcb, "CGcb")
)
write.csv(pl_all, file.path(out_dir, "pl_estimates.csv"), row.names = FALSE)
cat("PL estimates saved to", file.path(out_dir, "pl_estimates.csv"), "\n")

cat("\n========== PL ESTIMATION SUMMARY ==========\n")
summarize_pl(pl_CGaa, params, "CGa/a — correct spec")
summarize_pl(pl_CGbb, params, "CGb/b — correct spec")
summarize_pl(pl_CGcc, params, "CGc/c — correct spec  [reference]")
summarize_pl(pl_CGca, params, "CGc/a — misspec: over-interference")
summarize_pl(pl_CGcb, params, "CGc/b — misspec: under-dependence")

# ── 7. Parallel cluster setup (causal estimation only) ────────────────
n_cores <- max(1L, parallel::detectCores() - 1L)
cl      <- parallel::makeCluster(n_cores)
wd      <- getwd()
parallel::clusterExport(cl, c("wd", "m_star"), envir = environment())
invisible(parallel::clusterCall(cl, function() {
  setwd(wd)
  source("R/utils.R");   source("simulation/sim_utils.R")
  source("R/estimation.R"); source("R/causal.R")
}))
cat(sprintf("Parallel cluster ready: %d workers\n", n_cores))

cat("\nCalibrating causal estimation time...\n")
quick_causal_timing(reps_CGc, net_CGc, pl_CGcc, n_cores = n_cores,
                    alpha = alpha, K = K, burn = m_star)

# ── 8. Causal estimation — estimated params ────────────────────────────
cat("\n=== Causal estimation (estimated params) ===\n")
eff_CGaa <- run_causal_estimation(
  reps_CGa, net_CGa, pl_CGaa, alpha = alpha, K = K, burn = m_star,
  cl = cl, label = "causal est CGa/a"
)
eff_CGbb <- run_causal_estimation(
  reps_CGb, net_CGb, pl_CGbb, alpha = alpha, K = K, burn = m_star,
  cl = cl, label = "causal est CGb/b"
)
eff_CGcc <- run_causal_estimation(
  reps_CGc, net_CGc, pl_CGcc, alpha = alpha, K = K, burn = m_star,
  cl = cl, label = "causal est CGc/c"
)
eff_CGca <- run_causal_estimation(
  reps_CGc, net_CGa, pl_CGca, alpha = alpha, K = K, burn = m_star,
  cl = cl, label = "causal est CGc/a"
)
eff_CGcb <- run_causal_estimation(
  reps_CGc, net_CGb, pl_CGcb, alpha = alpha, K = K, burn = m_star,
  cl = cl, label = "causal est CGc/b"
)

# ── 9. Causal estimation — true params ───────────────────────────────
# Setting 1: each spec's own network + true params (model-specific benchmark).
# Setting 2: misspecified network + true params (oracle under wrong model).
#   Bias = PL estimation error alone; network misspecification is visible
#   as the difference between Setting 2 "true params" and Setting 1 truth.
cat("\n=== Causal estimation (true params) ===\n")
eff_CGaa_true <- run_causal_estimation(
  reps_CGa, net_CGa, params, alpha = alpha, K = K, burn = m_star,
  cl = cl, label = "causal true CGa/a"
)
eff_CGbb_true <- run_causal_estimation(
  reps_CGb, net_CGb, params, alpha = alpha, K = K, burn = m_star,
  cl = cl, label = "causal true CGb/b"
)
eff_CGcc_true <- run_causal_estimation(
  reps_CGc, net_CGc, params, alpha = alpha, K = K, burn = m_star,
  cl = cl, label = "causal true CGc/c"
)
# Misspec oracle: CGc data, wrong network, true params
eff_CGca_true <- run_causal_estimation(
  reps_CGc, net_CGa, params, alpha = alpha, K = K, burn = m_star,
  cl = cl, label = "causal true CGc/a"
)
eff_CGcb_true <- run_causal_estimation(
  reps_CGc, net_CGb, params, alpha = alpha, K = K, burn = m_star,
  cl = cl, label = "causal true CGc/b"
)

parallel::stopCluster(cl)

# ── 10. Save results to CSV ───────────────────────────────────────────
s1 <- list(
  CGa = list(est = eff_CGaa, true = eff_CGaa_true),
  CGb = list(est = eff_CGbb, true = eff_CGbb_true),
  CGc = list(est = eff_CGcc, true = eff_CGcc_true)
)
s2 <- list(
  CGa = list(est = eff_CGca, true = eff_CGca_true),
  CGb = list(est = eff_CGcb, true = eff_CGcb_true),
  CGc = list(est = eff_CGcc, true = eff_CGcc_true)
)
save_results(s1, s2, out_dir)

# ── 11. Results tables ───────────────────────────────────────────────
print_study2_combined_table(s1, s2)
