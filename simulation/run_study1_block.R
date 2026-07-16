# =====================================================================
# Study 1: N=20 block simulation — correct spec + misspecification
#
# Network: N=20 units as five independent 4-unit blocks.
# Ground truth: analytic closed-form via full enumeration (2^4=16 states).
# Estimator: Auto-G Gibbs with TRUE model parameters (no PL estimation).
#
# Setting 1 — correct specification
#   Data from CG(x), estimator uses CG(x) network  (x = a, b, c)
#
# Setting 2 — misspecification  (true DGP = CG(c))
#   Data always from CG(c); estimator uses CG(a), CG(b), or CG(c)
#   CG(c) row in Setting 2 reuses the Setting 1 CG(c) results.
#
# Coverage: simulation-only normal CI using empirical MC SD across replications.
#
# Outputs  (results/study1_N{N}_S{S}/)
#   study1_summary.csv       one row per scenario x estimand, both settings
#   study1_raw_s1.csv        per-rep point estimates, Setting 1
#   study1_raw_s2.csv        per-rep point estimates, Setting 2
#   block_networks.rds
#   reps_CGa.rds  reps_CGb.rds  reps_CGc.rds
#   study1_settings.rds
#
# Usage: Rscript simulation/run_study1_block.R
# =====================================================================

source("R/utils.R");         source("simulation/sim_utils.R")
source("R/estimation.R");    source("R/causal.R")

# ── §1. Configuration ─────────────────────────────────────────────────
cfg      <- make_study1_config()
params   <- make_parameters_N20()

N        <- cfg$N
S        <- cfg$S
K        <- cfg$K
m_star   <- cfg$m_star
alpha    <- cfg$alpha
seed     <- cfg$seed
dgp_burn     <- cfg$dgp_burn
dgp_parallel <- if (!is.null(cfg$dgp_parallel)) cfg$dgp_parallel else TRUE

scenarios <- c("CGa", "CGb", "CGc")
n_cores   <- max(1L, parallel::detectCores() - 1L)

format_tag_value <- function(x) {
  x <- format(unname(x), scientific = FALSE, trim = TRUE)
  x <- sub("^-", "m", x)
  x <- gsub("\\.", "p", x)
  x
}

make_study1_out_dir <- function(N, S, config, params) {
  key_params <- c(
    b3  = "beta3",
    b4  = "beta4",
    th  = "theta",
    g1  = "gamma1",
    eta = "eta",
    rho = "rho",
    om  = "omega"
  )
  param_tag <- paste(
    vapply(names(key_params), function(tag) {
      nm <- key_params[[tag]]
      sprintf("%s_%s", tag, format_tag_value(params[[nm]]))
    }, character(1L)),
    collapse = "_"
  )
  file.path(
    "results",
    sprintf("study1_N%d_S%d_K%d_burn%d_%s", N, S, config$K, config$dgp_burn, param_tag)
  )
}

out_dir <- make_study1_out_dir(N, S, cfg, params)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_run_settings <- function(out_dir, params, config) {
  lines <- c(
    "# MultiNetCausalRcpp Study 1 (N=20) run settings",
    sprintf("timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("out_dir: %s", normalizePath(out_dir, winslash = "/", mustWork = FALSE)),
    "",
    "[configuration]",
    sprintf("N: %s", config$N),
    sprintf("S: %s", config$S),
    sprintf("K: %s", config$K),
    sprintf("m_star: %s", config$m_star),
    sprintf("alpha: %s", config$alpha),
    sprintf("seed: %s", config$seed),
    sprintf("dgp_burn: %s", config$dgp_burn),
    sprintf("dgp_parallel: %s", if (!is.null(config$dgp_parallel)) config$dgp_parallel else TRUE),
    "",
    "[true_parameters]"
  )
  param_lines <- vapply(names(params), function(nm) {
    sprintf("%s: %s", nm, format(unname(params[[nm]]), scientific = FALSE, trim = TRUE))
  }, character(1L))
  lines <- c(lines, param_lines)
  writeLines(lines, file.path(out_dir, "run_settings.txt"))
}

save_run_settings(out_dir, params, cfg)
cat("Run settings saved to", file.path(out_dir, "run_settings.txt"), "\n")

start_time <- Sys.time()
cat("=== Study 1 (N=20): block simulation — correct spec + misspecification ===\n")
cat(sprintf("N=%d  S=%d  K=%d  m_star=%d  dgp_burn=%d  n_cores=%d\n",
            N, S, K, m_star, dgp_burn, n_cores))
cat("Output directory:", out_dir, "\n")

# ── §2. Build networks and compute analytic exact truths ──────────────
cat("\n=== Building networks and computing exact analytic truths ===\n")
networks <- setNames(
  lapply(scenarios, make_repeated_block_network, N = N),
  scenarios
)
truths <- setNames(
  lapply(scenarios, exact_truth_4unit, params = params, alpha = alpha),
  scenarios
)
for (sc in scenarios) {
  t <- truths[[sc]]
  cat(sprintf("  %s  exact truth:  DE=%+.4f  IE=%+.4f  ATE=%+.4f\n",
              sc, t$DE, t$IE, t$ATE))
}
saveRDS(networks, file.path(out_dir, "block_networks.rds"))

# block_id identifies the repeated 4-unit blocks and is used to report n_blocks.
block_id <- networks$CGc$block_id   # same for all three networks
n_blocks <- length(unique(block_id))
cat(sprintf("  Block structure: %d blocks x %d units\n", n_blocks, N / n_blocks))

# ── §3. Quick DGP timing + generate independent replications ──────────
cat(sprintf("\n=== Generating independent replications  (S=%d per DGP) ===\n", S))
cat("Each replication: fresh L/A/Y ~ Bernoulli(0.4), burn-in ",
    dgp_burn, " sweeps, save sweep ", dgp_burn + 1L, ".\n", sep = "")
cat("\nCalibrating Study 1 DGP generation time...\n")
quick_dgp_timing(networks$CGc, params, S_total = S, burn = dgp_burn,
                 n_cores = n_cores, n_groups = length(scenarios))

reps <- setNames(vector("list", length(scenarios)), scenarios)
for (j in seq_along(scenarios)) {
  sc <- scenarios[[j]]
  reps[[sc]] <- generate_replications(
    networks[[sc]], params, S = S,
    seed = seed + j, burn = dgp_burn,
    parallel = dgp_parallel, n_cores = n_cores,
    label = paste("Study 1 DGP", sc)
  )
  saveRDS(reps[[sc]], file.path(out_dir, paste0("reps_", sc, ".rds")))
}

# ── §4. Parallel cluster setup ────────────────────────────────────────
cl <- parallel::makeCluster(n_cores)
wd <- getwd()
parallel::clusterExport(cl, "wd", envir = environment())
invisible(parallel::clusterCall(cl, function() {
  setwd(wd)
  source("R/utils.R"); source("simulation/sim_utils.R")
  source("R/estimation.R"); source("R/causal.R")
}))
cat(sprintf("Parallel cluster ready: %d workers\n", n_cores))

# ── §5. Setting 1 — correct specification ────────────────────────────

cat("\n=== Setting 1: Correct specification ===\n")
eff_s1 <- setNames(vector("list", length(scenarios)), scenarios)
for (j in seq_along(scenarios)) {
  sc <- scenarios[[j]]
  cat(sprintf("  CG(%s) data + CG(%s) estimator ...", tolower(sub("CG","",sc)), tolower(sub("CG","",sc))))
  eff_s1[[sc]] <- run_causal_estimation(
    reps = reps[[sc]], network = networks[[sc]], params_source = params,
    alpha = alpha, K = K, burn = m_star, cl = cl,
    seed_offset = seed + 100000L * j,
    return_unit_level = TRUE
  )
  cat(" done\n")
}

# ── §6. Setting 2 — misspecification (true DGP = CG(c)) ──────────────
# CG(c) data throughout; estimator varies.
# CG(c) row reuses eff_s1$CGc (same data + same network = identical).
cat("\n=== Setting 2: Misspecification (true = CG(c)) ===\n")
eff_s2 <- list()
for (j in seq_along(c("CGa", "CGb"))) {
  sc <- c("CGa", "CGb")[[j]]
  cat(sprintf("  CG(c) data + CG(%s) estimator (wrong network) ...", tolower(sub("CG","",sc))))
  eff_s2[[sc]] <- run_causal_estimation(
    reps = reps$CGc, network = networks[[sc]], params_source = params,
    alpha = alpha, K = K, burn = m_star, cl = cl,
    seed_offset = seed + 300000L * j,
    return_unit_level = TRUE
  )
  cat(" done\n")
}

eff_s2$CGc <- eff_s1$CGc

if (!is.null(cl)) parallel::stopCluster(cl)

# ── §7. Save raw per-rep results ──────────────────────────────────────
save_raw <- function(eff_list) {
  do.call(rbind, lapply(names(eff_list), function(sc) {
    do.call(rbind, lapply(seq_along(eff_list[[sc]]), function(s) {
      e <- eff_list[[sc]][[s]]
      data.frame(scenario = sc, rep = s,
                 DE = e$DE, IE = e$IE, ATE = e$ATE,
                 stringsAsFactors = FALSE)
    }))
  }))
}
raw_s1 <- save_raw(eff_s1)
raw_s2 <- save_raw(eff_s2)
write.csv(raw_s1, file.path(out_dir, "study1_raw_s1.csv"), row.names = FALSE)
write.csv(raw_s2, file.path(out_dir, "study1_raw_s2.csv"), row.names = FALSE)

# ── §8. Summary table with MC-SD normal CI coverage ───────────────────
cat("\n=== Computing summary (MC-SD normal CI coverage) ===\n")
summary_df <- summary_study1_ci(
  eff_s1    = eff_s1,
  eff_s2    = eff_s2,
  truths    = truths,
  truth_cgc = truths$CGc,
  ci_method = "mc",
  ci_level  = 0.95
)
write.csv(summary_df, file.path(out_dir, "study1_summary.csv"), row.names = FALSE)

saveRDS(list(
  N = N, S = S, scenarios = scenarios, alpha = alpha,
  K = K, m_star = m_star,
  dgp_burn = dgp_burn, dgp_parallel = dgp_parallel,
  seed = seed, params = params,
  truths = truths, n_blocks = n_blocks
), file.path(out_dir, "study1_settings.rds"))

# ── §9. Print results tables (two separate tables) ────────────────────
cat(sprintf("\n  N=%d (%dx4-unit blocks)  S=%d  alpha=%.2f  K=%d  burn=%d\n",
            N, n_blocks, S, alpha, K, m_star))
print_study1_table(summary_df, ci_label = "MC-SD CI",
                    title = "Study 1: N=20, S=500, true parameters")

end_time <- Sys.time()
cat(sprintf("\nRuntime: %.2f minutes\n",
            as.numeric(difftime(end_time, start_time, units = "mins"))))
cat("Saved to:", out_dir, "\n")
cat("  study1_summary.csv\n  study1_raw_s1.csv\n  study1_raw_s2.csv\n")
cat("  block_networks.rds\n")
for (sc in scenarios) cat(sprintf("  reps_%s.rds\n", sc))

