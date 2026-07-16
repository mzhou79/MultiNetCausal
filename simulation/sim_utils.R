# =====================================================================
# Simulation utilities.
# This file contains all simulation/manuscript-only helpers for this repo:
# parameters, block-network helpers, reporting helpers, and plotting helpers.
# It intentionally lives under simulation/, not R/, because these are not part
# of the reusable estimation/effect engine.
# =====================================================================

# ── Parameter builders ──────────────────────────────────────────────────

make_parameters <- function() {
  list(
    beta0 = -1.5, beta1 = 0.8, beta2 = 0.3, beta3 = 0.4, beta4 = 0.3,
    theta = 0.5,
    lambda0 = -0.5, omega = 0.5,
    gamma0 = -0.2, rho = 0.4, gamma1 = 0.25, eta = 0.5
  )
}

# =====================================================================
# Parameters for Study 1 (N=20 block network simulation).
# Sourced by run_N20.R.
# =====================================================================

make_parameters_N20 <- function() {
  list(
    beta0 = -1.5, beta1 = 0.6, beta2 = 0.3, beta3 = 0.4, beta4 = 0.4,
    theta = 1.2,
    lambda0 = -0.5, omega = 1.2,
    gamma0 = -0.2, rho = 0.3, gamma1 = 0.3, eta = 1.2
  )
}

make_study1_config <- function() {
  list(
    N        = 20L,
    S        = 500L,
    K        = 50L,
    m_star   = 20L,
    alpha    = 0.5,
    seed     = 123L,
    dgp_burn = 1500L,
    dgp_parallel = TRUE
  )
}

# ── Study 1 block-network helpers ──────────────────────────────────────

binary_grid <- function(p) {
  as.matrix(expand.grid(rep(list(0:1), p)))
}

# =====================================================================
# Mo-original 4-unit block network utilities  (used by run_N20.R)
# =====================================================================

# Three chain-graph specifications on a 4-unit block.
# Convention: W_int and W_dep are both raw adjacency (not
#   row-standardised).
#   CGa: G_int = chain 1-2-3-4,  G_dep = chain 1-2-3-4
#   CGb: G_int = pairs (1-2, 3-4), G_dep = pairs (1-2, 3-4)
#   CGc: G_int = pairs (1-2, 3-4), G_dep = chain 1-2-3-4
.block_specs <- list(
  CGa = list(
    dep_edges     = matrix(c(1,2, 2,3, 3,4), ncol = 2L, byrow = TRUE),
    int_neighbors = list(c(2L), c(1L,3L), c(2L,4L), c(3L))
  ),
  CGb = list(
    dep_edges     = matrix(c(1,2, 3,4), ncol = 2L, byrow = TRUE),
    int_neighbors = list(c(2L), c(1L), c(4L), c(3L))
  ),
  CGc = list(
    dep_edges     = matrix(c(1,2, 2,3, 3,4), ncol = 2L, byrow = TRUE),
    int_neighbors = list(c(2L), c(1L), c(4L), c(3L))
  )
)

make_4unit_block_network <- function(scenario = c("CGa", "CGb", "CGc")) {
  scenario <- match.arg(scenario)
  spec     <- .block_specs[[scenario]]

  A_int <- matrix(0L, 4L, 4L)
  for (i in seq_len(4L)) {
    nb <- spec$int_neighbors[[i]]
    if (length(nb) > 0L) A_int[i, nb] <- 1L
  }

  A_dep <- matrix(0L, 4L, 4L)
  for (r in seq_len(nrow(spec$dep_edges))) {
    i <- spec$dep_edges[r, 1L]; j <- spec$dep_edges[r, 2L]
    A_dep[i, j] <- 1L; A_dep[j, i] <- 1L
  }
  stopifnot(isTRUE(all.equal(A_dep, t(A_dep))), all(diag(A_dep) == 0L))

  W_int <- A_int * 1.0
  W_dep <- A_dep * 1.0

  list(N = 4L, scenario = scenario,
       W_int = W_int, W_dep = W_dep, A_int = A_int, A_dep = A_dep,
       dep_edges = spec$dep_edges)
}

block_diagonal_repeat <- function(M, n_blocks) {
  bs  <- nrow(M)
  out <- matrix(0, bs * n_blocks, bs * n_blocks)
  for (b in seq_len(n_blocks)) {
    ix <- ((b - 1L) * bs + 1L):(b * bs)
    out[ix, ix] <- M
  }
  out
}

make_repeated_block_network <- function(scenario = c("CGa", "CGb", "CGc"), N = 20L) {
  scenario <- match.arg(scenario)
  if (N %% 4L != 0L) stop("N must be a multiple of 4 for repeated 4-unit blocks.")
  block          <- make_4unit_block_network(scenario)
  n_micro_blocks <- N / 4L
  W_int <- block_diagonal_repeat(block$W_int, n_micro_blocks)
  W_dep <- block_diagonal_repeat(block$W_dep, n_micro_blocks)
  list(
    N              = as.integer(N),
    scenario       = scenario,
    spec           = scenario,
    spec_label     = paste0(scenario, ": Mo-original repeated 4-unit blocks"),
    W_int          = W_int, W_dep = W_dep,
    A_int          = W_int, A_dep = W_dep,
    A_M            = W_int, A_A  = W_dep,
    n_micro_blocks = as.integer(n_micro_blocks),
    block_id       = rep(seq_len(n_micro_blocks), each = 4L)
  )
}

exact_truth_4unit <- function(scenario = c("CGa", "CGb", "CGc"),
                                      params,
                                      alpha = 0.5) {
  scenario <- match.arg(scenario)
  block    <- make_4unit_block_network(scenario)
  W_int    <- block$W_int
  edges    <- block$dep_edges

  params <- params_to_cpp(as_params(params))

  L_grid       <- binary_grid(4L)
  Y_grid       <- binary_grid(4L)
  A_minus_grid <- binary_grid(3L)

  U_L <- function(l) {
    params$lambda0[1] * sum(l) +
      params$nu_L[1, 1] * sum(apply(edges, 1L, function(e) l[e[1L]] * l[e[2L]]))
  }

  U_Y <- function(y, a, l) {
    node_part <- 0
    for (i in seq_len(4L)) {
      lin <- params$beta0 + params$beta1 * a[i] + params$beta2[1] * l[i] +
             params$beta3 * sum(W_int[i, ] * a) +
             params$beta4[1] * sum(W_int[i, ] * l)
      node_part <- node_part + y[i] * lin
    }
    dep_part <- params$theta *
      sum(apply(edges, 1L, function(e) y[e[1L]] * y[e[2L]]))
    node_part + dep_part
  }

  uL <- apply(L_grid, 1L, U_L)
  pL <- exp(uL - max(uL)); pL <- pL / sum(pL)

  exact_psi <- function(i, a) {
    out <- 0
    for (r in seq_len(nrow(L_grid))) {
      l  <- L_grid[r, ]
      uY <- apply(Y_grid, 1L, function(y) U_Y(y, a, l))
      pY <- exp(uY - max(uY)); pY <- pY / sum(pY)
      out <- out + pL[r] * sum(pY[Y_grid[, i] == 1L])
    }
    out
  }

  psi_zero <- vapply(seq_len(4L), function(i) exact_psi(i, rep(0L, 4L)), numeric(1L))

  DE_i <- IE_i <- numeric(4L)
  for (i in seq_len(4L)) {
    peers <- setdiff(seq_len(4L), i)
    for (r in seq_len(nrow(A_minus_grid))) {
      a_minus <- A_minus_grid[r, ]
      prob    <- alpha^sum(a_minus) * (1 - alpha)^(length(a_minus) - sum(a_minus))
      a1 <- a0 <- rep(0L, 4L)
      a1[peers] <- a_minus; a1[i] <- 1L
      a0[peers] <- a_minus; a0[i] <- 0L
      psi1 <- exact_psi(i, a1)
      psi0 <- exact_psi(i, a0)
      DE_i[i] <- DE_i[i] + prob * (psi1 - psi0)
      IE_i[i] <- IE_i[i] + prob * (psi0 - psi_zero[i])
    }
  }

  list(DE = mean(DE_i), IE = mean(IE_i), ATE = mean(DE_i + IE_i))
}

# ── Shared simulation runners/diagnostics ──────────────────────────────

quick_dgp_timing <- function(net, params, S_total, burn, n_check = 3,
                             n_cores = 1, n_groups = 3) {
  times <- numeric(n_check)
  for (i in seq_len(n_check)) {
    times[i] <- system.time(
      generate_replications(net, params, S = 1, seed = 9000L + i,
                            burn = burn, parallel = FALSE,
                            progress = FALSE)
    )["elapsed"]
  }
  cat(sprintf("\n--- DGP timing  [net: %s, burn=%d, n_check=%d] ---\n",
              net$spec, burn, n_check))
  cat(sprintf("  per rep:  %.3f s (mean)  %.3f s (min)  %.3f s (max)\n",
              mean(times), min(times), max(times)))
  cat(sprintf("  one S=%d rep set seq:  %.1f min\n",
              S_total, mean(times) * S_total / 60))
  if (n_cores > 1) {
    cat(sprintf("  one S=%d rep set par (%d cores): ~%.1f min\n",
                S_total, n_cores, mean(times) * S_total / 60 / n_cores))
  }
  cat(sprintf("  all %d DGP rep sets par estimate: ~%.1f min\n",
              n_groups, n_groups * mean(times) * S_total / 60 / max(1L, n_cores)))
  invisible(times)
}

quick_causal_check <- function(net, true_params, alpha = 0.5, K = 20,
                               burn = 5, seed = 25) {
  rep1 <- generate_replications(net, true_params, S = 1, seed = seed,
                                burn = 50)[[1]]
  est <- fit_pl(rep1$Y, rep1$A, rep1$L, as_network(net))
  net_eff <- as_network(net)
  ctrl <- list(seed_base = seed, seed_for_sampling = seed)
  eff_est <- estimate_effects(rep1$L, est, net_eff,
                              alpha = alpha, K = K, burn = burn,
                              control = ctrl)
  eff_true <- estimate_effects(rep1$L, true_params, net_eff,
                               alpha = alpha, K = K, burn = burn,
                               control = ctrl)
  cat(sprintf("\n=== Quick causal diagnostic  [net: %s, K=%d] ===\n",
              net$spec, K))
  cat(sprintf("  %-5s  %10s  %10s\n", "", "est.params", "true.params"))
  for (nm in c("DE", "IE", "ATE")) {
    cat(sprintf("  %-5s  %10.4f  %10.4f\n", nm, eff_est[[nm]], eff_true[[nm]]))
  }
  invisible(list(est = est, eff_est = eff_est, eff_true = eff_true))
}

run_causal_estimation <- function(reps, network, params_source,
                                  alpha = 0.5, K = 50, burn = 20,
                                  cl = NULL, progress = TRUE,
                                  progress_every = NULL,
                                  label = NULL,
                                  seed_offset = 0L,
                                  return_unit_level = FALSE,
                                  covariate_types = NULL) {
  force(reps); force(network); force(params_source); force(alpha)
  force(K); force(burn); force(seed_offset); force(return_unit_level)
  S_total <- length(reps)
  if (is.null(label)) label <- sprintf("causal %s", network$spec)
  if (is.null(progress_every)) progress_every <- max(1L, ceiling(S_total / 20L))
  progress_every <- max(1L, as.integer(progress_every))
  started_at <- Sys.time()
  network_eff <- as_network(network)

  report_progress <- function(done) {
    if (!isTRUE(progress)) return(invisible(NULL))
    elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
    pct <- 100 * done / S_total
    rate <- if (elapsed > 0) done / elapsed else NA_real_
    eta <- if (is.finite(rate) && rate > 0) (S_total - done) / rate else NA_real_
    cat(sprintf("[%s] %3.0f%% (%d/%d) elapsed %.1f min, ETA %.1f min\n",
                label, pct, done, S_total, elapsed / 60, eta / 60))
    flush.console()
  }

  use_single <- is.list(params_source) && "beta0" %in% names(params_source)
  worker <- function(s) {
    p <- if (use_single) params_source else params_source[[s]]
    seed_base <- as.integer(seed_offset + s)
    estimate_effects(
      L = reps[[s]]$L,
      params = p,
      network = network_eff,
      alpha = alpha,
      K = K,
      burn = burn,
      covariate_types = covariate_types,
      return_unit_level = return_unit_level,
      control = list(seed_base = seed_base, seed_for_sampling = seed_base)
    )
  }

  idx <- seq_along(reps)
  chunks <- split(idx, ceiling(seq_along(idx) / progress_every))
  out <- vector("list", S_total)
  for (chunk in chunks) {
    vals <- if (is.null(cl)) lapply(chunk, worker) else parallel::parLapply(cl, chunk, worker)
    out[chunk] <- vals
    report_progress(max(chunk))
  }
  out
}

quick_causal_timing <- function(reps, network, params_source,
                                n_check = 3, n_cores = 1,
                                alpha = 0.5, K = 50, burn = 20,
                                seed_offset = 0L,
                                covariate_types = NULL) {
  use_single <- is.list(params_source) && "beta0" %in% names(params_source)
  network_eff <- as_network(network)
  times <- numeric(n_check)
  for (i in seq_len(n_check)) {
    p <- if (use_single) params_source else params_source[[i]]
    seed_base <- as.integer(seed_offset + i)
    times[i] <- system.time(
      estimate_effects(
        L = reps[[i]]$L,
        params = p,
        network = network_eff,
        alpha = alpha,
        K = K,
        burn = burn,
        covariate_types = covariate_types,
        control = list(seed_base = seed_base, seed_for_sampling = seed_base)
      )
    )["elapsed"]
  }
  S_total <- length(reps)
  cat(sprintf("\n--- Causal timing  [net: %s, K=%d, n_check=%d] ---\n",
              network$spec, K, n_check))
  cat(sprintf("  per rep:  %.2f s (mean)  %.2f s (min)  %.2f s (max)\n",
              mean(times), min(times), max(times)))
  cat(sprintf("  S=%d seq:  %.1f min\n", S_total, mean(times) * S_total / 60))
  if (n_cores > 1) {
    cat(sprintf("  S=%d par (%d cores): ~%.1f min\n",
                S_total, n_cores, mean(times) * S_total / 60 / n_cores))
  }
  invisible(times)
}

# ── Reporting helpers ──────────────────────────────────────────────────

summary_pl_mc <- function(s1, s2, ci_level = 0.95) {
  specs     <- c("CGa", "CGb", "CGc")
  estimands <- c(DE = "DE", IE = "IE", ATE = "ATE")

  do.call(rbind, lapply(c(1, 2), function(setting) {
    slist <- if (setting == 1) s1 else s2
    do.call(rbind, lapply(specs, function(spec) {
      r <- slist[[spec]]
      do.call(rbind, lapply(names(estimands), function(key) {
        field <- estimands[[key]]
        truth <- if (setting == 1) mean(sapply(r$true, `[[`, field)) else mean(sapply(s2[["CGc"]]$true, `[[`, field))
        est   <- sapply(r$est, `[[`, field)
        ci    <- .mc_ci_summary(est, truth, level = ci_level)
        bias  <- mean(est) - truth
        data.frame(
          setting = setting,
          spec = spec,
          estimand = key,
          truth = truth,
          mean_est = mean(est),
          bias = bias,
          absolute_bias = abs(bias),
          mc_variance = ci$mc_variance,
          mc_sd = ci$mc_sd,
          ci_level = ci$ci_level,
          ci_half_width = ci$ci_half_width,
          ci_lower_mean = ci$ci_lower_mean,
          ci_upper_mean = ci$ci_upper_mean,
          coverage = ci$coverage,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
}

.mc_ci_summary <- function(est, truth, level = 0.95) {
  est <- as.numeric(est)
  mc_sd <- stats::sd(est, na.rm = TRUE)
  z <- stats::qnorm(1 - (1 - level) / 2)
  hw <- z * mc_sd
  lower <- est - hw
  upper <- est + hw
  list(
    mc_sd = mc_sd,
    mc_variance = mc_sd^2,
    ci_level = level,
    ci_half_width = hw,
    coverage = mean(lower <= truth & truth <= upper, na.rm = TRUE) * 100,
    ci_lower_mean = mean(lower, na.rm = TRUE),
    ci_upper_mean = mean(upper, na.rm = TRUE)
  )
}

.fmt_cov <- function(x) {
  if (is.na(x)) sprintf("%7s", "--") else sprintf("%6.1f%%", x)
}

.print_aligned_rows <- function(rows, align = NULL, indent = "  ", gap = "  ") {
  rows <- lapply(rows, function(x) {
    if (is.null(x)) return(NULL)
    as.character(x)
  })
  nonnull <- Filter(Negate(is.null), rows)
  if (!length(nonnull)) return(invisible(NULL))
  ncol <- max(vapply(nonnull, length, integer(1L)))
  rows <- lapply(rows, function(x) {
    if (is.null(x)) return(NULL)
    length(x) <- ncol
    x[is.na(x)] <- ""
    x
  })
  vals <- do.call(rbind, Filter(Negate(is.null), rows))
  widths <- apply(vals, 2L, function(z) max(nchar(z, type = "width"), na.rm = TRUE))
  if (is.null(align)) align <- rep("right", ncol)
  length(align) <- ncol
  align[is.na(align)] <- "right"
  fmt_cell <- function(x, w, a) {
    if (identical(a, "left")) sprintf("%-*s", w, x) else sprintf("%*s", w, x)
  }
  render <- function(x) {
    paste0(indent, paste(mapply(fmt_cell, x, widths, align, USE.NAMES = FALSE), collapse = gap))
  }
  for (i in seq_along(rows)) {
    row <- rows[[i]]
    if (is.null(row)) {
      cat("\n")
    } else if (identical(row[1], ".divider")) {
      cat(indent, strrep("-", nchar(render(rep("", ncol)), type = "width") - nchar(indent, type = "width")), "\n", sep = "")
    } else {
      cat(render(row), "\n", sep = "")
    }
  }
  invisible(widths)
}

.fmt_num <- function(x, digits = 4, signed = FALSE) {
  if (is.na(x)) return("--")
  if (signed) sprintf(paste0("%+.", digits, "f"), x) else sprintf(paste0("%.", digits, "f"), x)
}

.summary_combined <- function(s1, s2, ci_level = 0.95) {
  specs     <- c("CGa", "CGb", "CGc")
  estimands <- c(DE = "DE", IE = "IE", ATE = "ATE")
  do.call(rbind, lapply(specs, function(spec) {
    r1 <- s1[[spec]]
    r2 <- s2[[spec]]
    do.call(rbind, lapply(names(estimands), function(key) {
      field  <- estimands[[key]]
      truth1 <- mean(sapply(r1$true, `[[`, field))
      truth2 <- mean(sapply(s2[["CGc"]]$true, `[[`, field))
      pl1    <- sapply(r1$est, `[[`, field)
      pl2    <- sapply(r2$est, `[[`, field)
      ci1    <- .mc_ci_summary(pl1, truth1, level = ci_level)
      ci2    <- .mc_ci_summary(pl2, truth2, level = ci_level)
      data.frame(
        spec = spec,
        estimand = key,
        truth_s1 = truth1,
        truth_s2 = truth2,
        s1_pl = mean(pl1),
        s1_bias = mean(pl1) - truth1,
        s1_rmse = sqrt(mean((pl1 - truth1)^2)),
        s1_mc_var = ci1$mc_variance,
        s1_coverage = ci1$coverage,
        s2_pl = mean(pl2),
        s2_bias = mean(pl2) - truth2,
        s2_rmse = sqrt(mean((pl2 - truth2)^2)),
        s2_mc_var = ci2$mc_variance,
        s2_coverage = ci2$coverage,
        stringsAsFactors = FALSE
      )
    }))
  }))
}

print_study2_combined_table <- function(s1, s2, ci_level = 0.95,
    title = "Table 1: Correct Specification & Misspecification") {
  df <- .summary_combined(s1, s2, ci_level = ci_level)

  cat(sprintf("\n========== %s ==========\n", title))
  cat("  Bias = PL - Truth\n")
  cat(sprintf("  RMSE = sqrt(mean((PL - Truth)^2)); Coverage = mean_s[Truth in PL_s +/- %.2f x MC SD(PL_s)]\n",
              stats::qnorm(1 - (1 - ci_level) / 2)))
  cat("  S1 Truth = own-model oracle; S2 Truth = CG(c) oracle (true network + true params).\n")

  rows <- list(
    c("Model", "Est.", "S1 Truth", "S1 PL", "S1 Bias", "S1 RMSE", "S1 MC Var", "S1 95% CI",
      "S2 Truth", "S2 PL", "S2 Bias", "S2 RMSE", "S2 MC Var", "S2 95% CI"),
    c(".divider")
  )
  for (spec in c("CGa", "CGb", "CGc")) {
    lab <- tolower(sub("CG", "", spec))
    rows <- c(rows, list(NULL))
    for (key in c("DE", "IE", "ATE")) {
      z <- df[df$spec == spec & df$estimand == key, ]
      mdl <- if (key == "DE") sprintf("CG(%s)", lab) else ""
      rows <- c(rows, list(c(
        mdl, key,
        .fmt_num(z$truth_s1), .fmt_num(z$s1_pl), .fmt_num(z$s1_bias, signed = TRUE),
        .fmt_num(z$s1_rmse), .fmt_num(z$s1_mc_var), .fmt_cov(z$s1_coverage),
        .fmt_num(z$truth_s2), .fmt_num(z$s2_pl), .fmt_num(z$s2_bias, signed = TRUE),
        .fmt_num(z$s2_rmse), .fmt_num(z$s2_mc_var), .fmt_cov(z$s2_coverage)
      )))
    }
  }
  rows <- c(rows, list(c(".divider")))
  .print_aligned_rows(rows, align = c("left", "left", rep("right", 12)))
  invisible(df)
}

save_results <- function(s1, s2, out_dir = "results", ci_level = 0.95) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  specs     <- c("CGa", "CGb", "CGc")
  estimands <- c(DE = "DE", IE = "IE", ATE = "ATE")

  scalar_or_na <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_real_)
    as.numeric(x[[1]])
  }

  rep_rows <- lapply(c(1, 2), function(setting) {
    slist <- if (setting == 1) s1 else s2
    do.call(rbind, lapply(specs, function(spec) {
      do.call(rbind, lapply(c("est", "true"), function(ptype) {
        eff_list <- slist[[spec]][[ptype]]
        do.call(rbind, lapply(seq_along(eff_list), function(s) {
          e <- eff_list[[s]]
          data.frame(
            setting = setting,
            spec = spec,
            params_type = ptype,
            rep = s,
            DE = scalar_or_na(e$DE),
            IE = scalar_or_na(e$IE),
            ATE = scalar_or_na(e$ATE),
            stringsAsFactors = FALSE
          )
        }))
      }))
    }))
  })
  reps_df <- do.call(rbind, rep_rows)

  write.csv(reps_df, file.path(out_dir, "causal_reps.csv"), row.names = FALSE)

  sum_rows <- do.call(rbind, lapply(c(1, 2), function(setting) {
    slist <- if (setting == 1) s1 else s2
    do.call(rbind, lapply(specs, function(spec) {
      r <- slist[[spec]]
      do.call(rbind, lapply(names(estimands), function(key) {
        field <- estimands[[key]]
        truth <- if (setting == 1) {
          mean(sapply(r$true, `[[`, field))
        } else {
          mean(sapply(s2[["CGc"]]$true, `[[`, field))
        }
        est   <- sapply(r$est,  `[[`, field)
        ci    <- .mc_ci_summary(est, truth, level = ci_level)
        data.frame(
          setting = setting,
          spec = spec,
          estimand = key,
          truth = truth,
          mean_est = mean(est),
          bias = mean(est) - truth,
          mc_variance = ci$mc_variance,
          mc_sd = ci$mc_sd,
          ci_level = ci$ci_level,
          ci_half_width = ci$ci_half_width,
          ci_lower_mean = ci$ci_lower_mean,
          ci_upper_mean = ci$ci_upper_mean,
          coverage = ci$coverage,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
  write.csv(sum_rows, file.path(out_dir, "causal_summary.csv"), row.names = FALSE)

  pl_mc_rows <- summary_pl_mc(s1, s2, ci_level = ci_level)
  write.csv(pl_mc_rows, file.path(out_dir, "causal_summary_pl_mc.csv"), row.names = FALSE)

  cat("Saved causal_reps.csv, causal_summary.csv, causal_summary_pl_mc.csv to", out_dir, "\n")
  invisible(list(reps = reps_df, summary = sum_rows, pl_mc = pl_mc_rows))
}

summary_study1_ci <- function(eff_s1, eff_s2, truths, truth_cgc,
                              ci_method = "mc",
                              ci_level = 0.95) {
  if (!identical(ci_method, "mc")) {
    stop("Only ci_method = 'mc' is supported in the cleaned reporting helpers.", call. = FALSE)
  }
  estimands <- c(DE = "DE", IE = "IE", ATE = "ATE")
  scenarios <- c("CGa", "CGb", "CGc")

  do.call(rbind, lapply(seq_along(scenarios), function(si) {
    sc <- scenarios[si]
    do.call(rbind, lapply(names(estimands), function(key) {
      field <- estimands[[key]]
      r1 <- eff_s1[[sc]]
      r2 <- eff_s2[[sc]]
      ests1 <- sapply(r1, `[[`, field)
      ests2 <- sapply(r2, `[[`, field)
      t1 <- truths[[sc]][[field]]
      t2 <- truth_cgc[[field]]
      ci1 <- .mc_ci_summary(ests1, t1, level = ci_level)
      ci2 <- .mc_ci_summary(ests2, t2, level = ci_level)
      data.frame(
        scenario = sc,
        estimand = key,
        s1_truth = t1,
        s1_bias = mean(ests1) - t1,
        s1_rmse = sqrt(mean((ests1 - t1)^2)),
        s1_mc_se = ci1$mc_sd,
        s1_coverage = ci1$coverage,
        s2_truth = t2,
        s2_bias = mean(ests2) - t2,
        s2_rmse = sqrt(mean((ests2 - t2)^2)),
        s2_mc_se = ci2$mc_sd,
        s2_coverage = ci2$coverage,
        stringsAsFactors = FALSE
      )
    }))
  }))
}

print_study1_table <- function(df, ci_label = "95% CI",
                               title = "Study 1: N=20, S=500, true parameters") {
  scenarios <- c("CGa", "CGb", "CGc")
  estimands <- c("DE", "IE", "ATE")

  cat(sprintf("\n========== %s ==========\n", title))
  cat("  RMSE = sqrt(mean((est - truth)^2))\n")
  cat("  Setting 2 truth = CG(c).\n")

  rows <- list(
    c("Model", "Est.", "S1 Truth", "S1 Bias", "S1 RMSE", "S1 MC Var", paste("S1", ci_label),
      "S2 Truth", "S2 Bias", "S2 RMSE", "S2 MC Var", paste("S2", ci_label)),
    c(".divider")
  )
  for (sc in scenarios) {
    lab <- tolower(sub("CG", "", sc))
    rows <- c(rows, list(NULL))
    for (key in estimands) {
      r <- df[df$scenario == sc & df$estimand == key, ]
      mdl <- if (key == "DE") sprintf("CG(%s)", lab) else ""
      rows <- c(rows, list(c(
        mdl, key,
        .fmt_num(r$s1_truth), .fmt_num(r$s1_bias, signed = TRUE), .fmt_num(r$s1_rmse), .fmt_num(r$s1_mc_se^2, digits = 5), .fmt_cov(r$s1_coverage),
        .fmt_num(r$s2_truth), .fmt_num(r$s2_bias, signed = TRUE), .fmt_num(r$s2_rmse), .fmt_num(r$s2_mc_se^2, digits = 5), .fmt_cov(r$s2_coverage)
      )))
    }
  }
  rows <- c(rows, list(c(".divider")))
  .print_aligned_rows(rows, align = c("left", "left", rep("right", 10)))
  cat("  CG(c) Setting 2 row reuses Setting 1 results (same data and estimator network).\n")
  invisible(df)
}

# ── Plotting helpers ───────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(patchwork)
})

# ── internals ─────────────────────────────────────────────────────────

.net_layout <- function(net, seed = 20260622) {
  set.seed(seed)
  g <- if (!is.null(net$g_dep)) net$g_dep else
    igraph::graph_from_adjacency_matrix(net$A_A, mode = "undirected", diag = FALSE)
  igraph::layout_with_fr(g, niter = 1000L)
}

.make_tg <- function(A_mat, treat_vec, outcome_vec, xy) {
  g <- igraph::graph_from_adjacency_matrix(A_mat, mode = "undirected", diag = FALSE)
  tidygraph::as_tbl_graph(g) |>
    tidygraph::activate(nodes) |>
    dplyr::mutate(
      x       = xy[, 1L],
      y       = xy[, 2L],
      treat   = factor(ifelse(as.integer(treat_vec) == 1L, "Treated", "Control"),
                       levels = c("Control", "Treated")),
      outcome = factor(ifelse(as.integer(outcome_vec) == 1L, "Y = 1", "Y = 0"),
                       levels = c("Y = 0", "Y = 1"))
    )
}

.FILL  <- c("Control" = "#4393C3", "Treated" = "#D6604D")
.SHAPE <- c("Y = 0" = 21L, "Y = 1" = 22L)

.panel <- function(tg, title) {
  ggraph::ggraph(tg, layout = "manual", x = x, y = y) +
    ggraph::geom_edge_link(color = "gray60", alpha = 0.85, width = 0.4) +
    ggraph::geom_node_point(
      aes(fill = treat, shape = outcome),
      size = 1.6, stroke = 0.22, color = "gray25"
    ) +
    scale_fill_manual(name = "Treatment", values = .FILL,
                      guide = guide_legend(
                        override.aes = list(shape = 21L, size = 3.5,
                                            stroke = 0.4, color = "gray25"))) +
    scale_shape_manual(name = "Outcome", values = .SHAPE,
                       guide = guide_legend(
                         override.aes = list(fill = "gray60", size = 3.5,
                                             stroke = 0.4, color = "gray25"))) +
    labs(title = title) +
    ggraph::theme_graph(base_family = "", base_size = 11L) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 11.5, face = "bold"),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(4, 8, 4, 8)
    )
}

.assemble <- function(p1, p2, p3, caption) {
  fig <- (p1 + p2 + p3) +
    plot_layout(guides = "collect", ncol = 3L) &
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.text      = element_text(size = 9.5),
      legend.key.size  = unit(0.45, "cm"),
      legend.spacing.x = unit(0.4, "cm"),
      legend.margin    = margin(t = 2)
    )
  fig + plot_annotation(
    caption = caption,
    theme   = theme(
      plot.caption    = element_text(size = 8, hjust = 0, face = "italic",
                                     margin = margin(t = 6)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )
}

# ── exported functions ─────────────────────────────────────────────────

# Plot three panels for a network object from make_random_network().
# Displays by default; pass outfile to save.
plot_networks <- function(
    net, A, Y,
    titles  = c("Interference Network", "Auxiliary Network", "Augment Network"),
    caption = NULL,
    outfile = NULL,
    width = 13, height = 5.2, dpi = 300,
    layout_seed = 20260622
) {
  xy  <- .net_layout(net, seed = layout_seed)
  trt <- as.integer(A)
  out <- as.integer(Y)
  # if (is.null(caption))
  #   caption <- sprintf(
  #     "N = %d.  Color: blue = Control, red = Treated.  Shape: circle = Y=0, square = Y=1.",
  #     net$N
  #   )
  fig <- .assemble(
    .panel(.make_tg(net$A_M, trt, out, xy), titles[1L]),
    .panel(.make_tg(net$A_X, trt, out, xy), titles[2L]),
    .panel(.make_tg(net$A_A, trt, out, xy), titles[3L]),
    caption
  )
  if (!is.null(outfile)) {
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    ggsave(outfile, fig, width = width, height = height, dpi = dpi, bg = "white")
    message("Saved: ", outfile)
  }
  print(fig)
  invisible(fig)
}

# =====================================================================
# Simulation-only DGP / synthetic-data helpers.
# These are for synthetic-data generation and manuscript studies;
# practitioners analyzing observed data do NOT need this section.
# =====================================================================

# =====================================================================
# DGP: data-generating process (graph + joint chain over (L, A, Y)).
# Network-agnostic. Used for both simulation and specification studies.
# =====================================================================

# Erdos-Renyi two-network generator.
#   G_M: main/interference graph generated directly as an Erdos-Renyi random
#        graph with target mean degree int_deg.
#   G_X: auxiliary graph generated directly as an independent Erdos-Renyi
#        random graph with target mean degree dep_deg.
#   G_A: augmented/dependence graph, G_A = G_M union G_X.
#
# Here dep_deg means the requested auxiliary random-graph degree.  Because
# G_M and G_X are generated independently and then unioned, overlapping edges
# are counted once in G_A, so the final augmented/dependence degree is close to
# int_deg + dep_deg minus the finite-sample overlap.
#
# The returned object defaults to CG(c): W_int = G_M, W_dep = G_A, but it
# also stores raw adjacency matrices so CG(a)/CG(b)/CG(c) working networks
# can be built on the SAME underlying base graph.
make_random_network <- function(N, dep_deg, int_deg, seed = 12345) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("igraph is required: install.packages('igraph')")
  if (N < 2L) stop("N must be at least 2.")
  if (dep_deg < 0 || int_deg < 0) stop("dep_deg and int_deg must be non-negative.")
  set.seed(seed)

  # G_M: interference random graph with target mean degree int_deg.
  g_int <- igraph::sample_gnp(N, min(1, int_deg / (N - 1L)))

  # G_X: auxiliary random graph with target mean degree dep_deg.
  g_aux <- igraph::sample_gnp(N, min(1, dep_deg / (N - 1L)))

  A_M <- as.matrix(igraph::as_adjacency_matrix(g_int, sparse = FALSE))
  A_X_raw <- as.matrix(igraph::as_adjacency_matrix(g_aux, sparse = FALSE))

  # Store A_X as the auxiliary-only contribution beyond interference so that
  # A_A = A_M ∪ A_X and plots can show the added auxiliary edges cleanly.
  A_X <- pmax(A_X_raw - A_M, 0)
  A_A <- pmax(A_M, A_X_raw)

  g_aux_only <- igraph::graph_from_adjacency_matrix(A_X, mode = "undirected", diag = FALSE)
  g_dep <- igraph::graph_from_adjacency_matrix(A_A, mode = "undirected", diag = FALSE)

  net <- list(
    N = N, A_M = A_M, A_X = A_X, A_A = A_A,
    g_dep = g_dep, g_int = g_int, g_aux = g_aux_only,
    requested_aux_deg = dep_deg,
    requested_int_deg = int_deg,
    base_int_deg      = mean(rowSums(A_M)),
    base_aux_deg      = mean(rowSums(A_X_raw)),
    base_aux_only_deg = mean(rowSums(A_X)),
    base_dep_deg      = mean(rowSums(A_A)),
    frac_isolated_int = mean(rowSums(A_M) == 0)
  )
  make_network_spec(net, "CGc")
}

# Build the network specification (truth or working model) from a base_net
# produced by make_random_network.
#
#   CG(a): G_int = G_A, G_dep = G_A  (broad/over-interference working model)
#   CG(b): G_int = G_M, G_dep = G_M  (narrow/under-dependence working model)
#   CG(c): G_int = G_M, G_dep = G_A  (separated main + augmented model)
make_network_spec <- function(base_net, spec = c("CGa", "CGb", "CGc")) {
  spec <- match.arg(spec)
  if (is.null(base_net$A_M) || is.null(base_net$A_A)) {
    stop("base_net must contain raw A_M and A_A adjacency matrices.")
  }

  if (spec == "CGa") {
    A_int <- base_net$A_A
    A_dep <- base_net$A_A
    label <- "CG(a): G_int=G_A, G_dep=G_A"
  } else if (spec == "CGb") {
    A_int <- base_net$A_M
    A_dep <- base_net$A_M
    label <- "CG(b): G_int=G_M, G_dep=G_M"
  } else {
    A_int <- base_net$A_M
    A_dep <- base_net$A_A
    label <- "CG(c): G_int=G_M, G_dep=G_A"
  }

  # Both W_int and W_dep use raw adjacency (not row-standardized).
  W_int <- A_int * 1.0
  W_dep <- A_dep * 1.0

  list(N = base_net$N,
       spec = spec,
       spec_label = label,
       A_M = base_net$A_M,
       A_X = if (!is.null(base_net$A_X)) base_net$A_X else pmax(base_net$A_A - base_net$A_M, 0),
       A_A = base_net$A_A,
       A_int = A_int,
       A_dep = A_dep,
       W_int = W_int,
       W_dep = W_dep,
       int_deg = mean(rowSums(A_int)),
       dep_deg = mean(rowSums(A_dep)),
       requested_int_deg = if (!is.null(base_net$requested_int_deg)) base_net$requested_int_deg else NA_real_,
       requested_aux_deg = if (!is.null(base_net$requested_aux_deg)) base_net$requested_aux_deg else NA_real_,
       base_int_deg = mean(rowSums(base_net$A_M > 0)),
       base_aux_deg = if (!is.null(base_net$base_aux_deg)) base_net$base_aux_deg else NA_real_,
       base_aux_only_deg = if (!is.null(base_net$base_aux_only_deg)) base_net$base_aux_only_deg else NA_real_,
       base_dep_deg = mean(rowSums(base_net$A_A > 0)),
       frac_isolated_int = mean(rowSums(A_int) == 0))
}

# One joint Gibbs sweep over (L, A, Y) using the chain-graph conditionals.
# Interference terms use W_int; dependence terms use W_dep.
#   L_i: covariate dependence only, lambda0 + omega (W_dep L)_i.
#   A_i: covariate interference plus observed-treatment dependence,
#        gamma0 + gamma1 (W_int L)_i + eta (W_dep A)_i.
one_joint_sweep <- function(L, A, Y, W_int, W_dep, params) {
  N <- length(L)
  for (i in seq_len(N)) {
    lin_L <- params$lambda0 + params$omega * sum(W_dep[i, ] * L)
    L[i] <- rbinom(1, 1, logistic(lin_L))
    lin_A <- params$gamma0 +
      (if (!is.null(params$rho)) params$rho else 0) * L[i] +
      params$gamma1 * sum(W_int[i, ] * L) +
      params$eta * sum(W_dep[i, ] * A)
    A[i] <- rbinom(1, 1, logistic(lin_A))
    lin_Y <- params$beta0 + params$beta1 * A[i] + params$beta2 * L[i] +
      params$beta3 * sum(W_int[i, ] * A) + params$beta4 * sum(W_int[i, ] * L) +
      params$theta * sum(W_dep[i, ] * Y)
    Y[i] <- rbinom(1, 1, logistic(lin_Y))
  }
  list(L = L, A = A, Y = Y)
}

# Generate S approximately independent retained replications.
# Each replication starts a fresh joint Gibbs chain from an independent
# Bernoulli(init_prob) initialization for L, A, and Y, discards `burn` sweeps,
# and saves the next state.  With burn = 1500, this saves sweep 1501.
#
# The default path uses the C++ Gibbs backend when available.
generate_replications <- function(network, params, S, seed = 1,
                                   burn = 1500,
                                   init_prob = 0.4,
                                   parallel = FALSE, n_cores = NULL,
                                   progress = TRUE, progress_every = NULL,
                                   label = NULL, use_cpp = TRUE) {
  if (S < 1L) stop("S must be at least 1.")
  if (burn < 0L) stop("burn must be non-negative.")
  if (init_prob < 0 || init_prob > 1) stop("init_prob must be in [0, 1].")

  N <- network$N
  rep_seeds <- seed + seq_len(S) - 1L
  if (is.null(label)) label <- sprintf("DGP %s", network$spec)
  if (is.null(progress_every)) progress_every <- max(1L, ceiling(S / 20L))
  progress_every <- max(1L, as.integer(progress_every))
  started_at <- Sys.time()

  report_progress <- function(done) {
    if (!isTRUE(progress)) return(invisible(NULL))
    elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
    pct <- 100 * done / S
    rate <- if (elapsed > 0) done / elapsed else NA_real_
    eta <- if (is.finite(rate) && rate > 0) (S - done) / rate else NA_real_
    cat(sprintf("[%s] %3.0f%% (%d/%d) elapsed %.1f min, ETA %.1f min\n",
                label, pct, done, S, elapsed / 60, eta / 60))
    flush.console()
  }

  cpp_available <- FALSE
  if (isTRUE(use_cpp) && exists("load_gibbs_cpp", mode = "function", inherits = TRUE)) {
    cpp_available <- tryCatch({ load_gibbs_cpp(); exists("generate_replications_cpp", mode = "function", inherits = TRUE) },
                              error = function(e) FALSE)
  }

  run_cpp_chunk <- function(chunk) {
    generate_replications_cpp(network$W_int, network$W_dep, params,
                              S = length(chunk), seed = rep_seeds[[min(chunk)]],
                              burn = burn, init_prob = init_prob)
  }

  if (cpp_available) {
    idx <- seq_len(S)
    chunks <- split(idx, ceiling(seq_along(idx) / progress_every))
    out <- vector("list", S)

    if (!isTRUE(parallel) || S == 1L) {
      for (chunk in chunks) {
        out[chunk] <- run_cpp_chunk(chunk)
        report_progress(max(chunk))
      }
      return(out)
    }

    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("parallel package unavailable; falling back to sequential C++ DGP generation.")
      for (chunk in chunks) {
        out[chunk] <- run_cpp_chunk(chunk)
        report_progress(max(chunk))
      }
      return(out)
    }
    if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
    n_cores <- max(1L, min(as.integer(n_cores), length(chunks)))
    if (n_cores == 1L) {
      for (chunk in chunks) {
        out[chunk] <- run_cpp_chunk(chunk)
        report_progress(max(chunk))
      }
      return(out)
    }

    cpp_file <- file.path(getwd(), "src", "gibbs.cpp")
    W_int <- network$W_int
    W_dep <- network$W_dep
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl,
      varlist = c("cpp_file", "W_int", "W_dep", "params", "burn", "init_prob", "rep_seeds"),
      envir = environment())
    vals_by_chunk <- parallel::parLapply(cl, chunks, function(chunk) {
      if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required for C++ DGP generation.")
      Rcpp::sourceCpp(cpp_file, rebuild = FALSE, verbose = FALSE)
      generate_replications_cpp(W_int, W_dep, params,
                                S = length(chunk), seed = rep_seeds[[min(chunk)]],
                                burn = burn, init_prob = init_prob)
    })
    for (k in seq_along(chunks)) {
      out[chunks[[k]]] <- vals_by_chunk[[k]]
      report_progress(max(chunks[[k]]))
    }
    return(out)
  }

  # R fallback path.  Kept for transparency and for environments without Rcpp.
  logistic_fn <- get("logistic", mode = "function")
  logistic <- logistic_fn
  joint_sweep <- one_joint_sweep
  environment(joint_sweep) <- environment()

  draw_one <- function(s) {
    set.seed(rep_seeds[[s]])
    L <- rbinom(N, 1, init_prob)
    A <- rbinom(N, 1, init_prob)
    Y <- rbinom(N, 1, init_prob)

    for (m in seq_len(burn + 1L)) {
      state <- joint_sweep(L, A, Y, network$W_int, network$W_dep, params)
      L <- state$L; A <- state$A; Y <- state$Y
    }
    list(L = L, A = A, Y = Y)
  }

  out <- vector("list", S)
  for (s in seq_len(S)) {
    out[[s]] <- draw_one(s)
    if (s == 1L || s == S || s %% progress_every == 0L) report_progress(s)
  }
  out
}
