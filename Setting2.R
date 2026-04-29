###############################################################################
# SETTING 2 — MISSPECIFICATION STUDY, N = 20
# Truth generated from Chain Graph (c) (5 independent 4-unit blocks);
# Gibbs samplers fit CG(a) and CG(b).
###############################################################################

library(parallel)

## ============================================================================
## NETWORK SIZE
## ============================================================================

n_blocks   <- 5L
block_size <- 4L
N          <- n_blocks * block_size   # 20
in_block   <- matrix(1:N, nrow = block_size)  # columns = block IDs
local_pos  <- rep(1:block_size, times = n_blocks)

## ============================================================================
## SHARED PARAMETERS (paper Section 4.1)
## ============================================================================

beta0 <- -1.0; beta1 <- 0.5; beta2 <- 0.3; beta3 <- 0.85; beta4 <- 0.55
theta <- 0.42
eta   <- -0.5
omega <-  0.42
alpha <- 0.5

n_iter  <- 2000L   
burn_in <- 500L
n_runs  <- 200L
K_alloc <- 50L    

logistic <- function(x) 1 / (1 + exp(-x))

# All 2^4 block configs (used only for exact-truth computation)
block_cfg <- as.matrix(expand.grid(rep(list(0:1), block_size)))
colnames(block_cfg) <- paste0("v", 1:block_size)
other_cfg <- as.matrix(expand.grid(rep(list(0:1), block_size - 1)))

## ============================================================================
## WITHIN-BLOCK NETWORK SPECIFICATIONS
## ============================================================================

specs <- list(
  CGa = list(
    nbr_dep   = list(2L, c(1L, 3L), c(2L, 4L), 3L),
    edges_dep = list(c(1L, 2L), c(2L, 3L), c(3L, 4L)),
    nbr_int   = list(2L, c(1L, 3L), c(2L, 4L), 3L)
  ),
  CGb = list(
    nbr_dep   = list(2L, 1L, 4L, 3L),
    edges_dep = list(c(1L, 2L), c(3L, 4L)),
    nbr_int   = list(2L, 1L, 4L, 3L)
  ),
  CGc = list(
    nbr_dep   = list(2L, c(1L, 3L), c(2L, 4L), 3L),
    edges_dep = list(c(1L, 2L), c(2L, 3L), c(3L, 4L)),
    nbr_int   = list(2L, 1L, 4L, 3L)
  )
)

build_weights <- function(spec) {
  W_int <- matrix(0, block_size, block_size)
  for (i in 1:block_size) {
    nb <- spec$nbr_int[[i]]
    if (length(nb) > 0) W_int[i, nb] <- 1 / length(nb)
  }
  wY <- matrix(0, block_size, block_size)
  for (e in spec$edges_dep) {
    wY[e[1], e[2]] <- 0.5
    wY[e[2], e[1]] <- 0.5
  }
  list(W_int = W_int, wY = wY)
}

## ============================================================================
## EXACT (ANALYTICAL) TRUTH — PER BLOCK (unchanged from enumeration version)
## ============================================================================

make_exact_funcs <- function(spec, weights) {
  W_int <- weights$W_int; wY <- weights$wY

  U_L <- function(L) {
    s <- eta * sum(L)
    for (e in spec$edges_dep) s <- s + omega * L[e[1]] * L[e[2]]
    s
  }
  exact_fL <- function() {
    u <- apply(block_cfg, 1, U_L); u <- u - max(u)
    p <- exp(u); p / sum(p)
  }
  G_i <- function(i, a, l) {
    nb_i <- spec$nbr_int[[i]]
    beta0 + beta1 * a[i] + beta2 * l[i] +
      beta3 * sum(W_int[i, nb_i] * a[nb_i]) +
      beta4 * sum(W_int[i, nb_i] * l[nb_i])
  }
  U_Y <- function(y, a, l) {
    s <- 0
    for (i in 1:block_size) s <- s + y[i] * G_i(i, a, l)
    for (e in spec$edges_dep) {
      i <- e[1]; j <- e[2]
      s <- s + y[i] * y[j] * theta * wY[i, j]
    }
    s
  }
  exact_marg_Y <- function(i, a, l) {
    u <- apply(block_cfg, 1, function(y) U_Y(y, a, l))
    u <- u - max(u)
    p <- exp(u); p <- p / sum(p)
    sum(p[block_cfg[, i] == 1])
  }
  exact_psi <- function(i, a) {
    pL <- exact_fL(); out <- 0
    for (k in seq_len(nrow(block_cfg))) {
      l <- as.integer(block_cfg[k, ])
      out <- out + pL[k] * exact_marg_Y(i, a, l)
    }
    out
  }
  list(exact_psi = exact_psi)
}

compute_truth <- function(exact_psi) {
  psi_0 <- sapply(1:block_size, function(i) exact_psi(i, rep(0L, block_size)))
  DE_i <- numeric(block_size); IE_i <- numeric(block_size)
  for (i in 1:block_size) {
    for (k in seq_len(nrow(other_cfg))) {
      a_mi <- as.integer(other_cfg[k, ])
      pa <- alpha^sum(a_mi) * (1 - alpha)^(block_size - 1 - sum(a_mi))
      a1 <- integer(block_size); a0 <- integer(block_size); idx <- 1L
      for (j in 1:block_size) {
        if (j != i) { a1[j] <- a_mi[idx]; a0[j] <- a_mi[idx]; idx <- idx + 1L }
      }
      a1[i] <- 1L; a0[i] <- 0L
      psi_1 <- exact_psi(i, a1); psi_0c <- exact_psi(i, a0)
      DE_i[i] <- DE_i[i] + pa * (psi_1 - psi_0c)
      IE_i[i] <- IE_i[i] + pa * (psi_0c - psi_0[i])
    }
  }
  ATE_i <- DE_i + IE_i
  list(DE = mean(DE_i), IE = mean(IE_i), ATE = mean(ATE_i))
}

## ============================================================================
## GIBBS SAMPLER — takes a fixed allocation, returns psi_hat for all N units
## ============================================================================

make_gibbs_psi <- function(spec, weights) {
  W_int <- weights$W_int; wY <- weights$wY

  # Global-index neighbor lists
  global_nbr_dep <- vector("list", N)
  global_nbr_int <- vector("list", N)
  for (b in 1:n_blocks) {
    for (lp in 1:block_size) {
      u  <- in_block[lp, b]
      nd <- spec$nbr_dep[[lp]]
      ni <- spec$nbr_int[[lp]]
      global_nbr_dep[[u]] <- if (length(nd) > 0) in_block[nd, b] else integer(0)
      global_nbr_int[[u]] <- if (length(ni) > 0) in_block[ni, b] else integer(0)
    }
  }

  function(a_alloc, n_iter, burn_in) {
    L <- integer(N); Y <- integer(N)
    acc <- numeric(N); kept <- 0L
    for (m in 1:n_iter) {
      for (u in 1:N) {
        lp_u      <- local_pos[u]
        nb_dep    <- global_nbr_dep[[u]]
        nb_int    <- global_nbr_int[[u]]
        nb_dep_lp <- local_pos[nb_dep]
        nb_int_lp <- local_pos[nb_int]

        # L update
        p_L <- logistic(eta + omega * sum(L[nb_dep]))
        L[u] <- as.integer(runif(1) < p_L)

        # Y update
        lin <- beta0 + beta1 * a_alloc[u] + beta2 * L[u] +
          beta3 * sum(W_int[lp_u, nb_int_lp] * a_alloc[nb_int]) +
          beta4 * sum(W_int[lp_u, nb_int_lp] * L[nb_int]) +
          sum(wY[lp_u, nb_dep_lp] * theta * Y[nb_dep])
        p_Y <- logistic(lin)
        Y[u] <- as.integer(runif(1) < p_Y)
      }
      if (m > burn_in) { acc <- acc + Y; kept <- kept + 1L }
    }
    acc / kept
  }
}

## ============================================================================
## ONE REPLICATE — SAMPLED-ALLOCATION ESTIMATOR
## ============================================================================

make_one_replicate <- function(gibbs_psi) {
  function(run_id) {
    set.seed(10000L + run_id * 100L)

    # --- Baseline: psi_hat_i(a = 0) for IE, one Gibbs run ---
    psi_zero <- gibbs_psi(integer(N), n_iter, burn_in)  # length N

    # --- For each unit i, K sampled allocations ---
    DE_i <- numeric(N)
    IE_i <- numeric(N)

    for (i in 1:N) {
      diffs_DE <- numeric(K_alloc)
      psi_0_samples <- numeric(K_alloc)

      for (k in 1:K_alloc) {
        a_mi <- as.integer(rbinom(N - 1, size = 1, prob = alpha))

        # Build a^{(k,1)} and a^{(k,0)}
        a1 <- integer(N); a0 <- integer(N)
        idx <- 1L
        for (j in 1:N) {
          if (j != i) {
            a1[j] <- a_mi[idx]; a0[j] <- a_mi[idx]
            idx <- idx + 1L
          }
        }
        a1[i] <- 1L; a0[i] <- 0L

        psi_hat_1 <- gibbs_psi(a1, n_iter, burn_in)[i]
        psi_hat_0 <- gibbs_psi(a0, n_iter, burn_in)[i]

        diffs_DE[k]      <- psi_hat_1 - psi_hat_0
        psi_0_samples[k] <- psi_hat_0
      }

      DE_i[i] <- mean(diffs_DE)
      IE_i[i] <- mean(psi_0_samples) - psi_zero[i]
    }

    ATE_i <- DE_i + IE_i
    c(DE = mean(DE_i), IE = mean(IE_i), ATE = mean(ATE_i))
  }
}

## ============================================================================
## TRUTH UNDER CG(c) — computed once
## ============================================================================

cat("==================================================\n")
cat(sprintf("=== Truth under CG(c)  (N = %d, %d x %d blocks)\n",
            N, n_blocks, block_size))
cat("==================================================\n")
truth_spec    <- specs$CGc
truth_weights <- build_weights(truth_spec)
truth_exact   <- make_exact_funcs(truth_spec, truth_weights)
truth         <- compute_truth(truth_exact$exact_psi)
cat(sprintf("  DE_true  = %+.6f\n", truth$DE))
cat(sprintf("  IE_true  = %+.6f\n", truth$IE))
cat(sprintf("  ATE_true = %+.6f\n\n", truth$ATE))

cat(sprintf("Simulation settings:\n"))
cat(sprintf("  n_iter       = %d\n", n_iter))
cat(sprintf("  burn_in      = %d\n", burn_in))
cat(sprintf("  K_alloc      = %d\n", K_alloc))
cat(sprintf("  n_runs       = %d\n", n_runs))
cat(sprintf("  Gibbs runs per replicate = %d\n",
            2L * K_alloc * N + 1L))
cat(sprintf("  Total Gibbs runs (both misspec)  = %d\n\n",
            2L * n_runs * (2L * K_alloc * N + 1L)))

## ============================================================================
## RUN A MISSPECIFIED SAMPLER
## ============================================================================

run_misspec <- function(fit_model_name, truth) {
  cat("==================================================\n")
  cat(sprintf("=== Truth = CG(c), Gibbs sampler = %s\n", fit_model_name))
  cat("==================================================\n")

  spec    <- specs[[fit_model_name]]
  weights <- build_weights(spec)
  gibbs_psi     <- make_gibbs_psi(spec, weights)
  one_replicate <- make_one_replicate(gibbs_psi)

  n_cores <- max(1L, detectCores() - 1L)
  cat(sprintf("--- Running %d replicates on %d cores ---\n", n_runs, n_cores))

  t0 <- Sys.time()
  if (.Platform$OS.type != "windows") {
    results_list <- mclapply(
      1:n_runs, one_replicate,
      mc.cores = n_cores, mc.preschedule = TRUE, mc.set.seed = FALSE
    )
  } else {
    cl <- makeCluster(n_cores)
    clusterExport(cl,
                  varlist = c("N", "n_blocks", "block_size", "in_block", "local_pos",
                              "spec", "weights",
                              "beta0", "beta1", "beta2", "beta3", "beta4",
                              "theta", "eta", "omega", "alpha",
                              "n_iter", "burn_in", "K_alloc",
                              "logistic", "gibbs_psi", "one_replicate"),
                  envir = environment())
    results_list <- parLapply(cl, 1:n_runs, one_replicate)
    stopCluster(cl)
  }
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("Done in %.1f s (%.1f min)\n\n", elapsed, elapsed / 60))

  results_mat <- do.call(rbind, results_list)

  summarise_one <- function(name, hats, true_val) {
    mean_hat <- mean(hats); bias <- mean_hat - true_val; emp_SD <- sd(hats)
    covered  <- mean((hats - 1.96 * emp_SD <= true_val) &
                     (true_val <= hats + 1.96 * emp_SD))
    data.frame(fit_model = fit_model_name, estimand = name, truth = true_val,
               mean_est = mean_hat, bias = bias,
               emp_SE = emp_SD, coverage = covered)
  }

  summary_tbl <- rbind(
    summarise_one("DE",  results_mat[, "DE"],  truth$DE),
    summarise_one("IE",  results_mat[, "IE"],  truth$IE),
    summarise_one("ATE", results_mat[, "ATE"], truth$ATE)
  )

  cat(sprintf("--- Summary (truth CG(c), fit %s) ---\n", fit_model_name))
  print(summary_tbl, row.names = FALSE, digits = 4)
  cat("\n")

  list(fit_model = fit_model_name, results = results_mat, summary = summary_tbl)
}

## ============================================================================
## MAIN: run both misspecifications
## ============================================================================

all_results <- list()
for (m in c("CGa", "CGb")) {
  all_results[[m]] <- run_misspec(m, truth)
}

cat("==================================================\n")
cat("======  SETTING 2 (N=20) COMBINED SUMMARY  =======\n")
cat("===  Truth = CG(c); rows show fit model         ===\n")
cat("==================================================\n")
combined <- do.call(rbind, lapply(all_results, function(x) x$summary))
print(combined, row.names = FALSE, digits = 4)
