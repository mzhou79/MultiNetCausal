###############################################################################
# SETTING 1 — CORRECT SPECIFICATION
# Auto-g-computation simulation, N = 4 units.
#
# Runs three chain graphs sequentially, each with correctly specified Gibbs:
#   CG(a): chain dependence + chain interference (1—2—3—4)
#   CG(b): pair dependence + pair interference   ({1,2}, {3,4})
#   CG(c): chain dependence + pair interference  (the proposed structure)
#
# For each model we compute:
#   - Analytical truth (closed-form ψ_i(a) by enumeration over Y and L)
#   - Gibbs-based estimates of DE/IE/ATE across 500 independent runs
#   - Bias, empirical SE, and 95% Wald-CI coverage
###############################################################################

## ============================================================================
## SHARED PARAMETERS (paper Section 4.1)
## ============================================================================

N <- 4L

# Outcome model (eq 7): β = (β0, β1, β2, β3, β4), θ
beta0 <- -1.0; beta1 <- 0.5; beta2 <- 0.3; beta3 <- 1.2; beta4 <- 0.8
theta <- 0.6
# L model (eq 6)
eta   <- -0.5
omega <-  0.6
# Allocation regime
alpha <- 0.5

n_iter  <- 5000L
burn_in <- 1000L
n_runs  <- 500L

logistic <- function(x) 1 / (1 + exp(-x))

all_cfg <- as.matrix(expand.grid(rep(list(0:1), N)))
colnames(all_cfg) <- paste0("v", 1:N)

all_A     <- as.matrix(expand.grid(rep(list(0:1), N)))
all_A_key <- apply(all_A, 1, function(r) paste0(r, collapse = ""))
other_cfg <- as.matrix(expand.grid(rep(list(0:1), N - 1)))

## ============================================================================
## NETWORK SPECIFICATIONS — three chain graphs
## ============================================================================
#
# Each spec is a list with:
#   nbr_dep   : list of dependence neighbors per unit (governs L-L, Y-Y)
#   edges_dep : undirected dependence edges (i, j with i<j)
#   nbr_int   : list of interference neighbors per unit (governs β3, β4 terms)
#
# Once specified, weights W_int and wY are derived automatically.

specs <- list(
  CGa = list(
    nbr_dep   = list(2L, c(1L, 3L), c(2L, 4L), 3L),       # chain
    edges_dep = list(c(1L, 2L), c(2L, 3L), c(3L, 4L)),
    nbr_int   = list(2L, c(1L, 3L), c(2L, 4L), 3L)        # chain interference
  ),
  CGb = list(
    nbr_dep   = list(2L, 1L, 4L, 3L),                     # pairs
    edges_dep = list(c(1L, 2L), c(3L, 4L)),
    nbr_int   = list(2L, 1L, 4L, 3L)                      # pair interference
  ),
  CGc = list(
    nbr_dep   = list(2L, c(1L, 3L), c(2L, 4L), 3L),       # chain dependence
    edges_dep = list(c(1L, 2L), c(2L, 3L), c(3L, 4L)),
    nbr_int   = list(2L, 1L, 4L, 3L)                      # pair interference
  )
)

build_weights <- function(spec) {
  W_int <- matrix(0, N, N)
  for (i in 1:N) {
    nb <- spec$nbr_int[[i]]
    if (length(nb) > 0) W_int[i, nb] <- 1 / length(nb)
  }
  wY <- matrix(0, N, N)
  for (e in spec$edges_dep) {
    wY[e[1], e[2]] <- 0.5
    wY[e[2], e[1]] <- 0.5
  }
  list(W_int = W_int, wY = wY)
}

## ============================================================================
## EXACT (ANALYTICAL) COMPONENTS — parameterized by spec
## ============================================================================

make_exact_funcs <- function(spec, weights) {
  W_int <- weights$W_int
  wY    <- weights$wY

  U_L <- function(L) {
    s <- eta * sum(L)
    for (e in spec$edges_dep) s <- s + omega * L[e[1]] * L[e[2]]
    s
  }
  exact_fL <- function() {
    u <- apply(all_cfg, 1, U_L); u <- u - max(u)
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
    for (i in 1:N) s <- s + y[i] * G_i(i, a, l)
    for (e in spec$edges_dep) {
      i <- e[1]; j <- e[2]
      s <- s + y[i] * y[j] * theta * wY[i, j]
    }
    s
  }
  exact_marg_Y <- function(i, a, l) {
    u <- apply(all_cfg, 1, function(y) U_Y(y, a, l))
    u <- u - max(u)
    p <- exp(u); p <- p / sum(p)
    sum(p[all_cfg[, i] == 1])
  }
  exact_psi <- function(i, a) {
    pL <- exact_fL(); out <- 0
    for (k in seq_len(nrow(all_cfg))) {
      l <- as.integer(all_cfg[k, ])
      out <- out + pL[k] * exact_marg_Y(i, a, l)
    }
    out
  }
  list(exact_psi = exact_psi)
}

compute_truth <- function(exact_psi) {
  psi_0 <- sapply(1:N, function(i) exact_psi(i, rep(0L, N)))
  DE_i <- numeric(N); IE_i <- numeric(N)
  for (i in 1:N) {
    for (k in seq_len(nrow(other_cfg))) {
      a_mi <- as.integer(other_cfg[k, ])
      pa <- alpha^sum(a_mi) * (1 - alpha)^(N - 1 - sum(a_mi))
      a1 <- integer(N); a0 <- integer(N); idx <- 1L
      for (j in 1:N) {
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
## GIBBS SAMPLER — parameterized by spec
## ============================================================================

make_gibbs_psi <- function(spec, weights) {
  W_int <- weights$W_int
  wY    <- weights$wY

  function(a_alloc, n_iter, burn_in, seed) {
    set.seed(seed)
    L <- integer(N); Y <- integer(N)
    acc <- numeric(N); kept <- 0L
    for (m in 1:n_iter) {
      for (i in 1:N) {
        nb_dep <- spec$nbr_dep[[i]]
        nb_int <- spec$nbr_int[[i]]

        # L update (uses dependence neighbors)
        p_L <- logistic(eta + omega * sum(L[nb_dep]))
        L[i] <- as.integer(runif(1) < p_L)

        # Y update: interference uses nbr_int; Y-Y coupling uses nbr_dep with wY
        lin <- beta0 + beta1 * a_alloc[i] + beta2 * L[i] +
          beta3 * sum(W_int[i, nb_int] * a_alloc[nb_int]) +
          beta4 * sum(W_int[i, nb_int] * L[nb_int]) +
          sum(wY[i, nb_dep] * theta * Y[nb_dep])
        p_Y <- logistic(lin)
        Y[i] <- as.integer(runif(1) < p_Y)
      }
      if (m > burn_in) { acc <- acc + Y; kept <- kept + 1L }
    }
    acc / kept
  }
}

## ============================================================================
## ONE REPLICATE
## ============================================================================

make_one_replicate <- function(gibbs_psi) {
  function(run_id) {
    base_seed <- 10000L + run_id * 100L
    psi_tbl <- matrix(NA_real_, nrow = nrow(all_A), ncol = N)
    rownames(psi_tbl) <- all_A_key
    for (k in seq_len(nrow(all_A))) {
      a_vec <- as.integer(all_A[k, ])
      psi_tbl[k, ] <- gibbs_psi(a_vec, n_iter, burn_in, seed = base_seed + k)
    }
    psi0_key <- paste0(rep(0L, N), collapse = "")
    DE_i <- numeric(N); IE_i <- numeric(N)
    for (i in 1:N) {
      psi_base <- psi_tbl[psi0_key, i]
      for (k in seq_len(nrow(other_cfg))) {
        a_mi <- as.integer(other_cfg[k, ])
        pa <- alpha^sum(a_mi) * (1 - alpha)^(N - 1 - sum(a_mi))
        a1 <- integer(N); a0 <- integer(N); idx <- 1L
        for (j in 1:N) if (j != i) { a1[j] <- a_mi[idx]; a0[j] <- a_mi[idx]; idx <- idx + 1L }
        a1[i] <- 1L; a0[i] <- 0L
        k1 <- paste0(a1, collapse = ""); k0 <- paste0(a0, collapse = "")
        DE_i[i] <- DE_i[i] + pa * (psi_tbl[k1, i] - psi_tbl[k0, i])
        IE_i[i] <- IE_i[i] + pa * (psi_tbl[k0, i] - psi_base)
      }
    }
    ATE_i <- DE_i + IE_i
    c(DE = mean(DE_i), IE = mean(IE_i), ATE = mean(ATE_i))
  }
}

## ============================================================================
## RUN A FULL SIMULATION FOR ONE MODEL
## ============================================================================

run_model <- function(model_name) {
  cat("==================================================\n")
  cat(sprintf("===  Model: %s\n", model_name))
  cat("==================================================\n")

  spec    <- specs[[model_name]]
  weights <- build_weights(spec)
  exact   <- make_exact_funcs(spec, weights)
  gibbs_psi     <- make_gibbs_psi(spec, weights)
  one_replicate <- make_one_replicate(gibbs_psi)

  cat("--- Computing analytical truth ---\n")
  truth <- compute_truth(exact$exact_psi)
  cat(sprintf("  DE_true  = %+.6f\n", truth$DE))
  cat(sprintf("  IE_true  = %+.6f\n", truth$IE))
  cat(sprintf("  ATE_true = %+.6f\n\n", truth$ATE))

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
                  varlist = c("N", "spec", "weights",
                              "beta0", "beta1", "beta2", "beta3", "beta4",
                              "theta", "eta", "omega", "alpha",
                              "n_iter", "burn_in",
                              "all_A", "all_A_key", "other_cfg",
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
    data.frame(model = model_name, estimand = name, truth = true_val,
               mean_est = mean_hat, bias = bias,
               emp_SE = emp_SD, coverage = covered)
  }

  summary_tbl <- rbind(
    summarise_one("DE",  results_mat[, "DE"],  truth$DE),
    summarise_one("IE",  results_mat[, "IE"],  truth$IE),
    summarise_one("ATE", results_mat[, "ATE"], truth$ATE)
  )

  cat(sprintf("--- Summary (%s) ---\n", model_name))
  print(summary_tbl, row.names = FALSE, digits = 4)
  cat("\n")

  list(model = model_name, truth = truth, results = results_mat,
       summary = summary_tbl)
}

## ============================================================================
## MAIN: run all three models
## ============================================================================

all_results <- list()
for (m in c("CGa", "CGb", "CGc")) {
  all_results[[m]] <- run_model(m)
}

cat("==================================================\n")
cat("=========  SETTING 1 — COMBINED SUMMARY  =========\n")
cat("==================================================\n")
combined <- do.call(rbind, lapply(all_results, function(x) x$summary))
print(combined, row.names = FALSE, digits = 4)



