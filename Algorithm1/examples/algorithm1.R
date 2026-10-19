# Complete synthetic-data example. Run from the project folder:
# source("examples/algorithm1.R")
# For real data, replace Section 2 with the prepared Y, A, L and network.

source("R/utils.R")
source("R/estimation.R")
source("R/causal.R")
source("R/bootstrap.R")

# 1. Analysis settings.
alpha <- 0.5
B_alloc <- 1L          # Fixed treatment allocations per unit.
K <- 50L              # Retained Gibbs sweeps per allocation.
burn <- 20L
B_boot <- 5L          # Small runnable illustration, not final SE/CI inference.
bootstrap_burn <- 100L
n_cores <- 1L
seed <- 2810L

# 2. Synthetic observed data and networks (no participant data).
N <- 160L
set.seed(1801L)
W_M <- matrix(rbinom(N * N, 1L, 1 / (N - 1L)), N, N)
W_M[lower.tri(W_M)] <- t(W_M)[lower.tri(W_M)]
diag(W_M) <- 0
W_X <- matrix(rbinom(N * N, 1L, 4 / (N - 1L)), N, N)
W_X[lower.tri(W_X)] <- t(W_X)[lower.tri(W_X)]
diag(W_X) <- 0
W_A <- 1 * ((W_M + W_X) > 0)
network <- make_network(W_A, W_M, W_A, W_M, W_A)

# These parameters generate the example data only; estimation below refits them.
true_params <- as_params(list(
  lambda0 = -0.5, omega = 0.3,
  gamma0 = -0.2, rho = 0.4, gamma1 = 0.25, eta = 0.3,
  beta0 = -1.5, beta1 = 0.8, beta2 = 0.3,
  beta3 = 0.5, beta4 = 0.3, theta = 0.3
))
example_data <- sample_data(
  Y = rep(0, N), A = rep(0, N), L = rep(0, N),
  network = network, params = true_params,
  n_sweeps = 200L, seed = 1802L
)
Y <- example_data$Y
A <- example_data$A
L <- example_data$L

# 3. Fit the conditional models by pseudo-likelihood.
cat("\nFitting conditional models...\n")
fitted_params <- fit_pl(Y, A, L, network)
print(params_flat_vector(fitted_params))

# 4. Estimate effects using the fitted parameters and Algorithm 1.
cat("\nEstimating DE, IE and ATE...\n")
point_estimates <- estimate_effects(
  L = L, params = fitted_params, network = network,
  alpha = alpha, K = K, burn = burn,
  return_unit_level = TRUE,
  control = list(seed_base = seed, seed_for_sampling = seed),
  B = B_alloc, Y_start = Y
)
print(point_estimates[c("DE", "IE", "ATE")])

# 5. Parametric bootstrap: generate new data, refit and re-estimate effects.
cat("\nRunning parametric bootstrap...\n")
bootstrap_results <- bootstrap_effects(
  Y = Y, A = A, L = L, network = network,
  B = B_boot, bootstrap_burn = bootstrap_burn,
  alpha = alpha, K = K, burn = burn,
  seed = seed, n_cores = n_cores, B_alloc = B_alloc
)
print(bootstrap_results$summary)
cat("\nSuccessful bootstrap replications:", bootstrap_results$successes,
    "/", B_boot, "\n")

# Inspect fitted_params, point_estimates, bootstrap_results$draws and $summary.
# The small B_boot is only for demonstrating the pipeline, not reporting intervals.
