# Synthetic full-cohort / outcome-subnetwork example.
# One continuous and four binary covariates; all have network dependence.

source("R/utils.R")
source("R/estimation.R")
source("R/causal.R")
source("R/bootstrap.R")

# ---------------------------------------------------------------------
# 1. Full cohort and Study 2 random-network design.
# ---------------------------------------------------------------------
N <- 200L
N_outcome <- 150L
int_deg <- 1
aux_deg <- 4
network_seed <- 123L
outcome_seed <- 456L

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("igraph is required for this example.")
}

# Same construction as Study 2:
#   G_int: interference graph
#   G_aux: independently generated auxiliary graph
#   G_aug: union of G_int and G_aux; used as the dependence graph
set.seed(network_seed)
g_int <- igraph::sample_gnp(N, int_deg / (N - 1L), directed = FALSE, loops = FALSE)
g_aux <- igraph::sample_gnp(N, aux_deg / (N - 1L), directed = FALSE, loops = FALSE)
W_int_full <- as.matrix(igraph::as_adjacency_matrix(g_int, sparse = FALSE))
W_aux_raw_full <- as.matrix(igraph::as_adjacency_matrix(g_aux, sparse = FALSE))
W_aux_full <- pmax(W_aux_raw_full - W_int_full, 0)
W_aug_full <- pmax(W_int_full, W_aux_raw_full)
storage.mode(W_int_full) <- "double"
storage.mode(W_aux_full) <- "double"
storage.mode(W_aug_full) <- "double"

# Randomly select 150 of the 200 full-cohort units for the observed-outcome subnetwork.
set.seed(outcome_seed)
outcome_units <- sort(sample(seq_len(N), N_outcome, replace = FALSE))

# Induced outcome-subnetwork matrices embedded back into the full N x N index space.
W_int_Y <- matrix(0, N, N)
W_aux_Y <- matrix(0, N, N)
W_aug_Y <- matrix(0, N, N)
W_int_Y[outcome_units, outcome_units] <- W_int_full[outcome_units, outcome_units]
W_aux_Y[outcome_units, outcome_units] <- W_aux_full[outcome_units, outcome_units]
W_aug_Y[outcome_units, outcome_units] <- W_aug_full[outcome_units, outcome_units]

network <- make_network(
  W_dep_L = W_aug_full,
  W_int_A = W_int_full,
  W_dep_A = W_aug_full,
  W_int_Y = W_int_Y,
  W_dep_Y = W_aug_Y,
  outcome_units = outcome_units
)

network_summary <- data.frame(
  network = c("Interference", "Auxiliary only", "Augmented / dependence"),
  full_mean_degree = c(
    mean(rowSums(W_int_full)),
    mean(rowSums(W_aux_full)),
    mean(rowSums(W_aug_full))
  ),
  outcome_subnetwork_mean_degree = c(
    mean(rowSums(W_int_Y[outcome_units, outcome_units, drop = FALSE])),
    mean(rowSums(W_aux_Y[outcome_units, outcome_units, drop = FALSE])),
    mean(rowSums(W_aug_Y[outcome_units, outcome_units, drop = FALSE]))
  ),
  outcome_subnetwork_isolates = c(
    sum(rowSums(W_int_Y[outcome_units, outcome_units, drop = FALSE]) == 0),
    sum(rowSums(W_aux_Y[outcome_units, outcome_units, drop = FALSE]) == 0),
    sum(rowSums(W_aug_Y[outcome_units, outcome_units, drop = FALSE]) == 0)
  ),
  check.names = FALSE
)

# ---------------------------------------------------------------------
# 2. Strong-interference / weak-dependence parameters from Study 2.
# ---------------------------------------------------------------------
setting <- list(
  beta0 = -1.5, beta1 = 0.8, beta2 = 0.3,
  beta3 = 0.9, beta4 = 0.3, theta = 0.1,
  lambda0 = -0.5, omega = 0.1,
  gamma0 = -0.2, rho = 0.3, gamma1 = 0.1, eta = 0.1
)

covariate_types <- c("continuous", rep("binary", 4L))
p <- length(covariate_types)
cross_covariate <- matrix(0.1, p, p)
diag(cross_covariate) <- 0
true_params <- make_params(
  covariate_model = list(
    intercept = rep(setting$lambda0, p),
    within_node = cross_covariate,
    network_lag = diag(rep(setting$omega, p), p, p),
    sigma = rep(1, p)
  ),
  treatment_model = list(
    intercept = setting$gamma0,
    own_covariates = rep(setting$rho, p),
    network_covariates = rep(setting$gamma1, p),
    dependence = setting$eta
  ),
  outcome_model = list(
    intercept = setting$beta0,
    treatment = setting$beta1,
    own_covariates = rep(setting$beta2, p),
    network_treatment = setting$beta3,
    network_covariates = rep(setting$beta4, p),
    dependence = setting$theta
  ),
  meta = list(
    p = p,
    covariate_types = covariate_types,
    fit_status = list(ok = TRUE, issues = character(), model_status = list())
  )
)

# ---------------------------------------------------------------------
# 3. Generate one synthetic dataset from a fresh initial state.
# ---------------------------------------------------------------------
set.seed(789L)
L_start <- cbind(
  continuous = rnorm(N),
  binary_1 = rbinom(N, 1, 0.4),
  binary_2 = rbinom(N, 1, 0.4),
  binary_3 = rbinom(N, 1, 0.4),
  binary_4 = rbinom(N, 1, 0.4)
)
A_start <- rbinom(N, 1, 0.4)
Y_start <- rep(NA_real_, N)
Y_start[outcome_units] <- rbinom(N_outcome, 1, 0.4)

data <- sample_data(
  Y = Y_start,
  A = A_start,
  L = L_start,
  network = network,
  params = true_params,
  covariate_types = covariate_types,
  n_sweeps = 1501L,
  seed = 790L
)
L <- data$L
A <- data$A
Y <- data$Y

# ---------------------------------------------------------------------
# 4. Fit the observed-data models and estimate effects over 150 units.
# ---------------------------------------------------------------------
pl_fit <- fit_pl(Y, A, L, network, covariate_types)
if (!params_are_valid(pl_fit)) stop(params_problem_text(pl_fit))

effects <- estimate_effects(
  L = L,
  params = pl_fit,
  network = network,
  alpha = 0.5,
  K = 20L,
  burn = 10L,
  covariate_types = covariate_types,
  control = list(seed_base = 2001L, seed_for_sampling = 2001L)
)

# Small B keeps the submitted example quick; increase B for an analysis run.
bootstrap <- bootstrap_effects(
  Y = Y,
  A = A,
  L = L,
  network = network,
  covariate_types = covariate_types,
  B = 5L,
  bootstrap_burn = 100L,
  alpha = 0.5,
  K = 20L,
  burn = 10L,
  seed = 2001L,
  n_cores = 1L
)

cat("Full cohort N:", N, "\n")
cat("Outcome subnetwork N:", length(outcome_units), "\n\n")
cat("Realized network structure\n")
print(network_summary, row.names = FALSE)
cat("\nCovariates\n")
print(data.frame(
  covariate = colnames(L),
  type = covariate_types,
  mean = colMeans(L),
  sd = apply(L, 2L, stats::sd),
  check.names = FALSE
), row.names = FALSE)
cat("\nWithin-unit cross-covariate parameter matrix (rho_L)\n")
print(cross_covariate)
cat("\nEmpirical covariate correlation matrix\n")
print(round(stats::cor(L), 3))
cat("\nFitted within-unit cross-covariate coefficients (rho_L)\n")
print(round(pl_fit$covariate_model$within_node, 3))
cat("\nObserved-data effect estimates\n")
print(unlist(effects[c("DE", "IE", "ATE")]))
cat("\nSuccessful bootstrap replications:", bootstrap$successes, "/", bootstrap$settings$B, "\n\n")
print(bootstrap$summary, row.names = FALSE)
