# Auto-G-style causal-effect estimation. Scalar L is treated as p = 1.

gibbs_psi <- function(L, network, params, alpha, n_iter, burn,
                      seed = NULL, unit_fixed = NULL, treatment_value = NULL,
                      covariate_types = NULL) {
  load_gibbs_cpp()
  network <- as_network(network)
  L <- as_covariate_matrix(L, network$N)
  params <- as_params(params, covariate_types)
  if (!params_are_valid(params)) {
    stop("Cannot estimate effects with invalid parameters: ", params_problem_text(params), call. = FALSE)
  }
  covariate_types <- resolve_covariate_types(
    L, covariate_types %||% params$meta$covariate_types
  )
  gibbs_psi_cpp(
    L_start = L,
    W_dep_L = network$W_dep_L,
    W_int_Y = network$W_int_Y,
    W_dep_Y = network$W_dep_Y,
    params = params_to_cpp(params),
    alpha = as.numeric(alpha),
    n_iter = as.integer(n_iter),
    burn_in = as.integer(burn),
    seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
    unit_index = if (is.null(unit_fixed)) -1L else as.integer(unit_fixed),
    treatment_value = if (is.null(treatment_value)) 0L else as.integer(treatment_value),
    cov_types = as.integer(covariate_types == "continuous"),
    outcome_units = as.integer(network$outcome_units)
  )
}

estimate_effects <- function(L, params, network, alpha = 0.5, K = 50L, burn = 20L,
                             covariate_types = NULL, return_unit_level = FALSE,
                             control = NULL) {
  control <- control %||% list()
  network <- as_network(network)
  params <- as_params(params, covariate_types)
  if (!params_are_valid(params)) {
    stop("Cannot estimate effects with invalid parameters: ", params_problem_text(params), call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("alpha must be a single value in [0, 1].", call. = FALSE)
  }
  if (K < 1L || burn < 0L) stop("K must be positive and burn must be non-negative.", call. = FALSE)
  set.seed(control$seed_for_sampling %||% 1L)
  units <- network$outcome_units
  n_iter <- as.integer(K + burn)
  seed_base <- as.integer(control$seed_base %||% 1L)

  psi_zero <- gibbs_psi(
    L, network, params, alpha = 0, n_iter = n_iter, burn = burn,
    seed = seed_base + 99999L, covariate_types = covariate_types
  )
  DE_i <- IE_i <- numeric(length(units))
  for (index in seq_along(units)) {
    unit <- units[index]
    psi_one <- gibbs_psi(
      L, network, params, alpha = alpha, n_iter = n_iter, burn = burn,
      seed = seed_base + 100000L + index,
      unit_fixed = unit, treatment_value = 1L,
      covariate_types = covariate_types
    )[unit]
    psi_alpha_zero <- gibbs_psi(
      L, network, params, alpha = alpha, n_iter = n_iter, burn = burn,
      seed = seed_base + 200000L + index,
      unit_fixed = unit, treatment_value = 0L,
      covariate_types = covariate_types
    )[unit]
    DE_i[index] <- psi_one - psi_alpha_zero
    IE_i[index] <- psi_alpha_zero - psi_zero[unit]
  }

  result <- list(
    DE = mean(DE_i),
    IE = mean(IE_i),
    ATE = mean(DE_i + IE_i),
    meta = list(outcome_units = units)
  )
  if (isTRUE(return_unit_level)) result[c("DE_i", "IE_i")] <- list(DE_i, IE_i)
  result
}
