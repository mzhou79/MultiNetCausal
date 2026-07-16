# Parametric bootstrap for DE, IE, and ATE.

sample_data <- function(Y, A, L, network, params, covariate_types = NULL,
                        n_sweeps = 500L, seed = 1L) {
  network <- as_network(network)
  L <- as_covariate_matrix(L, network$N)
  covariate_types <- resolve_covariate_types(L, covariate_types)
  params <- as_params(params, covariate_types)
  if (!params_are_valid(params)) {
    stop("Cannot generate data with invalid parameters: ", params_problem_text(params), call. = FALSE)
  }
  A <- as.numeric(A)
  Y <- as.numeric(Y)
  if (length(A) != network$N || any(!is.finite(A)) || !all(A %in% c(0, 1))) {
    stop("A must be a complete binary vector of length N.", call. = FALSE)
  }
  if (length(Y) != network$N || any(!is.finite(Y[network$outcome_units])) ||
      !all(Y[network$outcome_units] %in% c(0, 1))) {
    stop("Y must be binary and observed for every outcome_units index.", call. = FALSE)
  }
  if (n_sweeps < 1L) stop("n_sweeps must be positive.", call. = FALSE)

  set.seed(seed)
  N <- network$N
  p <- ncol(L)
  L_current <- L
  A_current <- as.integer(A)
  Y_current <- integer(N)
  Y_current[network$outcome_units] <- as.integer(Y[network$outcome_units])

  for (sweep in seq_len(as.integer(n_sweeps))) {
    dependence_L <- network$W_dep_L %*% L_current
    for (k in seq_len(p)) {
      for (i in seq_len(N)) {
        linear <- params$covariate_model$intercept[k] +
          sum(params$covariate_model$within_node[k, ] * L_current[i, ]) +
          sum(params$covariate_model$network_lag[k, ] * dependence_L[i, ])
        old_value <- L_current[i, k]
        new_value <- if (covariate_types[k] == "continuous") {
          stats::rnorm(1L, linear, params$covariate_model$sigma[k])
        } else {
          stats::rbinom(1L, 1L, logistic(linear))
        }
        L_current[i, k] <- new_value
        change <- new_value - old_value
        if (change != 0) {
          dependence_L[, k] <- dependence_L[, k] + network$W_dep_L[, i] * change
        }
      }
    }

    interference_L_A <- network$W_int_A %*% L_current
    dependence_A <- as.numeric(network$W_dep_A %*% A_current)
    for (i in seq_len(N)) {
      linear <- params$treatment_model$intercept +
        sum(params$treatment_model$own_covariates * L_current[i, ]) +
        sum(params$treatment_model$network_covariates * interference_L_A[i, ]) +
        params$treatment_model$dependence * dependence_A[i]
      old_value <- A_current[i]
      A_current[i] <- stats::rbinom(1L, 1L, logistic(linear))
      change <- A_current[i] - old_value
      if (change != 0) dependence_A <- dependence_A + network$W_dep_A[, i] * change
    }

    interference_A_Y <- as.numeric(network$W_int_Y %*% A_current)
    interference_L_Y <- network$W_int_Y %*% L_current
    dependence_Y <- as.numeric(network$W_dep_Y %*% Y_current)
    for (i in network$outcome_units) {
      linear <- params$outcome_model$intercept +
        params$outcome_model$treatment * A_current[i] +
        sum(params$outcome_model$own_covariates * L_current[i, ]) +
        params$outcome_model$network_treatment * interference_A_Y[i] +
        sum(params$outcome_model$network_covariates * interference_L_Y[i, ]) +
        params$outcome_model$dependence * dependence_Y[i]
      old_value <- Y_current[i]
      Y_current[i] <- stats::rbinom(1L, 1L, logistic(linear))
      change <- Y_current[i] - old_value
      if (change != 0) dependence_Y <- dependence_Y + network$W_dep_Y[, i] * change
    }
  }

  Y_output <- as.numeric(Y_current)
  Y_output[-network$outcome_units] <- NA_real_
  colnames(L_current) <- colnames(L)
  list(L = L_current, A = A_current, Y = Y_output)
}

.bootstrap_one <- function(b, Y, A, L, network, observed_params, covariate_types,
                           bootstrap_burn, alpha, K, burn, seed) {
  tryCatch({
    generated <- sample_data(
      Y = Y, A = A, L = L, network = network,
      params = observed_params, covariate_types = covariate_types,
      n_sweeps = bootstrap_burn, seed = seed + b
    )
    fitted <- fit_pl(generated$Y, generated$A, generated$L, network, covariate_types)
    if (!params_are_valid(fitted)) stop(params_problem_text(fitted))
    effects <- estimate_effects(
      L = generated$L,
      params = fitted,
      network = network,
      alpha = alpha,
      K = K,
      burn = burn,
      covariate_types = covariate_types,
      control = list(seed_base = seed + 100000L + b,
                     seed_for_sampling = seed + 100000L + b)
    )
    data.frame(rep = b, success = TRUE, DE = effects$DE, IE = effects$IE,
               ATE = effects$ATE, error = "", stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(rep = b, success = FALSE, DE = NA_real_, IE = NA_real_,
               ATE = NA_real_, error = conditionMessage(e), stringsAsFactors = FALSE)
  })
}

bootstrap_summary <- function(point_estimates, draws) {
  estimands <- c("DE", "IE", "ATE")
  successful <- draws[draws$success, , drop = FALSE]
  if (nrow(successful) < 2L) stop("At least two successful bootstrap replications are required.", call. = FALSE)
  estimates <- unname(point_estimates[estimands])
  standard_errors <- vapply(estimands, function(name) stats::sd(successful[[name]]), numeric(1L))
  lower <- estimates - 1.96 * standard_errors
  upper <- estimates + 1.96 * standard_errors
  display <- data.frame(
    Estimand = estimands,
    Estimate = estimates,
    `Bootstrap SE` = standard_errors,
    `95% CI` = sprintf("[%.4f, %.4f]", lower, upper),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  list(summary = display,
       ci = data.frame(Estimand = estimands, Lower = lower, Upper = upper,
                       check.names = FALSE, stringsAsFactors = FALSE))
}

bootstrap_effects <- function(Y, A, L, network, covariate_types = NULL,
                              B = 200L, bootstrap_burn = 500L,
                              alpha = 0.5, K = 50L, burn = 20L,
                              seed = 20260713L, n_cores = 1L) {
  if (B < 2L) stop("B must be at least 2.", call. = FALSE)
  if (bootstrap_burn < 1L) stop("bootstrap_burn must be positive.", call. = FALSE)
  network <- as_network(network)
  L <- as_covariate_matrix(L, network$N)
  covariate_types <- resolve_covariate_types(L, covariate_types)
  observed_params <- fit_pl(Y, A, L, network, covariate_types)
  if (!params_are_valid(observed_params)) {
    stop("Observed pseudo-likelihood fit is invalid: ", params_problem_text(observed_params), call. = FALSE)
  }
  observed_effects <- estimate_effects(
    L, observed_params, network,
    alpha = alpha, K = K, burn = burn,
    covariate_types = covariate_types,
    control = list(seed_base = seed, seed_for_sampling = seed)
  )

  n_cores <- max(1L, min(as.integer(n_cores), as.integer(B)))
  indices <- seq_len(B)
  arguments <- list(
    Y = as.numeric(Y), A = as.numeric(A), L = L, network = network,
    observed_params = observed_params, covariate_types = covariate_types,
    bootstrap_burn = as.integer(bootstrap_burn), alpha = alpha,
    K = as.integer(K), burn = as.integer(burn), seed = as.integer(seed)
  )
  if (n_cores == 1L) {
    rows <- lapply(indices, function(b) do.call(.bootstrap_one, c(list(b = b), arguments)))
  } else {
    cluster <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    root <- getwd()
    parallel::clusterExport(cluster, "root", envir = environment())
    parallel::clusterCall(cluster, function() {
      setwd(root)
      source("R/utils.R")
      source("R/estimation.R")
      source("R/causal.R")
      source("R/bootstrap.R")
      load_gibbs_cpp()
      NULL
    })
    rows <- parallel::parLapply(cluster, indices, function(b, arguments) {
      do.call(.bootstrap_one, c(list(b = b), arguments))
    }, arguments = arguments)
  }
  draws <- do.call(rbind, rows)
  point_estimates <- unlist(observed_effects[c("DE", "IE", "ATE")], use.names = TRUE)
  reported <- bootstrap_summary(point_estimates, draws)
  list(
    estimates = point_estimates,
    summary = reported$summary,
    ci = reported$ci,
    draws = draws,
    successes = sum(draws$success),
    failures = sum(!draws$success),
    fit = observed_params,
    settings = list(B = B, bootstrap_burn = bootstrap_burn, alpha = alpha,
                    K = K, burn = burn, seed = seed, n_cores = n_cores)
  )
}
