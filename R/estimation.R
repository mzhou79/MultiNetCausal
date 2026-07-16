# Pseudo-likelihood estimation. Scalar L is treated as p = 1.

.fit_covariates <- function(L, covariate_types, W_dep_L) {
  p <- ncol(L)
  network_lag <- W_dep_L %*% L
  intercept <- rep(NA_real_, p)
  within_node <- matrix(0, p, p)
  network_coef <- matrix(NA_real_, p, p)
  sigma <- rep(1, p)
  statuses <- vector("list", p)
  issues <- character()

  for (k in seq_len(p)) {
    data <- data.frame(Lk = L[, k])
    other <- setdiff(seq_len(p), k)
    for (s in other) data[[paste0("L", s)]] <- L[, s]
    for (s in seq_len(p)) data[[paste0("WL", s)]] <- network_lag[, s]
    fit <- if (covariate_types[k] == "continuous") {
      capture_model_fit(stats::lm(Lk ~ ., data = data))
    } else {
      capture_model_fit(stats::glm(Lk ~ ., data = data, family = stats::binomial()))
    }
    status <- list(ok = fit$ok, warnings = fit$warnings, error = fit$error, issues = character())
    if (!fit$ok) {
      status$issues <- sprintf("L[%d] fit failed: %s", k, fit$error)
      statuses[[k]] <- status
      issues <- c(issues, status$issues)
      next
    }
    other_terms <- if (length(other)) paste0("L", other) else character()
    terms <- c("(Intercept)", other_terms, paste0("WL", seq_len(p)))
    extracted <- extract_named_coefficients(fit$fit, terms, sprintf("L[%d]", k))
    status$issues <- extracted$issues
    status$ok <- !length(extracted$issues)
    intercept[k] <- extracted$values[["(Intercept)"]]
    for (s in other) within_node[k, s] <- extracted$values[[paste0("L", s)]]
    for (s in seq_len(p)) network_coef[k, s] <- extracted$values[[paste0("WL", s)]]
    if (covariate_types[k] == "continuous") {
      sigma[k] <- summary(fit$fit)$sigma
      if (!is.finite(sigma[k]) || sigma[k] <= 0) {
        status$ok <- FALSE
        status$issues <- c(status$issues, sprintf("L[%d] has an invalid residual SD.", k))
      }
    }
    statuses[[k]] <- status
    issues <- c(issues, status$issues)
  }

  list(
    model = list(intercept = intercept, within_node = within_node,
                 network_lag = network_coef, sigma = sigma),
    status = statuses,
    issues = unique(issues)
  )
}

.fit_treatment <- function(A, L, network) {
  p <- ncol(L)
  W_int_L <- network$W_int_A %*% L
  data <- data.frame(A = A)
  for (k in seq_len(p)) data[[paste0("L", k)]] <- L[, k]
  for (k in seq_len(p)) data[[paste0("WL", k)]] <- W_int_L[, k]
  data$A_dep <- as.numeric(network$W_dep_A %*% A)
  fit <- capture_model_fit(stats::glm(A ~ ., data = data, family = stats::binomial()))
  status <- list(ok = fit$ok, warnings = fit$warnings, error = fit$error, issues = character())
  if (!fit$ok) {
    status$issues <- sprintf("A fit failed: %s", fit$error)
    return(list(model = list(intercept = NA_real_, own_covariates = rep(NA_real_, p),
                             network_covariates = rep(NA_real_, p), dependence = NA_real_), status = status))
  }
  terms <- c("(Intercept)", paste0("L", seq_len(p)), paste0("WL", seq_len(p)), "A_dep")
  extracted <- extract_named_coefficients(fit$fit, terms, "A")
  status$issues <- extracted$issues
  status$ok <- !length(extracted$issues)
  list(
    model = list(
      intercept = extracted$values[["(Intercept)"]],
      own_covariates = unname(extracted$values[paste0("L", seq_len(p))]),
      network_covariates = unname(extracted$values[paste0("WL", seq_len(p))]),
      dependence = extracted$values[["A_dep"]]
    ),
    status = status
  )
}

.fit_outcome <- function(Y, A, L, network) {
  p <- ncol(L)
  units <- network$outcome_units
  Y_state <- numeric(network$N)
  Y_state[units] <- Y[units]
  data <- data.frame(Y = Y[units], A = A[units])
  for (k in seq_len(p)) data[[paste0("L", k)]] <- L[units, k]
  data$A_int <- as.numeric(network$W_int_Y %*% A)[units]
  W_int_L <- network$W_int_Y %*% L
  for (k in seq_len(p)) data[[paste0("WL", k)]] <- W_int_L[units, k]
  data$Y_dep <- as.numeric(network$W_dep_Y %*% Y_state)[units]
  fit <- capture_model_fit(stats::glm(Y ~ ., data = data, family = stats::binomial()))
  status <- list(ok = fit$ok, warnings = fit$warnings, error = fit$error, issues = character())
  if (!fit$ok) {
    status$issues <- sprintf("Y fit failed: %s", fit$error)
    return(list(model = list(intercept = NA_real_, treatment = NA_real_,
                             own_covariates = rep(NA_real_, p), network_treatment = NA_real_,
                             network_covariates = rep(NA_real_, p), dependence = NA_real_), status = status))
  }
  terms <- c("(Intercept)", "A", paste0("L", seq_len(p)), "A_int",
             paste0("WL", seq_len(p)), "Y_dep")
  extracted <- extract_named_coefficients(fit$fit, terms, "Y")
  status$issues <- extracted$issues
  status$ok <- !length(extracted$issues)
  list(
    model = list(
      intercept = extracted$values[["(Intercept)"]],
      treatment = extracted$values[["A"]],
      own_covariates = unname(extracted$values[paste0("L", seq_len(p))]),
      network_treatment = extracted$values[["A_int"]],
      network_covariates = unname(extracted$values[paste0("WL", seq_len(p))]),
      dependence = extracted$values[["Y_dep"]]
    ),
    status = status
  )
}

fit_pl <- function(Y, A, L, network, covariate_types = NULL) {
  network <- as_network(network)
  L <- as_covariate_matrix(L, network$N)
  covariate_types <- resolve_covariate_types(L, covariate_types)
  A <- as.numeric(A)
  Y <- as.numeric(Y)
  if (length(A) != network$N || any(!is.finite(A)) || !all(A %in% c(0, 1))) {
    stop("A must be a complete binary vector of length N.", call. = FALSE)
  }
  if (length(Y) != network$N || any(!is.finite(Y[network$outcome_units])) ||
      !all(Y[network$outcome_units] %in% c(0, 1))) {
    stop("Y must be binary and observed for every outcome_units index.", call. = FALSE)
  }

  covariate_fit <- .fit_covariates(L, covariate_types, network$W_dep_L)
  treatment_fit <- .fit_treatment(A, L, network)
  outcome_fit <- .fit_outcome(Y, A, L, network)
  issues <- unique(c(covariate_fit$issues, treatment_fit$status$issues, outcome_fit$status$issues))
  warnings <- unique(c(
    unlist(lapply(covariate_fit$status, `[[`, "warnings"), use.names = FALSE),
    treatment_fit$status$warnings,
    outcome_fit$status$warnings
  ))
  fit_status <- list(
    ok = !length(issues),
    issues = issues,
    warnings = warnings,
    model_status = list(L = covariate_fit$status, A = treatment_fit$status, Y = outcome_fit$status)
  )
  make_params(
    covariate_model = covariate_fit$model,
    treatment_model = treatment_fit$model,
    outcome_model = outcome_fit$model,
    meta = list(
      p = ncol(L),
      covariate_types = covariate_types,
      fit_status = fit_status,
      outcome_units = network$outcome_units
    )
  )
}

summarize_pl <- function(pl_list, true_params, label) {
  flat <- lapply(pl_list, params_flat_vector)
  names_all <- unique(unlist(lapply(flat, names), use.names = FALSE))
  estimates <- do.call(rbind, lapply(flat, function(x) {
    row <- setNames(rep(NA_real_, length(names_all)), names_all)
    row[names(x)] <- x
    row
  }))
  truth <- params_flat_vector(as_params(true_params))
  cat(sprintf("\n--- PL estimates: %s  (S = %d) ---\n", label, length(pl_list)))
  cat(sprintf("  %-16s  %9s  %9s  %9s  %9s\n", "param", "true", "mean", "bias", "SE"))
  for (name in names_all) {
    truth_value <- truth[[name]] %||% NA_real_
    values <- estimates[, name]
    estimate <- mean(values, na.rm = TRUE)
    cat(sprintf("  %-16s  %9.3f  %9.3f  %+9.3f  %9.3f\n",
                name, truth_value, estimate, estimate - truth_value,
                stats::sd(values, na.rm = TRUE)))
  }
  invisible(estimates)
}
