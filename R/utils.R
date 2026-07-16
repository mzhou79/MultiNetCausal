logistic <- function(x) {
  x <- pmax(pmin(x, 35), -35)
  1 / (1 + exp(-x))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

normalize_outcome_units <- function(outcome_units, N) {
  if (is.null(outcome_units)) return(seq_len(N))
  if (!is.numeric(outcome_units) || any(!is.finite(outcome_units)) ||
      any(outcome_units != as.integer(outcome_units))) {
    stop("outcome_units must contain integer indices.", call. = FALSE)
  }
  outcome_units <- as.integer(outcome_units)
  if (length(outcome_units) == 0L || any(outcome_units < 1L) || any(outcome_units > N)) {
    stop(sprintf("outcome_units must contain indices in 1:%d.", N), call. = FALSE)
  }
  if (anyDuplicated(outcome_units)) {
    stop("outcome_units must contain unique indices.", call. = FALSE)
  }
  outcome_units
}

assert_square_matrix <- function(x, name) {
  if (is.null(x)) stop(sprintf("%s is required.", name), call. = FALSE)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) != ncol(x)) {
    stop(sprintf("%s must be square; got %d x %d.", name, nrow(x), ncol(x)), call. = FALSE)
  }
  if (any(!is.finite(x))) stop(sprintf("%s contains non-finite values.", name), call. = FALSE)
  x
}

make_network <- function(W_dep_L, W_int_A, W_dep_A, W_int_Y, W_dep_Y,
                         outcome_units = NULL) {
  matrices <- list(
    W_dep_L = assert_square_matrix(W_dep_L, "W_dep_L"),
    W_int_A = assert_square_matrix(W_int_A, "W_int_A"),
    W_dep_A = assert_square_matrix(W_dep_A, "W_dep_A"),
    W_int_Y = assert_square_matrix(W_int_Y, "W_int_Y"),
    W_dep_Y = assert_square_matrix(W_dep_Y, "W_dep_Y")
  )
  N <- nrow(matrices$W_dep_L)
  if (any(vapply(matrices, nrow, integer(1L)) != N)) {
    stop("All five network matrices must have the same dimension.", call. = FALSE)
  }
  structure(
    c(matrices, list(N = as.integer(N), outcome_units = normalize_outcome_units(outcome_units, N))),
    class = "multinet_network"
  )
}

as_network <- function(network, outcome_units = NULL) {
  if (inherits(network, "multinet_network")) {
    if (!is.null(outcome_units)) {
      network$outcome_units <- normalize_outcome_units(outcome_units, network$N)
    }
    return(network)
  }
  if (!is.list(network)) stop("network must be a list.", call. = FALSE)
  if (!is.null(network$W_dep_L) && !is.null(network$W_int_A) &&
      !is.null(network$W_dep_A) && !is.null(network$W_int_Y) &&
      !is.null(network$W_dep_Y)) {
    return(make_network(
      network$W_dep_L, network$W_int_A, network$W_dep_A,
      network$W_int_Y, network$W_dep_Y,
      outcome_units = outcome_units %||% network$outcome_units
    ))
  }
  if (is.null(network$W_int) || is.null(network$W_dep)) {
    stop("network must contain either the five network roles or legacy W_int and W_dep.", call. = FALSE)
  }
  make_network(
    W_dep_L = network$W_dep,
    W_int_A = network$W_int,
    W_dep_A = network$W_dep,
    W_int_Y = network$W_int,
    W_dep_Y = network$W_dep,
    outcome_units = outcome_units
  )
}

as_covariate_matrix <- function(L, N = NULL) {
  L <- if (is.matrix(L)) L else matrix(as.numeric(L), ncol = 1L)
  L <- as.matrix(L)
  storage.mode(L) <- "double"
  if (!is.null(N) && nrow(L) != N) {
    stop(sprintf("L has %d rows but the network has %d units.", nrow(L), N), call. = FALSE)
  }
  if (ncol(L) < 1L || any(!is.finite(L))) {
    stop("L must contain at least one finite covariate column and no missing values.", call. = FALSE)
  }
  L
}

resolve_covariate_types <- function(L, covariate_types = NULL) {
  p <- ncol(L)
  if (is.null(covariate_types)) {
    covariate_types <- vapply(seq_len(p), function(k) {
      if (all(L[, k] %in% c(0, 1))) "binary" else "continuous"
    }, character(1L))
  }
  if (length(covariate_types) != p || any(!covariate_types %in% c("binary", "continuous"))) {
    stop("covariate_types must contain only 'binary' or 'continuous', one value per column of L.", call. = FALSE)
  }
  for (k in which(covariate_types == "binary")) {
    if (!all(L[, k] %in% c(0, 1))) {
      stop(sprintf("Binary covariate column %d must contain only 0 and 1.", k), call. = FALSE)
    }
  }
  as.character(covariate_types)
}

make_params <- function(covariate_model, treatment_model, outcome_model, meta = list()) {
  if (is.null(meta$p)) meta$p <- length(covariate_model$intercept)
  structure(list(
    covariate_model = covariate_model,
    treatment_model = treatment_model,
    outcome_model = outcome_model,
    meta = meta
  ), class = c("multinet_params", "list"))
}

as_params <- function(params, covariate_types = NULL) {
  if (inherits(params, "multinet_params")) {
    if (!is.null(covariate_types)) params$meta$covariate_types <- covariate_types
    return(params)
  }
  if (!is.list(params)) stop("params must be a list or multinet_params object.", call. = FALSE)
  get_param <- function(name, default = NULL) {
    value <- params[[name, exact = TRUE]]
    if (is.null(value)) default else value
  }
  p <- max(
    1L,
    length(get_param("beta2", numeric())),
    length(get_param("beta4", numeric())),
    length(get_param("lambda0", numeric())),
    length(get_param("omega", numeric())),
    length(get_param("rho", numeric())),
    length(get_param("gamma1", numeric()))
  )
  covariate_types <- covariate_types %||% rep("binary", p)
  if (length(covariate_types) != p || any(!covariate_types %in% c("binary", "continuous"))) {
    stop("covariate_types must contain only 'binary' or 'continuous'.", call. = FALSE)
  }
  omega <- rep_len(as.numeric(get_param("omega", 0)), p)
  sigma <- rep_len(as.numeric(get_param("sigma", get_param("sigma_sq", 1))), p)
  rho_L <- as.matrix(get_param("rho_L", matrix(0, p, p)))
  nu_L <- as.matrix(get_param("nu_L", diag(omega, p, p)))
  if (!all(dim(rho_L) == c(p, p)) || !all(dim(nu_L) == c(p, p))) {
    stop("rho_L and nu_L must both be p x p matrices.", call. = FALSE)
  }
  make_params(
    covariate_model = list(
      intercept = rep_len(as.numeric(get_param("lambda0", 0)), p),
      within_node = rho_L,
      network_lag = nu_L,
      sigma = sigma
    ),
    treatment_model = list(
      intercept = as.numeric(get_param("gamma0", 0)),
      own_covariates = rep_len(as.numeric(get_param("rho", 0)), p),
      network_covariates = rep_len(as.numeric(get_param("gamma1", 0)), p),
      dependence = as.numeric(get_param("eta", 0))
    ),
    outcome_model = list(
      intercept = as.numeric(get_param("beta0", 0)),
      treatment = as.numeric(get_param("beta1", 0)),
      own_covariates = rep_len(as.numeric(get_param("beta2", 0)), p),
      network_treatment = as.numeric(get_param("beta3", 0)),
      network_covariates = rep_len(as.numeric(get_param("beta4", 0)), p),
      dependence = as.numeric(get_param("theta", 0))
    ),
    meta = list(
      p = as.integer(p),
      covariate_types = as.character(covariate_types),
      fit_status = get_param("fit_status", list(ok = TRUE, issues = character(), model_status = list())),
      legacy_input = TRUE
    )
  )
}

params_to_cpp <- function(params) {
  params <- as_params(params)
  list(
    lambda0 = as.numeric(params$covariate_model$intercept),
    rho_L = as.matrix(params$covariate_model$within_node),
    nu_L = as.matrix(params$covariate_model$network_lag),
    sigma = as.numeric(params$covariate_model$sigma),
    gamma0 = as.numeric(params$treatment_model$intercept),
    rho = as.numeric(params$treatment_model$own_covariates),
    gamma1 = as.numeric(params$treatment_model$network_covariates),
    eta = as.numeric(params$treatment_model$dependence),
    beta0 = as.numeric(params$outcome_model$intercept),
    beta1 = as.numeric(params$outcome_model$treatment),
    beta2 = as.numeric(params$outcome_model$own_covariates),
    beta3 = as.numeric(params$outcome_model$network_treatment),
    beta4 = as.numeric(params$outcome_model$network_covariates),
    theta = as.numeric(params$outcome_model$dependence)
  )
}

params_are_valid <- function(params) {
  params <- as_params(params)
  fit_status <- params$meta$fit_status %||% list(ok = TRUE)
  blocks <- params_to_cpp(params)
  isTRUE(fit_status$ok) && all(vapply(blocks, function(x) all(is.finite(as.numeric(x))), logical(1L)))
}

params_problem_text <- function(params) {
  params <- as_params(params)
  fit_status <- params$meta$fit_status %||% list(issues = character())
  issues <- unique(unlist(fit_status$issues, use.names = FALSE))
  if (!length(issues)) {
    blocks <- params_to_cpp(params)
    bad <- names(blocks)[!vapply(blocks, function(x) all(is.finite(as.numeric(x))), logical(1L))]
    if (length(bad)) issues <- paste("non-finite parameter block:", bad)
  }
  paste(issues, collapse = "; ")
}

params_flat_vector <- function(params) {
  params <- as_params(params)
  named_vector <- function(prefix, x) {
    x <- as.numeric(x)
    if (length(x) == 1L) setNames(x, prefix) else setNames(x, sprintf("%s[%d]", prefix, seq_along(x)))
  }
  out <- c(
    beta0 = params$outcome_model$intercept,
    beta1 = params$outcome_model$treatment,
    named_vector("beta2", params$outcome_model$own_covariates),
    beta3 = params$outcome_model$network_treatment,
    named_vector("beta4", params$outcome_model$network_covariates),
    theta = params$outcome_model$dependence,
    named_vector("lambda0", params$covariate_model$intercept),
    named_vector("omega", diag(as.matrix(params$covariate_model$network_lag))),
    gamma0 = params$treatment_model$intercept,
    named_vector("rho", params$treatment_model$own_covariates),
    named_vector("gamma1", params$treatment_model$network_covariates),
    eta = params$treatment_model$dependence
  )
  p <- params$meta$p
  if (p > 1L) {
    rho_L <- as.matrix(params$covariate_model$within_node)
    nu_L <- as.matrix(params$covariate_model$network_lag)
    for (k in seq_len(p)) for (s in seq_len(p)) {
      if (k != s) out[sprintf("rho_L[%d,%d]", k, s)] <- rho_L[k, s]
      out[sprintf("nu_L[%d,%d]", k, s)] <- nu_L[k, s]
    }
  }
  out
}

load_gibbs_cpp <- local({
  loaded <- FALSE
  function(root = getwd()) {
    needed <- c("gibbs_psi_cpp", "generate_replications_cpp")
    if (loaded && all(vapply(needed, exists, logical(1L), mode = "function", inherits = TRUE))) {
      return(invisible(TRUE))
    }
    if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.", call. = FALSE)
    cpp_file <- file.path(root, "src", "gibbs.cpp")
    if (!file.exists(cpp_file)) stop("Missing C++ Gibbs source: ", cpp_file, call. = FALSE)
    Rcpp::sourceCpp(cpp_file, rebuild = FALSE, verbose = FALSE)
    missing <- needed[!vapply(needed, exists, logical(1L), mode = "function", inherits = TRUE)]
    if (length(missing)) stop("C++ functions failed to load: ", paste(missing, collapse = ", "), call. = FALSE)
    loaded <<- TRUE
    invisible(TRUE)
  }
})

capture_model_fit <- function(expr) {
  warnings <- character()
  fit <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, fit = NULL, error = conditionMessage(fit), warnings = warnings))
  }
  list(ok = TRUE, fit = fit, error = NULL, warnings = warnings)
}

extract_named_coefficients <- function(fit, terms, model_name) {
  coefs <- stats::coef(fit)
  values <- setNames(rep(NA_real_, length(terms)), terms)
  problems <- character()
  for (term in terms) {
    if (!(term %in% names(coefs)) || length(coefs[[term]]) != 1L || !is.finite(coefs[[term]])) {
      problems <- c(problems, term)
    } else {
      values[[term]] <- unname(coefs[[term]])
    }
  }
  issues <- if (length(problems)) {
    sprintf("%s missing, aliased, or non-finite coefficients: %s", model_name, paste(problems, collapse = ", "))
  } else character()
  list(values = values, issues = issues)
}
