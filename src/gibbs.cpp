#include <Rcpp.h>
using namespace Rcpp;

static inline double clamp_logit(double x) {
  if (x > 35.0) return 35.0;
  if (x < -35.0) return -35.0;
  return x;
}

static inline double logistic_cpp(double x) {
  x = clamp_logit(x);
  return 1.0 / (1.0 + std::exp(-x));
}

static inline double param_double(const List& params, const char* name) {
  if (!params.containsElementNamed(name)) return 0.0;
  SEXP val_se = params[name];
  if (Rf_length(val_se) == 0) return 0.0;
  if (Rf_isReal(val_se)) {
    NumericVector val(val_se);
    if (NumericVector::is_na(val[0])) return 0.0;
    return val[0];
  } else if (Rf_isInteger(val_se)) {
    IntegerVector val(val_se);
    if (IntegerVector::is_na(val[0])) return 0.0;
    return static_cast<double>(val[0]);
  }
  return 0.0;
}

static inline double param_vec_elem(const List& params, const char* name, int k) {
  if (!params.containsElementNamed(name)) return 0.0;
  SEXP val_se = params[name];
  if (Rf_length(val_se) <= k) return 0.0;
  NumericVector val(val_se);
  if (NumericVector::is_na(val[k])) return 0.0;
  return val[k];
}

static inline double param_mat_elem(const List& params, const char* name,
                                    int r, int c, int nrows) {
  if (!params.containsElementNamed(name)) return 0.0;
  SEXP val_se = params[name];
  if (Rf_isMatrix(val_se)) {
    NumericMatrix M(val_se);
    if (r >= M.nrow() || c >= M.ncol()) return 0.0;
    if (NumericMatrix::is_na(M(r, c))) return 0.0;
    return M(r, c);
  }
  NumericVector val(val_se);
  int idx = r + c * nrows;
  if (idx >= Rf_length(val_se)) return 0.0;
  if (NumericVector::is_na(val[idx])) return 0.0;
  return val[idx];
}

static void set_seed_cpp(int seed) {
  if (seed == NA_INTEGER) return;
  Environment base = Environment::namespace_env("base");
  Function set_seed = base["set.seed"];
  set_seed(seed);
}

static NumericVector mat_vec_mult(const NumericMatrix& W, const IntegerVector& x) {
  int N = W.nrow();
  NumericVector out(N);
  for (int i = 0; i < N; ++i) {
    double s = 0.0;
    for (int j = 0; j < N; ++j) s += W(i, j) * x[j];
    out[i] = s;
  }
  return out;
}

static void add_col_delta(NumericVector& nb, const NumericMatrix& W, int col, double delta) {
  if (delta == 0.0) return;
  int N = W.nrow();
  for (int r = 0; r < N; ++r) {
    double w = W(r, col);
    if (w != 0.0) nb[r] += w * delta;
  }
}

static void add_col_delta_m(NumericMatrix& nb, const NumericMatrix& W,
                            int col, int k, double delta) {
  if (delta == 0.0) return;
  int N = W.nrow();
  for (int r = 0; r < N; ++r) {
    double w = W(r, col);
    if (w != 0.0) nb(r, k) += w * delta;
  }
}

// =====================================================================
// Scoped multi-covariate Bernoulli Gibbs.
// Separates L-dependence scope from Y interference/dependence scope.
// Scalar is handled as p = 1 by passing an N x 1 matrix.
// [[Rcpp::export]]
NumericVector gibbs_psi_cpp(NumericMatrix L_start,
                                            NumericMatrix W_dep_L,
                                            NumericMatrix W_int_Y,
                                            NumericMatrix W_dep_Y,
                                            List params,
                                            double alpha,
                                            int n_iter,
                                            int burn_in,
                                            int seed = NA_INTEGER,
                                            int unit_index = -1,
                                            int treatment_value = 0,
                                            IntegerVector cov_types = IntegerVector(),
                                            IntegerVector outcome_units = IntegerVector()) {
  RNGScope scope;
  set_seed_cpp(seed == NA_INTEGER ? NA_INTEGER : seed);

  int N = L_start.nrow();
  int p = L_start.ncol();
  if (W_dep_L.nrow() != N || W_dep_L.ncol() != N ||
      W_int_Y.nrow() != N || W_int_Y.ncol() != N ||
      W_dep_Y.nrow() != N || W_dep_Y.ncol() != N) {
    stop("Scoped Gibbs matrices must all be NxN.");
  }
  if (n_iter <= burn_in) stop("n_iter must be larger than burn_in.");

  if (cov_types.size() == 0) cov_types = IntegerVector(p, 0);
  if (cov_types.size() != p) stop("cov_types must have length p.");

  LogicalVector is_outcome(N, false);
  if (outcome_units.size() == 0) {
    for (int i = 0; i < N; ++i) is_outcome[i] = true;
  } else {
    for (int r = 0; r < outcome_units.size(); ++r) {
      int i = outcome_units[r] - 1;
      if (i < 0 || i >= N) stop("outcome_units contains an invalid index.");
      is_outcome[i] = true;
    }
  }

  int ui = (unit_index >= 1) ? unit_index - 1 : -1;
  if (ui >= N) stop("unit_index out of range.");

  double beta0 = param_double(params, "beta0");
  double beta1 = param_double(params, "beta1");
  double beta3 = param_double(params, "beta3");
  double theta = param_double(params, "theta");

  std::vector<double> beta2_v(p), beta4_v(p), lambda0_v(p), sigma_v(p);
  std::vector<double> rho_L(p * p, 0.0), nu_L(p * p, 0.0);
  for (int k = 0; k < p; ++k) {
    beta2_v[k] = param_vec_elem(params, "beta2", k);
    beta4_v[k] = param_vec_elem(params, "beta4", k);
    lambda0_v[k] = param_vec_elem(params, "lambda0", k);
    sigma_v[k] = param_vec_elem(params, "sigma", k);
    if (sigma_v[k] <= 0) sigma_v[k] = param_vec_elem(params, "sigma_sq", k);
    if (sigma_v[k] <= 0) sigma_v[k] = 1.0;
    for (int s = 0; s < p; ++s) {
      rho_L[k * p + s] = param_mat_elem(params, "rho_L", k, s, p);
      nu_L[k * p + s] = param_mat_elem(params, "nu_L", k, s, p);
    }
  }

  NumericMatrix L_cur = clone(L_start);
  IntegerVector Y(N, 0);
  IntegerVector A(N);
  for (int i = 0; i < N; ++i) A[i] = R::rbinom(1.0, alpha);
  if (ui >= 0) A[ui] = treatment_value;

  NumericVector int_nb_A(N);
  for (int i = 0; i < N; ++i) {
    double s = 0.0;
    for (int j = 0; j < N; ++j) s += W_int_Y(i, j) * A[j];
    int_nb_A[i] = s;
  }

  NumericMatrix dep_nb_L(N, p), int_nb_L(N, p);
  for (int k = 0; k < p; ++k) {
    for (int i = 0; i < N; ++i) {
      double s_dep = 0.0, s_int = 0.0;
      for (int j = 0; j < N; ++j) {
        double lv = L_cur(j, k);
        s_dep += W_dep_L(i, j) * lv;
        s_int += W_int_Y(i, j) * lv;
      }
      dep_nb_L(i, k) = s_dep;
      int_nb_L(i, k) = s_int;
    }
  }

  for (int i = 0; i < N; ++i) {
    double lin_Y = beta0 + beta1 * A[i] + beta3 * int_nb_A[i];
    for (int k = 0; k < p; ++k) {
      lin_Y += beta2_v[k] * L_cur(i, k) + beta4_v[k] * int_nb_L(i, k);
    }
    Y[i] = is_outcome[i] ? R::rbinom(1.0, logistic_cpp(lin_Y)) : 0;
  }

  NumericVector dep_nb_Y(N);
  for (int i = 0; i < N; ++i) {
    double s = 0.0;
    for (int j = 0; j < N; ++j) s += W_dep_Y(i, j) * Y[j];
    dep_nb_Y[i] = s;
  }

  NumericVector acc(N, 0.0);
  int kept = 0;
  for (int m = 1; m <= n_iter; ++m) {
    for (int i = 0; i < N; ++i) {
      int new_A = (ui >= 0 && i == ui) ? treatment_value : R::rbinom(1.0, alpha);
      int dA = new_A - A[i];
      if (dA != 0) {
        A[i] = new_A;
        add_col_delta(int_nb_A, W_int_Y, i, static_cast<double>(dA));
      }
    }

    for (int k = 0; k < p; ++k) {
      bool is_cont = (cov_types.size() > k && cov_types[k] == 1);
      for (int i = 0; i < N; ++i) {
        double lin_k = lambda0_v[k];
        for (int s = 0; s < p; ++s) {
          lin_k += rho_L[k * p + s] * L_cur(i, s);
          lin_k += nu_L[k * p + s] * dep_nb_L(i, s);
        }
        double old_Lk = L_cur(i, k);
        double new_Lk = is_cont ? R::rnorm(lin_k, sigma_v[k])
                                : R::rbinom(1.0, logistic_cpp(lin_k));
        L_cur(i, k) = new_Lk;
        double dL = new_Lk - old_Lk;
        if (dL != 0.0) {
          add_col_delta_m(dep_nb_L, W_dep_L, i, k, dL);
          add_col_delta_m(int_nb_L, W_int_Y, i, k, dL);
        }
      }
    }

    for (int i = 0; i < N; ++i) {
      if (!is_outcome[i]) continue;
      double lin_Y = beta0 + beta1 * A[i] + beta3 * int_nb_A[i] + theta * dep_nb_Y[i];
      for (int k = 0; k < p; ++k) {
        lin_Y += beta2_v[k] * L_cur(i, k) + beta4_v[k] * int_nb_L(i, k);
      }
      int old_Y = Y[i];
      Y[i] = R::rbinom(1.0, logistic_cpp(lin_Y));
      int dY = Y[i] - old_Y;
      if (dY != 0) add_col_delta(dep_nb_Y, W_dep_Y, i, static_cast<double>(dY));
    }

    if (m > burn_in) {
      for (int i = 0; i < N; ++i) acc[i] += Y[i];
      kept++;
    }
  }

  for (int i = 0; i < N; ++i) acc[i] /= kept;
  return acc;
}

// =====================================================================
// DGP generator: independent full-network (L,A,Y) chains.
// Each replication starts from fresh Bernoulli(init_prob) L/A/Y, runs
// burn + 1 sweeps, and returns the retained state.
// [[Rcpp::export]]
List generate_replications_cpp(NumericMatrix W_int,
                               NumericMatrix W_dep,
                               List params,
                               int S,
                               int seed,
                               int burn,
                               double init_prob = 0.4) {
  RNGScope scope;
  int N = W_int.nrow();
  if (S < 1) stop("S must be at least 1.");
  if (burn < 0) stop("burn must be non-negative.");
  if (init_prob < 0.0 || init_prob > 1.0) stop("init_prob must be in [0, 1].");
  if (W_int.ncol() != N || W_dep.nrow() != N || W_dep.ncol() != N) {
    stop("Incompatible dimensions in generate_replications_cpp.");
  }

  double beta0 = param_double(params, "beta0");
  double beta1 = param_double(params, "beta1");
  double beta2 = param_double(params, "beta2");
  double beta3 = param_double(params, "beta3");
  double beta4 = param_double(params, "beta4");
  double theta = param_double(params, "theta");
  double lambda0 = param_double(params, "lambda0");
  double omega = param_double(params, "omega");
  double gamma0 = param_double(params, "gamma0");
  double rho = param_double(params, "rho");
  double gamma1 = param_double(params, "gamma1");
  double eta = param_double(params, "eta");

  List reps(S);
  for (int s = 0; s < S; ++s) {
    set_seed_cpp(seed + s);
    IntegerVector L(N), A(N), Y(N);
    for (int i = 0; i < N; ++i) L[i] = R::rbinom(1.0, init_prob);
    for (int i = 0; i < N; ++i) A[i] = R::rbinom(1.0, init_prob);
    for (int i = 0; i < N; ++i) Y[i] = R::rbinom(1.0, init_prob);

    NumericVector dep_nb_L = mat_vec_mult(W_dep, L);
    NumericVector int_nb_L = mat_vec_mult(W_int, L);
    NumericVector dep_nb_A = mat_vec_mult(W_dep, A);
    NumericVector int_nb_A = mat_vec_mult(W_int, A);
    NumericVector dep_nb_Y = mat_vec_mult(W_dep, Y);

    for (int m = 0; m < burn + 1; ++m) {
      for (int i = 0; i < N; ++i) {
        int old_L = L[i];
        L[i] = R::rbinom(1.0, logistic_cpp(lambda0 + omega * dep_nb_L[i]));
        int dL = L[i] - old_L;
        if (dL != 0) {
          add_col_delta(dep_nb_L, W_dep, i, dL);
          add_col_delta(int_nb_L, W_int, i, dL);
        }

        int old_A = A[i];
        double lin_A = gamma0 + rho * L[i] + gamma1 * int_nb_L[i] + eta * dep_nb_A[i];
        A[i] = R::rbinom(1.0, logistic_cpp(lin_A));
        int dA = A[i] - old_A;
        if (dA != 0) {
          add_col_delta(dep_nb_A, W_dep, i, dA);
          add_col_delta(int_nb_A, W_int, i, dA);
        }

        int old_Y = Y[i];
        double lin_Y = beta0 + beta1 * A[i] + beta2 * L[i] +
          beta3 * int_nb_A[i] + beta4 * int_nb_L[i] + theta * dep_nb_Y[i];
        Y[i] = R::rbinom(1.0, logistic_cpp(lin_Y));
        int dY = Y[i] - old_Y;
        if (dY != 0) add_col_delta(dep_nb_Y, W_dep, i, dY);
      }
    }

    reps[s] = List::create(_["L"] = L, _["A"] = A, _["Y"] = Y);
  }
  return reps;
}