# Auto-g-computation with Multiple Network-induced Dependencies

This repository contains the code that can reproduce the simulation result of the manuscript.

## `R/`

- `utils.R` defines the five-network-role data structure and shared helpers.
- `estimation.R` implements the pseudo-likelihood fit of the covariate, treatment, and outcome full conditionals.
- `causal.R` implements the Gibbs-sampler auto-g-computation algorithm (Algorithm 1) that produces the DE/IE/ATE estimates.
- `bootstrap.R` implements the fixed-network parametric bootstrap used for uncertainty quantification.

## `simulation/`

- `run_study1_block.R` reproduces Study 1 (the N = 20 four-unit block network).
- `run_study2_random_network.R` reproduces Study 2 (the N = 200 random-network design).
- `sim_utils.R` holds the shared data-generating and evaluation utilities used by both.

## `results/`

- `final_study1_N20/` — output from the Study 1 (N = 20) validation setting.
- `final_strong_dep_N200/` — output from the Study 2 (N = 200) setting under
  the original interference/dependence parameters.
- `final_strong_inter_weak_dep_N200/` — output from the Study 2 stress
  setting with stronger interference and weaker within-variable dependence,
  including a `bootstrap_validation/` subfolder checking the bootstrap
  procedure against Monte Carlo variability.

## `examples/`

- `synthetic_full_cc.R` is a self-contained synthetic example illustrating the
full-cohort / outcome-complete-cohort setup used in the uConnect data
application, useful for confirming the pipeline runs end-to-end on
multivariate covariates without needing the real data.
