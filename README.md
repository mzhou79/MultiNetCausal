# MultiNetCausalRcpp

Submitted code for pseudo-likelihood estimation, Auto-G-style Gibbs effect estimation, and a parametric bootstrap.

## Structure

```text
R/
  utils.R       network, parameter, validation, and C++ loading helpers
  estimation.R  fit_pl()
  causal.R      gibbs_psi() and estimate_effects()
  bootstrap.R   sample_data() and bootstrap_effects()
src/
  gibbs.cpp     effect Gibbs sampler and original simulation DGP
examples/
  synthetic_full_cc.R
simulation/
  sim_utils.R
  run_study1_block.R
  run_study2_random_network.R
results/
  final_study1_N20/
  final_strong_dep_N200/
  final_strong_inter_weak_dep_N200/
```

The three final simulation result folders are included in this handoff. Each folder contains its run settings, saved networks/replications, and CSV summaries; the strong-interference/weak-dependence folder also contains the bootstrap validation outputs. Full workspace dumps are intentionally omitted because they retain stale function definitions and are not needed to reproduce or inspect the reported results.

There is no permanent `tests/` folder. The simulation scripts and the synthetic example are the executable checks of the submitted code.

## Five network roles

Use one network object:

```r
network <- make_network(
  W_dep_L,
  W_int_A,
  W_dep_A,
  W_int_Y,
  W_dep_Y,
  outcome_units = NULL
)
```

The five matrices must share the same N x N index space:

- `W_dep_L`: dependence network for the covariate model.
- `W_int_A`: covariate-interference network for the observed-treatment model.
- `W_dep_A`: dependence network for observed treatment.
- `W_int_Y`: treatment/covariate-interference network for the outcome model.
- `W_dep_Y`: outcome-dependence network.

`outcome_units = NULL` means all `1:N` units. When it is supplied:

- L and observed A are fit and updated on all N units;
- Y is fit and updated only on `outcome_units`;
- DE, IE, and ATE are calculated and averaged only over `outcome_units`.

Included isolates remain valid outcome units even when their network rows contain only zeros.

## Covariates

Only two covariate types are supported:

```r
covariate_types <- c("binary", "continuous")
```

Binary columns must contain only 0 and 1. Continuous columns must be finite numeric values. If `covariate_types` is omitted, a 0/1 column is inferred as binary and every other numeric column is inferred as continuous.

## Analysis

```r
source("R/utils.R")
source("R/estimation.R")
source("R/causal.R")
source("R/bootstrap.R")

pl_fit <- fit_pl(Y, A, L, network, covariate_types)

effects <- estimate_effects(
  L = L,
  params = pl_fit,
  network = network,
  alpha = 0.5,
  K = 50,
  burn = 20,
  covariate_types = covariate_types
)
```

Effect results use `DE`, `IE`, and `ATE`.

## Bootstrap report

```r
bootstrap <- bootstrap_effects(
  Y = Y,
  A = A,
  L = L,
  network = network,
  covariate_types = covariate_types,
  B = 200,
  bootstrap_burn = 500,
  alpha = 0.5,
  K = 50,
  burn = 20,
  seed = 20260713,
  n_cores = 1
)

print(bootstrap$summary, row.names = FALSE)
```

For each estimand, `Bootstrap SE` is the SD of the successful bootstrap point estimates. The reported interval is the estimate-centered 95% Wald interval:

```text
Estimate +/- 1.96 * Bootstrap SE
```

It is not a percentile interval. Failed replications remain in `bootstrap$draws` with an error message and are reported through `bootstrap$failures`.

## Synthetic full/CC example

Run from the folder root:

```r
source("examples/synthetic_full_cc.R")
```

The example uses the Study 2 random-network construction with N=200: an interference graph, an independently generated auxiliary graph, and their union as the augmented/dependence graph. It generates one dataset with one continuous and four binary covariates. Every covariate has augmented-network dependence, and the within-unit `rho_L` matrix includes nonzero cross-covariate terms as in the uConnect multicovariate model. L and A are observed for the full cohort; Y is observed on a randomly selected 150-unit induced subnetwork whose interference and augmented matrices are embedded in the full N x N index space.

## Manuscript simulation compatibility

The manuscript simulation remains the p = 1, all-unit special case:

```r
network <- make_network(
  W_dep_L = W_dep,
  W_int_A = W_int,
  W_dep_A = W_dep,
  W_int_Y = W_int,
  W_dep_Y = W_dep,
  outcome_units = NULL
)
```

Legacy simulation network lists containing only `W_int` and `W_dep` are mapped to these five roles by `as_network()`. The DGP, raw-adjacency convention, seeds, Gibbs allocation logic, and DE/IE/ATE definitions are unchanged.

Run the manuscript scripts from the folder root:

```r
source("simulation/run_study1_block.R")
source("simulation/run_study2_random_network.R")
```
