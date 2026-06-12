###############################################################################
# uConnect Real Data Application — Multiple Covariates (Binary + Continuous)
# Chain Graph (c): Interference = sexual network, Dependence = augmented network
###############################################################################

suppressPackageStartupMessages({
  library(parallel)
  library(igraph)
  library(dplyr)
})

set.seed(2026)

## ============================================================================
## 0. CONFIGURATION
## ============================================================================

# Column names in full_final
NODE_ID_COL    <- "su_id"
EGO_ID_COL     <- "ego_uniqueid1"   # bridge to edge files
TREATMENT_COL  <- "pepknow"
OUTCOME_COL    <- "sero"

# Edge file column names (same structure for all 6 edge files)
EDGE_FROM_COL  <- "uniqueid1_from"
EDGE_TO_COL    <- "uniqueid1_to"

# Covariates: names in full_final, and their types
COV_NAMES  <- c("agecalc", "numpars1", "healthcov", "jailever", "empstat")
COV_TYPES  <- c("continuous", "continuous", "binary", "binary", "binary")
# drugsmarij excluded: not in paper Table 2 covariate list
p <- length(COV_NAMES)

# Designated model parameters
# beta2 and beta4 are length-p vectors (one per covariate)
beta0  <- -1.00
beta1  <-  0.50                              # own treatment
beta2  <- c(0.10, 0.05, 0.20, 0.15, 0.15)  # own covariates (agecalc, numpars1, healthcov, jailever, empstat)
beta3  <-  0.85                              # neighbor treatment (interference)
beta4  <- c(0.10, 0.05, 0.15, 0.10, 0.10)  # neighbor covariates (interference)
theta  <-  0.42                              # Y-Y dependence (augmented)

# Covariate model parameters
eta       <- c(-0.50, -0.50, -0.50, -0.50, -0.50)  # length-p intercepts
sigma_cov <- c(1.0, 1.0, NA, NA, NA)                # SD for continuous (NA for binary)
# set to 1 because we scale below

# omega: p x p cross-unit covariate association matrix (augmented network)
# Diagonal = same covariate, off-diagonal = 0 (simplest designation)
omega_mat <- diag(c(0.20, 0.20, 0.42, 0.42, 0.42))

# rho: p x p within-unit covariate association (0 = independent within unit)
rho_mat <- matrix(0, p, p)

# Gibbs settings
alpha   <- 0.50   # Bernoulli treatment allocation probability
n_iter  <- 5000L
burn_in <- 1000L
n_runs  <- 40L    # replicates

logistic <- function(x) 1 / (1 + exp(-x))

## ============================================================================
## 1. LOAD AND PREPARE NODE DATA
## ============================================================================

cat("=== Preparing node data ===\n")

# ---- 1a. Full cohort N=618 (for covariate neighbor lookup, Section 5) ----
node618 <- full_final   # all 618 rows after MI, no missing
node618 <- node618[!duplicated(node618[[NODE_ID_COL]]), ]
cat(sprintf("  N=618 cohort: %d rows\n", nrow(node618)))

# Scale continuous covariates using N=618 means/SDs
# Store scale parameters so N=335 gets same scaling
scale_center <- sapply(c("agecalc","numpars1"), function(v) mean(node618[[v]], na.rm=TRUE))
scale_sd     <- sapply(c("agecalc","numpars1"), function(v) sd(node618[[v]],   na.rm=TRUE))

node618$agecalc_sc  <- (node618$agecalc  - scale_center["agecalc"])  / scale_sd["agecalc"]
node618$numpars1_sc <- (node618$numpars1 - scale_center["numpars1"]) / scale_sd["numpars1"]

# Scaled covariate names used in model matrix
COV_NAMES_SC <- c("agecalc_sc", "numpars1_sc", "healthcov", "jailever", "empstat")

# Build N=618 covariate matrix (used only for cross-unit neighbor terms in covariate update)
L_mat618 <- as.matrix(node618[, COV_NAMES_SC])
rownames(L_mat618) <- as.character(node618[[NODE_ID_COL]])

# ego_uniqueid1 -> su_id crosswalk for N=618
ego2suid_618 <- setNames(as.character(node618[[NODE_ID_COL]]),
                         as.character(node618[[EGO_ID_COL]]))

# ---- 1b. Selected cohort N=335 (at-risk, non-missing outcome) ----
# Identify N=335: HIV-negative at baseline AND outcome observed
# Adjust this filter to match your actual selection criteria
node335 <- full_final %>%
  filter(!is.na(.data[[OUTCOME_COL]]),
         !is.na(.data[[TREATMENT_COL]])) %>%
  # If you have a flag for at-risk, add: filter(atrisk == 1)
  distinct(.data[[NODE_ID_COL]], .keep_all = TRUE)

cat(sprintf("  N=335 cohort: %d rows\n", nrow(node335)))

node335$agecalc_sc  <- (node335$agecalc  - scale_center["agecalc"])  / scale_sd["agecalc"]
node335$numpars1_sc <- (node335$numpars1 - scale_center["numpars1"]) / scale_sd["numpars1"]

node335$A <- as.integer(node335[[TREATMENT_COL]])
node335$Y <- as.integer(node335[[OUTCOME_COL]])

# Drop any remaining NA in A or covariates (post-MI should be none)
keep335 <- complete.cases(node335[, c("A", COV_NAMES_SC)])
if (sum(!keep335) > 0)
  cat(sprintf("  Dropping %d rows with residual NA.\n", sum(!keep335)))
node335 <- node335[keep335, ]
N_total <- nrow(node335)
cat(sprintf("  Final N for Gibbs: %d\n", N_total))

# Integer index 1..N_total for Gibbs (indexed by su_id)
node335$idx <- seq_len(N_total)
suid2idx_335 <- setNames(node335$idx, as.character(node335[[NODE_ID_COL]]))
ego2idx_335  <- setNames(node335$idx, as.character(node335[[EGO_ID_COL]]))

# ego_uniqueid1 -> su_id crosswalk for N=335
ego2suid_335 <- setNames(as.character(node335[[NODE_ID_COL]]),
                         as.character(node335[[EGO_ID_COL]]))

# Observed covariate matrix for N=335 (starting point for Gibbs)
L_mat335_obs <- as.matrix(node335[, COV_NAMES_SC])  # N_total x p

## ============================================================================
## 2. LOAD AND CLEAN EDGE FILES
##    Key: edges use ego_uniqueid1 IDs, not su_id
##         -> translate via ego2idx_335 and ego2suid_618
## ============================================================================

# Generic cleaner: returns data.frame with cols i, j (integer indices)
# For N=335 edges: map ego_uniqueid1 -> idx in node335
# For N=618 edges: map ego_uniqueid1 -> su_id (used in L_mat618 lookup)

clean_edges_335 <- function(edge_df, label) {
  e <- data.frame(
    from = as.character(edge_df[[EDGE_FROM_COL]]),
    to   = as.character(edge_df[[EDGE_TO_COL]]),
    stringsAsFactors = FALSE
  )
  n0 <- nrow(e)
  e  <- e[e$from != e$to, ]                          # drop self-loops
  swap <- e$from > e$to                               # canonicalize undirected
  tmp  <- e$from[swap]; e$from[swap] <- e$to[swap]; e$to[swap] <- tmp
  e    <- unique(e)
  # Keep only edges where BOTH ego IDs map to a node335 unit
  in_335 <- e$from %in% names(ego2idx_335) & e$to %in% names(ego2idx_335)
  cat(sprintf("  [%s N=335] raw=%d, kept=%d (%d orphan dropped)\n",
              label, n0, sum(in_335), sum(!in_335)))
  e <- e[in_335, ]
  e$i <- unname(ego2idx_335[e$from])
  e$j <- unname(ego2idx_335[e$to])
  e[, c("i","j")]
}

clean_edges_618 <- function(edge_df, label) {
  e <- data.frame(
    from = as.character(edge_df[[EDGE_FROM_COL]]),
    to   = as.character(edge_df[[EDGE_TO_COL]]),
    stringsAsFactors = FALSE
  )
  n0 <- nrow(e)
  e  <- e[e$from != e$to, ]
  swap <- e$from > e$to
  tmp  <- e$from[swap]; e$from[swap] <- e$to[swap]; e$to[swap] <- tmp
  e    <- unique(e)
  # Keep edges where BOTH endpoints are in N=618
  in_618 <- e$from %in% names(ego2suid_618) & e$to %in% names(ego2suid_618)
  cat(sprintf("  [%s N=618] raw=%d, kept=%d (%d orphan dropped)\n",
              label, n0, sum(in_618), sum(!in_618)))
  e <- e[in_618, ]
  # Store as su_id strings (for L_mat618 rowname lookup)
  e$su_from <- unname(ego2suid_618[e$from])
  e$su_to   <- unname(ego2suid_618[e$to])
  e[, c("su_from","su_to")]
}

cat("\n=== Cleaning edge files ===\n")

# N=335 edges (used directly in Gibbs)
dep_edges_335 <- clean_edges_335(combined_edges_selected, "augmented")
int_edges_335 <- clean_edges_335(sex_all_plot_selected,   "sexual")

# N=618 augmented edges (for covariate cross-unit terms, Section 5)
dep_edges_618 <- clean_edges_618(combined_edges, "augmented")

# Sanity check: sexual edges should be subset of augmented
dep_set <- paste(dep_edges_335$i, dep_edges_335$j, sep="_")
int_set <- paste(int_edges_335$i, int_edges_335$j, sep="_")
n_miss  <- sum(!(int_set %in% dep_set))
if (n_miss > 0)
  cat(sprintf("  WARNING: %d sexual edges not found in augmented network.\n", n_miss))

## ============================================================================
## 3. BUILD ADJACENCY LISTS AND WEIGHTS
## ============================================================================

# ---- 3a. N=335 dependence neighbors (augmented) — used for Y-Y and own L update ----
nbr_dep <- vector("list", N_total)
for (i in seq_len(N_total)) nbr_dep[[i]] <- integer(0)
for (r in seq_len(nrow(dep_edges_335))) {
  i <- dep_edges_335$i[r]; j <- dep_edges_335$j[r]
  nbr_dep[[i]] <- c(nbr_dep[[i]], j)
  nbr_dep[[j]] <- c(nbr_dep[[j]], i)
}

# ---- 3b. N=335 interference neighbors (sexual) — used for A and L neighbor terms ----
nbr_int <- vector("list", N_total)
for (i in seq_len(N_total)) nbr_int[[i]] <- integer(0)
for (r in seq_len(nrow(int_edges_335))) {
  i <- int_edges_335$i[r]; j <- int_edges_335$j[r]
  nbr_int[[i]] <- c(nbr_int[[i]], j)
  nbr_int[[j]] <- c(nbr_int[[j]], i)
}

# ---- 3c. N=618 augmented neighbors FOR EACH N=335 unit (Section 5 covariate term) ----
# Returns su_id strings of N=618 neighbors (may include units outside N=335)
nbr_dep618_of_335 <- vector("list", N_total)
for (i in seq_len(N_total)) nbr_dep618_of_335[[i]] <- character(0)

suid335_set <- as.character(node335[[NODE_ID_COL]])

for (r in seq_len(nrow(dep_edges_618))) {
  su_f <- dep_edges_618$su_from[r]
  su_t <- dep_edges_618$su_to[r]
  # Only fill in for units that are in N=335
  if (su_f %in% suid335_set) {
    idx_f <- suid2idx_335[su_f]
    nbr_dep618_of_335[[idx_f]] <- c(nbr_dep618_of_335[[idx_f]], su_t)
  }
  if (su_t %in% suid335_set) {
    idx_t <- suid2idx_335[su_t]
    nbr_dep618_of_335[[idx_t]] <- c(nbr_dep618_of_335[[idx_t]], su_f)
  }
}
# Remove duplicates
nbr_dep618_of_335 <- lapply(nbr_dep618_of_335, unique)

# ---- 3d. Designated weights ----

# Interference weights w^a and w^l: uniform over sexual neighbors (row-normalized)
W_int_val <- lapply(nbr_int, function(nb) {
  if (length(nb) == 0L) numeric(0) else rep(1/length(nb), length(nb))
})

# Y-Y dependence weights w^y: symmetric 0.5 per paper simulation
wY_val <- lapply(nbr_dep, function(nb) {
  if (length(nb) == 0L) numeric(0) else rep(0.5, length(nb))
})

# Covariate cross-unit weights v_ij: uniform over N=618 augmented neighbors
v_618_val <- lapply(nbr_dep618_of_335, function(nb) {
  if (length(nb) == 0L) numeric(0) else rep(1/length(nb), length(nb))
})

## ============================================================================
## 4. COMPONENT DIAGNOSTICS (informational only — we run all N=335 jointly)
## ============================================================================

g_dep <- graph_from_data_frame(
  d        = dep_edges_335,
  vertices = data.frame(name = as.character(seq_len(N_total))),
  directed = FALSE
)
comps   <- components(g_dep)
cs_tab  <- table(as.integer(comps$csize))
cat("\n=== Augmented network component sizes (N=335, diagnostic only) ===\n")
print(data.frame(
  size      = as.integer(names(cs_tab)),
  n_comps   = as.integer(cs_tab),
  n_units   = as.integer(names(cs_tab)) * as.integer(cs_tab)
), row.names = FALSE)
cat(sprintf("  Isolates: %d  |  Largest component: %d\n",
            sum(comps$csize == 1), max(comps$csize)))

## ============================================================================
## 5. JOINT GIBBS SAMPLER — CG(c), Multiple Covariates, N=335 as one block
##
## Covariate update  : Eq. 5  — binary -> logistic, continuous -> Gaussian
##                     Cross-unit term uses N=618 augmented neighbors (Section 5)
## Outcome update    : Eq. 6  — beta2 and beta4 are length-p dot products
## ============================================================================

gibbs_psi_full <- function(a_alloc, n_iter, burn_in, seed) {
  set.seed(seed)
  N    <- N_total
  
  # Initialize from observed values
  L_cur <- L_mat335_obs          # N x p matrix
  Y_cur <- integer(N)
  
  acc  <- numeric(N)
  kept <- 0L
  
  for (m in 1:n_iter) {
    
    ## -- Covariate update: cycle through k = 1..p, then i = 1..N --
    ## Implements Algorithm 1 line 4: L^(m+1/k)_{-ki}
    for (k in 1:p) {
      for (i in 1:N) {
        
        # Within-unit term: sum_s!=k rho[k,s] * L_cur[i,s]  (zero if rho_mat = 0)
        within_term <- if (p > 1) sum(rho_mat[k, -k] * L_cur[i, -k]) else 0
        
        # Cross-unit term from N=618 augmented neighbors (Section 5, Eq. 4 expanded)
        # sum_j v_ij * sum_s omega[k,s] * L618[j, s]
        nb618 <- nbr_dep618_of_335[[i]]          # su_id strings
        cross_term <- if (length(nb618) > 0L) {
          v      <- v_618_val[[i]]               # uniform weights
          L_nb   <- L_mat618[nb618, , drop=FALSE] # neighbor rows from N=618 matrix
          # For each neighbor j: omega[k,] %*% L_nb[j,]
          # Then weight by v_ij and sum
          sum(v * as.numeric(L_nb %*% omega_mat[k, ]))
        } else 0
        
        lin_k <- eta[k] + within_term + cross_term
        
        if (COV_TYPES[k] == "binary") {
          L_cur[i, k] <- rbinom(1, 1, logistic(lin_k))
        } else {
          # Gaussian full conditional (Section 3.1, continuous H_ki)
          L_cur[i, k] <- rnorm(1, mean = lin_k, sd = sigma_cov[k])
        }
      }
    }
    
    ## -- Outcome update: Eq. 6 with vector beta2 and beta4 --
    for (i in 1:N) {
      nb_d <- nbr_dep[[i]]   # augmented dependence neighbors (N=335)
      nb_i <- nbr_int[[i]]   # sexual interference neighbors  (N=335)
      
      # Own covariate contribution: beta2 dot L_i  (Eq. 6: beta2^T * l_i)
      own_cov <- sum(beta2 * L_cur[i, ])
      
      # Neighbor treatment + covariate contributions (interference, sexual network)
      if (length(nb_i) > 0L) {
        wi       <- W_int_val[[i]]
        sumA_int <- sum(wi * a_alloc[nb_i])
        # wL_vec[s] = sum_j w^l_ij * L_cur[j, s]  for each s -> length-p vector
        wL_vec   <- colSums(wi * L_cur[nb_i, , drop=FALSE])
        sumL_int <- sum(beta4 * wL_vec)          # beta4^T * wL_vec
      } else {
        sumA_int <- 0
        sumL_int <- 0
      }
      
      # Y-Y dependence (augmented network, N=335)
      sumY_dep <- if (length(nb_d) > 0L)
        sum(wY_val[[i]] * theta * Y_cur[nb_d]) else 0
      
      lin <- beta0 + beta1 * a_alloc[i] + own_cov +
        beta3 * sumA_int + sumL_int + sumY_dep
      
      Y_cur[i] <- rbinom(1, 1, logistic(lin))
    }
    
    if (m > burn_in) { acc <- acc + Y_cur; kept <- kept + 1L }
  }
  
  acc / kept
}

## ============================================================================
## 6. BASELINE: psi at A = 0 for all N=335 (single run, shared across replicates)
## ============================================================================

cat("\n=== Baseline Gibbs run (A = 0 for all, single run) ===\n")
t0 <- Sys.time()
psi_baseline <- gibbs_psi_full(rep(0L, N_total), n_iter, burn_in, seed = 999L)
cat(sprintf("  Done in %.1f s  |  mean psi_baseline = %.4f\n",
            as.numeric(Sys.time() - t0, units = "secs"), mean(psi_baseline)))

## ============================================================================
## 7. ONE REPLICATE — 3 Gibbs runs
##
## Per replicate r (following paper Section 3.3 and Table 3 approach):
##   A_treat = 1 for all 335 -> psi_treat
##   A_ctrl  = 0 for all 335 -> psi_ctrl   (re-run with replicate seed)
##   DE_i^(r) = psi_treat_i - psi_ctrl_i
##   IE_i^(r) = psi_ctrl_i  - psi_baseline_i
##   TE_i^(r) = DE_i^(r) + IE_i^(r)
##   Network-level: mean over N_total units
## ============================================================================

one_replicate <- function(run_id) {
  base_seed <- 10000L + run_id * 1000L
  
  A_treat <- rep(1L, N_total)
  A_ctrl  <- rep(0L, N_total)
  
  psi_treat <- gibbs_psi_full(A_treat, n_iter, burn_in, seed = base_seed + 1L)
  psi_ctrl  <- gibbs_psi_full(A_ctrl,  n_iter, burn_in, seed = base_seed + 2L)
  
  DE_i <- psi_treat - psi_ctrl
  IE_i <- psi_ctrl  - psi_baseline
  TE_i <- DE_i + IE_i
  
  c(DE             = mean(DE_i),
    IE             = mean(IE_i),
    TE             = mean(TE_i),
    psi_treat_mean = mean(psi_treat),
    psi_ctrl_mean  = mean(psi_ctrl))
}

## ============================================================================
## 8. RUN 40 REPLICATES IN PARALLEL
## ============================================================================

n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("\n=== Running %d replicates (%d cores, 2 Gibbs runs each) ===\n",
            n_runs, n_cores))

t0 <- Sys.time()

if (.Platform$OS.type != "windows") {
  results_list <- mclapply(
    1:n_runs, one_replicate,
    mc.cores = n_cores, mc.preschedule = TRUE, mc.set.seed = FALSE
  )
} else {
  cl <- makeCluster(n_cores)
  clusterExport(cl, varlist = c(
    # data
    "N_total", "L_mat335_obs", "L_mat618",
    # neighbor lists and weights
    "nbr_dep", "nbr_int", "nbr_dep618_of_335",
    "W_int_val", "wY_val", "v_618_val",
    # parameters
    "beta0","beta1","beta2","beta3","beta4",
    "theta","eta","sigma_cov","omega_mat","rho_mat",
    "alpha","n_iter","burn_in",
    # covariate metadata
    "COV_TYPES","p",
    # functions and shared objects
    "logistic","gibbs_psi_full","psi_baseline"
  ), envir = environment())
  results_list <- parLapply(cl, 1:n_runs, one_replicate)
  stopCluster(cl)
}

elapsed <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("  Done in %.1f s (%.1f min)\n", elapsed, elapsed / 60))

results_mat <- do.call(rbind, results_list)

## ============================================================================
## 9. SUMMARY TABLE  (matches Table 3 in paper)
## ============================================================================

summarise_one <- function(name, vals) {
  m     <- mean(vals)
  emp_se <- sd(vals)
  data.frame(
    Estimand = name,
    Mean     = round(m, 6),
    Emp_SE   = round(emp_se, 6),
    CI_Lower = round(m - 1.96 * emp_se, 6),
    CI_Upper = round(m + 1.96 * emp_se, 6)
  )
}

summary_tbl <- rbind(
  summarise_one("DE", results_mat[, "DE"]),
  summarise_one("IE", results_mat[, "IE"]),
  summarise_one("TE", results_mat[, "TE"])
)

cat("\n============  uConnect DE / IE / TE  (N=335, 40 replicates)  ============\n")
print(summary_tbl, row.names = FALSE)

cat("\n--- Diagnostic: mean psi under each allocation ---\n")
diag_df <- data.frame(
  Allocation = c("A = 1 (treat all)", "A = 0 (control all)", "A = 0 (baseline, fixed)"),
  Mean_psi   = round(c(mean(results_mat[, "psi_treat_mean"]),
                       mean(results_mat[, "psi_ctrl_mean"]),
                       mean(psi_baseline)), 4)
)
print(diag_df, row.names = FALSE)

