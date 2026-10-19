# Mo uConnect: descriptive outcome-observed versus outcome-missing comparison.
# Source this file to load configuration/functions; then explicitly call
# run_missingness_check(). No packages, model fitting, imputation, or tests of MCAR.
# Keep inputs and generated aggregate outputs inside Stronghold.

DATA_DIR <- "E:/Data/Project_Data/R01_MSM_HIV_Networks/uConnect/URI's uConnect clean data"
RISK_DIR <- "E:/Projects/Team_Members/resources/HIV riskbehav stuff Mo"
NODE_FILE <- file.path(DATA_DIR, "ego_nodechars_w1.csv")
EDGE_FILE <- file.path(DATA_DIR, "uconnect_sociometric_edgelist.csv")
W1_FILE <- file.path(RISK_DIR, "wave1_hivrb.csv")
W3_FILE <- file.path(RISK_DIR, "wave3_hivrb.csv")
# If automatic resolution stops, inspect the Wave3 FILE and set its exact field.
W3_OUTCOME_COLUMN <- NULL

.clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}
.require_columns <- function(x, columns, label) {
  if (!is.data.frame(x) || anyDuplicated(names(x)))
    stop(label, ": expected a data frame with unique column names.", call. = FALSE)
  absent <- setdiff(columns, names(x))
  if (length(absent)) stop(label, ": missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
}
.unique_ids <- function(x, label) {
  if (anyNA(x)) stop(label, ": missing/blank IDs are not allowed.", call. = FALSE)
  if (anyDuplicated(x)) stop(label, ": duplicate IDs; resolve upstream, do not merge arbitrarily.", call. = FALSE)
}
.number <- function(x, label, binary = FALSE) {
  text <- .clean_text(x)
  value <- suppressWarnings(as.numeric(text))
  bad <- !is.na(text) & (is.na(value) | !is.finite(value))
  if (binary) bad <- bad | (!is.na(value) & !(value %in% c(0, 1)))
  if (any(bad)) stop(label, ": ", sum(bad), " invalid ",
                     if (binary) "binary (0/1)" else "numeric", " values.", call. = FALSE)
  value
}
.resolve_w3 <- function(risk_w3, column) {
  if (is.null(column)) {
    candidates <- intersect(c("hivriskbehav_v4_w3", "hivriskbehav_v4_w1"), names(risk_w3))
    if (length(candidates) != 1L)
      stop("Wave3 outcome schema is ambiguous or missing. Inspect wave3_hivrb.csv and edit W3_OUTCOME_COLUMN to its exact derived outcome column.", call. = FALSE)
    column <- candidates
  }
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column) || !column %in% names(risk_w3))
    stop("Edit W3_OUTCOME_COLUMN to one exact column present in the Wave3 FILE.", call. = FALSE)
  column
}
# Each simple undirected pair appears once; no group-induced graphs.
.simple_pairs <- function(from, to) {
  data.frame(from = pmin(from, to), to = pmax(from, to), stringsAsFactors = FALSE) |>
    unique()
}
.degree <- function(pairs, ids) {
  tabulate(match(c(pairs$from, pairs$to), ids), nbins = length(ids))
}
.abs_smd <- function(x, observed, binary = FALSE) {
  a <- x[observed]; b <- x[!observed]
  difference <- abs(mean(a) - mean(b))
  if (binary) {
    pooled <- sqrt((mean(a) * (1 - mean(a)) + mean(b) * (1 - mean(b))) / 2)
  } else {
    # Sample SD is undefined for a singleton group; do not invent a variance.
    if (length(a) < 2L || length(b) < 2L) return(NA_real_)
    pooled <- sqrt((stats::var(a) + stats::var(b)) / 2)
  }
  if (pooled == 0) return(if (difference == 0) 0 else Inf)
  difference / pooled
}

build_missingness_table <- function(node, edge, risk_w1, risk_w3,
                                    w3_outcome_column = NULL) {
  .require_columns(node, c("su_id", "UniqueID1", "pepknow", "agecalc", "numpars1",
                           "healthcov", "drugsmarij", "jailever", "empstat"), "Node")
  .require_columns(edge, c("uniqueid1_from", "uniqueid1_to", "origsource_num"), "Edge")
  .require_columns(risk_w1, c("su_id", "hivriskbehav_v4_w1"), "Wave1 risk")
  .require_columns(risk_w3, "su_id", "Wave3 risk")
  selected_w3 <- .resolve_w3(risk_w3, w3_outcome_column)
  raw_node_n <- nrow(node)
  raw_edge_n <- nrow(edge)
  node[] <- lapply(node, .clean_text)
  no_graph_id <- is.na(node$UniqueID1)
  node <- node[!no_graph_id, , drop = FALSE]
  # Keep-first is allowed only for identical duplicate graph-ID rows. Check ALL
  # supplied node fields (after trimming), not just the covariates used below.
  duplicate_ids <- unique(node$UniqueID1[duplicated(node$UniqueID1)])
  conflicts <- vapply(duplicate_ids, function(id) {
    nrow(unique(node[node$UniqueID1 == id, , drop = FALSE])) > 1L
  }, logical(1))
  if (any(conflicts)) stop("Node: conflicting duplicate UniqueID1 records; resolve upstream.", call. = FALSE)
  duplicate_rows <- sum(duplicated(node$UniqueID1))
  node <- node[!duplicated(node$UniqueID1), , drop = FALSE]
  .unique_ids(node$su_id, "Enrolled analysis su_id")
  risk_w1$su_id <- .clean_text(risk_w1$su_id)
  risk_w3$su_id <- .clean_text(risk_w3$su_id)
  .unique_ids(risk_w1$su_id, "Wave1 risk su_id")
  .unique_ids(risk_w3$su_id, "Wave3 risk su_id")
  r1 <- .number(risk_w1$hivriskbehav_v4_w1, "Wave1 risk", TRUE)
  r3 <- .number(risk_w3[[selected_w3]], "Wave3 risk", TRUE)
  m1 <- match(node$su_id, risk_w1$su_id)
  m3 <- match(node$su_id, risk_w3$su_id)
  age <- .number(node$agecalc, "Age")
  partners <- .number(node$numpars1, "Partners")
  healthcare <- .number(node$healthcov, "Healthcare", TRUE)
  marijuana <- .number(node$drugsmarij, "Marijuana", TRUE)
  jail <- .number(node$jailever, "Criminal justice", TRUE)
  employment <- rep(NA_real_, nrow(node))
  employment[node$empstat %in% c("1", "2")] <- 1
  employment[node$empstat %in% "3"] <- 0
  treatment <- .number(node$pepknow, "Treatment pepknow", TRUE)
  baseline <- data.frame(age, partners, healthcare, marijuana, jail, employment,
                         baseline_risk = r1[m1], treatment)
  full <- stats::complete.cases(baseline)
  # IMPORTANT: Y binds only to the actual Wave3 data frame, NEVER risk_w1.
  y <- r3[m3][full]
  observed <- !is.na(y)
  n_observed <- sum(observed); n_missing <- sum(!observed)
  if (!n_observed || !n_missing)
    stop("Full cohort must contain both outcome-observed and outcome-missing participants; one or both groups are empty.", call. = FALSE)
  x <- baseline[full, , drop = FALSE]
  full_ids <- node$UniqueID1[full]

  from <- .clean_text(edge$uniqueid1_from)
  to <- .clean_text(edge$uniqueid1_to)
  type <- .clean_text(edge$origsource_num)
  valid <- type %in% c("1", "2") & !is.na(from) & !is.na(to)
  valid <- valid & !is.na(from == to) & from != to
  valid[is.na(valid)] <- FALSE
  enrolled_edges <- valid & from %in% node$UniqueID1 & to %in% node$UniqueID1
  full_edges <- enrolled_edges & from %in% full_ids & to %in% full_ids
  augmented <- .simple_pairs(from[full_edges], to[full_edges])
  sexual_rows <- full_edges & type %in% "2"
  sexual <- .simple_pairs(from[sexual_rows], to[sexual_rows])
  x$sexual_degree <- .degree(sexual, full_ids)
  x$augmented_degree <- .degree(augmented, full_ids)

  count_percent <- function(z) sprintf("%d/%d (%.1f%%)", sum(z), length(z), 100 * mean(z))
  mean_sd <- function(z) sprintf("%.2f (%.2f)", mean(z), stats::sd(z))
  median_iqr <- function(z) {
    q <- stats::quantile(z, c(.25, .5, .75), names = FALSE)
    sprintf("%.2f [%.2f, %.2f]", q[2], q[1], q[3])
  }
  rows <- list()
  add <- function(variable, summary, a, b, smd = NA_real_) {
    data.frame(Variable = variable, Level_or_summary = summary,
               Outcome_observed = a, Outcome_missing = b, Absolute_SMD = smd,
               stringsAsFactors = FALSE)
  }
  rows[[1]] <- add("Group size", "N", as.character(n_observed), as.character(n_missing))
  quantitative <- c(age = "Age (years)", partners = "Sexual partners (past 6 months)",
                    sexual_degree = "Sexual degree", augmented_degree = "Augmented degree")
  for (key in names(quantitative)) {
    z <- x[[key]]; label <- quantitative[[key]]
    rows[[length(rows) + 1L]] <- add(label, "Mean (SD)", mean_sd(z[observed]),
                                    mean_sd(z[!observed]), .abs_smd(z, observed))
    rows[[length(rows) + 1L]] <- add(label, "Median [Q1, Q3] (IQR limits)",
                                    median_iqr(z[observed]), median_iqr(z[!observed]))
    if (key %in% c("sexual_degree", "augmented_degree"))
      rows[[length(rows) + 1L]] <- add(label, "Isolate (degree = 0), n/N (%)",
                                      count_percent(z[observed] == 0), count_percent(z[!observed] == 0),
                                      .abs_smd(as.numeric(z == 0), observed, TRUE))
  }
  binaries <- c(healthcare = "Healthcare coverage", marijuana = "Cannabis use history",
                jail = "Criminal justice involvement (jailever)", employment = "Employment (1/2 -> 1; 3 -> 0)",
                baseline_risk = "Baseline HIV risk behavior", treatment = "PrEP awareness")
  for (key in names(binaries)) {
    z <- x[[key]]
    for (level in c(0, 1))
      rows[[length(rows) + 1L]] <- add(binaries[[key]], paste0(level, ", n/N (%)"),
                                      count_percent(z[observed] == level), count_percent(z[!observed] == level),
                                      if (level == 1) .abs_smd(z, observed, TRUE) else NA_real_)
  }
  rows[[length(rows) + 1L]] <- add("Wave3 outcome positive", "1, n/N (%)",
                                  count_percent(y[observed] == 1), "Not observed")
  table <- do.call(rbind, rows)
  rownames(table) <- NULL

  actual <- c(raw_node_rows = raw_node_n,
              node_rows_missing_graph_id_removed = sum(no_graph_id),
              identical_duplicate_graph_rows_removed = duplicate_rows,
              enrolled = nrow(node), excluded_baseline_or_treatment = sum(!full),
              enrolled_missing_healthcare = sum(is.na(healthcare)),
              full = sum(full), outcome_observed = n_observed, outcome_missing = n_missing,
              wave1_raw_rows = nrow(risk_w1), wave3_raw_rows = nrow(risk_w3),
              wave1_raw_rows_unmatched_to_enrolled = sum(!risk_w1$su_id %in% node$su_id),
              wave3_raw_rows_unmatched_to_enrolled = sum(!risk_w3$su_id %in% node$su_id),
              wave1_raw_rows_unmatched_to_full = sum(!risk_w1$su_id %in% node$su_id[full]),
              wave3_raw_rows_unmatched_to_full = sum(!risk_w3$su_id %in% node$su_id[full]),
              enrolled_unmatched_wave1_record = sum(is.na(m1)),
              enrolled_unmatched_wave3_record = sum(is.na(m3)),
              enrolled_matched_wave3_missing_field = sum(!is.na(m3) & is.na(r3[m3])),
              full_unmatched_wave3_record = sum(is.na(m3[full])),
              full_matched_wave3_missing_field = sum(!is.na(m3[full]) & is.na(y)),
              raw_edge_rows = raw_edge_n, allowed_nonblank_nonself_edge_rows = sum(valid),
              enrolled_edge_rows_before_pair_dedup = sum(enrolled_edges),
              full_edge_rows_before_pair_dedup = sum(full_edges),
              full_sexual_undirected_edges = nrow(sexual),
              full_augmented_undirected_edges = nrow(augmented))
  expected <- c(enrolled = 618, full = 586, outcome_observed = 468, outcome_missing = 118)
  manuscript <- unname(expected[names(actual)])
  flow <- data.frame(Metric = names(actual), Actual = unname(actual), Manuscript_expected = manuscript,
                     Difference = unname(actual) - manuscript,
                     Status = ifelse(is.na(manuscript), "Not specified",
                                     ifelse(unname(actual) == manuscript, "Matches", "DEVIATION")),
                     stringsAsFactors = FALSE)
  notes <- c(
    paste("Wave3 selected column:", selected_w3),
    "Y uses the supplied Wave3 data frame only; baseline risk is never substituted into Y.",
    "The selected derived field is supplied by Mo. Its components and skip rules have NOT been reconstructed or verified.",
    "Enrollment: trim IDs; remove blank/missing UniqueID1; keep first identical graph-ID row, fail on conflicting duplicates; require unique nonmissing analysis IDs and risk IDs.",
    "Full cohort: complete age, partners, healthcare, marijuana, jail, recoded employment, baseline risk, and treatment. Missing healthcare is EXCLUDED.",
    "Employment: codes 1/2 -> 1, 3 -> 0; all other/missing codes -> NA. Age/partners retain original numeric scales.",
    "No baseline-HIV-negative restriction is applied. Outcome-observed means nonmissing actual Wave3 risk; missing includes unmatched records and matched missing fields.",
    "Unmatched Wave3 record versus matched missing field is a record-availability distinction, NOT a verified attrition cause.",
    "Raw risk unmatched-to-full counts include records belonging to excluded enrolled participants as well as records outside enrollment.",
    sprintf("Denominators: observed N=%d; missing N=%d. Binary rows show n/N (%%); baseline/treatment have no missing values in full.", n_observed, n_missing),
    "Outcome positivity denominator is the outcome-observed group only; missing group is Not observed, with no SMD.",
    "Sexual edges use type 2; augmented edges union types 1/2. Remove invalid endpoints/self-loops, restrict to enrolled then full, deduplicate unordered pairs, retain all full-cohort isolates.",
    "Both groups use degrees from the SAME full-cohort graph, including between-group edges. Full edge totals are cohort totals, NOT group edge counts.",
    "Continuous SMD = absolute mean difference / sqrt((sample variance observed + sample variance missing)/2). Binary SMD uses p*(1-p) variances; reported once at level 1.",
    "Zero pooled variance: equal means -> 0; different means -> Inf. Continuous SMD/SD with a singleton group is NA. Blank CSV SMD means not applicable/undefined.",
    "Median [Q1, Q3] reports the endpoints of the interquartile interval (R default quantile type 7).",
    "These are descriptive comparisons only: no fitting, imputation, significance tests/p-values, or proof/test of MCAR. SMD does not account for network dependence.",
    "Manuscript counts are reference checks only, never forced. Inspect every DEVIATION in cohort_flow.csv.",
    "Only aggregate outputs are returned/saved. No participant, edge, ID-level data, input snapshots, or original prep source are saved."
  )
  list(table = table, cohort_flow = flow, notes = notes)
}

run_missingness_check <- function(node_file = NODE_FILE, edge_file = EDGE_FILE,
                                  w1_file = W1_FILE, w3_file = W3_FILE,
                                  w3_outcome_column = W3_OUTCOME_COLUMN,
                                  output_dir = file.path(getwd(), "uconnect_missingness_check",
                                    format(Sys.time(), "%Y%m%d_%H%M%OS6"))) {
  # Never source the original prep or consume saved obs_y.
  paths <- c(node = node_file, edge = edge_file, wave1 = w1_file, wave3 = w3_file)
  if (any(!file.exists(paths))) stop("Input file(s) missing: ", paste(names(paths)[!file.exists(paths)], collapse = ", "), call. = FALSE)
  paths <- vapply(paths, normalizePath, character(1), winslash = "/", mustWork = TRUE)
  if (paths[["wave1"]] == paths[["wave3"]])
    stop("Wave1 and Wave3 must be distinct input files; refusing baseline substitution.", call. = FALSE)
  if (file.exists(output_dir) || dir.exists(output_dir))
    stop("Output directory already exists; refusing to overwrite: ", output_dir, call. = FALSE)
  read_input <- function(path) read.csv(path, colClasses = "character", stringsAsFactors = FALSE,
                                        check.names = FALSE, na.strings = c("", "NA"))
  node <- read_input(paths[["node"]]); edge <- read_input(paths[["edge"]])
  risk_w1 <- read_input(paths[["wave1"]]); risk_w3 <- read_input(paths[["wave3"]])
  selected <- .resolve_w3(risk_w3, w3_outcome_column)
  message("Wave3 outcome FILE: ", paths[["wave3"]], "\nWave3 outcome COLUMN: ", selected)
  result <- build_missingness_table(node, edge, risk_w1, risk_w3, selected)
  provenance <- c("Mo uConnect aggregate missingness comparison", paste("Run time:", Sys.time()),
                  paste("R:", R.version.string), paste(names(paths), "FILE:", paths),
                  paste("Wave3 outcome binding:", paths[["wave3"]], "::", selected),
                  result$notes, "", "Actual counts versus manuscript:",
                  capture.output(print(result$cohort_flow, row.names = FALSE)))
  # Create parent as needed, then atomically require a NEW final directory.
  dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(output_dir, showWarnings = FALSE))
    stop("Could not create a new output directory; nothing overwritten.", call. = FALSE)
  write.csv(result$table, file.path(output_dir, "table.csv"), row.names = FALSE, na = "")
  write.csv(result$cohort_flow, file.path(output_dir, "cohort_flow.csv"), row.names = FALSE, na = "")
  writeLines(enc2utf8(provenance), file.path(output_dir, "source_reconciliation.txt"), useBytes = TRUE)
  deviations <- result$cohort_flow[result$cohort_flow$Status == "DEVIATION", , drop = FALSE]
  if (nrow(deviations)) {
    message("MANUSCRIPT COUNT DEVIATIONS (not forced):")
    print(deviations, row.names = FALSE)
  }
  message("Aggregate outputs saved inside: ", normalizePath(output_dir, winslash = "/"))
  invisible(result)
}

message("Loaded functions/configuration only. Review paths and W3_OUTCOME_COLUMN, then run: run_missingness_check()")
