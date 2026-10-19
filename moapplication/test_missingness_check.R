# Synthetic-only acceptance tests. Run using Rscript /absolute/path/to/this/file.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[startsWith(args, "--file=")]
if (length(file_arg) != 1L) stop("Run this test file with Rscript.")
build_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
script <- file.path(build_dir, "missingness_check.R")
workspace <- tempfile("synthetic_tests_", tmpdir = build_dir)
dir.create(workspace)
old_wd <- getwd()
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAIL: ", label, call. = FALSE)
  checks <<- checks + 1L
  cat("PASS: ", label, "\n", sep = "")
}
expect_error <- function(expr, pattern, label) {
  text <- tryCatch({ force(expr); "NO ERROR" }, error = function(e) conditionMessage(e))
  check(grepl(pattern, text, fixed = TRUE), paste0(label, " [", pattern, "]"))
}

run_tests <- function() {
  on.exit({setwd(old_wd); unlink(workspace, recursive = TRUE)}, add = TRUE)
  setwd(workspace)
  before <- list.files(all.files = TRUE)
  source(script, local = .GlobalEnv)
  check(identical(before, list.files(all.files = TRUE)), "absolute-path source from different cwd writes nothing")
  node <- data.frame(
    su_id = sprintf("%03d", 1:8), UniqueID1 = paste0("g", 1:8),
    pepknow = c("0", "1", "0", "1", "0", "1", "0", "1"),
    agecalc = as.character(seq(20, 90, 10)),
    numpars1 = c(rep("1", 3), rep("0", 5)),
    healthcov = c(rep("1", 6), NA, "1"),
    drugsmarij = rep(c("0", "1"), 4),
    jailever = c(rep("1", 3), rep("0", 5)),
    empstat = c("1", "2", "3", "1", "2", "3", "1", "9"),
    hiv_positive = c("0", "0", "1", rep("0", 5)), stringsAsFactors = FALSE)
  node$su_id[1] <- " 001 "
  node$UniqueID1[1] <- " g1 "
  blank <- node[8, , drop = FALSE]; blank$UniqueID1 <- " "
  missing <- blank; missing$UniqueID1 <- NA
  node <- rbind(node, node[1, , drop = FALSE], blank, missing)
  risk_w1 <- data.frame(su_id = c(sprintf("%03d", 1:8), "999"),
                        hivriskbehav_v4_w1 = c("0", "1", "0", "1", "0", "1", "0", "1", "0"))
  risk_w3 <- data.frame(su_id = c("001", "002", "003", "005", "007", "999"),
                        hivriskbehav_v4_w3 = c("1", "0", "1", NA, "1", "0"))
  edge <- data.frame(
    uniqueid1_from = c("g1", "g4", "g1", "g2", "g1", "g4", "g1", "g3", "g1", "g1", "", "g3", "g3", "g2"),
    uniqueid1_to   = c("g4", "g1", "g4", "g5", "g2", "g5", "g5", "g3", "g7", "outside", "g3", "g6", "g6", "g5"),
    origsource_num = c("2", "2", "1", "2", "1", "1", "1", "2", "2", "2", "2", "3", NA, "2"))
  build <- function(n = node, e = edge, w1 = risk_w1, w3 = risk_w3, column = NULL)
    build_missingness_table(n, e, w1, w3, column)
  result <- build()
  metric <- function(name, r = result) r$cohort_flow$Actual[match(name, r$cohort_flow$Metric)]
  cell <- function(variable, summary, column, r = result) {
    tab <- r$table
    tab[tab$Variable == variable & tab$Level_or_summary == summary, column]
  }
  check(identical(names(result), c("table", "cohort_flow", "notes")), "builder returns aggregate tables/notes only")
  check(identical(names(result$table), c("Variable", "Level_or_summary", "Outcome_observed", "Outcome_missing", "Absolute_SMD")), "required table schema")
  check(metric("raw_node_rows") == 11 && metric("enrolled") == 8 &&
          metric("identical_duplicate_graph_rows_removed") == 1 &&
          metric("node_rows_missing_graph_id_removed") == 2, "exact raw/enrolled reconciliation and identical keep-first")
  check(metric("full") == 6 && metric("enrolled_missing_healthcare") == 1 &&
          metric("excluded_baseline_or_treatment") == 2, "missing healthcare and other employment excluded")
  check(metric("outcome_observed") == 3 && metric("outcome_missing") == 3, "actual counts retained, including HIV-positive participant")
  check(cell("Employment (1/2 -> 1; 3 -> 0)", "1, n/N (%)", "Outcome_observed") == "2/3 (66.7%)" &&
          cell("Employment (1/2 -> 1; 3 -> 0)", "0, n/N (%)", "Outcome_missing") == "1/3 (33.3%)", "employment codes 1/2 and 3 correctly recoded")
  check(cell("Wave3 outcome positive", "1, n/N (%)", "Outcome_observed") == "2/3 (66.7%)" &&
          cell("Baseline HIV risk behavior", "1, n/N (%)", "Outcome_observed") == "1/3 (33.3%)", "distinct W3 values, not W1, determine outcome")
  check(cell("Wave3 outcome positive", "1, n/N (%)", "Outcome_missing") == "Not observed" &&
          is.na(cell("Wave3 outcome positive", "1, n/N (%)", "Absolute_SMD")), "outcome missing label and no outcome SMD")
  check(metric("full_unmatched_wave3_record") == 2 && metric("full_matched_wave3_missing_field") == 1,
        "unmatched Wave3 versus matched missing field separated")
  check(metric("wave3_raw_rows_unmatched_to_enrolled") == 1 &&
          metric("wave3_raw_rows_unmatched_to_full") == 2 && metric("wave1_raw_rows_unmatched_to_full") == 3,
        "raw risk unmatched counts include outside-full records")
  check(metric("full_sexual_undirected_edges") == 2 && metric("full_augmented_undirected_edges") == 5,
        "graph bridge IDs, duplicates/opposite/type/endpoint/selfloop filtering")
  check(cell("Sexual degree", "Mean (SD)", "Outcome_observed") == "0.67 (0.58)" &&
          cell("Sexual degree", "Mean (SD)", "Outcome_missing") == "0.67 (0.58)",
        "sexual degrees include between-group edges on common full graph")
  check(cell("Augmented degree", "Mean (SD)", "Outcome_observed") == "1.67 (1.53)" &&
          cell("Augmented degree", "Mean (SD)", "Outcome_missing") == "1.67 (1.53)",
        "union degrees deduplicate cross-type pairs")
  check(cell("Augmented degree", "Isolate (degree = 0), n/N (%)", "Outcome_missing") == "1/3 (33.3%)" &&
          cell("Sexual degree", "Isolate (degree = 0), n/N (%)", "Outcome_observed") == "1/3 (33.3%)",
        "isolates preserved in both groups")
  check(cell("Age (years)", "Mean (SD)", "Outcome_observed") == "30.00 (10.00)" &&
          cell("Age (years)", "Median [Q1, Q3] (IQR limits)", "Outcome_missing") == "60.00 [55.00, 65.00]" &&
          cell("Age (years)", "Mean (SD)", "Absolute_SMD") == 3, "original numeric scale and descriptive continuous SMD")
  check(cell("Healthcare coverage", "1, n/N (%)", "Absolute_SMD") == 0 &&
          is.infinite(cell("Criminal justice involvement (jailever)", "1, n/N (%)", "Absolute_SMD")) &&
          is.infinite(cell("Sexual partners (past 6 months)", "Mean (SD)", "Absolute_SMD")) &&
          !any(is.nan(result$table$Absolute_SMD)), "zero-variance SMD is 0/Inf, never NaN")
  check(all(result$cohort_flow$Status[result$cohort_flow$Metric %in%
              c("enrolled", "full", "outcome_observed", "outcome_missing")] == "DEVIATION"), "manuscript deviations visible, not forced")
  check(!any(grepl("001|g1", unlist(result), fixed = FALSE)), "no fixture IDs in returned aggregates")

  legacy <- risk_w3; names(legacy)[2] <- "hivriskbehav_v4_w1"
  check(identical(build(w3 = legacy)$table, result$table), "legacy w1-suffix allowed ONLY within Wave3 file")
  ambiguous <- risk_w3; ambiguous$hivriskbehav_v4_w1 <- "0"
  expect_error(build(w3 = ambiguous), "ambiguous or missing", "ambiguous W3 schema stops")
  check(identical(build(w3 = ambiguous, column = "hivriskbehav_v4_w3")$table, result$table), "explicit Wave3 schema selection resolves ambiguity")
  absent <- risk_w3; names(absent)[2] <- "other_field"
  expect_error(build(w3 = absent), "ambiguous or missing", "missing W3 schema stops")
  expect_error(build(column = "typo"), "exact column", "invalid explicit W3 column stops")
  expect_error(build(w1 = rbind(risk_w1, risk_w1[1, ])), "duplicate IDs", "duplicate Wave1 IDs fail")
  expect_error(build(w3 = rbind(risk_w3, risk_w3[1, ])), "duplicate IDs", "duplicate Wave3 IDs fail")
  bad <- node[1:8, ]; bad$su_id[2] <- "001"
  expect_error(build(n = bad), "duplicate IDs", "duplicate analysis IDs fail after trimming")
  bad <- node; bad$agecalc[9] <- "999"
  expect_error(build(n = bad), "conflicting duplicate", "conflicting graph-ID duplicate fails")
  bad <- risk_w3; bad$hivriskbehav_v4_w3[1] <- "2"
  expect_error(build(w3 = bad), "invalid binary", "invalid Wave3 binary fails")
  bad <- risk_w1; bad$hivriskbehav_v4_w1[1] <- "9"
  expect_error(build(w1 = bad), "invalid binary", "invalid baseline risk fails")
  bad <- node[1:8, ]; bad$healthcov[1] <- "2"
  expect_error(build(n = bad), "invalid binary", "invalid baseline binary fails")
  bad <- node[1:8, ]; bad$pepknow[1] <- "yes"
  expect_error(build(n = bad), "invalid binary", "invalid treatment fails")
  bad <- node[1:8, ]; bad$agecalc[1] <- "not numeric"
  expect_error(build(n = bad), "invalid numeric", "malformed age fails without silent exclusion")
  bad <- risk_w3; bad$su_id[1] <- " "
  expect_error(build(w3 = bad), "missing/blank IDs", "blank risk ID fails")
  bad <- risk_w3; bad$hivriskbehav_v4_w3[] <- NA
  expect_error(build(w3 = bad), "groups are empty", "no observed outcome stops")
  all_observed <- data.frame(su_id = sprintf("%03d", 1:8), hivriskbehav_v4_w3 = "1")
  expect_error(build(w3 = all_observed), "groups are empty", "no missing outcome stops")
  bad <- node[1:8, ]; bad$pepknow[] <- NA
  expect_error(build(n = bad), "groups are empty", "empty full cohort stops")
  singleton <- risk_w3; singleton$hivriskbehav_v4_w3[2:3] <- NA
  check(is.na(cell("Age (years)", "Mean (SD)", "Absolute_SMD", build(w3 = singleton))), "singleton continuous variance explicitly undefined")
  empty_graph <- build(e = edge[FALSE, ])
  check(metric("full_augmented_undirected_edges", empty_graph) == 0 &&
          cell("Sexual degree", "Isolate (degree = 0), n/N (%)", "Outcome_observed", empty_graph) == "3/3 (100.0%)",
        "empty network retains full cohort as isolates")
  # CSV integration: leading-zero IDs must survive read.csv colClasses=character.
  paths <- file.path(workspace, c("ego_nodechars_w1.csv", "uconnect_sociometric_edgelist.csv", "wave1_hivrb.csv", "wave3_hivrb.csv"))
  fixtures <- list(node, edge, risk_w1, risk_w3)
  for (i in seq_along(paths)) write.csv(fixtures[[i]], paths[i], row.names = FALSE, na = "")
  before_hash <- tools::md5sum(paths)
  integrated <- run_missingness_check(paths[1], paths[2], paths[3], paths[4])
  check(identical(integrated$table, result$table), "actual CSV runner preserves IDs and matches pure aggregate builder")
  root <- file.path(workspace, "uconnect_missingness_check")
  dirs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  check(length(dirs) == 1L && grepl("^[0-9]{8}_[0-9]{6}[.][0-9]+$", basename(dirs)), "default timestamped output under getwd")
  output <- dirs[1]
  check(identical(sort(list.files(output)), sort(c("table.csv", "cohort_flow.csv", "source_reconciliation.txt"))), "only requested aggregate artifacts saved; no snapshots")
  saved <- read.csv(file.path(output, "table.csv"), check.names = FALSE)
  roundtrip <- all.equal(saved, result$table, tolerance = 1e-12)
  if (!isTRUE(roundtrip)) print(roundtrip)
  check(isTRUE(roundtrip), "saved table round-trips within 1e-12 including Inf and blank SMD")
  saved_flow <- read.csv(file.path(output, "cohort_flow.csv"))
  check(isTRUE(all.equal(saved_flow, result$cohort_flow)), "saved cohort flow matches builder")
  provenance <- readLines(file.path(output, "source_reconciliation.txt"))
  check(any(grepl(normalizePath(paths[4]), provenance, fixed = TRUE)) &&
          any(grepl("Wave3 selected column: hivriskbehav_v4_w3", provenance, fixed = TRUE)) &&
          any(grepl("DEVIATION", provenance, fixed = TRUE)) &&
          any(grepl("NOT been reconstructed or verified", provenance, fixed = TRUE)), "saved exact Wave3 binding, deviations and definition limitation")
  check(identical(before_hash, tools::md5sum(paths)), "input CSV fixtures not modified")
  hashes <- tools::md5sum(list.files(output, full.names = TRUE))
  expect_error(run_missingness_check(paths[1], paths[2], paths[3], paths[4], output_dir = output),
               "refusing to overwrite", "existing directory rejected")
  check(identical(hashes, tools::md5sum(list.files(output, full.names = TRUE))), "no-overwrite preserves all output bytes")
  expect_error(run_missingness_check(paths[1], paths[2], paths[3], paths[3]),
               "distinct input files", "runner rejects same W1/W3 file")
  invisible(TRUE)
}
run_tests()
cat(sprintf("RESULT: %d checks passed; synthetic-only; no real/enclave data accessed.\n", checks))
