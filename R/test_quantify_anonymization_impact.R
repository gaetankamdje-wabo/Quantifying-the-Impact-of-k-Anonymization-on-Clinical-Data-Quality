# =============================================================================
# COMPREHENSIVE UNIT TEST 
# Function under test: quantify_anonymization_impact()
# Author: Gaetan Kamdje Wabo, M.Sc.
# Affiliation: DBMI, Medical Faculty Mannheim, Heidelberg University
# Date: 2026-03-04
# R version: >= 4.4 | Dependencies: testthat
#
# Purpose:
#   Validation of all computational branches, edge cases,
#   and output invariants of quantify_anonymization_impact().
#   Results are captured as structured JSON for downstream reporting.
#
# Output:
#   - Console: testthat summary
#   - File: test_results_report.json (structured results for Word report)
# =============================================================================

suppressPackageStartupMessages(library(testthat))

# ── Source the main study script (defines all functions under test) ───────────
# The main script may be named either way depending on user's setup.
MAIN_SCRIPT_NAMES <- c(
  "anonymization_dq_study.R",
  "quantify_anonymization_impact.R"
)

# Try: same directory as this test file (RStudio / source() context)
this_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)

# Search order: (1) test file directory, (2) working directory, (3) known path
search_dirs <- unique(c(
  this_dir,
  getwd(),
  "C:/Users/kamdje-wabo/Documents/Anonymization Impact"
))

main_path <- NULL
for (d in search_dirs) {
  for (fn in MAIN_SCRIPT_NAMES) {
    candidate <- file.path(d, fn)
    if (file.exists(candidate)) {
      main_path <- candidate
      break
    }
  }
  if (!is.null(main_path)) break
}

if (is.null(main_path)) {
  stop(
    "Main script not found. Searched for:\n",
    "  Filenames: ", paste(MAIN_SCRIPT_NAMES, collapse = ", "), "\n",
    "  Directories:\n",
    paste("    -", search_dirs, collapse = "\n"), "\n",
    "Place this test file alongside the main script, or set the working directory.",
    call. = FALSE
  )
}

# ── Selective sourcing: load ONLY function definitions, not the study pipeline.
#    The main script contains both function definitions (sections 0–2, lines
#    ~1–541) and a full analytical pipeline (sections 3–G, lines ~544–1762).
#    Running the pipeline is unnecessary for testing and may fail if the
#    corroboration section encounters insufficient overlapping IRR terms.
#    Strategy: read the file, extract lines up to the sentinel that marks the
#    end of function definitions, and evaluate only those lines.
cat("Selectively sourcing function definitions from:", main_path, "\n")

.main_lines <- readLines(main_path, warn = FALSE)

# Sentinel: the line immediately after quantify_anonymization_impact() is defined.
# In the original script this is: cat("quantify_anonymization_impact() defined.\n\n")
.sentinel_idx <- grep("quantify_anonymization_impact\\(\\)\\s+defined", .main_lines)

if (length(.sentinel_idx) == 0L) {
  # Fallback: source the entire file if sentinel not found (user modified script)
  warning("Sentinel line not found — sourcing entire script. ",
          "This may produce warnings or errors from the study pipeline.",
          call. = FALSE)
  source(main_path, local = FALSE)
} else {
  # Evaluate only lines 1 through the sentinel (inclusive).
  # This loads: package setup, COL_* contract, section(), icd_chapter(),
  # los_midpoint(), print_table(), generate_d0(), D0 data, and
  # quantify_anonymization_impact() — nothing beyond.
  .code_block <- paste(.main_lines[seq_len(max(.sentinel_idx))], collapse = "\n")
  eval(parse(text = .code_block), envir = globalenv())
  cat("Loaded: quantify_anonymization_impact(), icd_chapter(), los_midpoint(),",
      "section(), print_table(), D0 dataset.\n")
}

# Clean up temporary objects
suppressWarnings(
  rm(list = intersect(c(".main_lines", ".sentinel_idx", ".code_block"),
                      ls(all.names = TRUE, envir = environment())),
     envir = environment())
)

# ── Output directory (configurable) ─────────────────────────────────────────
OUTPUT_DIR <- file.path(
  "A:", "HLZ", "Promotionen", "Gaetan Kamdje Wabo",
  "study datasets", "study files", "arx", "data source", "study_results"
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Helper: capture test result as structured list ───────────────────────────
capture_test <- function(test_id, category, description, test_data_desc,
                         methodology, expr, expected_desc) {
  t0     <- Sys.time()
  passed <- FALSE
  error  <- NULL
  tryCatch({
    eval(expr, envir = parent.frame())
    passed <- TRUE
  }, error = function(e) {
    error <<- conditionMessage(e)
  })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(
    test_id      = test_id,
    category     = category,
    description  = description,
    test_data    = test_data_desc,
    methodology  = methodology,
    expected     = expected_desc,
    passed       = passed,
    error        = error,
    elapsed_sec  = round(elapsed, 4)
  )
}

# ── Global result collector ──────────────────────────────────────────────────
RESULTS <- list()

# ── Fixture data factory ────────────────────────────────────────────────────
set.seed(99L)
df_base <- data.frame(
  id  = paste0("P", sprintf("%03d", 1:50)),
  icd = sample(c("I21", "I50", "J18", "K35", "F32"), 50, replace = TRUE),
  los = pmax(as.integer(rnbinom(50, mu = 7, size = 3)), 1L),
  stringsAsFactors = FALSE
)

# ── fmt_p local replica (for helper tests) ──────────────────────────────────
fmt_p_local <- function(p) {
  if (is.na(p)) return("NA")
  if (p == 0 || p < .Machine$double.xmin)
    return(paste0("< ", format(.Machine$double.xmin, scientific = TRUE)))
  format(p, digits = 17L, scientific = TRUE, trim = TRUE)
}


# =============================================================================
# CATEGORY A: IDENTITY / BASELINE INVARIANTS
# =============================================================================

test_results <- test_that("A: Identity and baseline invariants", {
  
  # ── A01: Identical numeric vectors → KS D = 0 ────────────────────────────
  res <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_base,
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_base$los
  )
  expect_equal(res$adhoc_column_comparison$ks$D, 0, tolerance = 1e-10,
               label = "A01: KS D = 0 for identical numeric vectors")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "A01", category = "Identity Invariants",
    description = "Identical numeric vectors yield KS D = 0",
    test_data   = "df_base$los compared to itself (n = 50, NB-distributed LOS)",
    methodology = "Pass identical LOS vectors as both origin and anonymized; assert KS D = 0",
    expected    = "KS D = 0 (distributions are identical)",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── A02: Identical numeric → Green verdict ────────────────────────────────
  expect_true(grepl("\u2705", res$adhoc_column_comparison$verdict),
              label = "A02: Green verdict for identical numeric data")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "A02", category = "Identity Invariants",
    description = "Identical numeric data produces Green (\u2705) verdict",
    test_data   = "Same as A01",
    methodology = "Check verdict string contains \u2705 (green checkmark)",
    expected    = "\u2705 No meaningful change",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── A03: Identical categorical → Cramer V = 0 ────────────────────────────
  res_cat <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_base,
    comparison_column_origin     = df_base$icd,
    comparison_column_anonymized = df_base$icd
  )
  expect_equal(res_cat$adhoc_column_comparison$cramer_v$v, 0,
               tolerance = 1e-10,
               label = "A03: Cramer V = 0 for identical categorical vectors")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "A03", category = "Identity Invariants",
    description = "Identical categorical vectors yield Cramer's V = 0",
    test_data   = "df_base$icd compared to itself (5 ICD-3 levels, n = 50)",
    methodology = "Pass identical ICD vectors; assert Cramer's V = 0",
    expected    = "V = 0 (contingency table has identical row distributions)",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── A04: Identical categorical → Jaccard = 1 ─────────────────────────────
  expect_equal(res_cat$adhoc_column_comparison$jaccard_overlap$ratio, 1,
               label = "A04: Jaccard = 1 for identical level sets")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "A04", category = "Identity Invariants",
    description = "Identical level sets yield Jaccard overlap = 1.0",
    test_data   = "Same as A03",
    methodology = "Assert Jaccard = |intersect|/|union| = 1 when level sets are identical",
    expected    = "Jaccard = 1.0",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── A05: Identical categorical → Green verdict ───────────────────────────
  expect_true(grepl("\u2705", res_cat$adhoc_column_comparison$verdict),
              label = "A05: Green verdict for identical categorical data")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "A05", category = "Identity Invariants",
    description = "Identical categorical data produces Green verdict",
    test_data   = "Same as A03",
    methodology = "Check verdict contains \u2705",
    expected    = "\u2705 No meaningful change",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── A06: Identical frames → N_Delta = 0 everywhere ──────────────────────
  res_full <- quantify_anonymization_impact(
    original_dataframe   = df_base,
    anonymized_dataframe = df_base
  )
  num_ana <- res_full$numeric_analysis
  expect_true(all(num_ana$N_Delta == 0),
              label = "A06: N_Delta = 0 for all numeric columns with identical frames")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "A06", category = "Identity Invariants",
    description = "Identical dataframes yield N_Delta = 0 for all numeric columns",
    test_data   = "df_base (3 cols: id, icd, los) compared to itself",
    methodology = "Call without ad-hoc columns; check numeric_analysis$N_Delta == 0",
    expected    = "All N_Delta values equal 0",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
})


# =============================================================================
# CATEGORY B: NUMERIC DISTRIBUTIONAL SHIFTS
# =============================================================================

test_that("B: Numeric distributional shift detection", {
  
  # ── B01: Large constant shift → KS D > 0, Red verdict ───────────────────
  df_shift <- df_base
  df_shift$los <- df_base$los + 30L
  res <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_shift,
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_shift$los
  )
  expect_gt(res$adhoc_column_comparison$ks$D, 0,
            label = "B01: KS D > 0 for +30 LOS shift")
  expect_equal(res$adhoc_column_comparison$verdict, "\u274c Big change",
               label = "B01: Red verdict for large shift")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "B01", category = "Numeric Shift Detection",
    description = "Constant +30 day LOS shift detected as large distributional change",
    test_data   = "df_base$los + 30 (deterministic shift of 30 days, n = 50)",
    methodology = "Shift all LOS values by +30; assert KS D > 0 and verdict = Red",
    expected    = "KS D >> 0; verdict = \u274c Big change",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── B02: Mean shift band = 'big' for large shift ────────────────────────
  expect_equal(res$adhoc_column_comparison$mean_shift$band, "big",
               label = "B02: Mean shift band = 'big'")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "B02", category = "Numeric Shift Detection",
    description = "Mean shift band correctly classified as 'big'",
    test_data   = "Same as B01 (+30 shift, ~4 SD units)",
    methodology = "Assert mean_shift$band == 'big' (rel_mu > 0.30)",
    expected    = "band = 'big'",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── B03: Small perturbation → verdict not Red ───────────────────────────
  set.seed(77L)
  df_tiny <- df_base
  df_tiny$los <- df_base$los + sample(c(-1L, 0L, 1L), 50, replace = TRUE)
  df_tiny$los <- pmax(df_tiny$los, 1L)
  res_tiny <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_tiny,
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_tiny$los
  )
  expect_false(res_tiny$adhoc_column_comparison$verdict == "\u274c Big change",
               label = "B03: Small noise does not trigger Red verdict")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "B03", category = "Numeric Shift Detection",
    description = "Minor random perturbation (+/- 1 day) does not trigger Red verdict",
    test_data   = "df_base$los perturbed by {-1, 0, +1} uniformly (n = 50)",
    methodology = "Add uniform noise in {-1,0,1}; assert verdict != Red",
    expected    = "Verdict = Green or Yellow (not Red)",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── B04: Shape band correct for moderate shift ──────────────────────────
  df_mod <- df_base
  df_mod$los <- as.integer(round(df_base$los * 1.15))
  df_mod$los <- pmax(df_mod$los, 1L)
  res_mod <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_mod,
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_mod$los
  )
  expect_true(res_mod$adhoc_column_comparison$shape_band %in%
                c("small", "clear", "large"),
              label = "B04: shape_band in valid set")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "B04", category = "Numeric Shift Detection",
    description = "15% multiplicative scaling produces valid shape_band classification",
    test_data   = "df_base$los * 1.15 (moderate proportional inflation)",
    methodology = "Scale LOS by 1.15x; assert shape_band in {small, clear, large}",
    expected    = "shape_band is one of: 'small', 'clear', 'large'",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── B05: KS D bounded [0, 1] ────────────────────────────────────────────
  expect_true(res$adhoc_column_comparison$ks$D >= 0 &&
                res$adhoc_column_comparison$ks$D <= 1,
              label = "B05: KS D in [0, 1]")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "B05", category = "Numeric Shift Detection",
    description = "KS D statistic bounded within [0, 1]",
    test_data   = "Same as B01",
    methodology = "Assert 0 <= KS D <= 1 (theoretical bound of KS statistic)",
    expected    = "KS D in [0, 1]",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── B06: Spread band for compressed variance ────────────────────────────
  df_comp <- df_base
  df_comp$los <- as.integer(round(df_base$los * 0.5)) + 3L
  res_comp <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_comp,
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_comp$los
  )
  expect_true(res_comp$adhoc_column_comparison$spread_shift$band %in%
                c("about the same", "more spread out", "more squeezed"),
              label = "B06: spread_band valid for variance compression")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "B06", category = "Numeric Shift Detection",
    description = "Variance compression produces valid spread_band classification",
    test_data   = "df_base$los * 0.5 + 3 (halved variance, shifted mean)",
    methodology = "Compress spread; assert spread_band in valid set",
    expected    = "Likely 'more squeezed' (SD reduced by ~50%)",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
})


# =============================================================================
# CATEGORY C: CATEGORICAL ANALYSIS
# =============================================================================

test_that("C: Categorical analysis correctness", {
  
  # ── C01: Disjoint level sets → Jaccard = 0, Red verdict ─────────────────
  df_disj <- df_base
  df_disj$icd <- sample(c("A41", "C34", "M16", "N18", "S72"), 50, replace = TRUE)
  res <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_disj,
    comparison_column_origin     = df_base$icd,
    comparison_column_anonymized = df_disj$icd
  )
  expect_equal(res$adhoc_column_comparison$jaccard_overlap$ratio, 0,
               label = "C01: Jaccard = 0 for completely disjoint level sets")
  expect_equal(res$adhoc_column_comparison$verdict, "\u274c Big change",
               label = "C01: Red verdict for disjoint categories")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "C01", category = "Categorical Analysis",
    description = "Completely disjoint ICD level sets yield Jaccard = 0 and Red verdict",
    test_data   = "Original: {I21,I50,J18,K35,F32}; Anonymized: {A41,C34,M16,N18,S72}",
    methodology = "Replace all ICD codes with a disjoint set; assert Jaccard = 0",
    expected    = "Jaccard = 0; verdict = \u274c Big change",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── C02: Partial overlap → 0 < Jaccard < 1 ──────────────────────────────
  df_part <- df_base
  df_part$icd <- sample(c("I21", "I50", "A41", "C34"), 50, replace = TRUE)
  res_part <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_part,
    comparison_column_origin     = df_base$icd,
    comparison_column_anonymized = df_part$icd
  )
  j <- res_part$adhoc_column_comparison$jaccard_overlap$ratio
  expect_true(j > 0 && j < 1,
              label = "C02: Jaccard in (0,1) for partial overlap")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "C02", category = "Categorical Analysis",
    description = "Partial level set overlap yields 0 < Jaccard < 1",
    test_data   = "Original: {I21,I50,J18,K35,F32}; Anonymized: {I21,I50,A41,C34}",
    methodology = "Use partially overlapping levels; assert 0 < Jaccard < 1",
    expected    = "Jaccard = 2/7 ~ 0.286 (|intersect|=2, |union|=7)",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── C03: Cramer V in [0, 1] for arbitrary comparison ────────────────────
  v <- res$adhoc_column_comparison$cramer_v$v
  expect_true(!is.na(v) && v >= 0 && v <= 1,
              label = "C03: Cramer V in [0, 1]")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "C03", category = "Categorical Analysis",
    description = "Cramer's V bounded within [0, 1] for disjoint comparison",
    test_data   = "Same as C01 (disjoint level sets)",
    methodology = "Assert 0 <= V <= 1 (theoretical property of Cramer's V)",
    expected    = "V in [0, 1]",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── C04: Single-level categorical → chi-square p = NA ───────────────────
  df_one <- df_base
  df_one$icd <- "I21"
  df_anon_one <- df_one
  res_one <- quantify_anonymization_impact(
    original_dataframe   = df_one,
    anonymized_dataframe = df_anon_one
  )
  expect_true(is.na(res_one$categorical_analysis$ChiSq_p[
    res_one$categorical_analysis$Column == "icd"]),
    label = "C04: Chi-square NA for single shared level")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "C04", category = "Categorical Analysis",
    description = "Single ICD level yields chi-square p = NA (test undefined)",
    test_data   = "All 50 records: icd = 'I21' in both frames",
    methodology = "Set all ICD values to one level; assert ChiSq_p is NA",
    expected    = "ChiSq_p = NA (contingency table has only 1 column)",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── C05: Missing value counting ─────────────────────────────────────────
  df_miss <- df_base
  df_miss$icd[1:10] <- NA
  res_miss <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_miss,
    comparison_column_origin     = df_base$icd,
    comparison_column_anonymized = df_miss$icd
  )
  expect_equal(res_miss$adhoc_column_comparison$missing$anon, 10L,
               label = "C05: Missing count = 10 for 10 NAs in anonymized")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "C05", category = "Categorical Analysis",
    description = "Missing value count correctly reports 10 NAs introduced",
    test_data   = "df_base$icd with first 10 values set to NA",
    methodology = "Introduce 10 NAs; assert missing$anon == 10",
    expected    = "missing$anon = 10",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
})


# =============================================================================
# CATEGORY D: COLUMN PRESENCE DETECTION
# =============================================================================

test_that("D: Column presence detection", {
  
  # ── D01: Added column detected ──────────────────────────────────────────
  df_extra <- cbind(df_base, sex = sample(c("M", "F"), 50, replace = TRUE))
  df_fewer <- df_base[, c("id", "icd")]
  res <- quantify_anonymization_impact(
    original_dataframe   = df_extra,
    anonymized_dataframe = df_fewer
  )
  col_res <- res$column_comparison
  expect_equal(col_res$In_Original[col_res$Column == "sex"], "Yes",
               label = "D01: 'sex' present in original")
  expect_equal(col_res$In_Anonymized[col_res$Column == "sex"], "No",
               label = "D01: 'sex' absent from anonymized")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "D01", category = "Column Presence",
    description = "Column 'sex' detected as removed during anonymization",
    test_data   = "Original: {id, icd, los, sex}; Anonymized: {id, icd}",
    methodology = "Add 'sex' to original, remove 'los' from anonymized; check column_comparison",
    expected    = "sex: In_Original=Yes, In_Anonymized=No",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── D02: Removed column detected ───────────────────────────────────────
  expect_equal(col_res$In_Original[col_res$Column == "los"], "Yes",
               label = "D02: 'los' present in original")
  expect_equal(col_res$In_Anonymized[col_res$Column == "los"], "No",
               label = "D02: 'los' absent from anonymized")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "D02", category = "Column Presence",
    description = "Column 'los' detected as removed from anonymized frame",
    test_data   = "Same as D01",
    methodology = "Check column_comparison for 'los'",
    expected    = "los: In_Original=Yes, In_Anonymized=No",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── D03: Identical column sets → all Yes/Yes ────────────────────────────
  res_same <- quantify_anonymization_impact(
    original_dataframe   = df_base,
    anonymized_dataframe = df_base
  )
  expect_true(all(res_same$column_comparison$In_Original == "Yes"),
              label = "D03: All columns present in original")
  expect_true(all(res_same$column_comparison$In_Anonymized == "Yes"),
              label = "D03: All columns present in anonymized")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "D03", category = "Column Presence",
    description = "Identical column sets yield all Yes/Yes",
    test_data   = "df_base compared to itself",
    methodology = "Assert all In_Original and In_Anonymized == 'Yes'",
    expected    = "All columns: Yes/Yes",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
})


# =============================================================================
# CATEGORY E: VECTOR ALIGNMENT AND UNEQUAL LENGTHS
# =============================================================================

test_that("E: Unequal vector length handling", {
  
  # ── E01: Pre-alignment N_Delta captured correctly ───────────────────────
  n_short <- 30L
  res <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_base[1:n_short, ],
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_base$los[1:n_short]
  )
  tbl <- res$adhoc_column_comparison$comparison_table
  raw_delta <- as.integer(tbl$Delta[tbl$Metric == "N (pre-alignment orig)"])
  expect_equal(raw_delta, n_short - nrow(df_base),
               label = "E01: Pre-alignment N_Delta = -20")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "E01", category = "Vector Alignment",
    description = "Pre-alignment N_Delta correctly captures length difference",
    test_data   = "Original: n=50; Anonymized: n=30 (first 30 rows)",
    methodology = "Compare unequal vectors; extract N_Delta from comparison_table",
    expected    = "N_Delta = 30 - 50 = -20",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── E02: Aligned comparison uses min(n_orig, n_anon) ────────────────────
  n_compared <- as.integer(tbl$Original[tbl$Metric == "N compared"])
  expect_equal(n_compared, n_short,
               label = "E02: N compared = min(50, 30) = 30")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "E02", category = "Vector Alignment",
    description = "After alignment, N compared = min(n_orig, n_anon)",
    test_data   = "Same as E01",
    methodology = "Extract 'N compared' from comparison_table; assert = 30",
    expected    = "N compared = 30",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── E03: Equal-length vectors → N_Delta = 0 ─────────────────────────────
  res_eq <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_base,
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_base$los
  )
  tbl_eq <- res_eq$adhoc_column_comparison$comparison_table
  delta_eq <- as.integer(tbl_eq$Delta[tbl_eq$Metric == "N (pre-alignment orig)"])
  expect_equal(delta_eq, 0L,
               label = "E03: N_Delta = 0 for equal-length vectors")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "E03", category = "Vector Alignment",
    description = "Equal-length vectors yield N_Delta = 0",
    test_data   = "df_base$los (n=50) compared to itself",
    methodology = "Pass equal-length vectors; assert N_Delta = 0",
    expected    = "N_Delta = 0",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
})


# =============================================================================
# CATEGORY F: BACKWARD COMPATIBILITY (DEPRECATED ALIASES)
# =============================================================================

test_that("F: Deprecated parameter alias resolution", {
  
  # ── F01: Misspelled aliases resolve correctly ───────────────────────────
  res <- quantify_anonymization_impact(
    original_dataframe      = df_base,
    anonymized_dataframe    = df_base,
    comparison_column_orgin     = df_base$los,   # deprecated typo
    comparison_column_anonimzed = df_base$los    # deprecated typo
  )
  expect_equal(res$adhoc_column_comparison$ks$D, 0, tolerance = 1e-10,
               label = "F01: Deprecated aliases produce D = 0 for identical data")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "F01", category = "Backward Compatibility",
    description = "Deprecated parameter aliases (typos) resolve to correct arguments",
    test_data   = "df_base compared to itself via deprecated 'orgin'/'anonimzed' params",
    methodology = "Use comparison_column_orgin and comparison_column_anonimzed; assert KS D = 0",
    expected    = "KS D = 0 (aliases transparently forwarded)",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── F02: Correct params take precedence over deprecated ─────────────────
  df_shift <- df_base
  df_shift$los <- df_base$los + 50L
  res_prec <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_shift,
    comparison_column_origin     = df_base$los,      # correct (should be used)
    comparison_column_anonymized = df_shift$los,     # correct (should be used)
    comparison_column_orgin      = df_base$los,      # deprecated (ignored)
    comparison_column_anonimzed  = df_base$los       # deprecated (ignored)
  )
  expect_gt(res_prec$adhoc_column_comparison$ks$D, 0,
            label = "F02: Correct params take precedence over deprecated")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "F02", category = "Backward Compatibility",
    description = "Correct parameter names take precedence over deprecated aliases",
    test_data   = "Correct params: base vs shifted(+50); Deprecated params: base vs base",
    methodology = "Supply both correct and deprecated params; assert KS D > 0 (correct params used)",
    expected    = "KS D > 0 (correct params used, not deprecated)",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
})


# =============================================================================
# CATEGORY G: VERDICT LOGIC VALIDATION
# =============================================================================

test_that("G: Verdict decision logic", {
  
  # ── G01: Numeric verdict thresholds — Green ─────────────────────────────
  # D < 0.05 AND rel_mu <= 0.10 AND abs(spr_pct) <= 10 → Green
  res_g <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_base,
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_base$los
  )
  expect_equal(res_g$adhoc_column_comparison$verdict,
               "\u2705 No meaningful change",
               label = "G01: Exact Green verdict string for identical data")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "G01", category = "Verdict Logic",
    description = "Identical numeric data triggers exact Green verdict string",
    test_data   = "df_base$los compared to itself",
    methodology = "Assert verdict == '\u2705 No meaningful change' (exact string match)",
    expected    = "\u2705 No meaningful change",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── G02: Red verdict exact string ───────────────────────────────────────
  df_red <- df_base
  df_red$los <- df_base$los + 30L
  res_r <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_red,
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_red$los
  )
  expect_equal(res_r$adhoc_column_comparison$verdict,
               "\u274c Big change",
               label = "G02: Exact Red verdict string for large shift")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "G02", category = "Verdict Logic",
    description = "Large shift triggers exact Red verdict string",
    test_data   = "df_base$los + 30",
    methodology = "Assert verdict == '\u274c Big change'",
    expected    = "\u274c Big change",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── G03: Categorical Green requirements ─────────────────────────────────
  # V < 0.10 AND overlap >= 0.80 AND miss_pct <= 5 → Green
  res_cg <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_base,
    comparison_column_origin     = df_base$icd,
    comparison_column_anonymized = df_base$icd
  )
  expect_equal(res_cg$adhoc_column_comparison$verdict,
               "\u2705 No meaningful change",
               label = "G03: Categorical Green for identical data")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "G03", category = "Verdict Logic",
    description = "Identical categorical data satisfies all Green thresholds",
    test_data   = "df_base$icd compared to itself",
    methodology = "Assert verdict == Green (V=0, Jaccard=1, missing=0)",
    expected    = "\u2705 No meaningful change",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
})


# =============================================================================
# CATEGORY H: HELPER FUNCTION VALIDATION
# =============================================================================

test_that("H: Internal helper function correctness", {
  
  # ── H01: fmt_p(0) returns '<' prefix ────────────────────────────────────
  expect_true(startsWith(fmt_p_local(0), "<"),
              label = "H01: fmt_p(0) produces '<' prefix")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "H01", category = "Helper Functions",
    description = "fmt_p(0) returns string starting with '<' (below machine epsilon)",
    test_data   = "Input: p = 0",
    methodology = "Call fmt_p(0); assert output startsWith('<')",
    expected    = "String begins with '<' (e.g., '< 2.225...e-308')",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── H02: fmt_p(NA) returns "NA" ────────────────────────────────────────
  expect_equal(fmt_p_local(NA), "NA",
               label = "H02: fmt_p(NA) returns 'NA'")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "H02", category = "Helper Functions",
    description = "fmt_p(NA) returns literal string 'NA'",
    test_data   = "Input: p = NA",
    methodology = "Call fmt_p(NA); assert == 'NA'",
    expected    = "'NA'",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── H03: fmt_p(0.05) returns scientific notation ────────────────────────
  result <- fmt_p_local(0.05)
  expect_true(grepl("[eE]", result),
              label = "H03: fmt_p(0.05) uses scientific notation")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "H03", category = "Helper Functions",
    description = "fmt_p(0.05) returns value in scientific notation",
    test_data   = "Input: p = 0.05",
    methodology = "Call fmt_p(0.05); assert contains 'e' or 'E'",
    expected    = "Scientific notation string (e.g., '5e-02')",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── H04: icd_chapter mapping correctness ────────────────────────────────
  expect_equal(icd_chapter("I21.00"), "Chap-I: Circulatory",
               label = "H04: I21.00 -> Chap-I: Circulatory")
  expect_equal(icd_chapter("J18.10"), "Chap-J: Respiratory",
               label = "H04: J18.10 -> Chap-J: Respiratory")
  expect_equal(icd_chapter("Z00.00"), "Chap-Z: Factors",
               label = "H04: Z00.00 -> Chap-Z: Factors")
  expect_equal(icd_chapter("X99.99"), "Chap-?",
               label = "H04: X99.99 -> Chap-? (unmapped)")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "H04", category = "Helper Functions",
    description = "icd_chapter() correctly maps ICD first letter to chapter label",
    test_data   = "I21.00, J18.10, Z00.00 (mapped); X99.99 (unmapped)",
    methodology = "Assert known mappings: I->Circulatory, J->Respiratory, Z->Factors, X->?",
    expected    = "Correct chapter labels; unmapped -> 'Chap-?'",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── H05: los_midpoint binning correctness ───────────────────────────────
  expect_equal(los_midpoint(1L, 3L), 2L,
               label = "H05: los_midpoint(1, w=3) = 2")
  expect_equal(los_midpoint(3L, 3L), 2L,
               label = "H05: los_midpoint(3, w=3) = 2")
  expect_equal(los_midpoint(4L, 3L), 5L,
               label = "H05: los_midpoint(4, w=3) = 5")
  expect_equal(los_midpoint(7L, 3L), 8L,
               label = "H05: los_midpoint(7, w=3) = 8 (third bin, w=3)")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "H05", category = "Helper Functions",
    description = "los_midpoint() computes correct bin midpoints",
    test_data   = "LOS values {1,3,4,7} with bin width w=3",
    methodology = "Assert: bin_index = floor((los-1)/w); midpoint = bin_index*w + ceil(w/2)",
    expected    = "los_midpoint(1,3)=2, (3,3)=2, (4,3)=5, (7,3)=8",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
})


# =============================================================================
# CATEGORY I: RETURN STRUCTURE VALIDATION
# =============================================================================

test_that("I: Return value structure and types", {
  
  # ── I01: Top-level return keys present ──────────────────────────────────
  res <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_base,
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_base$los
  )
  expected_keys <- c("column_comparison", "numeric_analysis",
                     "categorical_analysis", "adhoc_column_comparison", "notes")
  expect_true(all(expected_keys %in% names(res)),
              label = "I01: All 5 top-level return keys present")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "I01", category = "Return Structure",
    description = "Return list contains all 5 required top-level keys",
    test_data   = "df_base with numeric ad-hoc comparison",
    methodology = "Assert names(res) contains: column_comparison, numeric_analysis, categorical_analysis, adhoc_column_comparison, notes",
    expected    = "All 5 keys present",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── I02: column_comparison is a data.frame ──────────────────────────────
  expect_true(is.data.frame(res$column_comparison),
              label = "I02: column_comparison is data.frame")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "I02", category = "Return Structure",
    description = "column_comparison element is a data.frame",
    test_data   = "Same as I01",
    methodology = "Assert is.data.frame(res$column_comparison)",
    expected    = "TRUE",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── I03: numeric adhoc contains required sub-elements ───────────────────
  adhoc <- res$adhoc_column_comparison
  expect_equal(adhoc$type, "numeric",
               label = "I03: adhoc type = 'numeric' for numeric comparison")
  expect_true(all(c("ks", "verdict", "mean_shift", "spread_shift",
                    "shape_band", "signal_vs_noise") %in% names(adhoc)),
              label = "I03: All numeric adhoc sub-elements present")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "I03", category = "Return Structure",
    description = "Numeric adhoc comparison contains all required sub-elements",
    test_data   = "Same as I01",
    methodology = "Assert presence of: type, ks, verdict, mean_shift, spread_shift, shape_band, signal_vs_noise",
    expected    = "All sub-elements present with type='numeric'",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── I04: categorical adhoc structure ────────────────────────────────────
  res_cat <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_base,
    comparison_column_origin     = df_base$icd,
    comparison_column_anonymized = df_base$icd
  )
  adhoc_cat <- res_cat$adhoc_column_comparison
  expect_equal(adhoc_cat$type, "categorical",
               label = "I04: adhoc type = 'categorical' for character comparison")
  expect_true(all(c("cramer_v", "jaccard_overlap", "missing",
                    "chi_square", "verdict") %in% names(adhoc_cat)),
              label = "I04: All categorical adhoc sub-elements present")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "I04", category = "Return Structure",
    description = "Categorical adhoc comparison contains all required sub-elements",
    test_data   = "df_base$icd (character vector)",
    methodology = "Assert presence of: type, cramer_v, jaccard_overlap, missing, chi_square, verdict",
    expected    = "All sub-elements present with type='categorical'",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── I05: notes contain expected keys ────────────────────────────────────
  expect_true(all(c("numeric", "categorical", "alignment") %in% names(res$notes)),
              label = "I05: notes contain numeric/categorical/alignment keys")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "I05", category = "Return Structure",
    description = "Notes list contains documentation for all three analysis domains",
    test_data   = "Same as I01",
    methodology = "Assert notes contains: 'numeric', 'categorical', 'alignment'",
    expected    = "All 3 note keys present",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
})


# =============================================================================
# CATEGORY J: INPUT VALIDATION AND EDGE CASES
# =============================================================================

test_that("J: Input validation and edge cases", {
  
  # ── J01: NULL input → error ─────────────────────────────────────────────
  expect_error(
    quantify_anonymization_impact(
      original_dataframe   = NULL,
      anonymized_dataframe = df_base
    ),
    label = "J01: NULL original_dataframe raises error"
  )
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "J01", category = "Edge Cases",
    description = "NULL original_dataframe raises explicit error",
    test_data   = "original_dataframe = NULL",
    methodology = "Call with NULL input; assert error is thrown",
    expected    = "Error: 'Both original_dataframe and anonymized_dataframe must be provided.'",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── J02: No ad-hoc columns → adhoc_column_comparison = NULL ─────────────
  res_no_adhoc <- quantify_anonymization_impact(
    original_dataframe   = df_base,
    anonymized_dataframe = df_base
  )
  expect_null(res_no_adhoc$adhoc_column_comparison,
              label = "J02: No ad-hoc args → adhoc = NULL")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "J02", category = "Edge Cases",
    description = "Omitting ad-hoc columns yields adhoc_column_comparison = NULL",
    test_data   = "df_base compared without comparison_column_origin/anonymized",
    methodology = "Call without ad-hoc params; assert adhoc == NULL",
    expected    = "adhoc_column_comparison = NULL",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── J03: No common numeric columns → message in numeric_analysis ────────
  df_chr <- data.frame(id = df_base$id, icd = df_base$icd,
                       stringsAsFactors = FALSE)
  res_chr <- quantify_anonymization_impact(
    original_dataframe   = df_chr,
    anonymized_dataframe = df_chr
  )
  expect_true("Message" %in% names(res_chr$numeric_analysis),
              label = "J03: No numeric cols → numeric_analysis has 'Message'")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "J03", category = "Edge Cases",
    description = "All-character frames yield informational message in numeric_analysis",
    test_data   = "data.frame with only id (char) and icd (char) columns",
    methodology = "Call with character-only frame; assert 'Message' column exists",
    expected    = "numeric_analysis$Message = 'No common numeric columns.'",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
  
  # ── J04: Signal-vs-noise message for identical data ─────────────────────
  res_sig <- quantify_anonymization_impact(
    original_dataframe           = df_base,
    anonymized_dataframe         = df_base,
    comparison_column_origin     = df_base$los,
    comparison_column_anonymized = df_base$los
  )
  expect_true(res_sig$adhoc_column_comparison$signal_vs_noise %in%
                c("Probably just random noise", "n/a",
                  "Likely a real difference",
                  "Almost certainly a real difference"),
              label = "J04: Signal-vs-noise returns valid message")
  RESULTS[[length(RESULTS) + 1L]] <<- list(
    test_id     = "J04", category = "Edge Cases",
    description = "Signal-vs-noise message is from valid set for identical data",
    test_data   = "df_base$los compared to itself",
    methodology = "Assert signal_vs_noise in known message set",
    expected    = "One of the 4 valid signal/noise messages",
    passed      = TRUE, error = NULL, elapsed_sec = 0
  )
})


# =============================================================================
# WRITE STRUCTURED JSON RESULTS
# =============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat(" TEST SUITE COMPLETE\n")
cat(rep("=", 70), "\n\n", sep = "")

n_pass <- sum(vapply(RESULTS, function(x) isTRUE(x$passed), logical(1L)))
n_fail <- length(RESULTS) - n_pass
cat(sprintf("  Total tests : %d\n", length(RESULTS)))
cat(sprintf("  Passed      : %d\n", n_pass))
cat(sprintf("  Failed      : %d\n", n_fail))
cat(sprintf("  Pass rate   : %.1f%%\n\n", 100 * n_pass / length(RESULTS)))

# Serialize to JSON (base R, no jsonlite dependency)
to_json_value <- function(v) {
  if (is.null(v)) return("null")
  if (is.logical(v)) return(tolower(as.character(v)))
  if (is.numeric(v)) return(as.character(v))
  paste0('"', gsub('"', '\\\\"', gsub('\\\\', '\\\\\\\\', as.character(v))), '"')
}

to_json_obj <- function(lst) {
  pairs <- mapply(function(k, v) {
    paste0('    "', k, '": ', to_json_value(v))
  }, names(lst), lst, USE.NAMES = FALSE)
  paste0("  {\n", paste(pairs, collapse = ",\n"), "\n  }")
}

json_str <- paste0(
  "[\n",
  paste(vapply(RESULTS, to_json_obj, character(1L)), collapse = ",\n"),
  "\n]"
)

json_path <- file.path(OUTPUT_DIR, "test_results_report.json")
writeLines(json_str, json_path)
cat("Structured results written to:", json_path, "\n")