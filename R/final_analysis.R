# ============================================================================
# STUDY: Quantifying the Impact of k-Anonymization on Clinical Data Quality
# ============================================================================
#
# Author  : Gaetan Kamdje Wabo
# Affil.  : Dept. of Biomedical Informatics, 
#           Medical Faculty of Mannheim, University of Heidelberg
#           
# Date    : 2026-03-03
# R ver.  : >= 4.1.0 (requires: lme4 for Linear Mixed Models)
#
# PURPOSE
# -------
# This script performs a systematic, multi-level comparison of real-world
# clinical encounter data (encounter_id, main_diagnosis_icd, los_days)
# before and after k-anonymization at three privacy thresholds (k = 5,
# k = 10, k = 15) using the ARX Data Anonymization Tool (v3.9+).
#
# The analysis leverages quantify_anonymization_impact() to compute:
#   - Column presence / absence
#   - Numeric distributional analysis (KS D, sqrt(n_eff)*D, mean/SD shifts)
#   - Categorical analysis (Jaccard overlap, Cramer V, chi-square)
#   - Traffic-light verdicts with handling recommendations
#
#
# DATA LAYOUT
# -----------
# Source : A:\HLZ\Promotionen\Gaetan Kamdje Wabo\study datasets\
#          study files\arx\data source\
# Files  : D0.csv  (original, N approx 720,359 encounters)
#          D5.csv  (k = 5  anonymized)
#          D10.csv (k = 10 anonymized)
#          D15.csv (k = 15 anonymized)
# Schema : encounter_id ; main_diagnosis_icd ; los_days
# Sep    : semicolon (;)
# ============================================================================


# ============================================================================
# 0) CONFIGURATION
# ============================================================================

# -- Data source path --------------------------------------------------------
# Adjust this single path to point to the folder containing D0-D15 CSV files.
# Default assumes the working directory is the repository root.
DATA_DIR <- file.path("data")

# -- Output directory for plots and tables -----------------------------------
OUTPUT_DIR <- file.path("output")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# -- Plot settings for publication -------------------------------------------
PLOT_WIDTH  <- 10    # inches
PLOT_HEIGHT <- 7     # inches

# -- Anonymization levels to compare against D0 (original) -------------------
ANON_LEVELS <- c("D5", "D10", "D15")
K_VALUES    <- c(5, 10, 15)   # corresponding k-anonymity thresholds


# ============================================================================
# 1) SOURCE THE CORE FUNCTION
# ============================================================================
# The function quantify_anonymization_impact() must be in the same directory
# or sourced from its file. Adjust the path below if needed.

source(file.path("R", "quantify_anonymization_impact.R"))
cat("=== quantify_anonymization_impact() loaded successfully ===\n\n")



# 2) IMPORT REAL-WORLD DATA

# Each CSV uses semicolon separation and contains three columns:
#   encounter_id        (integer, identifying attribute)
#   main_diagnosis_icd  (string,  quasi-identifier, ICD-10-GM code)
#   los_days            (numeric, quasi-identifier, length of stay in days)
#
# D0.csv is the original (unanonymized) dataset.
# D5/D10/D15 are the outputs of ARX at increasing k-anonymity thresholds.

cat("--- Loading datasets ---\n")

# Helper: read a study CSV with validation
read_study_csv <- function(filename) {
  filepath <- file.path(DATA_DIR, filename)
  if (!file.exists(filepath)) stop(paste("File not found:", filepath))
  df <- read.csv(filepath, sep = ";", stringsAsFactors = FALSE,
                 header = TRUE, na.strings = c("", "NA", "*"))
  cat(sprintf("  %-10s : %d rows x %d cols\n", filename, nrow(df), ncol(df)))
  return(df)
}

# Import all four datasets
df_D0  <- read_study_csv("D0.csv")     # Original (unanonymized)
df_D5  <- read_study_csv("D5.csv")     # k = 5
df_D10 <- read_study_csv("D10.csv")    # k = 10
df_D15 <- read_study_csv("D15.csv")    # k = 15

# Store in a named list for iteration
datasets <- list(D0 = df_D0, D5 = df_D5, D10 = df_D10, D15 = df_D15)

cat("\n--- Dataset overview ---\n")
for (nm in names(datasets)) {
  df <- datasets[[nm]]
  cat(sprintf("  %s : N = %s | ICD distinct = %d | los_days range = [%s, %s]\n",
              nm,
              formatC(nrow(df), format = "d", big.mark = ","),
              length(unique(na.omit(df$main_diagnosis_icd))),
              ifelse(all(is.na(df$los_days)), "NA",
                     as.character(min(df$los_days, na.rm = TRUE))),
              ifelse(all(is.na(df$los_days)), "NA",
                     as.character(max(df$los_days, na.rm = TRUE)))))
}


# ============================================================================
# 3) RUN quantify_anonymization_impact() FOR EACH ANONYMIZATION LEVEL
# ============================================================================
# We compare D0 (original) against each of D5, D10, D15, running both:
#   (a) Full-dataset comparison (column presence, numeric, categorical)
#   (b) Ad-hoc deep-dive on los_days (numeric QI)
#   (c) Ad-hoc deep-dive on main_diagnosis_icd (categorical QI)

cat("\n=== Running anonymization impact analysis ===\n")

results_full <- list()   # full comparison (no ad-hoc)
results_los  <- list()   # ad-hoc on los_days
results_icd  <- list()   # ad-hoc on main_diagnosis_icd

for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]
  k   <- K_VALUES[i]
  df_anon <- datasets[[lvl]]
  
  cat(sprintf("\n--- D0 vs. %s (k = %d) ---\n", lvl, k))
  
  # (a) Full comparison without ad-hoc columns
  results_full[[lvl]] <- quantify_anonymization_impact(
    original_dataframe  = df_D0,
    anonymized_dataframe = df_anon
  )
  
  # (b) Ad-hoc: los_days (numeric quasi-identifier)
  results_los[[lvl]] <- quantify_anonymization_impact(
    original_dataframe           = df_D0,
    anonymized_dataframe         = df_anon,
    comparison_column_origin     = df_D0$los_days,
    comparison_column_anonymized = df_anon$los_days
  )
  
  # (c) Ad-hoc: main_diagnosis_icd (categorical quasi-identifier)
  results_icd[[lvl]] <- quantify_anonymization_impact(
    original_dataframe           = df_D0,
    anonymized_dataframe         = df_anon,
    comparison_column_origin     = df_D0$main_diagnosis_icd,
    comparison_column_anonymized = df_anon$main_diagnosis_icd
  )
  
  cat(sprintf("  los_days verdict  : %s\n",
              results_los[[lvl]]$adhoc_column_comparison$handling_recommendation$verdict))
  cat(sprintf("  ICD verdict       : %s\n",
              results_icd[[lvl]]$adhoc_column_comparison$handling_recommendation$verdict))
}


# ============================================================================
# 4) EXTRACT KEY METRICS INTO A SUMMARY TABLE
# ============================================================================
# Build a concise data.frame across anonymization levels.
# This table is the backbone of the scholarly plots and the final report.

cat("\n=== Building summary metrics table ===\n")

summary_metrics <- data.frame(
  k_level            = K_VALUES,
  label              = ANON_LEVELS,
  N_original         = nrow(df_D0),
  N_anonymized       = sapply(ANON_LEVELS, function(l) nrow(datasets[[l]])),
  records_lost       = sapply(ANON_LEVELS, function(l) nrow(df_D0) - nrow(datasets[[l]])),
  records_lost_pct   = sapply(ANON_LEVELS, function(l)
    100 * (nrow(df_D0) - nrow(datasets[[l]])) / nrow(df_D0)),
  los_KS_D           = sapply(ANON_LEVELS, function(l)
    as.numeric(results_los[[l]]$adhoc_column_comparison$meta$ks_D)),
  los_mean_shift     = sapply(ANON_LEVELS, function(l)
    as.numeric(results_los[[l]]$adhoc_column_comparison$meta$mean_anon) -
      as.numeric(results_los[[l]]$adhoc_column_comparison$meta$mean_orig)),
  los_sd_shift       = sapply(ANON_LEVELS, function(l)
    as.numeric(results_los[[l]]$adhoc_column_comparison$meta$sd_anon) -
      as.numeric(results_los[[l]]$adhoc_column_comparison$meta$sd_orig)),
  los_KS_z           = sapply(ANON_LEVELS, function(l)
    as.numeric(results_los[[l]]$adhoc_column_comparison$meta$ks_z)),
  los_verdict        = sapply(ANON_LEVELS, function(l)
    results_los[[l]]$adhoc_column_comparison$handling_recommendation$verdict),
  icd_jaccard        = sapply(ANON_LEVELS, function(l)
    as.numeric(results_icd[[l]]$adhoc_column_comparison$meta$overlap)),
  icd_cramer_v       = sapply(ANON_LEVELS, function(l)
    as.numeric(results_icd[[l]]$adhoc_column_comparison$meta$cramer_v)),
  icd_verdict        = sapply(ANON_LEVELS, function(l)
    results_icd[[l]]$adhoc_column_comparison$handling_recommendation$verdict),
  stringsAsFactors = FALSE, row.names = NULL
)

print(summary_metrics[, c("label", "k_level", "N_anonymized", "records_lost_pct",
                          "los_KS_D", "los_mean_shift",
                          "icd_jaccard", "icd_cramer_v")])


# ============================================================================
# 5) SCHOLARLY PLOTS
# ============================================================================
# All plots use base R graphics for zero-dependency reproducibility.
# Color palette: blue (#2E75B6) = mild, amber (#ED7D31) = moderate,
# red (#C00000) = strong anonymization.

cat("\n=== Generating publication-quality plots ===\n")

col_k5  <- "#2E75B6"
col_k10 <- "#ED7D31"
col_k15 <- "#C00000"
col_orig <- "#404040"
COLS <- c(col_k5, col_k10, col_k15)


# --------------------------------------------------------------------------
# PLOT 1: Record Suppression by k-Level
# --------------------------------------------------------------------------
# Shows how many records are lost (suppressed) at each anonymization level.
# This is a fundamental privacy-utility trade-off indicator.

pdf(file.path(OUTPUT_DIR, "Fig1_record_suppression.pdf"),
    width = PLOT_WIDTH, height = PLOT_HEIGHT)
par(mar = c(5, 6, 4, 2), family = "serif")

bp <- barplot(
  summary_metrics$records_lost_pct,
  names.arg = paste0("k = ", summary_metrics$k_level),
  col = COLS, border = NA,
  ylim = c(0, max(summary_metrics$records_lost_pct) * 1.25),
  ylab = "Records Suppressed (%)",
  xlab = "k-Anonymity Threshold",
  main = "Record Suppression Rate by Anonymization Level",
  cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1, cex.names = 1.2, las = 1
)
text(bp, summary_metrics$records_lost_pct,
     labels = sprintf("%.2f%%", summary_metrics$records_lost_pct),
     pos = 3, cex = 1.1, font = 2)
text(bp, summary_metrics$records_lost_pct / 2,
     labels = formatC(summary_metrics$records_lost, format = "d", big.mark = ","),
     cex = 0.9, col = "white", font = 2)
mtext(paste0("Original N = ", formatC(nrow(df_D0), format = "d", big.mark = ",")),
      side = 3, line = 0.3, cex = 0.95, font = 3)
dev.off()
cat("  Fig1_record_suppression.pdf\n")


# --------------------------------------------------------------------------
# PLOT 2: KS D Effect Size Across k-Levels (los_days)
# --------------------------------------------------------------------------
# KS D = maximum distance between empirical CDFs of original vs anonymized.
# D < 0.05 = negligible; D < 0.10 = small; D >= 0.10 = substantial.

pdf(file.path(OUTPUT_DIR, "Fig2_KS_D_los_days.pdf"),
    width = PLOT_WIDTH, height = PLOT_HEIGHT)
par(mar = c(5, 6, 4, 2), family = "serif")

plot(summary_metrics$k_level, summary_metrics$los_KS_D,
     type = "b", pch = 19, cex = 2.5, lwd = 3, col = col_orig,
     xlab = "k-Anonymity Threshold",
     ylab = "Kolmogorov-Smirnov D (Effect Size)",
     main = "Distributional Divergence of los_days\nOriginal vs. Anonymized",
     xlim = c(3, 17), ylim = c(0, max(summary_metrics$los_KS_D, 0.15) * 1.3),
     xaxt = "n", las = 1,
     cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)
axis(1, at = K_VALUES, labels = paste0("k = ", K_VALUES), cex.axis = 1.1)
points(summary_metrics$k_level, summary_metrics$los_KS_D,
       pch = 19, cex = 2.5, col = COLS)
abline(h = 0.05, lty = 2, col = "forestgreen", lwd = 1.5)
abline(h = 0.10, lty = 2, col = "darkorange",  lwd = 1.5)
text(16.5, 0.05, "D = 0.05\n(negligible)", cex = 0.85, col = "forestgreen", pos = 3)
text(16.5, 0.10, "D = 0.10\n(substantial)", cex = 0.85, col = "darkorange",  pos = 3)
text(summary_metrics$k_level, summary_metrics$los_KS_D,
     labels = sprintf("D = %.4f", summary_metrics$los_KS_D),
     pos = 3, cex = 1.0, font = 2, offset = 1.2)
dev.off()
cat("  Fig2_KS_D_los_days.pdf\n")


# --------------------------------------------------------------------------
# PLOT 3: ICD Code Jaccard Overlap and Cramer V (Dual-Axis)
# --------------------------------------------------------------------------
# Jaccard = label-space preservation (what fraction of ICD codes survive).
# Cramer V = distributional shift magnitude. Together they capture both
# structural and distributional aspects of categorical data quality loss.

pdf(file.path(OUTPUT_DIR, "Fig3_ICD_quality_metrics.pdf"),
    width = PLOT_WIDTH, height = PLOT_HEIGHT)
par(mar = c(5, 6, 4, 6), family = "serif")

plot(summary_metrics$k_level, summary_metrics$icd_jaccard,
     type = "b", pch = 17, cex = 2.5, lwd = 3, col = "#2E75B6",
     xlab = "k-Anonymity Threshold", ylab = "Jaccard Overlap Ratio",
     main = "ICD Code Quality: Label Preservation vs. Distributional Shift",
     xlim = c(3, 17), ylim = c(0, 1), xaxt = "n", las = 1,
     cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)
axis(1, at = K_VALUES, labels = paste0("k = ", K_VALUES), cex.axis = 1.1)
text(summary_metrics$k_level, summary_metrics$icd_jaccard,
     labels = sprintf("J = %.3f", summary_metrics$icd_jaccard),
     pos = 1, cex = 0.95, font = 2, col = "#2E75B6", offset = 1.2)

par(new = TRUE)
plot(summary_metrics$k_level, summary_metrics$icd_cramer_v,
     type = "b", pch = 15, cex = 2.5, lwd = 3, col = "#C00000",
     axes = FALSE, xlab = "", ylab = "", xlim = c(3, 17), ylim = c(0, 1))
axis(4, las = 1, col = "#C00000", col.axis = "#C00000", cex.axis = 1.1)
mtext("Cramer's V (Effect Size)", side = 4, line = 3.5, cex = 1.2, col = "#C00000")
text(summary_metrics$k_level, summary_metrics$icd_cramer_v,
     labels = sprintf("V = %.3f", summary_metrics$icd_cramer_v),
     pos = 3, cex = 0.95, font = 2, col = "#C00000", offset = 1.2)
abline(h = 0.10, lty = 3, col = "gray50")
abline(h = 0.30, lty = 3, col = "gray50")
text(3.5, 0.10, "V=0.10 (small)",  cex = 0.8, col = "gray50", pos = 3)
text(3.5, 0.30, "V=0.30 (medium)", cex = 0.8, col = "gray50", pos = 3)
legend("topright",
       legend = c("Jaccard Overlap (left axis)", "Cramer's V (right axis)"),
       pch = c(17, 15), col = c("#2E75B6", "#C00000"),
       lwd = 3, cex = 1.0, bg = "white", box.lwd = 0.5)
dev.off()
cat("  Fig3_ICD_quality_metrics.pdf\n")


# --------------------------------------------------------------------------
# PLOT 4: Empirical CDF Overlay for los_days
# --------------------------------------------------------------------------
# Visual equivalent of the KS test: max vertical distance between curves = D.

pdf(file.path(OUTPUT_DIR, "Fig4_ECDF_los_days.pdf"),
    width = PLOT_WIDTH, height = PLOT_HEIGHT)
par(mar = c(5, 6, 4, 2), family = "serif")

ecdf_orig <- ecdf(df_D0$los_days)
x_max <- quantile(df_D0$los_days, 0.99, na.rm = TRUE)
x_seq <- seq(0, x_max, length.out = 1000)

plot(x_seq, ecdf_orig(x_seq), type = "l", lwd = 3, col = col_orig,
     xlab = "Length of Stay (days)", ylab = expression(hat(F)(x)),
     main = "Empirical CDF: los_days\nOriginal vs. k-Anonymized Datasets",
     las = 1, cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)
for (i in seq_along(ANON_LEVELS)) {
  ecdf_anon <- ecdf(datasets[[ANON_LEVELS[i]]]$los_days)
  lines(x_seq, ecdf_anon(x_seq), lwd = 2.5, col = COLS[i], lty = i + 1)
}
legend("bottomright",
       legend = c("Original (D0)", paste0("k=", K_VALUES, " (", ANON_LEVELS, ")")),
       col = c(col_orig, COLS), lwd = c(3, rep(2.5, 3)),
       lty = c(1, 2, 3, 4), cex = 1.0, bg = "white", box.lwd = 0.5)
dev.off()
cat("  Fig4_ECDF_los_days.pdf\n")


# --------------------------------------------------------------------------
# PLOT 5: 2x2 Dashboard Panel
# --------------------------------------------------------------------------
# Single-page overview: (A) record loss, (B) KS D, (C) Jaccard, (D) Cramer V.

pdf(file.path(OUTPUT_DIR, "Fig5_dashboard_panel.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2), family = "serif")

barplot(summary_metrics$records_lost_pct,
        names.arg = paste0("k=", K_VALUES), col = COLS, border = NA,
        ylab = "Suppressed (%)", main = "(A) Record Suppression",
        cex.main = 1.4, cex.lab = 1.2, las = 1)

barplot(summary_metrics$los_KS_D,
        names.arg = paste0("k=", K_VALUES), col = COLS, border = NA,
        ylab = "KS D", main = "(B) Distributional Shift: los_days",
        cex.main = 1.4, cex.lab = 1.2, las = 1)
abline(h = 0.05, lty = 2, col = "forestgreen")
abline(h = 0.10, lty = 2, col = "darkorange")

barplot(summary_metrics$icd_jaccard,
        names.arg = paste0("k=", K_VALUES), col = COLS, border = NA,
        ylab = "Jaccard", main = "(C) ICD Label-Space Preservation",
        ylim = c(0, 1), cex.main = 1.4, cex.lab = 1.2, las = 1)
abline(h = 0.80, lty = 2, col = "forestgreen")
abline(h = 0.60, lty = 2, col = "darkorange")

barplot(summary_metrics$icd_cramer_v,
        names.arg = paste0("k=", K_VALUES), col = COLS, border = NA,
        ylab = "Cramer's V", main = "(D) ICD Distributional Distortion",
        ylim = c(0, max(summary_metrics$icd_cramer_v, 0.35) * 1.3),
        cex.main = 1.4, cex.lab = 1.2, las = 1)
abline(h = 0.10, lty = 2, col = "forestgreen")
abline(h = 0.30, lty = 2, col = "darkorange")

dev.off()
cat("  Fig5_dashboard_panel.pdf\n")


# --------------------------------------------------------------------------
# PLOT 6: Kernel Density of los_days
# --------------------------------------------------------------------------
# Reveals shape changes that summary statistics miss: mode shifts, tail
# compression, bimodality emergence/collapse.

pdf(file.path(OUTPUT_DIR, "Fig6_density_los_days.pdf"),
    width = PLOT_WIDTH, height = PLOT_HEIGHT)
par(mar = c(5, 6, 4, 2), family = "serif")

x_max_d <- quantile(df_D0$los_days, 0.99, na.rm = TRUE)
d_orig  <- density(df_D0$los_days[df_D0$los_days <= x_max_d], na.rm = TRUE)

y_max <- max(d_orig$y)
for (lvl in ANON_LEVELS) {
  vals <- datasets[[lvl]]$los_days
  vals <- vals[!is.na(vals) & vals <= x_max_d]
  if (length(vals) > 1) y_max <- max(y_max, max(density(vals)$y))
}

plot(d_orig, lwd = 3, col = col_orig,
     xlim = c(0, x_max_d), ylim = c(0, y_max * 1.1),
     xlab = "Length of Stay (days)", ylab = "Density",
     main = "Kernel Density: los_days\nOriginal vs. k-Anonymized Datasets",
     cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1, las = 1)
for (i in seq_along(ANON_LEVELS)) {
  vals <- datasets[[ANON_LEVELS[i]]]$los_days
  vals <- vals[!is.na(vals) & vals <= x_max_d]
  if (length(vals) > 1) lines(density(vals), lwd = 2.5, col = COLS[i], lty = i + 1)
}
legend("topright",
       legend = c("Original (D0)", paste0("k=", K_VALUES, " (", ANON_LEVELS, ")")),
       col = c(col_orig, COLS), lwd = c(3, rep(2.5, 3)),
       lty = c(1, 2, 3, 4), cex = 1.0, bg = "white", box.lwd = 0.5)
dev.off()
cat("  Fig6_density_los_days.pdf\n")


# --------------------------------------------------------------------------
# PLOT 7: ICD Code Rank-Frequency (Zipf Plot)
# --------------------------------------------------------------------------
# Clinical diagnosis codes follow a Zipf-like distribution. Anonymization
# via generalization alters the rank-frequency curve, collapsing the long
# tail of rare codes.

pdf(file.path(OUTPUT_DIR, "Fig7_ICD_zipf_rank.pdf"),
    width = PLOT_WIDTH, height = PLOT_HEIGHT)
par(mar = c(5, 6, 4, 2), family = "serif")

freq_orig <- sort(table(df_D0$main_diagnosis_icd), decreasing = TRUE)
n_codes   <- min(length(freq_orig), 100)

plot(seq_len(n_codes), as.numeric(freq_orig[seq_len(n_codes)]),
     type = "l", lwd = 3, col = col_orig, log = "y",
     xlab = "ICD Code Rank (most frequent first)",
     ylab = "Frequency (log scale)",
     main = "ICD Code Rank-Frequency (Zipf Plot)\nTop 100 Codes",
     cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1, las = 1)
for (i in seq_along(ANON_LEVELS)) {
  freq_anon <- sort(table(datasets[[ANON_LEVELS[i]]]$main_diagnosis_icd),
                    decreasing = TRUE)
  n_a <- min(length(freq_anon), 100)
  lines(seq_len(n_a), as.numeric(freq_anon[seq_len(n_a)]),
        lwd = 2.5, col = COLS[i], lty = i + 1)
}
legend("topright",
       legend = c("Original (D0)", paste0("k=", K_VALUES, " (", ANON_LEVELS, ")")),
       col = c(col_orig, COLS), lwd = c(3, rep(2.5, 3)),
       lty = c(1, 2, 3, 4), cex = 1.0, bg = "white", box.lwd = 0.5)
dev.off()
cat("  Fig7_ICD_zipf_rank.pdf\n")


# --------------------------------------------------------------------------
# PLOT 8: Mean and SD Shifts of los_days (Lollipop Chart)
# --------------------------------------------------------------------------
# Visualizes how much central tendency and spread drift after anonymization.

pdf(file.path(OUTPUT_DIR, "Fig8_los_mean_sd_shifts.pdf"),
    width = PLOT_WIDTH, height = PLOT_HEIGHT)
par(mar = c(5, 8, 4, 2), family = "serif")

y_pos <- c(1, 2, 3, 5, 6, 7)
vals  <- c(summary_metrics$los_mean_shift, summary_metrics$los_sd_shift)
cols  <- rep(COLS, 2)
labs  <- c(paste0("Mean shift (k=", K_VALUES, ")"),
           paste0("SD shift (k=", K_VALUES, ")"))

plot(vals, y_pos, pch = 19, cex = 2.5, col = cols,
     xlim = range(c(vals, 0)) * c(1.4, 1.4),
     yaxt = "n", xlab = "Shift (days)", ylab = "",
     main = "los_days: Mean and SD Shifts After Anonymization",
     cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1, las = 1)
segments(0, y_pos, vals, y_pos, lwd = 2.5, col = cols)
abline(v = 0, lty = 1, col = "gray30", lwd = 1.5)
axis(2, at = y_pos, labels = labs, las = 1, cex.axis = 0.95)
text(vals, y_pos, labels = sprintf("%+.3f", vals),
     pos = ifelse(vals >= 0, 4, 2), cex = 0.95, font = 2)
abline(h = 4, lty = 3, col = "gray70")
dev.off()
cat("  Fig8_los_mean_sd_shifts.pdf\n")


# ============================================================================
# 6) EXPORT SUMMARY TABLE AS CSV
# ============================================================================

write.csv(summary_metrics,
          file = file.path(OUTPUT_DIR, "Table1_summary_metrics.csv"),
          row.names = FALSE)
cat("\n  Table1_summary_metrics.csv saved\n")


# ============================================================================
# 7) CONSOLE REPORT: TRAFFIC-LIGHT VERDICTS
# ============================================================================

cat("\n")
cat("====================================================================\n")
cat("  TRAFFIC-LIGHT VERDICTS: Anonymization Impact Assessment\n")
cat("====================================================================\n\n")

for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]
  k   <- K_VALUES[i]
  
  cat(sprintf("  k = %d (%s)\n", k, lvl))
  cat(sprintf("    Records suppressed : %s (%.2f%%)\n",
              formatC(summary_metrics$records_lost[i], format="d", big.mark=","),
              summary_metrics$records_lost_pct[i]))
  cat(sprintf("    los_days           : %s\n",
              results_los[[lvl]]$adhoc_column_comparison$handling_recommendation$verdict))
  cat(sprintf("      KS D = %.4f | z = %.2f | mean shift = %+.3f | SD shift = %+.3f\n",
              summary_metrics$los_KS_D[i], summary_metrics$los_KS_z[i],
              summary_metrics$los_mean_shift[i], summary_metrics$los_sd_shift[i]))
  cat(sprintf("    main_diagnosis_icd : %s\n",
              results_icd[[lvl]]$adhoc_column_comparison$handling_recommendation$verdict))
  cat(sprintf("      Jaccard = %.3f | Cramer V = %.3f\n",
              summary_metrics$icd_jaccard[i], summary_metrics$icd_cramer_v[i]))
  cat("\n")
}

cat("====================================================================\n")
cat("  Decision Support:\n")
cat("    Green  : Safe for most analyses\n")
cat("    Yellow : Document shifts, run sensitivity analyses\n")
cat("    Red    : Avoid direct comparability, redesign anonymization\n")
cat("====================================================================\n")


# ============================================================================
# 8) SAVE FULL RESULTS AS .RDS FOR DOWNSTREAM USE
# ============================================================================

study_results <- list(
  summary_metrics = summary_metrics,
  results_full    = results_full,
  results_los     = results_los,
  results_icd     = results_icd,
  datasets        = datasets,
  metadata = list(
    date        = Sys.time(),
    R_version   = R.version.string,
    data_dir    = DATA_DIR,
    description = paste(
      "Multi-level k-anonymization impact study.",
      "Original (D0) vs. D5 (k=5), D10 (k=10), D15 (k=15).",
      "Clinical encounter data: encounter_id, main_diagnosis_icd, los_days.",
      "Anonymized using ARX with generalization + suppression."
    )
  )
)

saveRDS(study_results,
        file = file.path(OUTPUT_DIR, "study_results_full.rds"))
cat("\n  study_results_full.rds saved\n")

cat("\n=== Study complete. All outputs in: ===\n")
cat(paste0("  ", OUTPUT_DIR, "\n"))
cat("\nGenerated files:\n")
cat("  Fig1_record_suppression.pdf    - Suppression rates per k-level\n")
cat("  Fig2_KS_D_los_days.pdf         - KS D effect size trajectory\n")
cat("  Fig3_ICD_quality_metrics.pdf   - Jaccard + Cramer V dual-axis\n")
cat("  Fig4_ECDF_los_days.pdf         - Empirical CDF overlay\n")
cat("  Fig5_dashboard_panel.pdf       - 2x2 summary dashboard\n")
cat("  Fig6_density_los_days.pdf      - Kernel density overlay\n")
cat("  Fig7_ICD_zipf_rank.pdf         - Rank-frequency Zipf plot\n")
cat("  Fig8_los_mean_sd_shifts.pdf    - Mean/SD shift lollipop chart\n")
cat("  Table1_summary_metrics.csv     - Machine-readable summary table\n")
cat("  study_results_full.rds         - Full R object for reuse\n")

# ============================================================================
# 9) LINEAR MIXED MODEL: RANDOM INTERCEPTS FOR ICD-3 -> LOS
# ============================================================================
# Sections 2.5.3 (Model Specification) and 2.5.4 (Reproducibility Assessment)
#
# Rationale for LMM over GLM with fixed ICD effects:
#   ARX anonymized D0 -> D5/D10/D15 using cell-level suppression ONLY.
#   No ICD hierarchy / generalization was applied.  Individual values of
#   main_diagnosis_icd and/or los_days were replaced with NA ("*") to
#   satisfy k-anonymity; entire records that could not satisfy k were removed.
#
#   Consequence: main_diagnosis_icd (truncated to ICD-3) is a high-cardinality
#   categorical variable with potentially hundreds of levels.  Fitting each
#   level as a fixed dummy coefficient (as in NegBin GLM) is:
#     (a) computationally prohibitive for > ~300 categories,
#     (b) methodologically questionable (overfitting, unstable estimates
#         for rare categories, no borrowing of information across groups).
#
#   A Linear Mixed Model solves both problems by treating ICD-3 category
#   effects as random intercepts drawn from a common normal distribution:
#     - Rare categories are automatically shrunk toward the grand mean
#     - Variance components (sigma^2_u, sigma^2_e) are estimated efficiently
#     - The ICC directly quantifies ICD's explanatory share of LOS variance
#
# Model specification:
#   log(los_days + 1) ~ 1 + (1 | icd3)
#
#   Response:     log(LOS + 1)  -- log-transformation handles right-skew of
#                 LOS; +1 offset accommodates same-day cases (LOS = 0)
#   Fixed effect: Intercept (grand mean of log-transformed LOS)
#   Random effect: Per-ICD-3 intercept (deviation from grand mean)
#
#   Estimation: REML (Restricted Maximum Likelihood) -- preferred for
#               unbiased variance component estimation when fixed effects
#               are identical across models (Patterson & Thompson, 1971).
#
# Complete-case analysis: rows missing either ICD or LOS are excluded.
# NA counts are documented transparently per dataset (Section 9.1).
#
# Dependency: lme4 (standard R package for mixed-effects models)
# ============================================================================

cat("\n")
cat("####################################################################\n")
cat("#  SECTION 2.5.3 / 2.5.4 -- LMM Regression & Reproducibility     #\n")
cat("####################################################################\n\n")

library(lme4)


# ============================================================================
# 9.1) DOCUMENT ARX-INTRODUCED NAs
# ============================================================================
# Before any transformation, count exactly how many NAs ARX introduced
# per column per dataset.  This is essential context for interpreting
# the regression results: the modeling sample is a strict subset of the
# anonymized file, and that subset shrinks as k increases.

cat("--- 9.1  ARX-introduced NA documentation ---\n\n")

na_documentation <- data.frame(
  dataset          = character(0),
  k_level          = integer(0),
  N_rows           = integer(0),
  NA_icd           = integer(0),
  NA_icd_pct       = numeric(0),
  NA_los           = integer(0),
  NA_los_pct       = numeric(0),
  NA_either        = integer(0),
  NA_either_pct    = numeric(0),
  N_complete       = integer(0),
  N_complete_pct   = numeric(0),
  stringsAsFactors = FALSE
)

all_k <- c(0, K_VALUES)
all_names <- c("D0", ANON_LEVELS)

for (idx in seq_along(all_names)) {
  nm  <- all_names[idx]
  k   <- all_k[idx]
  df  <- datasets[[nm]]
  n   <- nrow(df)
  
  # Detect NAs in ICD: true NA, empty string, whitespace, or ARX suppression "*"
  icd_vals    <- trimws(as.character(df$main_diagnosis_icd))
  is_na_icd   <- is.na(icd_vals) | icd_vals == "" | icd_vals == "*"
  na_icd      <- sum(is_na_icd)
  
  # Detect NAs in los_days: true NA, or non-numeric after coercion
  los_numeric <- suppressWarnings(as.numeric(df$los_days))
  is_na_los   <- is.na(los_numeric)
  na_los      <- sum(is_na_los)
  
  # Either column missing
  is_na_either <- is_na_icd | is_na_los
  na_either    <- sum(is_na_either)
  n_complete   <- n - na_either
  
  na_documentation <- rbind(na_documentation, data.frame(
    dataset        = nm,
    k_level        = k,
    N_rows         = n,
    NA_icd         = na_icd,
    NA_icd_pct     = 100 * na_icd / n,
    NA_los         = na_los,
    NA_los_pct     = 100 * na_los / n,
    NA_either      = na_either,
    NA_either_pct  = 100 * na_either / n,
    N_complete     = n_complete,
    N_complete_pct = 100 * n_complete / n,
    stringsAsFactors = FALSE
  ))
}

# Print the NA documentation table
cat("  NA counts introduced by ARX suppression:\n\n")
print(na_documentation, row.names = FALSE)
cat("\n")

# Save as CSV
write.csv(na_documentation,
          file = file.path(OUTPUT_DIR, "Table2_NA_documentation.csv"),
          row.names = FALSE)
cat("  Table2_NA_documentation.csv saved\n")


# ============================================================================
# 9.2) PREPARE MODELING DATA
# ============================================================================
# Steps per dataset:
#   1. Coerce los_days to numeric (ARX may export as character after suppression)
#   2. Mark suppressed ICD values as NA
#   3. Extract ICD-3 = first 3 characters of ICD-10-GM code
#   4. Keep only complete cases (both los_days and icd3 non-missing, los >= 0)
#   5. Drop ICD-3 categories with < MIN_FREQ observations
#      (LMM handles sparse groups via shrinkage, so threshold is lower than GLM)
#   6. Compute log-transformed response: log(los_days + 1)
#   7. Set the reference factor ordering (most frequent ICD-3 in D0 first)
#
# Note on MIN_FREQ: LMM naturally shrinks random effects for sparse groups
# toward zero, so stability does not require the high threshold (30+) that
# fixed-effect GLMs need.  MIN_FREQ = 5 ensures each group contributes
# meaningfully to the between-group variance estimate while preserving
# more categories than the GLM approach.

cat("\n--- 9.2  Preparing modeling data ---\n")

MIN_FREQ <- 5   # minimum obs per ICD-3 category (LMM shrinkage allows lower)

prepare_model_data <- function(df, label) {
  
  # 1. Coerce los_days to numeric
  df$los_days <- suppressWarnings(as.numeric(df$los_days))
  
  # 2. Clean ICD: mark suppressed values as NA
  icd_raw   <- trimws(as.character(df$main_diagnosis_icd))
  icd_clean <- ifelse(is.na(icd_raw) | icd_raw == "" | icd_raw == "*",
                      NA_character_, icd_raw)
  
  # 3. Extract ICD-3 (first 3 characters, e.g. "I25.1" -> "I25")
  df$icd3 <- substr(icd_clean, 1, 3)
  
  # Validate: ICD-3 must be exactly 3 characters with letter + 2 digits
  valid_icd3 <- grepl("^[A-Z][0-9]{2}$", df$icd3)
  df$icd3[!valid_icd3] <- NA_character_
  
  # 4. Complete cases (both icd3 and los_days present, los >= 0)
  complete <- !is.na(df$icd3) & !is.na(df$los_days) & df$los_days >= 0
  df_cc    <- df[complete, c("icd3", "los_days")]
  
  # 5. Frequency filter
  freq       <- table(df_cc$icd3)
  keep_cats  <- names(freq[freq >= MIN_FREQ])
  df_model   <- df_cc[df_cc$icd3 %in% keep_cats, ]
  df_model$icd3 <- factor(df_model$icd3)
  
  # 6. Log-transformed response
  df_model$log_los <- log(df_model$los_days + 1)
  
  cat(sprintf("  %s: %d rows -> %d complete -> %d after freq filter (%d ICD-3 cats)\n",
              label, nrow(df), nrow(df_cc), nrow(df_model), length(keep_cats)))
  
  list(data = df_model, n_complete = nrow(df_cc),
       n_model = nrow(df_model), n_cats = length(keep_cats))
}

model_data <- list()
for (nm in all_names) {
  model_data[[nm]] <- prepare_model_data(datasets[[nm]], nm)
}

# Reference category: most frequent ICD-3 in D0
ref_cat <- names(which.max(table(model_data[["D0"]]$data$icd3)))
cat(sprintf("\n  Reference ICD-3 category (most frequent in D0): %s\n", ref_cat))

for (nm in all_names) {
  lvls <- levels(model_data[[nm]]$data$icd3)
  if (ref_cat %in% lvls) {
    model_data[[nm]]$data$icd3 <- relevel(model_data[[nm]]$data$icd3, ref = ref_cat)
  }
}


# ============================================================================
# 9.3) FIT LINEAR MIXED MODELS
# ============================================================================
# Model: log(los_days + 1) ~ 1 + (1 | icd3)
#
# lme4::lmer() estimates:
#   - Fixed intercept (grand mean of log-LOS across all ICD-3 groups)
#   - Random intercept variance sigma^2_u (between-group heterogeneity)
#   - Residual variance sigma^2_e (within-group variation)
#   - BLUPs (Best Linear Unbiased Predictions) per ICD-3 category
#
# REML = TRUE for unbiased variance component estimation.
# We also fit with ML (REML = FALSE) to obtain comparable AIC/BIC values.

cat("\n--- 9.3  Fitting Linear Mixed Models ---\n")

lmm_models_reml <- list()   # primary: REML for variance components
lmm_models_ml   <- list()   # secondary: ML for AIC/BIC comparison

for (nm in all_names) {
  df_m <- model_data[[nm]]$data
  cat(sprintf("\n  %s: N = %s, %d ICD-3 categories ...",
              nm, formatC(nrow(df_m), format = "d", big.mark = ","),
              model_data[[nm]]$n_cats))
  
  t0 <- proc.time()
  
  lmm_models_reml[[nm]] <- tryCatch(
    lmer(log_los ~ 1 + (1 | icd3), data = df_m, REML = TRUE),
    error = function(e) { cat(" [REML FAILED]"); NULL }
  )
  
  lmm_models_ml[[nm]] <- tryCatch(
    lmer(log_los ~ 1 + (1 | icd3), data = df_m, REML = FALSE),
    error = function(e) { cat(" [ML FAILED]"); NULL }
  )
  
  elapsed <- (proc.time() - t0)[3]
  lmm <- lmm_models_reml[[nm]]
  if (!is.null(lmm)) {
    vc <- as.data.frame(VarCorr(lmm))
    s2_u <- vc$vcov[vc$grp == "icd3"]
    s2_e <- vc$vcov[vc$grp == "Residual"]
    icc  <- s2_u / (s2_u + s2_e)
    cat(sprintf(" ICC=%.4f sigma2_u=%.4f sigma2_e=%.4f (%.1fs)\n",
                icc, s2_u, s2_e, elapsed))
  }
}


# ============================================================================
# 9.4) MODEL DIAGNOSTICS
# ============================================================================

cat("\n--- 9.4  Model diagnostics ---\n")


# ---- (a) Variance component extraction and ICC computation -----------------
# ICC = sigma^2_u / (sigma^2_u + sigma^2_e)
# Interpretation: proportion of total log-LOS variance attributable to
# differences between ICD-3 diagnosis groups.

cat("\n  (a) Variance components and ICC\n\n")

variance_components <- data.frame(
  dataset   = character(0),
  k_level   = integer(0),
  sigma2_u  = numeric(0),  # random intercept variance (between-group)
  sigma_u   = numeric(0),  # random intercept SD
  sigma2_e  = numeric(0),  # residual variance (within-group)
  sigma_e   = numeric(0),  # residual SD
  ICC       = numeric(0),
  stringsAsFactors = FALSE
)

for (idx in seq_along(all_names)) {
  nm  <- all_names[idx]
  k   <- all_k[idx]
  lmm <- lmm_models_reml[[nm]]
  if (is.null(lmm)) next
  
  vc   <- as.data.frame(VarCorr(lmm))
  s2_u <- vc$vcov[vc$grp == "icd3"]
  s2_e <- vc$vcov[vc$grp == "Residual"]
  icc  <- s2_u / (s2_u + s2_e)
  
  variance_components <- rbind(variance_components, data.frame(
    dataset   = nm,
    k_level   = k,
    sigma2_u  = s2_u,
    sigma_u   = sqrt(s2_u),
    sigma2_e  = s2_e,
    sigma_e   = sqrt(s2_e),
    ICC       = icc,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("    %s (k=%d): sigma2_u = %.6f | sigma2_e = %.6f | ICC = %.4f\n",
              nm, k, s2_u, s2_e, icc))
}


# ---- (b) Nakagawa & Schielzeth (2013) R-squared ----------------------------
# Marginal  R^2: variance explained by fixed effects / total
# Conditional R^2: variance explained by fixed + random / total
#
# For an intercept-only model (no varying fixed predictor):
#   Marginal R^2 = 0 (the intercept captures no cross-observation variance)
#   Conditional R^2 = ICC (= sigma^2_u / [sigma^2_u + sigma^2_e])
#
# We compute this manually to avoid extra package dependencies.

cat("\n  (b) Nakagawa & Schielzeth R-squared\n\n")

r_squared <- data.frame(
  dataset    = character(0),
  k_level    = integer(0),
  R2_marg    = numeric(0),
  R2_cond    = numeric(0),
  stringsAsFactors = FALSE
)

for (idx in seq_along(all_names)) {
  nm  <- all_names[idx]
  k   <- all_k[idx]
  lmm <- lmm_models_reml[[nm]]
  if (is.null(lmm)) next
  
  vc   <- as.data.frame(VarCorr(lmm))
  s2_u <- vc$vcov[vc$grp == "icd3"]
  s2_e <- vc$vcov[vc$grp == "Residual"]
  
  # Fixed-effect variance: Var(X * beta) -- for intercept-only, this is 0
  # because there is no predictor varying across observations.
  s2_f <- 0
  
  total <- s2_f + s2_u + s2_e
  r2_m  <- s2_f / total          # = 0 for intercept-only
  r2_c  <- (s2_f + s2_u) / total # = ICC for intercept-only
  
  r_squared <- rbind(r_squared, data.frame(
    dataset    = nm,
    k_level    = k,
    R2_marg    = r2_m,
    R2_cond    = r2_c,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("    %s: R2_marginal = %.4f | R2_conditional = %.4f\n",
              nm, r2_m, r2_c))
}


# ---- (c) Residual normality check (Shapiro-Wilk on subsample) -------------
# With N > 700k, a full Shapiro-Wilk test is not feasible.  We use a
# random subsample of 5000 residuals.  This is standard practice for
# large-N LMMs (Pinheiro & Bates, 2000).

cat("\n  (c) Residual normality (Shapiro-Wilk on 5000 subsample)\n\n")

set.seed(42)

normality_results <- list()
for (nm in all_names) {
  lmm <- lmm_models_reml[[nm]]
  if (is.null(lmm)) next
  
  res_all <- residuals(lmm)
  n_res   <- length(res_all)
  sub_idx <- sample.int(n_res, min(5000, n_res))
  sw      <- shapiro.test(res_all[sub_idx])
  
  normality_results[[nm]] <- list(W = sw$statistic, p = sw$p.value)
  cat(sprintf("    %s: W = %.6f, p = %.2e => %s\n",
              nm, sw$statistic, sw$p.value,
              ifelse(sw$p.value < 0.05,
                     "non-normal (expected with large N; CLT applies)",
                     "normal")))
}


# ---- (d) Residual QQ-plots (publication-quality) ---------------------------

cat("\n  (d) Residual QQ-plots\n")

for (nm in all_names) {
  lmm <- lmm_models_reml[[nm]]
  if (is.null(lmm)) next
  
  res_all <- residuals(lmm)
  n_res   <- length(res_all)
  # Subsample for visual clarity (max 10000 points)
  sub_idx <- sample.int(n_res, min(10000, n_res))
  
  pdf(file.path(OUTPUT_DIR, paste0("Fig9_QQ_residuals_", nm, ".pdf")),
      width = PLOT_WIDTH, height = PLOT_HEIGHT)
  par(mar = c(5, 6, 4, 2), family = "serif")
  
  qqnorm(res_all[sub_idx],
         main = paste0("Normal QQ-Plot of LMM Residuals: ", nm),
         xlab = "Theoretical Quantiles",
         ylab = "Sample Quantiles (log-LOS residuals)",
         pch = 16, cex = 0.3, col = adjustcolor("#404040", alpha.f = 0.3),
         cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1, las = 1)
  qqline(res_all[sub_idx], col = "#C00000", lwd = 2.5)
  
  vc   <- as.data.frame(VarCorr(lmm))
  s2_u <- vc$vcov[vc$grp == "icd3"]
  s2_e <- vc$vcov[vc$grp == "Residual"]
  icc  <- s2_u / (s2_u + s2_e)
  
  legend("topleft", bty = "n", cex = 1.0,
         legend = c(
           sprintf("N = %s (10k subsample)",
                   formatC(n_res, format = "d", big.mark = ",")),
           sprintf("ICC = %.4f", icc),
           sprintf("sigma_e = %.4f", sqrt(s2_e))))
  
  dev.off()
  cat(sprintf("    Fig9_QQ_residuals_%s.pdf\n", nm))
}


# ---- (e) Random effect normality (QQ-plot of BLUPs) -----------------------
# If the random effects are approximately normal, the BLUPs should align
# with the theoretical normal quantiles.

cat("\n  (e) Random effect (BLUP) QQ-plots\n")

for (nm in all_names) {
  lmm <- lmm_models_reml[[nm]]
  if (is.null(lmm)) next
  
  blups <- ranef(lmm)$icd3[, 1]
  
  pdf(file.path(OUTPUT_DIR, paste0("Fig9e_QQ_BLUPs_", nm, ".pdf")),
      width = PLOT_WIDTH, height = PLOT_HEIGHT)
  par(mar = c(5, 6, 4, 2), family = "serif")
  
  qqnorm(blups,
         main = paste0("Normal QQ-Plot of ICD-3 Random Intercepts: ", nm),
         xlab = "Theoretical Quantiles",
         ylab = "BLUP (random intercept deviation)",
         pch = 19, cex = 1.0, col = adjustcolor("#2E75B6", alpha.f = 0.6),
         cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1, las = 1)
  qqline(blups, col = "#C00000", lwd = 2.5)
  
  legend("topleft", bty = "n", cex = 1.0,
         legend = c(sprintf("%d ICD-3 groups", length(blups)),
                    sprintf("SD(BLUPs) = %.4f", sd(blups))))
  dev.off()
  cat(sprintf("    Fig9e_QQ_BLUPs_%s.pdf\n", nm))
}


# ---- (f) Prediction accuracy: RMSE and MAE --------------------------------
# Computed on the modeling data (in-sample) for each dataset.
# Back-transformed from log scale to original LOS days via exp(pred) - 1.

cat("\n  (f) Prediction accuracy (in-sample)\n\n")

prediction_accuracy <- data.frame(
  dataset        = character(0),
  k_level        = integer(0),
  RMSE_log       = numeric(0),
  MAE_log        = numeric(0),
  RMSE_days      = numeric(0),
  MAE_days       = numeric(0),
  stringsAsFactors = FALSE
)

for (idx in seq_along(all_names)) {
  nm  <- all_names[idx]
  k   <- all_k[idx]
  lmm <- lmm_models_reml[[nm]]
  if (is.null(lmm)) next
  
  df_m   <- model_data[[nm]]$data
  pred   <- fitted(lmm)                  # on log scale
  resid  <- df_m$log_los - pred
  rmse_l <- sqrt(mean(resid^2))
  mae_l  <- mean(abs(resid))
  
  # Back-transform to days
  pred_days  <- exp(pred) - 1
  resid_days <- df_m$los_days - pred_days
  rmse_d     <- sqrt(mean(resid_days^2))
  mae_d      <- mean(abs(resid_days))
  
  prediction_accuracy <- rbind(prediction_accuracy, data.frame(
    dataset        = nm,
    k_level        = k,
    RMSE_log       = rmse_l,
    MAE_log        = mae_l,
    RMSE_days      = rmse_d,
    MAE_days       = mae_d,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("    %s: RMSE = %.4f (log) / %.2f (days) | MAE = %.4f (log) / %.2f (days)\n",
              nm, rmse_l, rmse_d, mae_l, mae_d))
}


# ============================================================================
# 9.5) EXTRACT FIXED EFFECTS AND BLUPs
# ============================================================================
# Fixed effects: grand mean intercept and its SE/CI
# Random effects: BLUPs (conditional modes) for each ICD-3 category

cat("\n--- 9.5  Extracting fixed effects and BLUPs ---\n")

fixed_effects <- data.frame(
  dataset       = character(0),
  k_level       = integer(0),
  intercept     = numeric(0),
  intercept_se  = numeric(0),
  intercept_lo  = numeric(0),
  intercept_hi  = numeric(0),
  intercept_exp = numeric(0),  # back-transformed: exp(intercept) - 1 = grand mean LOS (days)
  stringsAsFactors = FALSE
)

blup_tables <- list()   # BLUPs per dataset

for (idx in seq_along(all_names)) {
  nm  <- all_names[idx]
  k   <- all_k[idx]
  lmm <- lmm_models_reml[[nm]]
  if (is.null(lmm)) next
  
  fe    <- fixef(lmm)
  fe_se <- sqrt(diag(vcov(lmm)))
  fe_lo <- fe - 1.96 * fe_se
  fe_hi <- fe + 1.96 * fe_se
  
  fixed_effects <- rbind(fixed_effects, data.frame(
    dataset       = nm,
    k_level       = k,
    intercept     = fe[["(Intercept)"]],
    intercept_se  = fe_se[["(Intercept)"]],
    intercept_lo  = fe_lo[["(Intercept)"]],
    intercept_hi  = fe_hi[["(Intercept)"]],
    intercept_exp = exp(fe[["(Intercept)"]]) - 1,
    stringsAsFactors = FALSE
  ))
  
  # Extract BLUPs
  re   <- ranef(lmm)$icd3
  blup <- data.frame(
    icd3       = rownames(re),
    blup       = re[, 1],
    fitted_log = fe[["(Intercept)"]] + re[, 1],     # group-specific mean (log)
    fitted_los = exp(fe[["(Intercept)"]] + re[, 1]) - 1,  # back-transformed (days)
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
  blup_tables[[nm]] <- blup
  
  cat(sprintf("  %s: intercept = %.4f (SE=%.4f) => grand mean LOS = %.2f days | %d BLUPs\n",
              nm, fe[["(Intercept)"]], fe_se[["(Intercept)"]],
              exp(fe[["(Intercept)"]]) - 1, nrow(blup)))
}


# ============================================================================
# 10) REPRODUCIBILITY ASSESSMENT (Section 2.5.4)
# ============================================================================

cat("\n")
cat("####################################################################\n")
cat("#  SECTION 2.5.4 -- Reproducibility Assessment (LMM)              #\n")
cat("####################################################################\n\n")

blup_d0 <- blup_tables[["D0"]]


# ---- 10.1  BLUP Concordance: Spearman rho ---------------------------------
# For shared ICD-3 categories between D0 and Dk, compute Spearman rank
# correlation of the BLUP vectors.  rho = 1.0 => perfect rank preservation
# of diagnosis-specific LOS effects.

cat("--- 10.1  BLUP concordance (Spearman rho) ---\n")

concordance <- data.frame(
  k_level      = integer(0),
  n_shared     = integer(0),
  n_D0_only    = integer(0),
  n_Dk_only    = integer(0),
  spearman_rho = numeric(0),
  rho_p        = numeric(0),
  pearson_r    = numeric(0),
  stringsAsFactors = FALSE)

for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]
  blup_dk <- blup_tables[[lvl]]
  if (is.null(blup_d0) || is.null(blup_dk)) next
  
  shared   <- intersect(blup_d0$icd3, blup_dk$icd3)
  d0_only  <- setdiff(blup_d0$icd3, blup_dk$icd3)
  dk_only  <- setdiff(blup_dk$icd3, blup_d0$icd3)
  
  if (length(shared) < 3) {
    cat(sprintf("  k=%d: only %d shared categories, skipping\n", k, length(shared)))
    next
  }
  
  v0 <- blup_d0$blup[match(shared, blup_d0$icd3)]
  vk <- blup_dk$blup[match(shared, blup_dk$icd3)]
  
  rt_sp <- cor.test(v0, vk, method = "spearman", exact = FALSE)
  rt_pe <- cor.test(v0, vk, method = "pearson")
  
  concordance <- rbind(concordance, data.frame(
    k_level      = k,
    n_shared     = length(shared),
    n_D0_only    = length(d0_only),
    n_Dk_only    = length(dk_only),
    spearman_rho = rt_sp$estimate,
    rho_p        = rt_sp$p.value,
    pearson_r    = rt_pe$estimate,
    stringsAsFactors = FALSE))
  
  cat(sprintf("  k=%d (%s): %d shared | rho_S = %.4f (p = %.2e) | r = %.4f | %d lost | %d new\n",
              k, lvl, length(shared), rt_sp$estimate, rt_sp$p.value,
              rt_pe$estimate, length(d0_only), length(dk_only)))
}


# ---- 10.2  BLUP Attenuation: BLUP_Dk / BLUP_D0 ---------------------------
# Ratio measures whether diagnosis-specific effects are preserved, shrunk,
# or inflated after anonymization.
# Ratio ~1 = preserved; <1 = shrinkage/attenuation; >1 = inflation.
# Only computed for categories where both D0 and Dk have nonzero BLUPs.

cat("\n--- 10.2  BLUP attenuation ---\n")

attenuation_data <- list()
for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]
  blup_dk <- blup_tables[[lvl]]
  if (is.null(blup_d0) || is.null(blup_dk)) next
  
  shared <- intersect(blup_d0$icd3, blup_dk$icd3)
  if (length(shared) == 0) next
  
  v0 <- blup_d0$blup[match(shared, blup_d0$icd3)]
  vk <- blup_dk$blup[match(shared, blup_dk$icd3)]
  
  # Attenuation ratio: only meaningful for categories with nonzero D0 effect
  # Threshold: exclude BLUPs very close to zero to avoid division artifacts
  nonzero <- abs(v0) > 1e-4
  v0_nz   <- v0[nonzero]
  vk_nz   <- vk[nonzero]
  ratio   <- vk_nz / v0_nz
  
  attenuation_data[[lvl]] <- data.frame(
    icd3 = shared[nonzero], blup_d0 = v0_nz, blup_dk = vk_nz,
    ratio = ratio, k_level = k,
    stringsAsFactors = FALSE)
  
  cat(sprintf("  k=%d: median ratio = %.4f | IQR [%.4f, %.4f] | %d categories\n",
              k, median(ratio), quantile(ratio, 0.25), quantile(ratio, 0.75),
              length(ratio)))
}


# ---- 10.3  ICC Trajectory -------------------------------------------------
# The ICC is the central metric: it directly measures how much of the total
# log-LOS variance is attributable to ICD-3 group membership.
# A decline in ICC with increasing k indicates anonymization is eroding
# the explanatory power of diagnosis codes.

cat("\n--- 10.3  ICC trajectory across anonymization levels ---\n\n")

icc_trajectory <- variance_components[, c("dataset", "k_level", "ICC",
                                          "sigma2_u", "sigma2_e")]
print(icc_trajectory, row.names = FALSE)

if (nrow(icc_trajectory) >= 2) {
  icc_d0 <- icc_trajectory$ICC[icc_trajectory$k_level == 0]
  if (length(icc_d0) == 1) {
    icc_trajectory$ICC_retention <- 100 * icc_trajectory$ICC / icc_d0
    cat(sprintf("\n  ICC retention (%%  of D0 ICC):\n"))
    for (j in seq_len(nrow(icc_trajectory))) {
      cat(sprintf("    %s: %.2f%%\n", icc_trajectory$dataset[j],
                  icc_trajectory$ICC_retention[j]))
    }
  }
}


# ---- 10.4  Fixed Effect (Grand Mean) Stability ----------------------------
# The fixed intercept represents the overall mean of log(LOS+1).
# Anonymization via suppression removes primarily rare/outlier combinations,
# which may systematically shift the grand mean.

cat("\n--- 10.4  Fixed effect (grand mean) stability ---\n\n")

print(fixed_effects[, c("dataset", "k_level", "intercept", "intercept_se",
                        "intercept_lo", "intercept_hi", "intercept_exp")],
      row.names = FALSE)

if (nrow(fixed_effects) >= 2) {
  int_d0 <- fixed_effects$intercept[fixed_effects$k_level == 0]
  if (length(int_d0) == 1) {
    cat("\n  Intercept shift from D0:\n")
    for (j in seq_len(nrow(fixed_effects))) {
      shift <- fixed_effects$intercept[j] - int_d0
      cat(sprintf("    %s: %+.6f (%.4f%% change)\n",
                  fixed_effects$dataset[j], shift,
                  100 * shift / abs(int_d0)))
    }
  }
}


# ---- 10.5  Model Fit Comparison (AIC, BIC, logLik) ------------------------
# Using ML-fitted models for comparable AIC/BIC (REML-based AIC/BIC are
# only comparable across models with identical fixed effects AND identical
# data -- which holds here, but ML is conventionally preferred for AIC).

cat("\n--- 10.5  Model fit comparison (ML-based AIC/BIC) ---\n\n")

fit_comparison <- data.frame(
  dataset   = character(0),
  k_level   = integer(0),
  AIC       = numeric(0),
  BIC       = numeric(0),
  logLik    = numeric(0),
  deviance  = numeric(0),
  N_obs     = integer(0),
  N_groups  = integer(0),
  stringsAsFactors = FALSE
)

for (idx in seq_along(all_names)) {
  nm  <- all_names[idx]
  k   <- all_k[idx]
  lmm_ml <- lmm_models_ml[[nm]]
  if (is.null(lmm_ml)) next
  
  ll   <- as.numeric(logLik(lmm_ml))
  smry <- summary(lmm_ml)
  
  fit_comparison <- rbind(fit_comparison, data.frame(
    dataset   = nm,
    k_level   = k,
    AIC       = AIC(lmm_ml),
    BIC       = BIC(lmm_ml),
    logLik    = ll,
    deviance  = -2 * ll,
    N_obs     = nobs(lmm_ml),
    N_groups  = smry$ngrps[["icd3"]],
    stringsAsFactors = FALSE
  ))
}

print(fit_comparison, row.names = FALSE)


# ---- 10.6  Corroboration with quantify_anonymization_impact() -------------
# Map KS D (los_days) and Jaccard (ICD) from Part 1 against BLUP concordance.
# This validates whether the distributional impact function correctly predicts
# the direction and magnitude of statistical reproducibility loss.

cat("\n--- 10.6  Corroboration: Part 1 metrics vs. BLUP concordance ---\n")

if (nrow(concordance) == length(ANON_LEVELS)) {
  corroboration <- data.frame(
    k_level          = concordance$k_level,
    spearman_rho     = concordance$spearman_rho,
    pearson_r        = concordance$pearson_r,
    los_KS_D         = summary_metrics$los_KS_D,
    icd_jaccard      = summary_metrics$icd_jaccard,
    icd_cramer_v     = summary_metrics$icd_cramer_v,
    records_lost_pct = summary_metrics$records_lost_pct,
    ICC              = variance_components$ICC[variance_components$k_level %in% K_VALUES],
    stringsAsFactors = FALSE)
  print(corroboration, row.names = FALSE)
  cat("\n  Expected pattern: KS D rises -> rho falls; Jaccard falls -> rho falls; ICC falls\n")
} else {
  cat("  Insufficient concordance data.\n")
}


# ============================================================================
# 11) SCHOLARLY PLOTS -- LMM Regression & Reproducibility
# ============================================================================

cat("\n=== Generating LMM regression & reproducibility plots ===\n")


# --- PLOT 9b: ICD Prevalence + LOS Distribution Panel (D0/D5/D10/D15) ------
# A 2x4 panel comparing all four datasets side by side:
#   Top row:    Top-15 ICD code prevalences (raw codes as exported by ARX)
#   Bottom row: LOS kernel density with summary statistics

cat("  Generating ICD prevalence + LOS distribution panel...\n")

pdf(file.path(OUTPUT_DIR, "Fig9b_prevalence_los_panel.pdf"),
    width = 16, height = 12)
par(mfrow = c(2, 4), family = "serif")

extract_icd_valid <- function(df) {
  icd_raw <- trimws(as.character(df$main_diagnosis_icd))
  icd_raw[is.na(icd_raw) | icd_raw == "" | icd_raw == "*"] <- NA_character_
  icd_raw
}

icd_by_ds <- lapply(datasets, extract_icd_valid)

icd_d0_valid    <- icd_by_ds[["D0"]][!is.na(icd_by_ds[["D0"]])]
freq_d0_icd     <- sort(table(icd_d0_valid), decreasing = TRUE)
top15           <- names(freq_d0_icd)[1:min(15, length(freq_d0_icd))]
top15_rev       <- rev(top15)

prev_xmax <- 100 * max(freq_d0_icd[top15]) / length(icd_d0_valid) * 1.3

panel_cols  <- c(col_orig, COLS)
panel_names <- c("D0 (Original)", paste0(ANON_LEVELS, " (k=", K_VALUES, ")"))

for (d in seq_along(all_names)) {
  nm <- all_names[d]
  par(mar = c(4, 7, 3, 1))
  
  icd_vec  <- icd_by_ds[[nm]]
  na_ct    <- sum(is.na(icd_vec))
  icd_valid <- icd_vec[!is.na(icd_vec)]
  n_valid  <- length(icd_valid)
  freq_ds  <- table(icd_valid)
  
  prev_pct <- numeric(length(top15_rev))
  names(prev_pct) <- top15_rev
  for (code in top15_rev) {
    if (code %in% names(freq_ds)) {
      prev_pct[code] <- 100 * freq_ds[code] / n_valid
    }
  }
  
  bp <- barplot(prev_pct, horiz = TRUE, las = 1,
                col = adjustcolor(panel_cols[d], alpha.f = 0.7), border = NA,
                xlab = "Prevalence (%)", main = panel_names[d],
                cex.main = 1.2, cex.lab = 1.0, cex.axis = 0.85,
                cex.names = 0.7, xlim = c(0, prev_xmax))
  text(prev_pct, bp, labels = sprintf("%.1f%%", prev_pct),
       pos = 4, cex = 0.65, font = 1)
  mtext(sprintf("N valid = %s | NA = %s",
                formatC(n_valid, format = "d", big.mark = ","),
                formatC(na_ct, format = "d", big.mark = ",")),
        side = 1, line = 2.8, cex = 0.7, font = 3)
}

# Bottom row: LOS density
los_xmax <- quantile(datasets[["D0"]]$los_days, 0.99, na.rm = TRUE)
los_ymax <- 0
for (nm in all_names) {
  lv <- suppressWarnings(as.numeric(datasets[[nm]]$los_days))
  lv <- lv[!is.na(lv) & lv >= 0 & lv <= los_xmax]
  if (length(lv) > 10) los_ymax <- max(los_ymax, max(density(lv)$y))
}

for (d in seq_along(all_names)) {
  nm <- all_names[d]
  par(mar = c(4, 4, 3, 1))
  
  los_raw <- suppressWarnings(as.numeric(datasets[[nm]]$los_days))
  na_ct   <- sum(is.na(los_raw))
  lv      <- los_raw[!is.na(los_raw) & los_raw >= 0 & los_raw <= los_xmax]
  
  if (length(lv) > 10) {
    dens <- density(lv)
    plot(dens, lwd = 3, col = panel_cols[d],
         xlim = c(0, los_xmax), ylim = c(0, los_ymax * 1.05),
         xlab = "Length of Stay (days)", ylab = "Density",
         main = panel_names[d],
         cex.main = 1.2, cex.lab = 1.0, cex.axis = 0.85, las = 1)
    polygon(c(dens$x, rev(dens$x)),
            c(dens$y, rep(0, length(dens$y))),
            col = adjustcolor(panel_cols[d], alpha.f = 0.15), border = NA)
    abline(v = median(lv), lty = 2, col = "gray40", lwd = 1.5)
    legend("topright", bty = "n", cex = 0.75,
           legend = c(
             sprintf("N = %s", formatC(length(lv), format = "d", big.mark = ",")),
             sprintf("NA = %s", formatC(na_ct, format = "d", big.mark = ",")),
             sprintf("Mean = %.1f", mean(lv)),
             sprintf("Median = %.0f", median(lv)),
             sprintf("SD = %.1f", sd(lv))))
  } else {
    plot.new()
    text(0.5, 0.5, paste0(panel_names[d], "\nInsufficient data"), cex = 1.2)
  }
}

dev.off()
cat("  Fig9b_prevalence_los_panel.pdf\n")


# --- PLOT 10: NA Introduction by ARX (Stacked Bar) -------------------------

pdf(file.path(OUTPUT_DIR, "Fig10_NA_introduction.pdf"),
    width = PLOT_WIDTH, height = PLOT_HEIGHT)
par(mar = c(5, 6, 4, 8), family = "serif", xpd = TRUE)

na_mat <- rbind(
  ICD_only = na_documentation$NA_icd_pct,
  LOS_only = na_documentation$NA_los_pct
)
colnames(na_mat) <- na_documentation$dataset

bp <- barplot(na_mat, beside = TRUE,
              col = c("#2E75B6", "#ED7D31"), border = NA,
              ylab = "Values Suppressed (%)",
              xlab = "Dataset",
              main = "ARX-Introduced NA Rates per Column",
              cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1,
              cex.names = 1.2, las = 1,
              ylim = c(0, max(na_mat) * 1.3))

text(bp, na_mat + max(na_mat) * 0.03,
     labels = sprintf("%.1f%%", na_mat),
     cex = 0.9, font = 2)

legend("topright", inset = c(-0.18, 0),
       legend = c("main_diagnosis_icd", "los_days"),
       fill = c("#2E75B6", "#ED7D31"), border = NA,
       cex = 1.0, box.lwd = 0.5)

dev.off()
cat("  Fig10_NA_introduction.pdf\n")


# --- PLOT 11: BLUP Concordance Scatter (D0 vs. Dk) --------------------------

pdf(file.path(OUTPUT_DIR, "Fig11_BLUP_concordance_scatter.pdf"),
    width = 14, height = 5)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2), family = "serif")

for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]
  ad <- attenuation_data[[lvl]]
  if (is.null(ad)) next
  
  lim <- range(c(ad$blup_d0, ad$blup_dk)) * 1.1
  
  plot(ad$blup_d0, ad$blup_dk,
       pch = 19, cex = 0.8, col = adjustcolor(COLS[i], alpha.f = 0.5),
       xlab = expression("BLUP"[D0]),
       ylab = bquote("BLUP"[.(lvl)]),
       main = bquote("D0 vs. " ~ .(lvl) ~ " (k=" ~ .(k) ~ ")"),
       xlim = lim, ylim = lim,
       cex.main = 1.2, cex.lab = 1.1, las = 1)
  abline(0, 1, lty = 2, col = "gray40", lwd = 2)
  
  rho_v <- concordance$spearman_rho[concordance$k_level == k]
  if (length(rho_v) == 1)
    legend("topleft", bty = "n", cex = 1.2, text.font = 2,
           legend = bquote(rho[S] == .(sprintf("%.4f", rho_v))))
}
dev.off()
cat("  Fig11_BLUP_concordance_scatter.pdf\n")


# --- PLOT 12: BLUP Attenuation Ratio Boxplot --------------------------------

if (length(attenuation_data) > 0) {
  pdf(file.path(OUTPUT_DIR, "Fig12_BLUP_attenuation_ratio.pdf"),
      width = PLOT_WIDTH, height = PLOT_HEIGHT)
  par(mar = c(5, 6, 4, 2), family = "serif")
  
  all_att <- do.call(rbind, attenuation_data)
  all_att$k_fac <- factor(all_att$k_level, levels = K_VALUES,
                          labels = paste0("k = ", K_VALUES))
  
  ql <- quantile(all_att$ratio, 0.01); qu <- quantile(all_att$ratio, 0.99)
  att_trim <- all_att[all_att$ratio >= ql & all_att$ratio <= qu, ]
  
  boxplot(ratio ~ k_fac, data = att_trim,
          col = adjustcolor(COLS, alpha.f = 0.4), border = COLS, lwd = 2,
          outline = FALSE,
          ylab = expression("BLUP"[Dk] / "BLUP"[D0]),
          xlab = "k-Anonymity Threshold",
          main = "BLUP Attenuation After Anonymization",
          cex.main = 1.3, cex.lab = 1.2, las = 1)
  abline(h = 1.0, lty = 2, col = "gray30", lwd = 2)
  text(3.3, 1.0, "Perfect preservation", cex = 0.85, col = "gray30", pos = 3)
  
  meds <- tapply(att_trim$ratio, att_trim$k_fac, median)
  text(seq_along(meds), meds, labels = sprintf("%.3f", meds),
       pos = 4, cex = 0.95, font = 2, col = COLS)
  
  dev.off()
  cat("  Fig12_BLUP_attenuation_ratio.pdf\n")
}


# --- PLOT 13: Theta/AIC replaced by Variance Component Trajectory -----------
# Dual panel: (A) ICC trajectory  (B) AIC trajectory

ok_idx <- which(!is.na(variance_components$ICC))

if (length(ok_idx) >= 2) {
  vc_ok  <- variance_components[ok_idx, ]
  fc_ok  <- fit_comparison[match(vc_ok$dataset, fit_comparison$dataset), ]
  col_ok <- c(col_orig, COLS)[ok_idx]
  
  pdf(file.path(OUTPUT_DIR, "Fig14_ICC_AIC_trajectory.pdf"),
      width = 12, height = 6)
  par(mfrow = c(1, 2), mar = c(5, 6, 4, 2), family = "serif")
  
  # Panel A: ICC trajectory
  plot(vc_ok$k_level, vc_ok$ICC,
       type = "b", pch = 19, cex = 2.5, lwd = 3, col = col_ok,
       xlab = "k-Anonymity (0 = Original)",
       ylab = "Intraclass Correlation (ICC)",
       main = "(A) Variance Explained by ICD-3 (ICC)",
       xaxt = "n", las = 1, cex.main = 1.2, cex.lab = 1.1,
       ylim = c(min(vc_ok$ICC) * 0.9, max(vc_ok$ICC) * 1.1))
  axis(1, at = vc_ok$k_level, labels = vc_ok$dataset)
  text(vc_ok$k_level, vc_ok$ICC,
       labels = sprintf("%.4f", vc_ok$ICC),
       pos = 3, cex = 0.9, font = 2, offset = 1.0)
  
  # Panel B: AIC (ML-based)
  if (!is.null(fc_ok) && nrow(fc_ok) >= 2 && !all(is.na(fc_ok$AIC))) {
    plot(fc_ok$k_level, fc_ok$AIC,
         type = "b", pch = 19, cex = 2.5, lwd = 3, col = col_ok,
         xlab = "k-Anonymity (0 = Original)",
         ylab = "AIC (ML)",
         main = "(B) Model AIC",
         xaxt = "n", las = 1, cex.main = 1.2, cex.lab = 1.1)
    axis(1, at = fc_ok$k_level, labels = fc_ok$dataset)
    text(fc_ok$k_level, fc_ok$AIC,
         labels = formatC(round(fc_ok$AIC), format = "d", big.mark = ","),
         pos = 3, cex = 0.9, font = 2, offset = 1.0)
  }
  
  dev.off()
  cat("  Fig14_ICC_AIC_trajectory.pdf\n")
} else {
  cat("  Fig14 SKIPPED: fewer than 2 models converged\n")
}


# --- PLOT 15: Corroboration (3-Panel) ----------------------------------------
# KS D, Jaccard, Cramer V  vs.  Spearman rho of BLUPs

if (nrow(concordance) == length(ANON_LEVELS)) {
  pdf(file.path(OUTPUT_DIR, "Fig15_corroboration.pdf"),
      width = 14, height = 5)
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 2), family = "serif")
  rho_v <- concordance$spearman_rho
  
  # (A) KS D vs rho
  plot(summary_metrics$los_KS_D, rho_v,
       pch = 19, cex = 3, col = COLS,
       xlab = "KS D (los_days)", ylab = expression(rho[S] ~ "(BLUPs)"),
       main = "(A) Distributional Shift vs. Concordance",
       cex.main = 1.1, cex.lab = 1.1, las = 1)
  text(summary_metrics$los_KS_D, rho_v, labels = paste0("k=", K_VALUES),
       pos = 1, cex = 1.0, font = 2, offset = 1.2)
  if (length(rho_v) >= 3) abline(lm(rho_v ~ summary_metrics$los_KS_D),
                                 lty = 2, col = "gray50", lwd = 1.5)
  
  # (B) Jaccard vs rho
  plot(summary_metrics$icd_jaccard, rho_v,
       pch = 19, cex = 3, col = COLS,
       xlab = "ICD Jaccard", ylab = expression(rho[S] ~ "(BLUPs)"),
       main = "(B) Label Preservation vs. Concordance",
       cex.main = 1.1, cex.lab = 1.1, las = 1)
  text(summary_metrics$icd_jaccard, rho_v, labels = paste0("k=", K_VALUES),
       pos = 1, cex = 1.0, font = 2, offset = 1.2)
  if (length(rho_v) >= 3) abline(lm(rho_v ~ summary_metrics$icd_jaccard),
                                 lty = 2, col = "gray50", lwd = 1.5)
  
  # (C) Cramer V vs rho
  plot(summary_metrics$icd_cramer_v, rho_v,
       pch = 19, cex = 3, col = COLS,
       xlab = expression("Cramer's V"), ylab = expression(rho[S] ~ "(BLUPs)"),
       main = "(C) Distortion vs. Concordance",
       cex.main = 1.1, cex.lab = 1.1, las = 1)
  text(summary_metrics$icd_cramer_v, rho_v, labels = paste0("k=", K_VALUES),
       pos = 1, cex = 1.0, font = 2, offset = 1.2)
  if (length(rho_v) >= 3) abline(lm(rho_v ~ summary_metrics$icd_cramer_v),
                                 lty = 2, col = "gray50", lwd = 1.5)
  dev.off()
  cat("  Fig15_corroboration.pdf\n")
}


# ===========================================================================
#  FIGURE: Anonymization Impact on Variance Architecture
# ===========================================================================
# A 2x2 publication-quality panel that tells the complete story of how
# k-anonymization erodes the variance structure of LOS explained by ICD-3.
#
# (A) Variance Decomposition:  Stacked bars -- ICD-3 (random) vs. residual
#     variance for each dataset, with ICC annotated.
# (B) BLUP Concordance Cloud:  All three Dk BLUP vectors plotted against
#     D0 BLUPs (color-coded), showing progressive fidelity loss.
# (C) BLUP Distribution Shift: Overlaid density curves of BLUPs for all
#     four datasets, revealing shrinkage of between-group heterogeneity.
# (D) ICC & Grand Mean LOS Trajectory: Dual-axis line plot of ICC (left)
#     and back-transformed grand mean LOS (right) across k-levels.
#
# This figure, because it integrates four complementary
# perspectives into a single visual narrative -- from macro (variance
# partition, ICC trend) to micro (individual category fidelity).

cat("\n  Generating FIGURE: Variance Architecture Impact...\n")

if (nrow(variance_components) >= 2 && !is.null(blup_d0)) {
  
  pdf(file.path(OUTPUT_DIR, "Fig_variance_architecture.pdf"),
      width = 16, height = 14)
  par(mfrow = c(2, 2), family = "serif")
  
  # === Panel (A): Variance Decomposition Stacked Bars ========================
  par(mar = c(5, 6, 4, 2))
  
  vc_plot <- variance_components
  var_mat <- rbind(
    "ICD-3 (random)" = vc_plot$sigma2_u,
    "Residual"       = vc_plot$sigma2_e
  )
  colnames(var_mat) <- vc_plot$dataset
  
  bp_a <- barplot(var_mat, col = c("#2E75B6", "#B0B0B0"), border = NA,
                  ylab = expression("Variance (" * sigma^2 * ")"),
                  xlab = "", main = "(A) Variance Decomposition: ICD-3 vs. Residual",
                  cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1,
                  cex.names = 1.2, las = 1)
  
  # Annotate ICC on each bar
  total_var <- vc_plot$sigma2_u + vc_plot$sigma2_e
  text(bp_a, total_var + max(total_var) * 0.03,
       labels = sprintf("ICC = %.4f", vc_plot$ICC),
       cex = 1.0, font = 2, col = "#2E75B6")
  
  legend("topright",
         legend = c(expression(sigma[u]^2 ~ "(ICD-3)"),
                    expression(sigma[e]^2 ~ "(Residual)")),
         fill = c("#2E75B6", "#B0B0B0"), border = NA,
         cex = 1.0, bg = "white", box.lwd = 0.5)
  
  
  # === Panel (B): BLUP Concordance Cloud (all k-levels overlaid) =============
  par(mar = c(5, 5, 4, 2))
  
  # Collect all shared BLUPs
  all_v0 <- numeric(0); all_vk <- numeric(0); all_col <- character(0)
  for (i in seq_along(ANON_LEVELS)) {
    lvl <- ANON_LEVELS[i]
    blup_dk <- blup_tables[[lvl]]
    if (is.null(blup_dk)) next
    shared <- intersect(blup_d0$icd3, blup_dk$icd3)
    if (length(shared) < 3) next
    v0 <- blup_d0$blup[match(shared, blup_d0$icd3)]
    vk <- blup_dk$blup[match(shared, blup_dk$icd3)]
    all_v0  <- c(all_v0, v0)
    all_vk  <- c(all_vk, vk)
    all_col <- c(all_col, rep(COLS[i], length(v0)))
  }
  
  if (length(all_v0) > 0) {
    lim <- range(c(all_v0, all_vk)) * 1.1
    plot(all_v0, all_vk,
         pch = 16, cex = 0.7, col = adjustcolor(all_col, alpha.f = 0.4),
         xlab = expression("BLUP"[D0] ~ "(random intercept)"),
         ylab = expression("BLUP"[Dk] ~ "(random intercept)"),
         main = "(B) BLUP Fidelity: D0 vs. All Anonymized",
         xlim = lim, ylim = lim,
         cex.main = 1.3, cex.lab = 1.1, las = 1)
    abline(0, 1, lty = 2, col = "gray30", lwd = 2.5)
    legend("topleft",
           legend = paste0("k=", K_VALUES, " (", ANON_LEVELS, ")"),
           pch = 16, col = COLS, pt.cex = 1.5,
           cex = 1.0, bg = "white", box.lwd = 0.5)
    
    # Annotate Spearman rho per k-level
    y_ann <- lim[1] + (lim[2] - lim[1]) * c(0.15, 0.10, 0.05)
    for (i in seq_along(ANON_LEVELS)) {
      rho_v <- concordance$spearman_rho[concordance$k_level == K_VALUES[i]]
      if (length(rho_v) == 1) {
        text(lim[1] + (lim[2] - lim[1]) * 0.02, y_ann[i],
             labels = bquote(rho[S]^{k==.(K_VALUES[i])} == .(sprintf("%.4f", rho_v))),
             adj = 0, cex = 0.95, col = COLS[i], font = 2)
      }
    }
  } else {
    plot.new()
    text(0.5, 0.5, "Insufficient shared BLUPs", cex = 1.2)
  }
  
  
  # === Panel (C): BLUP Distribution Shift (Density Overlay) ==================
  par(mar = c(5, 5, 4, 2))
  
  blup_densities <- list()
  y_max_blup <- 0
  x_range <- c(Inf, -Inf)
  
  for (nm in all_names) {
    if (is.null(blup_tables[[nm]])) next
    bvals <- blup_tables[[nm]]$blup
    if (length(bvals) >= 3) {
      dens <- density(bvals)
      blup_densities[[nm]] <- dens
      y_max_blup <- max(y_max_blup, max(dens$y))
      x_range[1] <- min(x_range[1], min(dens$x))
      x_range[2] <- max(x_range[2], max(dens$x))
    }
  }
  
  if (length(blup_densities) > 0) {
    plot(NULL, xlim = x_range, ylim = c(0, y_max_blup * 1.1),
         xlab = "Random Intercept (BLUP)", ylab = "Density",
         main = "(C) Distribution of ICD-3 Random Intercepts",
         cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.1, las = 1)
    
    all_plot_cols <- c(col_orig, COLS)
    all_ltys <- c(1, 2, 3, 4)
    for (d in seq_along(all_names)) {
      nm <- all_names[d]
      if (!is.null(blup_densities[[nm]])) {
        lines(blup_densities[[nm]], lwd = 2.5, col = all_plot_cols[d],
              lty = all_ltys[d])
      }
    }
    
    abline(v = 0, lty = 3, col = "gray50", lwd = 1)
    
    legend("topright",
           legend = c("D0 (Original)", paste0("k=", K_VALUES, " (", ANON_LEVELS, ")")),
           col = c(col_orig, COLS), lwd = 2.5, lty = c(1, 2, 3, 4),
           cex = 0.95, bg = "white", box.lwd = 0.5)
    
    # Annotate SD of BLUPs
    for (d in seq_along(all_names)) {
      nm <- all_names[d]
      if (!is.null(blup_tables[[nm]])) {
        sd_b <- sd(blup_tables[[nm]]$blup)
        mtext(sprintf("%s: SD=%.4f", nm, sd_b), side = 1,
              line = 2.5 + (d - 1) * 0.7, cex = 0.7, col = all_plot_cols[d],
              font = 2, adj = 0)
      }
    }
  }
  
  
  # === Panel (D): ICC & Grand Mean LOS Trajectory (Dual Axis) ================
  par(mar = c(5, 6, 4, 6))
  
  icc_vals <- variance_components$ICC
  gm_vals  <- fixed_effects$intercept_exp  # grand mean LOS in days
  k_axis   <- variance_components$k_level
  
  all_plot_cols_d <- c(col_orig, COLS)
  
  # Left axis: ICC
  plot(k_axis, icc_vals,
       type = "b", pch = 19, cex = 2.5, lwd = 3, col = "#2E75B6",
       xlab = "k-Anonymity Level (0 = Original)",
       ylab = "ICC (ICD-3 Variance Share)",
       main = "(D) ICC & Grand Mean LOS Across k-Levels",
       xaxt = "n", las = 1,
       ylim = c(min(icc_vals) * 0.9, max(icc_vals) * 1.1),
       cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.1)
  axis(1, at = k_axis, labels = variance_components$dataset, cex.axis = 1.1)
  
  text(k_axis, icc_vals, labels = sprintf("%.4f", icc_vals),
       pos = 3, cex = 0.85, font = 2, col = "#2E75B6", offset = 0.8)
  
  # Right axis: Grand Mean LOS (days)
  par(new = TRUE)
  plot(k_axis, gm_vals,
       type = "b", pch = 17, cex = 2.5, lwd = 3, col = "#C00000",
       axes = FALSE, xlab = "", ylab = "",
       xlim = range(k_axis),
       ylim = c(min(gm_vals) * 0.95, max(gm_vals) * 1.05))
  axis(4, las = 1, col = "#C00000", col.axis = "#C00000", cex.axis = 1.1)
  mtext("Grand Mean LOS (days)", side = 4, line = 3.5, cex = 1.1, col = "#C00000")
  
  text(k_axis, gm_vals, labels = sprintf("%.2f", gm_vals),
       pos = 1, cex = 0.85, font = 2, col = "#C00000", offset = 0.8)
  
  legend("top",
         legend = c("ICC (left axis)", "Grand Mean LOS (right axis)"),
         pch = c(19, 17), col = c("#2E75B6", "#C00000"),
         lwd = 3, cex = 1.0, bg = "white", box.lwd = 0.5, horiz = TRUE)
  
  dev.off()
  cat("  Fig_variance_architecture.pdf\n")
}


# ============================================================================
# 12) EXPORT LMM RESULTS
# ============================================================================

cat("\n--- Exporting LMM result tables ---\n")

# --- Table 3: Comprehensive LMM Summary (THE ONE SUMMARY TABLE) -----------
# One row per dataset, all key metrics that a reviewer would expect.

lmm_summary_table <- data.frame(
  dataset         = character(0),
  k_level         = integer(0),
  N_obs           = integer(0),
  N_groups        = integer(0),
  intercept       = numeric(0),
  intercept_SE    = numeric(0),
  intercept_CI95  = character(0),
  grand_mean_LOS  = numeric(0),
  sigma2_icd3     = numeric(0),
  sigma2_residual = numeric(0),
  ICC             = numeric(0),
  R2_conditional  = numeric(0),
  AIC_ML          = numeric(0),
  BIC_ML          = numeric(0),
  logLik_ML       = numeric(0),
  RMSE_log        = numeric(0),
  MAE_log         = numeric(0),
  RMSE_days       = numeric(0),
  MAE_days        = numeric(0),
  BLUP_rho_vs_D0  = numeric(0),
  BLUP_rho_p      = numeric(0),
  BLUP_atten_median = numeric(0),
  stringsAsFactors = FALSE
)

for (idx in seq_along(all_names)) {
  nm  <- all_names[idx]
  k   <- all_k[idx]
  
  # Gather all pieces
  fe_row   <- fixed_effects[fixed_effects$dataset == nm, ]
  vc_row   <- variance_components[variance_components$dataset == nm, ]
  r2_row   <- r_squared[r_squared$dataset == nm, ]
  fc_row   <- fit_comparison[fit_comparison$dataset == nm, ]
  pa_row   <- prediction_accuracy[prediction_accuracy$dataset == nm, ]
  co_row   <- concordance[concordance$k_level == k, ]
  at_row   <- attenuation_data[[nm]]
  
  ci_str <- if (nrow(fe_row) == 1)
    sprintf("[%.4f, %.4f]", fe_row$intercept_lo, fe_row$intercept_hi) else NA
  
  lmm_summary_table <- rbind(lmm_summary_table, data.frame(
    dataset         = nm,
    k_level         = k,
    N_obs           = if (nrow(fc_row) == 1) fc_row$N_obs else NA_integer_,
    N_groups        = if (nrow(fc_row) == 1) fc_row$N_groups else NA_integer_,
    intercept       = if (nrow(fe_row) == 1) fe_row$intercept else NA_real_,
    intercept_SE    = if (nrow(fe_row) == 1) fe_row$intercept_se else NA_real_,
    intercept_CI95  = if (!is.na(ci_str)) ci_str else NA_character_,
    grand_mean_LOS  = if (nrow(fe_row) == 1) fe_row$intercept_exp else NA_real_,
    sigma2_icd3     = if (nrow(vc_row) == 1) vc_row$sigma2_u else NA_real_,
    sigma2_residual = if (nrow(vc_row) == 1) vc_row$sigma2_e else NA_real_,
    ICC             = if (nrow(vc_row) == 1) vc_row$ICC else NA_real_,
    R2_conditional  = if (nrow(r2_row) == 1) r2_row$R2_cond else NA_real_,
    AIC_ML          = if (nrow(fc_row) == 1) fc_row$AIC else NA_real_,
    BIC_ML          = if (nrow(fc_row) == 1) fc_row$BIC else NA_real_,
    logLik_ML       = if (nrow(fc_row) == 1) fc_row$logLik else NA_real_,
    RMSE_log        = if (nrow(pa_row) == 1) pa_row$RMSE_log else NA_real_,
    MAE_log         = if (nrow(pa_row) == 1) pa_row$MAE_log else NA_real_,
    RMSE_days       = if (nrow(pa_row) == 1) pa_row$RMSE_days else NA_real_,
    MAE_days        = if (nrow(pa_row) == 1) pa_row$MAE_days else NA_real_,
    BLUP_rho_vs_D0  = if (nrow(co_row) == 1) co_row$spearman_rho else NA_real_,
    BLUP_rho_p      = if (nrow(co_row) == 1) co_row$rho_p else NA_real_,
    BLUP_atten_median = if (!is.null(at_row)) median(at_row$ratio) else NA_real_,
    stringsAsFactors = FALSE
  ))
}

write.csv(lmm_summary_table,
          file.path(OUTPUT_DIR, "Table3_LMM_summary.csv"),
          row.names = FALSE)
cat("  Table3_LMM_summary.csv\n")

# Print to console
cat("\n  === LMM SUMMARY TABLE ===\n\n")
print(lmm_summary_table, row.names = FALSE)

# --- Table 4: BLUP Concordance ---
write.csv(concordance,
          file.path(OUTPUT_DIR, "Table4_BLUP_concordance.csv"),
          row.names = FALSE)
cat("  Table4_BLUP_concordance.csv\n")

# --- Table 5: Variance Components ---
write.csv(variance_components,
          file.path(OUTPUT_DIR, "Table5_variance_components.csv"),
          row.names = FALSE)
cat("  Table5_variance_components.csv\n")

# --- Table 6: BLUP Attenuation ---
if (length(attenuation_data) > 0) {
  write.csv(do.call(rbind, attenuation_data),
            file.path(OUTPUT_DIR, "Table6_BLUP_attenuation.csv"),
            row.names = FALSE)
  cat("  Table6_BLUP_attenuation.csv\n")
}

# --- Table 7: Fixed Effects ---
write.csv(fixed_effects,
          file.path(OUTPUT_DIR, "Table7_fixed_effects.csv"),
          row.names = FALSE)
cat("  Table7_fixed_effects.csv\n")

# --- Table 8: Prediction Accuracy ---
write.csv(prediction_accuracy,
          file.path(OUTPUT_DIR, "Table8_prediction_accuracy.csv"),
          row.names = FALSE)
cat("  Table8_prediction_accuracy.csv\n")

# --- Table 9: Model Fit Comparison ---
write.csv(fit_comparison,
          file.path(OUTPUT_DIR, "Table9_model_fit_comparison.csv"),
          row.names = FALSE)
cat("  Table9_model_fit_comparison.csv\n")


# ============================================================================
# 13) UPDATE STUDY RESULTS AND RE-SAVE .RDS
# ============================================================================

study_results$na_documentation    <- na_documentation
study_results$lmm_models_reml    <- lmm_models_reml
study_results$lmm_models_ml      <- lmm_models_ml
study_results$variance_components <- variance_components
study_results$r_squared           <- r_squared
study_results$fixed_effects       <- fixed_effects
study_results$blup_tables         <- blup_tables
study_results$concordance         <- concordance
study_results$attenuation         <- attenuation_data
study_results$prediction_accuracy <- prediction_accuracy
study_results$fit_comparison      <- fit_comparison
study_results$normality_results   <- normality_results
study_results$lmm_summary_table  <- lmm_summary_table
study_results$model_info <- lapply(model_data, function(m) {
  m$data <- NULL; m })

saveRDS(study_results, file.path(OUTPUT_DIR, "study_results_full.rds"))
cat("\n  study_results_full.rds updated\n")


# ============================================================================
# 14) CONSOLE SUMMARY
# ============================================================================

cat("\n")
cat("====================================================================\n")
cat("  LMM REPRODUCIBILITY SUMMARY\n")
cat("====================================================================\n\n")
cat(sprintf("  Model  : log(los_days + 1) ~ 1 + (1 | icd3)\n"))
cat(sprintf("  Method : REML (lme4::lmer)\n"))
cat(sprintf("  Ref.   : %s | Min freq: >= %d\n\n", ref_cat, MIN_FREQ))

for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]
  cat(sprintf("  k = %d (%s)\n", k, lvl))
  
  # NA info
  na_row <- na_documentation[na_documentation$dataset == lvl, ]
  cat(sprintf("    NAs introduced : ICD %.1f%% | LOS %.1f%% | complete %.1f%%\n",
              na_row$NA_icd_pct, na_row$NA_los_pct, na_row$N_complete_pct))
  
  # Variance components
  vc_row <- variance_components[variance_components$dataset == lvl, ]
  if (nrow(vc_row) == 1)
    cat(sprintf("    sigma2_u = %.6f | sigma2_e = %.6f | ICC = %.4f\n",
                vc_row$sigma2_u, vc_row$sigma2_e, vc_row$ICC))
  
  # Fixed effect
  fe_row <- fixed_effects[fixed_effects$dataset == lvl, ]
  if (nrow(fe_row) == 1)
    cat(sprintf("    Grand mean LOS = %.2f days (intercept = %.4f)\n",
                fe_row$intercept_exp, fe_row$intercept))
  
  # Concordance
  cr <- concordance[concordance$k_level == k, ]
  if (nrow(cr) == 1)
    cat(sprintf("    Spearman rho (BLUPs) = %.4f (p = %.2e) | %d shared ICD-3\n",
                cr$spearman_rho, cr$rho_p, cr$n_shared))
  
  # Attenuation
  ad <- attenuation_data[[lvl]]
  if (!is.null(ad))
    cat(sprintf("    BLUP attenuation : median = %.4f\n", median(ad$ratio)))
  
  # Prediction accuracy
  pa_row <- prediction_accuracy[prediction_accuracy$dataset == lvl, ]
  if (nrow(pa_row) == 1)
    cat(sprintf("    RMSE = %.4f (log) / %.2f (days) | MAE = %.4f (log) / %.2f (days)\n",
                pa_row$RMSE_log, pa_row$RMSE_days, pa_row$MAE_log, pa_row$MAE_days))
  
  cat("\n")
}

cat("====================================================================\n")
cat("  INTERPRETATION KEY:\n")
cat("  rho ~ 1.0        : BLUP rank ordering preserved\n")
cat("  attenuation ~ 1.0: effect magnitudes preserved\n")
cat("  ICC stable        : variance architecture preserved\n")
cat("  RMSE stable       : predictive utility preserved\n")
cat("  grand mean stable : no systematic bias from suppression\n")
cat("====================================================================\n")
cat("\n=== Sections 2.5.3 / 2.5.4 complete ===\n")
