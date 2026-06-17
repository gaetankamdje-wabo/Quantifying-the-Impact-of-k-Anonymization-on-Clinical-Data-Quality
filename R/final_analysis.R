# ============================================================================
# STUDY: Quantifying the Impact of k-Anonymization on Clinical Data Quality
# Author : Gaetan Kamdje Wabo
# R ver. : >= 4.1.0 (requires: lme4)
#
# Schema (D0..D15): patient_id, encounter_id, admission_year,
#                   main_diagnosis_icd, los_days
# ARX suppression in THIS dataset is CELL-LEVEL: main_diagnosis_icd and
# los_days are replaced with the token "*" (jointly), no rows are deleted,
# admission_year is never suppressed. Suppression accounting is therefore
# done at the cell level, not by counting removed encounters.
# ============================================================================

# ---- 0) CONFIGURATION ------------------------------------------------------
DATA_DIR <- file.path(
  "A:", "HLZ", "Promotionen", "Gaetan Kamdje Wabo",
  "study datasets", "study files", "arx", "data source"
)
if (nzchar(Sys.getenv("KANON_DATA_DIR"))) DATA_DIR <- Sys.getenv("KANON_DATA_DIR")

OUTPUT_DIR <- file.path(DATA_DIR, "study_results_updated")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

PNG_PX <- 1200L
PNG_RES <- 150L
ANON_LEVELS <- c("D5", "D10", "D15")
K_VALUES    <- c(5, 10, 15)
NA_TOKENS   <- c("", "NA", "*")
set.seed(42)
N_BOOT <- 2000L

open_png <- function(name)
  grDevices::png(file.path(OUTPUT_DIR, name), width = PNG_PX, height = PNG_PX, res = PNG_RES)

# ---- 1) SOURCE THE CORE FUNCTION ------------------------------------------
source(file.path(getwd(), "quantify_anonymization_impact().R"))
cat("=== quantify_anonymization_impact() loaded ===\n\n")

# ---- 2) IMPORT DATA --------------------------------------------------------
# Read everything as character so the "*" token survives, then coerce the
# numeric columns. After coercion los_days is numeric (suppressed -> NA), which
# is REQUIRED for the impact function to take its numeric (KS) branch.
cat("--- Loading datasets ---\n")
read_study_csv <- function(filename) {
  fp <- file.path(DATA_DIR, filename)
  if (!file.exists(fp)) stop(paste("File not found:", fp))
  df <- read.csv(fp, sep = ",", stringsAsFactors = FALSE, header = TRUE,
                 na.strings = NA_TOKENS, colClasses = "character")
  df$los_days       <- suppressWarnings(as.numeric(df$los_days))
  df$admission_year <- suppressWarnings(as.integer(df$admission_year))
  cat(sprintf("  %-8s : %d rows x %d cols\n", filename, nrow(df), ncol(df)))
  df
}
df_D0  <- read_study_csv("D0.csv")
df_D5  <- read_study_csv("D5.csv")
df_D10 <- read_study_csv("D10.csv")
df_D15 <- read_study_csv("D15.csv")

datasets  <- list(D0 = df_D0, D5 = df_D5, D10 = df_D10, D15 = df_D15)
all_names <- c("D0", ANON_LEVELS)
all_k     <- c(0, K_VALUES)

cat("\n--- Dataset overview ---\n")
for (nm in names(datasets)) {
  df <- datasets[[nm]]
  cat(sprintf("  %s : N=%s | patients=%s | ICD distinct=%d | LOS range=[%s,%s]\n",
              nm, formatC(nrow(df), format="d", big.mark=","),
              formatC(length(unique(na.omit(df$patient_id))), format="d", big.mark=","),
              length(unique(na.omit(df$main_diagnosis_icd))),
              ifelse(all(is.na(df$los_days)), "NA", as.character(min(df$los_days, na.rm=TRUE))),
              ifelse(all(is.na(df$los_days)), "NA", as.character(max(df$los_days, na.rm=TRUE)))))
}
col_orig <- "#404040"; COLS <- c("#2E75B6", "#ED7D31", "#C00000")

# ---- 3) IMPACT ANALYSIS ----------------------------------------------------
cat("\n=== Running anonymization impact analysis ===\n")
results_los <- results_icd <- list()
for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]; df_anon <- datasets[[lvl]]
  cat(sprintf("\n--- D0 vs. %s (k=%d) ---\n", lvl, k))
  results_los[[lvl]] <- quantify_anonymization_impact(
    df_D0, df_anon, comparison_column_origin = df_D0$los_days,
    comparison_column_anonymized = df_anon$los_days)
  results_icd[[lvl]] <- quantify_anonymization_impact(
    df_D0, df_anon, comparison_column_origin = df_D0$main_diagnosis_icd,
    comparison_column_anonymized = df_anon$main_diagnosis_icd)
  cat(sprintf("  los_days verdict : %s | ks_z=%.2f | reason: %s\n",
              results_los[[lvl]]$adhoc_column_comparison$handling_recommendation$verdict,
              results_los[[lvl]]$adhoc_column_comparison$meta$ks_z,
              results_los[[lvl]]$adhoc_column_comparison$meta$ks_z_reason))
  cat(sprintf("  ICD verdict      : %s\n",
              results_icd[[lvl]]$adhoc_column_comparison$handling_recommendation$verdict))
}

# ---- 4) SUMMARY METRICS ----------------------------------------------------
cat("\n=== Building summary metrics table ===\n")
getm <- function(res, field) as.numeric(res$adhoc_column_comparison$meta[[field]])
summary_metrics <- data.frame(
  k_level = K_VALUES, label = ANON_LEVELS, N_original = nrow(df_D0),
  N_anonymized = sapply(ANON_LEVELS, function(l) nrow(datasets[[l]])),
  los_KS_D = sapply(ANON_LEVELS, function(l) getm(results_los[[l]], "ks_D")),
  los_KS_z = sapply(ANON_LEVELS, function(l) getm(results_los[[l]], "ks_z")),
  los_KS_z_reason = sapply(ANON_LEVELS, function(l)
    results_los[[l]]$adhoc_column_comparison$meta$ks_z_reason),
  los_mean_shift = sapply(ANON_LEVELS, function(l)
    getm(results_los[[l]], "mean_anon") - getm(results_los[[l]], "mean_orig")),
  los_sd_shift = sapply(ANON_LEVELS, function(l)
    getm(results_los[[l]], "sd_anon") - getm(results_los[[l]], "sd_orig")),
  los_verdict = sapply(ANON_LEVELS, function(l)
    results_los[[l]]$adhoc_column_comparison$handling_recommendation$verdict),
  icd_jaccard = sapply(ANON_LEVELS, function(l) getm(results_icd[[l]], "overlap")),
  icd_cramer_v = sapply(ANON_LEVELS, function(l) getm(results_icd[[l]], "cramer_v")),
  icd_verdict = sapply(ANON_LEVELS, function(l)
    results_icd[[l]]$adhoc_column_comparison$handling_recommendation$verdict),
  stringsAsFactors = FALSE, row.names = NULL)
print(summary_metrics[, c("label","k_level","los_KS_D","los_KS_z","icd_jaccard","icd_cramer_v")])
write.csv(summary_metrics, file.path(OUTPUT_DIR, "summary_metrics.csv"), row.names = FALSE)

# ---- 4b) DATA POINT REMOVAL / SUPPRESSION DOCUMENTATION -------------------
# Two distinct mechanisms, both reported:
#   encounters_removed : rows in D0 absent from Dk (row deletion; 0 here)
#   cells_*_suppressed : cells replaced by "*" within surviving rows
cat("\n=== Data point removal / cell suppression documentation ===\n")
id0 <- df_D0$encounter_id
removal_documentation <- data.frame(
  dataset = ANON_LEVELS, k_level = K_VALUES, N_D0 = nrow(df_D0),
  N_Dk = sapply(ANON_LEVELS, function(l) nrow(datasets[[l]])),
  encounters_removed = sapply(ANON_LEVELS, function(l)
    length(setdiff(id0, datasets[[l]]$encounter_id))),
  cells_icd_suppressed = sapply(ANON_LEVELS, function(l)
    sum(is.na(datasets[[l]]$main_diagnosis_icd))),
  cells_los_suppressed = sapply(ANON_LEVELS, function(l)
    sum(is.na(datasets[[l]]$los_days))),
  cells_year_suppressed = sapply(ANON_LEVELS, function(l)
    sum(is.na(datasets[[l]]$admission_year))),
  stringsAsFactors = FALSE)
removal_documentation$encounters_removed_pct <-
  100 * removal_documentation$encounters_removed / removal_documentation$N_D0
removal_documentation$cells_icd_suppressed_pct <-
  100 * removal_documentation$cells_icd_suppressed / removal_documentation$N_D0
removal_documentation$cells_los_suppressed_pct <-
  100 * removal_documentation$cells_los_suppressed / removal_documentation$N_D0
print(removal_documentation, row.names = FALSE)
write.csv(removal_documentation, file.path(OUTPUT_DIR, "data_point_removal.csv"), row.names = FALSE)

# ---- 4c) SUPPRESSION RATES BY ADMISSION YEAR ------------------------------
# Cell-level suppression. An encounter counts as "suppressed" in Dk if its
# main_diagnosis_icd OR its los_days was replaced by "*" (in this dataset the
# two are starred jointly), or the encounter is absent from Dk. admission_year
# is taken from D0 (never suppressed); Dk is aligned to D0 by encounter_id.
# Reported per year: count and % of that year's D0 encounters suppressed in Dk.
# (D0 itself has zero suppression, so this is not a per-year D0-vs-Dk test; D0
# is only the denominator.)
cat("\n=== Suppression by admission year (ICD or LOS starred) ===\n")
years <- sort(unique(na.omit(df_D0$admission_year)))
year_n_d0 <- sapply(years, function(y) sum(df_D0$admission_year == y, na.rm = TRUE))
yr_of <- df_D0$admission_year   # year per D0 encounter, in D0 row order

supp_by_year <- data.frame(admission_year = years, N_D0 = year_n_d0, stringsAsFactors = FALSE)
for (l in ANON_LEVELS) {
  dk  <- datasets[[l]]
  idx <- match(df_D0$encounter_id, dk$encounter_id)   # align Dk to D0
  icd_supp <- is.na(dk$main_diagnosis_icd[idx])       # "*" cell -> NA, or absent
  los_supp <- is.na(dk$los_days[idx])
  supp_flag <- icd_supp | los_supp
  cnt <- tapply(supp_flag, yr_of, sum)[as.character(years)]
  cnt[is.na(cnt)] <- 0
  supp_by_year[[paste0("suppressed_", l)]]     <- as.integer(cnt)
  supp_by_year[[paste0("suppressed_pct_", l)]] <- 100 * cnt / year_n_d0
}
cat(sprintf("  overall suppressed encounters: D5=%d (%.2f%%), D10=%d (%.2f%%), D15=%d (%.2f%%)\n",
            sum(supp_by_year$suppressed_D5),  100*sum(supp_by_year$suppressed_D5)/nrow(df_D0),
            sum(supp_by_year$suppressed_D10), 100*sum(supp_by_year$suppressed_D10)/nrow(df_D0),
            sum(supp_by_year$suppressed_D15), 100*sum(supp_by_year$suppressed_D15)/nrow(df_D0)))
print(supp_by_year, row.names = FALSE)
write.csv(supp_by_year, file.path(OUTPUT_DIR, "suppression_by_year.csv"), row.names = FALSE)

# ---- 4d) SUPPRESSED vs RETAINED COHORT CHARACTERIZATION -------------------
# "Suppressed" = encounters whose ICD/LOS was starred in Dk. Compared against
# the rest, on D0 (true) values, to expose selection bias.
cat("\n=== Suppressed vs retained cohort characterization ===\n")
describe_cohort <- function(df_sub) {
  los <- df_sub$los_days[is.finite(df_sub$los_days)]
  t <- sort(table(df_sub$main_diagnosis_icd[!is.na(df_sub$main_diagnosis_icd)]), decreasing = TRUE)
  data.frame(
    n = nrow(df_sub),
    n_patients = length(unique(na.omit(df_sub$patient_id))),
    los_mean = if (length(los)) mean(los) else NA_real_,
    los_median = if (length(los)) median(los) else NA_real_,
    los_sd = if (length(los) > 1) sd(los) else NA_real_,
    top_icd = if (length(t)) paste0(names(t)[1], " (", round(100*t[1]/sum(t),1), "%)") else NA_character_,
    stringsAsFactors = FALSE)
}
cohort_rows <- list()
for (i in seq_along(ANON_LEVELS)) {
  l <- ANON_LEVELS[i]; k <- K_VALUES[i]; dk <- datasets[[l]]
  idx <- match(id0, dk$encounter_id)
  is_supp <- is.na(idx) | is.na(dk$main_diagnosis_icd[idx]) | is.na(dk$los_days[idx])
  ret <- describe_cohort(df_D0[!is_supp, ]); ret$k_level <- k; ret$cohort <- "retained"
  sup <- describe_cohort(df_D0[ is_supp, ]); sup$k_level <- k; sup$cohort <- "suppressed"
  cohort_rows[[length(cohort_rows)+1]] <- ret
  cohort_rows[[length(cohort_rows)+1]] <- sup
}
cohort_summary <- do.call(rbind, cohort_rows)
cohort_summary <- cohort_summary[, c("k_level","cohort","n","n_patients",
                                     "los_mean","los_median","los_sd","top_icd")]
print(cohort_summary, row.names = FALSE)
write.csv(cohort_summary, file.path(OUTPUT_DIR, "cohort_characterization.csv"), row.names = FALSE)

# ============================================================================
# PART 2: DISTRIBUTIONAL FIGURES (PNG 1200x1200)
# ============================================================================
cat("\n=== Generating distributional figures ===\n")

open_png("Fig_KS_D_los_days.png")
par(mar = c(5,6,4,2), family = "serif")
plot(summary_metrics$k_level, summary_metrics$los_KS_D, type="b", pch=19, cex=2.5, lwd=3,
     col=col_orig, xlab="k-Anonymity Threshold", ylab="Kolmogorov-Smirnov D (Effect Size)",
     main="Distributional Divergence of los_days\nOriginal vs. Anonymized",
     xlim=c(3,17), ylim=c(0, max(c(summary_metrics$los_KS_D,0.15),na.rm=TRUE)*1.3),
     xaxt="n", las=1, cex.main=1.3, cex.lab=1.2, cex.axis=1.1)
axis(1, at=K_VALUES, labels=paste0("k = ",K_VALUES), cex.axis=1.1)
points(summary_metrics$k_level, summary_metrics$los_KS_D, pch=19, cex=2.5, col=COLS)
abline(h=c(0.05,0.10), lty=2, col=c("forestgreen","darkorange"), lwd=1.5)
text(summary_metrics$k_level, summary_metrics$los_KS_D, labels=sprintf("D=%.4f", summary_metrics$los_KS_D),
     pos=3, cex=1.0, font=2, offset=1.2)
dev.off(); cat("  Fig_KS_D_los_days.png\n")

open_png("Fig_ICD_quality_metrics.png")
par(mar=c(5,6,4,6), family="serif")
plot(summary_metrics$k_level, summary_metrics$icd_jaccard, type="b", pch=17, cex=2.5, lwd=3,
     col="#2E75B6", xlab="k-Anonymity Threshold", ylab="Jaccard Overlap Ratio",
     main="ICD Code Quality: Label Preservation vs. Distributional Shift",
     xlim=c(3,17), ylim=c(0,1), xaxt="n", las=1, cex.main=1.3, cex.lab=1.2, cex.axis=1.1)
axis(1, at=K_VALUES, labels=paste0("k = ",K_VALUES), cex.axis=1.1)
text(summary_metrics$k_level, summary_metrics$icd_jaccard, labels=sprintf("J=%.3f", summary_metrics$icd_jaccard),
     pos=1, cex=0.95, font=2, col="#2E75B6", offset=1.2)
par(new=TRUE)
plot(summary_metrics$k_level, summary_metrics$icd_cramer_v, type="b", pch=15, cex=2.5, lwd=3,
     col="#C00000", axes=FALSE, xlab="", ylab="", xlim=c(3,17), ylim=c(0,1))
axis(4, las=1, col="#C00000", col.axis="#C00000", cex.axis=1.1)
mtext("Cramer's V (Effect Size)", side=4, line=3.5, cex=1.2, col="#C00000")
text(summary_metrics$k_level, summary_metrics$icd_cramer_v, labels=sprintf("V=%.3f", summary_metrics$icd_cramer_v),
     pos=3, cex=0.95, font=2, col="#C00000", offset=1.2)
abline(h=c(0.10,0.30), lty=3, col="gray50")
legend("topright", legend=c("Jaccard (left)","Cramer's V (right)"), pch=c(17,15),
       col=c("#2E75B6","#C00000"), lwd=3, cex=1.0, bg="white", box.lwd=0.5)
dev.off(); cat("  Fig_ICD_quality_metrics.png\n")

open_png("Fig_ECDF_los_days.png")
par(mar=c(5,6,4,2), family="serif")
ecdf_orig <- ecdf(df_D0$los_days); x_max <- quantile(df_D0$los_days, .99, na.rm=TRUE)
x_seq <- seq(0, x_max, length.out=1000)
plot(x_seq, ecdf_orig(x_seq), type="l", lwd=3, col=col_orig, xlab="Length of Stay (days)",
     ylab=expression(hat(F)(x)), main="Empirical CDF: los_days\nOriginal vs. k-Anonymized",
     las=1, cex.main=1.3, cex.lab=1.2, cex.axis=1.1)
for (i in seq_along(ANON_LEVELS)) lines(x_seq, ecdf(datasets[[ANON_LEVELS[i]]]$los_days)(x_seq),
                                        lwd=2.5, col=COLS[i], lty=i+1)
legend("bottomright", legend=c("Original (D0)", paste0("k=",K_VALUES)),
       col=c(col_orig,COLS), lwd=c(3,rep(2.5,3)), lty=1:4, cex=1.0, bg="white", box.lwd=0.5)
dev.off(); cat("  Fig_ECDF_los_days.png\n")

open_png("Fig_density_los_days.png")
par(mar=c(5,6,4,2), family="serif")
x_max_d <- quantile(df_D0$los_days, .99, na.rm=TRUE)
d_orig <- density(df_D0$los_days[is.finite(df_D0$los_days) & df_D0$los_days<=x_max_d])
y_max <- max(d_orig$y)
for (lvl in ANON_LEVELS) { v <- datasets[[lvl]]$los_days; v <- v[is.finite(v)&v<=x_max_d]
if (length(v)>1) y_max <- max(y_max, max(density(v)$y)) }
plot(d_orig, lwd=3, col=col_orig, xlim=c(0,x_max_d), ylim=c(0,y_max*1.1),
     xlab="Length of Stay (days)", ylab="Density",
     main="Kernel Density: los_days\nOriginal vs. k-Anonymized", cex.main=1.3, cex.lab=1.2, cex.axis=1.1, las=1)
for (i in seq_along(ANON_LEVELS)) { v <- datasets[[ANON_LEVELS[i]]]$los_days; v <- v[is.finite(v)&v<=x_max_d]
if (length(v)>1) lines(density(v), lwd=2.5, col=COLS[i], lty=i+1) }
legend("topright", legend=c("Original (D0)", paste0("k=",K_VALUES)),
       col=c(col_orig,COLS), lwd=c(3,rep(2.5,3)), lty=1:4, cex=1.0, bg="white", box.lwd=0.5)
dev.off(); cat("  Fig_density_los_days.png\n")

open_png("Fig_ICD_zipf_rank.png")
par(mar=c(5,6,4,2), family="serif")
freq_orig <- sort(table(df_D0$main_diagnosis_icd), decreasing=TRUE); n_codes <- min(length(freq_orig),100)
plot(seq_len(n_codes), as.numeric(freq_orig[seq_len(n_codes)]), type="l", lwd=3, col=col_orig, log="y",
     xlab="ICD Code Rank (most frequent first)", ylab="Frequency (log scale)",
     main="ICD Code Rank-Frequency (Zipf Plot)\nTop 100 Codes", cex.main=1.3, cex.lab=1.2, cex.axis=1.1, las=1)
for (i in seq_along(ANON_LEVELS)) { fa <- sort(table(datasets[[ANON_LEVELS[i]]]$main_diagnosis_icd), decreasing=TRUE)
na_ <- min(length(fa),100); lines(seq_len(na_), as.numeric(fa[seq_len(na_)]), lwd=2.5, col=COLS[i], lty=i+1) }
legend("topright", legend=c("Original (D0)", paste0("k=",K_VALUES)),
       col=c(col_orig,COLS), lwd=c(3,rep(2.5,3)), lty=1:4, cex=1.0, bg="white", box.lwd=0.5)
dev.off(); cat("  Fig_ICD_zipf_rank.png\n")

open_png("Fig_los_mean_sd_shifts.png")
par(mar=c(5,9,4,2), family="serif")
y_pos <- c(1,2,3,5,6,7); vals <- c(summary_metrics$los_mean_shift, summary_metrics$los_sd_shift)
cols2 <- rep(COLS,2); labs <- c(paste0("Mean shift (k=",K_VALUES,")"), paste0("SD shift (k=",K_VALUES,")"))
plot(vals, y_pos, pch=19, cex=2.5, col=cols2, xlim=range(c(vals,0))*c(1.4,1.4), yaxt="n",
     xlab="Shift (days)", ylab="", main="los_days: Mean and SD Shifts After Anonymization",
     cex.main=1.3, cex.lab=1.2, cex.axis=1.1, las=1)
segments(0, y_pos, vals, y_pos, lwd=2.5, col=cols2); abline(v=0, lty=1, col="gray30", lwd=1.5)
axis(2, at=y_pos, labels=labs, las=1, cex.axis=0.95)
text(vals, y_pos, labels=sprintf("%+.3f", vals), pos=ifelse(vals>=0,4,2), cex=0.95, font=2)
abline(h=4, lty=3, col="gray70")
dev.off(); cat("  Fig_los_mean_sd_shifts.png\n")

open_png("Fig_prevalence_los_panel.png")
par(mfrow=c(2,4), family="serif")
icd_by_ds <- lapply(datasets, function(df) df$main_diagnosis_icd)
icd_d0_valid <- icd_by_ds[["D0"]][!is.na(icd_by_ds[["D0"]])]
freq_d0_icd <- sort(table(icd_d0_valid), decreasing=TRUE)
top15 <- names(freq_d0_icd)[1:min(15,length(freq_d0_icd))]; top15_rev <- rev(top15)
prev_xmax <- 100*max(freq_d0_icd[top15])/length(icd_d0_valid)*1.3
panel_cols <- c(col_orig,COLS); panel_names <- c("D0 (Original)", paste0(ANON_LEVELS," (k=",K_VALUES,")"))
for (d in seq_along(all_names)) {
  nm <- all_names[d]; par(mar=c(4,7,3,1))
  icd_vec <- icd_by_ds[[nm]]; na_ct <- sum(is.na(icd_vec)); iv <- icd_vec[!is.na(icd_vec)]; nv <- length(iv)
  fq <- table(iv); pv <- setNames(numeric(length(top15_rev)), top15_rev)
  for (code in top15_rev) if (code %in% names(fq)) pv[code] <- 100*fq[code]/nv
  bp <- barplot(pv, horiz=TRUE, las=1, col=adjustcolor(panel_cols[d],0.7), border=NA,
                xlab="Prevalence (%)", main=panel_names[d], cex.main=1.2, cex.lab=1.0,
                cex.axis=0.85, cex.names=0.7, xlim=c(0,prev_xmax))
  text(pv, bp, labels=sprintf("%.1f%%",pv), pos=4, cex=0.65)
  mtext(sprintf("N valid=%s | NA=%s", formatC(nv,format="d",big.mark=","),
                formatC(na_ct,format="d",big.mark=",")), side=1, line=2.8, cex=0.7, font=3)
}
los_xmax <- quantile(df_D0$los_days, .99, na.rm=TRUE); los_ymax <- 0
for (nm in all_names) { lv <- datasets[[nm]]$los_days; lv <- lv[is.finite(lv)&lv>=0&lv<=los_xmax]
if (length(lv)>10) los_ymax <- max(los_ymax, max(density(lv)$y)) }
for (d in seq_along(all_names)) {
  nm <- all_names[d]; par(mar=c(4,4,3,1)); lr <- datasets[[nm]]$los_days; na_ct <- sum(is.na(lr))
  lv <- lr[is.finite(lr)&lr>=0&lr<=los_xmax]
  if (length(lv)>10) {
    dn <- density(lv)
    plot(dn, lwd=3, col=panel_cols[d], xlim=c(0,los_xmax), ylim=c(0,los_ymax*1.05),
         xlab="Length of Stay (days)", ylab="Density", main=panel_names[d],
         cex.main=1.2, cex.lab=1.0, cex.axis=0.85, las=1)
    polygon(c(dn$x,rev(dn$x)), c(dn$y,rep(0,length(dn$y))), col=adjustcolor(panel_cols[d],0.15), border=NA)
    abline(v=median(lv), lty=2, col="gray40", lwd=1.5)
    legend("topright", bty="n", cex=0.75,
           legend=c(sprintf("N=%s",formatC(length(lv),format="d",big.mark=",")),
                    sprintf("NA=%s",formatC(na_ct,format="d",big.mark=",")),
                    sprintf("Mean=%.1f",mean(lv)), sprintf("Median=%.0f",median(lv)), sprintf("SD=%.1f",sd(lv))))
  } else { plot.new(); text(0.5,0.5,paste0(panel_names[d],"\nInsufficient data")) }
}
dev.off(); cat("  Fig_prevalence_los_panel.png\n")

# ============================================================================
# PART 3: LINEAR MIXED MODEL  log(los+1) ~ year + (1|icd3) + (1|patient_id)
# Each dataset fit ONCE with REML. AIC/BIC are REML-based and comparable here
# because the fixed-effects structure is identical across datasets. Fitting
# once (not REML+ML) roughly halves runtime.
# ============================================================================
cat("\n####################################################################\n")
cat("#  LMM Regression & Reproducibility                                 #\n")
cat("####################################################################\n\n")
suppressMessages(library(lme4))

# ARX replaces suppressed quasi-identifier cells with the token "*", which is
# read in as NA (see NA_TOKENS). In the anonymized datasets (D5/D10/D15) every
# NA in main_diagnosis_icd or los_days is therefore an ARX-suppressed cell, not
# natural missingness. admission_year is never a quasi-identifier and is never
# suppressed. D0 is the original dataset, so any NA there reflects true source
# missingness and serves only as the baseline. We label the counts as
# "suppressed" for k>0 and report D0 separately as the baseline.
cat("--- ARX suppression documentation (\"*\" cells read as NA) ---\n\n")
na_documentation <- data.frame()
for (idx in seq_along(all_names)) {
  nm <- all_names[idx]; df <- datasets[[nm]]; n <- nrow(df)
  supp_icd  <- sum(is.na(df$main_diagnosis_icd))
  supp_los  <- sum(is.na(df$los_days))
  supp_year <- sum(is.na(df$admission_year))
  supp_any  <- sum(is.na(df$main_diagnosis_icd)|is.na(df$los_days)|is.na(df$admission_year))
  na_documentation <- rbind(na_documentation, data.frame(
    dataset=nm, k_level=all_k[idx],
    source=ifelse(all_k[idx]==0, "baseline (true missingness)", "ARX suppression (*)"),
    N_rows=n,
    suppressed_icd=supp_icd,   suppressed_icd_pct=100*supp_icd/n,
    suppressed_los=supp_los,   suppressed_los_pct=100*supp_los/n,
    suppressed_year=supp_year, suppressed_year_pct=100*supp_year/n,
    suppressed_any=supp_any,   suppressed_any_pct=100*supp_any/n,
    N_complete=n-supp_any,     N_complete_pct=100*(n-supp_any)/n, stringsAsFactors=FALSE))
}
print(na_documentation, row.names = FALSE)
write.csv(na_documentation, file.path(OUTPUT_DIR, "NA_documentation.csv"), row.names = FALSE)

# Build per-dataset modeling frames once. Each is light (5 columns); the heavy
# object is the fitted model, which we create, extract from, and discard one at
# a time so memory never holds more than a single lmer fit (avoids OOM on the
# ~370k-patient random effect).
cat("\n--- Preparing modeling data (all encounters retained) ---\n")
prepare_model_data <- function(df, label) {
  icd3 <- substr(df$main_diagnosis_icd, 1, 3)
  icd3[!grepl("^[A-Z][0-9]{2}$", icd3)] <- NA_character_
  keep <- is.finite(df$los_days) & df$los_days >= 0 &
    !is.na(icd3) & !is.na(df$admission_year) & !is.na(df$patient_id)
  out <- data.frame(icd3=factor(icd3[keep]), patient_id=factor(df$patient_id[keep]),
                    year=factor(df$admission_year[keep]), los_days=df$los_days[keep],
                    log_los=log(df$los_days[keep]+1), stringsAsFactors=FALSE)
  cat(sprintf("  %s: %d rows -> %d modeled (%d ICD-3, %d patients, %d years)\n",
              label, nrow(df), nrow(out), nlevels(out$icd3), nlevels(out$patient_id), nlevels(out$year)))
  list(data=out, n_model=nrow(out), n_cats=nlevels(out$icd3), n_pat=nlevels(out$patient_id))
}
model_data <- lapply(all_names, function(nm) prepare_model_data(datasets[[nm]], nm))
names(model_data) <- all_names
ref_cat <- names(which.max(table(model_data[["D0"]]$data$icd3)))
cat(sprintf("\n  Reference ICD-3 (most frequent in D0): %s\n", ref_cat))
for (nm in all_names) if (ref_cat %in% levels(model_data[[nm]]$data$icd3))
  model_data[[nm]]$data$icd3 <- relevel(model_data[[nm]]$data$icd3, ref=ref_cat)

# ---- Single fit/extract/free loop -----------------------------------------
# All downstream tables are accumulated here from lightweight extracts. The
# model object is removed and gc() called before the next dataset is fit.
cat("\n--- Fitting LMMs (REML, one fit per dataset; freed after extraction) ---\n")
variance_components  <- data.frame()
r_squared            <- data.frame()
prediction_accuracy  <- data.frame()
fixed_effects        <- data.frame()
fit_comparison       <- data.frame()
normality_results    <- list()
blup_tables          <- list()

for (idx in seq_along(all_names)) {
  nm <- all_names[idx]; df_m <- model_data[[nm]]$data
  ckpt <- file.path(OUTPUT_DIR, paste0(".ckpt_", nm, ".rds"))
  if (file.exists(ckpt)) {                       # resume: reuse completed extract
    e <- readRDS(ckpt)
    variance_components <- rbind(variance_components, e$vc_row)
    r_squared           <- rbind(r_squared, e$r2_row)
    prediction_accuracy <- rbind(prediction_accuracy, e$pa_row)
    fixed_effects       <- rbind(fixed_effects, e$fe_row)
    fit_comparison      <- rbind(fit_comparison, e$fc_row)
    normality_results[[nm]] <- e$norm
    blup_tables[[nm]]       <- e$blup
    cat(sprintf("  %s: restored from checkpoint (ICC_icd=%.4f)\n", nm, e$vc_row$ICC_icd3))
    next
  }
  cat(sprintf("  %s: N=%s ...", nm, formatC(nrow(df_m), format="d", big.mark=",")))
  t0 <- proc.time()
  lmm <- tryCatch(
    lmer(log_los ~ year + (1|icd3) + (1|patient_id), data=df_m, REML=TRUE,
         control=lmerControl(optimizer="bobyqa", check.conv.singular=.makeCC("ignore",1e-4))),
    error=function(e){ cat(sprintf(" [FAIL: %s]", conditionMessage(e))); NULL })
  if (is.null(lmm)) { cat("\n"); next }
  
  vc <- as.data.frame(VarCorr(lmm))
  s2_icd <- vc$vcov[vc$grp=="icd3"]; s2_pat <- vc$vcov[vc$grp=="patient_id"]; s2_e <- vc$vcov[vc$grp=="Residual"]
  tot <- s2_icd+s2_pat+s2_e
  vc_row <- data.frame(dataset=nm, k_level=all_k[idx], sigma2_icd3=s2_icd, sigma2_patient=s2_pat,
                       sigma2_residual=s2_e, ICC_icd3=s2_icd/tot, ICC_patient=s2_pat/tot,
                       ICC_icd3_plus_patient=(s2_icd+s2_pat)/tot, stringsAsFactors=FALSE)
  
  s2_f <- var(as.vector(model.matrix(lmm) %*% fixef(lmm))); totf <- s2_f+s2_icd+s2_pat+s2_e
  r2_row <- data.frame(dataset=nm, k_level=all_k[idx], R2_marg=s2_f/totf,
                       R2_cond=(s2_f+s2_icd+s2_pat)/totf, stringsAsFactors=FALSE)
  
  res_all <- residuals(lmm)
  s5  <- res_all[sample.int(length(res_all), min(5000, length(res_all)))]
  s10 <- res_all[sample.int(length(res_all), min(10000, length(res_all)))]
  sw <- shapiro.test(s5); norm <- list(W=unname(sw$statistic), p=sw$p.value)
  
  pred <- fitted(lmm); resid_log <- df_m$log_los - pred
  pred_days <- exp(pred)-1; resid_days <- df_m$los_days - pred_days
  pa_row <- data.frame(dataset=nm, k_level=all_k[idx],
                       RMSE_log=sqrt(mean(resid_log^2)), MAE_log=mean(abs(resid_log)),
                       RMSE_days=sqrt(mean(resid_days^2)), MAE_days=mean(abs(resid_days)), stringsAsFactors=FALSE)
  
  fe <- fixef(lmm); fe_se <- sqrt(diag(as.matrix(vcov(lmm))))
  b0 <- fe[["(Intercept)"]]; se0 <- fe_se[["(Intercept)"]]
  fe_row <- data.frame(dataset=nm, k_level=all_k[idx], intercept=b0, intercept_se=se0,
                       intercept_lo=b0-1.96*se0, intercept_hi=b0+1.96*se0, intercept_exp=exp(b0)-1, stringsAsFactors=FALSE)
  re <- ranef(lmm)$icd3
  blup <- data.frame(icd3=rownames(re), blup=re[,1], fitted_los=exp(b0+re[,1])-1,
                     stringsAsFactors=FALSE, row.names=NULL)
  
  ll <- as.numeric(logLik(lmm)); ng <- summary(lmm)$ngrps
  fc_row <- data.frame(dataset=nm, k_level=all_k[idx], AIC=AIC(lmm), BIC=BIC(lmm),
                       logLik=ll, deviance=-2*ll, N_obs=nobs(lmm), N_icd3=ng[["icd3"]], N_patient=ng[["patient_id"]],
                       stringsAsFactors=FALSE)
  
  open_png(paste0("Fig_QQ_residuals_",nm,".png")); par(mar=c(5,6,4,2), family="serif")
  qqnorm(s10, main=paste0("Normal QQ-Plot of LMM Residuals: ",nm), xlab="Theoretical Quantiles",
         ylab="Sample Quantiles (log-LOS residuals)", pch=16, cex=0.3, col=adjustcolor("#404040",0.3),
         cex.main=1.3, cex.lab=1.2, cex.axis=1.1, las=1); qqline(s10, col="#C00000", lwd=2.5); dev.off()
  bl <- blup$blup
  open_png(paste0("Fig_QQ_BLUPs_",nm,".png")); par(mar=c(5,6,4,2), family="serif")
  qqnorm(bl, main=paste0("Normal QQ-Plot of ICD-3 Random Intercepts: ",nm), xlab="Theoretical Quantiles",
         ylab="BLUP (random intercept deviation)", pch=19, cex=1.0, col=adjustcolor("#2E75B6",0.6),
         cex.main=1.3, cex.lab=1.2, cex.axis=1.1, las=1); qqline(bl, col="#C00000", lwd=2.5)
  legend("topleft", bty="n", cex=1.0, legend=c(sprintf("%d ICD-3 groups",length(bl)), sprintf("SD(BLUPs)=%.4f",sd(bl)))); dev.off()
  
  # accumulate + checkpoint to disk so an interrupted run resumes here
  variance_components <- rbind(variance_components, vc_row)
  r_squared           <- rbind(r_squared, r2_row)
  prediction_accuracy <- rbind(prediction_accuracy, pa_row)
  fixed_effects       <- rbind(fixed_effects, fe_row)
  fit_comparison      <- rbind(fit_comparison, fc_row)
  normality_results[[nm]] <- norm
  blup_tables[[nm]]       <- blup
  saveRDS(list(vc_row=vc_row, r2_row=r2_row, pa_row=pa_row, fe_row=fe_row,
               fc_row=fc_row, norm=norm, blup=blup), ckpt)
  
  cat(sprintf(" ICC_icd=%.4f intercept=%.4f LOS=%.2f d (%.0fs)\n",
              s2_icd/tot, b0, exp(b0)-1, (proc.time()-t0)[3]))
  rm(lmm, res_all, s5, s10, pred, resid_log, pred_days, resid_days, re); gc(verbose=FALSE)
}

print(variance_components, row.names = FALSE)
write.csv(variance_components, file.path(OUTPUT_DIR, "variance_components.csv"), row.names = FALSE)
write.csv(prediction_accuracy, file.path(OUTPUT_DIR, "prediction_accuracy.csv"), row.names = FALSE)
write.csv(fixed_effects, file.path(OUTPUT_DIR, "fixed_effects.csv"), row.names = FALSE)
print(fit_comparison, row.names = FALSE)
write.csv(fit_comparison, file.path(OUTPUT_DIR, "model_fit_comparison.csv"), row.names = FALSE)

# ============================================================================
# PART 4: REPRODUCIBILITY ASSESSMENT
# ============================================================================
cat("\n####################################################################\n")
cat("#  Reproducibility Assessment                                       #\n")
cat("####################################################################\n\n")
blup_d0 <- blup_tables[["D0"]]
lin_ccc <- function(x,y){ mx<-mean(x);my<-mean(y);vx<-var(x);vy<-var(y);cxy<-cov(x,y); (2*cxy)/(vx+vy+(mx-my)^2) }
boot_ci <- function(x,y,stat_fun,B=N_BOOT){ n<-length(x); if(n<3) return(c(NA_real_,NA_real_))
bs<-numeric(B); for(b in seq_len(B)){ id<-sample.int(n,n,replace=TRUE); bs[b]<-stat_fun(x[id],y[id]) }
unname(quantile(bs, c(.025,.975), na.rm=TRUE)) }

cat("--- BLUP concordance (Spearman + Lin CCC, 95% bootstrap CI) ---\n")
concordance <- data.frame()
for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]; bdk <- blup_tables[[lvl]]
  if (is.null(blup_d0)||is.null(bdk)) next
  shared <- intersect(blup_d0$icd3, bdk$icd3); if (length(shared)<3){cat(sprintf("  k=%d: <3 shared\n",k)); next}
  v0 <- blup_d0$blup[match(shared,blup_d0$icd3)]; vk <- bdk$blup[match(shared,bdk$icd3)]
  rt_sp <- suppressWarnings(cor.test(v0,vk,method="spearman",exact=FALSE)); ccc <- lin_ccc(v0,vk)
  sp_ci <- boot_ci(v0,vk,function(a,b) suppressWarnings(cor(a,b,method="spearman")))
  ccc_ci <- boot_ci(v0,vk,lin_ccc)
  concordance <- rbind(concordance, data.frame(k_level=k, n_shared=length(shared),
                                               n_D0_only=length(setdiff(blup_d0$icd3,bdk$icd3)), n_Dk_only=length(setdiff(bdk$icd3,blup_d0$icd3)),
                                               spearman_rho=unname(rt_sp$estimate), rho_p=rt_sp$p.value, spearman_lo=sp_ci[1], spearman_hi=sp_ci[2],
                                               pearson_r=cor(v0,vk), lin_ccc=ccc, lin_ccc_lo=ccc_ci[1], lin_ccc_hi=ccc_ci[2], stringsAsFactors=FALSE))
  cat(sprintf("  k=%d: %d shared | rho=%.4f [%.4f,%.4f] | CCC=%.4f [%.4f,%.4f]\n",
              k, length(shared), rt_sp$estimate, sp_ci[1], sp_ci[2], ccc, ccc_ci[1], ccc_ci[2]))
}
write.csv(concordance, file.path(OUTPUT_DIR, "BLUP_concordance.csv"), row.names = FALSE)

cat("\n--- BLUP attenuation (ratio Dk/D0) ---\n")
attenuation_data <- list()
for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]; bdk <- blup_tables[[lvl]]
  if (is.null(blup_d0)||is.null(bdk)) next
  shared <- intersect(blup_d0$icd3, bdk$icd3); if (!length(shared)) next
  v0 <- blup_d0$blup[match(shared,blup_d0$icd3)]; vk <- bdk$blup[match(shared,bdk$icd3)]
  nz <- abs(v0)>1e-4
  attenuation_data[[lvl]] <- data.frame(icd3=shared[nz], blup_d0=v0[nz], blup_dk=vk[nz],
                                        ratio=vk[nz]/v0[nz], k_level=k, stringsAsFactors=FALSE)
  r <- attenuation_data[[lvl]]$ratio
  cat(sprintf("  k=%d: median ratio=%.4f IQR[%.4f,%.4f] n=%d\n", k, median(r), quantile(r,.25), quantile(r,.75), length(r)))
}
if (length(attenuation_data)) write.csv(do.call(rbind,attenuation_data), file.path(OUTPUT_DIR,"BLUP_attenuation.csv"), row.names=FALSE)

cat("\n--- Extreme individual distortion ---\n")
distortion_list <- list()
for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]; ad <- attenuation_data[[lvl]]; if (is.null(ad)) next
  ad$sign_reversal <- sign(ad$blup_dk)!=sign(ad$blup_d0); ad$pct_of_original <- 100*ad$ratio
  distortion_list[[lvl]] <- ad; n_rev <- sum(ad$sign_reversal)
  cat(sprintf("  k=%d: sign reversals=%d/%d (%.1f%%) | min=%.0f%% max=%.0f%% of original\n",
              k, n_rev, nrow(ad), 100*n_rev/nrow(ad), min(ad$pct_of_original), max(ad$pct_of_original)))
}
distortion_all <- if (length(distortion_list)) do.call(rbind, distortion_list) else NULL
if (!is.null(distortion_all)) write.csv(distortion_all, file.path(OUTPUT_DIR,"BLUP_extreme_distortion.csv"), row.names=FALSE)

cat("\n--- Variance inflation interrogation (between-diagnosis) ---\n")
s2_icd_d0 <- variance_components$sigma2_icd3[variance_components$dataset=="D0"]
var_inflation <- data.frame()
for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]
  s2_k <- variance_components$sigma2_icd3[variance_components$dataset==lvl]; if (!length(s2_k)) next
  n0 <- model_data[["D0"]]$n_cats; nk <- model_data[[lvl]]$n_cats
  ms0 <- median(as.numeric(table(model_data[["D0"]]$data$icd3)))
  msk <- median(as.numeric(table(model_data[[lvl]]$data$icd3)))
  var_inflation <- rbind(var_inflation, data.frame(k_level=k, sigma2_icd3_D0=s2_icd_d0,
                                                   sigma2_icd3_Dk=s2_k, abs_change=s2_k-s2_icd_d0, rel_change_pct=100*(s2_k-s2_icd_d0)/s2_icd_d0,
                                                   n_icd3_D0=n0, n_icd3_Dk=nk, icd3_dropped=n0-nk, median_group_size_D0=ms0,
                                                   median_group_size_Dk=msk, stringsAsFactors=FALSE))
  cat(sprintf("  k=%d: s2_icd %.4f -> %.4f (%+.1f%%) | groups %d -> %d (dropped %d) | med size %.0f -> %.0f\n",
              k, s2_icd_d0, s2_k, 100*(s2_k-s2_icd_d0)/s2_icd_d0, n0, nk, n0-nk, ms0, msk))
}
if (nrow(var_inflation)) write.csv(var_inflation, file.path(OUTPUT_DIR,"variance_inflation_interrogation.csv"), row.names=FALSE)

cat("\n--- Corroboration ---\n")
corroboration <- NULL
if (nrow(concordance)==length(ANON_LEVELS)) {
  corroboration <- data.frame(k_level=concordance$k_level, spearman_rho=concordance$spearman_rho,
                              lin_ccc=concordance$lin_ccc, los_KS_D=summary_metrics$los_KS_D, icd_jaccard=summary_metrics$icd_jaccard,
                              icd_cramer_v=summary_metrics$icd_cramer_v,
                              ICC_icd3=variance_components$ICC_icd3[variance_components$k_level %in% K_VALUES], stringsAsFactors=FALSE)
  print(corroboration, row.names = FALSE)
}

# ---- Reproducibility figures ----------------------------------------------
cat("\n=== Generating reproducibility figures ===\n")
open_png("Fig_BLUP_concordance_scatter.png"); par(mfrow=c(2,2), mar=c(5,5,4,2), family="serif")
for (i in seq_along(ANON_LEVELS)) {
  lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]; ad <- attenuation_data[[lvl]]
  if (is.null(ad)) { plot.new(); next }
  lim <- range(c(ad$blup_d0, ad$blup_dk))*1.1
  plot(ad$blup_d0, ad$blup_dk, pch=19, cex=0.8, col=adjustcolor(COLS[i],0.5),
       xlab=expression("BLUP"[D0]), ylab=bquote("BLUP"[.(lvl)]),
       main=bquote("D0 vs. "~.(lvl)~" (k="~.(k)~")"), xlim=lim, ylim=lim, cex.main=1.2, cex.lab=1.1, las=1)
  abline(0,1,lty=2,col="gray40",lwd=2)
  rho_v <- concordance$spearman_rho[concordance$k_level==k]; ccc_v <- concordance$lin_ccc[concordance$k_level==k]
  legend("topleft", bty="n", cex=1.0, text.font=2, legend=c(sprintf("rho = %.4f",rho_v), sprintf("CCC = %.4f",ccc_v)))
}
plot.new(); dev.off(); cat("  Fig_BLUP_concordance_scatter.png\n")

if (length(attenuation_data)) {
  open_png("Fig_BLUP_attenuation_ratio.png"); par(mar=c(5,6,4,2), family="serif")
  all_att <- do.call(rbind, attenuation_data)
  all_att$k_fac <- factor(all_att$k_level, levels=K_VALUES, labels=paste0("k=",K_VALUES))
  ql <- quantile(all_att$ratio,.01); qu <- quantile(all_att$ratio,.99)
  att_t <- all_att[all_att$ratio>=ql & all_att$ratio<=qu, ]
  boxplot(ratio~k_fac, data=att_t, col=adjustcolor(COLS,.4), border=COLS, lwd=2, outline=FALSE,
          ylab=expression("BLUP"[Dk]/"BLUP"[D0]), xlab="k-Anonymity Threshold",
          main="BLUP Attenuation After Anonymization", cex.main=1.3, cex.lab=1.2, las=1)
  abline(h=1.0, lty=2, col="gray30", lwd=2)
  meds <- tapply(att_t$ratio, att_t$k_fac, median)
  text(seq_along(meds), meds, labels=sprintf("%.3f",meds), pos=4, cex=0.95, font=2, col=COLS)
  dev.off(); cat("  Fig_BLUP_attenuation_ratio.png\n")
}

if (!is.null(distortion_all)) {
  open_png("Fig_extreme_distortion.png"); par(mfrow=c(2,2), mar=c(5,5,4,2), family="serif")
  for (i in seq_along(ANON_LEVELS)) {
    lvl <- ANON_LEVELS[i]; k <- K_VALUES[i]; ad <- distortion_list[[lvl]]; if (is.null(ad)){plot.new(); next}
    o <- order(abs(log(pmax(ad$ratio,1e-3))), decreasing=TRUE)
    top <- ad[o[seq_len(min(12,nrow(ad)))], ]; top <- top[order(top$pct_of_original), ]
    bcol <- ifelse(top$sign_reversal, "#C00000", COLS[i])
    bp <- barplot(top$pct_of_original, horiz=TRUE, names.arg=top$icd3, las=1, col=adjustcolor(bcol,.8),
                  border=NA, xlab="% of original BLUP", main=bquote("Extreme distortions, k="~.(k)),
                  cex.main=1.15, cex.names=0.8, cex.axis=0.9)
    abline(v=100, lty=2, col="gray30", lwd=1.5); abline(v=0, lty=1, col="gray50")
    text(top$pct_of_original, bp, labels=sprintf("%.0f%%",top$pct_of_original), pos=4, cex=0.65, font=2)
  }
  plot.new(); legend("center", bty="n", cex=1.1, legend=c("sign reversal","magnitude change only"),
                     fill=c("#C00000",COLS[1]), border=NA, title="Per-diagnosis BLUP distortion")
  dev.off(); cat("  Fig_extreme_distortion.png\n")
}

if (!is.null(corroboration)) {
  open_png("Fig_corroboration.png"); par(mfrow=c(2,2), mar=c(5,5,4,2), family="serif")
  rho_v <- concordance$spearman_rho
  plot(summary_metrics$los_KS_D, rho_v, pch=19, cex=3, col=COLS, xlab="KS D (los_days)",
       ylab=expression(rho[S]~"(BLUPs)"), main="(A) Distributional Shift vs. Concordance", cex.main=1.1, cex.lab=1.1, las=1)
  text(summary_metrics$los_KS_D, rho_v, labels=paste0("k=",K_VALUES), pos=1, cex=1.0, font=2)
  plot(summary_metrics$icd_jaccard, rho_v, pch=19, cex=3, col=COLS, xlab="ICD Jaccard",
       ylab=expression(rho[S]~"(BLUPs)"), main="(B) Label Preservation vs. Concordance", cex.main=1.1, cex.lab=1.1, las=1)
  text(summary_metrics$icd_jaccard, rho_v, labels=paste0("k=",K_VALUES), pos=1, cex=1.0, font=2)
  plot(summary_metrics$icd_cramer_v, rho_v, pch=19, cex=3, col=COLS, xlab="Cramer's V",
       ylab=expression(rho[S]~"(BLUPs)"), main="(C) Distortion vs. Concordance", cex.main=1.1, cex.lab=1.1, las=1)
  text(summary_metrics$icd_cramer_v, rho_v, labels=paste0("k=",K_VALUES), pos=1, cex=1.0, font=2)
  plot.new(); dev.off(); cat("  Fig_corroboration.png\n")
}

ok_idx <- which(!is.na(variance_components$ICC_icd3))
if (length(ok_idx)>=2) {
  vc_ok <- variance_components[ok_idx,]; fc_ok <- fit_comparison[match(vc_ok$dataset, fit_comparison$dataset),]
  col_ok <- c(col_orig,COLS)[ok_idx]
  open_png("Fig_ICC_AIC_trajectory.png"); par(mfrow=c(1,2), mar=c(5,6,4,2), family="serif")
  plot(vc_ok$k_level, vc_ok$ICC_icd3, type="b", pch=19, cex=2.5, lwd=3, col=col_ok,
       xlab="k-Anonymity (0=Original)", ylab="ICC (ICD-3 variance share)", main="(A) Variance Explained by ICD-3",
       xaxt="n", las=1, cex.main=1.2, cex.lab=1.1, ylim=c(min(vc_ok$ICC_icd3)*0.9, max(vc_ok$ICC_icd3)*1.1))
  axis(1, at=vc_ok$k_level, labels=vc_ok$dataset)
  text(vc_ok$k_level, vc_ok$ICC_icd3, labels=sprintf("%.4f",vc_ok$ICC_icd3), pos=3, cex=0.9, font=2, offset=1.0)
  if (!all(is.na(fc_ok$AIC))) {
    plot(fc_ok$k_level, fc_ok$AIC, type="b", pch=19, cex=2.5, lwd=3, col=col_ok,
         xlab="k-Anonymity (0=Original)", ylab="AIC (REML)", main="(B) Model AIC", xaxt="n", las=1, cex.main=1.2, cex.lab=1.1)
    axis(1, at=fc_ok$k_level, labels=fc_ok$dataset)
    text(fc_ok$k_level, fc_ok$AIC, labels=formatC(round(fc_ok$AIC), format="d", big.mark=","), pos=3, cex=0.9, font=2, offset=1.0)
  } else plot.new()
  dev.off(); cat("  Fig_ICC_AIC_trajectory.png\n")
}

if (nrow(variance_components)>=2 && !is.null(blup_d0)) {
  open_png("Fig_variance_architecture.png"); par(mfrow=c(2,2), family="serif")
  par(mar=c(5,6,4,2))
  vm <- rbind("ICD-3"=variance_components$sigma2_icd3, "Patient"=variance_components$sigma2_patient,
              "Residual"=variance_components$sigma2_residual); colnames(vm) <- variance_components$dataset
  barplot(vm, col=c("#2E75B6","#70AD47","#B0B0B0"), border=NA, ylab=expression("Variance ("*sigma^2*")"),
          main="(A) Variance Decomposition", cex.main=1.3, cex.lab=1.2, las=1)
  legend("topright", legend=c("ICD-3","Patient","Residual"), fill=c("#2E75B6","#70AD47","#B0B0B0"), border=NA, cex=1.0, bg="white")
  par(mar=c(5,5,4,2)); av0<-numeric(0); avk<-numeric(0); acol<-character(0)
  for (i in seq_along(ANON_LEVELS)) { lvl<-ANON_LEVELS[i]; bdk<-blup_tables[[lvl]]; if(is.null(bdk)) next
  sh<-intersect(blup_d0$icd3,bdk$icd3); if(length(sh)<3) next
  av0<-c(av0,blup_d0$blup[match(sh,blup_d0$icd3)]); avk<-c(avk,bdk$blup[match(sh,bdk$icd3)]); acol<-c(acol,rep(COLS[i],length(sh))) }
  if (length(av0)) { lim<-range(c(av0,avk))*1.1
  plot(av0,avk,pch=16,cex=0.7,col=adjustcolor(acol,.4), xlab=expression("BLUP"[D0]), ylab=expression("BLUP"[Dk]),
       main="(B) BLUP Fidelity", xlim=lim, ylim=lim, cex.main=1.3, cex.lab=1.1, las=1)
  abline(0,1,lty=2,col="gray30",lwd=2.5); legend("topleft", legend=paste0("k=",K_VALUES), pch=16, col=COLS, cex=1.0, bg="white")
  } else plot.new()
  par(mar=c(5,5,4,2))
  dens_list <- lapply(all_names, function(nm) if(!is.null(blup_tables[[nm]]) && nrow(blup_tables[[nm]])>=3) density(blup_tables[[nm]]$blup) else NULL)
  names(dens_list) <- all_names
  xr <- range(unlist(lapply(dens_list, function(d) if(!is.null(d)) range(d$x))))
  ym <- max(unlist(lapply(dens_list, function(d) if(!is.null(d)) max(d$y))))
  plot(NULL, xlim=xr, ylim=c(0,ym*1.1), xlab="Random Intercept (BLUP)", ylab="Density",
       main="(C) Distribution of ICD-3 BLUPs", cex.main=1.3, cex.lab=1.1, las=1)
  pc <- c(col_orig,COLS)
  for (d in seq_along(all_names)) if(!is.null(dens_list[[all_names[d]]])) lines(dens_list[[all_names[d]]], lwd=2.5, col=pc[d], lty=d)
  abline(v=0, lty=3, col="gray50"); legend("topright", legend=c("D0",paste0("k=",K_VALUES)), col=pc, lwd=2.5, lty=1:4, cex=0.95, bg="white")
  par(mar=c(5,6,4,6)); icc_v <- variance_components$ICC_icd3; gm_v <- fixed_effects$intercept_exp; ka <- variance_components$k_level
  plot(ka, icc_v, type="b", pch=19, cex=2.5, lwd=3, col="#2E75B6", xlab="k-Anonymity (0=Original)",
       ylab="ICC (ICD-3 share)", main="(D) ICC & Grand Mean LOS", xaxt="n", las=1,
       ylim=c(min(icc_v)*0.9, max(icc_v)*1.1), cex.main=1.3, cex.lab=1.1)
  axis(1, at=ka, labels=variance_components$dataset)
  text(ka, icc_v, labels=sprintf("%.4f",icc_v), pos=3, cex=0.85, font=2, col="#2E75B6", offset=0.8)
  par(new=TRUE)
  plot(ka, gm_v, type="b", pch=17, cex=2.5, lwd=3, col="#C00000", axes=FALSE, xlab="", ylab="", ylim=c(min(gm_v)*0.95, max(gm_v)*1.05))
  axis(4, las=1, col="#C00000", col.axis="#C00000"); mtext("Grand Mean LOS (days)", side=4, line=3.5, cex=1.1, col="#C00000")
  text(ka, gm_v, labels=sprintf("%.2f",gm_v), pos=1, cex=0.85, font=2, col="#C00000", offset=0.8)
  dev.off(); cat("  Fig_variance_architecture.png\n")
}

# ============================================================================
# PART 5: SUMMARY TABLE, JSON, RDS
# ============================================================================
cat("\n--- Exporting LMM summary table ---\n")
lmm_summary_table <- data.frame()
for (idx in seq_along(all_names)) {
  nm <- all_names[idx]; k <- all_k[idx]
  fe <- fixed_effects[fixed_effects$dataset==nm,]; vc <- variance_components[variance_components$dataset==nm,]
  r2 <- r_squared[r_squared$dataset==nm,]; fc <- fit_comparison[fit_comparison$dataset==nm,]
  pa <- prediction_accuracy[prediction_accuracy$dataset==nm,]; co <- concordance[concordance$k_level==k,]
  at <- attenuation_data[[nm]]; one <- function(d,col) if(nrow(d)==1) d[[col]] else NA
  lmm_summary_table <- rbind(lmm_summary_table, data.frame(
    dataset=nm, k_level=k, N_obs=one(fc,"N_obs"), N_icd3=one(fc,"N_icd3"), N_patient=one(fc,"N_patient"),
    intercept=one(fe,"intercept"), intercept_SE=one(fe,"intercept_se"),
    intercept_CI95=if(nrow(fe)==1) sprintf("[%.4f, %.4f]", fe$intercept_lo, fe$intercept_hi) else NA,
    grand_mean_LOS=one(fe,"intercept_exp"), sigma2_icd3=one(vc,"sigma2_icd3"),
    sigma2_patient=one(vc,"sigma2_patient"), sigma2_residual=one(vc,"sigma2_residual"),
    ICC_icd3=one(vc,"ICC_icd3"), ICC_patient=one(vc,"ICC_patient"),
    R2_marg=one(r2,"R2_marg"), R2_cond=one(r2,"R2_cond"), AIC=one(fc,"AIC"), BIC=one(fc,"BIC"), logLik=one(fc,"logLik"),
    RMSE_log=one(pa,"RMSE_log"), MAE_log=one(pa,"MAE_log"), RMSE_days=one(pa,"RMSE_days"), MAE_days=one(pa,"MAE_days"),
    BLUP_spearman=if(nrow(co)==1) co$spearman_rho else NA,
    BLUP_spearman_CI=if(nrow(co)==1) sprintf("[%.4f, %.4f]", co$spearman_lo, co$spearman_hi) else NA,
    BLUP_lin_ccc=if(nrow(co)==1) co$lin_ccc else NA,
    BLUP_lin_ccc_CI=if(nrow(co)==1) sprintf("[%.4f, %.4f]", co$lin_ccc_lo, co$lin_ccc_hi) else NA,
    BLUP_atten_median=if(!is.null(at)) median(at$ratio) else NA, stringsAsFactors=FALSE))
}
print(lmm_summary_table[, c("dataset","N_obs","ICC_icd3","ICC_patient","grand_mean_LOS","BLUP_spearman","BLUP_lin_ccc")], row.names=FALSE)
write.csv(lmm_summary_table, file.path(OUTPUT_DIR, "LMM_summary.csv"), row.names = FALSE)

cat("\n--- Writing test_results_report.json ---\n")
js <- function(x){ if(is.null(x)||(length(x)==1&&is.na(x))) return("null")
  if(is.logical(x)) return(tolower(as.character(x)))
  if(is.numeric(x)) return(ifelse(is.finite(x), formatC(x, format="g", digits=8), "null"))
  paste0("\"", gsub("\"","\\\\\"",as.character(x)), "\"") }
kv <- function(k,v) paste0("\"",k,"\": ",v); arr <- function(v) paste0("[", paste(v,collapse=", "), "]")
tests <- c(
  kv("schema_columns_ok", js(all(c("patient_id","encounter_id","admission_year","main_diagnosis_icd","los_days") %in% colnames(df_D0)))),
  kv("star_token_mapped_to_NA", js(
    !any(df_D5$main_diagnosis_icd == "*", na.rm=TRUE) &&    # no literal "*" survives the read
      sum(is.na(df_D5$main_diagnosis_icd)) > 0)),             # suppression now present as NA
  kv("ks_z_no_nan", js(!any(is.nan(summary_metrics$los_KS_z)))),
  kv("ks_z_all_finite", js(all(is.finite(summary_metrics$los_KS_z)))),
  kv("ks_z_reasons", arr(sapply(summary_metrics$los_KS_z_reason, js))),
  kv("ks_z_values", arr(sapply(summary_metrics$los_KS_z, js))),
  kv("cell_suppression_detected", js(
    all(removal_documentation$cells_icd_suppressed > 0) &&
      all(removal_documentation$encounters_removed == 0))),   # masking, not row deletion
  kv("lmm_has_patient_term", js(all(c("sigma2_patient") %in% colnames(variance_components)) && all(variance_components$sigma2_patient >= 0))),
  kv("lmm_has_year_fixed", js(nrow(r_squared) == length(all_names) && all(is.finite(r_squared$R2_marg)))),
  kv("concordance_has_CIs", js(all(c("spearman_lo","spearman_hi","lin_ccc_lo","lin_ccc_hi") %in% colnames(concordance)))),
  kv("n_models_converged", js(nrow(variance_components))),
  kv("between_diag_variance_D0", js(s2_icd_d0)),
  kv("between_diag_variance_by_k", arr(sapply(ANON_LEVELS, function(l) js(variance_components$sigma2_icd3[variance_components$dataset==l])))),
  kv("sign_reversal_counts_by_k", arr(sapply(ANON_LEVELS, function(l) js(if(!is.null(distortion_list[[l]])) sum(distortion_list[[l]]$sign_reversal) else NA))))
)
json <- paste0("{\n  ", paste(tests, collapse=",\n  "), "\n}\n")
writeLines(json, file.path(OUTPUT_DIR, "test_results_report.json")); cat(json)

study_results <- list(summary_metrics=summary_metrics, removal_documentation=removal_documentation,
                      suppression_by_year=supp_by_year, cohort_summary=cohort_summary, na_documentation=na_documentation,
                      variance_components=variance_components, r_squared=r_squared, fixed_effects=fixed_effects,
                      blup_tables=blup_tables, concordance=concordance, attenuation=attenuation_data, distortion=distortion_list,
                      variance_inflation=var_inflation, prediction_accuracy=prediction_accuracy, fit_comparison=fit_comparison,
                      normality_results=normality_results, lmm_summary_table=lmm_summary_table, corroboration=corroboration,
                      metadata=list(date=Sys.time(), R_version=R.version.string,
                                    model="log(los_days+1) ~ year + (1|icd3) + (1|patient_id)"))
saveRDS(study_results, file.path(OUTPUT_DIR, "study_results_full.rds"))
cat("\n  study_results_full.rds saved\n")
cat("\n=== Analysis complete. Outputs in:", OUTPUT_DIR, "===\n")