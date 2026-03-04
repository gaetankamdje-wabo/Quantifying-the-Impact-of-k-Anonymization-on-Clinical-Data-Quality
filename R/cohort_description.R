# ──────────────────────────────────────────────────────────────────────────────
# Descriptive cohort analysis of hospital encounter data
# Outputs:summary table
#
# Author : Gaetan Kamdje Wabo
# Context: Anonymization impact study – cohort characterization
# Date   : 04.03.2026

# ── 0. Dependencies (minimal set) ──────────────────────────────────────────────
library(data.table)   # Fast data manipulation
library(ggplot2)      # Publication-quality plots
library(flextable)    # Word-compatible tables
library(officer)      # DOCX export

# ── 1. Paths ───────────────────────────────────────────────────────────────────
# ------ CONFIGURE THESE PATHS -----------------------------------------------
# Point to the directory containing the source CSV files.
# Default assumes the working directory is the repository root.
path_input  <- file.path("data", "original data.csv")
path_output <- file.path("output")
# ----------------------------------------------------------------------------

# Create output directory if it does not exist
if (!dir.exists(path_output)) dir.create(path_output, recursive = TRUE)

# ── 2. Data Import ─────────────────────────────────────────────────────────────
dt <- fread(path_input, sep = ",", header = TRUE, na.strings = c("", "NA"))

# Parse dates; extract study year from admission_timestamp
dt[, admission_date := as.Date(admission_timestamp)]
dt[, discharge_date := as.Date(discharge_timestamp)]
dt[, study_year := year(admission_date)]

# Replace NA in los_days with 1 (same-day encounters)
dt[is.na(los_days), los_days := 1]

# Clean ICD column name for convenience
setnames(dt, "icd_code(Hauptdia)", "icd_hauptdia", skip_absent = TRUE)

cat("Data loaded:", nrow(dt), "rows |", uniqueN(dt$encounter_id), "unique encounters\n")

# ── 3. Study Period Overview ───────────────────────────────────────────────────
year_range <- range(dt$study_year, na.rm = TRUE)
cat("Study period:", year_range[1], "-", year_range[2], "\n")

# 4. ENCOUNTER ANALYSIS

# 4a. Total encounters
n_total <- uniqueN(dt$encounter_id)

# 4b. Encounters per year
enc_per_year <- dt[, .(n_encounters = uniqueN(encounter_id)), by = study_year]
setorder(enc_per_year, study_year)

enc_summary <- data.table(
  Metric = c("Total encounters", "Study period",
             "Encounters/year – Min", "Encounters/year – Max",
             "Encounters/year – Mean", "Encounters/year – Median"),
  Value  = c(
    format(n_total, big.mark = ","),
    paste0(year_range[1], " – ", year_range[2]),
    format(min(enc_per_year$n_encounters),  big.mark = ","),
    format(max(enc_per_year$n_encounters),  big.mark = ","),
    format(round(mean(enc_per_year$n_encounters), 1), big.mark = ","),
    format(median(enc_per_year$n_encounters), big.mark = ",")
  )
)

# 5. TOP-5 DIAGNOSIS PREVALENCE (YEARS-BASED)

# 5a. Identify top 5 ICD codes by overall frequency
icd_freq <- dt[!is.na(icd_hauptdia), .N, by = icd_hauptdia]
setorder(icd_freq, -N)
top5_codes <- icd_freq[1:5, icd_hauptdia]

# 5b. Yearly prevalence for top 5 (proportion of encounters per year)
yearly_totals <- dt[, .(total = uniqueN(encounter_id)), by = study_year]

top5_yearly <- dt[icd_hauptdia %in% top5_codes,
                  .(count = uniqueN(encounter_id)),
                  by = .(study_year, icd_hauptdia)]

top5_yearly <- merge(top5_yearly, yearly_totals, by = "study_year")
top5_yearly[, prevalence := round(count / total * 100, 2)]

# 5c. Mean prevalence across years per diagnosis
top5_mean_prev <- top5_yearly[, .(
  mean_prevalence_pct = round(mean(prevalence), 2),
  sd_prevalence_pct   = round(sd(prevalence), 2),
  total_count         = sum(count)
), by = icd_hauptdia]

setorder(top5_mean_prev, -total_count)

# Add overall proportion column
top5_mean_prev[, overall_pct := round(total_count / n_total * 100, 2)]

# 6. LENGTH OF STAY (LOS) ANALYSIS

# 6a. Overall LOS statistics (NA already replaced with 1)
los_overall <- data.table(
  Period  = "Entire study",
  Min     = min(dt$los_days, na.rm = TRUE),
  Max     = max(dt$los_days, na.rm = TRUE),
  Mean    = round(mean(dt$los_days, na.rm = TRUE), 2),
  Median  = median(dt$los_days, na.rm = TRUE),
  SD      = round(sd(dt$los_days, na.rm = TRUE), 2),
  IQR     = IQR(dt$los_days, na.rm = TRUE),
  N       = nrow(dt)
)

# 6b. LOS per year
los_per_year <- dt[, .(
  Min    = min(los_days, na.rm = TRUE),
  Max    = max(los_days, na.rm = TRUE),
  Mean   = round(mean(los_days, na.rm = TRUE), 2),
  Median = as.double(median(los_days, na.rm = TRUE)),
  SD     = round(sd(los_days, na.rm = TRUE), 2),
  IQR    = IQR(los_days, na.rm = TRUE),
  N      = .N
), by = study_year]

setorder(los_per_year, study_year)
los_per_year[, Period := as.character(study_year)]
los_per_year[, study_year := NULL]

# 6c. Combined LOS table (overall + per year)
los_combined <- rbindlist(list(los_overall, los_per_year), use.names = TRUE)

# 7. COMBINED SUMMARY TABLE (for the paper)

# Build a single cohort characteristics table
cohort_table <- data.table(
  Characteristic = c(
    "Study period",
    "Total hospital encounters, n",
    "",
    "Encounters per year",
    "  Min",
    "  Max",
    "  Mean (SD)",
    "  Median",
    "",
    "Length of stay (days) – Overall",
    "  Min",
    "  Max",
    "  Mean (SD)",
    "  Median (IQR)",
    "",
    "Top 5 principal diagnoses (ICD-10-GM)"
  ),
  Value = c(
    paste0(year_range[1], " – ", year_range[2]),
    format(n_total, big.mark = ","),
    "",
    "",
    format(min(enc_per_year$n_encounters), big.mark = ","),
    format(max(enc_per_year$n_encounters), big.mark = ","),
    paste0(format(round(mean(enc_per_year$n_encounters), 1), big.mark = ","),
           " (", format(round(sd(enc_per_year$n_encounters), 1), big.mark = ","), ")"),
    format(median(enc_per_year$n_encounters), big.mark = ","),
    "",
    "",
    as.character(los_overall$Min),
    format(los_overall$Max, big.mark = ","),
    paste0(los_overall$Mean, " (", los_overall$SD, ")"),
    paste0(los_overall$Median, " (", los_overall$IQR, ")"),
    "",
    ""
  )
)

# Append top 5 diagnoses rows
for (i in seq_len(nrow(top5_mean_prev))) {
  cohort_table <- rbindlist(list(
    cohort_table,
    data.table(
      Characteristic = paste0("  ", top5_mean_prev$icd_hauptdia[i]),
      Value = paste0(
        format(top5_mean_prev$total_count[i], big.mark = ","),
        " (", top5_mean_prev$overall_pct[i], "%) – ",
        "Mean annual prevalence: ",
        top5_mean_prev$mean_prevalence_pct[i], "% (",
        "\u00B1", top5_mean_prev$sd_prevalence_pct[i], "%)"
      )
    )
  ), use.names = TRUE)
}

# 8. TABLE EXPORT (DOCX)

# ── 8a. Cohort characteristics table ──
ft_cohort <- flextable(cohort_table) |>
  set_header_labels(Characteristic = "Characteristic", Value = "Value") |>
  bold(i = c(1, 2, 4, 10, 16), part = "body") |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Times New Roman", part = "all") |>
  autofit() |>
  set_caption("Table 1. Cohort Characteristics of Hospital Encounter Data")

doc <- read_docx() |>
  body_add_par("Table 1. Cohort Characteristics", style = "heading 1") |>
  body_add_flextable(ft_cohort) |>
  body_add_par("") |>
  body_add_par("") |>
  body_add_par("Table 2. Length of Stay by Year", style = "heading 1")

# ── 8b. LOS per year table ──
ft_los <- flextable(los_combined) |>
  set_header_labels(
    Period = "Period", Min = "Min", Max = "Max",
    Mean = "Mean", Median = "Median", SD = "SD", IQR = "IQR", N = "N"
  ) |>
  bold(i = 1, part = "body") |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Times New Roman", part = "all") |>
  autofit() |>
  set_caption("Table 2. Length of Stay Distribution – Overall and by Year")

doc <- doc |>
  body_add_flextable(ft_los) |>
  body_add_par("") |>
  body_add_par("") |>
  body_add_par("Table 3. Top-5 Diagnosis Yearly Prevalence", style = "heading 1")

# ── 8c. Top-5 prevalence per year (wide format) ──
top5_wide <- dcast(top5_yearly, study_year ~ icd_hauptdia,
                   value.var = "prevalence", fill = 0)
setnames(top5_wide, "study_year", "Year")

ft_prev <- flextable(top5_wide) |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Times New Roman", part = "all") |>
  autofit() |>
  set_caption("Table 3. Annual Prevalence (%) of Top 5 Principal Diagnoses")

doc <- doc |>
  body_add_flextable(ft_prev)

print(doc, target = file.path(path_output, "Fig16_cohort_summary_tables.docx"))
cat("Tables exported to DOCX.\n")

# 9. FIGURES

# Consistent theme for all figures
theme_pub <- theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    axis.title       = element_text(face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# ── 9a. Encounters per year ──
p1 <- ggplot(enc_per_year, aes(x = factor(study_year), y = n_encounters)) +
  geom_col(fill = "#2C5F8A", width = 0.7) +
  geom_text(aes(label = format(n_encounters, big.mark = ",")),
            vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = scales::comma_format()) +
  labs(title    = "Hospital Encounters per Year",
       subtitle = paste0("Study period: ", year_range[1], "–", year_range[2],
                         " | Total: N = ", format(n_total, big.mark = ",")),
       x = "Year", y = "Number of Encounters") +
  theme_pub

ggsave(file.path(path_output, "Fig16_encounters_per_year.png"),
       p1, width = 10, height = 5.5, dpi = 300)

# ── 9b. Top-5 diagnosis prevalence trends ──
p2 <- ggplot(top5_yearly, aes(x = study_year, y = prevalence,
                              color = icd_hauptdia, group = icd_hauptdia)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(year_range[1], year_range[2], by = 1)) +
  scale_color_brewer(palette = "Set1") +
  labs(title    = "Annual Prevalence of Top 5 Principal Diagnoses",
       subtitle = paste0("Proportion of total encounters per year (%) | ",
                         year_range[1], "–", year_range[2]),
       x = "Year", y = "Prevalence (%)", color = "ICD-10-GM") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(path_output, "Fig16_top5_diagnosis_prevalence.png"),
       p2, width = 10, height = 6, dpi = 300)

# ── 9c. LOS distribution (boxplot per year) ──
# Cap outliers visually at 99th percentile for readability
los_cap <- quantile(dt$los_days, 0.99, na.rm = TRUE)

p3 <- ggplot(dt, aes(x = factor(study_year), y = los_days)) +
  geom_boxplot(fill = "#5B9BD5", alpha = 0.7, outlier.size = 0.5,
               outlier.alpha = 0.3) +
  coord_cartesian(ylim = c(0, los_cap * 1.1)) +
  labs(title    = "Length of Stay Distribution by Year",
       subtitle = paste0("Overall median: ", los_overall$Median,
                         " days | Mean: ", los_overall$Mean,
                         " days (capped at 99th percentile for display)"),
       x = "Year", y = "Length of Stay (days)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(path_output, "Fig16_los_distribution_by_year.png"),
       p3, width = 10, height = 6, dpi = 300)

# ── 9d. LOS mean trend per year ──
p4 <- ggplot(los_per_year, aes(x = as.numeric(Period), y = Mean)) +
  geom_line(color = "#2C5F8A", linewidth = 1) +
  geom_point(color = "#2C5F8A", size = 2.5) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD),
              alpha = 0.15, fill = "#2C5F8A") +
  scale_x_continuous(breaks = seq(year_range[1], year_range[2], by = 1)) +
  labs(title    = "Mean Length of Stay Over Time",
       subtitle = "Shaded area = \u00B11 SD",
       x = "Year", y = "Mean LOS (days)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(path_output, "Fig16_los_mean_trend.png"),
       p4, width = 10, height = 5.5, dpi = 300)

# 10. SUMMARY
cat("\n══════════════════════════════════════════════════════\n")
cat("  COHORT ANALYSIS COMPLETE\n")
cat("══════════════════════════════════════════════════════\n")
cat("  Total encounters :", format(n_total, big.mark = ","), "\n")
cat("  Study period     :", year_range[1], "–", year_range[2], "\n")
cat("  Enc/year range   :", min(enc_per_year$n_encounters), "–",
    max(enc_per_year$n_encounters), "\n")
cat("  LOS overall      : Median =", los_overall$Median,
    "| Mean =", los_overall$Mean, "\n")
cat("  Top 5 diagnoses  :", paste(top5_codes, collapse = ", "), "\n")
cat("══════════════════════════════════════════════════════\n")
cat("\nOutputs saved to:\n ", path_output, "\n")
cat("  - Fig16_cohort_summary_tables.docx\n")
cat("  - Fig16_encounters_per_year.png\n")
cat("  - Fig16_top5_diagnosis_prevalence.png\n")
cat("  - Fig16_los_distribution_by_year.png\n")
cat("  - Fig16_los_mean_trend.png\n")
# END
