# ──────────────────────────────────────────────────────────────────────────────
# Revised Descriptive cohort analysis of hospital encounter data
# Output: Word document 
#


# ── 0. Dependencies 
library(data.table)
library(flextable)
library(officer)


# ── 1. Paths 
path_input <- "A:/xxx" #add local path for input files

path_output <- "A:/xxx" #add local path for output files

if (!dir.exists(path_output)) {
  dir.create(path_output, recursive = TRUE)
}

docx_output <- file.path(path_output, "Table1_cohort_characteristics.docx")


# ── 2. Study-period (just for verification purpose)
study_period_label <- "2010-01-01 – 2024-09-30"


# ── 3. Data import
dt <- fread(
  path_input,
  sep = ",",
  header = TRUE,
  na.strings = c("", "NA", "NaN", "NULL")
)

# Ensure expected data structure
required_cols <- c(
  "patient_id",
  "encounter_id",
  "admission_year",
  "main_diagnosis_icd",
  "los_days"
)

missing_cols <- setdiff(required_cols, names(dt))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "The following required columns are missing: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# Keep only required variables
dt <- dt[, ..required_cols]

# Basic type harmonization
dt[, patient_id := as.character(patient_id)]
dt[, encounter_id := as.character(encounter_id)]
dt[, admission_year := as.integer(admission_year)]
dt[, main_diagnosis_icd := as.character(main_diagnosis_icd)]
dt[, los_days := as.numeric(los_days)]

# Restrict to the intended study years.
# The exact upper bound 2024-09-30 was already applied during extraction,
# because only admission_year is available in this analysis dataset.
dt <- dt[
  admission_year >= 2010 &
    admission_year <= 2024
]

# Replace missing LOS values with 1 if same-day stays were encoded as missing.
# Remove this line if missing LOS should remain missing.
dt[is.na(los_days), los_days := 1]


# ── 4. Optional duplicate check ───────────────────────────────────────────────
# The dataset is expected to contain one row per hospital encounter.
# If exact duplicate encounter rows exist, they are removed.
dt <- unique(dt)

# Check whether one encounter_id appears more than once with different information
encounter_duplicate_check <- dt[
  ,
  .N,
  by = encounter_id
][N > 1]

if (nrow(encounter_duplicate_check) > 0) {
  warning(
    paste0(
      "There are ",
      nrow(encounter_duplicate_check),
      " encounter_id values appearing more than once. ",
      "Please verify whether the dataset truly contains one row per encounter."
    )
  )
}


# ── 5. Formatting  functions
fmt_n <- function(x) {
  format(
    x,
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

fmt_num <- function(x, digits = 1) {
  format(
    round(x, digits),
    big.mark = ",",
    decimal.mark = ".",
    scientific = FALSE,
    trim = TRUE,
    nsmall = digits
  )
}

fmt_pct <- function(x, digits = 2) {
  paste0(fmt_num(x, digits), "%")
}


# ── 6. Core cohort counts ─────────────────────────────────────────────────────
n_rows <- nrow(dt)

n_unique_patients <- uniqueN(dt$patient_id)
n_unique_encounters <- uniqueN(dt$encounter_id)

patient_encounter_counts <- unique(
  dt[, .(patient_id, encounter_id)]
)[
  ,
  .(n_encounters = uniqueN(encounter_id)),
  by = patient_id
]

n_patients_one_encounter <- patient_encounter_counts[
  n_encounters == 1,
  .N
]

n_patients_more_than_one_encounter <- patient_encounter_counts[
  n_encounters > 1,
  .N
]

n_encounters_from_patients_more_than_one <- patient_encounter_counts[
  n_encounters > 1,
  sum(n_encounters)
]

max_encounters_per_patient <- patient_encounter_counts[
  ,
  max(n_encounters, na.rm = TRUE)
]


# ── 7. Encounters per year ────────────────────────────────────────────────────
enc_per_year <- dt[
  ,
  .(n_encounters = uniqueN(encounter_id)),
  by = admission_year
]

setorder(enc_per_year, admission_year)

enc_year_min <- min(enc_per_year$n_encounters, na.rm = TRUE)
enc_year_max <- max(enc_per_year$n_encounters, na.rm = TRUE)
enc_year_mean <- mean(enc_per_year$n_encounters, na.rm = TRUE)
enc_year_sd <- sd(enc_per_year$n_encounters, na.rm = TRUE)
enc_year_median <- median(enc_per_year$n_encounters, na.rm = TRUE)


# ── 8. Length of stay summary ─────────────────────────────────────────────────
los_min <- min(dt$los_days, na.rm = TRUE)
los_max <- max(dt$los_days, na.rm = TRUE)
los_mean <- mean(dt$los_days, na.rm = TRUE)
los_sd <- sd(dt$los_days, na.rm = TRUE)
los_median <- median(dt$los_days, na.rm = TRUE)
los_iqr <- IQR(dt$los_days, na.rm = TRUE)


# ── 9. Top 5 principal diagnoses ──────────────────────────────────────────────
icd_freq <- dt[
  !is.na(main_diagnosis_icd) &
    main_diagnosis_icd != "",
  .(total_count = uniqueN(encounter_id)),
  by = main_diagnosis_icd
]

setorder(icd_freq, -total_count)

top5_codes <- icd_freq[
  1:min(5, .N),
  main_diagnosis_icd
]

top5_overall <- icd_freq[
  main_diagnosis_icd %in% top5_codes
]

top5_overall[
  ,
  overall_pct := total_count / n_unique_encounters * 100
]

# Annual totals
yearly_totals <- dt[
  ,
  .(year_total_encounters = uniqueN(encounter_id)),
  by = admission_year
]

# Annual counts for top 5 diagnoses
top5_yearly_counts <- dt[
  main_diagnosis_icd %in% top5_codes,
  .(diagnosis_count = uniqueN(encounter_id)),
  by = .(admission_year, main_diagnosis_icd)
]

# Complete all year-diagnosis combinations, including zero-count years
all_years <- sort(unique(dt$admission_year))

top5_grid <- CJ(
  admission_year = all_years,
  main_diagnosis_icd = top5_codes
)

top5_yearly <- merge(
  top5_grid,
  top5_yearly_counts,
  by = c("admission_year", "main_diagnosis_icd"),
  all.x = TRUE
)

top5_yearly[is.na(diagnosis_count), diagnosis_count := 0]

top5_yearly <- merge(
  top5_yearly,
  yearly_totals,
  by = "admission_year",
  all.x = TRUE
)

top5_yearly[
  ,
  annual_prevalence_pct := diagnosis_count / year_total_encounters * 100
]

top5_mean_prev <- top5_yearly[
  ,
  .(
    mean_annual_prevalence_pct = mean(annual_prevalence_pct, na.rm = TRUE),
    sd_annual_prevalence_pct = sd(annual_prevalence_pct, na.rm = TRUE)
  ),
  by = main_diagnosis_icd
]

top5_summary <- merge(
  top5_overall,
  top5_mean_prev,
  by = "main_diagnosis_icd",
  all.x = TRUE
)

setorder(top5_summary, -total_count)


# ── 10. Build Table 1 ─────────────────────────────────────────────────────────
cohort_table <- data.table(
  Characteristic = c(
    "Study period",
    "Total hospital encounters, n",
    "Unique encounters, n",
    "Unique patients, n",
    "",
    "Patient-level encounter structure",
    "  Patients with one encounter, n",
    "  Patients with more than one encounter, n",
    "  Encounters from patients with more than one encounter, n",
    "  Maximum encounters per patient, n",
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
    study_period_label,
    fmt_n(n_rows),
    fmt_n(n_unique_encounters),
    fmt_n(n_unique_patients),
    "",
    "",
    fmt_n(n_patients_one_encounter),
    fmt_n(n_patients_more_than_one_encounter),
    fmt_n(n_encounters_from_patients_more_than_one),
    fmt_n(max_encounters_per_patient),
    "",
    "",
    fmt_n(enc_year_min),
    fmt_n(enc_year_max),
    paste0(fmt_num(enc_year_mean, 1), " (", fmt_num(enc_year_sd, 1), ")"),
    fmt_n(enc_year_median),
    "",
    "",
    fmt_n(los_min),
    fmt_n(los_max),
    paste0(fmt_num(los_mean, 2), " (", fmt_num(los_sd, 2), ")"),
    paste0(fmt_num(los_median, 0), " (", fmt_num(los_iqr, 0), ")"),
    "",
    ""
  )
)

# Append top 5 diagnoses
for (i in seq_len(nrow(top5_summary))) {
  cohort_table <- rbindlist(
    list(
      cohort_table,
      data.table(
        Characteristic = paste0("  ", top5_summary$main_diagnosis_icd[i]),
        Value = paste0(
          fmt_n(top5_summary$total_count[i]),
          " (", fmt_pct(top5_summary$overall_pct[i], 2), ")",
          " – Mean annual prevalence: ",
          fmt_pct(top5_summary$mean_annual_prevalence_pct[i], 2),
          " (±",
          fmt_pct(top5_summary$sd_annual_prevalence_pct[i], 2),
          ")"
        )
      )
    ),
    use.names = TRUE
  )
}


# ── 11. Export Word document ──────────────────────────────────────────────────
section_rows <- which(
  cohort_table$Characteristic %in% c(
    "Study period",
    "Total hospital encounters, n",
    "Patient-level encounter structure",
    "Encounters per year",
    "Length of stay (days) – Overall",
    "Top 5 principal diagnoses (ICD-10-GM)"
  )
)

ft_cohort <- flextable(cohort_table) |>
  set_header_labels(
    Characteristic = "Characteristic",
    Value = "Value"
  ) |>
  bold(i = section_rows, part = "body") |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Times New Roman", part = "all") |>
  align(j = 1, align = "left", part = "all") |>
  align(j = 2, align = "left", part = "all") |>
  autofit() |>
  set_caption("Table 1. Cohort Characteristics of Hospital Encounter Data")

doc <- read_docx() |>
  body_add_par("Table 1. Cohort Characteristics", style = "heading 1") |>
  body_add_flextable(ft_cohort)

print(doc, target = docx_output)


# ── 12. Console summary ───────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════\n")
cat("  COHORT CHARACTERISTICS EXPORT COMPLETE\n")
cat("══════════════════════════════════════════════════════\n")
cat("  Study period                         :", study_period_label, "\n")
cat("  Total rows / hospital encounters     :", fmt_n(n_rows), "\n")
cat("  Unique encounters                    :", fmt_n(n_unique_encounters), "\n")
cat("  Unique patients                      :", fmt_n(n_unique_patients), "\n")
cat("  Patients with one encounter          :", fmt_n(n_patients_one_encounter), "\n")
cat("  Patients with >1 encounter           :", fmt_n(n_patients_more_than_one_encounter), "\n")
cat("  Output saved to                      :", docx_output, "\n")
cat("══════════════════════════════════════════════════════\n")