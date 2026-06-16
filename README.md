# Quantifying the Impact of k-Anonymization on Clinical Data Quality

**Companion code** for the manuscript:

> Kamdje Wabo, G. *Quantifying the Impact of Anonymization-Induced Clinical Data Quality Loss: A Methodological Case Study Using Primary Diagnosis Codes and Hospital Length of Stay.* JMIR Medical Informatics (under review), 2026.

---

## Overview

This repository contains the complete R analysis pipeline for a methodological case study that measures how **ARX k-anonymity** (k = 5, 10, 15) affects the analytical utility of two of the most frequently used elements in retrospective hospital research: **ICD-10-GM primary diagnosis codes** and **hospital length of stay (LOS)**.

Anonymization in this study uses **cell-level suppression combined with microaggregation**. No rows are deleted. ARX replaces individual quasi-identifier values with a missing token (`*`) until the k-anonymity threshold is met. The analysis quantifies the resulting distributional and inferential distortion and shows that meaningful distortion can occur that the anonymization tool itself does not report.

The pipeline proceeds in three stages.

| Stage | Script | Purpose |
|-------|--------|---------|
| **1. Core function** | `R/quantify_anonymization_impact.R` | Reusable, base-R function. Column-presence checks, numeric distributional analysis (KS D, scaled statistic z = sqrt(n_eff)*D, mean/SD shifts), categorical analysis (Jaccard overlap, Cramer V, chi-square over shared levels), and traffic-light verdicts. |
| **2. Cohort description** | `R/cohort_description.R` | Descriptive statistics and a publication-ready summary table for the original hospital encounter dataset. |
| **3. Main analysis** | `R/final_analysis.R` | Full impact analysis. Runs Stage 1 across k = 5/10/15, documents cell suppression, fits a three-level linear mixed model, assesses inferential reproducibility, and writes all figures, tables, and a results bundle. |

A standalone test suite (`R/test_quantify_anonymization_impact.R`) validates the core function.

---

## Repository structure

```
.
├── R/
│   ├── quantify_anonymization_impact.R        # Core comparison function (base R only)
│   ├── cohort_description.R                    # Cohort characterization (Stage 2)
│   ├── final_analysis.R                        # Main analysis pipeline (Stage 3)
│   └── test_quantify_anonymization_impact.R    # Unit-test script for the core function
├── tests/
│   └── test_results_report.json                # 42 unit-test results (all passed)
├── ARX certificates/
│   ├── Multimedia Appendix 1 - 5-anonymized data certificate.pdf
│   ├── Multimedia Appendix 2 - 10-anonymized data certificate.pdf
│   └── Multimedia Appendix 3 - 15-anonymized data certificate.pdf
├── ADAQI Checklist/
│   ├── Multimedia Appendix 6 - ADAQI Checklist Template.pdf      # Blank template
│   └── Multimedia Appendix 7 - ADAQI Checklist filled out.pdf    # Worked example (this study)
├── data/                                       # Place input CSVs here (not included; see Data)
├── output/                                     # Generated figures, tables, .rds (created on run)
├── LICENSE                                     # MIT
├── CITATION.cff
└── README.md
```

---

## Requirements

**R version:** >= 4.1.0

| Package | Used by | Purpose |
|---------|---------|---------|
| `lme4` | `final_analysis.R` | Linear mixed models |
| `data.table` | `cohort_description.R` | Fast data manipulation |
| `flextable` | `cohort_description.R` | Word-compatible summary table |
| `officer` | `cohort_description.R` | DOCX export |

```r
install.packages(c("lme4", "data.table", "flextable", "officer"))
```

The core function (`quantify_anonymization_impact.R`) and the impact and reproducibility logic in `final_analysis.R` use **base R only** for zero-dependency reproducibility. Only the mixed-model fit requires `lme4`, and only the cohort summary table requires `data.table`, `flextable`, and `officer`.

---

## Data

The study uses inpatient encounter data from **Universitatsmedizin Mannheim (UMM)** for the period **January 2010 to September 2024**, anonymized with the [ARX Data Anonymization Tool](https://arx.deidentifier.org/) version 3.9.2.

> **Data availability.** The raw data cannot be publicly shared due to data-protection requirements. Source data have been archived and access can be made available for individual requests based on approval of the Ethics and Use and Access Committees. Full project documentation, including R source code, function implementation, and anonymization certificates, is publicly available in this repository under the MIT License.

Because of these requirements, **no CSV data files are included** in this repository. To reproduce the analysis, place the following **comma-separated** CSV files in `data/`.

| File | Description |
|------|-------------|
| `D0.csv` | Original (unanonymized) dataset |
| `D5.csv` | k = 5 anonymized |
| `D10.csv` | k = 10 anonymized |
| `D15.csv` | k = 15 anonymized |

**Schema (all four files):**

```
patient_id, encounter_id, admission_year, main_diagnosis_icd, los_days
```

- `patient_id` — pseudonymous patient identifier (enables the patient-level random effect)
- `encounter_id` — unique encounter identifier (direct identifier, removed before anonymization)
- `admission_year` — year of admission (retained for the temporal-drift check, never anonymized)
- `main_diagnosis_icd` — primary ICD-10-GM diagnosis code (quasi-identifier)
- `los_days` — hospital length of stay in days (quasi-identifier and analytical outcome)

Suppressed quasi-identifier cells appear as the token `*` (ARX default). On import, `main_diagnosis_icd` and `los_days` are read as character so the token survives, then `los_days` is coerced to numeric (suppressed values become `NA`). Rows are never dropped.

For `cohort_description.R`, place a single file `original data.csv` in `data/` with the columns `patient_id`, `main_diagnosis_icd`, `los_days` (and admission and discharge timestamps if cohort-level temporal summaries are required).

---

## Quick start

### 1. Configure the data path

`final_analysis.R` reads its input directory from the `DATA_DIR` variable, which can also be set with the `KANON_DATA_DIR` environment variable without editing the script.

```r
# Option A: set an environment variable before launching R
Sys.setenv(KANON_DATA_DIR = "/path/to/your/data")

# Option B: edit DATA_DIR near the top of final_analysis.R
DATA_DIR <- file.path("data")
```

Outputs are written to `OUTPUT_DIR`, which defaults to a `study_results_updated/` subfolder of `DATA_DIR`. Set it to `output/` if you prefer the in-repository location.

### 2. Run the analysis

From the `R/` directory (so the core function is found via `getwd()`):

```r
setwd("R")

# Stage 2 (optional, standalone): cohort description
source("cohort_description.R")

# Stage 3: full impact + reproducibility analysis
# (this sources quantify_anonymization_impact.R automatically)
source("final_analysis.R")
```

To run the unit tests for the core function:

```r
source("R/test_quantify_anonymization_impact.R")
```

### 3. Inspect outputs

All results are written to `OUTPUT_DIR`.

**Figures (PNG, 1200 x 1200 px):**

| File | Content |
|------|---------|
| `Fig_KS_D_los_days.png` | KS D divergence of LOS across k |
| `Fig_ICD_quality_metrics.png` | ICD Jaccard overlap and Cramer V (dual axis) |
| `Fig_ECDF_los_days.png` | Empirical CDF overlay for LOS |
| `Fig_density_los_days.png` | Kernel density overlay for LOS |
| `Fig_ICD_zipf_rank.png` | ICD rank-frequency (Zipf) plot |
| `Fig_los_mean_sd_shifts.png` | Mean and SD shift chart |
| `Fig_prevalence_los_panel.png` | ICD prevalence and LOS panel per dataset |
| `Fig_QQ_residuals_D{0,5,10,15}.png` | Residual QQ-plots per dataset |
| `Fig_QQ_BLUPs_D{0,5,10,15}.png` | ICD-3 random-intercept (BLUP) QQ-plots per dataset |
| `Fig_ICC_AIC_trajectory.png` | ICC and AIC across k |
| `Fig_BLUP_concordance_scatter.png` | BLUP concordance, original vs anonymized |
| `Fig_BLUP_attenuation_ratio.png` | BLUP attenuation ratios across k |
| `Fig_variance_architecture.png` | Variance decomposition and BLUP fidelity |
| `Fig_extreme_distortion.png` | Extreme individual BLUP distortions and sign reversals |
| `Fig_corroboration.png` | Distributional and categorical metrics vs concordance |

**Tables (CSV):**

| File | Content |
|------|---------|
| `summary_metrics.csv` | KS D, KS z, mean/SD shifts, Jaccard, Cramer V, verdicts |
| `data_point_removal.csv` | Cell-suppression documentation (no rows deleted) |
| `suppression_by_year.csv` | Suppression distribution across admission years |
| `cohort_characterization.csv` | Suppressed vs retained cohort characteristics |
| `NA_documentation.csv` | ARX-introduced missing values per dataset |
| `variance_components.csv` | Three-level variance components and ICC |
| `variance_inflation_interrogation.csv` | Between-diagnosis variance interrogation |
| `fixed_effects.csv` | Fixed-effect estimates (intercept, grand-mean LOS) |
| `BLUP_concordance.csv` | Spearman rho and Lin CCC with 95% bootstrap CIs |
| `BLUP_attenuation.csv` | Per-diagnosis BLUP attenuation ratios |
| `BLUP_extreme_distortion.csv` | Sign reversals and extreme magnitude changes |
| `prediction_accuracy.csv` | RMSE and MAE (log and day scales) |
| `model_fit_comparison.csv` | AIC, BIC, log-likelihood, group counts |
| `LMM_summary.csv` | Consolidated mixed-model summary |
| `test_results_report.json` | Internal consistency checks for the run |

**R object:**

- `study_results_full.rds` — complete results bundle for downstream reuse.

---

## Methodology

### Core comparison function

`quantify_anonymization_impact()` accepts an original and an anonymized data frame (and optionally a single column pair for an ad-hoc comparison) and returns a named list with column-presence, numeric, categorical, and ad-hoc results.

- **Numeric analysis.** KS D statistic, effective sample size n_eff = n1*n2/(n1+n2), scaled statistic z = sqrt(n_eff)*D, and min/max/mean/SD shifts. Counts use doubles to avoid 32-bit integer overflow at large n, and a `z_reason` field records why z is `NA` whenever a precondition fails.
- **Categorical analysis.** Jaccard overlap ratio (|intersection| / |union|), Cramer V = sqrt(chi2 / (n * min(r-1, c-1))), and a chi-square p-value computed over shared levels only.
- **Traffic-light verdict.** A green, yellow, or red verdict with handling recommendations. Numeric: green if D < 0.05 and the normalized mean shift <= 0.10 and the SD change is within +/-10%; yellow if D < 0.10 and the normalized mean shift <= 0.30; otherwise red. Categorical: green if Cramer V < 0.10 and Jaccard >= 0.80 and missing-data drift <= 5 percentage points; yellow if Cramer V < 0.30 and Jaccard >= 0.60 and drift <= 10 points; otherwise red.

### Linear mixed model

```
log(los_days + 1) ~ year + (1 | icd3) + (1 | patient_id)
```

- **Response.** Log-transformed LOS, with a +1 offset for same-day encounters.
- **Fixed effect.** Admission year, retained to verify that suppression is not concentrated in particular years (temporal-drift check).
- **Random effects.** A random intercept for the three-character ICD category (`icd3`) and a random intercept for `patient_id`, the latter addressing the non-independence of repeated encounters from the same patient.
- **Estimation.** REML for variance components, with maximum likelihood used for comparable AIC and BIC.
- **Rationale.** The mixed model is the measurement instrument for reproducibility, not a predictive model. Its variance components and diagnosis-level best linear unbiased predictions (BLUPs) provide exactly the quantities needed to ask whether the diagnosis-to-LOS signal survives anonymization.

### Reproducibility assessment

The diagnosis-level signal recovered from anonymized data is compared with the original across complementary dimensions.

| Quantity | Interpretation |
|----------|----------------|
| ICC (diagnosis, patient) | Variance architecture |
| Spearman rho of BLUPs (95% CI) | Rank ordering of diagnosis-specific effects |
| Lin concordance correlation coefficient (95% CI) | Magnitude agreement, with precision (Pearson r) and accuracy (Cb) components |
| BLUP sign reversals and MAD | Individual-effect stability |
| RMSE, MAE, AIC, BIC | Model fit and prediction accuracy |
| Grand-mean LOS stability | Systematic bias from suppression |

Confidence intervals for the Spearman correlation and the Lin CCC are obtained by bootstrap (`N_BOOT = 2000`, `set.seed(42)`).

---

## Testing

The core function is validated by **42 unit tests** (all passing; see `tests/test_results_report.json`). The test script `R/test_quantify_anonymization_impact.R` sources only the function definition and exercises identity invariants, numeric shift detection, categorical analysis, column-presence handling, vector alignment, backward-compatible argument aliases, verdict logic, helper functions, the return structure, and edge cases.

```r
source("R/test_quantify_anonymization_impact.R")
```

---

## Reproducibility notes

- The pipeline is deterministic. The random seed is fixed (`set.seed(42)`) before bootstrapping.
- The mixed-model fit holds at most one model in memory at a time and checkpoints each dataset's extract to disk, so an interrupted run resumes without refitting completed datasets.
- Figures are written at 1200 x 1200 px and 150 dpi.

---

## Citation

If you use this code, please cite the accompanying paper:

```bibtex
@article{kamdjewabo2026anonymization,
  title   = {Quantifying the Impact of Anonymization-Induced Clinical Data Quality
             Loss: A Methodological Case Study Using Primary Diagnosis Codes
             and Hospital Length of Stay},
  author  = {Kamdje Wabo, Gaetan},
  journal = {JMIR Medical Informatics},
  year    = {2026},
  note    = {Under review}
}
```

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

**Gaetan Kamdje Wabo, M.Sc.**
Department of Biomedical Informatics,
Medical Faculty Mannheim, Heidelberg University
