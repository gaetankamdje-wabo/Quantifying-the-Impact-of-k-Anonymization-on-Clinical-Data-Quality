# Quantifying the Impact of k-Anonymization on Clinical Data Quality

**Companion code** for the manuscript:

> Kamdje Wabo, G. *Quantifying the Impact of Anonymization-Induced Clinical Data Quality Loss: An Advanced Methodological Study Using Primary Diagnosis Codes and Hospital Length of Stay.* JMIR Medical Informatics (under review), 2026.

---

## Overview

This repository provides the complete R analysis pipeline for a study that systematically measures how **ARX k-anonymity** (k = 5, 10, 15) with cell-level suppression affects the statistical utility of hospital encounter data — specifically **ICD-10-GM primary diagnosis codes** and **length of stay (LOS)**.

The analysis proceeds in three stages:

| Stage | Script | Purpose |
|-------|--------|---------|
| **1. Core function** | `R/quantify_anonymization_impact.R` | Reusable R function: column presence checks, numeric distributional analysis (KS D, mean/SD shifts), categorical analysis (Jaccard overlap, Cramér's V, χ²), and traffic-light verdicts |
| **2. Cohort description** | `R/cohort_description.R` | Descriptive statistics and publication-ready figures for the original hospital encounter dataset |
| **3. Main analysis** | `R/final_analysis.R` | Full impact analysis: runs Stage 1 across k = 5/10/15, generates 15+ scholarly figures, fits Linear Mixed Models (LMM), and assesses statistical reproducibility |

---

## Repository Structure

```
.
├── R/
│   ├── quantify_anonymization_impact.R   # Core comparison function (660 lines)
│   ├── cohort_description.R              # Cohort characterization (338 lines)
│   └── final_analysis.R                  # Main analysis pipeline (2106 lines)
├── tests/
│   └── test_results_report.json          # 30 unit tests (all passed)
├── data/                                 # Place your CSV files here (not included)
├── output/                               # Generated tables, figures, .rds files
├── LICENSE
├── CITATION.cff
└── README.md
```

---

## Requirements

**R version:** ≥ 4.1.0

**Packages:**

| Package | Purpose | Stage |
|---------|---------|-------|
| `data.table` | Fast data manipulation | 2 |
| `ggplot2` | Publication-quality plots | 2 |
| `flextable` | Word-compatible tables | 2 |
| `officer` | DOCX export | 2 |
| `lme4` | Linear Mixed Models | 3 |
| `scales` | Axis formatting | 2 |

Install all dependencies:

```r
install.packages(c("data.table", "ggplot2", "flextable", "officer", "lme4", "scales"))
```

No additional packages are needed. Stage 1 (`quantify_anonymization_impact.R`) and the core of Stage 3 use **base R only** to ensure zero-dependency reproducibility.

---

## Data

The study uses hospital encounter data from **Universitätsmedizin Mannheim** (2010–2024), anonymized with [ARX Data Anonymization Tool](https://arx.deidentifier.org/) v3.9+.

**Required input files** (semicolon-separated CSV, not included for privacy):

| File | Description | Columns |
|------|-------------|---------|
| `D0.csv` | Original dataset (~720,359 encounters) | `encounter_id; main_diagnosis_icd; los_days` |
| `D5.csv` | k = 5 anonymized | Same schema |
| `D10.csv` | k = 10 anonymized | Same schema |
| `D15.csv` | k = 15 anonymized | Same schema |

Suppressed values appear as `*` (ARX default). Place all four files in the `data/` directory.

---

## Quick Start

### 1. Configure data paths

In each R script, update the `DATA_DIR` / `path_input` variable to point to your `data/` folder:

```r
# In final_analysis.R (line 45–48):
DATA_DIR <- file.path("data")

# In cohort_description.R (line 16):
path_input <- file.path("data", "original data.csv")
```

### 2. Run the analysis

```r
# Step 1: Source the core function
source("R/quantify_anonymization_impact.R")

# Step 2: Cohort description (optional, standalone)
source("R/cohort_description.R")

# Step 3: Full analysis (sources Step 1 automatically)
source("R/final_analysis.R")
```

### 3. Inspect outputs

All results are saved to the `output/` directory:

**Figures (PDF, 300 DPI):**
- `Fig1` – Record suppression rates per k-level
- `Fig2` – KS D effect size trajectory (LOS)
- `Fig3` – ICD Jaccard overlap + Cramér's V (dual-axis)
- `Fig4` – Empirical CDF overlay (LOS)
- `Fig5` – 2×2 summary dashboard panel
- `Fig6` – Kernel density overlay (LOS)
- `Fig7` – ICD rank-frequency Zipf plot
- `Fig8` – Mean/SD shift lollipop chart
- `Fig9` – Residual QQ-plots (per dataset)
- `Fig9b` – ICD prevalence + LOS panel (D0/D5/D10/D15)
- `Fig9e` – BLUP QQ-plots
- `Fig10–15` – LMM reproducibility: ICC trajectory, BLUP concordance, attenuation, variance decomposition, prediction accuracy, grand mean stability

**Tables (CSV):**
- `Table1` – Summary metrics across k-levels
- `Table2` – NA documentation (ARX-introduced suppressions)
- `Table3` – LMM comprehensive summary
- `Table4` – BLUP concordance (Spearman ρ)
- `Table5` – Variance components
- `Table6` – BLUP attenuation ratios
- `Table7` – Fixed effects
- `Table8` – Prediction accuracy (RMSE/MAE)
- `Table9` – Model fit comparison (AIC/BIC)

**R objects:**
- `study_results_full.rds` – Complete results bundle for downstream reuse

---

## Methodology

### Core Comparison Function

`quantify_anonymization_impact()` accepts an original and an anonymized data frame and returns:

- **Column presence:** Detects added/removed columns
- **Numeric analysis:** KS D statistic, √(n_eff)·D signal gauge, min/max/mean/SD shifts
- **Categorical analysis:** Jaccard overlap ratio, Cramér's V, χ² p-value (shared levels)
- **Traffic-light verdict:** Green (✅), Yellow (⚠️), or Red (❌) with handling recommendations

### Linear Mixed Model (Method subsection)

```
log(los_days + 1) ~ 1 + (1 | icd3)
```

- **Response:** Log-transformed LOS (+1 offset for same-day cases)
- **Random effect:** Per-ICD-3 intercept (high-cardinality diagnosis groups)
- **Estimation:** REML for variance components; ML for AIC/BIC comparison
- **Rationale:** LMM handles hundreds of ICD-3 categories via shrinkage, avoiding the overfitting and instability of fixed-effect GLMs with high-cardinality categorical predictors

### Reproducibility Assessment (Method subsection)

Five dimensions of reproducibility are evaluated across k-levels:

| Metric | Interpretation |
|--------|---------------|
| Spearman ρ (BLUPs) | Rank ordering of diagnosis-specific effects preserved |
| BLUP attenuation ratio | Effect magnitudes preserved (ratio ≈ 1) |
| ICC trajectory | Variance architecture preserved |
| RMSE/MAE stability | Predictive utility preserved |
| Grand mean stability | No systematic bias from suppression |

---

## Testing

The core function is validated by **30 unit tests** across 10 categories:

| Category | Tests | Coverage |
|----------|-------|----------|
| Identity Invariants | A01–A06 | Identical data → D=0, V=0, Jaccard=1, Green |
| Numeric Shift Detection | B01–B06 | Constant shifts, perturbations, variance compression |
| Categorical Analysis | C01–C05 | Disjoint/partial overlap, single level, missing values |
| Column Presence | D01–D03 | Added/removed/identical columns |
| Alignment | E01 | Length-matched comparisons |
| Backward Compatibility | F01–F02 | Deprecated parameter aliases |
| Verdict Logic | G01–G03 | Exact verdict string matching |
| Helper Functions | H01–H05 | `fmt_p`, `icd_chapter`, `los_midpoint` |
| Return Structure | I01–I05 | All required output keys and sub-elements |
| Edge Cases | J01–J04 | NULL input, missing ad-hoc, character-only frames |

**Result: 30/30 passed** (see `tests/test_results_report.json`).

---

## Citation

If you use this code, please cite:

```bibtex
@article{kamdjewabo2026anonymization,
  title   = {Quantifying the Impact of Anonymization-Induced Clinical Data Quality
             Loss: An Advanced Methodological Study Using Primary Diagnosis Codes
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
