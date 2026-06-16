# data/

Place the input CSV files here. They are **not included** in the repository because the
raw data cannot be publicly shared due to data-protection requirements. Source data have
been archived and access can be made available for individual requests based on approval
of the Ethics and Use and Access Committees.

## Files required by `final_analysis.R` (comma-separated CSV)

- `D0.csv`   — original (unanonymized) dataset
- `D5.csv`   — k = 5 anonymized
- `D10.csv`  — k = 10 anonymized
- `D15.csv`  — k = 15 anonymized

**Schema (all four files):**

```
patient_id, encounter_id, admission_year, main_diagnosis_icd, los_days
```

Suppressed quasi-identifier cells appear as the token `*` (ARX default). `main_diagnosis_icd`
and `los_days` are read as character so the token survives, then `los_days` is coerced to
numeric (suppressed values become `NA`). Rows are never dropped.

## File required by `cohort_description.R`

- `original data.csv` — full hospital encounter dataset with columns `patient_id`,
  `main_diagnosis_icd`, `los_days` (plus admission and discharge timestamps if temporal
  summaries are needed).
