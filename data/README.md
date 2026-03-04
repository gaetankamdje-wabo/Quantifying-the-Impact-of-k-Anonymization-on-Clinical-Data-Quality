Place the following semicolon-separated CSV files in this directory:

- D0.csv   (original, unanonymized dataset)
- D5.csv   (k = 5 anonymized)
- D10.csv  (k = 10 anonymized)
- D15.csv  (k = 15 anonymized)

For cohort_description.R, additionally place:
- original data.csv  (full hospital encounter dataset with timestamps)

Schema: encounter_id ; main_diagnosis_icd ; los_days

These files are not included in the repository for patient privacy reasons.
