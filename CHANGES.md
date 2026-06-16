# Update log

This release aligns the companion code with the revised manuscript
("A Methodological Case Study ...", JMIR Medical Informatics, 2026).

## Code
- **Resolved duplicate core function.** Removed `quantify_anonymization_impact().R`
  (the parenthesised filename) and kept the newer implementation under the clean,
  shell- and git-safe name `quantify_anonymization_impact.R`. This version adds a
  large-n integer-overflow guard for `n_eff` and a `z_reason` field documenting why
  the scaled statistic z is `NA` in any given case. `final_analysis.R` and the test
  script now source the clean filename.
- **Exposed module-level helpers.** `fmt_p()`, `icd_chapter()`, and `los_midpoint()`
  are now defined at the top level of the core file, with the sentinel line the test
  suite uses for selective sourcing.
- **Corrected one test contract.** `los_midpoint(7, 3)` now asserts `8`
  (third bin) to match the documented binning formula
  `bin_index = floor((los-1)/w); midpoint = bin_index*w + ceil(w/2)`. The test script
  and `tests/test_results_report.json` were updated accordingly.

## Documentation
- **Rewrote `README.md`** to match the revised study: five-column comma-separated
  schema (`patient_id, encounter_id, admission_year, main_diagnosis_icd, los_days`),
  cell-level suppression with the `*` token and no row deletion, the three-level mixed
  model `log(los_days + 1) ~ year + (1|icd3) + (1|patient_id)`, the actual generated
  figure and table filenames, the `KANON_DATA_DIR` configuration, the 42-test suite,
  the ADAQI checklists, the renamed Multimedia Appendix certificates, and the data
  availability statement.
- **Updated `data/`, `output/`, and `tests/` READMEs** to the new schema, output
  filenames, and 42-test count.
- **Updated `CITATION.cff`** title to the case-study wording.
- **Extended `.gitignore`** to cover the default `study_results_updated/` output
  location and `.ckpt_*.rds` checkpoints.

## Note
The raw CSV data are not included. They cannot be shared publicly due to
data-protection requirements; access can be requested subject to Ethics and Use and
Access Committee approval (see README, Data availability).
