# output/

This directory is populated automatically when the analysis scripts run. By default
`final_analysis.R` writes to a `study_results_updated/` subfolder of `DATA_DIR`; set
`OUTPUT_DIR` to this `output/` folder to keep results in the repository instead.

Generated artifacts include:

- **Figures (PNG, 1200 x 1200 px):** `Fig_KS_D_los_days`, `Fig_ICD_quality_metrics`,
  `Fig_ECDF_los_days`, `Fig_density_los_days`, `Fig_ICD_zipf_rank`, `Fig_los_mean_sd_shifts`,
  `Fig_prevalence_los_panel`, `Fig_QQ_residuals_D{0,5,10,15}`, `Fig_QQ_BLUPs_D{0,5,10,15}`,
  `Fig_ICC_AIC_trajectory`, `Fig_BLUP_concordance_scatter`, `Fig_BLUP_attenuation_ratio`,
  `Fig_variance_architecture`, `Fig_extreme_distortion`, `Fig_corroboration`.
- **Tables (CSV):** `summary_metrics`, `data_point_removal`, `suppression_by_year`,
  `cohort_characterization`, `NA_documentation`, `variance_components`,
  `variance_inflation_interrogation`, `fixed_effects`, `BLUP_concordance`,
  `BLUP_attenuation`, `BLUP_extreme_distortion`, `prediction_accuracy`,
  `model_fit_comparison`, `LMM_summary`.
- **Run report:** `test_results_report.json` (internal consistency checks).
- **R object:** `study_results_full.rds` (complete results bundle).
