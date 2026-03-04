# Unit Tests

`test_results_report.json` contains the results of 30 unit tests validating the
core function `quantify_anonymization_impact()`.

All 30 tests pass. Categories covered:

- **A (Identity Invariants):** Identical data returns D=0, V=0, Jaccard=1, Green verdict
- **B (Numeric Shift Detection):** Constant shifts, perturbations, variance compression
- **C (Categorical Analysis):** Disjoint/partial overlap, single level, missing values
- **D (Column Presence):** Added, removed, and identical columns
- **E (Alignment):** Length-matched vector comparisons
- **F (Backward Compatibility):** Deprecated parameter aliases
- **G (Verdict Logic):** Exact Green/Yellow/Red string matching
- **H (Helper Functions):** `fmt_p`, `icd_chapter`, `los_midpoint`
- **I (Return Structure):** All 5 top-level keys and sub-elements validated
- **J (Edge Cases):** NULL input, missing ad-hoc, character-only frames

To re-run the tests, use the test script (not included; tests were executed
in an interactive R session using synthetic data).
