# tests/

`test_results_report.json` contains the results of the unit tests that validate the core
function `quantify_anonymization_impact()`. The test driver is `R/test_quantify_anonymization_impact.R`,
which sources only the function definition (up to its definition sentinel) and then runs
the suite.

**Result: 42/42 tests passed.**

Categories covered:

- **A — Identity invariants:** identical data returns KS D = 0, Cramer V = 0, Jaccard = 1, green verdict.
- **B — Numeric shift detection:** constant shifts, perturbations, variance compression.
- **C — Categorical analysis:** disjoint and partial overlap, single level, missing values.
- **D — Column presence:** added, removed, and identical columns.
- **E — Alignment:** length-matched vector comparisons.
- **F — Backward compatibility:** deprecated argument aliases.
- **G — Verdict logic:** exact green/yellow/red string matching.
- **H — Helper functions:** `fmt_p`, `icd_chapter`, `los_midpoint`.
- **I — Return structure:** all top-level keys and sub-elements present.
- **J — Edge cases:** NULL input, missing ad-hoc pair, character-only frames, large-n overflow guard.

To re-run:

```r
source("R/test_quantify_anonymization_impact.R")
```
