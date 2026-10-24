# Reproducibility and Analytical Criteria

The **R-ODAF-Shiny** interface follows the analytical criteria established in the R-ODAF framework.  
Core thresholds and filters are fixed to reduce user-induced analytical variability and support reproducible results across sessions.

### Parameters

| Criterion | Value / Rule | Description |
|------------|---------------|-------------|
| Minimum coverage | 5 million reads per sample | Ensures adequate sequencing depth |
| CPM threshold | 1 | Filters out very lowly expressed genes |
| FDR cutoff | 0.01 | Controls false discovery rate in DESeq2 results |
| Presence filter | ≥ 75% of samples within either comparison group | Removes sparsely detected genes |
| PCA outlier threshold | 20% variance | Applies the predefined PCA-based outlier criterion |

### Provenance and environment control

- The `Dockerfile` defines the version-pinned production software environment.
- `docs/sessionInfo.txt` records the R and package versions of the validated production environment.
- Each analysis generates a run summary containing the analysis settings, filtering outcomes, and run timing, and users can export `sessionInfo()` from the application.

### Validation

Numerical equivalence was evaluated against the reference R-ODAF differential-expression workflow using the bundled human colon organoid doxorubicin dataset.