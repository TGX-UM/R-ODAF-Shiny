R-ODAF-Shiny: Graphical Interface for the R-ODAF DEG Pipeline

This Shiny app provides a graphical interface for the R-ODAF DESeq2 pipeline.
It enables fixed-parameter differential gene expression (DEG) analysis with built-in visualization and reporting.
All analyses follow R-ODAF’s regulatory criteria for reproducibility and traceability.

Live Deployment

The containerized app is available at:
https://r-odaf.nl

Uploads are processed in-session and discarded automatically after session end.

Key Features

RNA-seq counts + metadata input

Filters for regulatory consistency:

Minimum coverage: 5M reads/sample

CPM threshold: 1

FDR cutoff: 0.01

PCA outlier variance threshold: 20%

Interactive exploration:

PCA (before/after outlier removal)

Final DEG tables (exportable)

Heatmaps, top-gene and average-expression plots

Exportable results and summary reports

TempO-Seq module is planned.

Run Locally using Docker

The validated R-ODAF-Shiny production environment is defined by the included Dockerfile. The Docker configuration pins the R base image, CRAN snapshot, Bioconductor release, and key package versions used by the deployed application.

See the Dockerfile and docs/sessionInfo.txt for details of the software environment.

Inputs expected

counts: genes × samples (CSV/TSV)

metadata: samples × attributes (must include a condition column)

Reproducibility

This repository includes:

Dockerfile – version-pinned production software environment

docs/sessionInfo.txt – R and package versions for the validated production environment

Fixed analysis thresholds to support reproducible results across runs

See also:

docs/PRIVACY.md – data handling

docs/REPRODUCIBILITY.md – fixed parameters and validation notes

Citation

Please cite the R-ODAF framework and this Shiny implementation:


R-ODAF-Shiny (this work)
Citation will follow soon!

R-ODAF (framework)
Verheijen M.C., Meier M.J., Asensio J.O., Gant T.W., Tong W., Yauk C.L., Caiment F. (2022).
R-ODAF: Omics Data Analysis Framework for Regulatory Application.
Regulatory Toxicology and Pharmacology, 131, 105143. https://doi.org/10.1016/j.yrtph.2022.105143

Credits

Developed by Saad Lodhi.
Part of the R-ODAF pipeline, Department of Translational Genomics at Maastricht University.
