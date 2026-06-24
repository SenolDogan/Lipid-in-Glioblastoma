# Data Availability

## Primary input data

Discovery cohort lipidomics, metabolomics, and clinical annotations: CPTAC/TCGA GBM 2021 (`gbm_cptac_2021`) via [cBioPortal](https://www.cbioportal.org/study/summary?id=gbm_cptac_2021).

External benchmarking: de-identified plasma metabolomics from Fontanilles et al. (2024), processed outputs under `external_validation/fontanilles_2024/`.

## Processed outputs in this repository

All analysis-ready files are listed in **Supplementary Table S17** (`Supplementary/tables/Supplementary_Table_S17.csv`).

Key files:

- `results/cluster_assignments.csv` — exploratory subtype labels (NMF Subtype A/B)
- `results/lipid_metabolite_correlations.csv` — correlation atlas edge list (1,676 edges)
- `results/differential_features.csv` — subtype contrasts with Cohen's d and FDR
- `results/feature_importance_rf.csv` — Random Forest benchmark rankings
- `Feature_Names_Mapping.xlsx` — feature display-name mapping

## Supplementary materials

- Figures S1–S8: `Supplementary/figures/` (correlation heatmap/network views are in main-text **Figure 3**)
- Tables S1–S22: `Supplementary/tables/`
- Consolidated bundles: `Supplementary/Supplementary_Figures.docx`, `Supplementary/Supplementary_Tables.xlsx`

## Code and quick verification

- Analysis and figure-generation scripts: `scripts/` and `requirements.txt`
- Reviewer quick-check: `python scripts/verify_release.py` (see root `README.md`)
