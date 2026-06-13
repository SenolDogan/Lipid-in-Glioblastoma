# Lipid-in-Glioblastoma

Harmonized, analysis-ready integrated lipidomics and metabolomics data for glioblastoma (GBM), with correlation and network annotations. This repository supports the *Scientific Data* Data Descriptor manuscript and provides processed tables, supplementary materials, analysis code, and figure-generation scripts.

**Manuscript:** `manuscript/LIPID.docx`  
**Repository URL:** https://github.com/SenolDogan/Lipid-in-Glioblastoma

## Quick start

1. Open **Supplementary Table S17** (`Supplementary/tables/Supplementary_Table_S17.csv`) as the file manifest.
2. Load `results/cluster_assignments.csv` to map sample identifiers to exploratory Subtype A/B labels.
3. Use `Feature_Names_Mapping.xlsx` for consistent feature display names in plots and tables.
4. Filter `results/lipid_metabolite_correlations.csv` (|ρ| ≥ 0.3, BH-FDR < 0.05) or `results/differential_features.csv` for custom reuse workflows.

## Repository layout

| Path | Description |
|------|-------------|
| `manuscript/` | Main manuscript (`LIPID.docx`), consolidated supplementary DOCX/XLSX |
| `Supplementary/figures/` | Supplementary Figure S1–S10 (PNG) |
| `Supplementary/tables/` | Supplementary Table S1–S22 (CSV) |
| `Supplementary/Supplementary_Figures.docx` | All supplementary figures with captions |
| `Supplementary/Supplementary_Tables.xlsx` | All supplementary tables in one workbook |
| `Supplementary/Supplementary_Table_19_Lipid_biological_source_inventory.xlsx` | S19 lipid biological source inventory (also exported as S19 CSV) |
| `results/` | Processed analysis-ready CSV outputs |
| `Feature_Names_Mapping.xlsx` | Feature ID to display-name mapping |
| `figures/main/` | Main manuscript figures (where available) |
| `figures/supplementary/` | Mirror of supplementary figure PNGs and plot sources |
| `scripts/` | Analysis, validation, and plotting scripts |
| `external_validation/fontanilles_2024/` | External plasma metabolomics benchmarking outputs |
| `docs/` | Supporting documentation |

## Main manuscript figures

| Figure | Description |
|--------|-------------|
| Figure 1 | Study workflow |
| Figure 2 | Integrated subtype discovery, feature prioritization, and subtype-associated molecular differences |

## Supplementary figures (S1–S10)

| ID | Description |
|----|-------------|
| S1 | Principal Component Analysis |
| S2 | Combined subtype-separating and network-based prioritization outputs |
| S3 | Dual-purpose candidate lipids (subtype discrimination + network centrality) |
| S4 | Metabolite category abundance by metabolic-lipidomic subtype |
| S5 | Survival feature importance (Random Forest regression) |
| S6 | Random Forest feature importance for subtype prediction |
| S7 | Top 20 lipid correlation heatmap with metabolites |
| S8 | Lipid–metabolite correlation network (top 100 FDR-significant associations) |
| S9 | Lipid-class remodeling and direction of subtype-associated features |
| S10 | Lipid class–metabolite category correlation matrix |

## Supplementary tables (S1–S22)

| ID | Description |
|----|-------------|
| S1 | Analysis-ready cohort and feature summary |
| S2 | Cluster assignments (k-means, hierarchical, NMF) with vital status |
| S3 | Bootstrap stability results |
| S4 | 1,676 significant lipid–metabolite correlations |
| S5 | Network statistics |
| S6 | Hub molecule analysis |
| S7 | Top 20 hub molecules |
| S8 | Novel/uncommon candidates from differential and hub analyses |
| S9 | Differential features between subtypes |
| S10 | Pathway enrichment results |
| S11 | PCA loadings (PC1–PC5) |
| S12 | Random Forest feature importance |
| S13 | Cross-validation results |
| S14 | Lipid-class vs metabolite-category correlations (supports S10) |
| S15 | Lipid class abundances by patient |
| S16 | Metabolite category abundances by patient |
| S17 | **File manifest** (start here) |
| S18 | External benchmarking summary (Fontanilles et al. 2024) |
| S19 | Biological source inventory for 581 lipid species |
| S20 | Top subtype-discriminating lipid features |
| S21 | Subtype-specific signatures and dual-purpose hub lipids |
| S22 | Significant lipid class–metabolite category correlations |

## Data policy

This repository includes **processed** result tables and derived outputs. Large raw CPTAC quantification matrices are not included; input data are obtained from [cBioPortal gbm_cptac_2021](https://www.cbioportal.org/study/summary?id=gbm_cptac_2021) as described in the manuscript.

## Reproducibility

- Python dependencies: `requirements.txt`
- Core analysis: `scripts/`
- Figure generation: `scripts/plots/` (Integration project) or `Lipidomics_Metabolomics_Integration/plots/` locally
- Fixed random seed: 42 (see manuscript Methods)

## Citation

If you use this resource, please cite the Data Descriptor manuscript and this repository (Reference 39 in the manuscript).
