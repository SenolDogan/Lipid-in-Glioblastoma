# Lipid-in-Glioblastoma

Harmonized, analysis-ready integrated lipidomics and metabolomics data for glioblastoma (GBM). This repository supports the *Scientific Data* Data Descriptor and is intended for **immediate reuse** without rebuilding preprocessing pipelines.

**Manuscript:** `manuscript/LIPID.docx`  
**Repository:** https://github.com/SenolDogan/Lipid-in-Glioblastoma

---

## Who is this for?

| You want to… | Start here | Time |
|--------------|------------|------|
| **Verify the release (reviewer / editor)** | [10-minute reviewer check](#10-minute-reviewer-check) | ~10 min |
| **Reuse processed tables** | [Two reuse workflows](#two-reuse-workflows) | ~5–15 min |
| **Re-run the full CPTAC pipeline from raw matrices** | [Full re-analysis (optional)](#full-re-analysis-optional) | hours + raw data download |

You do **not** need to read the manuscript first if you only want to verify or reuse the processed outputs.

---

## 10-minute reviewer check

Confirm that the repository is complete, readable, and consistent with the Data Descriptor **without** re-downloading CPTAC raw data.

### Step 0 — Requirements

- Python 3.9+ (3.10+ recommended)
- `git`
- ~200 MB free disk space

### Step 1 — Clone and install (~2 min)

```bash
git clone https://github.com/SenolDogan/Lipid-in-Glioblastoma.git
cd Lipid-in-Glioblastoma
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### Step 2 — Run the release verification script (~1 min)

```bash
python scripts/verify_release.py
```

**Expected:** all lines print `[PASS]` and the script exits with code `0`.

Numbers checked against the manuscript:

- **74** matched patients in `results/cluster_assignments.csv`
- **1,676** lipid–metabolite edges at |ρ| ≥ 0.3, BH-FDR < 0.05
- NMF subtype counts **62 / 12** (`nmf_subtype` 0 / 1)
- Network permutation empirical **p < 0.05** in `results/acceptance_network_permutation_summary.csv`
- **8** supplementary figure PNGs (S1–S8) and Supplementary Table S17 manifest

### Step 3 — Open the file manifest (~2 min)

```bash
head -n 8 Supplementary/tables/Supplementary_Table_S17.csv
```

Or open `Supplementary/tables/Supplementary_Table_S17.csv` in Excel or LibreOffice. This is the **master index** of every shared file.

### Step 4 — Spot-check Workflow 1 (subtype labels) (~2 min)

```python
import pandas as pd

df = pd.read_csv("results/cluster_assignments.csv")
print(df.shape)                     # expect (74, 5)
print(df["nmf_subtype"].value_counts())
```

Subtype A/B in the manuscript follow **`nmf_subtype`** (0 = Subtype A, 1 = Subtype B; Table 1; Figure 2).

### Step 5 — Spot-check Workflow 2 (correlation atlas) (~2 min)

```python
import pandas as pd

edges = pd.read_csv("results/lipid_metabolite_correlations.csv")
print(len(edges))                   # expect 1676
print(edges.head())
```

Figure 3 in the manuscript shows illustrative heatmap and network views of this atlas (edge list also in Supplementary Table S4).

### Step 6 — Optional: view supplementary figures (~1 min)

Open `Supplementary/Supplementary_Figures.docx` or any PNG under `Supplementary/figures/`.

---

## Two reuse workflows

These mirror the Data Overview in the manuscript.

### Workflow 1 — Subtype stratification (Figure 2; Table 1)

| Step | File | Action |
|------|------|--------|
| 1 | `results/cluster_assignments.csv` | Load NMF Subtype A/B (`nmf_subtype`) |
| 2 | `results/differential_features.csv` | Subtype contrasts (Cohen's d, FDR) |
| 3 | `Supplementary/tables/Supplementary_Table_S01.csv` – `S03.csv` | Cohort summary, assignments, bootstrap stability |
| 4 | `Feature_Names_Mapping.xlsx` | Map feature IDs to display names |

### Workflow 2 — Correlation screening (Figure 3)

| Step | File | Action |
|------|------|--------|
| 1 | `results/lipid_metabolite_correlations.csv` | Full atlas (1,676 edges) |
| 2 | `results/hub_molecule_analysis.csv` | Hub centrality scores |
| 3 | `Supplementary/tables/Supplementary_Table_S04.csv` – `S07.csv` | Edge list, network stats, top hubs |
| 4 | `results/lipid_class_abundance.csv`, `results/metabolite_category_abundance.csv` | Class-level summaries (Supplementary Figures S7–S8) |

---

## Repository layout

| Path | Description |
|------|-------------|
| `manuscript/` | Main manuscript (`LIPID.docx`) |
| `Supplementary/figures/` | Supplementary Figure S1–S8 (PNG) |
| `Supplementary/tables/` | Supplementary Table S1–S22 (CSV) |
| `Supplementary/Supplementary_Figures.docx` | All supplementary figures with captions and legends |
| `Supplementary/Supplementary_Tables.xlsx` | All supplementary tables in one workbook |
| `Supplementary/Supplementary_Table_19_Lipid_biological_source_inventory.xlsx` | S19 lipid biological source inventory |
| `results/` | Processed analysis-ready CSV outputs |
| `Feature_Names_Mapping.xlsx` | Feature ID to display-name mapping |
| `figures/main/` | Main manuscript figure sources |
| `figures/supplementary/` | Mirror of supplementary figure PNGs |
| `scripts/` | Analysis, validation, and plotting scripts |
| `scripts/verify_release.py` | **Reviewer quick-check script** |
| `external_validation/fontanilles_2024/` | External plasma metabolomics benchmarking outputs |
| `docs/` | Supporting documentation |

---

## Main manuscript figures

| Figure | Content | Primary files |
|--------|---------|---------------|
| **Figure 1** | Study workflow | `scripts/plots/figure1_study_workflow_flowchart.py` |
| **Figure 2** | Subtype discovery (t-SNE, RF importance, Cohen's d) | `results/cluster_assignments.csv`, `results/feature_importance_rf.csv`, `results/differential_features.csv` |
| **Figure 3** | Correlation atlas (heatmap + network panels) | `results/lipid_metabolite_correlations.csv`, `Supplementary/tables/Supplementary_Table_S04.csv` |

---

## Supplementary figures (S1–S8)

| ID | Description |
|----|-------------|
| S1 | Principal Component Analysis |
| S2 | Combined subtype-separating and network-based prioritization outputs |
| S3 | Dual-purpose candidate lipids (subtype discrimination + network centrality) |
| S4 | Metabolite category abundance by metabolic-lipidomic subtype |
| S5 | Survival feature importance (Random Forest regression) |
| S6 | Random Forest feature importance for subtype prediction |
| S7 | Lipid-class remodeling and direction of subtype-associated features |
| S8 | Lipid class vs metabolite category correlation matrix |

Captions and legends: `Supplementary/Supplementary_Figures.docx`

> **Note:** The correlation heatmap and network views that were previously archived as supplementary figures are now **Figure 3** in the main Data Descriptor.

---

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
| S14 | Lipid-class vs metabolite-category correlations (supports Supplementary Figure S8) |
| S15 | Lipid class abundances by patient |
| S16 | Metabolite category abundances by patient |
| S17 | **File manifest** (start here) |
| S18 | External benchmarking summary (Fontanilles et al. 2024) |
| S19 | Biological source inventory for 581 lipid species |
| S20 | Top subtype-discriminating lipid features |
| S21 | Subtype-specific signatures and dual-purpose hub lipids |
| S22 | Significant lipid class–metabolite category correlations |

---

## File format conventions

All tabular outputs in `results/` and `Supplementary/tables/` are:

- **Format:** CSV (comma-separated values)
- **Encoding:** UTF-8
- **Header row:** First row contains column names (stable across releases)
- **Delimiter:** comma (`,`)
- **Decimal separator:** period (`.`)

### Key output schemas

| File | Column names | Description |
|------|--------------|-------------|
| `results/cluster_assignments.csv` | `sample_id`, `kmeans_subtype`, `hierarchical_subtype`, `nmf_subtype`, `vital_status` | CPTAC sample IDs and subtype labels |
| `results/lipid_metabolite_correlations.csv` | `feature1`, `feature2`, `correlation`, `pvalue`, `fdr` | Significant lipid–metabolite edges |
| `results/differential_features.csv` | `cluster`, `feature`, `mean_cluster`, `mean_other`, `fold_change`, `pvalue`, `cohens_d`, `fdr` | Subtype differential statistics |
| `results/feature_importance_rf.csv` | `feature`, `importance` | Random Forest Gini importance |
| `results/hub_molecule_analysis.csv` | `molecule`, `degree_centrality`, `betweenness_centrality`, `closeness_centrality`, `eigenvector_centrality`, `degree`, `hub_score` | Network hub metrics |

See Supplementary Table S17 for the complete file manifest and field-level descriptions.

---

## Data policy

This repository includes **processed** result tables and derived outputs. Large raw CPTAC quantification matrices are not included; input data are obtained from [cBioPortal gbm_cptac_2021](https://www.cbioportal.org/study/summary?id=gbm_cptac_2021) as described in the manuscript.

---

## Full re-analysis (optional)

Re-running `scripts/lipidomics_metabolomics_analysis.py` from scratch requires:

1. CPTAC/TCGA GBM raw matrices from [cBioPortal](https://www.cbioportal.org/study/summary?id=gbm_cptac_2021)
2. Local path configuration in the analysis scripts (paths may need editing on your machine)

For most users and reviewers, the **processed `results/` tables are the primary distribution**. Full re-analysis is not required to validate the Data Descriptor.

Additional scripts:

- `scripts/acceptance_boost_analyses.py` — network permutation and ML benchmark summaries
- `scripts/validate_fontanilles_2024.py` — external plasma metabolomics benchmarking
- `scripts/plots/` — figure regeneration utilities

Python dependencies: `requirements.txt`  
Fixed random seed: **42** (see manuscript Methods)

---

## Citation

If you use this resource, please cite the Data Descriptor manuscript and this repository (Reference 39 in the manuscript).
