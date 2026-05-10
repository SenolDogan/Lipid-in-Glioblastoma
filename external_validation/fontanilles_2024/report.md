## External validation: Fontanilles et al. 2024 (Acta Neuropathologica Communications)

- **Source XLSX**: `external_validation/Fontanilles_2024_s40478-024-01861-5_MOESM1.xlsx`
- **Outputs**: `external_validation/fontanilles_2024/_validation`

### Data availability (this dataset works for validation A+B)
- **Has quantified metabolomics matrix**: yes (`Table_S3`, Metabolomics rows)
- **Has overall survival**: yes (`Table_S1`: `OS DELAY` + `OS`)

### QC summary (Inclusion timepoint)
- **Patients with Inclusion + OS**: 33 patients (33 inclusion rows)
- **Metabolomics variables tested (median-split Cox)**: 180
- **Nominally significant (p<0.05)**: 4 (2.2%)

### Reproducibility vs authors’ Supplementary Table_S9
- **Variables comparable (intersection)**: 148
- **HR direction match**: 51/148 (34.5%)
- **File**: `inclusion_vs_authors_table_s9.csv`

### Overlap vs our internal Cox features
- **Exact-name overlap count**: 0 (naming/platform differences expected).

### Artifacts
- `inclusion_univariate_cox.csv`
- `inclusion_vs_authors_table_s9.csv`
- `overlap_with_internal_cox.csv` (may be empty)
