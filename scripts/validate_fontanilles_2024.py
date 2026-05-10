#!/usr/bin/env python3
"""
External validation using Fontanilles et al. 2024 (s40478-024-01861-5) supplementary XLSX.

Goal
- Provide a real external dataset with quantified metabolomics/lipid-like metabolites + overall survival.
- Run basic integrity checks.
- Recompute univariate survival associations at Inclusion timepoint and compare to authors' Table_S9.
- Optionally compare overlapping feature names against our internal Cox results (if any overlap).

Inputs
- external_validation/Fontanilles_2024_s40478-024-01861-5_MOESM1.xlsx (downloaded already)

Outputs
- external_validation/fontanilles_2024/_validation/report.md
- external_validation/fontanilles_2024/_validation/inclusion_univariate_cox.csv
- external_validation/fontanilles_2024/_validation/inclusion_vs_authors_table_s9.csv
- external_validation/fontanilles_2024/_validation/overlap_with_internal_cox.csv
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
XLSX_PATH = HERE / "Fontanilles_2024_s40478-024-01861-5_MOESM1.xlsx"
OUT_DIR = HERE / "fontanilles_2024" / "_validation"

# Internal results (optional)
INTERNAL_COX = (
    Path(__file__).resolve().parents[1]
    / "Lipidomics_Metabolomics_Integration"
    / "results"
    / "cox_regression_results.csv"
)


def _norm_feature(x: str) -> str:
    x = (x or "").strip()
    # normalize unicode spaces and casing, keep punctuation used in lipid notation
    x = re.sub(r"\s+", " ", x)
    return x


def _patient_id_from_sample(sample: str) -> str:
    # e.g. 01_006_Inclusion -> 01_006
    parts = (sample or "").split("_")
    if len(parts) >= 2:
        return "_".join(parts[:2])
    return sample


def _timepoint_from_sample(sample: str) -> str:
    # e.g. 01_006_Inclusion -> Inclusion
    if sample and "_" in sample:
        return sample.split("_")[-1]
    return ""


@dataclass(frozen=True)
class SurvivalSpec:
    time_col: str
    event_col: str


def _detect_survival_cols(clin: pd.DataFrame) -> SurvivalSpec:
    cols = {c.lower(): c for c in clin.columns}
    # The sheet uses "OS DELAY" and "OS" (Yes/No-like) in the provided preview.
    time_col = cols.get("os delay") or cols.get("os_delay") or cols.get("os_delay_days") or cols.get("os_days")
    event_col = cols.get("os") or cols.get("death") or cols.get("status")
    if not time_col or not event_col:
        raise ValueError(f"Could not detect survival columns. Found columns include: {list(clin.columns)[:30]}")
    return SurvivalSpec(time_col=time_col, event_col=event_col)


def _event_to_int(v) -> int | float:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    s = str(v).strip().lower()
    if s in {"1", "dead", "deceased", "yes", "y", "true"}:
        return 1
    if s in {"0", "alive", "no", "n", "false"}:
        return 0
    return np.nan


def _cox_median_split(time: pd.Series, event: pd.Series, x: pd.Series) -> dict:
    """
    Simple 2-group Cox equivalent via log-rank approximation is possible,
    but we stick to CoxPH on a binary covariate for direct HR.
    """
    from lifelines import CoxPHFitter

    df = pd.DataFrame({"time": time, "event": event, "x": x})
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    if df["event"].sum() < 3:
        return {"n": len(df), "events": float(df["event"].sum()), "hr": np.nan, "p": np.nan}

    cph = CoxPHFitter()
    try:
        cph.fit(df, duration_col="time", event_col="event")
        hr = float(np.exp(cph.params_["x"]))
        p = float(cph.summary.loc["x", "p"])
        return {"n": len(df), "events": float(df["event"].sum()), "hr": hr, "p": p}
    except Exception:
        return {"n": len(df), "events": float(df["event"].sum()), "hr": np.nan, "p": np.nan}


def main() -> int:
    if not XLSX_PATH.exists():
        raise FileNotFoundError(f"Missing XLSX: {XLSX_PATH}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    xl = pd.ExcelFile(XLSX_PATH)
    required = {"Table_S1", "Table_S3", "Table_S9"}
    missing = required - set(xl.sheet_names)
    if missing:
        raise ValueError(f"Missing required sheets: {sorted(missing)}. Available: {xl.sheet_names}")

    clin = xl.parse("Table_S1")
    matrix = xl.parse("Table_S3")
    authors_s9 = xl.parse("Table_S9")

    # Clinical checks
    if "sample" not in [c.lower() for c in clin.columns]:
        raise ValueError("Clinical Table_S1 missing 'sample' column.")
    sample_col = [c for c in clin.columns if c.lower() == "sample"][0]
    clin = clin.copy()
    clin["patient_id"] = clin[sample_col].map(_patient_id_from_sample)
    clin["timepoint"] = clin[sample_col].map(_timepoint_from_sample)

    surv = _detect_survival_cols(clin)
    clin["os_time"] = pd.to_numeric(clin[surv.time_col], errors="coerce")
    clin["os_event"] = clin[surv.event_col].map(_event_to_int)

    # Use Inclusion timepoint for “baseline” survival analysis
    clin_incl = clin[clin["timepoint"].str.lower() == "inclusion"].copy()
    clin_incl = clin_incl.dropna(subset=["os_time", "os_event"])

    # Matrix: keep metabolomics only and Inclusion samples only
    if "omics_type" not in matrix.columns:
        raise ValueError("Matrix Table_S3 missing 'omics_type' column.")
    if "variable" not in matrix.columns:
        raise ValueError("Matrix Table_S3 missing 'variable' column.")

    sample_cols = [c for c in matrix.columns if re.match(r"^\d{2}_\d{3}_(Inclusion|W6|Progression)$", str(c))]
    if not sample_cols:
        # fallback: treat any column containing '_Inclusion' as sample
        sample_cols = [c for c in matrix.columns if str(c).endswith("_Inclusion")]

    matrix = matrix.copy()
    matrix["variable_norm"] = matrix["variable"].map(_norm_feature)

    metab = matrix[matrix["omics_type"].astype(str).str.lower() == "metabolomics"].copy()
    incl_cols = [c for c in sample_cols if str(c).endswith("_Inclusion")]
    metab_long = (
        metab[["variable_norm"] + incl_cols]
        .melt(id_vars=["variable_norm"], var_name="sample", value_name="value")
    )
    metab_long["patient_id"] = metab_long["sample"].map(_patient_id_from_sample)
    metab_long["value"] = pd.to_numeric(metab_long["value"], errors="coerce")

    # Merge with survival
    surv_map = clin_incl.set_index("patient_id")[["os_time", "os_event"]]
    metab_long = metab_long.join(surv_map, on="patient_id", how="inner")
    metab_long = metab_long.dropna(subset=["value", "os_time", "os_event"])

    # Compute median split Cox per metabolite
    rows = []
    for var, g in metab_long.groupby("variable_norm", sort=True):
        med = g["value"].median()
        x = (g["value"] >= med).astype(int)  # 1=High
        res = _cox_median_split(time=g["os_time"], event=g["os_event"], x=x)
        rows.append(
            {
                "variable": var,
                "median": float(med) if pd.notna(med) else np.nan,
                **res,
            }
        )
    out_cox = pd.DataFrame(rows).sort_values(["p", "variable"], na_position="last")
    out_cox.to_csv(OUT_DIR / "inclusion_univariate_cox.csv", index=False)

    # Compare with authors Table_S9 (Inclusion)
    a = authors_s9.copy()
    a.columns = [str(c).strip() for c in a.columns]
    a_incl = a[a["time"].astype(str).str.lower() == "inclusion"].copy() if "time" in a.columns else pd.DataFrame()
    if not a_incl.empty and "variable" in a_incl.columns:
        a_incl["variable_norm"] = a_incl["variable"].map(_norm_feature)
        # Authors Table_S9 is in two rows (High/Low) per variable; HR is typically on the "Low" row.
        label_col = "label" if "label" in a_incl.columns else None
        if label_col:
            a_low = a_incl[a_incl[label_col].astype(str).str.contains("low", case=False, na=False)].copy()
        else:
            a_low = a_incl.copy()

        a_low["HR"] = pd.to_numeric(a_low["HR"], errors="coerce")
        a_low["p.value"] = pd.to_numeric(a_low.get("p.value"), errors="coerce")
        a_hr = (
            a_low.dropna(subset=["HR"])
            .groupby("variable_norm", as_index=False)
            .agg(authors_hr=("HR", "first"), authors_p=("p.value", "first"))
        )
        merged = out_cox.merge(a_hr, left_on="variable", right_on="variable_norm", how="inner")
        merged["hr"] = pd.to_numeric(merged["hr"], errors="coerce")
        merged["hr_direction_match"] = np.sign(np.log(merged["hr"])) == np.sign(np.log(merged["authors_hr"]))
        merged = merged.drop(columns=["variable_norm"])
        merged.to_csv(OUT_DIR / "inclusion_vs_authors_table_s9.csv", index=False)
    else:
        merged = pd.DataFrame()

    # Optional overlap with our internal Cox features (likely limited, but we record it)
    overlap_df = pd.DataFrame()
    if INTERNAL_COX.exists():
        internal = pd.read_csv(INTERNAL_COX)
        # internal format: covariate, coef, exp(coef) ...
        if "covariate" in internal.columns and "exp(coef)" in internal.columns:
            internal["feature_norm"] = internal["covariate"].map(_norm_feature)
            # try exact overlaps against metabolomics variable names
            overlap = out_cox.merge(
                internal[["feature_norm", "exp(coef)", "p"]].rename(
                    columns={"feature_norm": "variable", "exp(coef)": "internal_hr", "p": "internal_p"}
                ),
                on="variable",
                how="inner",
            )
            if not overlap.empty:
                overlap["hr_direction_match"] = np.sign(np.log(overlap["hr"])) == np.sign(np.log(overlap["internal_hr"]))
            overlap_df = overlap.sort_values(["p", "variable"], na_position="last")
            overlap_df.to_csv(OUT_DIR / "overlap_with_internal_cox.csv", index=False)

    # Write report
    def pct(x: float) -> str:
        return f"{x*100:.1f}%"

    n_pat = clin_incl["patient_id"].nunique()
    n_incl_samples = clin_incl.shape[0]
    n_vars = out_cox.shape[0]
    n_sig = int((out_cox["p"] < 0.05).sum()) if "p" in out_cox.columns else 0
    n_comp = int(merged.shape[0]) if isinstance(merged, pd.DataFrame) else 0
    n_match = int(merged["hr_direction_match"].sum()) if n_comp and "hr_direction_match" in merged.columns else 0

    lines = []
    lines.append("## External validation: Fontanilles et al. 2024 (Acta Neuropathologica Communications)")
    lines.append("")
    lines.append(f"- **Source XLSX**: `{XLSX_PATH.relative_to(Path.cwd())}`")
    lines.append(f"- **Outputs**: `{OUT_DIR.relative_to(Path.cwd())}`")
    lines.append("")
    lines.append("### Data availability (this dataset works for validation A+B)")
    lines.append("- **Has quantified metabolomics matrix**: yes (`Table_S3`, Metabolomics rows)")
    lines.append("- **Has overall survival**: yes (`Table_S1`: `OS DELAY` + `OS`)")
    lines.append("")
    lines.append("### QC summary (Inclusion timepoint)")
    lines.append(f"- **Patients with Inclusion + OS**: {n_pat} patients ({n_incl_samples} inclusion rows)")
    lines.append(f"- **Metabolomics variables tested (median-split Cox)**: {n_vars}")
    lines.append(f"- **Nominally significant (p<0.05)**: {n_sig} ({pct(n_sig/max(n_vars,1))})")
    lines.append("")
    lines.append("### Reproducibility vs authors’ Supplementary Table_S9")
    if n_comp:
        lines.append(f"- **Variables comparable (intersection)**: {n_comp}")
        lines.append(f"- **HR direction match**: {n_match}/{n_comp} ({pct(n_match/max(n_comp,1))})")
        lines.append(f"- **File**: `inclusion_vs_authors_table_s9.csv`")
    else:
        lines.append("- **Comparable intersection**: not computed (Table_S9 format mismatch or no overlap detected).")
    lines.append("")
    if INTERNAL_COX.exists():
        lines.append("### Overlap vs our internal Cox features")
        if not overlap_df.empty:
            lines.append(f"- **Exact-name overlap count**: {len(overlap_df)}")
            lines.append("- **File**: `overlap_with_internal_cox.csv`")
        else:
            lines.append("- **Exact-name overlap count**: 0 (naming/platform differences expected).")
    lines.append("")
    lines.append("### Artifacts")
    lines.append(f"- `inclusion_univariate_cox.csv`")
    if n_comp:
        lines.append(f"- `inclusion_vs_authors_table_s9.csv`")
    if INTERNAL_COX.exists():
        lines.append(f"- `overlap_with_internal_cox.csv` (may be empty)")
    lines.append("")
    (OUT_DIR / "report.md").write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")

    print(f"Wrote: {OUT_DIR / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

