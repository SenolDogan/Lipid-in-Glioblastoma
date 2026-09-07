#!/usr/bin/env python3
"""
Quick release verification for reviewers (~30 seconds).

Run from the repository root:
    python scripts/verify_release.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
SUPP = ROOT / "Supplementary"


def check(name: str, ok: bool, detail: str = "") -> bool:
    status = "PASS" if ok else "FAIL"
    suffix = f" — {detail}" if detail else ""
    print(f"[{status}] {name}{suffix}")
    return ok


def main() -> int:
    ok = True

    manifest = SUPP / "tables" / "Supplementary_Table_S17.csv"
    ok &= check("Supplementary Table S17 manifest exists", manifest.is_file())

    ca_path = RESULTS / "cluster_assignments.csv"
    ok &= check("cluster_assignments.csv exists", ca_path.is_file())
    if ca_path.is_file():
        ca = pd.read_csv(ca_path)
        ok &= check("cluster_assignments row count", len(ca) == 74, f"n={len(ca)}")
        ok &= check(
            "cluster_assignments columns",
            {"sample_id", "kmeans_subtype", "hierarchical_subtype", "nmf_subtype", "vital_status"}.issubset(
                ca.columns
            ),
        )
        ok &= check(
            "nmf_subtype label counts",
            ca["nmf_subtype"].value_counts().to_dict() == {0: 62, 1: 12},
            "expect {0: 62, 1: 12}",
        )

    corr_path = RESULTS / "lipid_metabolite_correlations.csv"
    ok &= check("lipid_metabolite_correlations.csv exists", corr_path.is_file())
    if corr_path.is_file():
        corr = pd.read_csv(corr_path)
        ok &= check("correlation atlas row count", len(corr) == 1676, f"n={len(corr)}")
        ok &= check(
            "correlation atlas columns",
            {"feature1", "feature2", "correlation", "pvalue", "fdr"}.issubset(corr.columns),
        )
        if len(corr):
            within = (corr["correlation"].abs() >= 0.3) & (corr["fdr"] < 0.05)
            ok &= check(
                "correlation atlas thresholds",
                bool(within.all()),
                f"{int(within.sum())}/{len(corr)} rows meet |rho|>=0.3 and FDR<0.05",
            )

    perm_path = RESULTS / "acceptance_network_permutation_summary.csv"
    ok &= check("network permutation summary exists", perm_path.is_file())
    if perm_path.is_file():
        perm = pd.read_csv(perm_path)
        if "empirical_p_perm_ge_obs" in perm.columns and len(perm):
            pval = float(perm["empirical_p_perm_ge_obs"].iloc[0])
            ok &= check("network permutation empirical p-value", pval < 0.05, f"p={pval:.4f}")

    for rel in (
        "differential_features.csv",
        "feature_importance_rf.csv",
        "hub_molecule_analysis.csv",
        "lipid_class_abundance.csv",
        "metabolite_category_abundance.csv",
    ):
        ok &= check(f"results/{rel} exists", (RESULTS / rel).is_file())

    fig_dir = SUPP / "figures"
    if fig_dir.is_dir():
        pngs = sorted(fig_dir.glob("Supplementary_Figure_S*.png"))
        ok &= check("supplementary figure PNG count", len(pngs) == 8, f"found {len(pngs)}")
    else:
        ok &= check("Supplementary/figures directory exists", False)

    mapping = ROOT / "Feature_Names_Mapping.xlsx"
    ok &= check("Feature_Names_Mapping.xlsx exists", mapping.is_file())

    if ok:
        print("\nAll checks passed. Release is consistent with the manuscript.")
        return 0

    print("\nOne or more checks failed. See messages above.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
