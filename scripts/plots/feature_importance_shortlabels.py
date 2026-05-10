#!/usr/bin/env python3
"""
Create feature importance barplot with shortened lipid names for readability.
Reads results/feature_importance_rf.csv and writes plots/feature_importance_shortlabels.png
"""

from __future__ import annotations

import re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager


IN_CSV = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration/results/feature_importance_rf.csv"
OUT_PNG = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration/plots/feature_importance_shortlabels.png"
TIMES_FONT = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"


def shorten(name: str) -> str:
    s = str(name)
    # Keep first entry if multiple annotations separated by ';'
    if ";" in s:
        s = s.split(";")[0]
    # Remove duplicate markers / suffixes like _A, _dup#
    s = re.sub(r"_dup\d+$", "", s)
    s = re.sub(r"_A$", "", s)
    return s


def main() -> None:
    font_manager.fontManager.addfont(TIMES_FONT)
    plt.rcParams["font.family"] = "Times New Roman"

    df = pd.read_csv(IN_CSV)
    df = df.dropna(subset=["feature", "importance"]).copy()
    top = df.head(30).iloc[::-1].copy()  # reverse for barh
    top["short"] = top["feature"].map(shorten)

    plt.figure(figsize=(10.5, 11.5))
    plt.barh(range(len(top)), top["importance"].values, color="#ff6b81")
    plt.yticks(range(len(top)), top["short"].values, fontsize=14, fontweight="bold")
    plt.xlabel("Importance", fontsize=13, fontweight="bold")
    plt.title("Top 30 Features by Random Forest Importance", fontsize=16, fontweight="bold")
    plt.grid(axis="x", alpha=0.25)
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=600)
    plt.close()
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()

