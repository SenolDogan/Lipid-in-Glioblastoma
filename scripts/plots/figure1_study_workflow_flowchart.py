#!/usr/bin/env python3
"""
Generate manuscript-ready Figure 1 flowchart as a PNG.
No manuscript files are modified.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # headless backend for reliable file rendering
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


OUT_PNG = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration/plots/figure1_study_workflow_flowchart.png"
OUT_DPI = 600  # higher-resolution output for manuscript figures
FONT_FAMILY = "Times New Roman"
TEXT_SCALE = 1.25
TEXT_OFFSET = 5.0


def add_box(ax, x, y, w, h, text, fontsize=10.5, lw=1.1, fc="white", ec="#2b2b2b"):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.004,rounding_size=0.018",
        linewidth=lw,
        facecolor=fc,
        edgecolor=ec,
        mutation_aspect=1.0,
    )
    ax.add_patch(box)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize * TEXT_SCALE + TEXT_OFFSET,
        color="#111111",
        family=FONT_FAMILY,
        wrap=True,
        linespacing=0.75,
    )
    return box


def add_label(ax, x, y, label, fontsize=10):
    ax.text(
        x,
        y,
        label,
        ha="left",
        va="top",
        fontsize=fontsize * TEXT_SCALE + TEXT_OFFSET,
        color="#111111",
        family=FONT_FAMILY,
        weight="bold",
    )


def add_arrow(ax, x1, y1, x2, y2, lw=1.8, color="#1f1f1f", shrink=14):
    arr = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="-|>",
        mutation_scale=24,
        linewidth=lw + 0.7,
        color=color,
        shrinkA=shrink,
        shrinkB=shrink,
        connectionstyle="arc3,rad=0.0",
    )
    ax.add_patch(arr)


def main():
    # Taller canvas so enlarged boxes and arrows do not get compressed in Word.
    fig = plt.figure(figsize=(11.7, 14.8), dpi=OUT_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    # Add small margins to prevent stroke clipping at edges
    ax.set_xlim(-0.07, 1.07)
    # Tight vertical crop: keep only a small margin beyond the first/last box borders.
    ax.set_ylim(-3.32, 1.01)
    ax.axis("off")

    # Layout coordinates
    w_main, h_main = 0.74, 0.32
    x_main = 0.13
    # Slightly increase spacing between the first boxes for clearer arrows
    # Slightly increase vertical gaps to lengthen downward arrows
    yA, yB, yC, yD = 0.67, 0.17, -0.33, -0.83

    add_box(
        ax,
        x_main,
        yA,
        w_main,
        h_main,
        "Data sources\nCPTAC 2021 lipidome + metabolome\n+ TCGA GBM clinical annotations",
        fontsize=13,
    )
    add_box(
        ax,
        x_main,
        yB,
        w_main,
        h_main,
        "Cohort assembly\nMatch patient IDs across modalities\nRetain complete cases",
        fontsize=13,
    )
    add_box(
        ax,
        x_main,
        yC,
        w_main,
        h_main,
        "Preprocessing\nRemove features >50% missing\nMedian imputation\nScaling/normalization",
        fontsize=13,
    )
    add_box(
        ax,
        x_main,
        yD,
        w_main,
        h_main,
        "Integrated feature matrix\nLipids + metabolites (~700 features)",
        fontsize=13,
    )

    # Arrows main vertical
    add_arrow(ax, 0.5, yA, 0.5, yB + h_main, shrink=2)
    add_arrow(ax, 0.5, yB, 0.5, yC + h_main, shrink=2)
    add_arrow(ax, 0.5, yC, 0.5, yD + h_main, shrink=2)

    # Branches
    # Increase vertical gap between Box 4 and the A/B/C row
    y_branch_top = -1.56
    w_b, h_b = 0.325, 0.54
    # Glue A/B/C boxes closer together
    gap = 0.045
    xL = -0.005
    xM = xL + w_b + gap
    xR = xM + w_b + gap

    box5L = add_box(
        ax,
        xL,
        y_branch_top,
        w_b,
        h_b,
        "Subtype discovery\nUnsupervised clustering\n(K-means / Hierarchical / NMF)\nSelect k via silhouette & elbow",
        fontsize=11.6,
    )
    add_label(ax, xL + 0.012, y_branch_top + h_b - 0.01, "A", fontsize=11)

    box5M = add_box(
        ax,
        xM,
        y_branch_top,
        w_b,
        h_b,
        "Lipid–metabolite associations\nSpearman correlation\nBH-FDR control",
        fontsize=11.6,
    )
    add_label(ax, xM + 0.012, y_branch_top + h_b - 0.01, "B", fontsize=11)

    box5R = add_box(
        ax,
        xR,
        y_branch_top,
        w_b,
        h_b,
        "Differential & pathway analysis\nMann–Whitney U + BH-FDR\nEffect size (Cohen's d)\nClass/category enrichment",
        fontsize=11.6,
    )
    add_label(ax, xR + 0.012, y_branch_top + h_b - 0.01, "C", fontsize=11)

    # Left branch second box (subtypes)
    # Increase vertical gap between row 5 and row 6
    y6L = -2.14
    h6 = 0.40
    box6L = add_box(
        ax,
        xL,
        y6L,
        w_b,
        h6,
        "Metabolic–lipidomic subtypes\nSubtype 0 (TG-high)\nSubtype 1 (SM/Cer-high)",
        fontsize=11.6,
    )

    # Middle branch second box (network)
    y6M = y6L
    box6M = add_box(
        ax,
        xM,
        y6M,
        w_b,
        h6,
        "Network analysis\nConstruct correlation network\nHub prioritization (centrality)",
        fontsize=11.6,
    )

    # Connect integrated matrix to branches (diagonal arrows)
    # Increase gap between arrow tips and A/B/C boxes
    yD_bottom = yD
    # Diagonal arrows: add visible gap to box borders
    add_arrow(ax, 0.5, yD_bottom, xL + w_b / 2, y_branch_top + h_b, shrink=6)
    add_arrow(ax, xM + w_b / 2, yD_bottom, xM + w_b / 2, y_branch_top + h_b, shrink=2)
    add_arrow(ax, 0.5, yD_bottom, xR + w_b / 2, y_branch_top + h_b, shrink=6)

    # Connect 5L->6L, 5M->6M
    add_arrow(ax, xL + w_b / 2, y_branch_top, xL + w_b / 2, y6L + h6, shrink=6)
    add_arrow(ax, xM + w_b / 2, y_branch_top, xM + w_b / 2, y6M + h6, shrink=6)

    # Merge to clinical/prediction
    # Increase vertical gap between the 6th-row boxes and the Clinical box
    w7, h7 = 0.74, 0.48
    x7, y7 = 0.13, -2.80
    add_box(
        ax,
        x7,
        y7,
        w7,
        h7,
        "Clinical association & prediction\nKaplan–Meier / Cox\nCross-validated subtype prediction\nMinimal panel feasibility",
        fontsize=12.0,
    )

    # Merge arrows from 6L, 6M, 5R
    add_arrow(ax, xL + w_b / 2, y6L, xL + w_b / 2, y7 + h7, shrink=8)
    # Network analysis -> Clinical should point straight down
    add_arrow(ax, xM + w_b / 2, y6M, xM + w_b / 2, y7 + h7, shrink=8)
    add_arrow(ax, xR + w_b / 2, y_branch_top, 0.65, y7 + h7, shrink=8)

    # Final interpretation box
    # Increase vertical gap between Clinical box (row 7) and Interpretation box (row 8)
    w8, h8 = 0.74, 0.32
    x8, y8 = 0.13, -3.30
    add_box(
        ax,
        x8,
        y8,
        w8,
        h8,
        "Biological interpretation\nSubtype-specific metabolic programs\nDruggable pathway hypotheses",
        fontsize=12.2,
    )
    add_arrow(ax, 0.5, y7, 0.5, y8 + h8, shrink=1)

    fig.savefig(OUT_PNG, dpi=OUT_DPI, facecolor="white")
    plt.close(fig)
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()

