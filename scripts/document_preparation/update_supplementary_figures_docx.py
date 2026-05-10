#!/usr/bin/env python3
from __future__ import annotations

import datetime as _dt
import os
import re
import shutil
from typing import Optional

from docx import Document
from docx.oxml import OxmlElement
from docx.text.paragraph import Paragraph
from docx.shared import Inches


SUPP_DOCX = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration/Supplementary/Supplementary_Figures.docx"

PLOTS_DIR = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration/plots"

FIGS = {
    "1": {
        "img": os.path.join(PLOTS_DIR, "pca_analysis.png"),
        "caption": (
            "Supplementary Figure 1. Principal Component Analysis. PCA plot demonstrating separation between subtypes "
            "along PC1 and PC2; the first 5 principal components explain 50.1% of total variance."
        ),
    },
    "2": {
        "img": os.path.join(PLOTS_DIR, "figure5_dual_purpose_candidates.png"),
        "caption": (
            "Supplementary Figure 2. Dual-purpose candidate lipids prioritized by combining subtype discrimination "
            "(effect size) with network centrality. Shown are representative candidates that are both subtype-separating "
            "and network-central."
        ),
    },
    "3": {
        "img": os.path.join(PLOTS_DIR, "figure5_combined.png"),
        "caption": (
            "Supplementary Figure 3. Combined visualization of key subtype-separating and network-based prioritization "
            "outputs used for candidate nomination and downstream re-use."
        ),
    },
    "4": {
        "img": os.path.join(PLOTS_DIR, "survival_feature_importance.png"),
        "caption": (
            "Supplementary Figure 4. Feature importance for survival time prediction (Random Forest regression). "
            "Shown are the top features contributing to survival time prediction in exploratory analyses."
        ),
    },
    "5": {
        "img": os.path.join(PLOTS_DIR, "feature_importance_shortlabels.png"),
        "caption": (
            "Supplementary Figure 5. Random Forest feature importance for metabolic-lipidomic subtype prediction. "
            "Top features ranked by importance in a Random Forest classification model trained on integrated lipidomic "
            "and metabolomic profiles (n=74). Features are ranked by their contribution to subtype classification, with "
            "higher importance scores indicating greater discriminatory power between Subtype 0 and Subtype 1."
        ),
    },
}

CAP_RE = re.compile(r"^Supplementary Figure\s+(\d+)\.", flags=re.IGNORECASE)


def backup(path: str) -> str:
    ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = f"{path}.{ts}.bak"
    shutil.copy2(path, bak)
    return bak


def has_drawing(p: Paragraph) -> bool:
    return bool(p._p.xpath(".//w:drawing"))  # type: ignore[attr-defined]


def insert_paragraph_after(p: Paragraph) -> Paragraph:
    new_p_elm = OxmlElement("w:p")
    p._p.addnext(new_p_elm)  # type: ignore[attr-defined]
    return Paragraph(new_p_elm, p._parent)  # type: ignore[arg-type]


def delete_paragraph(paragraph: Paragraph) -> None:
    el = paragraph._element
    el.getparent().remove(el)
    paragraph._p = paragraph._element = None  # type: ignore[attr-defined]


def find_caption_idx(doc: Document, fig_id: str) -> Optional[int]:
    prefix = f"Supplementary Figure {fig_id}."
    for i, p in enumerate(doc.paragraphs):
        if (p.text or "").strip().startswith(prefix):
            return i
    return None


def delete_following_figure_block(doc: Document, cap_idx: int) -> None:
    """
    Remove drawing paragraphs and blank spacers immediately after a caption,
    stopping before the next supplementary caption.
    """
    i = cap_idx + 1
    while i < len(doc.paragraphs):
        t = (doc.paragraphs[i].text or "").strip()
        if CAP_RE.match(t):
            break
        if has_drawing(doc.paragraphs[i]) or t == "":
            delete_paragraph(doc.paragraphs[i])
            continue
        break


def insert_image_below_caption(doc: Document, cap_idx: int, img_path: str) -> None:
    pic_p = insert_paragraph_after(doc.paragraphs[cap_idx])
    run = pic_p.add_run()
    run.add_picture(img_path, width=Inches(6.8))


def ensure_figure(doc: Document, fig_id: str) -> None:
    fig = FIGS[fig_id]
    cap_idx = find_caption_idx(doc, fig_id)
    if cap_idx is None:
        # Append at end if missing
        if doc.paragraphs and (doc.paragraphs[-1].text or "").strip() != "":
            doc.add_paragraph("")
        doc.add_paragraph(fig["caption"])
        cap_idx = len(doc.paragraphs) - 1
    else:
        doc.paragraphs[cap_idx].text = fig["caption"]

    delete_following_figure_block(doc, cap_idx)
    insert_image_below_caption(doc, cap_idx, fig["img"])


def main() -> None:
    if not os.path.exists(SUPP_DOCX):
        raise SystemExit(f"Missing docx: {SUPP_DOCX}")
    for fig_id, info in FIGS.items():
        if not os.path.exists(info["img"]):
            raise SystemExit(f"Missing figure image for Supp Fig {fig_id}: {info['img']}")

    bak = backup(SUPP_DOCX)
    doc = Document(SUPP_DOCX)

    for fig_id in ["1", "2", "3", "4", "5"]:
        ensure_figure(doc, fig_id)

    doc.save(SUPP_DOCX)
    print(f"Updated: {SUPP_DOCX}")
    print(f"Backup : {bak}")


if __name__ == "__main__":
    main()

