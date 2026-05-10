#!/usr/bin/env python3
"""
Fix misplaced stacked citations under Results heading in LIPID.docx.

User requirement:
  - Do not place multiple Supplementary Table citations as a block directly under the Results title.
  - Instead, cite each Supplementary Table locally inside the relevant Results subsection.
  - Do not group supplementary citations as ranges or grouped lists.

Actions:
  1) Delete the stacked paragraphs:
       "Supporting data for this analysis are provided in Supplementary Table X."
       "Relevant results are summarized in Table 3."
  2) Insert/ensure local, single-item citations near the most relevant paragraphs:
       - Supp Table 9: cluster assignments -> near clustering/subtype assignment paragraph
       - Supp Table 6: bootstrap stability -> near stability mention
       - Supp Table 5: cross-validation -> near subtype prediction performance paragraph
       - Supp Table 3: differential features -> near differential features mention
       - Supp Table 10: network statistics -> near network section
       - Supp Table 13: summary statistics -> near dataset overview / characteristics paragraph

Creates a timestamped .bak backup and edits DOCX in place.
"""

from __future__ import annotations

import datetime as _dt
import os
import shutil
from typing import Optional

from docx import Document
from docx.text.paragraph import Paragraph


DOCX_PATH = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration/LIPID.docx"

STACK_PREFIX = "Supporting data for this analysis are provided in Supplementary Table "
STACK_SENTENCES = {
    "Supporting data for this analysis are provided in Supplementary Table 13.",
    "Supporting data for this analysis are provided in Supplementary Table 10.",
    "Supporting data for this analysis are provided in Supplementary Table 9.",
    "Supporting data for this analysis are provided in Supplementary Table 6.",
    "Supporting data for this analysis are provided in Supplementary Table 5.",
    "Supporting data for this analysis are provided in Supplementary Table 3.",
    "Relevant results are summarized in Table 3.",
}


def backup(path: str) -> str:
    ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = f"{path}.{ts}.bak"
    shutil.copy2(path, bak)
    return bak


def delete_paragraph(p: Paragraph) -> None:
    el = p._element
    el.getparent().remove(el)
    p._p = p._element = None  # type: ignore[attr-defined]


def find_first_contains(doc: Document, needle: str, start: int = 0) -> Optional[int]:
    for i in range(start, len(doc.paragraphs)):
        if needle in (doc.paragraphs[i].text or ""):
            return i
    return None


def add_parenthetical_if_missing(p: Paragraph, cite: str) -> bool:
    """
    Append ' (Supplementary Table X).' if not already present.
    Returns True if modified.
    """
    t = p.text or ""
    if cite in t:
        return False
    # avoid double periods
    if t.strip().endswith("."):
        p.text = t.rstrip()[:-1] + f" ({cite})."
    else:
        p.text = t + f" ({cite})."
    return True


def ensure_sentence_after(doc: Document, anchor_idx: int, sentence: str) -> None:
    """
    Insert sentence as its own paragraph after anchor_idx if not already present anywhere.
    """
    if any((p.text or "").strip() == sentence for p in doc.paragraphs):
        return
    nxt = doc.paragraphs[anchor_idx + 1] if anchor_idx + 1 < len(doc.paragraphs) else None
    if nxt is not None:
        nxt.insert_paragraph_before(sentence)
    else:
        doc.add_paragraph(sentence)


def main() -> None:
    if not os.path.exists(DOCX_PATH):
        raise SystemExit(f"Missing: {DOCX_PATH}")

    bak = backup(DOCX_PATH)
    doc = Document(DOCX_PATH)

    # 1) Delete stacked paragraphs under Results
    removed = 0
    for p in list(doc.paragraphs):
        if (p.text or "").strip() in STACK_SENTENCES or (p.text or "").startswith(STACK_PREFIX):
            delete_paragraph(p)
            removed += 1

    # 2) Insert local citations
    # Supp Table 13 (summary statistics): near the very first Results paragraph describing dataset size/features
    r1 = find_first_contains(doc, "We analyzed integrated lipidomic and metabolomic profiles", 0)
    if r1 is not None:
        add_parenthetical_if_missing(doc.paragraphs[r1], "Supplementary Table 13")

    # Supp Table 9 (cluster assignments): near paragraph describing concordant subtype assignments / methods
    clust = find_first_contains(doc, "Across k-means, hierarchical clustering, and non-negative matrix factorization", 0)
    if clust is not None:
        add_parenthetical_if_missing(doc.paragraphs[clust], "Supplementary Table 9")

    # Supp Table 6 (bootstrap stability): same paragraph contains stability mention
    stab = find_first_contains(doc, "bootstrap", 0) or clust
    if stab is not None:
        add_parenthetical_if_missing(doc.paragraphs[stab], "Supplementary Table 6")

    # Supp Table 3 (differential features): near first explicit mention of 728 differential features (pathway enrichment section)
    diff = find_first_contains(doc, "352/728", 0) or find_first_contains(doc, "728", 0)
    if diff is not None:
        add_parenthetical_if_missing(doc.paragraphs[diff], "Supplementary Table 3")

    # Supp Table 10 (network statistics): near network construction/statistics paragraph
    net = find_first_contains(doc, "correlation network", 0) or find_first_contains(doc, "network", 0)
    if net is not None:
        add_parenthetical_if_missing(doc.paragraphs[net], "Supplementary Table 10")

    # Supp Table 5 (cross-validation): near subtype prediction performance section
    cv = find_first_contains(doc, "Cross-validation yielded", 0) or find_first_contains(doc, "repeated stratified cross-validation", 0)
    if cv is not None:
        add_parenthetical_if_missing(doc.paragraphs[cv], "Supplementary Table 5")

    doc.save(DOCX_PATH)
    print(f"Updated: {DOCX_PATH}")
    print(f"Backup:  {bak}")
    print(f"Removed stacked Results-citation paragraphs: {removed}")


if __name__ == "__main__":
    main()

