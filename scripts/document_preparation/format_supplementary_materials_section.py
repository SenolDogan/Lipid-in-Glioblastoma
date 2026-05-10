#!/usr/bin/env python3
"""
In LIPID.docx, in the "Supplementary Materials" section, ensure that under the
"Figures" and "Tables" headings we list ONLY:
  - "Supplementary Figure X. <Title>"
  - "Supplementary Table X. <Title>"

No long explanations or additional supplementary text should remain under these headings.

Source of titles: the Supplementary Figure/Table caption paragraphs already present in the doc.
Handles duplicate numbering (e.g., two "Supplementary Figure 3") by converting to 3A/3B.
Creates a timestamped backup before saving.
"""

from __future__ import annotations

import datetime as _dt
import re
import shutil
from collections import defaultdict

from docx import Document
from docx.text.paragraph import Paragraph
from docx.oxml import OxmlElement


DOCX_PATH = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration/LIPID.docx"

SUPP_MAT_HEADING = "supplementary materials"

SUPP_FIG_CAP_RE = re.compile(r"^Supplementary\s+Figure\s+(\d+)\.\s*(.*)$", re.IGNORECASE)
SUPP_TAB_CAP_RE = re.compile(r"^Supplementary\s+Table\s+(\d+)\.\s*(.*)$", re.IGNORECASE)


def backup(path: str) -> str:
    ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = f"{path}.{ts}.bak"
    shutil.copy2(path, bak)
    return bak


def remove_paragraph(p: Paragraph) -> None:
    elm = p._p  # type: ignore[attr-defined]
    parent = elm.getparent()
    if parent is not None:
        parent.remove(elm)


def insert_paragraph_after(p: Paragraph, text: str) -> Paragraph:
    new_p_elm = OxmlElement("w:p")
    p._p.addnext(new_p_elm)  # type: ignore[attr-defined]
    new_p = Paragraph(new_p_elm, p._parent)  # type: ignore[arg-type]
    new_p.text = text
    return new_p


def collect_supplementary_titles(doc: Document) -> tuple[list[str], list[str]]:
    """
    Returns:
      - fig_lines: ["Supplementary Figure 1. Title", "Supplementary Figure 3A. ...", ...]
      - tab_lines: ["Supplementary Table 1. Title", ...]
    """
    fig_caps: list[tuple[int, str]] = []
    tab_caps: list[tuple[int, str]] = []

    for p in doc.paragraphs:
        t = (p.text or "").strip()
        m = SUPP_FIG_CAP_RE.match(t)
        if m:
            fig_caps.append((int(m.group(1)), (m.group(2) or "").strip()))
            continue
        m = SUPP_TAB_CAP_RE.match(t)
        if m:
            tab_caps.append((int(m.group(1)), (m.group(2) or "").strip()))

    # Sort by number but preserve within-number order for duplicates based on appearance
    fig_ordered = sorted(enumerate(fig_caps), key=lambda x: (x[1][0], x[0]))
    fig_counts = defaultdict(int)
    fig_total_per_n = defaultdict(int)
    for n, _title in fig_caps:
        fig_total_per_n[n] += 1

    fig_lines: list[str] = []
    for _idx, (n, title) in fig_ordered:
        fig_counts[n] += 1
        suffix = ""
        if fig_total_per_n[n] > 1:
            suffix = chr(ord("A") + fig_counts[n] - 1)
        label = f"Supplementary Figure {n}{suffix}."
        if title:
            fig_lines.append(f"{label} {title}")
        else:
            fig_lines.append(f"{label}")

    tab_caps_sorted = sorted(tab_caps, key=lambda x: x[0])
    tab_lines: list[str] = []
    for n, title in tab_caps_sorted:
        label = f"Supplementary Table {n}."
        tab_lines.append(f"{label} {title}".strip())

    return fig_lines, tab_lines


def find_section_idx(doc: Document, text: str) -> int | None:
    target = text.strip().lower()
    for i, p in enumerate(doc.paragraphs):
        if (p.text or "").strip().lower() == target:
            return i
    return None


def find_next_heading_idx(doc: Document, start_idx: int) -> int:
    """
    Find the next *top-level* section heading after start_idx.
    Heuristic: style name equals 'Heading 1' OR equals common major headings.
    """
    major = {"references", "discussion", "results", "materials and methods", "materials & methods"}
    for i in range(start_idx + 1, len(doc.paragraphs)):
        t = (doc.paragraphs[i].text or "").strip()
        if not t:
            continue
        style = (doc.paragraphs[i].style.name or "") if doc.paragraphs[i].style else ""
        if style.strip().lower() == "heading 1":
            return i
        if t.strip().lower() in major:
            return i
    return len(doc.paragraphs)


def replace_list_under_heading(doc: Document, heading_idx: int, lines: list[str], stop_idx: int) -> None:
    """
    Replace everything between heading_idx (exclusive) and stop_idx (exclusive)
    with paragraphs containing `lines`.
    """
    # Snapshot paragraphs to remove (indices will shift as we delete)
    to_remove = list(doc.paragraphs[heading_idx + 1 : stop_idx])
    for p in reversed(to_remove):
        remove_paragraph(p)

    # Insert new lines after heading, preserving order
    insert_after = doc.paragraphs[heading_idx]
    for line in lines:
        insert_after = insert_paragraph_after(insert_after, line)


def main() -> None:
    bak = backup(DOCX_PATH)
    doc = Document(DOCX_PATH)

    fig_lines, tab_lines = collect_supplementary_titles(doc)

    supp_idx = find_section_idx(doc, "Supplementary Materials") or find_section_idx(doc, "Supplementary material")
    if supp_idx is None:
        raise RuntimeError("Could not find 'Supplementary Materials' heading in document.")

    section_end = find_next_heading_idx(doc, supp_idx)

    # Find Supplementary Figures and Supplementary Tables subheadings within Supplementary Materials section
    fig_head = None
    tab_head = None
    for i in range(supp_idx + 1, section_end):
        t = (doc.paragraphs[i].text or "").strip().lower()
        if t in {"supplementary figures", "figures"}:
            fig_head = i
        elif t in {"supplementary tables", "tables"}:
            tab_head = i

    if fig_head is None or tab_head is None:
        raise RuntimeError("Could not find both 'Supplementary Figures' and 'Supplementary Tables' headings under Supplementary Materials.")

    # Replace later blocks first so earlier deletions don't invalidate indices.
    # Determine stop indices: each block stops at the next subheading or section end.
    blocks = []
    if fig_head < tab_head:
        blocks.append((tab_head, tab_lines, section_end))
        blocks.append((fig_head, fig_lines, tab_head))
    else:
        blocks.append((fig_head, fig_lines, section_end))
        blocks.append((tab_head, tab_lines, fig_head))

    for head_idx, lines, stop_idx in sorted(blocks, key=lambda x: x[0], reverse=True):
        replace_list_under_heading(doc, head_idx, lines, stop_idx)

    doc.save(DOCX_PATH)
    print(f"Updated: {DOCX_PATH}")
    print(f"Backup:  {bak}")


if __name__ == "__main__":
    main()

