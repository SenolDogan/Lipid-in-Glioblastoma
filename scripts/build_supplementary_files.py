#!/usr/bin/env python3
"""
Create supplementary deliverables inside `Supplementary/`:

1) `Supplementary_Figures.docx`
   - Contains all Supplementary Figures extracted from `LIPID.docx`
   - Each figure: caption/legend text + embedded image(s) immediately following the caption

2) `Supplementary_Tables.xlsx`
   - Contains all Supplementary Tables extracted from `LIPID.docx`
   - One worksheet per supplementary table (SuppTable1, SuppTable2, ...)
   - Each sheet starts with the caption, followed by the table content

Source of truth: `/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration/LIPID.docx`

Notes:
- This script is intentionally conservative: it extracts only images/tables that live
  directly under each Supplementary caption until the next caption.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Optional, Union

from docx import Document
from docx.oxml.table import CT_Tbl  # type: ignore
from docx.oxml.text.paragraph import CT_P  # type: ignore
from docx.table import Table
from docx.text.paragraph import Paragraph


ROOT = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration"
SRC_DOCX = os.path.join(ROOT, "LIPID.docx")
OUT_DIR = os.path.join(ROOT, "Supplementary")
OUT_FIG_DOCX = os.path.join(OUT_DIR, "Supplementary_Figures.docx")
OUT_TAB_XLSX = os.path.join(OUT_DIR, "Supplementary_Tables.xlsx")
RESULTS_DIR = os.path.join(ROOT, "results")


SUPP_FIG_CAP_RE = re.compile(r"^Supplementary\s+Figure\s+(\d+)\.\s*(.*)$", re.IGNORECASE)
SUPP_TAB_CAP_RE = re.compile(r"^Supplementary\s+Table\s+(\d+)\.\s*(.*)$", re.IGNORECASE)
ANY_CAP_RE = re.compile(r"^(Figure|Table|Supplementary\s+(Figure|Table))\s+\d+\.", re.IGNORECASE)


def iter_block_items(doc: Document):
    """Yield Paragraph and Table in document order."""
    body = doc.element.body
    for child in body.iterchildren():
        if isinstance(child, CT_P):
            yield Paragraph(child, doc)
        elif isinstance(child, CT_Tbl):
            yield Table(child, doc)


def paragraph_has_drawing(p: Paragraph) -> bool:
    return bool(p._p.xpath(".//w:drawing"))  # type: ignore[attr-defined]


@dataclass
class SuppFigure:
    number: int
    caption: str
    image_paragraphs: list[Paragraph]
    # Optional external fallback images when the source docx does not embed the image
    fallback_image_paths: list[str]


@dataclass
class SuppTable:
    number: int
    caption: str
    tables: list[Table]


def extract_supplementary_blocks(doc: Document) -> tuple[list[SuppFigure], list[SuppTable]]:
    supp_figs: list[SuppFigure] = []
    supp_tabs: list[SuppTable] = []

    current_fig: Optional[SuppFigure] = None
    current_tab: Optional[SuppTable] = None

    def close_open():
        nonlocal current_fig, current_tab
        if current_fig is not None:
            supp_figs.append(current_fig)
            current_fig = None
        if current_tab is not None:
            supp_tabs.append(current_tab)
            current_tab = None

    for block in iter_block_items(doc):
        if isinstance(block, Paragraph):
            text = (block.text or "").strip()

            mfig = SUPP_FIG_CAP_RE.match(text)
            mtab = SUPP_TAB_CAP_RE.match(text)

            if mfig:
                close_open()
                n = int(mfig.group(1))
                current_fig = SuppFigure(number=n, caption=text, image_paragraphs=[], fallback_image_paths=[])
                continue

            if mtab:
                close_open()
                n = int(mtab.group(1))
                current_tab = SuppTable(number=n, caption=text, tables=[])
                continue

            # Stop capture if we hit ANY caption while inside a block
            if (current_fig is not None or current_tab is not None) and ANY_CAP_RE.match(text):
                close_open()
                # do not continue; this paragraph will be processed by next loop iteration implicitly
                # (we've already seen it). So just start fresh.
                continue

            # Capture images under a Supplementary Figure
            if current_fig is not None and paragraph_has_drawing(block):
                current_fig.image_paragraphs.append(block)

        else:
            # Table block
            if current_tab is not None:
                current_tab.tables.append(block)

    close_open()
    supp_figs.sort(key=lambda x: x.number)
    supp_tabs.sort(key=lambda x: x.number)
    return supp_figs, supp_tabs


def infer_supp_figure_fallback_images(fig: SuppFigure) -> list[str]:
    """
    If a Supplementary Figure caption has no embedded image, infer the expected plot file
    from the repository `plots/` outputs based on caption keywords.
    """
    plots = os.path.join(ROOT, "plots")
    cap = fig.caption.lower()

    # PCA
    if "principal component" in cap or cap.startswith("supplementary figure 1"):
        return [os.path.join(plots, "pca_analysis.png")]

    # KM curves (multiple clustering methods)
    if "kaplan-meier" in cap or "survival curves" in cap:
        if "kmeans" in cap:
            return [os.path.join(plots, "survival_kmeans_subtype.png")]
        if "hierarchical" in cap:
            return [os.path.join(plots, "survival_hierarchical_subtype.png")]
        if "nmf" in cap:
            return [os.path.join(plots, "survival_nmf_subtype.png")]

    # Survival RF regression importance
    if "feature importance" in cap and "survival" in cap:
        return [os.path.join(plots, "survival_feature_importance.png")]

    return []


def write_supplementary_figures(figs: list[SuppFigure]) -> None:
    out = Document()
    out.add_paragraph("Supplementary Figures")

    # Handle duplicate numbering in source (e.g., two "Supplementary Figure 3" captions)
    seen_counts: dict[int, int] = {}

    for fig in figs:
        seen_counts[fig.number] = seen_counts.get(fig.number, 0) + 1
        occurrence = seen_counts[fig.number]
        label_suffix = ""
        if max(1, sum(1 for f in figs if f.number == fig.number)) > 1:
            # 1->A, 2->B, ...
            label_suffix = chr(ord("A") + occurrence - 1)

        out.add_paragraph("")  # spacing
        if label_suffix:
            out.add_paragraph(fig.caption.replace(f"Supplementary Figure {fig.number}.", f"Supplementary Figure {fig.number}{label_suffix}."))
        else:
            out.add_paragraph(fig.caption)

        # If source docx has no embedded image(s), fall back to known plot outputs
        fallback_paths = infer_supp_figure_fallback_images(fig)

        # Re-embed images by copying the underlying relationship target blobs.
        # python-docx doesn't expose direct "copy picture" API; easiest stable path:
        # use the existing runs' picture parts by re-inserting the same image bytes.
        for p in fig.image_paragraphs:
            # Find all blip embed rIds in this paragraph
            blips = p._p.xpath(".//a:blip")  # type: ignore[attr-defined]
            if not blips:
                continue
            for blip in blips:
                rId = blip.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed")
                if not rId:
                    continue
                part = p.part.related_parts[rId]
                img_bytes = part.blob
                # write to temp file (in-memory not supported by add_picture)
                tmp_path = os.path.join(OUT_DIR, f"._tmp_suppfig_{fig.number}_{rId}.png")
                with open(tmp_path, "wb") as f:
                    f.write(img_bytes)
                out.add_picture(tmp_path)
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass

        if not fig.image_paragraphs:
            if fallback_paths:
                for fp in fallback_paths:
                    if os.path.exists(fp):
                        out.add_picture(fp)
                    else:
                        out.add_paragraph(f"[Missing expected image file: {os.path.basename(fp)}]")
            else:
                out.add_paragraph("[Image not detected under caption in source document.]")

    out.save(OUT_FIG_DOCX)


def sanitize_sheet_name(name: str) -> str:
    # Excel sheet name constraints: max 31 chars; no []:*?/\\
    bad = r"[]:*?/\\"
    for ch in bad:
        name = name.replace(ch, "_")
    return name[:31]


def write_supplementary_tables_xlsx(tabs: list[SuppTable]) -> None:
    from openpyxl import Workbook
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    # Remove default sheet
    default = wb.active
    wb.remove(default)

    for st in tabs:
        ws = wb.create_sheet(sanitize_sheet_name(f"SuppTable{st.number}"))
        ws["A1"] = st.caption

        row_cursor = 3
        if not st.tables:
            # Fallback: load from results/*.csv where applicable
            csv_map = {
                2: "lipid_metabolite_correlations.csv",
                3: "differential_features.csv",
                4: "lipid_metabolite_category_correlation.csv",
                5: "cross_validation_results.csv",
                6: "bootstrap_stability_results.csv",
                7: "hub_molecule_analysis.csv",
                8: "top_20_hub_molecules.csv",
                9: "cluster_assignments.csv",
                10: "network_statistics.csv",
                11: "detailed_pathway_enrichment.csv",
                12: "pca_loadings.csv",
                13: "summary_statistics.csv",
                14: "lipid_class_abundance.csv",
                15: "metabolite_category_abundance.csv",
                16: "feature_importance_rf.csv",
            }
            csv_name = csv_map.get(st.number)
            if csv_name:
                csv_path = os.path.join(RESULTS_DIR, csv_name)
                if os.path.exists(csv_path):
                    import csv

                    with open(csv_path, newline="", encoding="utf-8") as f:
                        reader = csv.reader(f)
                        for r_i, row in enumerate(reader):
                            if r_i == 0:
                                # header
                                for c_i, v in enumerate(row, start=1):
                                    ws.cell(row=row_cursor, column=c_i, value=v)
                                row_cursor += 1
                                continue
                            for c_i, v in enumerate(row, start=1):
                                ws.cell(row=row_cursor, column=c_i, value=v)
                            row_cursor += 1
                else:
                    ws["A3"] = f"[Missing expected CSV source: {csv_name}]"
            else:
                ws["A3"] = "[Table not detected under caption in source document.]"
            continue

        for t_idx, table in enumerate(st.tables, start=1):
            if t_idx > 1:
                row_cursor += 2
                ws[f"A{row_cursor}"] = f"(continued {t_idx})"
                row_cursor += 1

            for r in table.rows:
                col_cursor = 1
                for c in r.cells:
                    text = (c.text or "").strip()
                    ws.cell(row=row_cursor, column=col_cursor, value=text)
                    col_cursor += 1
                row_cursor += 1

        # basic autosize
        max_col = ws.max_column
        for col in range(1, max_col + 1):
            letter = get_column_letter(col)
            # heuristic width
            ws.column_dimensions[letter].width = min(60, max(12, ws.column_dimensions[letter].width or 12))

    wb.save(OUT_TAB_XLSX)


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    doc = Document(SRC_DOCX)
    figs, tabs = extract_supplementary_blocks(doc)

    write_supplementary_figures(figs)
    write_supplementary_tables_xlsx(tabs)

    print("Wrote:", OUT_FIG_DOCX)
    print("Wrote:", OUT_TAB_XLSX)
    print("Supplementary Figures:", [f.number for f in figs])
    print("Supplementary Tables:", [t.number for t in tabs])


if __name__ == "__main__":
    main()

