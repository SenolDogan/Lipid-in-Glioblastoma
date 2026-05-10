#!/usr/bin/env python3
"""
Combine subtype discovery panels into a single multi-panel PNG.

Mapping used (based on available plot files):
  - Fig 2: clustering_tsne.png
  - Fig 5: feature_importance.png
  - Panel C: top_fdr_subtype_features.png

Output:
  - combined_F2_F5_F6.png

No manuscript files are modified.
"""

from __future__ import annotations

from PIL import Image, ImageDraw, ImageFont


PLOTS_DIR = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration/plots"
FIG2 = f"{PLOTS_DIR}/clustering_tsne.png"
FIG5 = f"{PLOTS_DIR}/feature_importance_shortlabels.png"
FIG6 = f"{PLOTS_DIR}/top_fdr_subtype_features.png"
TIMES_FONT = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"

OUT = f"{PLOTS_DIR}/combined_F2_F5_F6.png"


def load_rgb(path: str) -> Image.Image:
    im = Image.open(path)
    if im.mode != "RGB":
        im = im.convert("RGB")
    return im


def fit_to_width(im: Image.Image, width: int) -> Image.Image:
    w, h = im.size
    if w == width:
        return im
    scale = width / w
    return im.resize((width, int(h * scale)), Image.LANCZOS)


def fit_to_height(im: Image.Image, height: int) -> Image.Image:
    w, h = im.size
    if h == height:
        return im
    scale = height / h
    return im.resize((int(w * scale), height), Image.LANCZOS)


def center_crop_to_width(im: Image.Image, width: int) -> Image.Image:
    w, h = im.size
    if w == width:
        return im
    if w < width:
        # pad if needed (shouldn't happen often)
        canvas = Image.new("RGB", (width, h), "white")
        x0 = (width - w) // 2
        canvas.paste(im, (x0, 0))
        return canvas
    # crop
    x0 = (w - width) // 2
    return im.crop((x0, 0, x0 + width, h))


def pad_to_height(im: Image.Image, height: int, anchor: str = "top") -> Image.Image:
    """Pad image vertically to the requested height without resizing content."""
    w, h = im.size
    if h >= height:
        return im
    canvas = Image.new("RGB", (w, height), "white")
    y = 0 if anchor == "top" else (height - h) // 2
    canvas.paste(im, (0, y))
    return canvas


def add_equal_panel_gutters(im: Image.Image, n_panels: int = 3, gutter: int = 90) -> Image.Image:
    """Add white space between equally sized subpanels without resizing them."""
    w, h = im.size
    panel_w = w // n_panels
    panels = []
    for i in range(n_panels):
        x0 = i * panel_w
        x1 = w if i == n_panels - 1 else (i + 1) * panel_w
        panels.append(im.crop((x0, 0, x1, h)))

    canvas = Image.new("RGB", (w + gutter * (n_panels - 1), h), "white")
    x = 0
    for panel in panels:
        canvas.paste(panel, (x, 0))
        x += panel.size[0] + gutter
    return canvas


def add_panel_label(im: Image.Image, label: str, pad: int = 18) -> Image.Image:
    draw = ImageDraw.Draw(im)
    # Use a default font; if Times not available, fall back.
    try:
        font = ImageFont.truetype(TIMES_FONT, 48)
    except Exception:
        font = ImageFont.load_default()
    # White box behind label
    bbox = draw.textbbox((0, 0), label, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    x0, y0 = pad, pad
    margin = 10
    draw.rectangle([x0 - margin, y0 - margin, x0 + tw + margin, y0 + th + margin], fill="white", outline="black", width=2)
    draw.text((x0, y0), label, fill="black", font=font)
    return im


def main() -> None:
    fig2 = load_rgb(FIG2)
    fig5 = load_rgb(FIG5)
    fig6 = load_rgb(FIG6)

    # Layout: 2-column grid with Fig2 spanning full width on top,
    # then Fig5 (left) and Fig6 (right) on bottom row.
    target_w = 2400  # high-res combined width
    gutter = 40
    margin = 60

    fig2 = fit_to_width(fig2, target_w)
    fig2 = add_equal_panel_gutters(fig2, n_panels=3, gutter=190)

    half_w = (target_w - gutter) // 2
    panel_b_w = int(half_w * 1.38)
    panel_c_w = int(half_w * 0.95)
    fig5 = fit_to_width(fig5, panel_b_w)
    # Give panel B more horizontal room for clear labels/bars. Panel C remains
    # readable but is slightly narrower than before.
    fig6 = fit_to_width(fig6, panel_c_w)
    # Keep B visually larger than C, while aligning panel bottoms.
    fig6 = fit_to_height(fig6, fig5.size[1])

    bottom_w = fig5.size[0] + gutter + fig6.size[0]
    fig2 = fit_to_width(fig2, bottom_w)

    fig2 = add_panel_label(fig2, "A")
    fig5 = add_panel_label(fig5, "B")
    fig6 = add_panel_label(fig6, "C")

    # Compute canvas size
    top_h = fig2.size[1]
    bottom_h = max(fig5.size[1], fig6.size[1])
    canvas_w = margin * 2 + max(target_w, bottom_w)
    canvas_h = margin * 3 + top_h + bottom_h

    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
    x = margin
    y = margin
    canvas.paste(fig2, (x, y))

    y2 = y + top_h + margin
    canvas.paste(fig5, (x, y2))
    canvas.paste(fig6, (x + fig5.size[0] + gutter, y2))

    canvas.save(OUT, dpi=(600, 600))
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()

