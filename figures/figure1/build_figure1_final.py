"""Final Figure 1: panel e text, panel a permeabilisation label, panels b-c replaced by vector persistence diagrams."""
import numpy as np
import pymupdf

import os
from pathlib import Path as _P
_HERE = _P(__file__).resolve().parent
_B = _P(os.environ.get("BUNDLE_DIR", _HERE.parents[1] / "results" / "bundle"))
SRC = os.environ.get("FIG1_SOURCE_PDF", str(_HERE / "Fig1_authors_original.pdf"))   # the authors' original Figure 1
_OUTDIR = _P(os.environ.get("FIG_OUT", _HERE.parent / "output")); _OUTDIR.mkdir(parents=True, exist_ok=True)
OUT = str(_OUTDIR / "Figure1.pdf")
doc = pymupdf.open(SRC); page = doc[0]
rgb = lambda h: tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))

# 1 | panel e: replace text spans (text only; boxes and lines kept)
REPL = {"4,259 cells": "2,782 spots", "10 PCs": "50 PCs", "20×20 px": "6,498 px", "4,400-dim vector": "6,588-dim vector",
        "topology + HVG": "50 PCs + topology", "64-dim →n classes": "256 → 10 classes", "Adam · lr=0.001": "Adam · lr=0.005",
        "64-dim embed.": "256-dim embed.", "11 types · ~90% acc.": "10 types · 89% acc."}
spans = [sp for blk in page.get_text("dict")["blocks"] for line in blk.get("lines", []) for sp in line["spans"]
         if sp["text"].strip() in REPL]
assert len(spans) == len(REPL)
for sp in spans:
    r = pymupdf.Rect(sp["bbox"]); page.add_redact_annot(pymupdf.Rect(r.x0, r.y0 + 0.5, r.x1, r.y1 - 0.5))
page.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE, graphics=pymupdf.PDF_REDACT_LINE_ART_NONE)
for sp in spans:
    new = REPL[sp["text"].strip()]; r = pymupdf.Rect(sp["bbox"])
    if "→" in new:          # Helvetica has no arrow glyph: DejaVu Sans, measured with the font itself
        page.insert_font(fontname="dvs", fontfile="/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
        w = pymupdf.Font(fontfile="/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf").text_length(new, fontsize=sp["size"] * 0.93)
        page.insert_text(((r.x0 + r.x1) / 2 - w / 2, sp["origin"][1]), new, fontname="dvs", fontsize=sp["size"] * 0.93, color=rgb(f"{sp['color']:06x}"))
        continue
    w = pymupdf.get_text_length(new, fontname="helv", fontsize=sp["size"])
    page.insert_text(((r.x0 + r.x1) / 2 - w / 2, sp["origin"][1]), new, fontname="helv", fontsize=sp["size"], color=rgb(f"{sp['color']:06x}"))

# 2 | panel a: "30min 37°C" lives inside an embedded image -> white patch + vector text
z = 600 / 72
pix = page.get_pixmap(matrix=pymupdf.Matrix(z, z), clip=pymupdf.Rect(194, 324.8, 228, 331.8))
a = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)[..., :3].astype(int)
grey = (np.abs(a[..., 0] - a[..., 1]) < 25) & (np.abs(a[..., 1] - a[..., 2]) < 25) & (a.mean(-1) < 140)
ys, xs = np.where(grey)
perm = pymupdf.Rect(194 + xs.min() / z, 324.8 + ys.min() / z, 194 + (xs.max() + 1) / z, 324.8 + (ys.max() + 1) / z)
page.draw_rect(pymupdf.Rect(perm.x0 - 0.6, perm.y0 - 0.5, perm.x1 + 0.6, perm.y1 + 0.5), color=None, fill=(1, 1, 1), overlay=True)
size = (perm.y1 - perm.y0) / 0.72; txt = "40 min, 37 °C"
w = pymupdf.get_text_length(txt, fontname="helv", fontsize=size)
page.insert_text(((perm.x0 + perm.x1) / 2 - w / 2, perm.y1), txt, fontname="helv", fontsize=size, color=rgb("2b2b2b"))

# 3 | panels b-c: remove the old vector content, then place the new vector panels in the same footprint
RECT = {"TYLCV": pymupdf.Rect(428.0, 24.5, 570.0, 206.5), "mock": pymupdf.Rect(575.0, 24.5, 706.0, 206.5)}
page.add_redact_annot(pymupdf.Rect(428.0, 24.5, 706.0, 206.5), fill=False)      # both panels and the gap between them
page.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE, graphics=pymupdf.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                      text=pymupdf.PDF_REDACT_TEXT_REMOVE)
left = [d for d in page.get_drawings() if pymupdf.Rect(428.0, 24.5, 706.0, 206.5).contains(d["rect"])]
for grp, r in RECT.items():
    page.show_pdf_page(r, pymupdf.open(str(_HERE / "panels" / f"panel_{grp}.pdf")), 0)
# 4 | panel d title: the persistence images encode H0 and H1 together
page.draw_rect(pymupdf.Rect(726.5, 24.6, 786.6, 39.4), color=None, fill=(1, 1, 1), overlay=True)
page.insert_font(fontname="dvs", fontfile="/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
for i, line in enumerate(["Persistence image", "(H₀ and H₁)"]):
    w = pymupdf.get_text_length(line, fontname="helv", fontsize=6.3) * 1.04
    page.insert_text((756.55 - w / 2, 30.6 + i * 7.4), line, fontname="dvs", fontsize=6.3, color=(0, 0, 0))

# 5 | panels f-g and their legend: raster UMAPs of the previous model replaced by vector UMAPs of the current one
FG = {"f": pymupdf.Rect(450.5, 469.4, 578.1, 588.2), "g": pymupdf.Rect(611.3, 472.5, 734.2, 588.2),
      "legend": pymupdf.Rect(746.2, 456.7, 825.0, 588.2)}
for r in FG.values():
    page.add_redact_annot(r, fill=False)
page.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_REMOVE, graphics=pymupdf.PDF_REDACT_LINE_ART_NONE,
                      text=pymupdf.PDF_REDACT_TEXT_NONE)
for k, r in FG.items():
    page.show_pdf_page(r, pymupdf.open(str(_HERE / "panels" / f"panel_{k}.pdf")), 0)

doc.save(OUT, garbage=3, deflate=True)
out = pymupdf.open(OUT)[0]
out.get_pixmap(matrix=pymupdf.Matrix(300 / 72, 300 / 72)).save(OUT.replace(".pdf", ".png"))
print("written", OUT, "| old drawings left inside b-c footprints:", len(left), "| permeabilisation font", round(size, 2))
