#!/usr/bin/env python3
"""
Figure 3l — KRT7 IHC (KPCA.B, day-26 tumour): add a calibrated scale bar to a QuPath export
============================================================================================
PURPOSE
    The 3l image is a QuPath export of a Zeiss 20x/0.8 scan (full resolution 0.2200804 um/px,
    from the CZI metadata / QuPath pixelCalibration). The bar length depends on the export
    downsample, so this script:
      1. checks each export's pixel size against the full scan sizes in the QuPath K7 project
         (downsample d is plausible if export_width * d matches a scan width, or is within it
         for a cropped region) and prints the implied field size for each candidate downsample;
      2. --downsample D: burns the final bar in and writes the panel.

INPUTS
    - one or more QuPath-exported images (png/tif/jpg, or QuPath SVG snapshots wrapping one embedded PNG,
      e.g. data/20260909_KPCA_KRT7 images/krt7_*.svg — 976 x 933 px viewer snapshots)
    - $DATA_ROOT/20260606_murine_secb_ihc/output/qupath_k7/data/*/server.json (scan sizes)
OUTPUTS
    - figures_dir/figure3/Fig3_l_KRT7_IHC_<name>.{png,tif,svg}       (--downsample D); the SVG holds the
      image + a VECTOR bar in the same coordinate frame, so it stays true when resized in Affinity

MANUSCRIPT PANEL(S): Fig 3l (run with --downsample 2 --um 100 on krt7_d26_t8.svg). Bar = black
100 um by default (--um), bottom-right, label in caption (or --label to draw "100 µm"). RUNTIME TIER: fast.

usage:
    python figures/figure3/06_fig3l_krt7_ihc_scale_bar.py --downsample 2 --um 100 img1.png [...]
"""
from __future__ import annotations

import argparse, base64, glob, io, json, os, sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

FULL_UM_PER_PX = 0.2200804394006009
QUPATH_K7 = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "20260606_murine_secb_ihc/output/qupath_k7/data")
CANDIDATES = {2: "10x export (downsample 2)", 10: "downsample 10"}
Image.MAX_IMAGE_PIXELS = None


def load_image(fp):
    """Raster file, or the single embedded raster of a QuPath SVG export."""
    if fp.lower().endswith(".svg"):
        ims = list(ET.parse(fp).getroot().iter("{http://www.w3.org/2000/svg}image"))
        if len(ims) != 1: raise SystemExit(f"{fp}: expected 1 embedded image, found {len(ims)}")
        href = ims[0].get("{http://www.w3.org/1999/xlink}href")
        return Image.open(io.BytesIO(base64.b64decode(href.split(",", 1)[1])))
    return Image.open(fp)


def scan_sizes():
    out = []
    for p in glob.glob(f"{QUPATH_K7}/*/server.json"):
        md = json.load(open(p)).get("metadata", {})
        if md.get("width"): out.append((md["width"], md["height"]))
    return out


def plausibility(w, h, d, scans):
    """Implied full-res size of the export; 'whole' if it matches a scan (±2%), 'crop' if it fits inside one."""
    W, H = w * d, h * d
    for sw, sh in scans:
        if abs(W - sw) / sw < 0.02 and abs(H - sh) / sh < 0.02: return "matches a whole scan"
    if any(W <= sw * 1.02 and H <= sh * 1.02 for sw, sh in scans): return "fits inside a scan (crop)"
    return "LARGER than every scan — impossible"


def draw_bar(im, um, um_per_px, label=None, margin_frac=0.04):
    im = im.convert("RGB").copy(); dr = ImageDraw.Draw(im)
    w, h = im.size; L = um / um_per_px; t = max(2, round(h * 0.012))
    x1 = w - w * margin_frac; x0 = x1 - L; y1 = h - h * margin_frac; y0 = y1 - t
    dr.rectangle([x0, y0, x1, y1], fill=(0, 0, 0))
    if label:
        try: f = ImageFont.truetype("Arial.ttf", max(10, round(h * 0.035)))
        except OSError: f = ImageFont.load_default()
        tw = dr.textlength(label, font=f); dr.text((x0 + (L - tw) / 2, y0 - f.size * 1.25), label, fill=(0, 0, 0), font=f)
    return im, L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("images", nargs="+")
    ap.add_argument("--downsample", type=float)
    ap.add_argument("--um", type=float, default=100)
    ap.add_argument("--label", action="store_true", help='draw "<um> µm" above the bar')
    a = ap.parse_args()
    scans = scan_sizes()
    for fp in a.images:
        im = load_image(fp); w, h = im.size; name = Path(fp).stem
        print(f"\n{fp}: {w} x {h} px")
        for d, desc in CANDIDATES.items():
            print(f"  if {desc:28s}: {FULL_UM_PER_PX*d:.3f} um/px, {a.um:.0f} um = {a.um/(FULL_UM_PER_PX*d):6.1f} px; "
                  f"field {w*FULL_UM_PER_PX*d/1000:.2f} x {h*FULL_UM_PER_PX*d/1000:.2f} mm -> {plausibility(w, h, d, scans)}")
        if a.downsample:
            p, L = draw_bar(im, a.um, FULL_UM_PER_PX * a.downsample, label=(f"{a.um:.0f} µm" if a.label else None))
            for ext in ("png", "tif"):
                op = path("figures_dir", "figure3", f"Fig3_l_KRT7_IHC_{name}.{ext}"); p.save(op, dpi=im.info.get("dpi", (300, 300)))
            buf = io.BytesIO(); im.convert("RGB").save(buf, "PNG"); b64 = base64.b64encode(buf.getvalue()).decode()
            L = a.um / (FULL_UM_PER_PX * a.downsample); t = max(2, round(h * 0.012))
            x1 = w - w * 0.04; y1 = h - h * 0.04
            lab = (f'<text x="{x1 - L / 2:.1f}" y="{y1 - t - 6:.1f}" font-family="Arial" font-size="{round(h * 0.035)}" '
                   f'text-anchor="middle">{a.um:.0f} µm</text>') if a.label else ""
            svg = (f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{w}" height="{h}" viewBox="0 0 {w} {h}">'
                   f'<image width="{w}" height="{h}" xlink:href="data:image/png;base64,{b64}"/>'
                   f'<rect id="scalebar_{a.um:.0f}um" x="{x1 - L:.2f}" y="{y1 - t:.2f}" width="{L:.2f}" height="{t}" fill="#000"/>{lab}</svg>')
            op_svg = path("figures_dir", "figure3", f"Fig3_l_KRT7_IHC_{name}.svg"); open(op_svg, "w").write(svg)
            print(f"  final ({a.downsample:g}x downsample, {a.um:.0f} um = {L:.1f} px) -> {op[:-4]}.{{png,tif}}")


if __name__ == "__main__":
    main()
