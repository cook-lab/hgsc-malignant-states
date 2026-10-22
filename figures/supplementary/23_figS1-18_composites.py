#!/usr/bin/env python3
"""Full-page composites of every Supplementary Figure S1-S18 (one lettered PNG per figure).

PURPOSE
    Assembles the promoted standalone panels in final_panels/supplementary/figures/ into
    final_panels/supplementary/figures/composites/SuppFig<N>_composite.png (7.5-in page width, 300 dpi).
    Each row is scaled by ONE factor so the row fills the page width (panel fonts stay consistent within a row);
    panels are top-aligned; bold panel letters (A, B, ...) are drawn at each panel's top-left. Figures whose
    single file already carries its letters (S16-S18) or that have a single panel (S5, S7, S8, S11, S13, S14)
    get no extra letters. S15 is one file holding two panels: A at its top-left and B over the pooled violins
    (INNER). S14 is the canonical-marker dot plot only.
INPUTS
    the promoted SuppFig<N>*.png files named in SPEC, in the folder set on the D line
OUTPUTS
    D/composites/SuppFig<N>_composite.png
MANUSCRIPT PANEL(S)
    Figures S1-S18 (lettered pages).
RUNTIME TIER
    fast.
Usage: python3 23_figS1-18_composites.py [N ...]   (default: all figures whose panels exist)
"""
import os, sys
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = None

D = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "final_panels/supplementary/figures")
OUT = f"{D}/composites"; os.makedirs(OUT, exist_ok=True)
DPI, PAGE_W, MARGIN, GAP_X, GAP_Y = 300, 7.5, 0.1, 0.15, 0.18
SRC_DPI = 450
FONT = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
LETTER_PT = 10

# figure -> rows -> [(file stem, letter or None)]
L = lambda *a: list(a)
SPEC = {
    1: [[("SuppFig1a_atlas_qc_total_counts_by_study", "A")], [("SuppFig1b_atlas_qc_genes_by_study", "B")],
        [("SuppFig1c_atlas_qc_doublet_score_by_study", "C")]],
    2: [[("SuppFig2a_atlas_umap_by_study", "A"), ("SuppFig2b_atlas_study_contribution_level1", "B")],
        [("SuppFig2c_atlas_per_study_umaps", "C")]],
    3: [[("SuppFig3a_atlas_umap_treatment_status", "A"), ("SuppFig3b_atlas_umap_anatomic_site", "B"), ("SuppFig3c_atlas_umap_metastatic_site", "C")],
        [("SuppFig3d_atlas_umap_TP53_status", "D"), ("SuppFig3e_atlas_umap_HRD_status", "E"), ("SuppFig3f_atlas_umap_BRCA_status", "F")]],
    4: [[("SuppFig4a_atlas_scanvi_prediction_probability", "A")], [("SuppFig4b_atlas_canonical_markers_dotplot", "B")],
        [("SuppFig4c_atlas_cnv_aneuploid_epithelial_vs_reference", "C")]],
    5: [[("SuppFig5_atlas_nmf_factor_umaps", None)]],
    6: [[("SuppFig6a_atlas_SecA_signature_umaps", "A")], [("SuppFig6b_atlas_SecB_signature_umaps", "B")]],
    7: [[("SuppFig7_atlas_clone_epitype_composition", None)]],
    8: [[("SuppFig8_atlas_functional_heatmap", None)]],
    9: [[("SuppFig9a_atlas_epitype_by_treatment_status", "A"), ("SuppFig9b_atlas_epitype_by_anatomic_site", "B")],
        [("SuppFig9c_atlas_epitype_by_TP53_status", "C"), ("SuppFig9d_atlas_epitype_by_HRD_status", "D"), ("SuppFig9e_atlas_epitype_by_BRCA_status", "E")],
        [("SuppFig9f_atlas_epitype_by_vazquez_garcia_annotation", "F"), ("SuppFig9g_atlas_epitype_by_geistlinger_annotation", "G")]],
    10: [[("SuppFig10a_kpca_ortholog_expression", "A"), ("SuppFig10b_kpca_expressed_ortholog_umaps", "B")],
         [("SuppFig10c_kpca_vs_human_factor_alignment", "C")], [("SuppFig10d_kpca_nmf_program_umaps", "D")]],
    11: [[("SuppFig11_xenium_whole_tissue_celltype", None)], [("SuppFig11_celltype_legend", None)]],
    12: [[("SuppFig12a_xenium_TMA_A_paired_cores_celltype", "A")],
         [("SuppFig12b_xenium_FTE_cores_celltype", "B"), ("SuppFig12_celltype_legend", None)]],
    13: [[("SuppFig13_xenium_TMA_B_cores_celltype", None)], [("SuppFig13_celltype_legend", None)]],
    14: [[("SuppFig14b_xenium_canonical_markers_dotplot", None)]],
    15: [[("SuppFig15_xenium_polarization_distribution", "A")]],
    16: [[("SuppFig16_per_sample_gams_environment", None)]],   # legend is already inside the panel file
    17: [[("SuppFig17_per_sample_gams_epithelial", None)]],
    18: [[("SuppFig18_per_sample_gams_macrophage", None)]],
}
# a legend drawn alone is kept at its native physical size (not stretched to the page width)
# optional per-row max width (inches) to keep tall figures on a letter page
ROW_MAXW = {(2, 1): 5.4, (8, 0): 6.5}
# extra letters inside one file: stem -> [(letter, x fraction, y fraction of the panel)]
INNER = {"SuppFig15_xenium_polarization_distribution": [("B", 0.79, 0.0)]}
NATIVE = {"SuppFig11_celltype_legend", "SuppFig12_celltype_legend", "SuppFig13_celltype_legend", "SuppFig16-18_sample_legend"}


def build(n):
    rows = SPEC[n]
    stems = [s for r in rows for s, _ in r]
    missing = [s for s in stems if not os.path.exists(f"{D}/{s}.png")]
    if missing:
        print(f"SF{n}: skipped (missing {', '.join(missing)})"); return
    px = lambda inch: int(round(inch * DPI))
    usable = px(PAGE_W - 2 * MARGIN)
    placed, y = [], px(MARGIN)
    R = []
    for ri, r in enumerate(rows):
        ims = [(Image.open(f"{D}/{s}.png").convert("RGB"), s, let) for s, let in r]
        nat = [(im.width * DPI / SRC_DPI, im.height * DPI / SRC_DPI) for im, _, _ in ims]   # native size at the page dpi
        R.append((ri, r, ims, nat, px(GAP_X) * (len(ims) - 1)))
    # ONE scale for the whole figure (consistent font size across rows): the largest factor <= 1 at which every
    # non-legend row fits its allowed width
    F = min([1.0] + [((px(ROW_MAXW[(n, ri)]) if (n, ri) in ROW_MAXW else usable) - g) / sum(w for w, _ in nat)
                     for ri, r, ims, nat, g in R if not (len(ims) == 1 and r[0][0] in NATIVE)])
    for ri, r, ims, nat, gaps in R:
        f = min(F, usable / nat[0][0]) if (len(ims) == 1 and r[0][0] in NATIVE) else F   # stand-alone legends shrink with the figure
        roww = int(sum(w for w, _ in nat) * f) + gaps
        x = px(MARGIN) + max(0, (usable - roww) // 2)
        rh = 0
        for (im, s, let), (w, h) in zip(ims, nat):
            W, H = max(1, int(w * f)), max(1, int(h * f))
            placed.append((im.resize((W, H), Image.LANCZOS), x, y, let))
            placed += [(None, x + int(fx * W), y + int(fy * H), l2) for l2, fx, fy in INNER.get(s, [])]
            x += W + px(GAP_X); rh = max(rh, H)
        y += rh + px(GAP_Y)
    page = Image.new("RGB", (px(PAGE_W), y - px(GAP_Y) + px(MARGIN)), "white")
    dr = ImageDraw.Draw(page); font = ImageFont.truetype(FONT, int(LETTER_PT / 72 * DPI))
    for im, x, yy, let in placed:
        if im is not None:
            page.paste(im, (x, yy))
        if let:
            dr.text((x, yy), let, fill="black", font=font)
    fp = f"{OUT}/SuppFig{n}_composite.png"
    page.save(fp, dpi=(DPI, DPI), optimize=True)
    print(f"SF{n}: {page.width / DPI:.1f} x {page.height / DPI:.1f} in -> {os.path.basename(fp)}")


for n in ([int(a) for a in sys.argv[1:]] or sorted(SPEC)):
    build(n)
