#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/scripts/00_setup.py, the script that ran. Kept: the paths,
# gene lists, parameters and ortholog helpers that the kept steps import. Removed: the plot style, palettes and
# figure helper, unused paths, gene sets and helpers, and the configuration check. The untrimmed script, which the
# docstring below describes, is in the lab project folder.

"""
Step 00 — Setup: paths, style, gene sets, shared helpers
=========================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

Imported by every downstream script:

    import sys; sys.path.insert(0, os.path.dirname(__file__))
    from importlib import import_module
    S = import_module("00_setup")

Running it directly prints a configuration check.
"""

import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

BASE      = "/Volumes/CookLab/Sarah/202608_kpca_secb_mouse"
OUTPUT    = os.path.join(BASE, "output")
EXTERNAL  = os.path.join(BASE, "data", "external")

# --- Mouse source data (KPCA-B_B2, one sample of 3-plex FlexPool2) -----------
KPCA_DIR   = "/Volumes/CookLab/Sarah/2025_kpca/KPCA-B"
KPCA_CR    = os.path.join(KPCA_DIR, "sample_filtered_feature_bc_matrix.h5")
KPCA_CB    = os.path.join(KPCA_DIR, "cellbender",
                          "KPCA-B_background_removed_filtered.h5")
KPCA_ANNOT = os.path.join(KPCA_DIR, "cell_types.csv")

# CellRanger called exactly this many cells for KPCA-B_B2. CellBender ran on the
# per-sample RAW matrix and over-called (45,585 barcodes, its own prior was
# 20,336, convergence_indicator 0.505). Its barcode set is a strict SUPERSET of
# CellRanger's, so we intersect down to these and never trust its cell calls.
N_CELLS_CELLRANGER = 28973

# --- Human atlas (what we transfer from) ------------------------------------
ATLAS      = "/Volumes/CookLab/Sarah/202601_final_atlas"
H_GENECLASS = os.path.join(ATLAS, "output", "11e_nmf_characterization",
                           "11e_gene_classification.csv")
H_LOADINGS  = os.path.join(ATLAS, "output", "11d_epithelial_nmf",
                           "11d_nmf_loadings.csv")

# 11d factor identities, established in 11d_factor_scores.csv:
#   Factor_2 -> SecB (7/8 canonical SecB genes in top-100)
#   Factor_3 -> SecA (8/9 canonical SecA genes in top-100)
H_SECB_FACTOR = "Factor_2"
H_SECA_FACTOR = "Factor_3"

# --- External resources -----------------------------------------------------
MGI_URL = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
MGI_RPT = os.path.join(EXTERNAL, "HOM_MouseHumanSequence.rpt")
ORTHO_MAP = os.path.join(OUTPUT, "02_orthologs", "ortholog_map.csv")

# ============================================================================
# GENE SIGNATURES  (human symbols; ortholog-mapped downstream)
# ============================================================================

SECA_GENES = ["BCAM", "MECOM", "FBXO21", "LGR5", "LPAR3", "PBX1", "SOX17",
              "WT1", "RCN2"]
SECB_GENES = ["KRT17", "KRT19", "KRT7", "LCN2", "MMP7", "PRSS22", "SLPI",
              "TACSTD2"]

# ============================================================================
# MOUSE COMPARTMENT MARKERS  (for de novo annotation in 01)
# ============================================================================

MOUSE_COMPARTMENT_MARKERS = {
    "Epithelial":  ["Epcam", "Krt8", "Krt18", "Krt7", "Cdh1", "Krt19"],
    "Fibroblast":  ["Col1a1", "Col1a2", "Dcn", "Pdgfrb", "Lum", "Fn1"],
    "Endothelial": ["Pecam1", "Cdh5", "Kdr", "Egfl7"],
    "Myeloid":     ["Lyz2", "Adgre1", "Csf1r", "Itgam", "C1qa"],
    "T/NK":        ["Cd3e", "Cd3d", "Trbc2", "Nkg7", "Klrb1c"],
    "B":           ["Cd79a", "Cd79b", "Ms4a1"],
    "Mast":        ["Cma1", "Cpa3", "Ms4a2"],
    "Immune_pan":  ["Ptprc"],
}

# ============================================================================
# PARAMETERS
# ============================================================================

SEED    = 42
N_HVG   = 3000
NMF_KS  = list(range(4, 13))   # k = 4..12
NMF_SEEDS = 20                 # random inits per k
N_PERM  = 10_000               # gene-label permutation null in 05

# ============================================================================
# HELPERS
# ============================================================================


def human_to_mouse_naive(sym):
    """Case-conversion fallback: KRT17 -> Krt17. Used only where MGI has no
    entry; 02_orthologs.py reports how often this fires."""
    if not isinstance(sym, str) or not sym:
        return None
    return sym[0].upper() + sym[1:].lower() if sym[0].isalpha() else sym


def load_ortholog_map():
    """human_symbol -> mouse_symbol (1:1 preferred). Requires 02_orthologs.py."""
    if not os.path.exists(ORTHO_MAP):
        raise FileNotFoundError(
            f"{ORTHO_MAP} missing — run scripts/02_orthologs.py first")
    m = pd.read_csv(ORTHO_MAP)
    return m


def map_genes(human_genes, ortho, restrict_to=None):
    """Map a list of human symbols to mouse symbols.

    Returns (mapped_dict, report_df). `restrict_to` is the set of symbols
    actually present in the mouse matrix, so the report can distinguish
    'no ortholog' from 'ortholog exists but not in this matrix'.
    """
    lut = dict(zip(ortho["human_symbol"], ortho["mouse_symbol"]))
    rows, mapped = [], {}
    for g in human_genes:
        mm = lut.get(g)
        src = "mgi"
        if mm is None or (isinstance(mm, float) and np.isnan(mm)):
            mm, src = human_to_mouse_naive(g), "case_conversion"
        present = restrict_to is None or (mm in restrict_to)
        if mm is not None and present:
            mapped[g] = mm
        rows.append({"human_symbol": g, "mouse_symbol": mm, "source": src,
                     "in_matrix": bool(present),
                     "status": ("mapped" if (mm and present) else
                                ("no_ortholog" if not mm else "off_panel"))})
    return mapped, pd.DataFrame(rows)
