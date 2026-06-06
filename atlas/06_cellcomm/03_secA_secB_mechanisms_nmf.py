#!/usr/bin/env python3
"""
SecA vs SecB communication mechanisms (NMF labels)
==================================================
HGSC malignant-states atlas backend.

Deeper analysis of the global LIANA+ results (NMF epithelial labels) from
01_cellcomm_nmf.py: (1) pathway-level aggregation of L-R pairs, (2) level2 TME
subtype partners, and (3) cross-reference of SecA-vs-SecB DEGs with the LIANA
L-R pairs to flag communication-active differentially expressed genes.

Migration note: this consolidates the original 17d wrapper, which exec-patched
the (Leiden-based, non-canonical) 16c script at runtime. The patches -- NMF
label schema, 17b LIANA input, 17d output prefixes -- are baked in here as a
standalone script. "Transitioning" is renamed to "Intermediate" per convention.

INPUTS:
  - output_root/06_cellcomm/tables/17b_liana_global.csv  (from 01_cellcomm_nmf.py)
  - <data_root>/2026_final_atlas/output/fig_secretory_polarization/data/
      panel_i_deg_results.csv  (SecA-vs-SecB Wilcoxon DEGs)

OUTPUTS (output_root/06_cellcomm/17d_secA_secB_mechanisms_nmf/):
  - tables/17d_*.csv, figs/17d_*.svg/pdf, 17d_secA_secB_mechanisms_nmf.html

MANUSCRIPT PANELS: supporting Fig 5F mechanism analysis / Supp Data 7.

RUNTIME TIER: fast (post-hoc analysis of LIANA + DEG tables).

SEEDING: deterministic (no stochastic step).

Usage:
    python 03_secA_secB_mechanisms_nmf.py
"""

import base64
import io
import os
import warnings

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

warnings.filterwarnings("ignore")

import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402

# ── Cook Lab style v1.2 ─────────────────────────────────────
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6,
    "figure.dpi":      150,
    "savefig.dpi":     450,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
    "savefig.bbox":    "tight",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

# ── Paths ────────────────────────────────────────────────────
LIANA_CSV   = path("output_root", "06_cellcomm", "tables", "17b_liana_global.csv")
DEG_CSV     = path("data_root", "2026_final_atlas", "output",
                   "fig_secretory_polarization", "data", "panel_i_deg_results.csv")
OUT_DIR     = path("output_root", "06_cellcomm", "17d_secA_secB_mechanisms_nmf")
FIG_DIR     = os.path.join(OUT_DIR, "figs")
TABLE_DIR   = os.path.join(OUT_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# ── v2 epithelial schema (from 16b) ─────────────────────────
V2_MAP = {
    "SecA epithelium":          "SecA",
    "Intermediate epithelium":  "Intermediate",
    "SecB epithelium":          "SecB",
    "Ciliated epithelial cell": "Ciliated",
}

V2_PALETTE = {
    "SecA":            "#E6A141",
    "SecB":            "#B8741A",
    "Stress-response": "#D9C5A2",
    "Intermediate":    "#C49A5E",
    "Ciliated":        "#E05A2C",
}

# ── Level-1 compartment mapping (from 16b) ──────────────────
LEVEL1_MAP = {
    # Epithelial (NMF labels)
    "SecA epithelium":          "Epithelial",
    "Intermediate epithelium":  "Epithelial",
    "SecB epithelium":          "Epithelial",
    "Ciliated epithelial cell": "Epithelial",
    # Mesothelial
    "Mesothelial cell":                          "Mesothelial",
    "Hypoxic mesothelial cell":                  "Mesothelial",
    # Fibroblast
    "Activated fibroblast":                      "Fibroblast",
    "Myo-fibroblastic cancer-associated fibroblast": "Fibroblast",
    "PI16+ universal fibroblast":                "Fibroblast",
    "Cycling fibroblast":                        "Fibroblast",
    "Inflammatory fibroblast":                   "Fibroblast",
    "Hypoxic inflammatory cancer-associated fibroblast": "Fibroblast",
    "Ovarian stromal cell":                      "Fibroblast",
    "Ovarian steroidogenic cell":                "Fibroblast",
    "Schwann cell":                              "Fibroblast",
    # Smooth muscle
    "Contractile smooth muscle cell":            "Smooth muscle",
    "Stress-response smooth muscle cell":        "Smooth muscle",
    "Inflammatory fibroblast-like smooth muscle cell": "Smooth muscle",
    # Pericyte
    "Pericyte":                                  "Pericyte",
    # Endothelial
    "Angiogenic endothelial cell":               "Endothelial",
    "Venous endothelial cell":                   "Endothelial",
    "Arterial endothelial cell":                 "Endothelial",
    "Lymphatic endothelial cell":                "Endothelial",
    "Cycling endothelial cell":                  "Endothelial",
    # T/NK
    "CD8 effector/exhausted T cell":             "T/NK cell",
    "CD8 effector T cell":                       "T/NK cell",
    "CD8 tissue-resident memory T cell":         "T/NK cell",
    "CD4 naive T cell":                          "T/NK cell",
    "CD4 Regulatory T cell":                     "T/NK cell",
    "CD4 regulatory T cell":                     "T/NK cell",
    "Cycling T/NK cell":                         "T/NK cell",
    "Quiescent T cell":                          "T/NK cell",
    "MAIT cell":                                 "T/NK cell",
    "NK cell":                                   "T/NK cell",
    "CD56bright NK cell":                        "T/NK cell",
    "CD56dim NK cell":                           "T/NK cell",
    "Metallothionein-high stress-response T cell": "T/NK cell",
    "Metallothionein-stress T cell":             "T/NK cell",
    "gdT cell":                                  "T/NK cell",
    # B cell
    "IFN-activated B cell":                      "B cell",
    "Activated B cell":                          "B cell",
    "Cycling B cell":                            "B cell",
    # Plasma
    "Plasma cell":                               "Plasma cell",
    "Cycling plasma cell":                       "Plasma cell",
    # Macrophage
    "C1Q tissue-resident macrophage":            "Macrophage",
    "Cycling C1Q+ tissue-resident macrophage":   "Macrophage",
    "Monocyte-derived macrophage":               "Macrophage",
    "Inflammatory macrophage":                   "Macrophage",
    "Hypoxic macrophage":                        "Macrophage",
    "Cycling macrophage":                        "Macrophage",
    "Monocyte":                                  "Macrophage",
    "Classical monocyte":                        "Macrophage",
    # DC
    "Type 1 DC":                                 "DC",
    "Type 2 DC":                                 "DC",
    "Conventional dendritic cell type 1":        "DC",
    "Conventional dendritic cell type 2":        "DC",
    "Mature DC":                                 "DC",
    "Mature dendritic cell":                     "DC",
    "Plasmacytoid DC":                           "DC",
    "Plasmacytoid dendritic cell":               "DC",
    # Neutrophil / Mast
    "Neutrophil":                                "Neutrophil",
    "Mast cell":                                 "Mast cell",
    "Cycling mast cell":                         "Mast cell",
    # Other
    "Hematopoietic stem cell":                   "Other",
}

COMPARTMENT_ORDER = [
    "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte", "Endothelial",
    "T/NK cell", "B cell", "Plasma cell", "Macrophage", "DC",
    "Neutrophil", "Mast cell",
]

SIG_THRESH = 0.05  # magnitude_rank <= 0.05

# ── Curated L-R to pathway mapping ──────────────────────────
# Keys: gene symbols that appear as ligand_complex or receptor_complex
# (or as components in underscore-delimited complexes).
# Each gene maps to one primary pathway. For genes in multiple pathways,
# the most canonical assignment is used.

GENE_TO_PATHWAY = {
    # --- WNT ---
    "WNT5A": "WNT", "WNT5B": "WNT", "WNT3A": "WNT", "WNT2": "WNT",
    "WNT4": "WNT", "WNT7A": "WNT", "WNT7B": "WNT", "WNT11": "WNT",
    "DKK3": "WNT", "SFRP1": "WNT", "SFRP2": "WNT",
    "FZD1": "WNT", "FZD2": "WNT", "FZD3": "WNT", "FZD4": "WNT",
    "FZD5": "WNT", "FZD6": "WNT", "FZD7": "WNT", "FZD8": "WNT",
    "LRP5": "WNT", "LRP6": "WNT", "RYK": "WNT", "PTK7": "WNT",
    "KREMEN1": "WNT",

    # --- NOTCH ---
    "JAG1": "NOTCH", "JAG2": "NOTCH", "DLL1": "NOTCH", "DLL3": "NOTCH",
    "DLL4": "NOTCH", "MFNG": "NOTCH",
    "NOTCH1": "NOTCH", "NOTCH2": "NOTCH", "NOTCH3": "NOTCH",
    "NOTCH4": "NOTCH", "MAML2": "NOTCH", "NCSTN": "NOTCH",

    # --- TGFb ---
    "TGFB1": "TGFb", "TGFB2": "TGFb", "TGFB3": "TGFb",
    "GDF11": "TGFb", "INHBA": "TGFb",
    "TGFBR1": "TGFb", "TGFBR2": "TGFb", "TGFBR3": "TGFb",
    "ACVR1": "TGFb", "ACVR1B": "TGFb", "ACVR2A": "TGFb",
    "ACVRL1": "TGFb", "ENG": "TGFb",
    "LTBP1": "TGFb", "LTBP3": "TGFb",

    # --- BMP ---
    "BMP1": "BMP", "BMP2": "BMP", "BMP4": "BMP", "BMP7": "BMP",
    "RGMB": "BMP",
    "BMPR1A": "BMP", "BMPR2": "BMP",

    # --- Hedgehog ---
    "SHH": "Hedgehog", "IHH": "Hedgehog", "DHH": "Hedgehog",
    "PTCH1": "Hedgehog", "PTCH2": "Hedgehog", "SMO": "Hedgehog",

    # --- JAK-STAT / Interleukin ---
    "IL1B": "JAK-STAT/Interleukin", "IL6": "JAK-STAT/Interleukin",
    "IL15RA": "JAK-STAT/Interleukin", "IL16": "JAK-STAT/Interleukin",
    "IL18": "JAK-STAT/Interleukin", "OSM": "JAK-STAT/Interleukin",
    "CLCF1": "JAK-STAT/Interleukin", "LIF": "JAK-STAT/Interleukin",
    "CSF1": "JAK-STAT/Interleukin",
    "IL1R1": "JAK-STAT/Interleukin", "IL1RAP": "JAK-STAT/Interleukin",
    "IL6R": "JAK-STAT/Interleukin", "IL6ST": "JAK-STAT/Interleukin",
    "IL2RG": "JAK-STAT/Interleukin", "IL17RC": "JAK-STAT/Interleukin",
    "IL18BP": "JAK-STAT/Interleukin", "LIFR": "JAK-STAT/Interleukin",
    "OSMR": "JAK-STAT/Interleukin", "CSF3R": "JAK-STAT/Interleukin",
    "LEPR": "JAK-STAT/Interleukin",

    # --- TNFa ---
    "TNF": "TNFa", "TNFSF10": "TNFa", "TNFSF12": "TNFa",
    "TNFSF13B": "TNFa", "TNFSF14": "TNFa", "TNFSF9": "TNFa",
    "LTB": "TNFa", "FADD": "TNFa",
    "TNFRSF1A": "TNFa", "TNFRSF1B": "TNFa", "TNFRSF10B": "TNFa",
    "TNFRSF10D": "TNFa", "TNFRSF12A": "TNFa", "TNFRSF13C": "TNFa",
    "TNFRSF14": "TNFa", "TNFRSF21": "TNFa", "TNFRSF25": "TNFa",
    "TNFRSF9": "TNFa", "LTBR": "TNFa", "FAS": "TNFa",
    "TRADD": "TNFa", "TRAF2": "TNFa", "RIPK1": "TNFa",

    # --- EGFR ---
    "AREG": "EGFR", "HBEGF": "EGFR", "EGF": "EGFR", "EREG": "EGFR",
    "EGFR": "EGFR", "ERBB2": "EGFR", "ERBB3": "EGFR",

    # --- PDGF ---
    "PDGFA": "PDGF", "PDGFB": "PDGF", "PDGFC": "PDGF", "PDGFD": "PDGF",
    "PDGFRA": "PDGF", "PDGFRB": "PDGF",

    # --- VEGF ---
    "VEGFA": "VEGF", "VEGFB": "VEGF", "VEGFC": "VEGF",
    "PGF": "VEGF", "PIGF": "VEGF",
    "KDR": "VEGF", "FLT1": "VEGF",
    "NRP1": "VEGF", "NRP2": "VEGF",

    # --- FGF ---
    "FGF7": "FGF", "FGF23": "FGF", "FGF1": "FGF", "FGF2": "FGF",
    "FGFR1": "FGF", "FGFR2": "FGF", "FGFR3": "FGF",

    # --- HGF/MET ---
    "HGF": "HGF/MET", "MET": "HGF/MET",

    # --- Semaphorin ---
    "SEMA3C": "Semaphorin", "SEMA4A": "Semaphorin", "SEMA4B": "Semaphorin",
    "SEMA4C": "Semaphorin", "SEMA4D": "Semaphorin", "SEMA5A": "Semaphorin",
    "SEMA7A": "Semaphorin",
    "PLXNA1": "Semaphorin", "PLXNA2": "Semaphorin", "PLXNA3": "Semaphorin",
    "PLXNB1": "Semaphorin", "PLXNB2": "Semaphorin", "PLXNC1": "Semaphorin",
    "PLXND1": "Semaphorin",

    # --- Ephrin ---
    "EFNA1": "Ephrin", "EFNA4": "Ephrin", "EFNA5": "Ephrin",
    "EFNB1": "Ephrin", "EFNB2": "Ephrin",
    "EPHA2": "Ephrin", "EPHB2": "Ephrin", "EPHB4": "Ephrin",
    "EPHB6": "Ephrin",

    # --- Chemokine ---
    "CCL2": "Chemokine", "CCL3": "Chemokine", "CCL4": "Chemokine",
    "CCL5": "Chemokine", "CCL20": "Chemokine",
    "CXCL2": "Chemokine", "CXCL10": "Chemokine", "CXCL12": "Chemokine",
    "CXCL14": "Chemokine", "CKLF": "Chemokine",
    "CCR1": "Chemokine", "CCR5": "Chemokine", "CCRL2": "Chemokine",
    "CXCR3": "Chemokine", "CXCR4": "Chemokine", "ACKR3": "Chemokine",

    # --- Complement ---
    "C1QA": "Complement", "C1QB": "Complement", "C3": "Complement",
    "CFH": "Complement", "SERPING1": "Complement",
    "C1QBP": "Complement", "C3AR1": "Complement", "C5AR1": "Complement",
    "CD46": "Complement",

    # --- MHC / Antigen presentation ---
    "B2M": "MHC/Antigen presentation", "HLA-A": "MHC/Antigen presentation",
    "HLA-B": "MHC/Antigen presentation", "HLA-C": "MHC/Antigen presentation",
    "HLA-E": "MHC/Antigen presentation", "HLA-F": "MHC/Antigen presentation",
    "CD74": "MHC/Antigen presentation", "CALR": "MHC/Antigen presentation",
    "CD4": "MHC/Antigen presentation", "CD2": "MHC/Antigen presentation",
    "CD6": "MHC/Antigen presentation",

    # --- Integrin ---
    "ITGA1": "Integrin", "ITGA2": "Integrin", "ITGA3": "Integrin",
    "ITGA4": "Integrin", "ITGA5": "Integrin", "ITGA6": "Integrin",
    "ITGA11": "Integrin", "ITGAE": "Integrin",
    "ITGAL": "Integrin", "ITGAM": "Integrin", "ITGAV": "Integrin",
    "ITGAX": "Integrin", "ITGB1": "Integrin", "ITGB2": "Integrin",
    "ITGB3BP": "Integrin", "ITGB4": "Integrin", "ITGB5": "Integrin",
    "ITGB7": "Integrin", "ITGB8": "Integrin",
    # Integrin ligands
    "FN1": "Integrin", "ICAM1": "Integrin", "ICAM2": "Integrin",
    "ICAM3": "Integrin", "VCAM1": "Integrin",
    "COL1A1": "Collagen/ECM", "COL1A2": "Collagen/ECM",
    "COL3A1": "Collagen/ECM", "COL4A1": "Collagen/ECM",
    "COL4A2": "Collagen/ECM", "COL4A5": "Collagen/ECM",
    "COL5A1": "Collagen/ECM", "COL5A2": "Collagen/ECM",
    "COL5A3": "Collagen/ECM", "COL6A1": "Collagen/ECM",
    "COL6A2": "Collagen/ECM", "COL6A3": "Collagen/ECM",
    "COL8A1": "Collagen/ECM", "COL12A1": "Collagen/ECM",
    "COL14A1": "Collagen/ECM", "COL15A1": "Collagen/ECM",
    "COL16A1": "Collagen/ECM", "COL18A1": "Collagen/ECM",
    "COL27A1": "Collagen/ECM",
    "DDR1": "Collagen/ECM", "DDR2": "Collagen/ECM",
    "LAMA4": "Collagen/ECM", "LAMA5": "Collagen/ECM",
    "LAMB1": "Collagen/ECM", "LAMB2": "Collagen/ECM",
    "LAMC1": "Collagen/ECM",

    # --- Galectin ---
    "LGALS1": "Galectin", "LGALS3": "Galectin", "LGALS3BP": "Galectin",
    "LGALS8": "Galectin", "LGALS9": "Galectin",

    # --- GAS6/AXL/MERTK (TAM receptor) ---
    "GAS6": "GAS6/TAM", "PROS1": "GAS6/TAM",
    "AXL": "GAS6/TAM", "MERTK": "GAS6/TAM",

    # --- Midkine/Pleiotrophin ---
    "MDK": "Midkine/PTN", "PTN": "Midkine/PTN",

    # --- SPP1/Osteopontin ---
    "SPP1": "SPP1/Osteopontin",

    # --- Thrombospondin ---
    "THBS1": "Thrombospondin", "THBS2": "Thrombospondin",
    "THBS3": "Thrombospondin",
    "CD47": "Thrombospondin", "SIRPA": "Thrombospondin",

    # --- MIF ---
    "MIF": "MIF",

    # --- SLIT/ROBO ---
    "SLIT2": "SLIT/ROBO", "SLIT3": "SLIT/ROBO",
    "ROBO1": "SLIT/ROBO", "ROBO3": "SLIT/ROBO",

    # --- ANGPT ---
    "ANGPT2": "Angiopoietin", "ANGPTL2": "Angiopoietin",
    "ANGPTL4": "Angiopoietin",
    "TIE1": "Angiopoietin",

    # --- Protease/TIMP ---
    "ADAM10": "Protease/TIMP", "ADAM12": "Protease/TIMP",
    "ADAM15": "Protease/TIMP", "ADAM17": "Protease/TIMP",
    "ADAM28": "Protease/TIMP", "ADAM9": "Protease/TIMP",
    "MMP2": "Protease/TIMP", "MMP7": "Protease/TIMP",
    "PLAU": "Protease/TIMP", "PLAT": "Protease/TIMP",
    "PLAUR": "Protease/TIMP",
    "TIMP1": "Protease/TIMP", "TIMP2": "Protease/TIMP",
    "TIMP3": "Protease/TIMP",
    "SERPINA1": "Protease/TIMP", "SERPINE1": "Protease/TIMP",
    "SERPINE2": "Protease/TIMP", "SERPINF1": "Protease/TIMP",

    # --- Immune checkpoint ---
    "TIGIT": "Immune checkpoint", "CD96": "Immune checkpoint",
    "PVR": "Immune checkpoint", "LAG3": "Immune checkpoint",
    "HAVCR2": "Immune checkpoint", "CD40": "Immune checkpoint",
    "TNFSF13B": "Immune checkpoint",
    "LILRB2": "Immune checkpoint", "LILRB3": "Immune checkpoint",
    "LILRB4": "Immune checkpoint",

    # --- Tetraspanin ---
    "CD9": "Tetraspanin", "CD63": "Tetraspanin", "CD81": "Tetraspanin",
    "CD82": "Tetraspanin", "CD151": "Tetraspanin",
    "TSPAN1": "Tetraspanin", "TSPAN5": "Tetraspanin",
    "TSPAN12": "Tetraspanin", "TSPAN14": "Tetraspanin",
    "TSPAN15": "Tetraspanin", "TSPAN17": "Tetraspanin",

    # --- IGF ---
    "IGF1": "IGF", "IGFBP4": "IGF",
    "IGF1R": "IGF", "IGF2R": "IGF", "INSR": "IGF",

    # --- Endothelin ---
    "EDN1": "Endothelin", "EDN3": "Endothelin",
    "EDNRA": "Endothelin", "EDNRB": "Endothelin",

    # --- S100 / Alarmin ---
    "S100A1": "S100/Alarmin", "S100A4": "S100/Alarmin",
    "S100A8": "S100/Alarmin", "S100A9": "S100/Alarmin",

    # --- Selectin ---
    "SELL": "Selectin", "SELPLG": "Selectin",

    # --- Nectin/CD ---
    "NECTIN1": "Nectin", "NECTIN2": "Nectin", "NECTIN3": "Nectin",

    # --- APP ---
    "APP": "APP",

    # --- CD44 ---
    "CD44": "CD44/Hyaluronan", "HAS2": "CD44/Hyaluronan",

    # --- BSG/CD147 ---
    "BSG": "BSG/CD147",

    # --- ANXA ---
    "ANXA1": "Annexin", "ANXA2": "Annexin",

    # --- GRN ---
    "GRN": "Granulin", "SORT1": "Granulin",

    # --- POSTN ---
    "POSTN": "Periostin",

    # --- NAMPT/PBEF ---
    "NAMPT": "NAMPT/PBEF",

    # --- ADM ---
    "ADM": "Adrenomedullin", "CALCRL": "Adrenomedullin",
    "RAMP1": "Adrenomedullin", "RAMP2": "Adrenomedullin",
}

PATHWAY_ORDER = [
    "WNT", "NOTCH", "TGFb", "BMP", "Hedgehog",
    "JAK-STAT/Interleukin", "TNFa", "EGFR", "PDGF", "VEGF", "FGF",
    "HGF/MET", "IGF", "Semaphorin", "Ephrin", "SLIT/ROBO",
    "Chemokine", "Complement", "MHC/Antigen presentation",
    "Immune checkpoint", "Integrin", "Collagen/ECM",
    "Galectin", "Thrombospondin", "Tetraspanin",
    "GAS6/TAM", "Midkine/PTN", "SPP1/Osteopontin", "MIF",
    "Protease/TIMP", "S100/Alarmin", "Angiopoietin",
    "CD44/Hyaluronan", "Annexin", "Adrenomedullin",
]


# ============================================================================
# UTILITY
# ============================================================================

def get_genes_from_complex(complex_name):
    """Split underscore-delimited complex into individual gene symbols."""
    if pd.isna(complex_name):
        return []
    return str(complex_name).split("_")


def assign_pathway(ligand_complex, receptor_complex):
    """Return pathway name for an L-R pair, or None if unmapped."""
    genes = get_genes_from_complex(ligand_complex) + \
            get_genes_from_complex(receptor_complex)
    pathways = set()
    for g in genes:
        if g in GENE_TO_PATHWAY:
            pathways.add(GENE_TO_PATHWAY[g])
    if len(pathways) == 1:
        return pathways.pop()
    elif len(pathways) > 1:
        # Prefer the ligand's pathway assignment
        for g in get_genes_from_complex(ligand_complex):
            if g in GENE_TO_PATHWAY:
                return GENE_TO_PATHWAY[g]
        return sorted(pathways)[0]
    return None


def fig_to_base64(fig, dpi=150):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return b64


# ============================================================================
# 1. LOAD & ANNOTATE
# ============================================================================

print("[1] Loading global LIANA results...")
df = pd.read_csv(LIANA_CSV)
print(f"    {len(df):,} total interactions")

# Add v2 and level1 annotations
df["source_v2"]     = df["source"].map(V2_MAP)
df["target_v2"]     = df["target"].map(V2_MAP)
df["source_level1"] = df["source"].map(LEVEL1_MAP).fillna("Unknown")
df["target_level1"] = df["target"].map(LEVEL1_MAP).fillna("Unknown")

# Filter to significant
sig = df[df["magnitude_rank"] <= SIG_THRESH].copy()
print(f"    {len(sig):,} significant (magnitude_rank <= {SIG_THRESH})")

# L-R pair label
sig["lr_pair"] = sig["ligand_complex"] + " -> " + sig["receptor_complex"]

# Assign pathways
sig["pathway"] = sig.apply(
    lambda r: assign_pathway(r["ligand_complex"], r["receptor_complex"]),
    axis=1,
)
n_mapped = sig["pathway"].notna().sum()
print(f"    {n_mapped:,} / {len(sig):,} significant interactions mapped "
      f"to a pathway ({100*n_mapped/len(sig):.1f}%)")

fig_data = {}

# ============================================================================
# 2. ANALYSIS 1: PATHWAY-LEVEL AGGREGATION
# ============================================================================

print("\n[2] Analysis 1: Pathway-level aggregation...")

pseudo = 1  # pseudocount for log2FC

pathway_rows = []
for direction in ["incoming", "outgoing"]:
    if direction == "incoming":
        sub_a = sig[sig["target_v2"] == "SecA"]
        sub_b = sig[sig["target_v2"] == "SecB"]
    else:
        sub_a = sig[sig["source_v2"] == "SecA"]
        sub_b = sig[sig["source_v2"] == "SecB"]

    # Count by pathway
    count_a = sub_a["pathway"].value_counts().rename("n_SecA")
    count_b = sub_b["pathway"].value_counts().rename("n_SecB")
    merged = pd.concat([count_a, count_b], axis=1).fillna(0).astype(int)
    merged["log2fc"] = np.log2((merged["n_SecB"] + pseudo) /
                                (merged["n_SecA"] + pseudo))
    merged["direction"] = direction
    merged.index.name = "pathway"
    pathway_rows.append(merged.reset_index())

pathway_df = pd.concat(pathway_rows, ignore_index=True)
pathway_df = pathway_df[pathway_df["pathway"].notna()].copy()
pathway_df.to_csv(os.path.join(TABLE_DIR,
                  "17d_pathway_aggregation.csv"), index=False)
print(f"    Saved pathway aggregation: {len(pathway_df)} entries")

# ── Figure 2a: Heatmap of pathway activity per pole ──────────
for direction in ["incoming", "outgoing"]:
    sub = pathway_df[pathway_df["direction"] == direction].copy()
    # Order by pathway order, keeping only those present
    ordered = [p for p in PATHWAY_ORDER if p in sub["pathway"].values]
    unmapped_pathways = [p for p in sub["pathway"].values
                         if p not in PATHWAY_ORDER]
    ordered = ordered + sorted(unmapped_pathways)
    sub = sub.set_index("pathway").reindex(ordered).dropna(subset=["n_SecA"])

    if len(sub) == 0:
        continue

    fig, axes = plt.subplots(1, 3, figsize=(10, max(6, len(sub) * 0.3)),
                             gridspec_kw={"width_ratios": [1, 1, 1.2]})

    # Panel 1: SecA counts
    axes[0].barh(range(len(sub)), sub["n_SecA"].values,
                 color=V2_PALETTE["SecA"], edgecolor="#333333", linewidth=0.3)
    axes[0].set_yticks(range(len(sub)))
    axes[0].set_yticklabels(sub.index, fontsize=6)
    axes[0].set_xlabel("Count", fontsize=7)
    axes[0].set_title("SecA", fontsize=9, fontweight="bold",
                       color=V2_PALETTE["SecA"])
    axes[0].invert_xaxis()
    axes[0].invert_yaxis()

    # Panel 2: SecB counts
    axes[1].barh(range(len(sub)), sub["n_SecB"].values,
                 color=V2_PALETTE["SecB"], edgecolor="#333333", linewidth=0.3)
    axes[1].set_yticks(range(len(sub)))
    axes[1].set_yticklabels([])
    axes[1].set_xlabel("Count", fontsize=7)
    axes[1].set_title("SecB", fontsize=9, fontweight="bold",
                       color=V2_PALETTE["SecB"])
    axes[1].invert_yaxis()

    # Panel 3: log2FC bar
    colors = [V2_PALETTE["SecA"] if v < 0 else V2_PALETTE["SecB"]
              for v in sub["log2fc"].values]
    axes[2].barh(range(len(sub)), sub["log2fc"].values,
                 color=colors, edgecolor="#333333", linewidth=0.3)
    axes[2].axvline(0, color="black", linewidth=0.5)
    axes[2].set_yticks(range(len(sub)))
    axes[2].set_yticklabels([])
    axes[2].set_xlabel("log2FC (SecB/SecA)", fontsize=7)
    axes[2].set_title("Differential", fontsize=9, fontweight="bold")
    axes[2].invert_yaxis()

    fig.suptitle(f"Pathway-level communication — {direction}",
                 fontsize=11, fontweight="bold", y=1.02)
    plt.tight_layout()

    fig.savefig(os.path.join(FIG_DIR,
                f"17d_pathway_{direction}.svg"),
                format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR,
                f"17d_pathway_{direction}.pdf"),
                format="pdf", bbox_inches="tight")
    fig_data[f"pathway_{direction}"] = fig_to_base64(fig)

# ── Figure 2b: Pathway log2FC summary (combined directions) ──
pivot = pathway_df.pivot(index="pathway", columns="direction",
                          values="log2fc")
# Order
ordered = [p for p in PATHWAY_ORDER if p in pivot.index]
unordered = [p for p in pivot.index if p not in PATHWAY_ORDER]
pivot = pivot.reindex(ordered + sorted(unordered))
pivot = pivot.dropna(how="all")

if len(pivot) > 0:
    fig, ax = plt.subplots(figsize=(6, max(5, len(pivot) * 0.28)))
    cmap = sns.diverging_palette(30, 220, s=80, l=55, as_cmap=True)
    vmax = max(abs(pivot.min().min()), abs(pivot.max().max()), 1)
    sns.heatmap(pivot, cmap=cmap, center=0, vmin=-vmax, vmax=vmax,
                linewidths=0.5, linecolor="white",
                cbar_kws={"label": "log2FC (SecB/SecA)", "shrink": 0.6},
                ax=ax, annot=True, fmt=".1f", annot_kws={"fontsize": 5})
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Pathway log2FC: SecB vs SecA", fontsize=10,
                  fontweight="bold")
    ax.tick_params(axis="y", labelsize=6)
    ax.tick_params(axis="x", labelsize=7)
    plt.tight_layout()

    fig.savefig(os.path.join(FIG_DIR, "17d_pathway_heatmap.svg"),
                format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, "17d_pathway_heatmap.pdf"),
                format="pdf", bbox_inches="tight")
    fig_data["pathway_heatmap"] = fig_to_base64(fig)


# ============================================================================
# 3. ANALYSIS 2: TME SUBTYPE-LEVEL (LEVEL2) PARTNERS
# ============================================================================

print("\n[3] Analysis 2: TME level2 subtype partners...")

# All non-epithelial level2 cell types in the LIANA data
epi_names = set(V2_MAP.keys())

level2_rows = []
for direction in ["incoming", "outgoing"]:
    if direction == "incoming":
        sub_a = sig[(sig["target_v2"] == "SecA") &
                    (~sig["source"].isin(epi_names))]
        sub_b = sig[(sig["target_v2"] == "SecB") &
                    (~sig["source"].isin(epi_names))]
        partner_col = "source"
    else:
        sub_a = sig[(sig["source_v2"] == "SecA") &
                    (~sig["target"].isin(epi_names))]
        sub_b = sig[(sig["source_v2"] == "SecB") &
                    (~sig["target"].isin(epi_names))]
        partner_col = "target"

    count_a = sub_a[partner_col].value_counts().rename("n_SecA")
    count_b = sub_b[partner_col].value_counts().rename("n_SecB")
    merged = pd.concat([count_a, count_b], axis=1).fillna(0).astype(int)

    merged["total"] = merged["n_SecA"] + merged["n_SecB"]
    merged["prop_SecA"] = merged["n_SecA"] / max(merged["n_SecA"].sum(), 1) * 100
    merged["prop_SecB"] = merged["n_SecB"] / max(merged["n_SecB"].sum(), 1) * 100
    merged["delta_prop"] = merged["prop_SecB"] - merged["prop_SecA"]
    merged["log2fc"] = np.log2((merged["n_SecB"] + pseudo) /
                                (merged["n_SecA"] + pseudo))
    # Map to level1 for coloring
    merged["level1"] = merged.index.map(
        lambda x: LEVEL1_MAP.get(x, "Unknown"))
    merged["direction"] = direction
    merged.index.name = "level2_celltype"
    level2_rows.append(merged.reset_index())

level2_df = pd.concat(level2_rows, ignore_index=True)
level2_df = level2_df[level2_df["level1"] != "Epithelial"].copy()
level2_df.to_csv(os.path.join(TABLE_DIR,
                 "17d_level2_tme_partners.csv"), index=False)
print(f"    Saved: {len(level2_df)} level2 celltype x direction entries")

# ── Figure 3: Level2 TME partners (dot plot style) ──────────
COMPARTMENT_PALETTE = {
    "Epithelial":    "#E6A141",
    "Mesothelial":   "#D4A574",
    "Fibroblast":    "#C4B9A8",
    "Smooth muscle": "#D14E6C",
    "Pericyte":      "#B87A7A",
    "Endothelial":   "#7D4E4E",
    "T/NK cell":     "#87CEFA",
    "B cell":        "#5665B6",
    "Plasma cell":   "#8A5DAF",
    "Macrophage":    "#8FBC8F",
    "DC":            "#2E8B57",
    "Neutrophil":    "#6B8E23",
    "Mast cell":     "#8B9B6B",
}

for direction in ["incoming", "outgoing"]:
    sub = level2_df[level2_df["direction"] == direction].copy()
    # Filter to cell types with at least 5 total interactions
    sub = sub[sub["total"] >= 5].copy()
    sub = sub.sort_values("log2fc")

    if len(sub) == 0:
        continue

    fig, ax = plt.subplots(figsize=(7, max(5, len(sub) * 0.3)))
    colors = [COMPARTMENT_PALETTE.get(r["level1"], "#999999")
              for _, r in sub.iterrows()]

    ax.barh(range(len(sub)), sub["log2fc"].values,
            color=colors, edgecolor="#333333", linewidth=0.3)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(sub["level2_celltype"].values, fontsize=5.5)
    ax.set_xlabel("log2FC (SecB / SecA)", fontsize=8)
    ax.set_title(f"Level2 TME partner differential — {direction}",
                 fontsize=10, fontweight="bold")

    # Add count annotations
    for i, (_, row) in enumerate(sub.iterrows()):
        label = f"A:{int(row['n_SecA'])} B:{int(row['n_SecB'])}"
        x_pos = row["log2fc"]
        ha = "left" if x_pos >= 0 else "right"
        offset = 0.03 if x_pos >= 0 else -0.03
        ax.text(x_pos + offset, i, label, va="center", ha=ha,
                fontsize=4, color="#666666")

    # Legend for compartment colors
    seen = []
    handles = []
    for _, r in sub.iterrows():
        if r["level1"] not in seen:
            seen.append(r["level1"])
            handles.append(Patch(
                facecolor=COMPARTMENT_PALETTE.get(r["level1"], "#999999"),
                label=r["level1"]))
    ax.legend(handles=handles, fontsize=5, frameon=False,
              loc="lower right", ncol=2)

    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR,
                f"17d_level2_partners_{direction}.svg"),
                format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR,
                f"17d_level2_partners_{direction}.pdf"),
                format="pdf", bbox_inches="tight")
    fig_data[f"level2_partners_{direction}"] = fig_to_base64(fig)

# ── Top 10 most differential level2 partners (each direction) ──
for direction in ["incoming", "outgoing"]:
    sub = level2_df[level2_df["direction"] == direction].copy()
    sub = sub[sub["total"] >= 5]
    sub["abs_log2fc"] = sub["log2fc"].abs()
    top = sub.nlargest(10, "abs_log2fc")
    top.to_csv(os.path.join(TABLE_DIR,
               f"17d_top10_diff_level2_{direction}.csv"), index=False)


# ============================================================================
# 4. ANALYSIS 3: DEG CROSS-REFERENCE
# ============================================================================

print("\n[4] Analysis 3: DEG cross-reference...")

deg = pd.read_csv(DEG_CSV)
print(f"    {len(deg):,} DEGs loaded")

# Identify significant DEGs
DEG_PVAL_THRESH = 0.05
DEG_LOG2FC_THRESH = 0.25
deg_sig = deg[
    (deg["pval_adj"] < DEG_PVAL_THRESH) &
    (deg["log2fc"].abs() > DEG_LOG2FC_THRESH)
].copy()
print(f"    {len(deg_sig):,} significant DEGs "
      f"(|log2fc| > {DEG_LOG2FC_THRESH}, pval_adj < {DEG_PVAL_THRESH})")

# Collect all unique ligand and receptor gene symbols from LIANA significant
all_ligand_genes = set()
all_receptor_genes = set()
for _, row in sig.iterrows():
    for g in get_genes_from_complex(row["ligand_complex"]):
        all_ligand_genes.add(g)
    for g in get_genes_from_complex(row["receptor_complex"]):
        all_receptor_genes.add(g)

all_lr_genes = all_ligand_genes | all_receptor_genes

# Cross-reference: which DEGs appear in L-R pairs?
deg_sig["is_ligand"] = deg_sig["gene"].isin(all_ligand_genes)
deg_sig["is_receptor"] = deg_sig["gene"].isin(all_receptor_genes)
deg_sig["is_lr_gene"] = deg_sig["is_ligand"] | deg_sig["is_receptor"]

deg_lr = deg_sig[deg_sig["is_lr_gene"]].copy()

n_lig = deg_lr["is_ligand"].sum()
n_rec = deg_lr["is_receptor"].sum()
n_both = (deg_lr["is_ligand"] & deg_lr["is_receptor"]).sum()
print(f"    {len(deg_lr)} DEGs are communication-active "
      f"({n_lig} ligands, {n_rec} receptors, {n_both} both)")

# For each DEG-LR gene, compute communication log2FC
# (number of significant interactions involving that gene in SecB vs SecA)
deg_comm_rows = []
for _, drow in deg_lr.iterrows():
    gene = drow["gene"]

    # Find interactions where this gene is part of ligand or receptor
    mask_lig = sig["ligand_complex"].apply(
        lambda x: gene in get_genes_from_complex(x))
    mask_rec = sig["receptor_complex"].apply(
        lambda x: gene in get_genes_from_complex(x))
    gene_sig = sig[mask_lig | mask_rec]

    # Count SecA vs SecB involvement (incoming + outgoing)
    n_seca = len(gene_sig[
        (gene_sig["source_v2"] == "SecA") | (gene_sig["target_v2"] == "SecA")])
    n_secb = len(gene_sig[
        (gene_sig["source_v2"] == "SecB") | (gene_sig["target_v2"] == "SecB")])

    comm_log2fc = np.log2((n_secb + pseudo) / (n_seca + pseudo))

    deg_comm_rows.append({
        "gene": gene,
        "deg_log2fc": drow["log2fc"],
        "deg_pval_adj": drow["pval_adj"],
        "deg_score": drow["score"],
        "is_ligand": drow["is_ligand"],
        "is_receptor": drow["is_receptor"],
        "n_interactions_SecA": n_seca,
        "n_interactions_SecB": n_secb,
        "comm_log2fc": comm_log2fc,
        "pathway": GENE_TO_PATHWAY.get(gene, "Other/Unknown"),
    })

deg_comm = pd.DataFrame(deg_comm_rows)
deg_comm = deg_comm.sort_values("deg_score", ascending=False)
deg_comm.to_csv(os.path.join(TABLE_DIR,
                "17d_deg_lr_crossref.csv"), index=False)
print(f"    Saved DEG-LR cross-reference: {len(deg_comm)} genes")

# Summary table
summary_counts = {
    "Total significant DEGs": len(deg_sig),
    "DEGs that are ligands": int(n_lig),
    "DEGs that are receptors": int(n_rec),
    "DEGs that are both L and R": int(n_both),
    "Total DEGs in L-R pairs": len(deg_lr),
    "SecB-enriched DEG-communicators (deg_log2fc>0, comm_log2fc>0)":
        int(((deg_comm["deg_log2fc"] > 0) &
             (deg_comm["comm_log2fc"] > 0)).sum()),
    "SecA-enriched DEG-communicators (deg_log2fc<0, comm_log2fc<0)":
        int(((deg_comm["deg_log2fc"] < 0) &
             (deg_comm["comm_log2fc"] < 0)).sum()),
}
summary_df = pd.DataFrame(list(summary_counts.items()),
                           columns=["Metric", "Value"])
summary_df.to_csv(os.path.join(TABLE_DIR,
                  "17d_deg_lr_summary.csv"), index=False)

# ── Figure 4a: Scatter of DEG log2fc vs communication log2fc ──
if len(deg_comm) > 0:
    fig, ax = plt.subplots(figsize=(6, 5))

    # Color by role
    for _, r in deg_comm.iterrows():
        if r["is_ligand"] and r["is_receptor"]:
            c, marker = "#8856a7", "D"  # purple diamond
        elif r["is_ligand"]:
            c, marker = "#e08214", "o"  # orange circle
        else:
            c, marker = "#1b7837", "s"  # green square
        ax.scatter(r["deg_log2fc"], r["comm_log2fc"],
                   c=c, marker=marker, s=25, alpha=0.7,
                   edgecolors="#333333", linewidths=0.3)

    # Label top genes (high |deg_score| and communication-active)
    top_label = deg_comm.nlargest(25, "deg_score")
    for _, r in top_label.iterrows():
        ax.annotate(r["gene"],
                    (r["deg_log2fc"], r["comm_log2fc"]),
                    fontsize=4.5, alpha=0.8,
                    xytext=(3, 3), textcoords="offset points")

    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.set_xlabel("DEG log2FC (SecB vs SecA)", fontsize=8)
    ax.set_ylabel("Communication log2FC (SecB vs SecA)", fontsize=8)
    ax.set_title("DEG expression vs communication activity",
                 fontsize=10, fontweight="bold")

    handles = [
        plt.scatter([], [], c="#e08214", marker="o", s=25, label="Ligand"),
        plt.scatter([], [], c="#1b7837", marker="s", s=25, label="Receptor"),
        plt.scatter([], [], c="#8856a7", marker="D", s=25,
                    label="Both L+R"),
    ]
    ax.legend(handles=handles, fontsize=6, frameon=False, loc="best")

    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "17d_deg_vs_comm.svg"),
                format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, "17d_deg_vs_comm.pdf"),
                format="pdf", bbox_inches="tight")
    fig_data["deg_vs_comm"] = fig_to_base64(fig)

# ── Figure 4b: Top DEG-communicators bar chart ──────────────
# Show top 30 DEGs by |deg_score| that are L-R genes
top_deg_comm = deg_comm.nlargest(30, "deg_score").copy()
if len(top_deg_comm) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(10, 7),
                             gridspec_kw={"width_ratios": [1, 1]})

    # Left: DEG log2FC
    top_sorted = top_deg_comm.sort_values("deg_log2fc")
    colors_deg = [V2_PALETTE["SecB"] if v > 0 else V2_PALETTE["SecA"]
                  for v in top_sorted["deg_log2fc"]]
    axes[0].barh(range(len(top_sorted)), top_sorted["deg_log2fc"].values,
                 color=colors_deg, edgecolor="#333333", linewidth=0.3)
    axes[0].set_yticks(range(len(top_sorted)))
    labels = []
    for _, r in top_sorted.iterrows():
        tag = ""
        if r["is_ligand"] and r["is_receptor"]:
            tag = " [L+R]"
        elif r["is_ligand"]:
            tag = " [L]"
        else:
            tag = " [R]"
        labels.append(f"{r['gene']}{tag}")
    axes[0].set_yticklabels(labels, fontsize=5.5)
    axes[0].axvline(0, color="black", linewidth=0.5)
    axes[0].set_xlabel("DEG log2FC", fontsize=7)
    axes[0].set_title("Expression (SecB vs SecA)", fontsize=9,
                       fontweight="bold")

    # Right: Communication log2FC
    colors_comm = [V2_PALETTE["SecB"] if v > 0 else V2_PALETTE["SecA"]
                   for v in top_sorted["comm_log2fc"]]
    axes[1].barh(range(len(top_sorted)), top_sorted["comm_log2fc"].values,
                 color=colors_comm, edgecolor="#333333", linewidth=0.3)
    axes[1].set_yticks(range(len(top_sorted)))
    axes[1].set_yticklabels([])
    axes[1].axvline(0, color="black", linewidth=0.5)
    axes[1].set_xlabel("Communication log2FC", fontsize=7)
    axes[1].set_title("Communication (SecB vs SecA)", fontsize=9,
                       fontweight="bold")

    # Add count annotations
    for i, (_, row) in enumerate(top_sorted.iterrows()):
        label = f"A:{int(row['n_interactions_SecA'])} " \
                f"B:{int(row['n_interactions_SecB'])}"
        x_pos = row["comm_log2fc"]
        ha = "left" if x_pos >= 0 else "right"
        offset = 0.02 if x_pos >= 0 else -0.02
        axes[1].text(x_pos + offset, i, label, va="center", ha=ha,
                     fontsize=4, color="#666666")

    fig.suptitle("Top DEG-communicators: expression vs communication",
                 fontsize=11, fontweight="bold", y=1.01)
    plt.tight_layout()

    fig.savefig(os.path.join(FIG_DIR, "17d_top_deg_communicators.svg"),
                format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, "17d_top_deg_communicators.pdf"),
                format="pdf", bbox_inches="tight")
    fig_data["top_deg_comm"] = fig_to_base64(fig)

# ── Figure 4c: Interaction counts per pole for top DEG-LR genes ──
top15 = deg_comm.nlargest(15, "deg_score")
if len(top15) > 0:
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(top15))
    w = 0.35
    ax.bar(x - w/2, top15["n_interactions_SecA"].values, w,
           color=V2_PALETTE["SecA"], edgecolor="#333333", linewidth=0.3,
           label="SecA")
    ax.bar(x + w/2, top15["n_interactions_SecB"].values, w,
           color=V2_PALETTE["SecB"], edgecolor="#333333", linewidth=0.3,
           label="SecB")
    ax.set_xticks(x)
    gene_labels = []
    for _, r in top15.iterrows():
        fc_str = f"({r['deg_log2fc']:+.1f})"
        gene_labels.append(f"{r['gene']} {fc_str}")
    ax.set_xticklabels(gene_labels, rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("Significant interactions", fontsize=8)
    ax.set_title("Interaction counts for top DEG-communicators",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, frameon=False)
    plt.tight_layout()

    fig.savefig(os.path.join(FIG_DIR, "17d_deg_interaction_counts.svg"),
                format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, "17d_deg_interaction_counts.pdf"),
                format="pdf", bbox_inches="tight")
    fig_data["deg_interaction_counts"] = fig_to_base64(fig)


# ============================================================================
# 4b. CLINICAL CONTEXT FIGURES
# ============================================================================

print("\n[4b] Clinical context figures...")

# ── Figure 5a: Niche-dependency pathway classification ──────────
# Classify pathways as TME-dependent (require stromal/immune partners)
# vs autonomous/autocrine (can function without TME)
niche_dependent = [
    "Collagen/ECM", "TNFa", "NOTCH", "FGF", "TGFb", "Chemokine",
    "PDGF", "VEGF", "Complement", "MHC/Antigen presentation",
    "Immune checkpoint", "Semaphorin", "HGF/MET",
]
niche_independent = [
    "WNT", "SLIT/ROBO", "JAK-STAT/Interleukin", "Ephrin",
    "Integrin", "EGFR", "Galectin", "Tetraspanin",
]

# Count interactions per pole for each category
dep_rows = []
for cat_name, cat_pathways in [("TME-dependent", niche_dependent),
                                ("Autonomous/autocrine", niche_independent)]:
    for direction in ["incoming", "outgoing"]:
        sub_pw = pathway_df[
            (pathway_df["direction"] == direction) &
            (pathway_df["pathway"].isin(cat_pathways))
        ]
        n_a = sub_pw["n_SecA"].sum()
        n_b = sub_pw["n_SecB"].sum()
        dep_rows.append({"category": cat_name, "direction": direction,
                         "SecA": n_a, "SecB": n_b})

dep_summary = pd.DataFrame(dep_rows)

fig, axes = plt.subplots(1, 2, figsize=(9, 4))

for i, direction in enumerate(["incoming", "outgoing"]):
    sub = dep_summary[dep_summary["direction"] == direction]
    x = np.arange(len(sub))
    w = 0.35
    axes[i].bar(x - w/2, sub["SecA"].values, w,
                color=V2_PALETTE["SecA"], edgecolor="#333333", linewidth=0.3,
                label="SecA")
    axes[i].bar(x + w/2, sub["SecB"].values, w,
                color=V2_PALETTE["SecB"], edgecolor="#333333", linewidth=0.3,
                label="SecB")
    axes[i].set_xticks(x)
    axes[i].set_xticklabels(sub["category"].values, fontsize=7)
    axes[i].set_ylabel("Significant interactions", fontsize=7)
    axes[i].set_title(f"{direction.capitalize()} pathways", fontsize=9,
                      fontweight="bold")
    axes[i].legend(fontsize=6, frameon=False)

    # Add count labels on bars
    for j, (_, row) in enumerate(sub.iterrows()):
        axes[i].text(j - w/2, row["SecA"] + 5, str(int(row["SecA"])),
                     ha="center", fontsize=5.5, color=V2_PALETTE["SecA"])
        axes[i].text(j + w/2, row["SecB"] + 5, str(int(row["SecB"])),
                     ha="center", fontsize=5.5, color=V2_PALETTE["SecB"])

fig.suptitle("Communication by niche dependency", fontsize=11,
             fontweight="bold", y=1.02)
plt.tight_layout()

fig.savefig(os.path.join(FIG_DIR, "17d_niche_dependency.svg"),
            format="svg", bbox_inches="tight")
fig.savefig(os.path.join(FIG_DIR, "17d_niche_dependency.pdf"),
            format="pdf", bbox_inches="tight")
fig_data["niche_dependency"] = fig_to_base64(fig)

# ── Figure 5b: SecB autonomous communication profile ──────────
# Show SecB-enriched pathway log2FCs highlighting autonomous pathways
outgoing = pathway_df[pathway_df["direction"] == "outgoing"].copy()
outgoing = outgoing[outgoing["pathway"].notna()].set_index("pathway")
# Select key pathways to highlight
highlight_pathways = ["WNT", "SLIT/ROBO", "JAK-STAT/Interleukin", "Ephrin",
                      "Collagen/ECM", "TNFa", "NOTCH", "FGF", "TGFb",
                      "Chemokine", "Integrin", "EGFR"]
sub_out = outgoing.reindex([p for p in highlight_pathways
                             if p in outgoing.index]).copy()

if len(sub_out) > 0:
    sub_out = sub_out.sort_values("log2fc")
    fig, ax = plt.subplots(figsize=(7, max(4, len(sub_out) * 0.4)))

    bar_colors = []
    for p in sub_out.index:
        if p in niche_independent:
            bar_colors.append("#4393C3")  # blue for autonomous
        else:
            bar_colors.append("#D6604D")  # red for TME-dependent

    ax.barh(range(len(sub_out)), sub_out["log2fc"].values,
            color=bar_colors, edgecolor="#333333", linewidth=0.3)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_yticks(range(len(sub_out)))
    ylabels = []
    for p in sub_out.index:
        tag = " *" if p in niche_independent else ""
        ylabels.append(f"{p}{tag}")
    ax.set_yticklabels(ylabels, fontsize=6.5)
    ax.set_xlabel("log2FC outgoing (SecB / SecA)", fontsize=7)
    ax.set_title("Outgoing pathway profile: autonomous vs TME-dependent",
                 fontsize=9, fontweight="bold")

    handles = [
        Patch(facecolor="#4393C3", label="Autonomous/autocrine"),
        Patch(facecolor="#D6604D", label="TME-dependent"),
    ]
    ax.legend(handles=handles, fontsize=6, frameon=False, loc="lower right")

    # Count annotations
    for i, (p, row) in enumerate(sub_out.iterrows()):
        label = f"A:{int(row['n_SecA'])} B:{int(row['n_SecB'])}"
        x_pos = row["log2fc"]
        ha = "left" if x_pos >= 0 else "right"
        offset = 0.05 if x_pos >= 0 else -0.05
        ax.text(x_pos + offset, i, label, va="center", ha=ha,
                fontsize=4.5, color="#666666")

    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "17d_autonomous_profile.svg"),
                format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, "17d_autonomous_profile.pdf"),
                format="pdf", bbox_inches="tight")
    fig_data["autonomous_profile"] = fig_to_base64(fig)

# ── Figure 5c: Concordance/discordance pie ──────────────────────
if len(deg_comm) > 0:
    concordant_secb = ((deg_comm["deg_log2fc"] > 0) &
                       (deg_comm["comm_log2fc"] > 0)).sum()
    concordant_seca = ((deg_comm["deg_log2fc"] < 0) &
                       (deg_comm["comm_log2fc"] < 0)).sum()
    discordant_up_down = ((deg_comm["deg_log2fc"] > 0) &
                          (deg_comm["comm_log2fc"] < 0)).sum()
    discordant_down_up = ((deg_comm["deg_log2fc"] < 0) &
                          (deg_comm["comm_log2fc"] > 0)).sum()

    cat_counts = {
        f"Concordant SecB\n(expr+ comm+)\nn={concordant_secb}": concordant_secb,
        f"Concordant SecA\n(expr\u2212 comm\u2212)\nn={concordant_seca}": concordant_seca,
        f"Discordant\n(SecB expr, SecA comm)\nn={discordant_up_down}": discordant_up_down,
        f"Discordant\n(SecA expr, SecB comm)\nn={discordant_down_up}": discordant_down_up,
    }

    fig, ax = plt.subplots(figsize=(5, 5))
    cat_colors = ["#B8741A", "#E6A141", "#7570B3", "#A6CEE3"]
    vals = list(cat_counts.values())
    labels = list(cat_counts.keys())
    wedges, texts, autotexts = ax.pie(
        vals, labels=labels, autopct="%1.0f%%", colors=cat_colors,
        startangle=90, textprops={"fontsize": 6},
        wedgeprops={"edgecolor": "#333333", "linewidth": 0.5})
    for t in autotexts:
        t.set_fontsize(7)
        t.set_fontweight("bold")
    ax.set_title("DEG-communicator concordance", fontsize=9,
                 fontweight="bold")
    plt.tight_layout()

    fig.savefig(os.path.join(FIG_DIR, "17d_concordance_pie.svg"),
                format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, "17d_concordance_pie.pdf"),
                format="pdf", bbox_inches="tight")
    fig_data["concordance_pie"] = fig_to_base64(fig)

# ── Figure 5d: TME engagement ratio ──────────────────────────────
# Total significant interactions per pole across all TME compartments
tme_incoming_a = len(sig[(sig["target_v2"] == "SecA") &
                          (sig["source_level1"] != "Epithelial")])
tme_incoming_b = len(sig[(sig["target_v2"] == "SecB") &
                          (sig["source_level1"] != "Epithelial")])
tme_outgoing_a = len(sig[(sig["source_v2"] == "SecA") &
                          (sig["target_level1"] != "Epithelial")])
tme_outgoing_b = len(sig[(sig["source_v2"] == "SecB") &
                          (sig["target_level1"] != "Epithelial")])

fig, ax = plt.subplots(figsize=(5, 4))
categories = ["Incoming\n(TME → Epi)", "Outgoing\n(Epi → TME)"]
seca_vals = [tme_incoming_a, tme_outgoing_a]
secb_vals = [tme_incoming_b, tme_outgoing_b]
x = np.arange(len(categories))
w = 0.35

bars_a = ax.bar(x - w/2, seca_vals, w, color=V2_PALETTE["SecA"],
                edgecolor="#333333", linewidth=0.3, label="SecA")
bars_b = ax.bar(x + w/2, secb_vals, w, color=V2_PALETTE["SecB"],
                edgecolor="#333333", linewidth=0.3, label="SecB")
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=8)
ax.set_ylabel("Significant TME interactions", fontsize=7)
ax.set_title("TME engagement: SecA vs SecB", fontsize=9, fontweight="bold")
ax.legend(fontsize=7, frameon=False)

# Add ratio annotations
for j in range(len(categories)):
    ratio = seca_vals[j] / max(secb_vals[j], 1)
    ax.text(x[j], max(seca_vals[j], secb_vals[j]) + 20,
            f"{ratio:.1f}x", ha="center", fontsize=7, fontweight="bold",
            color="#555555")

# Add count labels
for bar in bars_a:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 3,
            str(int(bar.get_height())), ha="center", fontsize=5.5,
            color=V2_PALETTE["SecA"])
for bar in bars_b:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 3,
            str(int(bar.get_height())), ha="center", fontsize=5.5,
            color=V2_PALETTE["SecB"])

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "17d_tme_engagement.svg"),
            format="svg", bbox_inches="tight")
fig.savefig(os.path.join(FIG_DIR, "17d_tme_engagement.pdf"),
            format="pdf", bbox_inches="tight")
fig_data["tme_engagement"] = fig_to_base64(fig)


# ============================================================================
# 5. HTML REPORT
# ============================================================================

print("\n[5] Generating HTML report...")


def img_tag(key, width="100%"):
    if key in fig_data:
        return (f'<img src="data:image/png;base64,{fig_data[key]}" '
                f'style="width:{width};">')
    return "<p><em>Figure not generated</em></p>"


def df_to_html(dataframe, max_rows=40, float_fmt=".3f"):
    """Convert a dataframe to a styled HTML table."""
    if len(dataframe) > max_rows:
        dataframe = dataframe.head(max_rows)
    rows_html = ""
    cols = dataframe.columns.tolist()
    header = "".join(f"<th>{c}</th>" for c in cols)

    for _, r in dataframe.iterrows():
        cells = ""
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                if abs(v) < 0.001 and v != 0:
                    cells += f"<td>{v:.2e}</td>"
                else:
                    cells += f"<td>{v:{float_fmt}}</td>"
            elif isinstance(v, bool):
                cells += f"<td>{'Yes' if v else ''}</td>"
            else:
                cells += f"<td>{v}</td>"
        rows_html += f"<tr>{cells}</tr>"

    return f"<table><tr>{header}</tr>{rows_html}</table>"


# Pathway table for HTML
pathway_html_sections = ""
for direction in ["incoming", "outgoing"]:
    sub = pathway_df[pathway_df["direction"] == direction].copy()
    sub = sub.sort_values("log2fc", ascending=False)
    sub_display = sub[["pathway", "n_SecA", "n_SecB", "log2fc"]].copy()
    sub_display.columns = ["Pathway", "n SecA", "n SecB", "log2FC"]
    pathway_html_sections += f"""
    <h3>Pathway counts — {direction}</h3>
    {df_to_html(sub_display, float_fmt=".2f")}
    """

# Level2 tables for HTML
level2_html_sections = ""
for direction in ["incoming", "outgoing"]:
    sub = level2_df[level2_df["direction"] == direction].copy()
    sub = sub[sub["total"] >= 5].sort_values("log2fc", ascending=False)
    sub_display = sub[["level2_celltype", "level1", "n_SecA", "n_SecB",
                        "log2fc"]].copy()
    sub_display.columns = ["Level2 cell type", "Compartment",
                           "n SecA", "n SecB", "log2FC"]
    level2_html_sections += f"""
    <h3>Level2 TME partners — {direction}</h3>
    {df_to_html(sub_display, float_fmt=".2f")}
    """

# DEG-LR table for HTML
deg_table_display = deg_comm.head(40)[
    ["gene", "deg_log2fc", "deg_pval_adj", "is_ligand", "is_receptor",
     "n_interactions_SecA", "n_interactions_SecB", "comm_log2fc",
     "pathway"]
].copy()
deg_table_display.columns = [
    "Gene", "DEG log2FC", "DEG p_adj", "Ligand?", "Receptor?",
    "n SecA", "n SecB", "Comm log2FC", "Pathway"
]

# Summary stats HTML
summary_html = ""
for _, r in summary_df.iterrows():
    summary_html += f"<li><strong>{r['Metric']}:</strong> {r['Value']}</li>\n"

html = f"""<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>17d — SecA vs SecB Communication Mechanisms</title>
<style>
    body {{ font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto;
           padding: 20px; background: #fafafa; color: #333; }}
    h1 {{ color: #333; border-bottom: 2px solid #E6A141; padding-bottom: 8px; }}
    h2 {{ color: #555; margin-top: 30px; border-bottom: 1px solid #ddd;
          padding-bottom: 5px; }}
    h3 {{ color: #777; }}
    .fig {{ text-align: center; margin: 15px 0; }}
    .fig img {{ border: 1px solid #ddd; border-radius: 4px; }}
    table {{ border-collapse: collapse; width: 100%; margin: 10px 0;
             font-size: 11px; }}
    th, td {{ border: 1px solid #ddd; padding: 4px 6px; text-align: left; }}
    th {{ background: #f0f0f0; font-weight: bold; }}
    tr:nth-child(even) {{ background: #f9f9f9; }}
    .schema {{ background: #fff; padding: 12px; border-left: 4px solid #E6A141;
               margin: 10px 0; }}
    .two-col {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
    ul {{ margin: 5px 0; }}
</style>
</head><body>

<h1>17d — SecA vs SecB Communication Mechanisms</h1>
<p>Deeper analyses of global LIANA+ results through the v2 epithelial schema:
pathway-level aggregation, level2 TME partner breakdown, and DEG cross-reference.
Significance threshold: magnitude_rank &le; {SIG_THRESH}.</p>

<div class="schema">
<strong>v2 Schema:</strong>
SecA = NMF Factor_2 &lt; p50 | Intermediate = p50-p75 |
SecB = Adaptive secretory epithelial |
Stress-response and Ciliated kept separate
</div>

<hr>

<h2>Analysis 1: Pathway-Level Aggregation</h2>
<p>L-R pairs mapped to canonical signalling pathways. log2FC computed as
log2((n_SecB + 1) / (n_SecA + 1)). Positive = SecB-enriched, negative = SecA-enriched.</p>

<h3>1a. Pathway activity per pole (incoming)</h3>
<div class="fig">{img_tag("pathway_incoming", "85%")}</div>

<h3>1b. Pathway activity per pole (outgoing)</h3>
<div class="fig">{img_tag("pathway_outgoing", "85%")}</div>

<h3>1c. Pathway heatmap (combined)</h3>
<div class="fig">{img_tag("pathway_heatmap", "55%")}</div>

{pathway_html_sections}

<hr>

<h2>Analysis 2: TME Level2 Subtype Partners</h2>
<p>Communication broken down to individual level2 TME cell types rather than
level1 compartments. Minimum 5 interactions required for display.</p>

<h3>2a. Level2 partner differential (incoming)</h3>
<div class="fig">{img_tag("level2_partners_incoming", "70%")}</div>

<h3>2b. Level2 partner differential (outgoing)</h3>
<div class="fig">{img_tag("level2_partners_outgoing", "70%")}</div>

{level2_html_sections}

<hr>

<h2>Analysis 3: DEG Cross-Reference</h2>
<p>Cross-reference of SecA-vs-SecB differentially expressed genes with LIANA L-R pairs.
DEG significance: |log2FC| &gt; {DEG_LOG2FC_THRESH}, adjusted p-value &lt; {DEG_PVAL_THRESH}.
Positive DEG log2FC = SecB-enriched expression.</p>

<h3>Summary</h3>
<ul>
{summary_html}
</ul>

<h3>3a. DEG expression vs communication activity</h3>
<div class="fig">{img_tag("deg_vs_comm", "55%")}</div>
<p><em>Each point is a significant DEG that also participates in L-R signalling.
x-axis: differential expression (SecB vs SecA); y-axis: differential communication
activity (SecB vs SecA). Concordant genes (same quadrant) suggest expression
changes translate to communication changes.</em></p>

<h3>3b. Top DEG-communicators</h3>
<div class="fig">{img_tag("top_deg_comm", "85%")}</div>

<h3>3c. Interaction counts for top DEG-communicators</h3>
<div class="fig">{img_tag("deg_interaction_counts", "70%")}</div>

<h3>3d. DEG-communicator table (top 40 by score)</h3>
{df_to_html(deg_table_display, float_fmt=".3f")}

<hr>

<h2>Integrated Interpretation</h2>

<div class="schema" style="border-left-color: #B8741A;">
<strong>Overview:</strong> SecA and SecB epithelial cells occupy fundamentally different
communication niches within the HGSC tumour microenvironment. SecA is a high-connectivity,
outward-signalling state that actively remodels its surroundings, while SecB is a lower-connectivity
state with selective, specialised signalling axes.
</div>

<h3>I. SecA: The microenvironment-engaged progenitor state</h3>
<p>SecA cells (NMF Factor_2 &lt; p50) dominate both incoming and outgoing
communication, with ~3&times; more significant interactions than SecB. This is not simply a cell
number artefact &mdash; it reflects broad engagement across multiple signalling modalities:</p>
<ul>
    <li><strong>ECM / Collagen axis:</strong> SecA is heavily enriched for collagen- and ECM-based
    interactions (Collagen/ECM, SPP1/Osteopontin), consistent with its role in matrix remodelling
    and its association with mesenchymal-like features.</li>
    <li><strong>TNF&alpha; and NOTCH (incoming):</strong> SecA receives TNF&alpha; and NOTCH signals
    from the TME while SecB receives essentially none. This suggests SecA is the target of
    inflammatory and differentiation cues from the stroma &mdash; potentially maintaining its
    progenitor-like state or licensing proliferation.</li>
    <li><strong>FGF and TGF&beta; (outgoing):</strong> SecA sends FGF and TGF&beta; signals outward;
    SecB sends very little. These are classic stromal remodelling pathways, consistent with SecA
    actively shaping the fibroblast and immune compartments.</li>
    <li><strong>Chemokine signalling:</strong> CXCR4, a SecA-upregulated gene, is one of the
    strongest concordant DEG-communicators (DEG log2FC = &minus;1.02, comm log2FC = &minus;1.31).
    SecA cells are both expressing and using CXCR4 for active chemokine communication, potentially
    driving migration and homing behaviours.</li>
</ul>

<h3>II. SecB: Selective signalling in a differentiated state</h3>
<p>SecB (Adaptive secretory epithelial) communicates through a narrower but functionally distinct
set of pathways:</p>
<ul>
    <li><strong>WNT (outgoing):</strong> The most striking SecB-enriched pathway. SecB sends WNT
    signals at ~40&times; the rate of SecA (78 vs 2 interactions), driven by receptors FZD1 and
    FZD7 which are both SecB-upregulated DEGs. This is a highly concordant expression-communication
    axis and may reflect autocrine WNT maintenance of the differentiated state or paracrine
    signalling to the stroma.</li>
    <li><strong>JAK-STAT/Interleukin (outgoing):</strong> SecB sends more JAK-STAT signals,
    including via IL15RA, suggesting a role in immune modulation distinct from SecA&rsquo;s
    TNF&alpha;-driven immune engagement.</li>
    <li><strong>SLIT/ROBO (outgoing):</strong> SecB uniquely sends SLIT/ROBO signals (7 vs 0).
    SLIT/ROBO is a guidance cue pathway involved in axon repulsion and, in cancer, tumour cell
    migration and immune cell guidance.</li>
    <li><strong>Ephrin (incoming):</strong> SecB receives more ephrin signals than SecA, suggesting
    cell-contact-dependent signalling that may maintain spatial organisation within tumour
    architecture.</li>
</ul>

<h3>III. The expression-communication discordance</h3>
<p>A key finding is the large number of <strong>discordant</strong> genes: DEGs that are
upregulated in SecB but whose communication activity is enriched in SecA. Of 86
DEG-communicators, the majority (65%) are discordant. The most striking examples:</p>
<ul>
    <li><strong>SLPI</strong> (SecB log2FC = +1.01, comm log2FC = &minus;2.63): SecB expresses
    SLPI at high levels but the ligand-receptor interaction scores are far weaker. SLPI in SecB
    may function primarily as a protease inhibitor (autocrine/protective) rather than as a
    paracrine signalling ligand.</li>
    <li><strong>TNFSF9</strong> (SecB log2FC = +1.23, comm = SecA-only): A TNF superfamily member
    expressed by SecB but only functionally communicating through SecA. SecB may lack the
    appropriate receptor context on its TME partners.</li>
    <li><strong>EGFR, ITGA2, ITGA3, ITGA5, ITGB1</strong> (all SecB-upregulated, all SecA-enriched
    communication): The integrin receptors are broadly upregulated in SecB but far more communication-active
    in SecA. This dissociation suggests that integrin signalling depends on the <em>partner cell</em>
    context &mdash; SecA&rsquo;s ECM-remodelling partners provide the ligands that activate these
    receptors, while SecB&rsquo;s TME niche does not.</li>
</ul>

<h3>IV. TME partner specificity</h3>
<p>At the level1 compartment level, SecA and SecB engage the same TME partners in roughly similar
proportions. The biological specificity emerges at level2:</p>
<ul>
    <li>All TME subtypes show SecA-preferring log2FC (negative), reflecting SecA&rsquo;s higher
    overall connectivity. The <em>relative</em> enrichment within each pole&rsquo;s interaction
    budget is more informative.</li>
    <li>The lack of strongly SecB-preferring TME subtypes reinforces SecB&rsquo;s characterisation
    as a more self-contained, less TME-engaged state.</li>
</ul>

<h3>V. Synthesis</h3>
<div class="schema" style="border-left-color: #E6A141;">
<p><strong>SecA</strong> functions as a <em>microenvironment architect</em> &mdash; a high-connectivity
progenitor state that sends ECM/FGF/TGF&beta; signals to remodel the stroma, receives TNF&alpha;/NOTCH
cues to maintain its phenotype, and actively engages CXCR4-chemokine axes for migration. Its broad
signalling repertoire is consistent with a plastic, stem-like state that depends on and shapes its niche.</p>

<p><strong>SecB</strong> functions as a <em>specialised signaller</em> &mdash; a lower-connectivity
differentiated state with selective WNT, SLIT/ROBO, and JAK-STAT outgoing axes. The widespread
expression-communication discordance (high receptor expression, low interaction scores) suggests SecB
has decoupled from many TME signalling circuits. Its high expression of integrins, SLPI, and EGFR
may serve cell-autonomous functions (adhesion, protease defence, survival signalling) rather than
paracrine communication.</p>

<p>This communication architecture directly supports the transcriptomic and functional characterisation
from Figure 2: SecA&rsquo;s pathway enrichment for proliferation, OXPHOS, and MYC targets aligns with
its high TME engagement and receipt of growth/differentiation signals, while SecB&rsquo;s Hallmark
enrichment for EMT and hypoxia aligns with its selective WNT signalling and relative communication
independence.</p>
</div>

<hr>

<h2>VI. Clinical Context: Niche Independence and Treatment Response</h2>

<div class="schema" style="border-left-color: #B8741A;">
<strong>Key insight:</strong> SecB&rsquo;s communication architecture &mdash; low TME engagement,
autocrine WNT signalling, and cell-autonomous receptor functions &mdash; directly explains two
clinically observed phenotypes: enrichment in ascites and maintenance post-chemotherapy.
</div>

<h3>6a. TME engagement: SecA vs SecB</h3>
<div class="fig">{img_tag("tme_engagement", "45%")}</div>
<p>SecA engages the TME at several-fold the rate of SecB in both signalling directions.
This quantifies the fundamental difference: SecA is TME-coupled, SecB is TME-decoupled.</p>

<h3>6b. Niche-dependency classification</h3>
<div class="fig">{img_tag("niche_dependency", "75%")}</div>
<p>Pathways classified as TME-dependent (requiring stromal/immune partners: Collagen/ECM, TNF&alpha;,
NOTCH, FGF, TGF&beta;, Chemokine, PDGF, VEGF, Complement, MHC, Immune checkpoint, Semaphorin, HGF/MET)
vs autonomous/autocrine (can function without TME: WNT, SLIT/ROBO, JAK-STAT, Ephrin, Integrin, EGFR,
Galectin, Tetraspanin). SecA dominates TME-dependent signalling; SecB&rsquo;s communication budget is
more balanced between the two categories.</p>

<h3>6c. Outgoing pathway profile</h3>
<div class="fig">{img_tag("autonomous_profile", "60%")}</div>
<p>SecB-enriched outgoing pathways (positive log2FC) are concentrated in autonomous/autocrine categories
(blue), while SecA-enriched pathways (negative log2FC) are concentrated in TME-dependent categories (red).
This asymmetry is the molecular basis for niche independence.</p>

<h3>6d. Expression&ndash;communication discordance</h3>
<div class="fig">{img_tag("concordance_pie", "35%")}</div>
<p>The majority of DEG-communicators are discordant: genes upregulated in SecB but with communication
activity enriched in SecA. This means SecB expresses the receptors (EGFR, integrins, SLPI) but does not
engage them in paracrine signalling &mdash; they serve cell-autonomous functions instead.</p>

<h3>Connecting to ascites enrichment</h3>
<div class="schema" style="border-left-color: #4393C3;">
<p>Ascites is a niche-depleted environment: free-floating cells detached from stroma, ECM, and
organised tissue architecture. SecA&rsquo;s communication repertoire &mdash; Collagen/ECM remodelling,
CXCR4-mediated chemotaxis, FGF/TGF&beta; stromal signalling, and TNF&alpha;/NOTCH reception &mdash;
all require tissue-resident partners that are absent in ascites. Without these inputs, SecA cells
cannot maintain their signalling network.</p>

<p>SecB, by contrast, is pre-adapted to niche-depleted conditions:</p>
<ul>
    <li><strong>Autocrine WNT signalling</strong> (FZD1, FZD7, WNT ligands) can sustain stemness
    and survival without stromal co-stimulation.</li>
    <li><strong>SLIT/ROBO outgoing signals</strong> function in cell-autonomous migration guidance.</li>
    <li><strong>Cell-autonomous receptor usage</strong> (integrins for adhesion/survival signalling,
    SLPI for protease defence, EGFR for survival) provides protective functions independent of
    paracrine ligand availability.</li>
    <li><strong>Low TME dependency</strong> means SecB does not require the fibroblast/macrophage/
    endothelial interactions that are unavailable in ascites fluid.</li>
</ul>
<p>This communication profile predicts selective survival of SecB in ascites &mdash; not because SecB
actively thrives there, but because SecA cannot maintain its signalling requirements.</p>
</div>

<h3>Connecting to post-chemotherapy maintenance</h3>
<div class="schema" style="border-left-color: #D6604D;">
<p>Chemotherapy disrupts the TME through multiple mechanisms: direct cytotoxicity to stromal and immune
cells, destruction of vascular architecture, and ECM damage. This effectively recapitulates the
niche-depleted state of ascites in solid tumour sites. The same communication-independence that
advantages SecB in ascites also advantages it post-chemotherapy:</p>
<ul>
    <li><strong>Stromal destruction</strong> removes the fibroblast/macrophage partners that SecA
    depends on for Collagen/ECM, FGF, and TGF&beta; signalling.</li>
    <li><strong>Loss of TNF&alpha;/NOTCH inputs</strong> from damaged immune cells deprives SecA of
    the inflammatory and differentiation cues it requires.</li>
    <li><strong>SecB&rsquo;s WNT autocrine loop</strong> is resilient to niche destruction &mdash;
    it only requires the cancer cell itself.</li>
    <li><strong>Cell-autonomous survival</strong> via SLPI (protease defence), integrins (anchorage-
    independent survival), and EGFR (growth factor-independent activation) provides SecB with
    intrinsic resistance mechanisms.</li>
    <li>SecB&rsquo;s <strong>Hypoxia and EMT Hallmark enrichment</strong> (from Figure 2 functional
    characterisation) further supports survival in the damaged, hypoxic post-chemo microenvironment.</li>
</ul>
<p>The model is thus: chemotherapy selectively eliminates SecA by destroying the TME niche it depends
on, while SecB persists through communication-autonomous survival programs. This is not classical drug
resistance (efflux pumps, target mutations) but rather <em>ecological resistance</em> &mdash; survival
through pre-existing independence from the niche that chemotherapy destroys.</p>
</div>

<h3>Model summary</h3>
<div class="schema" style="border-left-color: #333333;">
<table style="border: none; font-size: 12px; width: 100%;">
<tr style="background: #f0f0f0;">
    <th style="width: 25%;">Feature</th>
    <th style="width: 37.5%; color: #E6A141;">SecA (niche-dependent)</th>
    <th style="width: 37.5%; color: #B8741A;">SecB (niche-independent)</th>
</tr>
<tr><td><strong>TME connectivity</strong></td>
    <td>High (~3&times; more sig. interactions)</td>
    <td>Low (selective pathways only)</td></tr>
<tr><td><strong>Key outgoing</strong></td>
    <td>ECM, FGF, TGF&beta; (stromal remodelling)</td>
    <td>WNT, SLIT/ROBO, JAK-STAT (autocrine/selective)</td></tr>
<tr><td><strong>Key incoming</strong></td>
    <td>TNF&alpha;, NOTCH (inflammatory/differentiation)</td>
    <td>Ephrin (cell-contact only)</td></tr>
<tr><td><strong>Receptor usage</strong></td>
    <td>Paracrine (ligand-dependent)</td>
    <td>Cell-autonomous (ligand-independent)</td></tr>
<tr><td><strong>Ascites survival</strong></td>
    <td>Disadvantaged (niche absent)</td>
    <td>Advantaged (niche not required)</td></tr>
<tr><td><strong>Post-chemo survival</strong></td>
    <td>Disadvantaged (niche destroyed)</td>
    <td>Advantaged (autonomous programs)</td></tr>
<tr><td><strong>Resistance type</strong></td>
    <td>&mdash;</td>
    <td>Ecological (niche-independent persistence)</td></tr>
</table>
</div>

<hr>
<p style="font-size:11px; color:#999;">
Generated by 17d_secA_secB_mechanisms_nmf.py &middot;
HGSC Atlas &middot; Cook Lab 2026
</p>
</body></html>
"""

html_path = os.path.join(OUT_DIR,
                         "17d_secA_secB_communication_mechanisms.html")
with open(html_path, "w") as f:
    f.write(html)
print(f"    Saved: {html_path}")

print("\nDone!")
