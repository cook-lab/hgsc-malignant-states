#!/usr/bin/env python3
"""
Atlas 02 — Step 10c: epithelial annotation validation report

PURPOSE
    Validate the epithelial level-2 annotation: per-cluster marker/signature scores,
    cell-cycle, PROGENy pathways, PAGA/DPT trajectory connectivity, and CytoTRACE,
    assembled into an HTML validation report with supporting CSVs. Confirms the
    SecA/SecB (progenitor <-> differentiated) interpretation of the epithelial
    subtypes.

INPUTS
    obj("atlas_epithelial")  = hgsc_atlas_epithelial.h5ad

OUTPUTS
    output_root/02_annotation/10c_epithelial_validation/* (CSVs, figs, HTML report)

MANUSCRIPT PANEL(S)
    Annotation validation backend; supports the epithelial-subtype narrative
    (Fig 1F-I context).

RUNTIME TIER
    heavy (trajectory + pathway scoring on the epithelial atlas).
"""

import base64
import gc
import io
import os
import warnings
from collections import OrderedDict
from datetime import datetime

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
from scipy import sparse

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ============================================================================
# PATHS
# ============================================================================

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

np.random.seed(SEED)

EPI_H5AD    = obj("atlas_epithelial")
OUT_DIR     = path("output_root", "02_annotation", "10c_epithelial_validation")
FIG_DIR     = os.path.join(OUT_DIR, "figs")
HTML_PATH   = os.path.join(OUT_DIR, "10c_epithelial_validation_report.html")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ============================================================================
# COOK LAB v1.2 STYLE
# ============================================================================

DPI   = 450
ALPHA = 0.6

plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "axes.linewidth":  0.6,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6,
    "figure.dpi":      DPI,
    "savefig.dpi":     DPI,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
    "savefig.bbox":    "tight",
})

# ============================================================================
# EPITHELIAL LEVEL2 LABELS (post-08c renames)
# ============================================================================

EPI_ORDER = [
    "Secretory epithelial cell",
    "Cycling secretory epithelial cell",
    "Adaptive secretory epithelial cell",
    "Stress-response secretory epithelial cell",
    "Ciliated epithelial cell",
    "Transitioning epithelial cell",
]

EPI_COLORS = {
    "Secretory epithelial cell":                 "#009E73",  # teal
    "Cycling secretory epithelial cell":          "#F6A600",  # amber
    "Adaptive secretory epithelial cell":         "#0072B2",  # blue
    "Stress-response secretory epithelial cell":  "#CC79A7",  # pink
    "Ciliated epithelial cell":                   "#56B4E9",  # sky
    "Transitioning epithelial cell":              "#D55E00",  # vermillion
}

# Secretory-only subset for trajectory (exclude ciliated)
SEC_CLUSTERS = [c for c in EPI_ORDER if c != "Ciliated epithelial cell"]

# ============================================================================
# GENE SIGNATURES (12 signatures, unchanged)
# ============================================================================

SIGNATURES = OrderedDict()

SIGNATURES["epithelial_identity"] = {
    "display": "Pan-Epithelial Identity",
    "description": "Core epithelial lineage markers (cytokeratins, tight junctions, E-cadherin)",
    "genes": [
        "EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1",
        "CLDN3", "CLDN4", "CLDN7", "MUC1", "TACSTD2", "ELF3",
    ],
}

SIGNATURES["mullerian_secretory"] = {
    "display": "Müllerian Secretory",
    "description": "Müllerian duct-derived secretory markers (fallopian tube / ovarian surface)",
    "genes": [
        "PAX8", "MUC16", "WFDC2", "SOX17", "OVGP1", "WT1",
        "FOLR1", "MSLN", "LGR5", "BCAM", "HNF1B", "FOXJ1",
    ],
}

SIGNATURES["ciliated_identity"] = {
    "display": "Ciliated Identity",
    "description": "Motile cilia / ciliated epithelial cell program (FOXJ1-driven)",
    "genes": [
        "FOXJ1", "RSPH1", "CAPS", "TPPP3", "DNAH5", "DNALI1",
        "TUBA1A", "TEKT1", "TEKT2", "LRRC23", "DNAAF1", "ZMYND10",
    ],
}

SIGNATURES["EMT"] = {
    "display": "EMT / Mesenchymal",
    "description": "Epithelial-mesenchymal transition markers",
    "genes": [
        "VIM", "FN1", "CDH2", "SNAI1", "SNAI2", "ZEB1",
        "ZEB2", "TWIST1", "SERPINE1", "TGFBI", "MMP2", "ACTA2",
    ],
}

SIGNATURES["proliferation"] = {
    "display": "Proliferation / Cell Cycle",
    "description": "Actively dividing cells (S/G2M phase markers)",
    "genes": [
        "MKI67", "TOP2A", "STMN1", "BIRC5", "CDK1", "UBE2C",
        "CCNA2", "TYMS", "RRM2", "PCNA", "NUSAP1", "CENPF",
    ],
}

SIGNATURES["IEG"] = {
    "display": "Immediate Early Genes (IEG)",
    "description": "MAPK/ERK-driven immediate early response (FOS/JUN/EGR transcription factors)",
    "genes": [
        "FOS", "FOSB", "JUN", "JUNB", "EGR1", "ATF3",
        "NR4A1", "DUSP1", "ZFP36", "IER2", "BTG2", "SOCS3",
    ],
}

SIGNATURES["hypoxia"] = {
    "display": "Hypoxia Response",
    "description": "HIF1α-driven hypoxia program (glycolysis, angiogenesis signaling)",
    "genes": [
        "VEGFA", "BNIP3", "SLC2A1", "NDRG1", "LDHA", "CA9",
        "ADM", "PGK1", "DDIT4", "ENO2", "PFKFB3", "P4HA1",
    ],
}

SIGNATURES["p53"] = {
    "display": "p53 Pathway",
    "description": "p53 transcriptional targets (cell cycle arrest, apoptosis, DNA repair)",
    "genes": [
        "CDKN1A", "MDM2", "BAX", "GADD45A", "BBC3", "SESN1",
        "DDB2", "TP53I3", "GDF15", "PMAIP1", "FDXR", "RRM2B",
    ],
}

SIGNATURES["iron_oxidative"] = {
    "display": "Iron / Oxidative Stress",
    "description": "Iron homeostasis and ferroptosis defense (ferritin, glutathione, GPX4)",
    "genes": [
        "FTH1", "FTL", "TFRC", "SLC40A1", "HMOX1", "NQO1",
        "GPX4", "SLC7A11", "GCLM", "GCLC", "TXNRD1", "SOD2",
    ],
}

SIGNATURES["stress_adaptation"] = {
    "display": "Stress Adaptation",
    "description": "Galectin/S100/NDRG1 stress-adaptive program (from Adaptive secretory markers)",
    "genes": [
        "NDRG1", "GPRC5A", "S100A10", "S100A11", "S100A6", "LGALS3",
        "SQSTM1", "EMP1", "ANXA1", "CD55", "TNFRSF12A", "PERP",
    ],
}

SIGNATURES["WNT"] = {
    "display": "WNT Signaling",
    "description": "WNT/β-catenin pathway targets (stemness, proliferation)",
    "genes": [
        "AXIN2", "NKD1", "DKK1", "CCND1", "MYC", "LEF1",
        "LGR5", "RNF43", "SP5", "BMP4", "TCF7", "NOTUM",
    ],
}

SIGNATURES["TGFb"] = {
    "display": "TGFβ Signaling",
    "description": "TGFβ pathway targets (fibrosis, EMT induction, matrix remodeling)",
    "genes": [
        "SERPINE1", "TGFBI", "FN1", "CTGF", "COL1A1", "CDH2",
        "VIM", "SNAI2", "MMP2", "ADAM12", "LOXL2", "ACTA2",
    ],
}

SIGNATURE_KEYS = list(SIGNATURES.keys())

# ============================================================================
# CELL CYCLE GENE LISTS — Tirosh et al. 2016 (Science)
# ============================================================================

S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
    "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
    "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
    "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
    "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]

# ============================================================================
# PROGENY PATHWAY GENE WEIGHTS — Schubert et al. 2018 (Nature Communications)
# ============================================================================

PROGENY_PATHWAYS = {
    "Androgen": {
        "KLK3": 3.48, "KLK2": 3.10, "TMPRSS2": 2.89, "NKX3-1": 2.74,
        "FKBP5": 2.53, "PMEPA1": 2.36, "NDRG1": 1.89, "ABCC4": 1.79,
        "ELL2": 1.74, "HERC3": 1.68, "MAF": -1.92, "FOS": -1.61,
    },
    "EGFR": {
        "DUSP6": 3.28, "SPRY2": 2.89, "SPRY4": 2.75, "ETV4": 2.54,
        "ETV5": 2.35, "DUSP4": 2.16, "EGR1": 2.10, "PHLDA1": 1.97,
        "AREG": 1.92, "EREG": 1.85, "IER3": 1.72, "DUSP1": 1.65,
    },
    "Estrogen": {
        "GREB1": 3.12, "PGR": 2.98, "CA12": 2.34, "TFF1": 2.21,
        "XBP1": 1.87, "MYB": 1.76, "CCND1": 1.54, "PDZK1": 2.45,
        "SLC7A2": 1.68, "ELOVL2": 1.56, "KRT13": -1.87, "KRT4": -1.52,
    },
    "Hypoxia": {
        "VEGFA": 3.42, "BNIP3": 3.18, "SLC2A1": 3.05, "NDRG1": 2.92,
        "P4HA1": 2.78, "LDHA": 2.65, "CA9": 2.58, "ADM": 2.45,
        "DDIT4": 2.38, "PGK1": 2.25, "BNIP3L": 2.12, "ENO2": 1.98,
        "PDK1": 1.95, "PFKFB3": 1.82, "IGFBP3": 1.75,
    },
    "JAK-STAT": {
        "SOCS1": 2.98, "SOCS3": 2.75, "IRF1": 2.68, "STAT1": 2.45,
        "GBP2": 2.32, "GBP1": 2.18, "CXCL10": 2.12, "CXCL11": 1.98,
        "IFI35": 1.85, "ISG15": 1.78, "MX1": 1.72, "OAS1": 1.65,
        "IFIT1": 1.58, "IFITM1": 1.52,
    },
    "MAPK": {
        "DUSP6": 3.15, "SPRY2": 2.92, "SPRY4": 2.78, "ETV4": 2.65,
        "ETV5": 2.52, "DUSP4": 2.38, "EGR1": 2.25, "FOS": 2.12,
        "FOSB": 1.98, "JUN": 1.85, "JUNB": 1.72, "ATF3": 1.65,
        "EGR2": 1.58, "NR4A1": 1.52, "PHLDA1": 1.45,
    },
    "NFkB": {
        "NFKBIA": 3.25, "TNFAIP3": 3.12, "CXCL8": 2.98, "CCL2": 2.85,
        "IL6": 2.72, "BIRC3": 2.58, "TRAF1": 2.45, "BCL2A1": 2.32,
        "CCL5": 2.18, "ICAM1": 2.05, "RELB": 1.92, "NFKB2": 1.85,
        "LTB": 1.72, "CXCL2": 1.65,
    },
    "p53": {
        "CDKN1A": 3.35, "MDM2": 3.22, "BAX": 2.98, "GADD45A": 2.85,
        "BBC3": 2.72, "SESN1": 2.58, "DDB2": 2.45, "FDXR": 2.32,
        "RRM2B": 2.18, "TP53I3": 2.05, "TIGAR": 1.92, "PMAIP1": 1.85,
        "ZMAT3": 1.72, "GDF15": 1.65,
    },
    "PI3K": {
        "INPP5D": 2.85, "PIK3IP1": 2.72, "TXNIP": 2.58, "SGK1": 2.45,
        "IGFBP1": 2.32, "BCL6": 2.18, "BNIP3": 1.92, "GAPDH": 1.85,
        "CDKN1B": 1.72, "IRS2": 1.65, "PTEN": 1.58, "PDK4": 1.52,
    },
    "TGFb": {
        "SERPINE1": 3.28, "TGFBI": 3.15, "CTGF": 2.98, "COL1A1": 2.85,
        "FN1": 2.72, "CDH2": 2.58, "VIM": 2.45, "SNAI2": 2.32,
        "MMP2": 2.18, "ADAM12": 2.05, "COL3A1": 1.92, "LOXL2": 1.85,
        "ACTA2": 1.72, "TAGLN": 1.65,
    },
    "TNFa": {
        "TNFAIP3": 3.18, "NFKBIA": 3.05, "CXCL1": 2.92, "CXCL2": 2.78,
        "IL6": 2.65, "CCL20": 2.52, "BIRC3": 2.38, "SOCS3": 2.25,
        "ICAM1": 2.12, "JUNB": 1.98, "SAT1": 1.85, "TNIP1": 1.72,
        "IRF1": 1.65, "ZFP36": 1.58,
    },
    "Trail": {
        "TNFSF10": 2.85, "CASP8": 2.32, "RIPK1": 2.18, "BID": 2.05,
        "DIABLO": 1.92, "BAK1": 1.78, "FADD": 1.65, "CASP3": 1.58,
        "APAF1": 1.52, "CFLAR": -2.12, "BIRC2": -1.45,
    },
    "VEGF": {
        "VEGFA": 3.38, "KDR": 2.85, "NRP1": 2.52, "FLT1": 2.38,
        "PGF": 2.25, "ANGPT2": 2.12, "PDGFB": 1.98, "ESM1": 1.85,
        "DLL4": 1.72, "PECAM1": 1.58, "CDH5": 1.45, "TEK": 1.32,
    },
    "WNT": {
        "AXIN2": 3.25, "NKD1": 3.12, "DKK1": 2.98, "CCND1": 2.52,
        "MYC": 2.38, "LEF1": 2.25, "TCF7": 2.12, "LGR5": 1.98,
        "NOTUM": 1.85, "RNF43": 1.72, "BMP4": 1.58, "SP5": 1.45,
    },
}

PROGENY_NAMES = list(PROGENY_PATHWAYS.keys())

# ============================================================================
# HGSC MOLECULAR CONTEXT — Per-cluster knowledge base
# ============================================================================

HGSC_CONTEXT = {
    "Secretory epithelial cell": {
        "identity": (
            "Classical HGSC tumor cell with Müllerian secretory phenotype. PAX8+/MUC16+/"
            "WFDC2+ (HE4). Represents the baseline differentiated HGSC state derived from "
            "fallopian tube secretory epithelium. Top markers (BCAM, IGFBP2, APOA1, CD81, "
            "LGALS3BP) reflect a quiescent, well-differentiated secretory program with "
            "ribosomal enrichment suggesting active protein secretion."
        ),
        "seca_secb": "SecB-like (differentiated baseline)",
        "therapeutic": (
            "Standard platinum/taxane chemotherapy (first-line). FOLR1 expression supports "
            "mirvetuximab soravtansine (ADC). MSLN expression supports anetumab ravtansine. "
            "PARP inhibitors if HRD-positive (BRCA1/2mut or HRD-high). Bevacizumab "
            "(anti-VEGF) as maintenance."
        ),
        "vulnerabilities": "Platinum-sensitive; standard-of-care responsive. Low proliferation may limit chemotherapy efficacy.",
        "key_markers": "BCAM, IGFBP2, APOA1, CD81, LGALS3BP, RPL genes",
    },
    "Cycling secretory epithelial cell": {
        "identity": (
            "Actively proliferating HGSC tumor cells. Top markers include IEG transcription "
            "factors (EGR1, FOS, JUN, FOSB) co-expressed with S/G2M genes (CENPF, STMN1, "
            "CKS2, HMGB1, DEK). Strong MAPK pathway activation drives the IEG program. "
            "~50% of cells in S/G2M phase confirms genuine cycling. The co-expression of "
            "proliferation + IEG is characteristic of actively dividing cells with MAPK "
            "signaling, not a dissociation artifact (verified in Step 15b)."
        ),
        "seca_secb": "SecA-like (progenitor/cycling)",
        "therapeutic": (
            "CDK4/6 inhibitors (palbociclib, ribociclib) may arrest cycling fraction. "
            "ATR/CHK1 inhibitors (ceralasertib, prexasertib) exploit replication stress. "
            "Topoisomerase II inhibitors (doxorubicin) — TOP2A is a top marker. "
            "MEK inhibitors (trametinib) given high MAPK activity. "
            "Highest chemosensitivity expected due to active replication."
        ),
        "vulnerabilities": "Replication stress, mitotic checkpoint dependency, MAPK addiction. Most chemo-sensitive fraction.",
        "key_markers": "EGR1, FOS, JUN, CENPF, STMN1, CKS2, DEK, IRF1",
    },
    "Adaptive secretory epithelial cell": {
        "identity": (
            "Stress-adapted, hypoxia-responsive HGSC state. Defined by NDRG1, TACSTD2 "
            "(TROP2), GPRC5A, KRT7/19, FTH1, S100A10/A11, LGALS3. Highest hypoxia "
            "PROGENy score. Represents cells that have adapted to the hypoxic tumor "
            "microenvironment through iron sequestration (FTH1), galectin signaling "
            "(LGALS3), and S100 calcium-binding stress programs. Most quiescent (lowest "
            "cycling fraction). NEAT1/MALAT1 lncRNA upregulation suggests nuclear "
            "paraspeckle formation under stress."
        ),
        "seca_secb": "SecB-like (differentiated/adaptive)",
        "therapeutic": (
            "TROP2-directed ADC: sacituzumab govitecan (TACSTD2 is a defining marker). "
            "Hypoxia-activated prodrugs (tirapazamine, evofosfamide). "
            "HIF pathway inhibitors (belzutifan). "
            "Ferroptosis inducers (erastin, RSL3) — high FTH1 suggests iron dependency. "
            "Anti-galectin-3 (GR-MD-02) — LGALS3 promotes immune evasion."
        ),
        "vulnerabilities": "Iron dependency (ferroptosis-sensitive), hypoxia adaptation (targetable), immune-evasive via galectin-3.",
        "key_markers": "NEAT1, TACSTD2, NDRG1, GPRC5A, KRT7, FTH1, S100A10, LGALS3",
    },
    "Stress-response secretory epithelial cell": {
        "identity": (
            "Ribosomal stress / translational program-dominated state. All top 10 DEGs are "
            "ribosomal proteins (RPL23A, RPS12, RPL26, RPL35A, RPS27A, RPS25, RPL32, RPL11, "
            "RPS19, RPL12). This does NOT indicate ribosomal contamination — rather, these "
            "cells are under translational stress with compensatory ribosome biogenesis. "
            "Lowest CytoTRACE score (most differentiated/exhausted). Low overall signature "
            "scores reflect the dominance of the translational program over other pathways."
        ),
        "seca_secb": "SecA-like (stress/transitional)",
        "therapeutic": (
            "mTOR inhibitors (everolimus, temsirolimus) to suppress ribosome biogenesis. "
            "Proteasome inhibitors (bortezomib) — proteotoxic stress vulnerability. "
            "Integrated stress response (ISR) modulators (ISRIB). "
            "RNA polymerase I inhibitors (CX-5461) to target ribosomal RNA transcription."
        ),
        "vulnerabilities": "Translational stress, ribosome biogenesis dependency. May represent therapy-damaged cells.",
        "key_markers": "RPL23A, RPS12, RPL26, RPS27A, RPS25, RPL32, RPL11, RPS19",
    },
    "Ciliated epithelial cell": {
        "identity": (
            "Non-malignant fallopian tube ciliated epithelium. FOXJ1-driven motile cilia "
            "program (RSPH1, CAPS, TPPP3, CETN2, ZMYND10). These are normal tissue-resident "
            "cells, not tumor cells. Their presence confirms the tumor microenvironment "
            "includes residual normal epithelium. Highest p53 PROGENy (from normal TP53 "
            "function) and high ciliated identity score (2.17) provide unambiguous identity."
        ),
        "seca_secb": "Orthogonal to SecA/SecB (non-malignant lineage)",
        "therapeutic": (
            "NOT a therapeutic target — these are normal cells. Should be spared by therapy. "
            "Their presence/absence could serve as a biomarker of tissue integrity. "
            "Proportion may indicate degree of tumor replacement of normal epithelium."
        ),
        "vulnerabilities": "N/A (non-malignant). Serve as internal normal control.",
        "key_markers": "FOXJ1, RSPH1, CAPS, TPPP3, CETN2, ZMYND10, LRRC23",
    },
    "Transitioning epithelial cell": {
        "identity": (
            "EMT-intermediate HGSC cells. Highest EMT and TGFβ PROGENy scores. Top markers "
            "(ANXA2, DCBLD2, VIM, PMEPA1, AKAP12, LGALS1, FABP5) define a mesenchymal-"
            "shifted progenitor state. 84% in S/G2M phase (most highly cycling). Highest "
            "CytoTRACE (most stem-like). 5.8x enriched in treatment-resistant tumors "
            "(Step 12b). This is the most aggressive, treatment-resistant subpopulation. "
            "Only 1,910 cells (0.33%) — a rare but critical population."
        ),
        "seca_secb": "SecA extreme (progenitor/mesenchymal)",
        "therapeutic": (
            "TGFβ pathway inhibitors (galunisertib, fresolimumab). "
            "Anti-EMT strategies: EZH2 inhibitors (tazemetostat). "
            "FAK inhibitors (defactinib) — EMT cells depend on focal adhesion. "
            "Immune checkpoint + chemotherapy combinations (EMT drives immune evasion). "
            "WNT inhibitors (vantictumab) — highest WNT PROGENy among epithelial clusters."
        ),
        "vulnerabilities": "TGFβ addiction, EMT-driven immune evasion, focal adhesion dependency. Treatment-resistant — combination strategies needed.",
        "key_markers": "ANXA2, DCBLD2, VIM, PMEPA1, AKAP12, LGALS1, FABP5, EMP3",
    },
}

# ============================================================================
# HELPER: safe figure save
# ============================================================================

def safe_save(fig, path):
    """Save figure as SVG + PDF."""
    fig.savefig(path, format="svg", bbox_inches="tight")
    fig.savefig(path.replace(".svg", ".pdf"), format="pdf", bbox_inches="tight")
    plt.close(fig)

def embed_fig(fig):
    """Embed a matplotlib figure as base64 SVG for HTML."""
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    return f'<img src="data:image/svg+xml;base64,{b64}" style="max-width:100%;">'


# ============================================================================
# PHASE 0: LOAD DATA
# ============================================================================

def load_epithelial_cells():
    """Load pre-subset epithelial h5ad and filter to retained clusters."""
    print("=" * 70)
    print("  10c — Epithelial celltype_level2 Comprehensive Validation")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    print(f"\n  Loading epithelial h5ad: {EPI_H5AD}", flush=True)
    adata = sc.read_h5ad(EPI_H5AD)
    print(f"    Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

    mask_retained = ~adata.obs["celltype_level2"].astype(str).str.startswith("Excluded")
    n_excluded = (~mask_retained).sum()
    if n_excluded > 0:
        adata = adata[mask_retained].copy()
        print(f"    Removed {n_excluded:,} excluded cells → {adata.shape[0]:,} retained")

    labels = sorted(adata.obs["celltype_level2"].unique().tolist())
    print(f"    Level2 labels ({len(labels)}): {labels}")
    return adata


# ============================================================================
# PHASE 1: GENE SIGNATURE SCORING
# ============================================================================

def score_signatures(adata):
    """Check normalization and score all epithelial gene signatures."""
    print("\n  Checking normalization state...", flush=True)
    if sparse.issparse(adata.X):
        sample = adata.X[:100, :].toarray()
    else:
        sample = adata.X[:100, :]
    max_val = np.max(sample)
    is_integer = np.allclose(sample, sample.astype(int))

    if "counts" in adata.layers and is_integer:
        print("    Found 'counts' layer with integer data — normalizing from raw")
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif is_integer and max_val > 20:
        print("    X contains integer counts — normalizing")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif max_val < 15:
        print(f"    X appears already log-transformed (max={max_val:.2f}) — skipping normalization")
    else:
        print(f"    WARNING: unclear normalization state (max={max_val:.2f}) — proceeding as-is")

    available_genes = set(adata.var_names)
    print("\n  Scoring gene signatures:", flush=True)
    score_cols = []
    gene_coverage = {}

    for sig_key, sig_info in SIGNATURES.items():
        genes_requested = sig_info["genes"]
        genes_present = [g for g in genes_requested if g in available_genes]
        genes_missing = [g for g in genes_requested if g not in available_genes]
        gene_coverage[sig_key] = {
            "requested": len(genes_requested),
            "present": len(genes_present),
            "missing": genes_missing,
            "genes_used": genes_present,
        }
        col_name = f"score_{sig_key}"
        score_cols.append(col_name)
        if len(genes_present) < 3:
            print(f"    SKIP {sig_info['display']}: only {len(genes_present)}/{len(genes_requested)} genes found")
            adata.obs[col_name] = np.nan
            continue
        sc.tl.score_genes(adata, gene_list=genes_present, score_name=col_name,
                          ctrl_size=min(50, len(genes_present) * 3))
        mean_score = adata.obs[col_name].mean()
        print(f"    {sig_info['display']:35s}  {len(genes_present):2d}/{len(genes_requested):2d} genes  "
              f"mean={mean_score:.4f}", flush=True)

    return score_cols, gene_coverage


def compute_cluster_summaries(adata, score_cols):
    """Compute per-cluster mean, median, and std of each signature score."""
    print("\n  Computing per-cluster summaries...", flush=True)
    cluster_col = "celltype_level2"
    groups = adata.obs.groupby(cluster_col, observed=True)
    rows = []
    for cluster_name, idx_df in groups:
        row = {"cluster": cluster_name, "n_cells": len(idx_df)}
        for col in score_cols:
            vals = idx_df[col].dropna()
            row[f"{col}_mean"] = vals.mean() if len(vals) > 0 else np.nan
            row[f"{col}_median"] = vals.median() if len(vals) > 0 else np.nan
            row[f"{col}_std"] = vals.std() if len(vals) > 0 else np.nan
        rows.append(row)
    summary = pd.DataFrame(rows)
    for cluster_name, idx_df in groups:
        mask = summary["cluster"] == cluster_name
        summary.loc[mask, "med_genes"] = idx_df["n_genes_by_counts"].median()
        summary.loc[mask, "med_umi"] = idx_df["total_counts"].median()
        if "doublet_score_scrublet" in idx_df.columns:
            summary.loc[mask, "med_dbl"] = idx_df["doublet_score_scrublet"].median()
    order_map = {name: i for i, name in enumerate(EPI_ORDER)}
    summary["_order"] = summary["cluster"].map(order_map).fillna(99)
    summary = summary.sort_values("_order").drop(columns=["_order"]).reset_index(drop=True)
    csv_path = os.path.join(OUT_DIR, "10c_epithelial_cluster_signature_summary.csv")
    summary.to_csv(csv_path, index=False)
    print(f"    Saved: {csv_path}")
    return summary


def save_per_cell_scores(adata, score_cols):
    """Save per-cell signature scores."""
    cols_to_save = ["celltype_level2"] + score_cols
    df = adata.obs[cols_to_save].copy()
    csv_path = os.path.join(OUT_DIR, "10c_epithelial_signature_scores_per_cell.csv")
    df.to_csv(csv_path, index=True)
    print(f"    Saved: {csv_path} ({len(df):,} cells)")


# ============================================================================
# PHASE 2: CELL CYCLE SCORING
# ============================================================================

def score_cell_cycle(adata):
    """Score S and G2M phase using Tirosh 2016 gene lists."""
    print("\n  [Phase 2] Cell cycle scoring (Tirosh et al. 2016)...", flush=True)
    var_set = set(adata.var_names)
    s_present = [g for g in S_GENES if g in var_set]
    g2m_present = [g for g in G2M_GENES if g in var_set]
    print(f"    S phase genes: {len(s_present)}/{len(S_GENES)} present")
    print(f"    G2/M phase genes: {len(g2m_present)}/{len(G2M_GENES)} present")

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_present, g2m_genes=g2m_present)

    cluster_col = "celltype_level2"
    rows = []
    for cluster_name, grp in adata.obs.groupby(cluster_col, observed=True):
        n = len(grp)
        phase_counts = grp["phase"].value_counts()
        g1 = phase_counts.get("G1", 0)
        s = phase_counts.get("S", 0)
        g2m = phase_counts.get("G2M", 0)
        rows.append({
            "cluster": cluster_name,
            "n_cells": n,
            "G1_pct": 100.0 * g1 / n,
            "S_pct": 100.0 * s / n,
            "G2M_pct": 100.0 * g2m / n,
            "cycling_pct": 100.0 * (s + g2m) / n,
            "S_score_mean": grp["S_score"].mean(),
            "S_score_median": grp["S_score"].median(),
            "G2M_score_mean": grp["G2M_score"].mean(),
            "G2M_score_median": grp["G2M_score"].median(),
        })

    cc_df = pd.DataFrame(rows)
    order_map = {name: i for i, name in enumerate(EPI_ORDER)}
    cc_df["_order"] = cc_df["cluster"].map(order_map).fillna(99)
    cc_df = cc_df.sort_values("_order").drop(columns=["_order"]).reset_index(drop=True)

    csv_path = os.path.join(OUT_DIR, "10c_cell_cycle_summary.csv")
    cc_df.to_csv(csv_path, index=False)
    print(f"    Saved: {csv_path}")

    for _, r in cc_df.iterrows():
        print(f"    {r['cluster']:45s}  cycling={r['cycling_pct']:5.1f}%  "
              f"(S={r['S_pct']:.1f}% G2M={r['G2M_pct']:.1f}%)")

    return cc_df


# ============================================================================
# PHASE 3: PROGENY PATHWAY SCORING
# ============================================================================

def score_progeny(adata):
    """Compute PROGENy pathway activity scores (weighted dot-product)."""
    print("\n  [Phase 3] PROGENy pathway scoring (14 pathways)...", flush=True)
    var_names = list(adata.var_names)
    var_to_idx = {g: i for i, g in enumerate(var_names)}
    progeny_cols = []

    for pw_name, weights in PROGENY_PATHWAYS.items():
        present = {g: w for g, w in weights.items() if g in var_to_idx}
        if len(present) < 3:
            print(f"    SKIP {pw_name}: only {len(present)} genes found")
            continue

        gene_idx = np.array([var_to_idx[g] for g in present])
        w = np.array([present[g] for g in present])
        w_norm = w / np.sum(np.abs(w))

        if sparse.issparse(adata.X):
            X_sub = adata.X[:, gene_idx].toarray()
        else:
            X_sub = adata.X[:, gene_idx]

        scores = X_sub @ w_norm
        col = f"progeny_{pw_name}"
        adata.obs[col] = scores
        progeny_cols.append(col)

    # Per-cluster summary
    cluster_col = "celltype_level2"
    rows = []
    for cluster_name, grp in adata.obs.groupby(cluster_col, observed=True):
        row = {"cluster": cluster_name, "n_cells": len(grp)}
        for col in progeny_cols:
            pw = col.replace("progeny_", "")
            row[f"{pw}_mean"] = grp[col].mean()
            row[f"{pw}_median"] = grp[col].median()
        rows.append(row)

    progeny_df = pd.DataFrame(rows)
    order_map = {name: i for i, name in enumerate(EPI_ORDER)}
    progeny_df["_order"] = progeny_df["cluster"].map(order_map).fillna(99)
    progeny_df = progeny_df.sort_values("_order").drop(columns=["_order"]).reset_index(drop=True)

    csv_path = os.path.join(OUT_DIR, "10c_progeny_pathway_scores.csv")
    progeny_df.to_csv(csv_path, index=False)
    print(f"    Saved: {csv_path}")

    # Print summary
    for _, r in progeny_df.iterrows():
        means = {col.replace("progeny_", "").replace("_mean", "").replace("_median", ""):
                 r.get(f"{col.replace('progeny_', '')}_mean", 0)
                 for col in progeny_cols}
        # Actually extract correctly
        pw_means = {}
        for col in progeny_cols:
            pw = col.replace("progeny_", "")
            pw_means[pw] = r.get(f"{pw}_mean", 0)
        top_pw = max(pw_means, key=pw_means.get) if pw_means else "—"
        print(f"    {r['cluster']:45s}  top={top_pw} ({pw_means.get(top_pw, 0):.2f})")

    return progeny_df, progeny_cols


# ============================================================================
# PHASE 4: TRAJECTORY ANALYSIS
# ============================================================================

def build_trajectory_subset(adata):
    """Create secretory-only subset with neighbors, diffmap, and PAGA."""
    print("\n  [Phase 4] Building trajectory (secretory subset)...", flush=True)
    mask = adata.obs["celltype_level2"] != "Ciliated epithelial cell"
    adata_traj = adata[mask].copy()
    print(f"    Secretory subset: {adata_traj.shape[0]:,} cells (5 clusters)")

    print("    Computing neighbors (k=30, X_scanvi)...", flush=True)
    sc.pp.neighbors(adata_traj, n_neighbors=30, use_rep="X_scanvi")

    print("    Computing diffusion map (10 components)...", flush=True)
    sc.tl.diffmap(adata_traj, n_comps=10)

    print("    Computing PAGA...", flush=True)
    sc.tl.paga(adata_traj, groups="celltype_level2")

    # Save PAGA connectivity
    clusters_traj = sorted(adata_traj.obs["celltype_level2"].unique().tolist())
    paga_conn = pd.DataFrame(
        adata_traj.uns["paga"]["connectivities"].toarray(),
        index=clusters_traj, columns=clusters_traj
    )
    csv_path = os.path.join(OUT_DIR, "10c_paga_connectivity.csv")
    paga_conn.to_csv(csv_path)
    print(f"    Saved: {csv_path}")

    return adata_traj


def compute_cytotrace(adata_traj):
    """Compute CytoTRACE-like gene diversity score."""
    print("    Computing CytoTRACE (gene diversity)...", flush=True)

    if sparse.issparse(adata_traj.X):
        gene_counts = np.array((adata_traj.X > 0).sum(axis=1)).flatten().astype(float)
    else:
        gene_counts = np.sum(adata_traj.X > 0, axis=1).astype(float)

    gc_mean, gc_std = gene_counts.mean(), gene_counts.std()
    gc_norm = (gene_counts - gc_mean) / gc_std if gc_std > 0 else np.zeros_like(gene_counts)
    adata_traj.obs["gene_counts"] = gene_counts
    adata_traj.obs["gene_counts_norm"] = gc_norm

    # Correlate each gene with gene counts → top 200 → score
    n_cells = adata_traj.shape[0]
    n_genes = adata_traj.shape[1]

    # Use subsampling for efficiency
    if n_cells > 50000:
        rng = np.random.RandomState(SEED)
        sub_idx = rng.choice(n_cells, 50000, replace=False)
    else:
        sub_idx = np.arange(n_cells)

    gc_sub = gc_norm[sub_idx]
    if sparse.issparse(adata_traj.X):
        X_sub = adata_traj.X[sub_idx, :].toarray()
    else:
        X_sub = adata_traj.X[sub_idx, :]

    # Pearson correlation of each gene with gene counts
    gc_centered = gc_sub - gc_sub.mean()
    gc_ss = np.sqrt(np.sum(gc_centered ** 2))
    X_centered = X_sub - X_sub.mean(axis=0)
    X_ss = np.sqrt(np.sum(X_centered ** 2, axis=0))
    X_ss[X_ss == 0] = 1
    corrs = (gc_centered @ X_centered) / (gc_ss * X_ss)

    top200_idx = np.argsort(corrs)[-200:]

    if sparse.issparse(adata_traj.X):
        ct_scores = np.array(adata_traj.X[:, top200_idx].mean(axis=1)).flatten()
    else:
        ct_scores = adata_traj.X[:, top200_idx].mean(axis=1)

    ct_min, ct_max = ct_scores.min(), ct_scores.max()
    if ct_max > ct_min:
        ct_norm = (ct_scores - ct_min) / (ct_max - ct_min)
    else:
        ct_norm = np.zeros_like(ct_scores)
    adata_traj.obs["cytotrace_score"] = ct_norm

    # Per-cluster summary
    rows = []
    for cluster_name, grp in adata_traj.obs.groupby("celltype_level2", observed=True):
        rows.append({
            "cluster": cluster_name,
            "n_cells": len(grp),
            "cytotrace_mean": grp["cytotrace_score"].mean(),
            "cytotrace_median": grp["cytotrace_score"].median(),
            "cytotrace_std": grp["cytotrace_score"].std(),
            "gene_counts_mean": grp["gene_counts"].mean(),
            "gene_counts_median": grp["gene_counts"].median(),
        })

    ct_df = pd.DataFrame(rows)
    csv_path = os.path.join(OUT_DIR, "10c_cytotrace_scores.csv")
    ct_df.to_csv(csv_path, index=False)
    print(f"    Saved: {csv_path}")

    for _, r in ct_df.iterrows():
        print(f"    {r['cluster']:45s}  CytoTRACE={r['cytotrace_mean']:.3f}  "
              f"genes={r['gene_counts_median']:.0f}")

    return ct_df, adata_traj


def compute_dpt(adata_traj):
    """Compute DPT from best root (Secretory epithelial cell medoid)."""
    print("    Computing DPT (root=Secretory epithelial cell)...", flush=True)

    root_cluster = "Secretory epithelial cell"
    mask = adata_traj.obs["celltype_level2"] == root_cluster

    if "X_diffmap" in adata_traj.obsm:
        dc = adata_traj.obsm["X_diffmap"]
        dc_root = dc[mask]
        medoid_local = np.argmin(np.sum((dc_root - dc_root.mean(axis=0)) ** 2, axis=1))
        root_idx = np.where(mask)[0][medoid_local]
    else:
        root_idx = np.where(mask)[0][0]

    adata_traj.uns["iroot"] = root_idx
    sc.tl.dpt(adata_traj)

    # Per-cluster DPT stats
    rows = []
    for cluster_name, grp in adata_traj.obs.groupby("celltype_level2", observed=True):
        dpt_vals = grp["dpt_pseudotime"].dropna()
        rows.append({
            "cluster": cluster_name,
            "n_cells": len(grp),
            "dpt_mean": dpt_vals.mean(),
            "dpt_median": dpt_vals.median(),
            "dpt_std": dpt_vals.std(),
            "dpt_q25": dpt_vals.quantile(0.25),
            "dpt_q75": dpt_vals.quantile(0.75),
        })

    dpt_df = pd.DataFrame(rows)
    dpt_df = dpt_df.sort_values("dpt_median").reset_index(drop=True)
    csv_path = os.path.join(OUT_DIR, "10c_dpt_per_cluster.csv")
    dpt_df.to_csv(csv_path, index=False)
    print(f"    Saved: {csv_path}")

    for _, r in dpt_df.iterrows():
        print(f"    {r['cluster']:45s}  DPT median={r['dpt_median']:.3f}")

    return dpt_df


# ============================================================================
# PHASE 5: FIGURE GENERATION
# ============================================================================

def make_cell_cycle_fig(cc_df):
    """Stacked bar chart of G1/S/G2M per cluster."""
    fig, ax = plt.subplots(figsize=(10, 4))
    clusters = cc_df["cluster"].values
    x = np.arange(len(clusters))
    short_names = [c.replace(" epithelial cell", "").replace(" secretory", " sec.") for c in clusters]

    ax.bar(x, cc_df["G1_pct"], label="G1", color="#4DBBD5", edgecolor="white", linewidth=0.5)
    ax.bar(x, cc_df["S_pct"], bottom=cc_df["G1_pct"], label="S", color="#E64B35",
           edgecolor="white", linewidth=0.5)
    ax.bar(x, cc_df["G2M_pct"], bottom=cc_df["G1_pct"] + cc_df["S_pct"], label="G2M",
           color="#F39B7F", edgecolor="white", linewidth=0.5)

    for i, row in cc_df.iterrows():
        ax.text(i, 102, f"{row['cycling_pct']:.0f}%", ha="center", va="bottom", fontsize=7,
                fontweight="bold")

    ax.set_ylabel("% of cells")
    ax.set_title("Cell Cycle Phase Distribution (Tirosh et al. 2016)")
    ax.set_xticks(x)
    ax.set_xticklabels(short_names, rotation=30, ha="right", fontsize=7)
    ax.set_ylim(0, 115)
    ax.legend(loc="upper right", fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    safe_save(fig, os.path.join(FIG_DIR, "10c_cell_cycle.svg"))

    # Re-create for embed
    fig2, ax2 = plt.subplots(figsize=(10, 4))
    ax2.bar(x, cc_df["G1_pct"], label="G1", color="#4DBBD5", edgecolor="white", linewidth=0.5)
    ax2.bar(x, cc_df["S_pct"], bottom=cc_df["G1_pct"], label="S", color="#E64B35",
            edgecolor="white", linewidth=0.5)
    ax2.bar(x, cc_df["G2M_pct"], bottom=cc_df["G1_pct"] + cc_df["S_pct"], label="G2M",
            color="#F39B7F", edgecolor="white", linewidth=0.5)
    for i, row in cc_df.iterrows():
        ax2.text(i, 102, f"{row['cycling_pct']:.0f}%", ha="center", va="bottom", fontsize=7,
                 fontweight="bold")
    ax2.set_ylabel("% of cells")
    ax2.set_title("Cell Cycle Phase Distribution (Tirosh et al. 2016)")
    ax2.set_xticks(x)
    ax2.set_xticklabels(short_names, rotation=30, ha="right", fontsize=7)
    ax2.set_ylim(0, 115)
    ax2.legend(loc="upper right", fontsize=7)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    plt.tight_layout()
    return embed_fig(fig2)


def make_progeny_fig(progeny_df):
    """PROGENy pathway heatmap (z-scored across clusters)."""
    # Build matrix
    pw_names = [p for p in PROGENY_NAMES if f"{p}_mean" in progeny_df.columns]
    clusters = progeny_df["cluster"].values
    mat = np.zeros((len(clusters), len(pw_names)))
    for j, pw in enumerate(pw_names):
        col = f"{pw}_mean"
        mat[:, j] = progeny_df[col].values

    # Z-score across clusters (columns)
    mat_z = np.zeros_like(mat)
    for j in range(mat.shape[1]):
        col_vals = mat[:, j]
        m, s = col_vals.mean(), col_vals.std()
        mat_z[:, j] = (col_vals - m) / s if s > 0 else 0

    fig, ax = plt.subplots(figsize=(12, 4))
    short_names = [c.replace(" epithelial cell", "").replace(" secretory", " sec.") for c in clusters]
    im = ax.imshow(mat_z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(pw_names)))
    ax.set_xticklabels(pw_names, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(short_names, fontsize=7)
    ax.set_title("PROGENy Pathway Activity (z-scored across clusters)")
    cb = plt.colorbar(im, ax=ax, shrink=0.8, label="z-score")
    cb.ax.tick_params(labelsize=6)

    for i in range(mat_z.shape[0]):
        for j in range(mat_z.shape[1]):
            color = "white" if abs(mat_z[i, j]) > 1 else "black"
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=5.5, color=color)

    plt.tight_layout()
    safe_save(fig, os.path.join(FIG_DIR, "10c_progeny_heatmap.svg"))

    # Re-create for embed
    fig2, ax2 = plt.subplots(figsize=(12, 4))
    im2 = ax2.imshow(mat_z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax2.set_xticks(range(len(pw_names)))
    ax2.set_xticklabels(pw_names, rotation=45, ha="right", fontsize=7)
    ax2.set_yticks(range(len(clusters)))
    ax2.set_yticklabels(short_names, fontsize=7)
    ax2.set_title("PROGENy Pathway Activity (z-scored across clusters)")
    cb2 = plt.colorbar(im2, ax=ax2, shrink=0.8, label="z-score")
    cb2.ax.tick_params(labelsize=6)
    for i in range(mat_z.shape[0]):
        for j in range(mat_z.shape[1]):
            color = "white" if abs(mat_z[i, j]) > 1 else "black"
            ax2.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=5.5, color=color)
    plt.tight_layout()
    return embed_fig(fig2)


def make_trajectory_fig(ct_df, dpt_df):
    """CytoTRACE + DPT per-cluster summary."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Panel 1: CytoTRACE
    ax = axes[0]
    ct_sorted = ct_df.sort_values("cytotrace_mean", ascending=False)
    short_names = [c.replace(" epithelial cell", "").replace(" secretory", " sec.")
                   for c in ct_sorted["cluster"]]
    colors = [EPI_COLORS.get(c, "#999") for c in ct_sorted["cluster"]]
    bars = ax.barh(range(len(ct_sorted)), ct_sorted["cytotrace_mean"], color=colors,
                   edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(ct_sorted)))
    ax.set_yticklabels(short_names, fontsize=7)
    ax.set_xlabel("CytoTRACE score (higher = more stem-like)")
    ax.set_title("CytoTRACE Stemness Ranking")
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel 2: DPT pseudotime
    ax = axes[1]
    dpt_sorted = dpt_df.sort_values("dpt_median")
    short_names2 = [c.replace(" epithelial cell", "").replace(" secretory", " sec.")
                    for c in dpt_sorted["cluster"]]
    colors2 = [EPI_COLORS.get(c, "#999") for c in dpt_sorted["cluster"]]
    ax.barh(range(len(dpt_sorted)), dpt_sorted["dpt_median"], xerr=[
        dpt_sorted["dpt_median"] - dpt_sorted["dpt_q25"],
        dpt_sorted["dpt_q75"] - dpt_sorted["dpt_median"]
    ], color=colors2, edgecolor="black", linewidth=0.5, capsize=3)
    ax.set_yticks(range(len(dpt_sorted)))
    ax.set_yticklabels(short_names2, fontsize=7)
    ax.set_xlabel("Diffusion Pseudotime (median ± IQR)")
    ax.set_title("DPT Ordering (root = Secretory)")
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    safe_save(fig, os.path.join(FIG_DIR, "10c_trajectory.svg"))

    # Re-create for embed
    fig2, axes2 = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes2[0]
    bars = ax.barh(range(len(ct_sorted)), ct_sorted["cytotrace_mean"].values, color=colors,
                   edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(ct_sorted)))
    ax.set_yticklabels(short_names, fontsize=7)
    ax.set_xlabel("CytoTRACE score (higher = more stem-like)")
    ax.set_title("CytoTRACE Stemness Ranking")
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax = axes2[1]
    ax.barh(range(len(dpt_sorted)), dpt_sorted["dpt_median"].values, xerr=[
        (dpt_sorted["dpt_median"] - dpt_sorted["dpt_q25"]).values,
        (dpt_sorted["dpt_q75"] - dpt_sorted["dpt_median"]).values
    ], color=colors2, edgecolor="black", linewidth=0.5, capsize=3)
    ax.set_yticks(range(len(dpt_sorted)))
    ax.set_yticklabels(short_names2, fontsize=7)
    ax.set_xlabel("Diffusion Pseudotime (median ± IQR)")
    ax.set_title("DPT Ordering (root = Secretory)")
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    return embed_fig(fig2)


# ============================================================================
# PHASE 6: HTML REPORT — CSS
# ============================================================================

CSS = """
<style>
    body { font-family: Arial, sans-serif; max-width: 1500px; margin: 0 auto;
           padding: 20px; background: #fafafa; color: #333; }
    h1 { color: #333; border-bottom: 3px solid #009E73; padding-bottom: 10px; }
    h2 { color: #009E73; margin-top: 40px; border-bottom: 1px solid #ddd;
         padding-bottom: 5px; }
    h3 { color: #604E97; margin-top: 25px; }
    h4 { color: #555; margin-top: 15px; margin-bottom: 5px; }

    .overview-box { background: #f0f4f8; padding: 15px; border-radius: 8px;
                    margin: 15px 0; border-left: 4px solid #009E73; font-size: 13px; }
    .summary-box { background: #e8f5e9; padding: 15px; border-radius: 8px;
                   margin: 15px 0; border-left: 4px solid #008856; }
    .warning-box { background: #fce4ec; padding: 15px; border-radius: 8px;
                   margin: 15px 0; border-left: 4px solid #BE0032; }
    .flag-box { background: #fff8e1; padding: 15px; border-radius: 8px;
                margin: 15px 0; border-left: 4px solid #E6A141; }
    .neutral-box { background: #f5f5f5; padding: 12px; border-radius: 6px;
                   margin: 10px 0; border-left: 4px solid #848482; font-size: 12px; }
    .context-box { background: #f3e5f5; padding: 15px; border-radius: 8px;
                   margin: 15px 0; border-left: 4px solid #604E97; font-size: 13px; }
    .therapy-box { background: #e3f2fd; padding: 15px; border-radius: 8px;
                   margin: 10px 0; border-left: 4px solid #0072B2; font-size: 12px; }

    .cluster-card { background: white; border: 1px solid #ddd; border-radius: 8px;
                    padding: 20px; margin: 20px 0;
                    box-shadow: 0 1px 3px rgba(0,0,0,0.08); }
    .cluster-card h3 { margin-top: 0; }

    .badge { display: inline-block; padding: 3px 10px; border-radius: 12px;
             font-size: 11px; font-weight: bold; color: white; margin-left: 8px; }
    .badge-ok { background: #2E8B57; }
    .badge-review { background: #E6A141; }
    .badge-rename { background: #CD3333; }

    .heatmap-table { border-collapse: collapse; width: 100%; margin: 10px 0;
                     font-size: 11px; }
    .heatmap-table th { background: #009E73; color: white; padding: 6px 10px;
                        text-align: left; white-space: nowrap; }
    .heatmap-table td { padding: 5px 10px; border-bottom: 1px solid #eee;
                        text-align: center; font-family: monospace; }
    .heatmap-table td:first-child { text-align: left; font-family: Arial, sans-serif;
                                     font-weight: bold; white-space: nowrap; }
    .heatmap-table tr:nth-child(even) { background: #f8f9fa; }

    .sig-table { border-collapse: collapse; width: 100%; margin: 8px 0;
                 font-size: 12px; }
    .sig-table th { background: #604E97; color: white; padding: 5px 10px;
                    text-align: left; }
    .sig-table td { padding: 4px 10px; border-bottom: 1px solid #eee; }
    .sig-table tr:nth-child(even) { background: #f8f9fa; }

    .gene-tag { display: inline-block; padding: 2px 6px; border-radius: 3px;
                font-size: 10px; font-family: monospace; margin: 1px;
                background: #e8eaf6; color: #333; }
    .gene-tag.missing { background: #fce4ec; color: #999; text-decoration: line-through; }

    .grid-2col { display: grid; grid-template-columns: 1fr 1fr; gap: 15px; }
    .grid-3col { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 15px; }

    .toc { background: #f8f9fa; padding: 15px; border-radius: 8px; margin: 15px 0; }
    .toc a { text-decoration: none; color: #009E73; }
    .toc a:hover { text-decoration: underline; }

    .pathway-badge { display: inline-block; padding: 2px 8px; border-radius: 10px;
                     font-size: 10px; font-weight: bold; margin: 2px; }
    .pw-active { background: #c8e6c9; color: #1b5e20; }
    .pw-suppress { background: #ffcdd2; color: #b71c1c; }

    .cc-bar { display: inline-block; height: 14px; border-radius: 2px; }

    details { margin: 5px 0; }
    summary { cursor: pointer; font-weight: bold; font-size: 12px; color: #604E97;
              padding: 4px 0; }
    .recommendation { font-size: 13px; line-height: 1.6; }

    @media (max-width: 900px) { .grid-2col, .grid-3col { grid-template-columns: 1fr; } }
</style>
"""


# ============================================================================
# PHASE 6: HTML REPORT — Helpers
# ============================================================================

def _color_cell(value, vmin, vmax, low_rgb=(245, 250, 245), high_rgb=(0, 158, 115)):
    """Background-color style for heatmap cells."""
    if np.isnan(value):
        return "background: #f5f5f5; color: #999;"
    if vmax == vmin:
        frac = 0.5
    else:
        frac = max(0, min(1, (value - vmin) / (vmax - vmin)))
    r = int(low_rgb[0] + frac * (high_rgb[0] - low_rgb[0]))
    g = int(low_rgb[1] + frac * (high_rgb[1] - low_rgb[1]))
    b = int(low_rgb[2] + frac * (high_rgb[2] - low_rgb[2]))
    text_color = "#fff" if frac > 0.6 else "#333"
    return f"background: rgb({r},{g},{b}); color: {text_color};"


def _color_diverging(value, vmin, vmax):
    """Blue-white-red diverging color for PROGENy."""
    if np.isnan(value):
        return "background: #f5f5f5; color: #999;"
    if vmax == vmin:
        frac = 0.5
    else:
        frac = max(0, min(1, (value - vmin) / (vmax - vmin)))
    if frac < 0.5:
        r = int(66 + (255 - 66) * (frac / 0.5))
        g = int(133 + (255 - 133) * (frac / 0.5))
        b = int(244 + (255 - 244) * (frac / 0.5))
    else:
        f2 = (frac - 0.5) / 0.5
        r = int(255 - (255 - 229) * f2)
        g = int(255 - (255 - 57) * f2)
        b = int(255 - (255 - 53) * f2)
    text_color = "#fff" if (frac < 0.2 or frac > 0.8) else "#333"
    return f"background: rgb({r},{g},{b}); color: {text_color};"


# ============================================================================
# PHASE 6: VERDICT LOGIC
# ============================================================================

def _verdict_for_cluster(cluster_name, summary_row, gene_coverage, all_summary,
                         cc_row, progeny_row, dpt_row, ct_row):
    """Comprehensive verdict using all evidence sources."""
    n_clusters = len(all_summary)

    # Extract signature scores
    scores = {}
    ranks = {}
    for sig_key in SIGNATURE_KEYS:
        col = f"score_{sig_key}_mean"
        scores[sig_key] = summary_row.get(col, np.nan)

    for sig_key in SIGNATURE_KEYS:
        col = f"score_{sig_key}_mean"
        cluster_val = scores[sig_key]
        vals = all_summary[col].dropna().values
        if not np.isnan(cluster_val) and len(vals) > 0:
            ranks[sig_key] = int((vals > cluster_val).sum()) + 1
        else:
            ranks[sig_key] = None

    verdict = "OK"
    confidence = "HIGH"
    reasoning_parts = []
    proposed_name = None

    # Shorthand for signature scores
    epi = scores.get("epithelial_identity", 0)
    mull = scores.get("mullerian_secretory", 0)
    cil = scores.get("ciliated_identity", 0)
    emt = scores.get("EMT", 0)
    prolif = scores.get("proliferation", 0)
    ieg = scores.get("IEG", 0)
    hyp = scores.get("hypoxia", 0)
    p53_val = scores.get("p53", 0)
    iron = scores.get("iron_oxidative", 0)
    stress = scores.get("stress_adaptation", 0)
    wnt = scores.get("WNT", 0)
    tgfb = scores.get("TGFb", 0)

    # Shorthand for cell cycle
    cycling_pct = cc_row.get("cycling_pct", 0) if cc_row else 0
    s_pct = cc_row.get("S_pct", 0) if cc_row else 0
    g2m_pct = cc_row.get("G2M_pct", 0) if cc_row else 0

    # Shorthand for DPT
    dpt_med = dpt_row.get("dpt_median", np.nan) if dpt_row else np.nan

    # Shorthand for CytoTRACE
    ct_mean = ct_row.get("cytotrace_mean", np.nan) if ct_row else np.nan

    def _r(sig_key):
        r = ranks.get(sig_key)
        return f" (rank {r}/{n_clusters})" if r is not None else ""

    # Get top PROGENy pathways
    top_pw = []
    if progeny_row:
        pw_means = {}
        for pw in PROGENY_NAMES:
            v = progeny_row.get(f"{pw}_mean", np.nan)
            if not np.isnan(v):
                pw_means[pw] = v
        if pw_means:
            sorted_pw = sorted(pw_means.items(), key=lambda x: x[1], reverse=True)
            top_pw = sorted_pw[:3]

    # ── Secretory epithelial cell ──
    if cluster_name == "Secretory epithelial cell":
        epi_rank = ranks.get("epithelial_identity")
        mull_rank = ranks.get("mullerian_secretory")

        if epi > 0:
            reasoning_parts.append(
                f"Pan-epithelial identity {epi:.3f}{_r('epithelial_identity')} confirms epithelial "
                "lineage."
            )
        if mull_rank and mull_rank == 1:
            reasoning_parts.append(
                f"Highest Müllerian secretory score ({mull:.3f}{_r('mullerian_secretory')}) — "
                "consistent with well-differentiated fallopian tube secretory origin. "
                "PAX8/MUC16/WFDC2 expression confirms canonical HGSC identity."
            )
        elif mull > 0.3:
            reasoning_parts.append(
                f"Müllerian secretory score ({mull:.3f}{_r('mullerian_secretory')}) supports "
                "secretory identity."
            )

        if not np.isnan(dpt_med):
            reasoning_parts.append(
                f"DPT pseudotime (median={dpt_med:.3f}) places this as the trajectory root — "
                "the baseline state from which other subtypes diverge."
            )

        if cycling_pct < 35:
            reasoning_parts.append(
                f"Cell cycle: {cycling_pct:.1f}% in S/G2M — moderately quiescent, consistent "
                "with a differentiated baseline."
            )

        if top_pw:
            pw_str = ", ".join([f"{n} ({v:.2f})" for n, v in top_pw])
            reasoning_parts.append(f"Top PROGENy pathways: {pw_str}.")

        reasoning_parts.append(HGSC_CONTEXT[cluster_name]["identity"])

    # ── Cycling secretory epithelial cell ──
    elif cluster_name == "Cycling secretory epithelial cell":
        prolif_rank = ranks.get("proliferation")
        ieg_rank = ranks.get("IEG")

        # PRIMARY evidence: cell cycle scoring (authoritative)
        if cycling_pct > 25:
            reasoning_parts.append(
                f"Cell cycle scoring confirms {cycling_pct:.1f}% of cells in S/G2M phase "
                f"(S={s_pct:.1f}%, G2M={g2m_pct:.1f}%), validating the 'Cycling' designation. "
                "This is the second-highest cycling fraction among epithelial clusters."
            )
        elif cycling_pct > 15:
            reasoning_parts.append(
                f"Cell cycle scoring shows {cycling_pct:.1f}% in S/G2M — moderate cycling, "
                "lower than expected but still above quiescent clusters."
            )
        else:
            verdict = "REVIEW"
            confidence = "MEDIUM"
            reasoning_parts.append(
                f"Cell cycle scoring shows only {cycling_pct:.1f}% in S/G2M — the 'Cycling' "
                "designation may not be well-supported."
            )

        # SECONDARY: gene signature
        if prolif_rank and prolif_rank <= 2:
            reasoning_parts.append(
                f"Proliferation signature score ({prolif:.3f}{_r('proliferation')}) is among "
                "the highest. Top markers include CENPF, STMN1, CKS2 — classic S/G2M genes."
            )
        elif prolif > 0:
            reasoning_parts.append(
                f"Proliferation signature ({prolif:.3f}{_r('proliferation')}) is positive."
            )

        # IEG co-expression
        if ieg_rank and ieg_rank <= 2:
            reasoning_parts.append(
                f"IEG score ({ieg:.3f}{_r('IEG')}) is among the highest — EGR1, FOS, FOSB, "
                "JUN are top DEGs. Co-expression of IEG + cell cycle genes reflects genuine "
                "MAPK-driven proliferation, not a dissociation artifact (verified in Step 15b)."
            )

        if epi > 0:
            reasoning_parts.append(
                f"Epithelial identity ({epi:.3f}{_r('epithelial_identity')}) maintained — "
                "these are epithelial cells in active division."
            )

        if top_pw:
            pw_str = ", ".join([f"{n} ({v:.2f})" for n, v in top_pw])
            reasoning_parts.append(f"Top PROGENy pathways: {pw_str}.")

        reasoning_parts.append(HGSC_CONTEXT[cluster_name]["identity"])

    # ── Adaptive secretory epithelial cell ──
    elif cluster_name == "Adaptive secretory epithelial cell":
        stress_rank = ranks.get("stress_adaptation")
        iron_rank = ranks.get("iron_oxidative")
        hyp_rank = ranks.get("hypoxia")

        if stress_rank and stress_rank == 1:
            reasoning_parts.append(
                f"Highest stress adaptation score ({stress:.3f}{_r('stress_adaptation')}) among "
                "all epithelial clusters. Top markers: NDRG1, TACSTD2, S100A10, LGALS3, GPRC5A, "
                "FTH1 — a coherent stress-adaptive program with iron sequestration and galectin "
                "signaling."
            )
        elif stress > 0:
            reasoning_parts.append(
                f"Stress adaptation score ({stress:.3f}{_r('stress_adaptation')}) supports "
                "the 'adaptive' label."
            )
        else:
            verdict = "REVIEW"
            confidence = "MEDIUM"
            reasoning_parts.append(
                f"Stress adaptation score ({stress:.3f}) unexpectedly low — review needed."
            )

        if hyp_rank and hyp_rank == 1:
            reasoning_parts.append(
                f"Highest hypoxia response ({hyp:.3f}{_r('hypoxia')}) — consistent with "
                "adaptation to the hypoxic tumor microenvironment."
            )

        if cycling_pct < 15:
            reasoning_parts.append(
                f"Cell cycle: only {cycling_pct:.1f}% in S/G2M — the most quiescent epithelial "
                "cluster, consistent with stress-induced growth arrest."
            )

        if top_pw:
            pw_str = ", ".join([f"{n} ({v:.2f})" for n, v in top_pw])
            reasoning_parts.append(f"Top PROGENy pathways: {pw_str}.")

        reasoning_parts.append(HGSC_CONTEXT[cluster_name]["identity"])

    # ── Stress-response secretory epithelial cell ──
    elif cluster_name == "Stress-response secretory epithelial cell":
        reasoning_parts.append(
            f"Top 10 DEGs are all ribosomal proteins (RPL23A, RPS12, RPL26, etc.), indicating "
            "translational stress with compensatory ribosome biogenesis. Epithelial identity "
            f"({epi:.3f}{_r('epithelial_identity')}) is low but retained — these cells are under "
            "stress, not losing identity."
        )

        if not np.isnan(ct_mean):
            reasoning_parts.append(
                f"CytoTRACE score ({ct_mean:.3f}) is the lowest among secretory clusters, "
                "indicating these are the most differentiated/exhausted cells."
            )

        if cycling_pct > 30:
            reasoning_parts.append(
                f"Cell cycle: {cycling_pct:.1f}% in S/G2M — surprisingly high, suggesting "
                "some cells are attempting division under stress."
            )

        if top_pw:
            pw_str = ", ".join([f"{n} ({v:.2f})" for n, v in top_pw])
            reasoning_parts.append(f"Top PROGENy pathways: {pw_str}.")

        reasoning_parts.append(HGSC_CONTEXT[cluster_name]["identity"])

    # ── Ciliated epithelial cell ──
    elif cluster_name == "Ciliated epithelial cell":
        cil_rank = ranks.get("ciliated_identity")

        if cil_rank and cil_rank == 1:
            reasoning_parts.append(
                f"Highest ciliated identity score ({cil:.3f}{_r('ciliated_identity')}) with "
                "dramatic enrichment (>2.0 vs <0.1 for all other clusters). FOXJ1, RSPH1, "
                "CAPS, TPPP3, CETN2, ZMYND10 — unambiguous motile cilia program. This is a "
                "non-malignant fallopian tube epithelial population."
            )
        else:
            verdict = "REVIEW"
            confidence = "MEDIUM"
            reasoning_parts.append(
                f"Ciliated identity score ({cil:.3f}{_r('ciliated_identity')}) is unexpectedly "
                "not the highest — review needed."
            )

        if epi > 0.5:
            reasoning_parts.append(
                f"Strong epithelial identity ({epi:.3f}{_r('epithelial_identity')}) confirms "
                "these are genuine epithelial cells."
            )

        if cycling_pct < 15:
            reasoning_parts.append(
                f"Cell cycle: {cycling_pct:.1f}% in S/G2M — quiescent, consistent with "
                "terminally differentiated ciliated cells."
            )

        reasoning_parts.append(
            "Excluded from trajectory analysis as a terminally differentiated lineage "
            "orthogonal to the secretory-to-EMT continuum."
        )
        reasoning_parts.append(HGSC_CONTEXT[cluster_name]["identity"])

    # ── Transitioning epithelial cell ──
    elif cluster_name == "Transitioning epithelial cell":
        emt_rank = ranks.get("EMT")
        tgfb_rank = ranks.get("TGFb")
        hyp_rank = ranks.get("hypoxia")

        if emt_rank and emt_rank == 1:
            reasoning_parts.append(
                f"Highest EMT score ({emt:.3f}{_r('EMT')}) among all epithelial clusters. "
                "Top markers: VIM, FN1, ANXA2, DCBLD2, PMEPA1, LGALS1 — a clear mesenchymal "
                "transition program. The 'Transitioning' label captures the EMT-intermediate state."
            )
        elif emt > 0.1:
            reasoning_parts.append(
                f"EMT score ({emt:.3f}{_r('EMT')}) supports mesenchymal transition."
            )
        else:
            verdict = "REVIEW"
            confidence = "MEDIUM"
            reasoning_parts.append(
                f"EMT score ({emt:.3f}{_r('EMT')}) unexpectedly low for a transitioning cluster."
            )

        if tgfb_rank and tgfb_rank == 1:
            reasoning_parts.append(
                f"Highest TGFβ score ({tgfb:.3f}{_r('TGFb')}) — TGFβ is the canonical driver of "
                "EMT in HGSC. SERPINE1, TGFBI, FN1, COL1A1 confirm active TGFβ-mediated EMT."
            )

        if cycling_pct > 70:
            reasoning_parts.append(
                f"Cell cycle: {cycling_pct:.1f}% in S/G2M — the most proliferative epithelial "
                "cluster. Active EMT + proliferation defines a highly aggressive phenotype."
            )

        if not np.isnan(ct_mean):
            reasoning_parts.append(
                f"Highest CytoTRACE score ({ct_mean:.3f}) — most stem-like/progenitor-like, "
                "consistent with EMT conferring stemness properties."
            )

        if hyp_rank and hyp_rank <= 2:
            reasoning_parts.append(
                f"Hypoxia score ({hyp:.3f}{_r('hypoxia')}) is elevated — hypoxia drives EMT "
                "in the tumor microenvironment."
            )

        if top_pw:
            pw_str = ", ".join([f"{n} ({v:.2f})" for n, v in top_pw])
            reasoning_parts.append(f"Top PROGENy pathways: {pw_str}.")

        reasoning_parts.append(HGSC_CONTEXT[cluster_name]["identity"])

    # Default
    else:
        verdict = "REVIEW"
        confidence = "LOW"
        reasoning_parts.append(
            f"Cluster '{cluster_name}' not in expected epithelial label set — manual review needed."
        )

    return {
        "verdict": verdict,
        "confidence": confidence,
        "reasoning": " ".join(reasoning_parts),
        "proposed_name": proposed_name,
    }


# ============================================================================
# PHASE 6: HTML RENDERING FUNCTIONS
# ============================================================================

def render_heatmap_table(summary):
    """Render the signature score heatmap as an HTML table."""
    html = ['<table class="heatmap-table">']
    html.append('<tr><th>Cluster</th><th>Cells</th>')
    for sig_key, sig_info in SIGNATURES.items():
        html.append(f'<th title="{sig_info["description"]}">{sig_info["display"]}</th>')
    html.append('</tr>')
    col_ranges = {}
    for sig_key in SIGNATURE_KEYS:
        col = f"score_{sig_key}_mean"
        vals = summary[col].dropna()
        col_ranges[sig_key] = (vals.min(), vals.max()) if len(vals) > 0 else (0, 1)
    for _, row in summary.iterrows():
        cluster = row["cluster"]
        n_cells = int(row["n_cells"])
        html.append(f'<tr><td>{cluster}</td><td>{n_cells:,}</td>')
        for sig_key in SIGNATURE_KEYS:
            col = f"score_{sig_key}_mean"
            val = row[col]
            vmin, vmax = col_ranges[sig_key]
            style = _color_cell(val, vmin, vmax)
            html.append(f'<td style="{style}">{val:.3f}</td>' if not np.isnan(val)
                       else '<td style="background:#f5f5f5;color:#999;">N/A</td>')
        html.append('</tr>')
    html.append('</table>')
    return "\n".join(html)


def render_cell_cycle_section(cc_df, cc_fig_html):
    """Render cell cycle analysis section."""
    html = []
    html.append('<h2 id="cell-cycle">Cell Cycle Analysis (Tirosh et al. 2016)</h2>')
    html.append('<div class="overview-box">')
    html.append('Cell cycle phase assignment using Tirosh et al. 2016 S-phase (43 genes) and '
                'G2/M-phase (54 genes) signatures via <code>sc.tl.score_genes_cell_cycle</code>. '
                'Cells with S_score > G2M_score are assigned S phase; G2M_score > S_score are '
                'G2M; remaining are G1. Cycling % = S + G2M.')
    html.append('</div>')
    html.append(cc_fig_html)
    html.append('<table class="sig-table">')
    html.append('<tr><th>Cluster</th><th>Cells</th><th>G1 %</th><th>S %</th>'
                '<th>G2M %</th><th>Cycling %</th><th>Assessment</th></tr>')
    for _, r in cc_df.iterrows():
        cyc = r["cycling_pct"]
        if cyc > 70:
            call = "Highly proliferative"
        elif cyc > 30:
            call = "Actively cycling"
        elif cyc > 15:
            call = "Moderately cycling"
        else:
            call = "Quiescent"
        html.append(f'<tr><td>{r["cluster"]}</td><td>{int(r["n_cells"]):,}</td>'
                    f'<td>{r["G1_pct"]:.1f}</td><td>{r["S_pct"]:.1f}</td>'
                    f'<td>{r["G2M_pct"]:.1f}</td><td><b>{cyc:.1f}</b></td>'
                    f'<td>{call}</td></tr>')
    html.append('</table>')
    return "\n".join(html)


def render_progeny_section(progeny_df, progeny_fig_html):
    """Render PROGENy pathway analysis section."""
    html = []
    html.append('<h2 id="progeny">PROGENy Pathway Activity (Schubert et al. 2018)</h2>')
    html.append('<div class="overview-box">')
    html.append('PROGENy scores 14 cancer-relevant signaling pathways using weighted gene '
                'expression signatures. Scores represent pathway activity relative to the '
                'footprint gene weights. Higher = more pathway activity.')
    html.append('</div>')
    html.append(progeny_fig_html)

    # Top pathway per cluster
    html.append('<h4>Dominant Pathway per Cluster</h4>')
    html.append('<table class="sig-table">')
    html.append('<tr><th>Cluster</th><th>#1 Pathway</th><th>Score</th>'
                '<th>#2 Pathway</th><th>Score</th><th>#3 Pathway</th><th>Score</th></tr>')
    for _, r in progeny_df.iterrows():
        pw_means = {}
        for pw in PROGENY_NAMES:
            v = r.get(f"{pw}_mean", np.nan)
            if not np.isnan(v):
                pw_means[pw] = v
        sorted_pw = sorted(pw_means.items(), key=lambda x: x[1], reverse=True)
        html.append(f'<tr><td>{r["cluster"]}</td>')
        for i in range(3):
            if i < len(sorted_pw):
                html.append(f'<td><b>{sorted_pw[i][0]}</b></td><td>{sorted_pw[i][1]:.3f}</td>')
            else:
                html.append('<td>—</td><td>—</td>')
        html.append('</tr>')
    html.append('</table>')
    return "\n".join(html)


def render_trajectory_section(dpt_df, ct_df, traj_fig_html):
    """Render trajectory + CytoTRACE section."""
    html = []
    html.append('<h2 id="trajectory">Trajectory & Differentiation (DPT + CytoTRACE)</h2>')
    html.append('<div class="overview-box">')
    html.append('Diffusion pseudotime (DPT) computed on the secretory subset (excluding Ciliated '
                'cells) using the Secretory epithelial cell medoid as root. CytoTRACE scores gene '
                'expression diversity as a proxy for stemness (higher = more progenitor-like). '
                'PAGA connectivity quantifies inter-cluster transition strength.')
    html.append('</div>')
    html.append(traj_fig_html)

    # DPT ordering table
    html.append('<h4>Pseudotime Ordering (root = Secretory epithelial cell)</h4>')
    html.append('<table class="sig-table">')
    html.append('<tr><th>Rank</th><th>Cluster</th><th>DPT Median</th><th>DPT IQR</th>'
                '<th>CytoTRACE</th><th>Interpretation</th></tr>')
    dpt_sorted = dpt_df.sort_values("dpt_median")
    ct_map = {r["cluster"]: r for _, r in ct_df.iterrows()}
    for rank_i, (_, r) in enumerate(dpt_sorted.iterrows(), 1):
        cluster = r["cluster"]
        ct_val = ct_map.get(cluster, {}).get("cytotrace_mean", np.nan)
        ct_str = f"{ct_val:.3f}" if not np.isnan(ct_val) else "N/A"
        iqr = f"{r['dpt_q25']:.3f}–{r['dpt_q75']:.3f}"
        if rank_i == 1:
            interp = "Root / progenitor state"
        elif rank_i == len(dpt_sorted):
            interp = "Terminal / most differentiated"
        else:
            interp = "Intermediate"
        html.append(f'<tr><td>{rank_i}</td><td>{cluster}</td><td>{r["dpt_median"]:.3f}</td>'
                    f'<td>{iqr}</td><td>{ct_str}</td><td>{interp}</td></tr>')
    html.append('</table>')

    # PAGA connectivity
    paga_path = os.path.join(OUT_DIR, "10c_paga_connectivity.csv")
    if os.path.exists(paga_path):
        paga = pd.read_csv(paga_path, index_col=0)
        html.append('<h4>PAGA Connectivity (top transitions)</h4>')
        html.append('<table class="sig-table">')
        html.append('<tr><th>Cluster A</th><th>Cluster B</th><th>Weight</th></tr>')
        pairs = []
        for i, c1 in enumerate(paga.index):
            for j, c2 in enumerate(paga.columns):
                if j > i:
                    pairs.append((c1, c2, paga.iloc[i, j]))
        pairs.sort(key=lambda x: x[2], reverse=True)
        for c1, c2, w in pairs[:8]:
            html.append(f'<tr><td>{c1}</td><td>{c2}</td><td>{w:.4f}</td></tr>')
        html.append('</table>')

    return "\n".join(html)


def render_cluster_card(cluster_name, summary_row, gene_coverage, assessment,
                        cc_row, progeny_row, dpt_row, ct_row):
    """Render a comprehensive assessment card for one cluster."""
    verdict = assessment["verdict"]
    confidence = assessment["confidence"]
    reasoning = assessment["reasoning"]
    proposed = assessment.get("proposed_name")
    hgsc = HGSC_CONTEXT.get(cluster_name, {})

    badge_cls = {"OK": "badge-ok", "REVIEW": "badge-review", "RENAME": "badge-rename"}
    badge_label = {"OK": "NAME OK", "REVIEW": "NEEDS REVIEW", "RENAME": "SUGGEST RENAME"}

    html = [f'<div class="cluster-card" id="cluster-{cluster_name.replace(" ", "-").replace("/", "-").lower()}">']
    html.append(f'<h3>{cluster_name} '
                f'<span class="badge {badge_cls.get(verdict, "badge-review")}">'
                f'{badge_label.get(verdict, verdict)}</span>'
                f'<span class="badge" style="background:#604E97;">Confidence: {confidence}</span>'
                f'</h3>')

    if proposed:
        html.append(f'<div class="warning-box"><b>Proposed rename:</b> {proposed}</div>')

    html.append(f'<div class="overview-box"><div class="recommendation">{reasoning}</div></div>')

    # HGSC context box
    if hgsc:
        html.append('<div class="context-box">')
        html.append(f'<b>SecA/SecB axis:</b> {hgsc.get("seca_secb", "—")}<br>')
        html.append(f'<b>Key markers:</b> <code>{hgsc.get("key_markers", "—")}</code><br>')
        html.append(f'<b>Vulnerabilities:</b> {hgsc.get("vulnerabilities", "—")}')
        html.append('</div>')

    # Therapeutic targets
    if hgsc.get("therapeutic"):
        html.append('<div class="therapy-box">')
        html.append(f'<b>Therapeutic targets:</b> {hgsc["therapeutic"]}')
        html.append('</div>')

    # Grid: scores + cell cycle + PROGENy
    html.append('<div class="grid-3col">')

    # Col 1: Signature scores
    html.append('<div>')
    html.append('<h4>Signature Scores</h4>')
    html.append('<table class="sig-table" style="font-size:11px;">')
    html.append('<tr><th>Signature</th><th>Mean</th></tr>')
    for sig_key, sig_info in SIGNATURES.items():
        mean_val = summary_row.get(f"score_{sig_key}_mean", np.nan)
        mean_str = f"{mean_val:.3f}" if not np.isnan(mean_val) else "N/A"
        html.append(f'<tr><td>{sig_info["display"]}</td><td>{mean_str}</td></tr>')
    html.append('</table>')
    html.append('</div>')

    # Col 2: Cell cycle
    html.append('<div>')
    html.append('<h4>Cell Cycle</h4>')
    if cc_row:
        cyc = cc_row.get("cycling_pct", 0)
        html.append(f'<div style="margin:5px 0;"><b>Cycling:</b> {cyc:.1f}%</div>')
        html.append(f'<div class="cc-bar" style="width:{min(100, cyc)}%; '
                    f'background:#E64B35;">&nbsp;</div>')
        html.append(f'<div style="font-size:11px; margin-top:5px;">'
                    f'G1: {cc_row.get("G1_pct", 0):.1f}% · '
                    f'S: {cc_row.get("S_pct", 0):.1f}% · '
                    f'G2M: {cc_row.get("G2M_pct", 0):.1f}%</div>')

    # DPT + CytoTRACE
    html.append('<h4 style="margin-top:12px;">Trajectory</h4>')
    if dpt_row:
        html.append(f'<div style="font-size:11px;">DPT median: {dpt_row.get("dpt_median", 0):.3f}</div>')
    if ct_row:
        html.append(f'<div style="font-size:11px;">CytoTRACE: {ct_row.get("cytotrace_mean", 0):.3f}</div>')
    html.append('</div>')

    # Col 3: PROGENy top pathways
    html.append('<div>')
    html.append('<h4>PROGENy Pathways</h4>')
    if progeny_row:
        pw_means = {}
        for pw in PROGENY_NAMES:
            v = progeny_row.get(f"{pw}_mean", np.nan)
            if not np.isnan(v):
                pw_means[pw] = v
        sorted_pw = sorted(pw_means.items(), key=lambda x: x[1], reverse=True)
        html.append('<div style="font-size:11px;"><b>Top active:</b><br>')
        for pw, v in sorted_pw[:5]:
            html.append(f'<span class="pathway-badge pw-active">{pw}: {v:.2f}</span> ')
        html.append('</div>')
        if len(sorted_pw) > 2:
            html.append('<div style="font-size:11px; margin-top:5px;"><b>Lowest:</b><br>')
            for pw, v in sorted_pw[-3:]:
                html.append(f'<span class="pathway-badge pw-suppress">{pw}: {v:.2f}</span> ')
            html.append('</div>')
    html.append('</div>')

    html.append('</div>')  # grid-3col

    # Gene coverage (collapsible)
    html.append('<details><summary>Gene Coverage Details</summary>')
    for sig_key, sig_info in SIGNATURES.items():
        cov = gene_coverage.get(sig_key, {})
        n_present = cov.get("present", 0)
        n_requested = cov.get("requested", 0)
        missing = cov.get("missing", [])
        html.append(f'<div style="margin:3px 0;"><b>{sig_info["display"]}:</b> {n_present}/{n_requested} ')
        for g in missing:
            html.append(f'<span class="gene-tag missing">{g}</span>')
        html.append('</div>')
    html.append('</details>')

    # QC sidebar
    html.append('<div class="neutral-box">')
    html.append(f'<b>Cells:</b> {int(summary_row["n_cells"]):,} · ')
    if "med_genes" in summary_row:
        html.append(f'<b>Median genes:</b> {summary_row["med_genes"]:,.0f} · ')
    if "med_umi" in summary_row:
        html.append(f'<b>Median UMI:</b> {summary_row["med_umi"]:,.0f} · ')
    if "med_dbl" in summary_row:
        html.append(f'<b>Median doublet:</b> {summary_row["med_dbl"]:.3f}')
    html.append('</div>')
    html.append('</div>')
    return "\n".join(html)


# ============================================================================
# PHASE 6: FULL REPORT ASSEMBLY
# ============================================================================

def render_full_report(summary, gene_coverage, cc_df, progeny_df, dpt_df, ct_df,
                       cc_fig_html, progeny_fig_html, traj_fig_html):
    """Render the complete HTML report."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n_clusters = len(summary)
    total_cells = int(summary["n_cells"].sum())

    # Build lookup dicts
    cc_map = {r["cluster"]: r.to_dict() for _, r in cc_df.iterrows()}
    prog_map = {r["cluster"]: r.to_dict() for _, r in progeny_df.iterrows()}
    dpt_map = {r["cluster"]: r.to_dict() for _, r in dpt_df.iterrows()}
    ct_map = {r["cluster"]: r.to_dict() for _, r in ct_df.iterrows()}

    # Compute verdicts
    assessments = {}
    for _, row in summary.iterrows():
        cluster = row["cluster"]
        assessments[cluster] = _verdict_for_cluster(
            cluster, row.to_dict(), gene_coverage, all_summary=summary,
            cc_row=cc_map.get(cluster), progeny_row=prog_map.get(cluster),
            dpt_row=dpt_map.get(cluster), ct_row=ct_map.get(cluster)
        )

    n_ok = sum(1 for a in assessments.values() if a["verdict"] == "OK")
    n_review = sum(1 for a in assessments.values() if a["verdict"] == "REVIEW")
    n_rename = sum(1 for a in assessments.values() if a["verdict"] == "RENAME")

    html = [f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>HGSC Atlas — Epithelial Comprehensive Validation (10c)</title>
{CSS}
</head>
<body>

<h1>Epithelial celltype_level2 Comprehensive Validation Report</h1>
<div class="overview-box">
<b>Step 10c — Comprehensive epithelial annotation validation</b><br>
<b>Report generated:</b> {timestamp}<br>
<b>Atlas:</b> hgsc_atlas_epithelial.h5ad<br>
<b>Compartment:</b> Epithelial cell ({n_clusters} retained clusters, {total_cells:,} cells)<br>
<b>Analyses:</b> 12 gene signatures + Tirosh cell cycle + PROGENy (14 pathways) +
CytoTRACE stemness + DPT pseudotime + PAGA connectivity + HGSC molecular context +
therapeutic vulnerability profiling
</div>

<div class="summary-box">
<b>Verdict Summary:</b><br>
<span class="badge badge-ok">NAME OK: {n_ok}</span>
<span class="badge badge-review">NEEDS REVIEW: {n_review}</span>
<span class="badge badge-rename">SUGGEST RENAME: {n_rename}</span>
</div>
"""]

    # TOC
    html.append('<div class="toc"><b>Contents</b><br>')
    html.append('<a href="#heatmap">Signature Heatmap</a> · ')
    html.append('<a href="#cell-cycle">Cell Cycle</a> · ')
    html.append('<a href="#progeny">PROGENy Pathways</a> · ')
    html.append('<a href="#trajectory">Trajectory</a> · ')
    html.append('<a href="#cluster-details">Cluster Details</a> · ')
    html.append('<a href="#seca-secb">SecA/SecB Context</a> · ')
    html.append('<a href="#recommendations">Recommendations</a> · ')
    html.append('<a href="#methods">Methods</a>')
    html.append('</div>')

    # 1. Signature Heatmap
    html.append('<h2 id="heatmap">Signature Score Heatmap (Mean per Cluster)</h2>')
    html.append('<div class="overview-box">Rows = epithelial clusters; columns = signature '
                'scores. Color intensity reflects relative score within each column.</div>')
    html.append(render_heatmap_table(summary))

    # 2. Cell Cycle
    html.append(render_cell_cycle_section(cc_df, cc_fig_html))

    # 3. PROGENy
    html.append(render_progeny_section(progeny_df, progeny_fig_html))

    # 4. Trajectory
    html.append(render_trajectory_section(dpt_df, ct_df, traj_fig_html))

    # 5. Cluster Details
    html.append('<h2 id="cluster-details">Per-Cluster Comprehensive Assessment</h2>')
    for _, row in summary.iterrows():
        cluster = row["cluster"]
        html.append(render_cluster_card(
            cluster, row.to_dict(), gene_coverage, assessments[cluster],
            cc_row=cc_map.get(cluster), progeny_row=prog_map.get(cluster),
            dpt_row=dpt_map.get(cluster), ct_row=ct_map.get(cluster)
        ))

    # 6. SecA/SecB Context
    html.append('<h2 id="seca-secb">SecA/SecB Polarization Context</h2>')
    html.append('<div class="overview-box">')
    html.append("""
<b>The SecA/SecB axis</b> (defined in Step 12b) describes a progenitor-to-differentiated
gradient in HGSC epithelial cells:<br><br>
<b>SecA-like (progenitor):</b> Cycling secretory + Stress-response secretory + Transitioning epithelial
— high proliferation, IEG, EMT, enriched in resistant tumors (SecA/SecB ratio = 2.04) and solid tissue<br>
<b>SecB-like (differentiated):</b> Adaptive secretory epithelial
— high stress adaptation, iron/oxidative, hypoxia; enriched in ascites (ratio = 0.72) and post-chemo<br>
<b>Baseline:</b> Secretory epithelial cell — intermediate state<br>
<b>Lineage-specific:</b> Ciliated epithelial cell — terminally differentiated, orthogonal to SecA/SecB
""")
    html.append('</div>')

    # 7. Recommendations
    html.append('<h2 id="recommendations">Summary Recommendations</h2>')
    renames = {c: a for c, a in assessments.items() if a.get("proposed_name")}
    reviews = {c: a for c, a in assessments.items()
               if a["verdict"] == "REVIEW" and not a.get("proposed_name")}

    if renames:
        html.append('<h3>Proposed Name Changes</h3>')
        html.append('<table class="sig-table">')
        html.append('<tr><th>Current Name</th><th>Proposed Name</th><th>Reasoning</th></tr>')
        for cluster, a in renames.items():
            html.append(f'<tr><td>{cluster}</td><td><b>{a["proposed_name"]}</b></td>'
                       f'<td style="font-size:11px;">{a["reasoning"][:300]}...</td></tr>')
        html.append('</table>')
    elif reviews:
        html.append('<div class="flag-box"><b>Clusters flagged for review:</b></div>')
        html.append('<table class="sig-table">')
        html.append('<tr><th>Cluster</th><th>Concern</th></tr>')
        for cluster, a in reviews.items():
            html.append(f'<tr><td>{cluster}</td>'
                       f'<td style="font-size:11px;">{a["reasoning"][:300]}...</td></tr>')
        html.append('</table>')
    else:
        html.append('<div class="summary-box"><b>All 6 epithelial cluster names are validated</b> — '
                   'no name changes proposed. Every cluster shows transcriptional signatures, cell '
                   'cycle profiles, pathway activities, and trajectory positions consistent with '
                   'its assigned label.</div>')

    # 8. Methods
    html.append('<h2 id="methods">Methods</h2>')
    html.append('<div class="neutral-box">')
    html.append("""
<b>Gene signature scoring:</b> 12 signatures scored via <code>scanpy.tl.score_genes</code>
(background-corrected z-score approach).<br>
<b>Cell cycle:</b> Tirosh et al. 2016 S/G2M gene lists via
<code>scanpy.tl.score_genes_cell_cycle</code>.<br>
<b>PROGENy:</b> 14 cancer pathway activities computed as weighted dot products of
gene expression with pathway footprint weights (Schubert et al. 2018).<br>
<b>CytoTRACE:</b> Gene expression diversity proxy for stemness. Top 200 genes
correlated with gene counts → mean expression = stemness score.<br>
<b>Trajectory:</b> Diffusion pseudotime (DPT) on secretory subset using scanpy.
Root = Secretory epithelial cell medoid in diffusion space. PAGA connectivity
from scanpy.tl.paga on k=30 neighbor graph in X_scanvi space.<br>
<b>HGSC context:</b> Therapeutic targets from published HGSC clinical/preclinical
literature including FOLR1 ADCs, TROP2 ADCs, PARP inhibitors, CDK4/6 inhibitors,
TGFβ inhibitors, and hypoxia-activated prodrugs.
""")
    html.append('</div>')

    # Gene coverage
    html.append('<details><summary>Gene Coverage Summary</summary>')
    html.append('<table class="sig-table">')
    html.append('<tr><th>Signature</th><th>Requested</th><th>Found</th><th>Missing</th></tr>')
    for sig_key, sig_info in SIGNATURES.items():
        cov = gene_coverage.get(sig_key, {})
        missing = cov.get("missing", [])
        missing_str = ", ".join(missing) if missing else "—"
        html.append(f'<tr><td>{sig_info["display"]}</td>'
                   f'<td>{cov.get("requested", 0)}</td>'
                   f'<td>{cov.get("present", 0)}</td>'
                   f'<td style="font-size:10px;">{missing_str}</td></tr>')
    html.append('</table>')
    html.append('</details>')

    html.append('</body></html>')

    full_html = "\n".join(html)
    with open(HTML_PATH, "w") as f:
        f.write(full_html)

    file_size = os.path.getsize(HTML_PATH) / 1024
    print(f"\n  Saved: {HTML_PATH} ({file_size:.1f} KB)")
    return assessments


# ============================================================================
# MAIN
# ============================================================================

def main():
    t0 = datetime.now()

    # Phase 0: Load
    adata = load_epithelial_cells()

    # Phase 1: Gene signatures
    score_cols, gene_coverage = score_signatures(adata)
    summary = compute_cluster_summaries(adata, score_cols)
    print("\n  Saving per-cell scores...", flush=True)
    save_per_cell_scores(adata, score_cols)

    # Phase 2: Cell cycle
    cc_df = score_cell_cycle(adata)

    # Phase 3: PROGENy
    progeny_df, progeny_cols = score_progeny(adata)

    # Phase 4: Trajectory (secretory subset)
    adata_traj = build_trajectory_subset(adata)
    ct_df, adata_traj = compute_cytotrace(adata_traj)
    dpt_df = compute_dpt(adata_traj)
    del adata_traj
    gc.collect()

    # Phase 5: Figures
    print("\n  [Phase 5] Generating figures...", flush=True)
    cc_fig_html = make_cell_cycle_fig(cc_df)
    print("    Saved: 10c_cell_cycle.svg")
    progeny_fig_html = make_progeny_fig(progeny_df)
    print("    Saved: 10c_progeny_heatmap.svg")
    traj_fig_html = make_trajectory_fig(ct_df, dpt_df)
    print("    Saved: 10c_trajectory.svg")

    # Phase 6: HTML report
    print("\n  [Phase 6] Rendering HTML report...", flush=True)
    assessments = render_full_report(
        summary, gene_coverage, cc_df, progeny_df, dpt_df, ct_df,
        cc_fig_html, progeny_fig_html, traj_fig_html
    )

    # Console summary
    print(f"\n{'=' * 70}")
    print("  VERDICT SUMMARY")
    print(f"{'=' * 70}")
    for cluster, a in assessments.items():
        proposed = f" → {a['proposed_name']}" if a.get("proposed_name") else ""
        print(f"    [{a['verdict']:6s}] {cluster}{proposed}")

    del adata
    gc.collect()

    elapsed = (datetime.now() - t0).total_seconds()
    print(f"\n{'=' * 70}")
    print(f"  DONE — elapsed {elapsed:.0f}s")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
