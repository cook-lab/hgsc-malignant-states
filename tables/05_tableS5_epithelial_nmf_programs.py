#!/usr/bin/env python3
"""Supplementary Table 5 — epithelial NMF programs (k = 10), SecA/SecB gene classification, consensus metaprograms.

Rebuilds the whole workbook from the atlas NMF outputs that the final figures read (read-only):
  ATLAS/11d_epithelial_nmf/11d_nmf_loadings.csv     H matrix (10 factors x 3,007 genes), canonical step 11d
                                                    (scripts/11d_epithelial_nmf.py, 2026-02-27). The same run's
                                                    11d_nmf_usage.csv feeds Figures 1h-m, S5 and S6.
  ATLAS/11d_epithelial_nmf/11d_factor_scores.csv    signature-gene overlap per factor (step 11d)
  ATLAS/11e_nmf_characterization/11e_gene_classification.csv   90th-percentile SecA/SecB classification (step 11e);
                                                    recomputed here from 11d and asserted identical
  ATLAS/11f_consensus_nmf/11f_mp_summary.csv, 11f_mp_top_genes.csv, 11f_patient_program_counts.csv
                                                    per-patient consensus NMF (k = 5-10, 30 runs) metaprograms
  ATLAS/11g_secAB_consensus_nmf/11g_mp_summary.csv, 11g_mp_top_genes.csv
                                                    consensus NMF restricted to the 354 SecA/SecB-specific genes
  figures/supplementary/11_figS5_nmf_factor_umaps.py   FACTOR_NAMES (the labels printed on Figure S5)
  figures/figure1/07_fig1j_program_factor_umaps.py         PANELS (Fig 1j factor labels)
OUTPUTS: SuppTable5_epithelial_nmf_programs.xlsx (+ one CSV per sheet) via supp_xlsx.write_table.
MANUSCRIPT ITEM(S): Table S5 (k = 10, 3,007 genes; 177 / 177 / 124 classification; consensus NMF).
Run:  python 05_tableS5_epithelial_nmf_programs.py
"""
import ast
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from supp_xlsx import write_table  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, ".."))
ATLAS = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "2026_final_atlas/output")
F_LOAD = f"{ATLAS}/11d_epithelial_nmf/11d_nmf_loadings.csv"
F_SCORES = f"{ATLAS}/11d_epithelial_nmf/11d_factor_scores.csv"
F_CLASS = f"{ATLAS}/11e_nmf_characterization/11e_gene_classification.csv"
F_11F = f"{ATLAS}/11f_consensus_nmf"
F_11G = f"{ATLAS}/11g_secAB_consensus_nmf"
F_SF5 = f"{REPO}/figures/supplementary/11_figS5_nmf_factor_umaps.py"
F_FIG1 = f"{REPO}/figures/figure1/07_fig1j_program_factor_umaps.py"

SECB, SECA = "Factor_2", "Factor_3"
# Signature genes the canonical 11d run forced into the HVG set (scripts/11d_epithelial_nmf.py lines 88-89)
CANON_11D_SECA = ["BCAM", "MECOM", "FBXO21", "LGR5", "LPAR3", "PBX1", "SOX17", "WT1", "RCN2"]
CANON_11D_SECB = ["KRT17", "KRT19", "KRT7", "LCN2", "MMP7", "PRSS22", "SLPI", "TACSTD2"]

# Labels written on the atlas step-11d annotated overview (11d_all_factors_umap_annotated.png, 2026-03-24),
# transcribed from the image because no script in the atlas tree writes that file.
ATLAS_11D_LABEL = {
    "Factor_1": "Pan-epithelial / metabolic", "Factor_2": "SecB pole (keratins, differentiation)",
    "Factor_3": "SecA pole (stemness, TFs)", "Factor_4": "Interferon response",
    "Factor_5": "Cell cycle / proliferation", "Factor_6": "IEG / immediate early genes",
    "Factor_7": "Progenitor (IGFBP2, BCAM, SOX4)", "Factor_8": "Secretory / CLU / stress",
    "Factor_9": "EMT-like (VIM, FN1, LGALS1)", "Factor_10": "Metallothionein stress",
}
FIG_PANEL = {"Factor_2": "Fig 1i; SF5; SF6b", "Factor_3": "Fig 1h; SF5; SF6a", "Factor_5": "Fig 1j; SF5",
             "Factor_6": "Fig 1k; SF5", "Factor_4": "Fig 1l; SF5", "Factor_9": "Fig 1m; SF5"}


def _literal(path, name):
    """Read a module-level literal (dict/list) from a generator without executing it."""
    tree = ast.parse(open(path).read())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise KeyError(f"{name} not found in {path}")


def main():
    L = pd.read_csv(F_LOAD, index_col=0)                     # factors x genes
    factors = list(L.index)
    genes = list(L.columns)
    assert len(factors) == 10 and len(genes) == 3007, L.shape
    spec = L / L.sum(axis=0).replace(0, np.nan)              # share of each gene's total loading held by the factor
    norm = L.div(L.max(axis=1), axis=0)                      # loading / factor max (0-1)
    score_spec = (L * spec).fillna(0)                        # loading x specificity: ranks factor-distinctive genes

    sf5 = _literal(F_SF5, "FACTOR_NAMES")
    fig1 = {p[1]: p[2] for p in _literal(F_FIG1, "PANELS")}
    fs = pd.read_csv(F_SCORES, index_col=0)

    # ---------------- Factor_annotation ----------------
    rows = []
    for f in factors:
        name, markers = (sf5[f].split("\n") + [""])[:2]
        rows.append({
            "factor": f,
            "program_SF5_label": name.replace("—", "-"),
            "SF5_marker_genes_shown": markers.strip("()"),
            "Fig1_label": fig1.get(f, {"Factor_2": "SecB", "Factor_3": "SecA"}.get(f, "")),
            "atlas_11d_annotation": ATLAS_11D_LABEL[f],
            "figure_panels": FIG_PANEL.get(f, "SF5"),
            "top10_genes_by_loading": ", ".join(L.loc[f].nlargest(10).index),
            "top10_genes_by_specificity_weighted_loading": ", ".join(score_spec.loc[f].nlargest(10).index),
            "n_genes_nonzero_loading": int((L.loc[f] > 0).sum()),
            "SecA_sig_genes_in_top100 (11d 9-gene set)": int(fs.loc[f, "secA_top100"]),
            "SecB_sig_genes_in_top100 (11d 8-gene set)": int(fs.loc[f, "secB_top100"]),
        })
    annot = pd.DataFrame(rows)

    # ---------------- Loadings_long (all 3,007 genes x 10 factors) ----------------
    long = (L.reset_index().rename(columns={"index": "factor"})
            .melt(id_vars="factor", var_name="gene", value_name="loading"))
    long["rank_in_factor"] = long.groupby("factor")["loading"].rank(ascending=False, method="first").astype(int)
    long = long.merge(norm.reset_index().rename(columns={"index": "factor"})
                      .melt(id_vars="factor", var_name="gene", value_name="loading_norm_to_factor_max"),
                      on=["factor", "gene"])
    long = long.merge(spec.fillna(0).reset_index().rename(columns={"index": "factor"})
                      .melt(id_vars="factor", var_name="gene", value_name="specificity_share_of_gene_total"),
                      on=["factor", "gene"])
    long["factor"] = pd.Categorical(long["factor"], categories=factors, ordered=True)
    long = long.sort_values(["factor", "rank_in_factor"])[
        ["factor", "rank_in_factor", "gene", "loading", "loading_norm_to_factor_max",
         "specificity_share_of_gene_total"]].reset_index(drop=True)
    long["factor"] = long["factor"].astype(str)
    for c in ["loading", "loading_norm_to_factor_max", "specificity_share_of_gene_total"]:
        long[c] = long[c].round(6)

    # ---------------- Top50_wide ----------------
    top50 = pd.DataFrame({f: L.loc[f].nlargest(50).index.tolist() for f in factors})
    top50.insert(0, "rank", range(1, 51))
    top50s = pd.DataFrame({f: score_spec.loc[f].nlargest(50).index.tolist() for f in factors})
    top50s.insert(0, "rank", range(1, 51))

    # ---------------- SecA_SecB_classification (recomputed, then checked against 11e) ----------------
    f2, f3 = L.loc[SECB], L.loc[SECA]
    t2, t3 = np.percentile(f2.values, 90), np.percentile(f3.values, 90)
    cls = np.where((f3 > t3) & (f2 <= t2), "SecA-specific",
          np.where((f2 > t2) & (f3 <= t3), "SecB-specific",
          np.where((f2 > t2) & (f3 > t3), "Shared", "Background")))
    rank3 = f3.rank(ascending=False, method="first").astype(int)
    rank2 = f2.rank(ascending=False, method="first").astype(int)
    C = pd.DataFrame({
        "gene": genes, "class": cls,
        "F3_SecA_loading": f3.values.round(6), "F3_rank": rank3.values,
        "F3_percentile": (f3.rank(pct=True, method="max") * 100).round(2).values,
        "F2_SecB_loading": f2.values.round(6), "F2_rank": rank2.values,
        "F2_percentile": (f2.rank(pct=True, method="max") * 100).round(2).values,
    })
    C["in_forced_11d_signature"] = np.where(C.gene.isin(CANON_11D_SECA), "SecA (11d)",
                                   np.where(C.gene.isin(CANON_11D_SECB), "SecB (11d)", ""))
    ref = pd.read_csv(F_CLASS, index_col=0).set_index("gene")
    chk = C.set_index("gene")
    assert (ref.loc[chk.index, "class"] == chk["class"]).all(), "classification differs from 11e output"
    assert np.allclose(ref.loc[chk.index, "secA_loading"], f3.values, atol=1e-5)
    counts = C["class"].value_counts()
    assert (counts["SecA-specific"], counts["SecB-specific"], counts["Shared"]) == (177, 177, 124), counts
    order = {"SecA-specific": 0, "SecB-specific": 1, "Shared": 2, "Background": 3}
    C["_o"] = C["class"].map(order)
    C["_s"] = np.where(C["class"] == "SecB-specific", -C.F2_SecB_loading, -C.F3_SecA_loading)
    C = C.sort_values(["_o", "_s"]).drop(columns=["_o", "_s"]).reset_index(drop=True)

    # ---------------- Consensus metaprograms (11f; text) + restricted rerun (11g) ----------------
    s11f = pd.read_csv(f"{F_11F}/11f_mp_summary.csv")
    tg11f = pd.read_csv(f"{F_11F}/11f_mp_top_genes.csv")
    s11g = pd.read_csv(f"{F_11G}/11g_mp_summary.csv")
    tg11g = pd.read_csv(f"{F_11G}/11g_mp_top_genes.csv")
    pp = pd.read_csv(f"{F_11F}/11f_patient_program_counts.csv")
    s11f.insert(1, "description", s11f.metaprogram.map({"MP1": "General secretory program",
                                                          "MP2": "Metallothionein response program"}))
    cons = s11f.drop(columns=["top10_genes"])
    cons_genes = tg11f.rename(columns={"weight": "gene_weight"})
    cons_genes["gene_weight"] = cons_genes["gene_weight"].round(6)
    s11g = s11g.drop(columns=["top10_genes"])
    g11g = tg11g.rename(columns={"weight": "gene_weight", "gene_class": "11e_class"})
    g11g["gene_weight"] = g11g["gene_weight"].round(6)
    g11g["11e_class"] = g11g["11e_class"].map({"SecA": "SecA-specific", "SecB": "SecB-specific"})

    n_nz = int((L > 0).sum().sum())
    readme = [
        "Source (all read-only): atlas NMF outputs in 202601_final_atlas/output/ - 11d_epithelial_nmf/11d_nmf_loadings.csv "
        "(canonical step 11d, run 2026-02-27; the same run's 11d_nmf_usage.csv is what Fig 1h-m, SF5 and SF6 plot), "
        "11d_factor_scores.csv, 11e_nmf_characterization/11e_gene_classification.csv, 11f_consensus_nmf/, 11g_secAB_consensus_nmf/.",
        "Method (step 11d, scripts/11d_epithelial_nmf.py): all 575,366 atlas epithelial cells (hgsc_atlas_final_epithelial.h5ad, "
        "secretory + 6,678 ciliated), counts normalized to 10,000 and log1p; 3,000 HVGs (scanpy seurat_v3 on counts) plus the "
        "11d SecA/SecB signature genes forced in (BCAM, MECOM, FBXO21, LGR5, LPAR3, PBX1, SOX17, WT1, RCN2 / KRT17, KRT19, KRT7, "
        "LCN2, MMP7, PRSS22, SLPI, TACSTD2) = 3,007 genes. sklearn NMF k = 10 (init nndsvda, max_iter 500, seed 42) fitted on "
        "a 50,000-cell training set stratified by level-2 epithelial subtype, then all cells projected; usage normalized per cell.",
        "Factor identity: Factor_3 = SecA (progenitor-like) and Factor_2 = SecB (adaptive), as in Fig 1h/i, SF5 and SF6 and config "
        "polarization.factor = Factor_2.",
        "Sheet Factor_annotation: SF5 label and marker genes (as printed on SF5, from 11_figS5_nmf_factor_umaps.py), Fig 1 label "
        "(07_fig1j_program_factor_umaps.py), the atlas 11d annotation, the figure panels showing each factor, top 10 genes by raw "
        "loading and by specificity-weighted loading, and signature-gene overlap in each factor's top 100 (11d_factor_scores.csv).",
        f"Sheet Loadings_long: ALL 3,007 genes for every factor (30,070 rows; {n_nz:,} non-zero). Columns: factor, rank_in_factor "
        "(1 = highest loading), gene, loading (raw H-matrix weight), loading_norm_to_factor_max (loading / the factor's maximum "
        "loading, 0-1), specificity_share_of_gene_total (loading / sum of that gene's loadings over the 10 factors). The previous "
        "Supp Data 4 CSV (supplemental/supplementary_table_NMF_factor_genes.csv, 2026-05-22) is the same matrix restricted to "
        "the 17,083 non-zero loadings; values agree to 5e-7.",
        "Sheets Top50_wide / Top50_specific_wide: row 1 = factor, then the 50 highest-ranked genes per factor, ranked by raw "
        "loading and by loading x specificity respectively. Raw loadings rank abundant genes (ribosomal, FTH1, GAPDH, S100A6) "
        "highly in several factors; the specificity ranking shows the genes that set each factor apart.",
        "Sheet SecA_SecB_classification: all 3,007 genes. A gene is SecA-specific if its Factor_3 loading is above the 90th "
        "percentile of Factor_3 loadings and its Factor_2 loading is not above the Factor_2 90th percentile; SecB-specific is the "
        "mirror; Shared = above both; Background = neither. Recomputed from 11d_nmf_loadings.csv and checked to be identical to "
        "11e_gene_classification.csv: 177 SecA-specific, 177 SecB-specific, 124 Shared, 2,529 Background (matches methods). "
        "F3_percentile / F2_percentile = % of the 3,007 genes whose loading in that factor is <= the gene's loading (tied values, including zero loadings, share the upper value). These 177 + 177 lists are the SecA-/SecB-specific "
        "lists used for the TCGA-OV signatures (see Supp Table 6).",
        "Sheets Consensus_metaprograms / Consensus_MP_genes (step 11f, scripts/11f_consensus_nmf.py): per-patient NMF on the same "
        "3,007 genes for patients with >= 200 epithelial cells, k = 5-10, 30 runs per k, per-patient consensus programs by "
        "average-linkage clustering of cosine distance; programs pooled across patients and clustered (cosine threshold 0.3, "
        f">= 5 programs from >= 3 patients) into metaprograms. {len(pp)} patients, {int(pp.n_programs.sum())} programs, 2 metaprograms: "
        "MP1 general secretory (410 programs, 138 patients, 13 studies) and MP2 metallothionein (17 programs, 17 patients). "
        "Top 50 genes per metaprogram by weight.",
        "Sheets SecAB_restricted_MPs / SecAB_restricted_MP_genes (step 11g, scripts/11g_secAB_consensus_nmf.py): the same "
        "consensus procedure restricted to the 354 SecA- plus SecB-specific genes (k = 2-6). Three metaprograms (two SecA-leaning, "
        "one SecB, 12 patients). This analysis is not described in the manuscript text and is included for completeness.",
        "Signature genes in this classification: WT1, FBXO21, RCN2, PBX1, SOX17, MECOM and LPAR3 are SecA-specific; KRT7, LCN2, TACSTD2, SLPI, KRT17, MMP7 and PRSS22 are SecB-specific; KRT19 and BCAM are Shared; LGR5 is Background (Factor_3 rank 595).",
        "Notes: the NMF used 3,007 genes (3,000 highly variable genes plus 7 forced-in signature genes) and all epithelial "
        "cells, including ciliated; it was fitted on a 50,000-cell stratified training set and projected onto all 575,366 "
        "cells. The forced-in genes were the 9 + 8 11d sets (with BCAM, WT1 and MMP7), not the 7 + 7 scoring signature in "
        "shared/signatures.yml.",
    ]
    write_table(5, "epithelial_nmf_programs",
                title="Supplementary Table 5. Epithelial NMF programs (k = 10): factor annotation, gene loadings, "
                      "SecA/SecB gene classification and consensus metaprograms",
                readme=readme,
                sheets=[("Factor_annotation", annot), ("Loadings_long", long), ("Top50_wide", top50),
                        ("Top50_specific_wide", top50s), ("SecA_SecB_classification", C),
                        ("Consensus_metaprograms", cons), ("Consensus_MP_genes", cons_genes),
                        ("SecAB_restricted_MPs", s11g), ("SecAB_restricted_MP_genes", g11g)])
    print(counts.to_string())


if __name__ == "__main__":
    main()
