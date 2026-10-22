#!/usr/bin/env python3
"""
Fig 2b prep — cell-cycle phase of every epithelial cell, from the final atlas
=============================================================================
PURPOSE
    Assigns a cell-cycle phase to every epithelial cell of the final atlas and
    summarises the phases per epithelial state (Fig 2b; text: SecA 37.8% S/G2M,
    SecB 80.6% G1).

      * Expression: raw counts (layers/counts) of all 575,366 epithelial cells,
        normalize_total(target_sum=1e4) -> log1p.
      * scanpy.tl.score_genes_cell_cycle with the Tirosh et al. 2016 S (43) and
        G2/M (54) lists (the regev_lab_cell_cycle_genes list used in the scanpy
        tutorials; identical to Seurat cc.genes). Genes absent from the
        atlas var are dropped (reported). scanpy defaults (n_bins=25,
        ctrl_size=len(list), random_state=0); phase = scanpy default rule
        (G1 if both scores < 0, else the larger of S / G2M).
      * State from load_epithelial()['state'] (SecA / Intermediate / SecB / Ciliated).

INPUTS
    data_root/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
        (layers/counts only; its obs metadata are NOT used)
    figures/_prep/atlas_final_tables.py  (load_epithelial)

OUTPUTS (output_root/figure2/fig2b_cell_cycle/)
    cell_cycle_per_cell.parquet          barcode, S_score, G2M_score, phase, state
    cell_cycle_phase_proportions.csv     state x {G1, S, G2M} fractions (input of
                                         figures/figure2/02_fig2bcde_cellcycle_composition_paired.py)
    cell_cycle_phase_counts.csv          state x phase cell counts

MANUSCRIPT PANEL(S): Fig 2b (input).
RUNTIME TIER: moderate (loads the epithelial counts layer, ~10-20 min on SMB, ~25 GB RAM).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402
from figures._prep.atlas_final_tables import STATE_ORDER, load_epithelial  # noqa: E402

EPI_H5AD = path("data_root", "2026_final_atlas", "celltype_h5ad", "hgsc_atlas_final_epithelial.h5ad")
OUT_DIR = path("output_root", "figure2", "fig2b_cell_cycle")
os.makedirs(OUT_DIR, exist_ok=True)

# Tirosh et al. 2016 (regev_lab_cell_cycle_genes.txt: first 43 = S, remaining 54 = G2/M)
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP",
    "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
    "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
    "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B",
    "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
    "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN",
    "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]
assert len(S_GENES) == 43 and len(G2M_GENES) == 54

import anndata as ad  # noqa: E402
import h5py  # noqa: E402
import scanpy as sc  # noqa: E402
from anndata.io import read_elem  # noqa: E402

epi = load_epithelial(columns=["state"])
print("[1] reading epithelial counts layer ...", flush=True)
with h5py.File(EPI_H5AD, "r") as f:
    idx = pd.Index(np.asarray(read_elem(f["obs"][f["obs"].attrs["_index"]])).astype(str))
    genes = np.asarray(read_elem(f["var"][f["var"].attrs["_index"]])).astype(str)
    X = read_elem(f["layers/counts"])
assert idx.equals(epi.index), "epithelial h5ad order != load_epithelial() order"
X = X.astype(np.float32)
adata = ad.AnnData(X=X, obs=pd.DataFrame(index=idx), var=pd.DataFrame(index=genes))
print(f"    {adata.n_obs:,} cells x {adata.n_vars:,} genes; nnz = {X.nnz:,}", flush=True)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
s_use = [g for g in S_GENES if g in adata.var_names]
g_use = [g for g in G2M_GENES if g in adata.var_names]
print(f"    S genes present {len(s_use)}/43 (missing {sorted(set(S_GENES) - set(s_use))}); "
      f"G2M present {len(g_use)}/54 (missing {sorted(set(G2M_GENES) - set(g_use))})", flush=True)
print("[2] score_genes_cell_cycle ...", flush=True)
sc.tl.score_genes_cell_cycle(adata, s_genes=s_use, g2m_genes=g_use)

cc = adata.obs[["S_score", "G2M_score", "phase"]].copy()
cc["state"] = epi.loc[cc.index, "state"].values
cc.index.name = "barcode"
cc.to_parquet(os.path.join(OUT_DIR, "cell_cycle_per_cell.parquet"))
counts = pd.crosstab(cc.state, cc.phase).reindex(index=STATE_ORDER, columns=["G1", "S", "G2M"])
props = counts.div(counts.sum(axis=1), axis=0)
counts.to_csv(os.path.join(OUT_DIR, "cell_cycle_phase_counts.csv"))
props.to_csv(os.path.join(OUT_DIR, "cell_cycle_phase_proportions.csv"))
print((props * 100).round(1))
