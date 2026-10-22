#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/scripts/27_final_annotation_split.py, the script that ran.
# Kept: kpca_annotated.h5ad (final cell types, 20,598 cells; read by 31, 32 and 33). Removed: the label table and
# the annotation figure. The untrimmed script, which the docstring below describes, is in the lab project folder.

"""
Step 27 — Final annotation with full immune split (DEG-based)
=============================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

Top-DEG review (step 26 follow-up) corrected the coarse immune calls:
  cl0  Neutrophil            (S100a8/9, Csf3r, Il1b, Acod1)
  cl1  Neutrophil (cycling)  (Camp/Ngp/Ltf + Mki67/Top2a) — n=36, under review
  cl3  Macrophage            (C1qa/b/c, Mrc1, Csf1r, Trem2)
  cl4  Monocyte              (Ccr2, Plac8, Itgam, Vcan)
  cl13 Dendritic cell        (Zbtb46, Clec9a, Xcr1, Flt3, Ccr7) — was mislabelled "Mast"
Tumour = 5,6,8,9,11,12,15 ; Fibroblast = 10 ; Endothelial = 2 ; T/NK = 7 ; B = 14.

cl1 is drawn distinctly + Mki67 shown so Sarah can decide keep-vs-merge on the UMAP.

Input:   output/21_working/kpca_working.h5ad, output/24_whole_umap/24_allcell_ucell.csv
Output:  output/26_final_annotation/  (kpca_annotated.h5ad updated, split figure)
"""
import os, sys
from importlib import import_module
import numpy as np, pandas as pd, scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")
OUT = os.path.join(S.OUTPUT, "26_final_annotation")

# cl1 (36 cycling granulocytes) DROPPED per Sarah 2026-09-09 — too few to stand alone.
DROP = {"1"}
CELLTYPE = {
    "5": "Tumour", "6": "Tumour", "8": "Tumour", "9": "Tumour",
    "11": "Tumour", "12": "Tumour", "15": "Tumour",
    "10": "Fibroblast", "2": "Endothelial",
    "3": "Macrophage", "4": "Monocyte", "0": "Neutrophil",
    "13": "Dendritic", "7": "T/NK", "14": "B",
}
ORDER = ["Tumour", "Fibroblast", "Endothelial", "Macrophage", "Monocyte",
         "Neutrophil", "Dendritic", "T/NK", "B"]


def main():
    a = sc.read_h5ad(os.path.join(S.OUTPUT, "21_working", "kpca_working.h5ad"))
    a.obs["leiden"] = a.obs["leiden"].astype(str)
    uc = pd.read_csv(os.path.join(S.OUTPUT, "24_whole_umap", "24_allcell_ucell.csv"),
                     index_col=0)
    for c in ["SecA_UCell", "SecB_UCell", "polarization"]:
        a.obs[c] = uc.loc[a.obs_names, c].values
    n0 = a.n_obs
    a = a[~a.obs["leiden"].isin(DROP)].copy()
    print(f"dropped cl{sorted(DROP)}: {n0} -> {a.n_obs} cells")
    a.obs["celltype"] = pd.Categorical(a.obs["leiden"].map(CELLTYPE), categories=ORDER)
    a.obs["is_tumour"] = (a.obs["celltype"] == "Tumour").values
    print(a.obs["celltype"].value_counts().reindex(ORDER).to_string())

    a.write(os.path.join(OUT, "kpca_annotated.h5ad"))

    print("saved", OUT)


if __name__ == "__main__":
    main()
