#!/usr/bin/env python3
"""
scFEA metabolic flux on ALL epithelial cells (Fig 2a flux radar, SF8 flux block)
================================================================================
PURPOSE
    Runs scFEA (Alghamdi et al., Genome Res 2021; human M168 model) on all
    575,366 epithelial cells of the final atlas. The scFEA implementation of
    atlas/04_functional/as_run/21_epitype_functional_characterization.py does not
    scale to all cells, so this script re-implements it.

    The network, loss and hyperparameters are those of that implementation (100 epochs,
    lr 0.008, lambdas 1/1/1/0.01, torch seed 16, full-batch Adam; input = the
    epithelial h5ad X, i.e. log-normalized expression restricted to the 550 M168
    genes present in the atlas). Two implementation changes make it scale, and
    neither changes the mathematics:
      1. No module-blocked matrix. The original fed each module's encoder the full
         gene vector with non-module genes set to 0 (a cells x 168*550 dense
         matrix, ~270 GB here). Because the first layer has no bias, zeroed inputs
         contribute nothing and receive zero gradient, so each encoder is evaluated
         on its own module genes only: x[:, idx] @ W[:, idx].T. The full weight
         matrices are still created (same initialization and RNG stream), and the
         unused columns get zero gradient, so Adam leaves them unchanged exactly as
         in the original.
      2. The per-cell Pearson term and per-cell prediction loops are vectorized
         (same formula, including the +1e-8 in the denominator).
    `--validate` (optional, not part of the figure run) reruns the 20,000 cells whose
    flux the as_run implementation saved (output/21_epitype_functional_characterization/
    21_flux_per_cell.parquet) and compares per-module flux with that output.

INPUTS   data_root/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad (X)
         data_root/2026_final_atlas/tools/scFEA/data/{module_gene_m168,cmMat_c70_m168,
         cName_c70_m168}.csv;  figures/_prep/atlas_final_tables.load_epithelial()
OUTPUTS  output_root/figure2/fig2a_scfea_all/21_flux_per_cell.parquet (575,366 x 168)
         21_balance_per_cell.parquet (scFEA compound balance, written with the flux; not read
         by any figure); run_log.json (run record)
         (--validate: validate_vs_2026-05_subset.csv instead)
MANUSCRIPT PANEL(S): Fig 2a (flux radar), Supplementary Figure 8 (flux block).
RUNTIME TIER: moderate (CPU; minutes per 10 epochs on all cells).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

EPI_H5AD = path("data_root", "2026_final_atlas", "celltype_h5ad", "hgsc_atlas_final_epithelial.h5ad")
SCFEA = path("data_root", "2026_final_atlas", "tools", "scFEA", "data")
OLD_FLUX = path("data_root", "2026_final_atlas", "output", "21_epitype_functional_characterization",
                "21_flux_per_cell.parquet")
OUT_DIR = path("output_root", "figure2", "fig2a_scfea_all")
os.makedirs(OUT_DIR, exist_ok=True)
EPOCH, LR = 100, 0.008
LAMB_BA, LAMB_NG, LAMB_CELL, LAMB_MOD = 1.0, 1.0, 1.0, 1e-2
TORCH_SEED = 16

ap = argparse.ArgumentParser()
ap.add_argument("--validate", action="store_true",
                help="rerun the 20,000 cells of the 2026-05 subset and compare with its output")
ap.add_argument("--threads", type=int, default=12)
args = ap.parse_args()
torch.set_num_threads(args.threads)

# ---------- 1. Model files and genes ----------
moduleGene = pd.read_csv(os.path.join(SCFEA, "module_gene_m168.csv"), index_col=0)
cmMat = pd.read_csv(os.path.join(SCFEA, "cmMat_c70_m168.csv"), header=None).values
cName = pd.read_csv(os.path.join(SCFEA, "cName_c70_m168.csv")).columns
module_genes = [[g for g in moduleGene.iloc[i, :].values if pd.notna(g)] for i in range(moduleGene.shape[0])]
moduleLen = np.array([len(g) for g in module_genes])
all_mod_genes = {g for gs in module_genes for g in gs}

import h5py  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from scipy import sparse  # noqa: E402

with h5py.File(EPI_H5AD, "r") as f:
    var = np.asarray(read_elem(f["var"][f["var"].attrs["_index"]])).astype(str)
    obs_idx = np.asarray(read_elem(f["obs"][f["obs"].attrs["_index"]])).astype(str)
gene_names = sorted(all_mod_genes & set(var))          # same gene order as the as_run implementation
col_of = {g: i for i, g in enumerate(var)}
cols = np.array([col_of[g] for g in gene_names])
n_modules, n_genes, n_comps = len(module_genes), len(gene_names), cmMat.shape[0]
gpos = {g: i for i, g in enumerate(gene_names)}
mod_idx = [torch.tensor([gpos[g] for g in gs if g in gpos], dtype=torch.long) for gs in module_genes]
print(f"[1] M168: {n_modules} modules, {len(all_mod_genes)} genes, {n_genes} in the atlas", flush=True)

# ---------- 2. Cells ----------
epi = load_epithelial(columns=["state"])
assert (epi.index.values == obs_idx).all(), "epithelial h5ad row order != loader order"
if args.validate:
    old = pd.read_parquet(OLD_FLUX)
    cells = old.index.astype(str).values
    rows = pd.Index(obs_idx).get_indexer(cells)
    assert (rows >= 0).all()
else:
    cells = obs_idx
    rows = np.arange(len(obs_idx))
n_cells = len(cells)

# ---------- 3. Expression (X restricted to the module genes), chunked CSR read ----------
t0 = time.time()
X = np.zeros((n_cells, n_genes), dtype=np.float32)
order = np.argsort(rows)
srows = rows[order]
with h5py.File(EPI_H5AD, "r") as f:
    indptr = f["X/indptr"][:]
    CH = 50_000
    for a in range(0, len(srows), CH):
        r = srows[a:a + CH]
        lo, hi = int(r.min()), int(r.max()) + 1
        p0, p1 = indptr[lo], indptr[hi]
        blk = sparse.csr_matrix((f["X/data"][p0:p1], f["X/indices"][p0:p1], indptr[lo:hi + 1] - p0),
                                shape=(hi - lo, len(var)))
        X[order[a:a + CH]] = blk[r - lo][:, cols].toarray()
assert X.max() <= 50, "X does not look log-normalized"
print(f"[3] expression {X.shape} read in {time.time() - t0:.0f}s", flush=True)

x = torch.from_numpy(X)
geneExprSum = X.sum(axis=1)
geneExprScale = torch.from_numpy((geneExprSum / geneExprSum.mean()).astype(np.float32))
module_scale = np.column_stack([X[:, mi.numpy()].sum(axis=1) if len(mi) else np.zeros(n_cells, np.float32)
                                for mi in mod_idx]) / moduleLen
module_scale = torch.from_numpy(module_scale.astype(np.float32))
cm = torch.from_numpy(cmMat.astype(np.float32))
del X


# ---------- 4. Network (scFEA architecture; module-restricted forward) ----------
class FLUX(nn.Module):
    def __init__(self, n_modules, f_in, f_out=1):
        super().__init__()
        self.m_encoder = nn.ModuleList([
            nn.Sequential(nn.Linear(f_in, 8, bias=False), nn.Tanhshrink(),
                          nn.Linear(8, f_out), nn.Tanhshrink())
            for _ in range(n_modules)])

    def forward(self, x):
        outs = []
        for i, enc in enumerate(self.m_encoder):
            lin1, act1, lin2, act2 = enc
            h = x[:, mod_idx[i]] @ lin1.weight[:, mod_idx[i]].T   # == lin1(masked full vector)
            outs.append(act2(lin2(act1(h))))
        m = torch.cat(outs, 1)
        c = m @ cm.T                                              # == updateC loop
        return m, c


def scfea_loss(m, c, geneScale, moduleScale):
    total1 = torch.sum(c ** 2, dim=1)
    total2 = torch.sum(torch.abs(m) - m, dim=1)
    diff = (torch.sum(m, dim=1) - geneScale) ** 2
    total3 = torch.sqrt(diff) if bool((diff > 0).all()) else diff
    xm = m - m.mean(dim=1, keepdim=True)
    ym = moduleScale - moduleScale.mean(dim=1, keepdim=True)
    corr = (xm * ym).sum(1) / (torch.sqrt((xm ** 2).sum(1) * (ym ** 2).sum(1)) + 1e-8)
    total4 = 1 - torch.abs(corr)
    return (LAMB_BA * total1.sum() + LAMB_NG * total2.sum() + LAMB_CELL * total3.sum()
            + LAMB_MOD * total4.sum())


torch.manual_seed(TORCH_SEED)
net = FLUX(n_modules, f_in=n_genes)
opt = torch.optim.Adam(net.parameters(), lr=LR)
t1 = time.time()
net.train()
for ep in range(EPOCH):
    m, c = net(x)
    loss = scfea_loss(m, c, geneExprScale, module_scale)
    opt.zero_grad()
    loss.backward()
    opt.step()
    if ep % 10 == 0 or ep == EPOCH - 1:
        print(f"    epoch {ep + 1}/{EPOCH}  loss {loss.item():.4g}  ({time.time() - t1:.0f}s)", flush=True)

net.eval()
with torch.no_grad():
    parts_m, parts_c = [], []
    for a in range(0, n_cells, 100_000):
        mm, cc = net(x[a:a + 100_000])
        parts_m.append(mm.numpy()); parts_c.append(cc.numpy())
flux = pd.DataFrame(np.vstack(parts_m), index=cells, columns=moduleGene.index)
balance = pd.DataFrame(np.vstack(parts_c), index=cells, columns=cName)

# ---------- 5. Output ----------
if args.validate:
    old = old.loc[cells]
    r = [np.corrcoef(flux[k], old[k])[0, 1] if old[k].std() > 0 and flux[k].std() > 0 else np.nan
         for k in flux.columns]
    cmp_ = pd.DataFrame({"module": flux.columns, "pearson_r": r,
                         "max_abs_diff": (flux - old).abs().max().values})
    cmp_.to_csv(os.path.join(OUT_DIR, "validate_vs_2026-05_subset.csv"), index=False)
    st = epi.loc[cells, "state"].astype(str)
    zn = flux.groupby(st.values).mean(); zo = old.groupby(st.values).mean()
    print(f"[validate] per-module r: median {np.nanmedian(r):.6f}, min {np.nanmin(r):.6f}; "
          f"max |Δ| {cmp_.max_abs_diff.max():.3g}; per-state means max |Δ| {(zn - zo).abs().max().max():.3g}")
else:
    flux.to_parquet(os.path.join(OUT_DIR, "21_flux_per_cell.parquet"))
    balance.to_parquet(os.path.join(OUT_DIR, "21_balance_per_cell.parquet"))
    json.dump({"n_cells": int(n_cells), "n_modules": n_modules, "n_genes": n_genes, "epochs": EPOCH,
               "lr": LR, "torch_seed": TORCH_SEED, "final_loss": float(loss.item()),
               "train_seconds": round(time.time() - t1)}, open(os.path.join(OUT_DIR, "run_log.json"), "w"),
              indent=1)
    print(f"[5] flux {flux.shape} -> {OUT_DIR}", flush=True)
