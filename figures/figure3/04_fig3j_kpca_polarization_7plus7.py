#!/usr/bin/env python3
"""
Figure 3j -- KPCA.B mouse secretory polarization rescored with the 7 + 7 human signature orthologues
=====================================================================================================
PURPOSE
    Scores the KPCA.B tumour epithelium with the mouse orthologues of the 7 + 7 human scoring signature
    (shared/signatures.yml) and draws Fig 3j. The hand-off object also stores an earlier 5 + 4-gene score
    (SecA_UCell_epi / SecB_UCell_epi / polarization_epi: Bcam, Mecom, Pbx1, Sox17, Rcn2 vs Krt19, Krt7, Lcn2,
    Prss22). It is used here only to validate the UCell implementation and for the agreement rows of the stats
    table, which Table S6 quotes (tables/06_tableS6_secA_secB_signatures.py).

    * Genes: shared/signatures.yml SecA / SecB (7 + 7) -> MGI 1:1 orthologues (02_orthologs/ortholog_map.csv; all 14
      are 1:1). Every orthologue is scored, however sparsely detected (detection is reported).
    * Score: UCell (Andreatta & Carmona 2021), as implemented in the KPCA.B project (scripts/24_secb_whole_umap.py),
      which made the stored scores: per cell, all genes of the log-normalised matrix ranked by expression (descending;
      ties kept in gene order, a stable sort), ranks above maxRank = 1500 set to 1501; score = 1 - (sum of signature
      ranks - n(n+1)/2) / (n x maxRank). Polarization = SecB - SecA.
    * Cells: the 12,429 tumour (epithelial) cells of Fig 3j,k (celltype == "Tumour", the cells with the
      dedicated epithelial UMAP epi_umap_1/2), ranked over all genes of the object.
    * Validation (asserted): this implementation reproduces the stored SecA_UCell_epi and SecB_UCell_epi exactly from
      their own 5 + 4 genes, so the only difference from the stored score is the gene set.
    * Panel: 1.8 in UMAP square as for Fig 3i (models/kpca_mouse/scripts/41b_export_1p8in.py), gold -> tan -> brown
      diverging map symmetric about 0 at the 99th percentile of |polarization|, polarized cells on top (seed 0),
      colour bar under the UMAP; the layout is shared with Fig 3k in figures/figure3/_kpca_umap_panel.py.

INPUTS (read-only)
    $DATA_ROOT/202608_kpca_secb_mouse/output/33_handoff/kpca_final_annotated.h5ad   ($KPCA_H5AD)
    $DATA_ROOT/202608_kpca_secb_mouse/output/02_orthologs/ortholog_map.csv           ($KPCA_ORTHO)
    shared/signatures.yml

OUTPUTS
    output_root/figure3/kpca_polarization_7plus7_cells.csv   per cell: SecA/SecB UCell (7 + 7), polarization, and
                                                             the stored polarization_epi for comparison
    output_root/figure3/kpca_polarization_7plus7_stats.csv   gene detection; agreement with the stored 5 + 4 score
    figures_dir/figure3/Fig3_j_mouse_polarization_7plus7.{svg,png}

MANUSCRIPT PANEL(S): Fig 3j.
RUNTIME TIER: fast (~1 min; reads the 20,598-cell object).
USAGE (repo root): OUTPUT_ROOT=... python figures/figure3/04_fig3j_kpca_polarization_7plus7.py
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
import matplotlib
matplotlib.use("Agg")
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import pearsonr, spearmanr

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO)
from config.config import path  # noqa: E402

MOUSE = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "202608_kpca_secb_mouse/output")
H5AD = os.environ.get("KPCA_H5AD", f"{MOUSE}/33_handoff/kpca_final_annotated.h5ad")
ORTHO = os.environ.get("KPCA_ORTHO", f"{MOUSE}/02_orthologs/ortholog_map.csv")
SIGS = os.path.join(REPO, "shared", "signatures.yml")
OUT = path("output_root", "figure3")
FIG = path("figures_dir", "figure3")
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
MAX_RANK = 1500
OLD_SECA = ["Bcam", "Mecom", "Pbx1", "Sox17", "Rcn2"]          # genes of the stored SecA_UCell_epi
OLD_SECB = ["Krt19", "Krt7", "Lcn2", "Prss22"]                 # genes of the stored SecB_UCell_epi


def ucell(X, gene_idx, max_rank=MAX_RANK, step=1000):
    """UCell scores of one gene set for every row of the cells x genes matrix X (KPCA.B project convention:
    descending stable sort over all genes, so ties keep gene order; ranks above max_rank -> max_rank + 1)."""
    n = len(gene_idx)
    out = np.empty(X.shape[0])
    for s in range(0, X.shape[0], step):
        blk = X[s:s + step]
        blk = blk.toarray() if sp.issparse(blk) else np.asarray(blk)
        order = np.argsort(-blk, axis=1, kind="stable")
        ranks = np.empty_like(order)
        ranks[np.arange(blk.shape[0])[:, None], order] = np.arange(1, blk.shape[1] + 1)[None, :]
        sig = ranks[:, gene_idx].astype(np.float64)
        sig[sig > max_rank] = max_rank + 1
        out[s:s + step] = 1 - (sig.sum(1) - n * (n + 1) / 2) / (n * max_rank)
    return out


# ---------------------------------------------------------------- load
a = ad.read_h5ad(H5AD)
tum = a.obs["epi_umap_1"].notna().values
assert tum.sum() == 12_429 and (a.obs.loc[tum, "celltype"].astype(str) == "Tumour").all()
t = a[tum]
X = t.X   # log-normalised (uns["log1p"]), as for the stored scores
genes = pd.Index(t.var_names)

# ---------------------------------------------------------------- validate against the stored scores
for gs, col in [(OLD_SECA, "SecA_UCell_epi"), (OLD_SECB, "SecB_UCell_epi")]:
    s = ucell(X, genes.get_indexer(gs))
    d = float(np.abs(s - t.obs[col].to_numpy(float)).max())
    print(f"[validate] {col} from {', '.join(gs)}: max |diff| = {d:.2e}")
    assert d < 1e-6, f"UCell implementation does not reproduce {col}"

# ---------------------------------------------------------------- 7 + 7 orthologues
sig = yaml.safe_load(open(SIGS))
om = pd.read_csv(ORTHO)
om = om[om["human_symbol"].isin(sig["SecA"] + sig["SecB"])].set_index("human_symbol")
assert len(om) == 14 and (om["ortholog_class"] == "1:1").all(), "not all 7 + 7 genes have a 1:1 orthologue"
SECA = [om.at[h, "mouse_symbol"] for h in sig["SecA"]]
SECB = [om.at[h, "mouse_symbol"] for h in sig["SecB"]]
assert genes.get_indexer(SECA + SECB).min() >= 0, "orthologue missing from the object"
det = np.asarray((X[:, genes.get_indexer(SECA + SECB)] > 0).mean(0)).ravel() * 100
seca = ucell(X, genes.get_indexer(SECA))
secb = ucell(X, genes.get_indexer(SECB))
pol = secb - seca
old = t.obs["polarization_epi"].to_numpy(float)
rs, rp = spearmanr(pol, old).correlation, pearsonr(pol, old)[0]
flip = float(np.mean(np.sign(pol) != np.sign(old)) * 100)
print(f"[7+7] SecA {SECA}\n      SecB {SECB}")
print(f"[7+7] polarization median {np.median(pol):.3f} (stored {np.median(old):.3f}); "
      f"Spearman {rs:.3f}, Pearson {rp:.3f} with the stored score; sign differs in {flip:.1f}% of cells")

cells = pd.DataFrame({"cell": t.obs_names, "SecA_UCell_7": seca, "SecB_UCell_7": secb,
                      "polarization_7plus7": pol, "polarization_epi_stored": old})
cells.to_csv(f"{OUT}/kpca_polarization_7plus7_cells.csv", index=False)
stats = pd.DataFrame(
    [{"item": f"detected (% of tumour cells): {g} ({'SecA' if g in SECA else 'SecB'})", "value": round(v, 1)}
     for g, v in zip(SECA + SECB, det)]
    + [{"item": "cells", "value": int(tum.sum())},
       {"item": "polarization median (7 + 7)", "value": round(float(np.median(pol)), 4)},
       {"item": "polarization median (stored 5 + 4)", "value": round(float(np.median(old)), 4)},
       {"item": "Spearman r, 7 + 7 vs stored", "value": round(float(rs), 4)},
       {"item": "Pearson r, 7 + 7 vs stored", "value": round(float(rp), 4)},
       {"item": "% cells with a different sign", "value": round(flip, 2)},
       {"item": "SecA 7 vs stored SecA Spearman", "value": round(float(spearmanr(seca, t.obs["SecA_UCell_epi"]).correlation), 4)},
       {"item": "SecB 7 vs stored SecB Spearman", "value": round(float(spearmanr(secb, t.obs["SecB_UCell_epi"]).correlation), 4)}])
stats.to_csv(f"{OUT}/kpca_polarization_7plus7_stats.csv", index=False)
print(stats.to_string(index=False))

# ---------------------------------------------------------------- panel (Fig 3i square; colour bar under the UMAP)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _kpca_umap_panel import epi_umap_panel  # noqa: E402

SEED = 0
DIV_CMAP = LinearSegmentedColormap.from_list("seca_secb_gold_brown",
                                             [(0.00, "#E6A141"), (0.50, "#C08E48"), (1.00, "#5A3E2B")], N=256)
ex = t.obs["epi_umap_1"].to_numpy(float); ey = t.obs["epi_umap_2"].to_numpy(float)
eo = np.random.default_rng(SEED).permutation(t.n_obs)
vabs = float(np.quantile(np.abs(pol), 0.99))
stem = os.path.join(FIG, "Fig3_j_mouse_polarization_7plus7")
epi_umap_panel(ex, ey, pol, eo[np.argsort(np.abs(pol[eo]))], DIV_CMAP, Normalize(vmin=-vabs, vmax=vabs),
               "Polarization", "(SecB - SecA)", [-vabs, 0.0, vabs], [f"-{vabs:.2f}", "0", f"+{vabs:.2f}"], stem,
               end_labels=("SecA", "SecB"))   # polarized cells drawn on top
print(f"[save] {stem}.svg/.png (colour limits +/-{vabs:.3f})")
