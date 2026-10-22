#!/usr/bin/env python3
"""
Figure S10 — KPCA.B mouse tumour NMF programs are not cleanly aligned with the human SecA / SecB programs
=========================================================================================================
PURPOSE
  (a) dot plot of the 7 + 7 human SecA / SecB scoring genes and their mouse orthologs: % of human secretory
      cells vs % of KPCA.B tumour cells expressing each; (b) UMAPs of the orthologs expressed in >= 20% of
      tumour cells; (c) Spearman correlation of the gene loadings of the four KPCA.B NMF programs with the human
      NMF factors and the SecB - SecA axis (no mouse program specifically matches Factor 2 = SecB or Factor 3 =
      SecA; ns = FDR > 0.05); (d) per-cell usage of each KPCA.B NMF program on the tumour UMAP.

Mouse cells: the 12,429-cell `celltype == "Tumour"` manifold of the final KPCA.B object (the epithelial set
of main Fig 3i-k; NOT the `compartment` column, which mislabels tumour via Col1a1 ambient). This set is
identical (same barcodes, same log1p-CP10k X) to the cells the final k = 4 NMF (mouse step 31, 2026-09-09)
was fitted on — asserted below. Every cell is used; no subsampling.

INPUTS (read-only)
  Mouse project 202608_kpca_secb_mouse/output/
    33_handoff/kpca_final_annotated.h5ad      (2026-09-10) celltype, epi_umap_1/2 (Fig 3i-k embedding), X
    31_fulltumour_nmf/31_program_usage{,_raw}.csv, 31_program_loadings.csv, 31_program_annotation.csv,
                       31_program_summary.csv, kpca_fulltumour_programs.h5ad (cell list check)   (2026-09-09)
    32_secb_align/32_alignment_stats.csv, 32_verdict.csv                                        (2026-09-09)
    02_orthologs/ortholog_map.csv (MGI human-mouse homology), 02_coverage_summary.csv          (2026-08-26)
  Human
    202601_final_atlas/output/11d_epithelial_nmf/11d_nmf_loadings.csv (10 factors x 3,007 genes)
    shared/signatures.yml (7 + 7 scoring genes)
    final_panels/supplementary/tables/SuppTable6_..._{Scoring_signatures,Extended_lists}.csv
        (177/177 SecA-/SecB-specific lists; per-gene % expressing in atlas SecA / SecB cells)

METHOD NOTES
  * Mouse loadings over ALL genes are re-derived exactly as mouse step 32 (`extend_all_genes`: NNLS per gene
    with the step-31 raw usage W fixed, rows L2-normalised), then restricted to the shared 1:1-preferred
    ortholog space of the human 3,007 NMF genes (step 05 `align_gene_space`). The script asserts that its
    Spearman matrix reproduces 32_alignment_stats.csv; panel c plots the 32 values (final, Sep 9), with the
    step-32 permutation FDR (BH over the 44 tests).

OUTPUTS
  final_panels/supplementary/figures/SuppFig10{a,b,c,d}_*.{svg,png (450 dpi)}
  final_panels/supplementary/tables/SuppFig10_kpca_{ortholog_expression,program_alignment}.csv (source data of a, c)
  (set SF10_OUT=<dir> to render elsewhere, e.g. a scratch dir)

MANUSCRIPT PANEL(S)
  Figure S10a-d.

RUNTIME ~2-3 min (NNLS over 19,377 genes x 12,429 cells).
"""
import os
import numpy as np
import pandas as pd
import yaml
import anndata as ad
import scipy.sparse as sp
from scipy.optimize import nnls
from scipy.stats import spearmanr, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# ============================================================================ paths
MOUSE = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "202608_kpca_secb_mouse/output")
H5AD = f"{MOUSE}/33_handoff/kpca_final_annotated.h5ad"
D31, D32, D02 = f"{MOUSE}/31_fulltumour_nmf", f"{MOUSE}/32_secb_align", f"{MOUSE}/02_orthologs"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
H_LOAD = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "2026_final_atlas/output/11d_epithelial_nmf/11d_nmf_loadings.csv")
SIGS = f"{REPO}/shared/signatures.yml"
FIN = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "final_panels/supplementary")
ST6 = f"{FIN}/tables/SuppTable6_secA_secB_scoring_signatures"
OUT = os.environ.get("SF10_OUT", f"{FIN}/figures")
TAB = os.environ.get("SF10_OUT", f"{FIN}/tables")
os.makedirs(OUT, exist_ok=True); os.makedirs(TAB, exist_ok=True)
STEMS = {"a": "SuppFig10a_kpca_ortholog_expression", "b": "SuppFig10b_kpca_expressed_ortholog_umaps",
         "c": "SuppFig10c_kpca_vs_human_factor_alignment", "d": "SuppFig10d_kpca_nmf_program_umaps"}   # Figure S10: a dot plot, b expressed-ortholog UMAPs, c heatmap, d NMF program UMAPs
EXPR_MIN_PCT = 20.0   # panel B shows scoring-gene orthologs detected in >= 20% of KPCA.B tumour cells
SEED = 0          # shuffle seed of the Fig 3i-k export (41b_export_1p8in.py)

# ============================================================================ style
FA, FK, FS = 6.5, 5.5, 5
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FK, "axes.labelsize": FK, "axes.titlesize": FA, "xtick.labelsize": FK, "ytick.labelsize": FK,
    "legend.fontsize": FS, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "svg.fonttype": "none", "pdf.fonttype": 42,
    "figure.dpi": 450, "savefig.dpi": 450,
})
SECA, SECB, GREY = "#E6A141", "#9A7D55", "#8C8C8C"
USAGE_CMAP = LinearSegmentedColormap.from_list(      # SF5 (human NMF usage) colourmap
    "beige_burgundy", [(0.00, "#DDD5CA"), (0.48, "#DDD5CA"), (0.65, "#C8878A"), (0.85, "#A04858"), (1.00, "#6B2838")], N=256)
EXPR_CMAP = LinearSegmentedColormap.from_list("grey_seq", ["#EDEDED", "#9E9E9E", "#262626"], N=256)
RHO_CMAP = LinearSegmentedColormap.from_list(
    "blue_white_burgundy", [(0.0, "#3E6A93"), (0.5, "#FFFFFF"), (0.75, "#C8878A"), (1.0, "#6B2838")], N=256)

FACTOR_SHORT = {"Factor_1": "F1 Housekeeping", "Factor_2": "F2 SecB", "Factor_3": "F3 SecA", "Factor_4": "F4 Interferon",
                "Factor_5": "F5 Cell cycle", "Factor_6": "F6 AP-1 / stress", "Factor_7": "F7 Progenitor B",
                "Factor_8": "F8 Secretory mat.", "Factor_9": "F9 EMT", "Factor_10": "F10 Metallothionein",
                "Polarization_SecB_minus_SecA": "SecB − SecA axis"}
PROGS = ["mP1", "mP2", "mP3", "mP4"]
# display annotation of the step-31 top-loaded genes (31_program_annotation.csv); mP1 = ambient-flagged there
PROG_NAME = {"mP1": "ECM / collagen", "mP2": "Secretory (SecB-like)", "mP3": "ECM / mesenchymal TF",
             "mP4": "ECM / wound response"}
PROG_GENES = {"mP1": "Col1a1, Col3a1, Col1a2, Fn1", "mP2": "Clu, Lcn2, Ltf, Krt7, Cldn4",
              "mP3": "Fn1, Ccn2, Nr2f1, Meis2, Esr1", "mP4": "Dst, Col7a1, Serpine1, Plau"}


def save(fig, key, tight=True):
    for ext in ("svg", "png"):
        fig.savefig(f"{OUT}/{STEMS[key]}.{ext}", format=ext, dpi=450, bbox_inches="tight" if tight else None,
                    pad_inches=0.02)
    plt.close(fig)
    print(f"[save] {STEMS[key]}.svg/.png")


def dense_col(X, j):
    c = X[:, j]
    return (c.toarray() if sp.issparse(c) else np.asarray(c)).ravel()


# ============================================================================ load
print("[load] mouse final object")
A = ad.read_h5ad(H5AD)
tum = (A.obs["celltype"].astype(str) == "Tumour").to_numpy()
T = A[tum].copy()
print(f"  celltype == Tumour: {T.n_obs:,} cells x {T.n_vars:,} genes")
assert T.n_obs == 12429
assert T.obs["epi_umap_1"].notna().all() and A.obs.loc[~tum, "epi_umap_1"].isna().all()   # = Fig 3i-k embedding set

usage = pd.read_csv(f"{D31}/31_program_usage.csv", index_col=0)[PROGS]
usage_raw = pd.read_csv(f"{D31}/31_program_usage_raw.csv", index_col=0)[PROGS]
ann = pd.read_csv(f"{D31}/31_program_annotation.csv").set_index("program")
summ31 = pd.read_csv(f"{D31}/31_program_summary.csv").set_index("program")
fit_cells = ad.read_h5ad(f"{D31}/kpca_fulltumour_programs.h5ad", backed="r").obs_names
assert set(fit_cells) == set(T.obs_names) == set(usage.index), "NMF fit cells != 12,429 Tumour cells"
print(f"  step-31 NMF fit cells == Tumour cells: True ({len(fit_cells):,}); k = {len(PROGS)}, "
      f"recurrence {summ31['recurrence'].min():.2f}, {int(summ31['n_seeds'].min())} seeds")
T = T[usage.index].copy()
assert np.allclose(T.obs[[f"NMF_{p}" for p in PROGS]].to_numpy(float), usage.to_numpy(float))

ortho = pd.read_csv(f"{D02}/ortholog_map.csv")
LUT = dict(zip(ortho["human_symbol"], ortho["mouse_symbol"]))
CLS = dict(zip(ortho["human_symbol"], ortho["ortholog_class"]))
mouse_genes = set(T.var_names)


def to_mouse(g):
    """step 00 `map_genes`: MGI first, case-conversion fallback."""
    m = LUT.get(g)
    if not isinstance(m, str):
        m = g[0].upper() + g[1:].lower()
    return m


sig = yaml.safe_load(open(SIGS))
SECA_7, SECB_7 = sig["SecA"], sig["SecB"]
st6 = pd.read_csv(f"{ST6}_Scoring_signatures.csv").set_index("gene")
ext = pd.read_csv(f"{ST6}_Extended_lists.csv")
LISTS = {c: ext.loc[ext["list"] == c, "gene"].tolist() for c in ["SecA-specific", "SecB-specific"]}
HL = pd.read_csv(H_LOAD, index_col=0)

# ============================================================================ mouse expression stats
X = T.X.tocsc() if sp.issparse(T.X) else T.X
vidx = {g: i for i, g in enumerate(T.var_names)}


def mstats(m):
    if m not in vidx:
        return np.nan, np.nan
    v = dense_col(X, vidx[m])
    return 100 * float((v > 0).mean()), float(v.mean())


rows = []
for cls, genes in [("SecA", SECA_7), ("SecB", SECB_7)]:
    for g in genes:
        m = to_mouse(g); pct, mu = mstats(m)
        rows.append({"set": "scoring_7", "class": cls, "human_gene": g, "mouse_ortholog": m,
                     "ortholog_class": CLS.get(g, "case-conversion"), "on_mouse_panel": m in mouse_genes,
                     "human_pct_SecA_cells": st6.loc[g, "pct_expr_atlas_SecA"],
                     "human_pct_SecB_cells": st6.loc[g, "pct_expr_atlas_SecB"],
                     "mouse_pct_expressing": pct, "mouse_mean_log1p_cp10k": mu})
for cls, genes in LISTS.items():
    e = ext.set_index("gene")
    for g in genes:
        m = to_mouse(g); pct, mu = mstats(m)
        rows.append({"set": "specific_177", "class": cls.split("-")[0], "human_gene": g, "mouse_ortholog": m,
                     "ortholog_class": CLS.get(g, "case-conversion"), "on_mouse_panel": m in mouse_genes,
                     "human_pct_SecA_cells": e.loc[g, "pct_expr_atlas_SecA"],
                     "human_pct_SecB_cells": e.loc[g, "pct_expr_atlas_SecB"],
                     "mouse_pct_expressing": pct, "mouse_mean_log1p_cp10k": mu})
EX = pd.DataFrame(rows)
cov = EX[EX["set"] == "specific_177"].groupby("class")["on_mouse_panel"].agg(["sum", "size"])
cov_ref = pd.read_csv(f"{D02}/02_coverage_summary.csv").set_index("signature")
assert cov.loc["SecA", "sum"] == cov_ref.loc["SecA_specific_177", "mapped"] == 171
assert cov.loc["SecB", "sum"] == cov_ref.loc["SecB_specific_177", "mapped"] == 172
print(f"  177-gene coverage on the mouse panel: SecA {cov.loc['SecA','sum']}/177, SecB {cov.loc['SecB','sum']}/177")
EX.to_csv(f"{TAB}/SuppFig10_kpca_ortholog_expression.csv", index=False)

# ============================================================================ extended loadings (= step 32)
print("[nnls] extending mouse loadings to all genes (W fixed) ...", flush=True)
W = usage_raw.loc[T.obs_names].to_numpy(float)
Hx = np.zeros((W.shape[1], T.n_vars))
for j in range(T.n_vars):
    Hx[:, j], _ = nnls(W, np.clip(dense_col(X, j), 0, None))
ML = pd.DataFrame(Hx, index=PROGS, columns=T.var_names)
ML = ML.div(np.linalg.norm(ML.to_numpy(), axis=1, keepdims=True) + 1e-12)

# shared ortholog space (= step 05 align_gene_space)
seen, hk, mk = set(), [], []
for g in HL.columns:
    m = to_mouse(g)
    if m in mouse_genes and m not in seen:
        seen.add(m); hk.append(g); mk.append(m)
Hm = HL[hk].copy(); Mm = ML[mk].copy(); Mm.columns = hk          # both indexed by human symbol
Hm.loc["Polarization_SecB_minus_SecA"] = Hm.loc["Factor_2"] - Hm.loc["Factor_3"]
print(f"  shared orthologs: {len(hk):,}")
FACS = [f"Factor_{i}" for i in range(1, 11)] + ["Polarization_SecB_minus_SecA"]
rho = pd.DataFrame([[spearmanr(Mm.loc[p], Hm.loc[f]).correlation for f in FACS] for p in PROGS], index=PROGS, columns=FACS)

st32 = pd.read_csv(f"{D32}/32_alignment_stats.csv")
R32 = st32.pivot(index="mouse_program", columns="human_factor", values="spearman").loc[PROGS, FACS]
P32 = st32.pivot(index="mouse_program", columns="human_factor", values="p_adj").loc[PROGS, FACS]
maxdiff = float(np.abs(rho.to_numpy() - R32.to_numpy()).max())
print(f"  reproduces 32_alignment_stats.csv: max |d rho| = {maxdiff:.2e}")
assert maxdiff < 1e-6
V = pd.read_csv(f"{D32}/32_verdict.csv").iloc[0]
FLOOR, CEIL, FRAC, GAP = float(V["floor"]), float(V["ceiling"]), float(V["frac_floor_to_ceiling"]), float(V["gap_to_second"])

al = st32.copy(); al["verdict_floor"] = FLOOR; al["verdict_ceiling"] = CEIL
al.to_csv(f"{TAB}/SuppFig10_kpca_program_alignment.csv", index=False)

# ============================================================================ PANEL A — human secretory vs mouse tumour (one dot plot)
print("[A]")
S7 = EX[EX["set"] == "scoring_7"].reset_index(drop=True)
# human: all secretory epithelial cells pooled (cell-weighted SecA / Intermediate / SecB; Supp Table 6 state sizes)
N_STATE = {"SecA": 284344, "Intermediate": 142172, "SecB": 142172}
hp = st6.loc[S7["human_gene"], ["pct_expr_atlas_SecA", "pct_expr_atlas_Intermediate", "pct_expr_atlas_SecB"]].to_numpy(float)
w = np.array([N_STATE["SecA"], N_STATE["Intermediate"], N_STATE["SecB"]], float)
S7["human_pct_secretory"] = (hp * w).sum(1) / w.sum()
EX = EX.merge(S7[["human_gene", "human_pct_secretory"]], on="human_gene", how="left")
EX.to_csv(f"{TAB}/SuppFig10_kpca_ortholog_expression.csv", index=False)
n = len(S7)
fig = plt.figure(figsize=(2.3, 3.05))
axd = fig.add_axes([0.50, 0.15, 0.36, 0.72])
COLX = [0, 1.35]
SMAX = 70.0                                   # marker area (pt^2) at 100%
yy = np.arange(n)[::-1].astype(float); yy[7:] -= 0.6       # gap between SecA and SecB blocks
for i, r in S7.iterrows():
    c = SECA if r["class"] == "SecA" else SECB
    for x, v in [(COLX[0], r.human_pct_secretory), (COLX[1], r.mouse_pct_expressing)]:
        axd.scatter(x, yy[i], s=SMAX * v / 100, c=c, edgecolors="#333333", linewidths=0.3, zorder=3)
        axd.text(x + 0.30, yy[i], f"{v:.0f}%", va="center", ha="left", fontsize=FS, color="#333333")
axd.set_xlim(-0.45, 2.15); axd.set_ylim(yy.min() - 0.7, yy.max() + 0.7)
axd.set_xticks(COLX); axd.set_xticklabels(["Human\nsecretory", "KPCA.B\ntumour"], fontsize=FK)
axd.xaxis.tick_top()
axd.set_yticks(yy); axd.set_yticklabels([])
for sname in ("left", "right", "bottom", "top"):
    axd.spines[sname].set_visible(False)
axd.tick_params(length=0)
for i, r in S7.iterrows():
    c = SECA if r["class"] == "SecA" else SECB
    axd.text(-2.75, yy[i], r.human_gene, ha="left", va="center", fontsize=FK, fontstyle="italic", color=c, fontweight="bold")
    axd.text(-0.55, yy[i], r.mouse_ortholog, ha="right", va="center", fontsize=FK, fontstyle="italic", color="#333333")
for cls, sl, c in [("SecA", slice(0, 7), SECA), ("SecB", slice(7, 14), SECB)]:
    y0, y1 = yy[sl].min(), yy[sl].max()
    axd.plot([-2.92, -2.92], [y0 - 0.25, y1 + 0.25], color=c, lw=1.2, clip_on=False)
    axd.text(-3.05, (y0 + y1) / 2, cls, rotation=90, ha="right", va="center", fontsize=FA, fontweight="bold", color=c)
lg = fig.add_axes([0.05, 0.02, 0.85, 0.08]); lg.axis("off")
lg.text(0.0, 0.5, "% cells expressing", fontsize=FS, va="center")
for k, pv in enumerate([10, 40, 80]):
    lg.scatter(0.47 + 0.18 * k, 0.5, s=SMAX * pv / 100, c="#BBBBBB", edgecolors="#333333", linewidths=0.3)
    lg.text(0.47 + 0.18 * k + 0.05, 0.5, f"{pv}%", va="center", fontsize=FS)
lg.set_xlim(0, 1); lg.set_ylim(0, 1)
save(fig, "a")

# ============================================================================ PANEL B — expressed ortholog UMAPs
print("[B]")
ex_ = T.obs["epi_umap_1"].to_numpy(float); ey_ = T.obs["epi_umap_2"].to_numpy(float)
xr, yr = ex_.max() - ex_.min(), ey_.max() - ey_.min(); px, py = 0.04 * xr, 0.04 * yr
GUT = 0.30 * xr   # empty strip right of the cells so colour bars never sit on the UMAP
# same draw order as Fig 3i-k: permutation of the epithelial subset in the final-object order; expressing cells drawn on top
order_names = A.obs_names[tum]
perm = np.random.default_rng(SEED).permutation(len(order_names))
pos = pd.Series(np.arange(T.n_obs), index=T.obs_names).loc[order_names[perm]].to_numpy()
EXPd = S7[S7["mouse_pct_expressing"] >= EXPR_MIN_PCT]
rowsB = [EXPd[EXPd["class"] == "SecA"], EXPd[EXPd["class"] == "SecB"]]
ncol = max(len(r) for r in rowsB)
print("  expressed (>= %.0f%%): SecA %s | SecB %s" % (EXPR_MIN_PCT, list(rowsB[0].mouse_ortholog), list(rowsB[1].mouse_ortholog)))
SIDE = 1.15
fig, axes = plt.subplots(2, ncol, figsize=(ncol * SIDE + 0.22, 2 * SIDE + 0.36), squeeze=False)
fig.subplots_adjust(left=0.22 / (ncol * SIDE + 0.22), right=0.995, top=1 - 0.20 / (2 * SIDE + 0.36), bottom=0.005,
                    wspace=0.04, hspace=0.14)
for ri, (R, c) in enumerate(zip(rowsB, [SECA, SECB])):
    cmap = LinearSegmentedColormap.from_list(f"g{ri}", ["#E4E0DA", c, "#3A2A1A" if ri else "#8A5A10"], N=256)
    for ci in range(ncol):
        ax = axes[ri, ci]; ax.set_xticks([]); ax.set_yticks([])
        if ci >= len(R):
            ax.axis("off"); continue
        r = R.iloc[ci]; v = dense_col(X, vidx[r.mouse_ortholog])
        vmax = float(np.quantile(v[v > 0], 0.95)) if (v > 0).any() else 1.0; norm = Normalize(0, vmax)
        o = pos[np.argsort(v[pos] > 0, kind="stable")]                       # zeros first, expressing on top
        ax.scatter(ex_[o], ey_[o], c=v[o], cmap=cmap, norm=norm, s=0.6, linewidths=0, rasterized=True)
        ax.set_xlim(ex_.min() - px, ex_.max() + px + GUT); ax.set_ylim(ey_.min() - py, ey_.max() + py); ax.set_aspect("equal")   # GUT = blank strip for the colour bar
        for sname in ("top", "right"):
            ax.spines[sname].set_visible(False)
        ax.spines["bottom"].set_bounds(ex_.min() - px, ex_.min() - px + 0.18 * xr)
        ax.spines["left"].set_bounds(ey_.min() - py, ey_.min() - py + 0.18 * yr)
        ax.spines["bottom"].set_linewidth(0.4); ax.spines["left"].set_linewidth(0.4)
        ax.set_title(f"$\\bf{{\\it{{{r.mouse_ortholog}}}}}$  {r.mouse_pct_expressing:.0f}%", fontsize=FA, color=c, loc="left", pad=2)
        cax = ax.inset_axes([0.83, 0.06, 0.04, 0.30])
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="vertical")
        cb.outline.set_linewidth(0.3); cax.yaxis.set_ticks_position("right")
        cb.ax.tick_params(length=1.5, width=0.3, pad=1, labelsize=FS - 0.5)
        cb.set_ticks([0, vmax]); cb.set_ticklabels(["0", f"{vmax:.1f}"])
save(fig, "b", tight=False)

# ============================================================================ PANEL C — alignment heatmap
print("[C]")
fig = plt.figure(figsize=(3.55, 1.7))
ax = fig.add_axes([0.23, 0.05, 0.66, 0.56])
M = R32.to_numpy(float); VA = 0.45
ax.imshow(M, cmap=RHO_CMAP, vmin=-VA, vmax=VA, aspect="auto")
colbest = M.argmax(axis=0)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        bold = (j in (1, 2)) and colbest[j] == i
        txt = f"{M[i, j]:.2f}" + ("" if P32.iloc[i, j] < 0.05 else "ⁿˢ")
        ax.text(j, i, txt, ha="center", va="center", fontsize=4.6,
                fontweight="bold" if bold else "normal", color="white" if abs(M[i, j]) > 0.33 else "#222222")
for j, c in [(1, SECB), (2, SECA)]:
    ax.add_patch(Rectangle((j - 0.5, -0.5), 1, M.shape[0], fill=False, ec=c, lw=1.3, zorder=5, clip_on=False))
ax.axvline(9.5, color="white", lw=1.5)
ax.set_xticks(range(len(FACS))); ax.set_xticklabels([FACTOR_SHORT[f] for f in FACS], rotation=50, ha="left", fontsize=FS)
ax.xaxis.tick_top(); ax.tick_params(length=0)
for t, f in zip(ax.get_xticklabels(), FACS):
    if f in ("Factor_2", "Factor_3"):
        t.set_color(SECB if f == "Factor_2" else SECA); t.set_fontweight("bold")
ax.set_yticks(range(len(PROGS)))
ax.set_yticklabels([f"{p} {PROG_NAME[p].split(' (')[0]}" for p in PROGS], fontsize=FS)
for s in ax.spines.values():
    s.set_visible(False)
cax = fig.add_axes([0.905, 0.12, 0.018, 0.42])
cb = fig.colorbar(ScalarMappable(norm=Normalize(-VA, VA), cmap=RHO_CMAP), cax=cax)
cb.outline.set_linewidth(0.3); cb.set_ticks([-0.4, 0, 0.4]); cb.ax.tick_params(length=1.5, width=0.3, pad=1, labelsize=FS)
cb.set_label("Spearman ρ (loadings)", fontsize=FS, labelpad=2)
save(fig, "c")

# ============================================================================ PANEL D — NMF program usage UMAPs
print("[D]")
SIDE = 1.15
fig, axes = plt.subplots(1, 4, figsize=(4 * SIDE + 0.22, SIDE + 0.22))   # one row (author 2026-10-01)
# fixed canvas (no tight crop): 0.12 in left / 0.06 in top margin leave room for the composite panel letter
fig.subplots_adjust(left=0.22 / (4 * SIDE + 0.22), right=0.995, top=1 - 0.20 / (SIDE + 0.22), bottom=0.01,
                    wspace=0.04)
for ax, p in zip(axes.ravel(), PROGS):
    v = usage[p].to_numpy(float)
    vmax = float(np.quantile(v, 0.99)); norm = Normalize(0, vmax)
    ax.scatter(ex_[pos], ey_[pos], c=v[pos], cmap=USAGE_CMAP, norm=norm, s=0.6, marker="o", linewidths=0,
               rasterized=True, alpha=0.95)
    ax.set_xlim(ex_.min() - px, ex_.max() + px + GUT); ax.set_ylim(ey_.min() - py, ey_.max() + py); ax.set_aspect("equal")   # GUT = blank strip for the colour bar
    ax.set_xticks([]); ax.set_yticks([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_bounds(ex_.min() - px, ex_.min() - px + 0.18 * xr)
    ax.spines["left"].set_bounds(ey_.min() - py, ey_.min() - py + 0.18 * yr)
    ax.spines["bottom"].set_linewidth(0.4); ax.spines["left"].set_linewidth(0.4)
    ax.set_title(f"{p}  {PROG_NAME[p]}" + ("*" if bool(ann.loc[p, "flag_ambient"]) else ""), fontsize=FK, fontweight="bold", loc="left", pad=2)
    cax = ax.inset_axes([0.83, 0.06, 0.04, 0.30])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=USAGE_CMAP), cax=cax, orientation="vertical")
    cb.outline.set_linewidth(0.3); cax.yaxis.set_ticks_position("right")
    cb.ax.tick_params(length=1.5, width=0.3, pad=1, labelsize=FS)
    cb.set_ticks([0, vmax]); cb.set_ticklabels(["0.00", f"{vmax:.2f}"])
save(fig, "d", tight=False)

print("Done.")
