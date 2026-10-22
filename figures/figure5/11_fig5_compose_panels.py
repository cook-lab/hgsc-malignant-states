#!/usr/bin/env python3
"""
Figure 5 — subpanels sized from a full-page layout (letter, 0.5 in margins).

PURPOSE
    Every panel is drawn natively on ONE inch-based page grid so text sizes are true at
    print (labels 7 pt, ticks 6 pt) and panels align exactly. The page itself is not saved:
    each subpanel is cut from it as its own file (one component group at a time, cropped).
    Also runs the per-tissue tests of 5f and 5i and the paired tests of 5b and 5j.

INPUTS (output_root/figures/figure5/ and output_root/split_recompute/)
    gam_curves_5a.csv (02); fig5b_morphometrics_presplit.csv (03); roi_data/ (05);
    fig5c_gam_persample.csv (07), fig5c_tme_percell_presplit.csv (06); fig5e_mac_gam_curves.csv (09),
    fig5e_macrophage_percell_presplit.csv (08); fig5e_mac_morphology_unitmedians_presplit.csv (10);
    config/study_codes.py (title of panel e)

OUTPUTS (output_root/figures/figure5/)
    panels/<file>.{svg,png} (ROI maps PNG only); Figure5_composite_stats_eh.csv (per-tissue tests;
    rows "e" = panel f, rows "h" = panel i). File and component names keep earlier letters:
      a  Fig5a_epithelial_programs                 epithelial program GAMs
      b  Fig5b_epithelial_nuclear_area             epithelial nuclear area
      c  Fig5c_ROI_celltype (+ Fig5c_celltype_legend, also the key of e)   WT01 ROI 06 cell types
      d  Fig5c_ROI_{proliferation,nfkb,glycolysis,LDHA}                    same ROI, scores
      e  Fig5d_TMA_core_celltype                   TMA A core 2 cell types
      f  Fig5e_immune_composition_hypoxia          immune composition GAMs (hypoxia)
      g  Fig5f_TMA_core_hypoxia                    TMA core hypoxia
      h  Fig5g_TMA_core_highlights                 SecB + immune highlights; each immune
                                                   core sits under its GAM
      i  Fig5h_macrophage_programs_polarization    macrophage program GAMs (styled as a)
      j  Fig5i_macrophage_cell_area                macrophage cell area (as b)

MANUSCRIPT PANEL(S): Fig 5a–j
RUNTIME TIER: moderate
All inputs are SingleR / raw-count outputs (see figures/figure5/).
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from scipy.stats import wilcoxon

import os
R = os.environ.get("OUTPUT_ROOT", "output"); FIG = f"{R}/figures/figure5"; SR = f"{R}/split_recompute"
ROI = f"{FIG}/roi_data"
FA, FK, FL, FP = 7, 6, 7, 10          # axis label, tick, curve label, panel letter (pt)
LW = 0.5
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": LW,
                     "xtick.major.width": LW, "ytick.major.width": LW,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                     "xtick.labelsize": FK, "ytick.labelsize": FK})
PAL = {"SecA epithelium": "#E6A141", "Intermediate epithelium": "#C08E48", "SecB epithelium": "#9A7D55",
       "Ciliated epithelium": "#E07850", "Mesothelial": "#D4A574", "Fibroblast": "#DDD5CA",
       "Smooth muscle": "#D14E6C", "Pericyte": "#B87A7A", "Endothelial": "#7D4E4E",
       "Macrophage": "#8FBC8F", "Mast cell": "#8B9B6B", "Neutrophil": "#6B8E23",
       "Conventional dendritic cell": "#2E8B57", "Plasmacytoid dendritic cell": "#A0A0A0",
       "NK cell": "#56AFC4", "T cell": "#87CEFA", "Plasma cell": "#8A5DAF", "B cell": "#5665B6"}
SHORT = {"SecA epithelium": "SecA epi", "Intermediate epithelium": "Int. epi", "SecB epithelium": "SecB epi",
         "Ciliated epithelium": "Ciliated epi", "Conventional dendritic cell": "cDC",
         "Plasmacytoid dendritic cell": "pDC"}
GRAD = LinearSegmentedColormap.from_list("grad", list(zip([0, .28, .55, .80, 1],
                                         ["#F6EFE5", "#ECDDD0", "#D89E97", "#A03A4A", "#3A111A"])))
GREY, GREY_EDGE = "#E8E8E8", "#C4C4C4"

# ---------------------------------------------------------------- page grid --
W, H = 7.5, 10.0
fig = plt.figure(figsize=(W, H))
COMP, _CUR = {}, [None]          # component name -> its axes (for per-component export)
def comp(name): _CUR[0] = name
def ax_in(x, y, w, h, **kw):
    """axes from inches, (x, y) = top-left corner measured from the page top-left."""
    a = fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)
    COMP.setdefault(_CUR[0], []).append(a)
    return a
def style(a):
    for s in ("top", "right"): a.spines[s].set_visible(False)
    a.tick_params(labelsize=FK, width=LW, length=2.5, pad=1.5)
def fmt_p(p):
    return "p < 0.0001" if p < 1e-4 else (f"p = {p:.4f}" if p < 1e-3 else f"p = {p:.3f}")
def bracket(ax, x1, x2, y, txt, bold=True):
    yr = ax.get_ylim(); h = (yr[1] - yr[0]) * 0.04
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=0.4, clip_on=False)
    ax.text((x1 + x2) / 2, y + h * 1.3, txt, ha="center", va="bottom", fontsize=FK, clip_on=False,
            fontweight="bold" if bold else "normal")

def repel(ys, gap):
    pos = list(ys)
    for _ in range(300):
        moved = False
        for i in range(1, len(pos)):
            if pos[i] - pos[i - 1] < gap:
                m = (pos[i] + pos[i - 1]) / 2; pos[i - 1] = m - gap / 2; pos[i] = m + gap / 2; moved = True
        if not moved: break
    return pos

def gam_panel(ax, curves, xcol, order, colors, xlab, ylab, gap, lw=1.4, stats=None):
    ax.axhline(0, ls=":", color="#999999", lw=0.4, zorder=0)
    ends = []
    for f in order:
        d = curves[curves.feature == f].sort_values(xcol)
        c = colors[f]
        ax.fill_between(d[xcol], d.lower, d.upper, color=c, alpha=0.13, lw=0, zorder=1)
        ax.plot(d[xcol], d.fitted, color=c, lw=lw, zorder=2, solid_capstyle="round")
        ends.append((d.fitted.iloc[-1], f, c))
    ends.sort(); xmax = curves[xcol].max(); xmin = curves[xcol].min()
    for (y, f, c), yp in zip(ends, repel([e[0] for e in ends], gap)):
        xt = xmax + (xmax - xmin) * 0.015
        if stats is None:
            ax.text(xt, yp, f, color=c, fontsize=FL, fontweight="bold", va="center", ha="left", clip_on=False)
        else:   # label + per-tissue statistic underneath
            ax.text(xt, yp, f, color=c, fontsize=FL, fontweight="bold", va="bottom", ha="left", clip_on=False)
            n, pv = stats[f]
            ax.text(xt, yp, f"{n}/8 tissues, {fmt_p(pv)}", fontsize=5.5, va="top", ha="left", clip_on=False,
                    color="black" if pv < 0.05 else "#666666", fontweight="bold" if pv < 0.05 else "normal")
    ax.set_xlim(xmin - (xmax - xmin) * 0.01, xmax + (xmax - xmin) * 0.005)
    ax.set_xlabel(xlab, fontsize=FA, labelpad=2); ax.set_ylabel(ylab, fontsize=FA, labelpad=2); style(ax)

def paired_box(ax_wt, ax_tma, wt, tma, cats, labels, cols, ylab, wt_title, tma_title):
    """b / i : WT paired traces (all categories) + TMA box of first vs last category."""
    x = np.arange(len(cats)) + 1
    for row in wt[cats].values:
        ax_wt.plot(x, row, color="grey", lw=0.6, alpha=0.6, zorder=2)
    for j, c in enumerate(cats):
        ax_wt.scatter(np.full(len(wt), x[j]), wt[c], s=12, color=cols[j], edgecolor="white", lw=0.3, zorder=3)
    p = wilcoxon(wt[cats[0]], wt[cats[-1]]).pvalue
    ax_wt.set_xticks(x); ax_wt.set_xticklabels(labels, fontsize=FK); ax_wt.set_xlim(0.5, len(cats) + 0.5)
    ax_wt.set_ylabel(ylab, fontsize=FA, labelpad=2); ax_wt.set_title(wt_title, fontsize=FA, pad=12)
    style(ax_wt); bracket(ax_wt, 1, len(cats), wt.values.max() * 1.01, fmt_p(p), p < 0.05)
    pair = [cats[0], cats[-1]]; rng = np.random.default_rng(42)
    bp = ax_tma.boxplot([tma[c].values for c in pair], positions=[1, 2], widths=0.55, showfliers=False,
                        patch_artist=True, medianprops=dict(color="black", lw=0.8),
                        whiskerprops=dict(color="black", lw=0.5), capprops=dict(color="black", lw=0.5),
                        boxprops=dict(lw=0.5))
    for patch, c in zip(bp["boxes"], (cols[0], cols[-1])):
        patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor("black")
    for j, c in enumerate(pair):
        v = tma[c].values
        ax_tma.scatter(1 + j + rng.uniform(-0.14, 0.14, len(v)), v, s=1.5, color=(cols[0], cols[-1])[j],
                       alpha=0.3, edgecolors="none", zorder=1)
    p2 = wilcoxon(tma[pair[0]], tma[pair[1]]).pvalue
    ax_tma.set_xticks([1, 2]); ax_tma.set_xticklabels([labels[0], labels[-1]], fontsize=FK)
    ax_tma.set_xlim(0.4, 2.6); ax_tma.set_title(tma_title, fontsize=FA, pad=12); style(ax_tma)
    top = max(c.get_ydata().max() for c in bp["caps"])
    bracket(ax_tma, 1, 2, top * 1.01, fmt_p(p2), p2 < 0.05)
    return p, p2

# ---------------------------------------------------------------- polygons --
def load_poly(stem):
    v = pd.read_csv(f"{ROI}/{stem}_poly.csv.gz"); a = pd.read_csv(f"{ROI}/{stem}_attr.csv")
    polys, ids = [], []
    for (cell, part), g in v.groupby(["cell", "part"], sort=False):
        polys.append(g[["x", "y"]].values); ids.append(cell)
    a = a.set_index("barcode")
    return polys, np.array(ids), a
def draw_cells(ax, polys, faces, edge=GREY_EDGE, lw=0.08, order=None):
    # rasterized: in the SVG/PDF the cell polygons are one 600-dpi image (2026-10-21, author); titles,
    # colour bars, boxes and scale bars stay vector. PNG output is unchanged.
    idx = np.arange(len(polys)) if order is None else order
    ax.add_collection(PolyCollection([polys[i] for i in idx], facecolors=[faces[i] for i in idx],
                                     edgecolors=edge, linewidths=lw, rasterized=True))
def frame(ax, polys, pad=8):
    allv = np.vstack(polys)
    x0, y0 = allv.min(0) - pad; x1, y1 = allv.max(0) + pad
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal"); ax.axis("off")
    return x0, x1, y0, y1
def scalebar(ax, um, x0, x1, y0, y1, label=True):
    L = um; xs = x1 - (x1 - x0) * 0.05 - L; ys = y0 + (y1 - y0) * 0.04
    ax.plot([xs, xs + L], [ys, ys], color="black", lw=1.2, solid_capstyle="butt")
    if label: ax.text(xs + L / 2, ys + (y1 - y0) * 0.025, f"{um} µm", ha="center", va="bottom", fontsize=5.5)

def value_map(ax, polys, ids, a, col, epi_only, title, italic=False, tloc="center"):
    vals = a.loc[ids, col].values.astype(float)
    lab = a.loc[ids, "lab"].values
    show = np.isin(lab, ["SecA epithelium", "Intermediate epithelium", "SecB epithelium",
                         "Ciliated epithelium"]) if epi_only else np.ones(len(ids), bool)
    lo, hi = np.nanpercentile(vals[show], [2, 98])
    norm = Normalize(lo, hi)
    faces = np.array([GRAD(norm(v)) if s and np.isfinite(v) else matplotlib.colors.to_rgba(GREY)
                      for v, s in zip(vals, show)], dtype=object)
    order = np.lexsort((np.where(show, vals, -np.inf),))
    draw_cells(ax, polys, faces, order=order)
    ext = frame(ax, polys)
    ax.set_title(title, fontsize=FA, pad=3 if tloc == "left" else 2, fontweight="bold",
                 style="italic" if italic else "normal", loc=tloc)
    return norm, ext

def cbar(fig_ax, norm, label=None, ticks_fmt="{:.2f}"):
    sm = plt.cm.ScalarMappable(norm=norm, cmap=GRAD)
    cb = fig.colorbar(sm, cax=fig_ax, orientation="horizontal")
    cb.outline.set_linewidth(0.4); cb.ax.tick_params(labelsize=5.5, width=0.4, length=1.5, pad=1)
    cb.set_ticks([norm.vmin, norm.vmax]); cb.set_ticklabels(["low", "high"])
    if label: cb.set_label(label, fontsize=5.5, labelpad=1)

from matplotlib.ticker import MaxNLocator
def cbar_top(fig_ax, norm, dec=2):
    """small horizontal colour bar, numeric ticks labelled ABOVE the bar"""
    sm = plt.cm.ScalarMappable(norm=norm, cmap=GRAD)
    cb = fig.colorbar(sm, cax=fig_ax, orientation="horizontal")
    cb.outline.set_linewidth(0.4)
    t = [v for v in MaxNLocator(nbins=3, steps=[1, 2, 2.5, 5, 10]).tick_values(norm.vmin, norm.vmax)
         if norm.vmin <= v <= norm.vmax]
    cb.set_ticks(t); cb.set_ticklabels([f"{v:.{dec}f}".rstrip("0").rstrip(".") if dec else f"{v:.0f}" for v in t])
    fig_ax.xaxis.set_ticks_position("top"); fig_ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=5, width=0.4, length=1.5, pad=0.8)

# ========================================================================= a
cur = pd.read_csv(f"{FIG}/gam_curves_5a.csv")
ORD_A = ["NF-κB", "Glycolysis", "LDHA", "TNFA", "Oncogenic growth", "DNA damage repair",
         "Epigenetic silencing", "Proliferation"]
COL_A = {"NF-κB": "#8A5DAF", "Glycolysis": "#D14E6C", "LDHA": "#B87A7A", "TNFA": "#C04040",
         "Oncogenic growth": "#5665B6", "DNA damage repair": "#8FBC8F",
         "Epigenetic silencing": "#56AFC4", "Proliferation": "#E6A141"}
ROW1_Y, ROW1_H = 0.22, 1.9
comp("a_epithelial_programs")
ax = ax_in(0.55, ROW1_Y, 3.35, ROW1_H)
gam_panel(ax, cur, "polarization", ORD_A, COL_A,
          "Secretory cell polarization score (SecA → SecB)",
          "Secretory epithelial cell-intrinsic\nprogram score (z)", gap=0.23)
ax.set_xticks(np.arange(-0.5, 0.51, 0.25))

# ========================================================================= b
mb = pd.read_csv(f"{SR}/fig5b_morphometrics_presplit.csv")
EPI = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]
def wideb(src):
    return mb[mb.source == src].pivot_table(index="unit", columns="cell_label_split",
                                           values="nucleus_area").reindex(columns=EPI).dropna()
wtb, tmb = wideb("WT"), wideb("TMA")
comp("b_nuclear_area_whole_tissue"); axw = ax_in(5.55, ROW1_Y + 0.25, 0.95, ROW1_H - 0.25)
comp("b_nuclear_area_TMA");          axt = ax_in(6.72, ROW1_Y + 0.25, 0.72, ROW1_H - 0.25)
pb = paired_box(axw, axt, wtb, tmb, EPI, ["SecA", "Int.", "SecB"], [PAL[c] for c in EPI],
                "Nuclear area (µm²)", f"Whole tissue (n={len(wtb)})", f"TMA (n={len(tmb)})")

# ========================================================================= c / d (WT01 ROI 06)
ROW2_Y, MAP = 2.62, 1.12
polys, ids, a = load_poly("WT01_roi06")
LEG_X = 0.0
# cell-type legend (shared by panels c and e), left column
present = [k for k in PAL if k in set(a.lab) | set(pd.read_csv(f"{ROI}/tma_core2_attr.csv").lab)]
comp("cd_celltype_legend")
lax = ax_in(LEG_X, ROW2_Y + 0.05, 1.0, 3.0); lax.axis("off")
lax.text(0, 1, "Cell type", fontsize=FA, fontweight="bold", va="top", transform=lax.transAxes)
for i, k in enumerate(present):
    yy = 1 - 0.075 - i * 0.052
    lax.add_patch(Rectangle((0, yy - 0.018), 0.09, 0.036, color=PAL[k], transform=lax.transAxes, clip_on=False))
    lax.text(0.13, yy, SHORT.get(k, k), fontsize=FK, va="center", transform=lax.transAxes)
comp("c_ROI_celltype")
cax_ = ax_in(1.05, ROW2_Y + 0.12, MAP, MAP)
faces = np.array([matplotlib.colors.to_rgba(PAL.get(l, GREY)) for l in a.loc[ids, "lab"]], dtype=object)
draw_cells(cax_, polys, faces, edge="#8A8A8A", lw=0.06)
ext = frame(cax_, polys); scalebar(cax_, 100, *ext)
MAPS_C = [("proliferation", "Proliferation", False), ("nfkb", "NF-κB", False),
          ("glycolysis", "Glycolysis", False), ("LDHA", "LDHA", True)]
x0c, stepc = 2.55, 1.24
for j, (col, title, it) in enumerate(MAPS_C):
    comp(f"c_ROI_{col}")
    mx = ax_in(x0c + j * stepc, ROW2_Y + 0.12, MAP, MAP)
    norm, _ = value_map(mx, polys, ids, a, col, True, title, italic=it, tloc="left")
    cbar_top(ax_in(x0c + j * stepc + MAP * 0.56, ROW2_Y + 0.12 - 0.075, MAP * 0.44, 0.045), norm)

# ========================================================================= e (TMA core) / f (immune composition)
ROW3_Y = 4.30
cp, cid, ca = load_poly("tma_core2")
CORE = 1.12
comp("d_core_celltype")
dax = ax_in(1.05, ROW3_Y + 0.2, CORE, CORE)
faces = np.array([matplotlib.colors.to_rgba(PAL.get(l, GREY)) for l in ca.loc[cid, "lab"]], dtype=object)
draw_cells(dax, cp, faces, edge="#8A8A8A", lw=0.06)
cext = frame(dax, cp); scalebar(dax, 100, *cext)
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from config.study_codes import code_of
_code = code_of("TMA A", ca.patient_id.iloc[0])   # CHUM = TMA A; private key, never written out
dax.set_title(f"TMA core ({_code})", fontsize=FA, pad=2, fontweight="bold")   # study code (de-identified)

ps = pd.read_csv(f"{FIG}/fig5c_gam_persample.csv"); ps = ps[ps.cohort == "WT"]
pc = pd.read_csv(f"{SR}/fig5c_tme_percell_presplit.csv"); pc = pc[(pc.cohort == "WT") & np.isfinite(pc.x_hyp)]
TYPES = ["Macrophage", "T cell", "NK cell", "B cell", "Plasma cell"]
GX0, GSTEP, GW = 2.72, 0.955, 0.80                      # shared column grid of panels f and h
E_H = 1.38
stats_e = {}
for j, ct in enumerate(TYPES):
    comp("e_immune_" + ct.replace(" ", "_"))
    ax = ax_in(GX0 + j * GSTEP + 0.08, ROW3_Y + 0.22, GW - 0.08, E_H)
    c = PAL[ct]; sub = ps[ps.cell_type == ct]; grid = np.sort(sub.x.unique()); M = []
    for smp, g in sub.groupby("sample"):
        g = g.sort_values("x")
        ax.plot(g.x, 100 * g.fit, color=c, lw=0.4, ls=(0, (1.6, 1.4)), alpha=0.65, zorder=1)
        M.append(np.interp(grid, g.x, 100 * g.fit))
    M = np.vstack(M); mu = M.mean(0); se = M.std(0, ddof=1) / np.sqrt(M.shape[0])
    ax.fill_between(grid, mu - 1.96 * se, mu + 1.96 * se, color=c, alpha=0.2, lw=0, zorder=2)
    ax.plot(grid, mu, color=c, lw=1.5, zorder=3, solid_capstyle="round")
    lo_, hi_ = [], []
    for s in sorted(ps["sample"].unique()):          # test on ALL 8 tissues (raw bin fractions)
        g = pc[pc["sample"] == s]; q1, q2 = np.percentile(g.x_hyp, [20, 80])
        lo_.append(100 * (g[g.x_hyp <= q1].cell_label == ct).mean())
        hi_.append(100 * (g[g.x_hyp >= q2].cell_label == ct).mean())
    p = wilcoxon(lo_, hi_).pvalue; stats_e[ct] = (mu[0], mu[-1], p)
    ax.set_title(ct, fontsize=FA, color=c, fontweight="bold", pad=2)
    nfit = sub["sample"].nunique()
    if nfit < 8:   # tissues with < 50 cells of this type have no GAM fit
        ax.text(0.96, 0.02, f"curve: {nfit}/8 tissues", transform=ax.transAxes, ha="right",
                va="bottom", fontsize=5, color="#666666")
    ax.text(0.04, 0.98, f"{mu[0]:.1f}→{mu[-1]:.1f}%\n{fmt_p(p)}", transform=ax.transAxes, va="top",
            fontsize=5.5, color="black" if p < 0.05 else "#666666", fontweight="bold" if p < 0.05 else "normal")
    ax.set_xticks([-4, -2, 0, 2]); style(ax)
    if j == 0: ax.set_ylabel("% of local niche\n(mean of tissues)", fontsize=FA, labelpad=2)
    if j == 2: ax.set_xlabel("Local epithelial hypoxia (z)", fontsize=FA, labelpad=2)

# ========================================================================= g (TMA core hypoxia) / h (highlights)
ROW4_Y = 6.22
comp("f_core_hypoxia")
fax = ax_in(0.05, ROW4_Y + 0.2, CORE, CORE)
norm, fext = value_map(fax, cp, cid, ca, "hypoxia", False, "Hypoxia")
cbar(ax_in(0.31, ROW4_Y + 0.2 + CORE + 0.05, CORE - 0.52, 0.06), norm)
_lab = ca.loc[cid, "lab"].values
_cen = np.array([p_.mean(0) for p_ in cp])
_sb = _cen[_lab == "SecB epithelium"]
BOX = (*np.percentile(_sb, 5, axis=0) - 25, *np.percentile(_sb, 95, axis=0) + 25)   # x0, y0, x1, y1
def roi_box(ax):
    ax.add_patch(Rectangle((BOX[0], BOX[1]), BOX[2] - BOX[0], BOX[3] - BOX[1], fill=False,
                           edgecolor="black", lw=0.7, zorder=5))
roi_box(fax)
def highlight(ax, ct, title):
    lab = ca.loc[cid, "lab"].values; hit = lab == ct
    faces = np.array([matplotlib.colors.to_rgba(PAL[ct]) if h else matplotlib.colors.to_rgba(GREY)
                      for h in hit], dtype=object)
    order = np.argsort(hit, kind="stable")
    draw_cells(ax, cp, faces, order=order); frame(ax, cp); roi_box(ax)
    ax.set_title(f"{title} (n={hit.sum()})", fontsize=FK, color=PAL[ct], fontweight="bold", pad=2)
comp("g_core_SecB"); g0 = ax_in(1.25, ROW4_Y + 0.2, CORE, CORE); highlight(g0, "SecB epithelium", "SecB epi")
for j, ct in enumerate(TYPES):
    comp("g_core_" + ct.replace(" ", "_"))
    highlight(ax_in(GX0 + j * GSTEP + 0.08 + (GW - 0.08 - CORE * 0.86) / 2, ROW4_Y + 0.2 + 0.08,
                    CORE * 0.86, CORE * 0.86), ct, ct)

# ========================================================================= i (macrophage programs) / j (cell area)
ROW5_Y = 7.72
mac = pd.read_csv(f"{FIG}/fig5e_mac_gam_curves.csv"); mac = mac[mac.axis == "nb_pol"]
ORD_H = ["Glycolytic", "NF-κB activated", "Chemokine", "Tissue-resident identity"]
COL_H = {"Glycolytic": "#D14E6C", "NF-κB activated": "#8A5DAF", "Chemokine": "#4E9C81",
         "Tissue-resident identity": "#8FBC8F"}
# per-tissue test for panel i (same design as panel f): mean score in the top vs bottom 20% of
# local epithelial polarization within each tissue, paired Wilcoxon across the 8 tissues
mpc = pd.read_csv(f"{SR}/fig5e_macrophage_percell_presplit.csv"); mpc = mpc[np.isfinite(mpc.nb_pol)]
HCOL = {"Glycolytic": "glycolytic", "NF-κB activated": "nfkb_activated", "Chemokine": "chemokine",
        "Tissue-resident identity": "tissue_resident"}
stats_h = {}
for f, col in HCOL.items():
    lo_, hi_ = [], []
    for smp, g in mpc.groupby("sample"):
        q1, q2 = np.percentile(g.nb_pol, [20, 80])
        lo_.append(g[g.nb_pol <= q1][col].mean()); hi_.append(g[g.nb_pol >= q2][col].mean())
    d_ = np.array(hi_) - np.array(lo_); up = np.median(d_) > 0
    stats_h[f] = (int(np.sum(d_ > 0) if up else np.sum(d_ < 0)), wilcoxon(hi_, lo_).pvalue)
comp("h_macrophage_programs")
ax = ax_in(0.55, ROW5_Y + 0.22, 3.35, ROW1_H - 0.15)
gam_panel(ax, mac, "x", ORD_H, COL_H, "Local epithelial polarization (SecA → SecB, z)",
          "Macrophage-intrinsic\nprogram score (z)", gap=0.05)   # descriptive: stats_h -> table only
um = pd.read_csv(f"{SR}/fig5e_mac_morphology_unitmedians_presplit.csv")
def widei(co):
    return um[um.cohort == co].pivot_table(index="unit", columns="tertile", values="cell_area")[["low", "high"]].dropna()
wti, tmi = widei("WT"), widei("TMA")
comp("i_macrophage_area_whole_tissue"); axw = ax_in(5.55, ROW5_Y + 0.47, 0.95, ROW1_H - 0.45)
comp("i_macrophage_area_TMA");          axt = ax_in(6.72, ROW5_Y + 0.47, 0.72, ROW1_H - 0.45)
pi_ = paired_box(axw, axt, wti, tmi, ["low", "high"], ["SecA\nniche", "SecB\nniche"],
                 [PAL["SecA epithelium"], PAL["SecB epithelium"]], "Macrophage cell area (µm²)",
                 f"Whole tissue (n={len(wti)})", f"TMA (n={len(tmi)})")

# ---- one file per subpanel: the page above is not saved. Each subpanel hides every other
#      component and is cropped to what remains, so it keeps its page size and fonts;
#      hidden artists are not written.
import os
CDIR = f"{FIG}/panels"; os.makedirs(CDIR, exist_ok=True)
all_axes = [a_ for v in COMP.values() for a_ in v]
from matplotlib.transforms import Bbox
fig.canvas.draw(); rend = fig.canvas.get_renderer()
def comp_bbox(axs, pad=0.02):
    """matplotlib's own tight box for the visible component, extended to cover any
    axis label / title / text that is wider than its axes (e.g. the shared x-label of panel f)."""
    tb = fig.get_tightbbox(rend)                                   # inches, visible artists only
    ext = []
    for a_ in axs:
        for t in (a_.xaxis.label, a_.yaxis.label, a_.title, *a_.texts):
            if t.get_visible() and t.get_text():
                ext.append(t.get_window_extent(rend).transformed(fig.dpi_scale_trans.inverted()))
    bb = Bbox.union([tb] + ext)
    return Bbox.from_extents(bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad)
PANELS = {   # subpanel file -> the components drawn for it (file letters -> panels: see the module docstring)
    "Fig5a_epithelial_programs":            ["a_epithelial_programs"],
    "Fig5b_epithelial_nuclear_area":        ["b_nuclear_area_whole_tissue", "b_nuclear_area_TMA"],
    "Fig5c_celltype_legend":                ["cd_celltype_legend"],
    "Fig5c_ROI_celltype":                   ["c_ROI_celltype"],
    "Fig5c_ROI_proliferation":              ["c_ROI_proliferation"],
    "Fig5c_ROI_nfkb":                       ["c_ROI_nfkb"],
    "Fig5c_ROI_glycolysis":                 ["c_ROI_glycolysis"],
    "Fig5c_ROI_LDHA":                       ["c_ROI_LDHA"],
    "Fig5d_TMA_core_celltype":              ["d_core_celltype"],
    "Fig5e_immune_composition_hypoxia":     [k for k in COMP if k.startswith("e_immune_")],
    "Fig5f_TMA_core_hypoxia":               ["f_core_hypoxia"],
    "Fig5g_TMA_core_highlights":            [k for k in COMP if k.startswith("g_core_")],
    "Fig5h_macrophage_programs_polarization": ["h_macrophage_programs"],
    "Fig5i_macrophage_cell_area":           ["i_macrophage_area_whole_tissue", "i_macrophage_area_TMA"],
}
assert sorted(sum(PANELS.values(), [])) == sorted(COMP), "every component must belong to one subpanel"
for name, comps in PANELS.items():
    axs = [a_ for c_ in comps for a_ in COMP[c_]]
    for a_ in all_axes: a_.set_visible(a_ in axs)
    fig.canvas.draw(); bb = comp_bbox(axs)
    for ext in (("png",) if name.startswith("Fig5c_ROI_") else ("svg", "png")):   # ROI maps: PNG, as Fig 4
        fig.savefig(f"{CDIR}/{name}.{ext}", bbox_inches=bb,
                    dpi=450 if ext == "png" else 600, facecolor="white")   # svg: dpi of the rasterized cells
    print(f"  {name:42s} {bb.width:.2f} x {bb.height:.2f} in")
for a_ in all_axes: a_.set_visible(True)
print(f"subpanels: {len(PANELS)} -> {CDIR}")
print("b  WT p=%.4g TMA p=%.4g (n=%d/%d)" % (pb[0], pb[1], len(wtb), len(tmb)))
print("i  WT p=%.4g TMA p=%.4g (n=%d/%d)" % (pi_[0], pi_[1], len(wti), len(tmi)))
for k, v in stats_e.items(): print(f"e  {k:12s} {v[0]:.1f} -> {v[1]:.1f}%  p={v[2]:.4f}")
for k, v in stats_h.items(): print(f"h  {k:24s} {v[0]}/8  p={v[1]:.4f}")
pd.DataFrame([dict(panel="e", feature=k, pct_low=v[0], pct_high=v[1], p_paired_wilcoxon=v[2], n_tissues=8)
              for k, v in stats_e.items()] +
             [dict(panel="h", feature=k, n_concordant=v[0], p_paired_wilcoxon=v[1], n_tissues=8)
              for k, v in stats_h.items()]).to_csv(f"{FIG}/Figure5_composite_stats_eh.csv", index=False)
print("saved", f"{FIG}/Figure5_composite_stats_eh.csv")
