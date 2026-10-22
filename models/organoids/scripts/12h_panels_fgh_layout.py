#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202603_organoids/scripts/12h_panels_fgh_layout.py, the script that ran. Kept:
# fig3F_marker_vs_function (Fig 3f) and fig3GH_shared (Fig 3g, h). Removed: the alternate single-panel G and H
# layouts. The untrimmed script, which the header below describes, is in the lab project folder.

# 12h — Letter-page LAYOUT revision of the SecB marker-vs-function panels F/G/H.
# Preserves 12g_split_panels_fgh.py (originals) — this is a new script.
# Data / scoring / points / thresholds are IDENTICAL to 12g/12c
# (12_scores_all_sources.csv; Function = z(HIF) - z(Proliferation)). Layout only.
#
# Changes requested by the Figure-3 assembly:
#   F  fig3F_marker_vs_function  — REMOVE atlas reference points (plot + legend);
#        instead SHADE the top-right "real SecB niche" rectangle (from the niche
#        crosshair to the top-right corner) in SecB brown (#5A3E2B, alpha 0.18)
#        and label it "SecB niche". Keep all in-vitro systems + their legend.
#   G/H individual  fig3G_HIF_induction / fig3H_proliferation — unchanged content,
#        redrawn at true footprint (no tight-bbox crop) so 8 pt renders true.
#   GH combined  fig3GH_shared — NEW: G (HIF) on top, H (proliferation) on bottom,
#        two vertically-stacked axes sharing the systems x-axis (labels drawn once
#        at the very bottom). Each keeps its own y-axis + atlas-SecB reference line.
#        Sized so combined G+H ~ height of F (~3.0 in tall, ~3.6 in wide).
#
# All panels: fonts 8 pt text / 6 pt ticks; drawn at true footprint (no tight
# crop) so the SVG placed 1:1 in Affinity renders 8 pt as true 8 pt.
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgba

# ---- data (identical to 12c/12g) ----
SRC = "/Volumes/CookLab/Sarah/202603_organoids/output/12_GSE284507_spheroid_organoid"
OUT = "/Users/sarahnersesian/202609_spatial_rerun/figures/figure3"
os.makedirs(OUT, exist_ok=True)
r = pd.read_csv(f"{SRC}/12_scores_all_sources.csv")
r["hgsc"] = r["hgsc"].astype(str).str.upper() == "TRUE"
def z(x): return (x - x.mean()) / (x.std() + 1e-9)
r["Function"] = z(r["Hypoxia_HIF"]) - z(r["Proliferation"])   # same definition as 12c
# Quiescence = -Proliferation: plot the negated proliferation so that, like HIF,
# UP = more SecB-like (more quiescent) and the atlas-SecB line sits near the TOP.
# No scoring changes — the plotted value is simply multiplied by -1.
r["Quiescence"] = -r["Proliferation"]

# ---- style: 8 pt text / 6 pt tick+category labels ----
FA, FK = 8, 6
plt.rcParams.update({"font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FA, "axes.labelsize": FA, "axes.titlesize": 0,
    "xtick.labelsize": FK, "ytick.labelsize": FK, "legend.fontsize": FA,
    "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 450, "figure.dpi": 150})

GOLD, AMBER, BROWN = "#E6A141", "#C08E48", "#6B4D2E"   # atlas Epithelial / (unused) / SecB brown
# Atlas-lineage palette recolour (system -> atlas level-1 colour):
ASC   = "#5665B6"   # Our PDO ascites -> atlas B cell (dark blue/indigo)
PRI   = "#87CEFA"   # Our PDO primary -> atlas T/NK cell (light blue)
HG    = "#8A5DAF"   # DepMap HGSC     -> atlas Plasma cell (purple)
DLP_P = "#D4A574"   # Cold co-cx physoxia/normoxia -> atlas Mesothelial (beige)
DLP_H = "#D14E6C"   # Cold co-cx hypoxia          -> atlas Smooth muscle (pink/red)
G_ORG = "#8FBC8F"   # GSE284507 organoid -> atlas Macrophage (light green)
G_SPH = "#2E8B57"   # GSE284507 spheroid -> atlas DC (darker green)
NICHE_FILL = "#5A3E2B"   # SecB-brown shading for the "real SecB niche" region
def _origin_colors(samp):   # per-point colour for Our-PDO by ascites/primary origin
    return [ASC if str(s).startswith("OCAD") else PRI for s in samp]
def amed(sub, col): return r.loc[(r.grp == "atlas") & (r["sub"] == sub), col].median()


# ============================================================ Panel F: scatter
figF = plt.figure(figsize=(4.4, 3.0), constrained_layout=True)   # widened; legend sits to the right
axA = figF.add_subplot(111)
# in-vitro systems only (atlas reference points removed per request)
# EVERY system is a CIRCLE at the SAME uniform size (MS). Fill is slightly
# transparent (alpha in the RGBA facecolour) so dense overlaps stay readable,
# while the thin black edge is drawn at full opacity so points stay crisp.
MS = 18          # uniform marker size for ALL systems
MA = 0.70        # marker fill alpha
def _fc(col):    # translucent facecolour, single colour or list of colours
    return [to_rgba(c, MA) for c in col] if isinstance(col, list) else to_rgba(col, MA)
mo = r.grp == "our"
oc = _origin_colors(r.loc[mo, "sample"])
axA.scatter(r.SecB_marker[mo], r.Function[mo], s=MS, facecolor=_fc(oc), edgecolor="black", lw=0.4, marker="o", zorder=4)
mh = r.hgsc
axA.scatter(r.SecB_marker[mh], r.Function[mh], s=MS, facecolor=_fc(HG), edgecolor="black", lw=0.4, marker="o", zorder=4)
for cc, col in [("physoxia", DLP_P), ("hypoxia", DLP_H)]:
    m = (r.grp == "delapuente") & (r["sub"] == cc)
    axA.scatter(r.SecB_marker[m], r.Function[m], s=MS, facecolor=_fc(col), edgecolor="black", lw=0.4, marker="o", zorder=6)
mgo = (r.grp == "gse284507") & (r["sub"] == "organoid")
mgs = (r.grp == "gse284507") & (r["sub"] == "spheroid")
axA.scatter(r.SecB_marker[mgo], r.Function[mgo], s=MS, facecolor=_fc(G_ORG), edgecolor="black", lw=0.4, marker="o", zorder=5)
axA.scatter(r.SecB_marker[mgs], r.Function[mgs], s=MS, facecolor=_fc(G_SPH), edgecolor="black", lw=0.4, marker="o", zorder=5)

# niche crosshair = atlas SecB median marker (x) & function (y)
sx, sy = amed("SecB epithelium", "SecB_marker"), amed("SecB epithelium", "Function")
# freeze limits from the in-vitro points; the atlas-SecB marker threshold (sx) sits
# at the far-right edge of the in-vitro range, so extend the x-axis a little past it
# to give the top-right niche block visible width (that it stays nearly empty is the
# biological message: in-vitro systems rarely reach the real SecB niche).
xlim, ylim = axA.get_xlim(), axA.get_ylim()
xmax = max(xlim[1], sx + 0.07)
xlim = (xlim[0], xmax)
# darker, near-solid SecB-brown niche block so the empty niche stands out and the
# distance of the in-vitro points from it reads at a glance
axA.add_patch(Rectangle((sx, sy), xmax - sx, ylim[1] - sy,
                        facecolor=NICHE_FILL, alpha=0.50, edgecolor=NICHE_FILL,
                        lw=0.6, zorder=0))
axA.axhline(sy, color=BROWN, ls=(0, (4, 2)), lw=0.7, alpha=.8, zorder=1)
axA.axvline(sx, color=BROWN, ls=(0, (4, 2)), lw=0.7, alpha=.8, zorder=1)
axA.set_xlim(xlim); axA.set_ylim(ylim)
# label sits just LEFT of the (narrow) dark box, in SecB brown on the white margin,
# so the full "SecB niche" reads (white-on-white clipped when placed over the box)
axA.text(sx - 0.006, ylim[1] - 0.06, "SecB niche", fontsize=FK, color=NICHE_FILL,
         va="top", ha="right", fontweight="bold")
axA.set_xlabel("SecB MARKER score  (appearance)")
axA.set_ylabel("SecB FUNCTION  (HIF-induction − proliferation)")
for s in ("top", "right"): axA.spines[s].set_visible(False)
for s in ("left", "bottom"): axA.spines[s].set_linewidth(0.4)
axA.tick_params(width=0.4, length=2)
def _lh(fc, lab):   # circle legend handle in the new atlas colours
    return Line2D([0], [0], marker='o', color='none', markerfacecolor=fc,
                  markeredgecolor='k', markeredgewidth=0.3, markersize=4, ls='none', label=lab)
# Legend pulled OUT of the plotting area, to the RIGHT of the axes as a single
# vertical column (no dot overlap); in assembly it sits to the right of this panel.
axA.legend(handles=[
    _lh(ASC, 'PDO ascites'), _lh(PRI, 'PDO primary'),
    _lh(HG, 'PDO DepMap'), _lh(DLP_P, 'PDO physoxia'),
    _lh(DLP_H, 'PDO hypoxia'), _lh(G_ORG, 'Organoid'),
    _lh(G_SPH, 'Spheroid')],
    loc="center left", bbox_to_anchor=(1.01, 0.5), ncol=1, frameon=False,
    fontsize=FA, handletextpad=0.25, labelspacing=0.4, borderpad=0.1)
figF.savefig(f"{OUT}/fig3F_marker_vs_function.svg")
figF.savefig(f"{OUT}/fig3F_marker_vs_function.png", dpi=450)
plt.close(figF)

# ============================================ Panels G & H: function-arm strips
# Atlas SecA/SecB dot clouds REMOVED — only the atlas-SecB dashed reference line
# (drawn from amed(...) via `ref`) is kept, so the focus stays on model systems.
groups = [("Our PDO", r.grp == "our", PRI),
          ("DepMap HGSC", r.hgsc, HG),
          ("Cold physoxia", (r.grp == "delapuente") & (r["sub"] == "physoxia"), DLP_P),
          ("Cold hypoxia", (r.grp == "delapuente") & (r["sub"] == "hypoxia"), DLP_H),
          ("GSE284507 org", (r.grp == "gse284507") & (r["sub"] == "organoid"), G_ORG),
          ("GSE284507 sph", (r.grp == "gse284507") & (r["sub"] == "spheroid"), G_SPH)]
# compact HORIZONTAL x labels (two-line where needed) so they fit without angling
DISPLAY = {"Our PDO": "Our PDO", "DepMap HGSC": "DepMap",
           "Cold physoxia": "Cold\nphysoxia", "Cold hypoxia": "Cold\nhypoxia",
           "GSE284507 org": "GSE org", "GSE284507 sph": "GSE sph"}


def strip(ax, col, ylab, ref, refc, show_xticklabels=True):
    x = np.arange(len(groups)) + 1
    for i, (lab, m, c) in enumerate(groups):
        v = r.loc[m, col].values
        if len(v) >= 20:
            vp = ax.violinplot([v], positions=[x[i]], widths=0.8, showextrema=False)
            for b in vp["bodies"]:
                b.set_facecolor(c); b.set_edgecolor("none"); b.set_alpha(0.4)
        if len(v) >= 5:
            ax.boxplot([v], positions=[x[i]], widths=0.16, showfliers=False, patch_artist=True,
                       boxprops=dict(facecolor="white", edgecolor="black", lw=0.4),
                       medianprops=dict(color="black", lw=0.8),
                       whiskerprops=dict(color="black", lw=0.35), capprops=dict(color="black", lw=0.35))
        jit = (np.random.default_rng(1).uniform(-0.13, 0.13, len(v)) if len(v) > 1 else np.array([0.0]))
        # Our PDO points coloured by ascites(dark blue)/primary(light blue) origin;
        # every other system uses its single atlas-lineage colour. Circle glyphs,
        # UNIFORM size across systems, translucent fill + full-opacity black edge.
        pt_c = _origin_colors(r.loc[m, "sample"].values) if lab == "Our PDO" else c
        ax.scatter(np.full(len(v), x[i]) + jit, v, s=14, facecolor=_fc(pt_c),
                   edgecolor="black", lw=0.25, marker="o", zorder=5)
    ax.axhline(ref, color=refc, ls=(0, (4, 2)), lw=0.7, alpha=.8)
    ax.set_xticks(x)
    if show_xticklabels:
        ax.set_xticklabels([DISPLAY[g] for g, _, _ in groups], rotation=0,
                           ha="center", fontsize=FK)
    else:
        ax.set_xticklabels([])
    ax.set_ylabel(ylab)
    ax.set_xlim(0.5, len(groups) + 1.5)   # right margin holds the atlas-SecB label
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_linewidth(0.4)
    ax.tick_params(width=0.4, length=2)


# ---- Combined G + H with shared x-axis (G top, H bottom) ----
figGH, (axg, axh) = plt.subplots(2, 1, sharex=True, figsize=(3.6, 3.0),
                                 constrained_layout=True)
strip(axg, "Hypoxia_HIF", "HIF induction", amed("SecB epithelium", "Hypoxia_HIF"),
      BROWN, show_xticklabels=False)
axg.text(len(groups) + 0.15, amed("SecB epithelium", "Hypoxia_HIF"), "atlas SecB",
         fontsize=FK, color=BROWN, va="bottom", ha="left")
strip(axh, "Quiescence", "Quiescence", amed("SecB epithelium", "Quiescence"),
      BROWN, show_xticklabels=True)
axh.text(len(groups) + 0.15, amed("SecB epithelium", "Quiescence"), "atlas SecB",
         fontsize=FK, color=BROWN, va="bottom", ha="left")
axh.set_xticklabels([DISPLAY[g] for g, _, _ in groups], rotation=0, ha="center", fontsize=FK)
figGH.savefig(f"{OUT}/fig3GH_shared.svg")
figGH.savefig(f"{OUT}/fig3GH_shared.png", dpi=450)
plt.close(figGH)

print(f"niche crosshair: SecB_marker={sx:.3f}  Function={sy:.3f}")
