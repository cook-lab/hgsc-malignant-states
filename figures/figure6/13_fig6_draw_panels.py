#!/usr/bin/env python3
"""
Figure 6 — every panel drawn on one letter page grid (7.5 x 10 in), house style of 11_fig5_compose_panels.py.

PURPOSE
    Draws panels a–i. Each panel file is drawn alone (only that panel's objects) on the same page grid and cropped
    to its region, so every panel keeps its page size; the full page is drawn once to record the regions and is
    not saved. Panel keys inside draw() are earlier letters, mapped to the published panels by NAMES below:
      a  Fig6a_visium_section_SecA_SecB           Visium HD WT01 whole section (SecA-high / SecB-high bins)
      b  Fig6b_visium_ROI_SecA_SecB               ROI: SecA-high / SecB-high bins with territory contours
      c  Fig6c_visium_ROI_CA9                     ROI: CA9 (3 x 3-bin mean), territories outlined
      m  Fig6d_visium_ROI_macrophage_niche        ROI: SecB-niche vs non-SecB macrophage bins       (panel d)
      d  Fig6e_visium_niche_macrophage_volcano    Visium niche vs distal macrophage volcano         (panel e)
      e  Fig6f_atlas_hypoxic_macrophage_dotplot   atlas macrophage dot plot + Xenium bridge row     (panel f)
      f  Fig6g_hypoxic_macrophage_SecB_neighbours SecB neighbours, whole tissues + TMA              (panel g)
      g  Fig6h_crossplatform_matrix               cross-platform matrix (Xenium / Visium WT01 / atlas) (panel h)
      h  Fig6i_ROI_macrophage_VEGFA_CTSL_CD44_ICAM1  Xenium ROI, cell type + macrophage genes        (panel i, top row)
      k  Fig6j_ROI_SecB_CD55_IL32_MMP7_MCT4          same ROI, epithelial genes                     (panel i, bottom row)
    Every tile of rows h and k is also written alone (Fig6i1-i5, Fig6j1-j4).

INPUTS
    <fig6_singler>/: visium_n2_bins_WT01.csv.gz (01, 01b), visiumWT01_niche_vs_distal_volcano.csv (02),
    visium_candidates_2sections.csv (03), reciprocal_summary.csv (05), fig6_colocal_units.csv (06),
    xenium_pct_summary.csv (07), atlas_pb/ (08), atlas_hypoxic_vs_other_paired.csv (09), atlas_epi_pct.csv (10),
    fig6_concordance_table.csv (11), atlas_subtype_dotmatrix.csv (08b), fig6_rois.json (ROI picker);
    the Fig 2g DESeq2 table; WT02 cell polygons and counts through fig6_roi_picker/roi_server.py.

OUTPUTS
    <FINAL>/<name>.{pdf,png,svg} for each panel and tile (PNG of rows h, k transparent);
    <fig6_singler>/fig6h_matrix_standardized.csv (the plotted values of panel h).

Usage: python 13_fig6_draw_panels.py ["ROI name"]   (default ROI = "example 2 (700 um)" from fig6_rois.json)
MANUSCRIPT PANEL(S): Fig 6a–i. RUNTIME TIER: moderate.
"""
import sys, os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize, TwoSlopeNorm, to_rgba
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox
from scipy.stats import wilcoxon
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig6_roi_picker"))
import roi_server as RS

R = os.environ.get("OUTPUT_ROOT", "output"); O = f"{R}/fig6_singler"
FINAL = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "final_panels/figure6")
PAN = FINAL; SHOW_LETTERS = True; TILE_ONLY = None; TILEREG = {}; SUPP = False
LETMAP = {"a": "a", "b": "b", "m": "d", "c": "c", "d": "e", "e": "f", "f": "g", "g": "h", "h": "i", "i": "x", "j": "x", "k": "j"}
ROI = {r["name"]: r for r in json.load(open(f"{O}/fig6_rois.json"))}[sys.argv[1] if len(sys.argv) > 1 else "example 2 (700 um)"]
FA, FK, FL, FP, LW = 7, 6, 7, 10, 0.5
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "svg.fonttype": "none",
                     "pdf.fonttype": 42, "axes.linewidth": LW, "xtick.major.width": LW, "ytick.major.width": LW,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.labelsize": FK, "ytick.labelsize": FK})
PAL, GRAD, DIV, GREY, GREY_EDGE = RS.PAL, RS.GRAD, RS.DIV, RS.GREY, RS.GREY_EDGE
HYP, RES, NICHE, SECB, SECA = "#A03A4A", "#8FBC8F", "#2F6B45", PAL["SecB epithelium"], PAL["SecA epithelium"]
W, H = 7.5, 10.0
fig = None
def ax_in(x, y, w, h, **kw): return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)
def letter(x, y, s):
    if SHOW_LETTERS: fig.text(x / W, 1 - y / H, LETMAP.get(s, s), fontsize=FP, fontweight="bold", va="top", ha="left")
def style(a):
    for s in ("top", "right"): a.spines[s].set_visible(False)
    a.tick_params(labelsize=FK, width=LW, length=2.5, pad=1.5)
def fmt_p(p): return "p < 0.0001" if p < 1e-4 else (f"p = {p:.4f}" if p < 1e-3 else f"p = {p:.3f}")
def bracket(ax, x1, x2, y, txt, bold=True):
    yr = ax.get_ylim(); h = (yr[1] - yr[0]) * 0.04
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=0.4, clip_on=False)
    ax.text((x1 + x2) / 2, y + h * 1.3, txt, ha="center", va="bottom", fontsize=5.5, clip_on=False, fontweight="bold" if bold else "normal")
def cbar(cax, norm, cmap=GRAD, lo="low", hi="high"):
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.4); cb.ax.tick_params(labelsize=5, width=0.4, length=1.5, pad=1); cb.set_ticks([norm.vmin, norm.vmax]); cb.set_ticklabels([lo, hi])
from matplotlib.ticker import MaxNLocator
def cbar_num(fig_ax, norm, cmap=None, top=True, dec=1):
    """Fig 5 style small horizontal colour bar with numeric ticks (above if top else below)."""
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap or GRAD), cax=fig_ax, orientation="horizontal"); cb.outline.set_linewidth(0.4)
    t = [v for v in MaxNLocator(nbins=3, steps=[1, 2, 2.5, 5, 10]).tick_values(norm.vmin, norm.vmax) if norm.vmin - 1e-9 <= v <= norm.vmax + 1e-9]
    cb.set_ticks(t); cb.set_ticklabels([f"{v:.{dec}f}".rstrip("0").rstrip(".") if dec else f"{v:.0f}" for v in t])
    if top: fig_ax.xaxis.set_ticks_position("top"); fig_ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=5, width=0.4, length=1.5, pad=0.8); return cb
def stars(p): return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
def legend_sq(x, y, items):
    for i, (t, c) in enumerate(items):
        xx = x + (i % 2) * 1.0; yy = y + (i // 2) * 0.13
        fig.patches.append(Rectangle((xx / W, 1 - (yy + 0.035) / H), 0.07 / W, 0.07 / H, color=c, transform=fig.transFigure))
        fig.text((xx + 0.1) / W, 1 - yy / H, t, fontsize=5.5, va="center")
REG = {}


XR = pd.read_csv(f"{O}/reciprocal_summary.csv")

def draw(only=None):
    global REG, XR
    REG = {}
    want = lambda k: only is None or only == k

    Y1 = 0.12; IM = 1.62; WB, HB = 1.62, 1.62          # row 1: a whole section + three square ROI tiles
    B = pd.read_csv(f"{O}/visium_n2_bins_WT01.csv.gz")
    nm = B[B.cls == "mac (SecB niche)"]; Wx, Wy = 150.0, 90.0   # 16 um bins -> 2.4 x 1.44 mm window, same aspect as the tiles
    P_ = nm[["col", "row"]].values
    cnt = [int(((np.abs(P_[:, 0] - q[0]) <= Wx / 2) & (np.abs(P_[:, 1] - q[1]) <= Wy / 2)).sum()) for q in P_]
    zc, zr = P_[int(np.argmax(cnt))]; ZX = (zc - Wx / 2, zc + Wx / 2); ZY = (zr - Wy / 2, zr + Wy / 2)
    VR = [r for r in json.load(open(f"{O}/fig6_rois.json")) if r.get("kind") == "visium"]   # Visium window saved in the ROI picker (the last entry is used)
    if VR:
        v_ = VR[-1]; Wx = Wy = float(v_["size"]) / 16.0; zc, zr = float(v_["x"]), float(v_["y"])   # square, side = saved window width
        ZX = (zc - Wx / 2, zc + Wx / 2); ZY = (zr - Wy / 2, zr + Wy / 2); print("Visium window:", v_["name"], v_["size"], "um")
    if want('a'):
        letter(0.0, Y1 - 0.07, "a")
        ax = ax_in(0.08, Y1 + 0.18, 1.75, IM)
        CA = {"other": GREY, "mac (other)": GREY, "mac (distal)": GREY, "mac (SecB niche)": GREY, "epithelium": "#EFE6D8", "SecA-hi": SECA, "SecB-hi": SECB}
        for c in ["other", "mac (other)", "mac (distal)", "mac (SecB niche)", "epithelium", "SecA-hi", "SecB-hi"]:
            q = B[B.cls == c]; ax.scatter(q.col, q.row, s=0.12, c=CA[c], marker="s", lw=0, rasterized=True)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
        ax.add_patch(Rectangle((ZX[0], ZY[0]), Wx, Wy, fill=False, ec="black", lw=0.7))
        xl, yl = ax.get_xlim(), ax.get_ylim(); Lb = 1000 / 16
        ax.plot([xl[1] - Lb - 5, xl[1] - 5], [yl[0] - 3, yl[0] - 3], color="k", lw=1.2, solid_capstyle="butt"); ax.text(xl[1] - Lb / 2 - 5, yl[0] - 6, "1 mm", ha="center", va="bottom", fontsize=5.5)
        legend_sq(1.2, Y1 + 0.02, [("SecA-high", SECA)]); legend_sq(1.2, Y1 + 0.15, [("SecB-high", SECB)])   # top right

    REG["a"] = (0.0, 0.0, 1.9, 2.1)
    Bz = B[(B.col.between(*ZX)) & (B.row.between(*ZY))].copy()
    # smoothed SecB (solid brown) / SecA (dashed orange) territory contours on panels b-d
    from scipy.ndimage import gaussian_filter
    nxT, nyT = int(np.ceil(ZX[1] - ZX[0])) + 1, int(np.ceil(ZY[1] - ZY[0])) + 1
    giT = ((Bz.row - ZY[0]).astype(int).clip(0, nyT - 1), (Bz.col - ZX[0]).astype(int).clip(0, nxT - 1))
    def terr(mask):
        F_ = np.zeros((nyT, nxT)); E_ = np.zeros((nyT, nxT)); np.add.at(F_, giT, mask)
        np.add.at(E_, giT, Bz.cls.isin(["SecA-hi", "SecB-hi", "epithelium"]).values.astype(float))
        return gaussian_filter(F_, 2.5) / (gaussian_filter(E_, 2.5) + 0.05)
    fB_, fA_ = terr((Bz.cls == "SecB-hi").values.astype(float)), terr((Bz.cls == "SecA-hi").values.astype(float))
    def territory(a_):
        gx_, gy_ = np.arange(nxT) + ZX[0], np.arange(nyT) + ZY[0]
        a_.contour(gx_, gy_, fB_, levels=[0.45], colors=SECB, linewidths=0.9, zorder=5)
        a_.contour(gx_, gy_, fA_, levels=[0.45], colors=SECA, linewidths=0.9, linestyles="--", zorder=5)
    def roi_ax(x):
        a_ = ax_in(x, Y1 + 0.18, WB, HB); a_.set_aspect("equal"); a_.set_xlim(*ZX); a_.set_ylim(ZY[1], ZY[0]); a_.axis("off"); return a_
    S_BIN = (WB * 72 / Wx) ** 2 * 0.95
    if want('b'):
        letter(1.9, Y1 - 0.07, "b")
        bx = roi_ax(2.02)
        CB = {"other": GREY, "mac (other)": GREY, "epithelium": "#EFE6D8", "SecA-hi": SECA, "SecB-hi": SECB, "mac (distal)": GREY, "mac (SecB niche)": GREY}
        for c in ["other", "epithelium", "mac (other)", "SecA-hi", "SecB-hi", "mac (distal)", "mac (SecB niche)"]:
            q = Bz[Bz.cls == c]; bx.scatter(q.col, q.row, s=S_BIN, c=CB[c], marker="s", lw=0)
        territory(bx)
        Lb = 250 / 16; bx.plot([ZX[1] - Lb - 3, ZX[1] - 3], [ZY[1] - 3, ZY[1] - 3], color="k", lw=1.2); bx.text(ZX[1] - Lb / 2 - 3, ZY[1] - 5, "250 µm", fontsize=5, ha="center")
    REG["b"] = (1.9, 0.0, 3.75, 2.1)
    if want('m'):   # macrophage bins only: SecB-niche vs non-SecB (distal + other) macrophage bins
        letter(5.6, Y1 - 0.07, "m")
        mx_ = roi_ax(5.72)
        CM = {c_: GREY for c_ in ["other", "epithelium", "SecA-hi", "SecB-hi"]}; CM.update({"mac (other)": RES, "mac (distal)": RES, "mac (SecB niche)": NICHE})
        for c in ["other", "epithelium", "SecA-hi", "SecB-hi", "mac (other)", "mac (distal)", "mac (SecB niche)"]:
            q = Bz[Bz.cls == c]; mx_.scatter(q.col, q.row, s=S_BIN, c=CM[c], marker="s", lw=0)
        territory(mx_)
        legend_sq(5.75, Y1 + 1.86, [("Mac, SecB niche", NICHE), ("Mac, non-SecB", RES)])
    REG["m"] = (5.6, 0.0, 7.5, 2.1)
    if want('c'):
        letter(3.75, Y1 - 0.07, "c")
        cx_ = roi_ax(3.87)
        nb = cKDTree(Bz[["col", "row"]].values).query_ball_point(Bz[["col", "row"]].values, r=1.0, p=np.inf)
        v = np.array([Bz["CA9"].values[i].mean() for i in nb]); nn = Normalize(0, np.quantile(v, 0.98)); o_ = np.argsort(v)
        cx_.scatter(Bz.col.values[o_], Bz.row.values[o_], s=S_BIN, marker="s", lw=0, c=[GRAD(nn(x)) for x in v[o_]])
        territory(cx_)
        fig.text(3.87 / W, 1 - (Y1 + 0.18 + HB + 0.06) / H, "CA9", fontsize=FA, fontweight="bold", fontstyle="italic", ha="left", va="top")   # bottom left
        cbar_num(ax_in(3.87 + WB * 0.5, Y1 + 0.18 + HB + 0.07, WB * 0.44, 0.045), nn, top=False)                                          # bottom right
    REG["c"] = (3.75, 0.0, 5.6, 2.25)

    # ======================================================================= row 2: d, e, f
    Y2 = 2.12
    if want('d'):
        letter(0.0, Y2, "d")
        V = pd.read_csv(f"{O}/visiumWT01_niche_vs_distal_volcano.csv")
        PANMAC = ["CD68", "CSF1R", "AIF1", "CD14", "FCGR3A", "TYROBP", "LAPTM5"]
        V["x"] = V.log2FC - V[V.gene.isin(PANMAC)].log2FC.median(); V["y"] = -np.log10(V.fdr_bin.clip(lower=1e-45)); sig = V.fdr_bin < 0.05
        ax = ax_in(0.42, Y2 + 0.28, 1.3, 1.2)
        STROMAL = V.gene.str.startswith("COL") | V.gene.isin(["SPARC", "FN1", "DCN", "LUM", "POSTN", "TIMP3", "BGN", "MGP", "IGFBP7", "AEBP1", "CTHRC1", "VCAN", "TMSB4X", "VIM"])
        STR = "#8C7B6B"; ms, co = sig & V.mac_sourced & ~STROMAL, sig & ~V.mac_sourced; st = sig & STROMAL & (V.x < 0)
        ax.scatter(V.x[~sig], V.y[~sig], s=1, c="#D9D9D9", lw=0, rasterized=True)
        ax.scatter(V.x[co], V.y[co], s=2, facecolors="none", edgecolors=SECB, lw=0.3, alpha=0.5, rasterized=True)
        ax.scatter(V.x[ms & (V.x > 0)], V.y[ms & (V.x > 0)], s=4, c=NICHE, lw=0, zorder=3); ax.scatter(V.x[ms & (V.x < 0)], V.y[ms & (V.x < 0)], s=4, c=RES, lw=0, zorder=3)
        ax.scatter(V.x[st], V.y[st], s=6, marker="^", c=STR, lw=0, zorder=4)
        mt = V.gene.str.startswith("MT-")
        up = V[ms & (V.x > 0.25) & ~mt].sort_values("fdr_bin").head(8); dn = V[ms & (V.x < -0.25) & ~mt].sort_values("fdr_bin").head(3)
        stl = V[st].sort_values("fdr_bin").head(3)
        for dd, col, side in ((up, NICHE, 1), (dn, RES, -1), (stl, STR, -2)):
            for k, (_, r) in enumerate(dd.sort_values("y", ascending=False).iterrows()):
                tx = V.x.max() * 1.03 if side > 0 else V.x.min() * 0.97; ty = V.y.max() * ((0.97 - 0.1 * k) if side > 0 else ((0.5 - 0.09 * k) if side == -1 else (0.97 - 0.09 * k)))
                ax.annotate(r.gene, (r.x, r.y), xytext=(tx, ty), fontsize=5.2, fontstyle="italic", fontweight="bold", color=col, ha="left",
                            va="center", arrowprops=dict(arrowstyle="-", lw=0.25, color="#999"))
        ax.axvline(0, color="#999", lw=0.4, ls=":"); style(ax)
        ax.set_xlabel("log$_2$FC (content-adjusted)\ndistal ← mac bins → SecB niche", fontsize=6, labelpad=1); ax.set_ylabel("−log$_{10}$ adj. p (bins)", fontsize=6, labelpad=1)
    REG["d"] = (0.0, 2.1, 2.5, 3.95)
    if want('e'):
        letter(3.1, Y2, "e")
        DM = pd.read_csv(f"{O}/atlas_subtype_dotmatrix.csv")
        GD = ["SPP1", "MMP9", "VEGFA", "SLC2A1", "NDRG1", "LDHA", "CD44", "CTSL", "TGM2", "C1QC", "MRC1", "A2M"]
        CTD = [("Hypoxic macrophage", "Hypoxic"), ("Cycling macrophage", "Cycling"), ("Monocyte-derived macrophage", "Mono-derived"),
               ("C1Q tissue-resident macrophage", "C1Q resident"), ("Cycling C1Q+ tissue-resident macrophage", "Cycling C1Q+"), ("Classical monocyte", "Monocyte")]
        ax = ax_in(3.35, Y2 + 0.4, 2.3, 1.02)
        for j, g in enumerate(GD):
            d = DM[DM.gene == g].set_index("celltype"); vv = d["mean"]; sc_ = (vv - vv.min()) / (vv.max() - vv.min() + 1e-9)
            for i, (ct, _) in enumerate(CTD): ax.scatter(j, i, s=1 + 28 * d.loc[ct, "pct"] / 100, c=[GRAD(0.15 + 0.85 * sc_[ct])], edgecolors="#555", linewidths=0.2)
        XRm = XR[XR.comparison == "mac_niche_vs_ctrl"].set_index("feature"); yx = len(CTD) + 0.5
        nXd = TwoSlopeNorm(0, -1.5, 1.5)
        for j, g in enumerate(GD):
            if g not in XRm.index: ax.text(j, yx, "–", ha="center", va="center", fontsize=5, color="#aaa"); continue
            r = XRm.loc[g]; nd = r.n_up if r.mean_effect > 0 else r.n_sec - r.n_up
            ax.scatter(j, yx, s=3 + 30 * (nd / r.n_sec) ** 2, c=[DIV(nXd(np.clip(r.mean_effect, -1.5, 1.5)))], edgecolors="k" if r.p < 0.05 else "none", linewidths=0.5)
        ax.axhline(len(CTD) - 0.05, color="#999", lw=0.4); ax.axvline(8.5, color="#bbb", lw=0.4)
        ax.set_xticks(range(len(GD))); ax.set_xticklabels(GD, rotation=45, ha="right", rotation_mode="anchor", fontsize=5.8, fontstyle="italic")
        ax.set_yticks(list(range(len(CTD))) + [yx]); ax.set_yticklabels([c[1] for c in CTD] + ["Xenium niche\nvs control"], fontsize=5.8)
        ax.get_yticklabels()[0].set_color(HYP); ax.get_yticklabels()[0].set_fontweight("bold")
        ax.set_xlim(-0.6, len(GD) - 0.4); ax.set_ylim(len(CTD) + 1.0, -0.6); ax.tick_params(length=0, pad=1)
        for s_ in ax.spines.values(): s_.set_visible(False)
        cb_ = cbar_num(ax_in(3.4, Y2 + 0.13, 0.6, 0.04), Normalize(0, 1), top=True); cb_.set_ticks([0, 1]); cb_.set_ticklabels(["low", "high"])
        fig.text(3.7 / W, 1 - (Y2 + 0.2) / H, "mean expression (scaled)", fontsize=4.8, ha="center", va="top")
        lk = ax_in(4.2, Y2 + 0.05, 0.6, 0.16); lk.axis("off"); lk.set_xlim(-0.5, 2.5); lk.set_ylim(-1, 1)
        for kk, pc in enumerate([25, 50, 100]): lk.scatter(kk, 0.3, s=1 + 28 * pc / 100, c="#999", clip_on=False); lk.text(kk, -1.2, f"{pc}%", ha="center", fontsize=4.8)
        cbar_num(ax_in(4.9, Y2 + 0.13, 0.55, 0.04), nXd, cmap=DIV, top=True)
        fig.text(5.175 / W, 1 - (Y2 + 0.2) / H, "Xenium log$_2$FC", fontsize=4.8, ha="center", va="top")
    REG["e"] = (2.5, 2.1, 5.75, 3.95)
    if want('f'):
        letter(5.6, Y2, "f")
        U = pd.read_csv(f"{O}/fig6_colocal_units.csv"); wt, tm = U[U.cohort == "WT"], U[U.cohort == "TMA"]
        a1 = ax_in(6.22, Y2 + 0.4, 0.45, 1.05); a2 = ax_in(6.98, Y2 + 0.4, 0.45, 1.05); cats = ["secb_resident", "secb_hypoxic"]; cols = [RES, HYP]
        for row in wt[cats].values: a1.plot([1, 2], row, color="grey", lw=0.5, alpha=0.6)
        for j, c in enumerate(cats): a1.scatter(np.full(len(wt), j + 1), wt[c], s=9, color=cols[j], edgecolor="white", lw=0.3, zorder=3)
        p1 = wilcoxon(wt[cats[0]], wt[cats[1]]).pvalue; a1.set_xticks([1, 2]); a1.set_xticklabels(["Res.", "Hyp."], fontsize=5.5); a1.set_xlim(0.5, 2.5); style(a1)
        a1.set_ylabel("SecB within 50 µm\n(% of epithelium)", fontsize=6, labelpad=1)
        bracket(a1, 1, 2, wt[cats].values.max() * 1.02, stars(p1), p1 < 0.05)
        bp = a2.boxplot([tm[c] for c in cats], positions=[1, 2], widths=0.55, showfliers=False, patch_artist=True, medianprops=dict(color="k", lw=0.8),
                        whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5), boxprops=dict(lw=0.5))
        for pt, c in zip(bp["boxes"], cols): pt.set_facecolor(c); pt.set_alpha(0.55)
        p2 = wilcoxon(tm[cats[0]], tm[cats[1]]).pvalue; a2.set_xticks([1, 2]); a2.set_xticklabels(["Res.", "Hyp."], fontsize=5.5); a2.set_xlim(0.4, 2.6); style(a2)
        top_ = max(c.get_ydata().max() for c in bp["caps"]); bracket(a2, 1, 2, top_ * 1.02, stars(p2), p2 < 0.05)
    REG["f"] = (5.75, 2.1, 7.5, 3.95)

    # ======================================================================= row 3: g matrix
    Y3 = 3.95
    if want('g'):
        letter(0.0, Y3, "g")
        VC2 = pd.read_csv(f"{O}/visium_candidates_2sections.csv"); AT = pd.read_csv(f"{O}/atlas_hypoxic_vs_other_paired.csv")
        D2 = pd.read_csv(f"{R}/figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv").set_index("gene"); ATX = pd.read_csv(f"{O}/fig6_concordance_table.csv").set_index("gene")
        BEATS = [("Survival", [("VEGFA", "M"), ("BIRC3", "S"), ("NFKB2", "S")]),
                 ("Self-protection &\nimmunoregulation", [("CD55", "S"), ("CD47", "S"), ("LGALS9", "S"), ("IL32", "S")]),
                 ("Glycolysis & lactate", [("LDHA", "B"), ("SLC2A1", "B"), ("NDRG1", "B"), ("SLC16A3", "S")]),
                 ("Adhesion & matrix", [("SPP1", "M"), ("CD44", "M"), ("ICAM1", "S"), ("CTSL", "B"), ("MMP7", "S"), ("TGM2", "B"), ("TIMP1", "S")])]
        ROWS = [("secb", "Xenium, SecB at mac contact", "xc")] if SUPP else \
               [("mac", "Xenium", "xm"), ("mac", "Visium", "vm"), ("mac", "Atlas", "am"),
                ("contact", "Xenium", "xc"),   # SecB <=30 um of a niche macrophage vs SecB with no macrophage within 30 um
                ("secb", "Xenium", "xs"), ("secb", "Visium", "vs"), ("secb", "Atlas", "as")]
        MI = [i for i, r in enumerate(ROWS) if r[0] == "mac"]; SI = [i for i, r in enumerate(ROWS) if r[0] == "secb"]; CI = [i for i, r in enumerate(ROWS) if r[0] == "contact"]
        genes, xs, xpos = [], [], 0
        for bt, gl in BEATS:
            xs.append((bt, xpos, xpos + len(gl) - 1))
            for g, side in gl: genes.append((g, side, xpos)); xpos += 1
            xpos += 0.6
        MY = Y3 + 0.42; MH = 0.93 * len(ROWS) / 6; ax = ax_in(1.3, MY, 5.35, MH)
        nX = TwoSlopeNorm(0, -1.5, 1.5); nV = TwoSlopeNorm(0, -0.6, 0.6); nA = TwoSlopeNorm(0, -0.6, 0.6)
        def dash(x, y): ax.text(x, y, "–", ha="center", va="center", fontsize=5, color="#bbb")
        for g, side, xg in genes:
            for blk, rr in (("M", MI), ("C", CI), ("S", SI)):
                if not rr: continue
                y0, y1 = min(rr), max(rr)
                if side in ({"C": "S"}.get(blk, blk), "B"): ax.add_patch(Rectangle((xg - 0.47, y0 - 0.47), 0.94, y1 - y0 + 0.94, color="#F1ECE3", lw=0, zorder=0))
        # size = % of cells (bins) expressing the gene in the focal group; colour = effect; ring = significant
        XP = pd.read_csv(f"{O}/xenium_pct_summary.csv").pivot(index="gene", columns="group", values="pct")
        EP = pd.read_csv(f"{O}/atlas_epi_pct.csv").set_index("gene").pct_SecB
        Gp = pd.read_csv(f"{O}/atlas_pb/groups.csv"); gnm = open(f"{O}/atlas_pb/genes.txt").read().split("\n"); NZp = np.load(f"{O}/atlas_pb/nnz.npy")
        hk = (Gp.celltype_level2 == "Hypoxic macrophage").values; AMP = {g_: 100 * NZp[hk, gnm.index(g_)].sum() / Gp.n_cells[hk].sum() for g_, _, _ in genes if g_ in gnm}
        DOTS = []
        def dot(xg, i, pct, eff, norm, sig): DOTS.append((xg, i, pct, eff, sig))
        for i, (blk, lab, kind) in enumerate(ROWS):
            for g, side, xg in genes:
                if kind in ("xm", "xs", "xc"):
                    r = XR[(XR.comparison == {"xm": "mac_niche_vs_ctrl", "xs": "secb_vs_seca", "xc": "secb_contact"}[kind]) & (XR.feature == g)]
                    col = {"xm": "mac_niche", "xs": "secb", "xc": "secb_contact"}[kind]
                    if not len(r) or g not in XP.index: dash(xg, i); continue
                    r = r.iloc[0]; dot(xg, i, XP.loc[g, col], r.mean_effect, nX, r.p < 0.05)
                elif kind in ("vm", "vs"):
                    sd = "mac" if kind == "vm" else "secb"; r = VC2[(VC2.side == sd) & (VC2["sample"] == "WT01") & (VC2.gene == g)]
                    if not len(r): dash(xg, i); continue
                    r = r.iloc[0]
                    if sd == "mac" and r.effect > 0 and not bool(r.mac_sourced): ax.text(xg, i, "×", ha="center", va="center", fontsize=6, color=SECB); continue
                    dot(xg, i, 100 * r.det, r.effect, nV, r.p < 0.05 and abs(r.effect) >= 0.1)
                elif kind == "am":
                    r = AT[AT.feature == g]
                    if g not in AMP: dash(xg, i); continue
                    if len(r): r = r.iloc[0]; eff, sg = r.mean_effect, (r.p < 0.05 and abs(r.mean_effect) >= 0.05)
                    elif g in ATX.index and isinstance(ATX.loc[g, "atlas_mac_v"], str): eff, sg = float(ATX.loc[g, "atlas_mac_v"].split()[0]), ATX.loc[g, "atlas_mac"] in ("up", "down")
                    else: dash(xg, i); continue
                    dot(xg, i, AMP[g], eff, nA, sg)
                else:
                    if g not in D2.index or g not in EP.index: dash(xg, i); continue
                    r = D2.loc[g]; dot(xg, i, EP[g], r.log2FoldChange, nX, r.padj < 0.05 and abs(r.log2FoldChange) >= 0.25)
        nZ = TwoSlopeNorm(0, -2, 2); DT = pd.DataFrame(DOTS, columns=["x", "i", "pct", "eff", "sig"])
        DT["z"] = DT.eff / DT.groupby("i").eff.transform(lambda v: np.sqrt(np.mean(v ** 2)) + 1e-9)   # standardized per comparison (RMS = 1), 0 kept at 0
        for _, r in DT.iterrows():
            ax.scatter(r.x, r.i, s=2 + 44 * np.clip(r.pct, 0, 100) / 100, c=[DIV(nZ(np.clip(r.z, -2, 2)))], edgecolors="k" if r.sig else "none", linewidths=0.6, zorder=3)
        DT.to_csv(f"{O}/fig6h_matrix_standardized.csv", index=False)
        ax.set_xlim(-0.6, xpos - 0.4); ax.set_ylim(len(ROWS) - 0.5, -0.55); ax.set_xticks([x for _, _, x in genes])
        ax.set_xticklabels([g for g, _, _ in genes], fontsize=5.8, fontstyle="italic", rotation=45, ha="right", rotation_mode="anchor")   # genes along the bottom
        ax.tick_params(axis="x", length=0, pad=1.5); ax.set_yticks(range(len(ROWS))); ax.set_yticklabels([r[1] for r in ROWS], fontsize=5.8); ax.tick_params(axis="y", length=0, pad=2)
        for s_ in ax.spines.values(): s_.set_visible(False)
        if MI and SI: ax.axhline(max(MI) + 0.5, color="#999", lw=0.5)
        if CI: ax.axhline(max(CI) + 0.5, color="#999", lw=0.5)
        for bt, a_, b_ in xs:
            ax.text((a_ + b_) / 2, -0.95, bt, ha="center", va="bottom", fontsize=6, fontweight="bold"); ax.plot([a_ - 0.3, b_ + 0.3], [-0.75, -0.75], color="k", lw=0.6, clip_on=False)
        for lab_, rr, col_ in [("Macro-\nphage", MI, NICHE), ("SecB at mac\ncontact", CI, SECB), ("SecB\nepithelium", SI, SECB)]:
            if not rr or SUPP: continue
            y0, y1 = min(rr), max(rr)
            top_in = MY + (y0 - 0.5 + 0.55) / (len(ROWS) + 0.05) * MH; bot_in = MY + (y1 + 0.5 + 0.55) / (len(ROWS) + 0.05) * MH
            fig.lines.append(plt.Line2D([0.88 / W, 0.88 / W], [1 - (top_in + 0.02) / H, 1 - (bot_in - 0.02) / H], transform=fig.transFigure, color=col_, lw=1.2,
                                        ls=(0, (1.5, 1)) if rr is CI else "-"))
            one = len(rr) == 1   # single-row block: label horizontal (a rotated label would overrun the row)
            fig.text(0.84 / W, 1 - ((top_in + bot_in) / 2) / H, lab_, ha="right", va="center", rotation=0 if one else 90, multialignment="right" if one else "center",
                     fontsize=5.2 if one else 5.8, fontweight="bold", color=col_, linespacing=1.0)
        ly = MY
        cax_ = ax_in(6.85, ly + 0.05, 0.05, min(MH * 0.55, 0.55)); cz = fig.colorbar(plt.cm.ScalarMappable(nZ, DIV), cax=cax_, orientation="vertical")
        cz.outline.set_linewidth(0.4); cz.set_ticks([-2, 0, 2]); cz.ax.tick_params(labelsize=5, width=0.4, length=1.5, pad=1)
        fig.text(6.72 / W, 1 - (ly - 0.02) / H, "standardized\neffect (z)", fontsize=5, ha="left", va="bottom")
        lk = ax_in(6.75, ly + min(MH * 0.55, 0.55) + 0.2, 0.2, 0.55); lk.axis("off"); lk.set_xlim(-1, 1); lk.set_ylim(2.6, -0.6)
        for kk, pc in enumerate([25, 50, 100]):
            lk.scatter(0, kk, s=2 + 44 * pc / 100, c="#999", clip_on=False); lk.text(0.9, kk, f"{pc}%", ha="left", va="center", fontsize=4.8)
        lk.text(-0.5, -0.75, "% expressing", ha="left", va="bottom", fontsize=5)

    REG["g"] = (0.0, 3.95, 7.5, 5.9)

    # ======================================================================= ROI rows (Xenium)
    def roi_row(y, letter_, tiles, roi, start_col=0, ncol=5, caption=None):
        D = RS.load(roi["sample"]); x0, y0, S = roi["x"], roi["y"], roi["size"]
        idx = RS.roi_idx(D, x0, y0, S); d = D["d"].iloc[idx]; pl = RS.polys(D, idx); lab = d.cell_label_singler.values
        letter(0.0, y, letter_); n = ncol; TW = 1.2; gap = (7.3 - n * TW) / max(n - 1, 1)

        for j, (t, title) in enumerate(tiles):
            xx = 0.12 + (j + start_col) * (TW + gap)
            TILEREG[(letter_, j)] = ((xx - 0.2, y + 0.04, xx + TW + 0.2, y + 0.22 + TW + 0.24), t.split("@")[0].replace("*", "x"))
            if TILE_ONLY is not None and TILE_ONLY != j: continue
            ax = ax_in(xx, y + 0.22, TW, TW); sub = ""
            if t == "celltype":
                faces = [to_rgba(PAL.get(l, GREY)) for l in lab]; order = np.arange(len(pl)); edge = "#8A8A8A"
            else:
                gpart, _, tk = t.partition("@"); types = RS.TYPES.get(tk or "all"); show = np.ones(len(lab), bool) if types is None else np.isin(lab, types)
                raw = RS.gene_vals(D, idx, gpart); lo, hi = np.nanpercentile(raw[show], [2, 98]); nrm_ = Normalize(lo, hi if hi > lo else lo + 1)
                faces = [GRAD(nrm_(v_)) if sh else to_rgba(GREY) for v_, sh in zip(raw, show)]; order = np.lexsort((np.where(show, raw, -1),)); edge = GREY_EDGE
                cbar_num(ax_in(xx + TW * 0.56, y + 0.22 - 0.075, TW * 0.44, 0.045), nrm_, top=True)      # Fig 5: scale top right
            ax.add_collection(PolyCollection([pl[i] for i in order], facecolors=[faces[i] for i in order], edgecolors=edge, linewidths=0.05))
            ax.set_xlim(x0 - S / 2, x0 + S / 2); ax.set_ylim(y0 - S / 2, y0 + S / 2); ax.set_aspect("equal"); ax.axis("off")
            if t != "celltype": ax.set_title(title, fontsize=FA, fontweight="bold", fontstyle="italic", pad=3, loc="left")   # Fig 5: gene top left
            if j == 0:
                L_ = 100; xs_ = x0 + S / 2 - S * 0.05 - L_; ys_ = y0 - S / 2 + S * 0.04
                ax.plot([xs_, xs_ + L_], [ys_, ys_], color="k", lw=1.2, solid_capstyle="butt"); ax.text(xs_ + L_ / 2, ys_ + S * 0.025, "100 µm", ha="center", fontsize=5)
    Y4 = 5.92
    if want('h'):
        roi_row(Y4, "h", [("celltype", "Cell type"), ("VEGFA@mac", "VEGFA"), ("CTSL@mac", "CTSL"), ("CD44@mac", "CD44"), ("ICAM1@mac", "ICAM1")], ROI,
                caption=f"Xenium {ROI['sample']}, {ROI['size']:g} µm ROI; macrophage genes on macrophages (others grey)")
    REG["h"] = (0.0, 5.9, 7.5, 7.78)

    Y6 = 7.8
    if want('k'):
        roi_row(Y6, "k", [("CD55@epi", "CD55"), ("IL32@epi", "IL32"), ("MMP7@epi", "MMP7"), ("SLC16A3@epi", "MCT-4")], ROI, start_col=1,
                caption="Same ROI; SecB / epithelial genes on epithelium (others grey)")
    REG["k"] = (0.0, 7.78, 7.5, 9.7)



def page():
    global fig
    fig = plt.figure(figsize=(W, H)); return fig
# the full page is drawn once, not saved: it records every panel's crop region (REG)
page(); draw()
plt.close(fig); regions = dict(REG)
NAMES = {"a": "Fig6a_visium_section_SecA_SecB", "b": "Fig6b_visium_ROI_SecA_SecB", "m": "Fig6d_visium_ROI_macrophage_niche",
         "c": "Fig6c_visium_ROI_CA9", "d": "Fig6e_visium_niche_macrophage_volcano", "e": "Fig6f_atlas_hypoxic_macrophage_dotplot",
         "f": "Fig6g_hypoxic_macrophage_SecB_neighbours", "g": "Fig6h_crossplatform_matrix", "h": "Fig6i_ROI_macrophage_VEGFA_CTSL_CD44_ICAM1",
         "k": "Fig6j_ROI_SecB_CD55_IL32_MMP7_MCT4"}
SHOW_LETTERS = False
TRANSPARENT = {"h", "k"}   # the two ROI rows of panel 6i (files Fig6i*, Fig6j*): transparent-background PNGs
for k, (x0, y0, x1, y1) in regions.items():
    page(); draw(only=k); bb = Bbox([[x0, H - y1], [x1, H - y0]])
    for ext in ("png", "pdf", "svg"):
        tr_ = (k in TRANSPARENT and ext == "png")
        fig.savefig(f"{PAN}/{NAMES[k]}.{ext}", dpi=450 if ext == "png" else None, facecolor="none" if tr_ else "white", transparent=tr_, bbox_inches=bb)
    plt.close(fig)
# each ROI tile (keys h and k, the two rows of panel 6i) also as its own file, drawn alone
for (lt, j), (bx_, gname) in sorted(TILEREG.items()):
    TILE_ONLY = j; page(); draw(only=lt); x0, y0, x1, y1 = bx_; bb = Bbox([[x0, H - y1], [x1, H - y0]])
    nm_ = f"{NAMES[lt].split('_')[0]}{j + 1}_ROI_{'MCT4_SLC16A3' if gname == 'SLC16A3' else gname}"
    for ext in ("png", "pdf", "svg"):
        tr_ = (lt in TRANSPARENT and ext == "png")
        fig.savefig(f"{PAN}/{nm_}.{ext}", dpi=450 if ext == "png" else None, facecolor="none" if tr_ else "white", transparent=tr_, bbox_inches=bb)
    plt.close(fig)
TILE_ONLY = None
print("saved panels", sorted(regions), "+ ROI tiles", len(TILEREG), "| ROI:", ROI["name"], ROI["sample"], ROI["size"])
