#!/usr/bin/env python3
"""Figures S16 / S17 / S18 -- per-tissue GAM small multiples.

PURPOSE
    Plots the output of 21_figS16-18_per_sample_gams.R, whose internal figure ids SF14 / SF15 / SF16 are
    Figures S16 / S17 / S18 (this script writes the final numbers). Thin coloured lines = individual whole
    tissues (WT01-WT08; solid = treatment-naive, dashed = post-chemotherapy), black dashed = pooled GAM
    (identical to the main-figure fit) with its 95% CI. Corner: n tissues whose curve moves in the pooled
    direction / 8, and the paired Wilcoxon p (top vs bottom 20% of x within tissue). The tissue legend is
    drawn inside each figure.

INPUTS
    OUTPUT_ROOT/supplementary/sf14_16_per_sample_gams/{sample_meta, curves_pooled, curves_per_sample,
    per_sample_consistency, pooled_vs_main}.csv (from 21_figS16-18_per_sample_gams.R)

OUTPUTS
    <out_dir>/SuppFig16_per_sample_gams_environment.{svg,png}, SuppFig17_per_sample_gams_epithelial.{svg,png},
    SuppFig18_per_sample_gams_macrophage.{svg,png} (out_dir = argv[1], default .../sf14_16_per_sample_gams/figures)
    .../sf14_16_per_sample_gams/SuppFig16-18_per_sample_consistency.csv (per-panel statistics)

MANUSCRIPT PANEL(S)
    Figures S16a-d, S17a-n, S18a-h.

RUNTIME TIER
    fast (reads the curve tables).
"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R = os.environ.get("OUTPUT_ROOT", "output")
IN = f"{R}/supplementary/sf14_16_per_sample_gams"
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{IN}/figures"
os.makedirs(OUT, exist_ok=True)

FA, FK, FT, FS, FP = 6.5, 5.5, 6.5, 5.2, 8      # axis, tick, title, stat, panel letter (pt)
LW = 0.5
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": LW,
                     "xtick.major.width": LW, "ytick.major.width": LW, "xtick.major.size": 2.2,
                     "ytick.major.size": 2.2, "xtick.labelsize": FK, "ytick.labelsize": FK})

meta = pd.read_csv(f"{IN}/sample_meta.csv").sort_values("code")
SCOL = dict(zip(meta.code, ["#D14E6C", "#E07850", "#E6A141", "#8A5DAF",
                            "#5665B6", "#56AFC4", "#4E9C81", "#7D4E4E"]))
SLS = {c: ("-" if t.lower().startswith("treatment") else (0, (3.0, 1.6))) for c, t in zip(meta.code, meta.treatment)}
TREAT = dict(zip(meta.code, meta.treatment))

pool = pd.read_csv(f"{IN}/curves_pooled.csv")
samp = pd.read_csv(f"{IN}/curves_per_sample.csv")
st = pd.read_csv(f"{IN}/per_sample_consistency.csv")

# Figure S17 (internal SF15) shows only genes named in the text or main figures: TGFB1 (text names only the
# TGFβ pathway), INHBA, C7, CDH2 and ITGB5, also fitted by 21_..., are dropped here and the panels re-lettered.
DROP = {"SF15": {"TGFB1", "INHBA", "C7", "CDH2", "ITGB5"}}
for fid, genes in DROP.items():
    keep = [pn for pn, ft in pool[pool.figure == fid][["panel", "feature"]].drop_duplicates().itertuples(index=False) if ft not in genes]
    relab = {old: chr(65 + i) for i, old in enumerate(keep)}
    def _apply(df):
        m = df.figure == fid
        df = df[~(m & df.feature.isin(genes))].copy()
        m = df.figure == fid; df.loc[m, "panel"] = df.loc[m, "panel"].map(relab); return df
    pool, samp, st = _apply(pool), _apply(samp), _apply(st)

FIGS = {
    "SF14": ("SuppFig16_per_sample_gams_environment", "Focal polarization (SecA → SecB)", "Feature (z)"),
    "SF15": ("SuppFig17_per_sample_gams_epithelial", "Focal polarization (SecA → SecB)", "Epithelial expression (z)"),
    "SF16": ("SuppFig18_per_sample_gams_macrophage", "Local epithelial polarization\n(SecA → SecB, z)",
             "Macrophage expression (z)"),
}
SUFFIX = {  # value type shown under the title
    "Cell density": "cells within 50 µm", "Distance to nearest vascular cell (µm)": "µm",
    "Matrix remodelling": "neighbourhood UCell", "Hypoxia": "neighbourhood UCell",
    "Glycolytic": "UCell", "NF-κB activated": "UCell"}
UCELL_EPI = {"Proliferation", "Myc", "PI3K", "Wnt", "DNA damage repair", "NF-κB", "Glycolysis"}


def fmt_p(p):
    return "p < 0.0001" if p < 1e-4 else (f"p = {p:.4f}" if p < 1e-3 else f"p = {p:.3f}")


def legend_handles():
    h = [Line2D([], [], color=SCOL[c], lw=1.1, ls=SLS[c],
                label=f"{c} ({'naive' if TREAT[c].lower().startswith('treatment') else 'post-chemo'})")
         for c in meta.code]
    h.append(Line2D([], [], color="black", lw=1.5, ls=(0, (3.5, 1.8)), label="Pooled GAM (95% CI)"))
    return h


def draw(fig_id):
    stem, xlab, ylab = FIGS[fig_id]
    P = pool[pool.figure == fig_id]; S = samp[samp.figure == fig_id]; T = st[st.figure == fig_id]
    panels = list(dict.fromkeys(P.panel))
    ncol = 4; nrow = int(np.ceil(len(panels) / ncol))
    pw, ph = 1.72, 1.42
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * pw, nrow * ph + 0.55), squeeze=False)
    for k, ax in enumerate(axes.flat):
        if k >= len(panels):
            ax.axis("off"); continue
        pan = panels[k]; p = P[P.panel == pan]; feat = p.feature.iloc[0]
        ax.axhline(0, ls=":", color="#999999", lw=0.4, zorder=0)
        for c in meta.code:
            s = S[(S.panel == pan) & (S.code == c)]
            ax.plot(s.x, s.fitted, color=SCOL[c], lw=0.8, ls=SLS[c], alpha=0.95, zorder=2)
        ax.fill_between(p.x, p.lower, p.upper, color="black", alpha=0.12, lw=0, zorder=3)
        ax.plot(p.x, p.fitted, color="black", lw=1.5, ls=(0, (3.5, 1.8)), zorder=4)
        t = T[T.panel == pan].iloc[0]
        up = t.pooled_delta > 0; kk = t.k_rise if up else t.k_fall
        sig = t.wilcoxon_p < 0.05
        ax.text(0.03, 0.97, f"{kk}/8 {'rise' if up else 'fall'}\n{fmt_p(t.wilcoxon_p)}", transform=ax.transAxes,
                ha="left", va="top", fontsize=FS, color="black" if sig else "#666666",
                fontweight="bold" if sig else "normal", zorder=6,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.4))
        sub = SUFFIX.get(feat, "UCell" if feat in UCELL_EPI else "log1p CP10k")
        ttl = "Distance to nearest\nvascular cell (µm)" if feat.startswith("Distance") else (
              f"{feat} ({sub})" if len(feat) + len(sub) < 24 else f"{feat}\n({sub})")
        ax.set_title(ttl,
                     fontsize=FT, pad=2.5, fontweight="bold")
        ax.text(-0.2, 1.08, pan.lower(), transform=ax.transAxes, fontsize=FP, fontweight="bold", va="bottom", ha="left")
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.tick_params(pad=1.2)
        ymin = min(S[S.panel == pan].fitted.min(), p.lower.min()); ymax = max(S[S.panel == pan].fitted.max(), p.upper.max())
        pad = 0.06 * (ymax - ymin); ax.set_ylim(ymin - pad, ymax + pad + 0.22 * (ymax - ymin))   # headroom for stats text
        if k % ncol == 0: ax.set_ylabel(ylab, fontsize=FA, labelpad=1.5)
        if k >= len(panels) - ncol: ax.set_xlabel(xlab, fontsize=FA, labelpad=1.5)
    # legend: in the first empty slot if there is one, otherwise underneath
    empty = [ax for k, ax in enumerate(axes.flat) if k >= len(panels)]
    if empty:
        empty[0].legend(handles=legend_handles(), loc="center left", frameon=False, fontsize=FK, handlelength=2.6,
                        borderaxespad=0, labelspacing=0.45)
        fig.tight_layout(pad=0.4, h_pad=0.9, w_pad=0.7)
    else:
        H = fig.get_size_inches()[1]
        fig.tight_layout(pad=0.4, h_pad=0.9, w_pad=0.7, rect=(0, 0.42 / H, 1, 1))
        fig.legend(handles=legend_handles(), loc="lower center", ncol=5, frameon=False, fontsize=FK,
                   handlelength=2.6, columnspacing=1.2, bbox_to_anchor=(0.5, 0.0))
    for ext in ("svg", "png"):
        fig.savefig(f"{OUT}/{stem}.{ext}", dpi=450, bbox_inches="tight", facecolor="white")
    plt.close(fig); print("saved", stem, len(panels), "panels")


for f in FIGS: draw(f)

# ---- final stats table (-> final_panels/supplementary/tables/SuppFig16-18_per_sample_consistency.csv after review)
c = pd.read_csv(f"{IN}/pooled_vs_main.csv")
MAINFIG = {"SF14": "Fig 4f", "SF15": "Fig 5a", "SF16": "Fig 5h (fig5e_mac_gam_curves nb_pol)"}
c["main_figure_curve"] = c.figure.map(MAINFIG) + ": " + c.main_feature
t = st.merge(c[["figure", "feature", "main_figure_curve", "max_abs_dfit"]]
             .rename(columns={"max_abs_dfit": "pooled_vs_main_max_abs_diff"}), on=["figure", "feature"], how="left")
t["pooled_direction"] = np.where(t.pooled_delta > 0, "rise", "fall")
t["k_consistent"] = np.where(t.pooled_delta > 0, t.k_rise, t.k_fall)
t["figure"] = t.figure.map({"SF14": "SuppFig16", "SF15": "SuppFig17", "SF16": "SuppFig18"})   # internal id -> published figure
t = t[["figure", "panel", "feature", "source", "x_axis", "n_cells", "pooled_left", "pooled_right", "pooled_delta",
       "pooled_direction", "k_rise", "k_fall", "k_consistent", "k_top_gt_bottom", "median_top_minus_bottom",
       "wilcoxon_p", "per_sample_delta", "main_figure_curve", "pooled_vs_main_max_abs_diff"]]
for k in ["pooled_left", "pooled_right", "pooled_delta", "median_top_minus_bottom"]: t[k] = t[k].round(4)
t.to_csv(f"{IN}/SuppFig16-18_per_sample_consistency.csv", index=False); print("saved stats table")
