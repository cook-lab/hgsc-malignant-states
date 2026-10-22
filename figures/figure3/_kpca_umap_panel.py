"""
Shared layout of the KPCA.B epithelial UMAP panels, Fig 3j (polarization) and Fig 3k (Krt7)
===========================================================================================
The 1.8 in UMAP square of Fig 3i (models/kpca_mouse/scripts/41b_export_1p8in.py: same axes box, 4% padding,
L-shaped axis stubs, bold title and italic subtitle top left, point size 1.4, alpha 0.95), with the colour bar as a
horizontal bar centred UNDER the square. The canvas is 1.8 in wide and 1.8 in + CB_STRIP tall, so the UMAP squares
line up with Fig 3i.

Used by figures/figure3/04_fig3j_kpca_polarization_7plus7.py (3j) and 05_fig3k_kpca_krt7_umap.py (3k).
"""
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

SIDE = 1.8        # in, UMAP square (as Fig 3i)
CB_STRIP = 0.34   # in, strip under the square that holds the colour bar and its labels
CB_W, CB_H = 1.0, 0.065   # in, colour bar
CB_Y = 0.17       # in, bottom of the colour bar above the canvas edge
FA, FK = 8, 6

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": FA, "svg.fonttype": "none", "pdf.fonttype": 42, "figure.dpi": 450,
                     "savefig.dpi": 450})


def epi_umap_panel(ex, ey, vals, order, cmap, norm, title, subtitle, ticks, ticklabels, stem, end_labels=None):
    """Draw cells (ex, ey) coloured by vals, in the plotting order `order`, and save stem.{svg,png}.
    end_labels = (left, right) text beside the ends of the colour bar, e.g. ("SecA", "SecB")."""
    W, H = SIDE, SIDE + CB_STRIP
    fig = plt.figure(figsize=(W, H))
    # same box as the 1.8 x 1.8 in Fig 3i panel (subplots_adjust left 0.03, right 0.99, bottom 0.03, top 0.99)
    ax = fig.add_axes([0.03, (CB_STRIP + 0.03 * SIDE) / H, 0.96, 0.96 * SIDE / H])
    ax.scatter(ex[order], ey[order], c=vals[order], cmap=cmap, norm=norm, s=1.4, marker="o", linewidths=0,
               rasterized=True, alpha=0.95)
    exr, eyr = ex.max() - ex.min(), ey.max() - ey.min()
    epx, epy = 0.04 * exr, 0.04 * eyr
    ax.set_xlim(ex.min() - epx, ex.max() + epx); ax.set_ylim(ey.min() - epy, ey.max() + epy)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_bounds(ex.min() - epx, ex.min() - epx + 0.18 * exr)
    ax.spines["left"].set_bounds(ey.min() - epy, ey.min() - epy + 0.18 * eyr)
    ax.spines["bottom"].set_linewidth(0.4); ax.spines["left"].set_linewidth(0.4)
    ax.text(0.04, 0.97, title, transform=ax.transAxes, fontsize=FA, fontweight="bold", ha="left", va="top")
    if subtitle:
        ax.text(0.04, 0.875, subtitle, transform=ax.transAxes, fontsize=FK, fontstyle="italic", ha="left",
                va="top", color="#444444")
    cax = fig.add_axes([(W - CB_W) / 2 / W, CB_Y / H, CB_W / W, CB_H / H])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.3)
    cb.ax.tick_params(length=1.5, width=0.3, pad=1, labelsize=FK)
    cb.set_ticks(ticks); cb.set_ticklabels(ticklabels)
    if end_labels:
        yc = (CB_Y + CB_H / 2) / H
        fig.text((W - CB_W) / 2 / W - 0.03 / W, yc, end_labels[0], fontsize=FK, ha="right", va="center")
        fig.text((W + CB_W) / 2 / W + 0.03 / W, yc, end_labels[1], fontsize=FK, ha="left", va="center")
    fig.savefig(stem + ".svg", format="svg", pad_inches=0.0)
    fig.savefig(stem + ".png", format="png", pad_inches=0.0)
    plt.close(fig)
