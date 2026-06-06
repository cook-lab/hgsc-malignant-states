#!/usr/bin/env python3
"""
Atlas 08 — Xenium reference bridge (2/2): downsample + markers + coverage report

PURPOSE
    1. Downsample hgsc_atlas_xenium.h5ad to 1000 cells per xenium_celltype.
    2. Run Wilcoxon one-vs-rest DEG markers per type.
    3. Cross-reference markers against the Xenium gene panel.
    4. Produce an HTML panel-coverage report.
    The downsampled object is the canonical SingleR reference for Xenium annotation.

INPUTS
    DATA_ROOT/2026_final_atlas/hgsc_atlas_xenium.h5ad     (from add_xenium_celltype.py)
    DATA_ROOT/2026_final_atlas/output/xenium_panel_genes.txt

OUTPUTS
    obj("xenium_ref")  = output/xenium_celltype/xenium_celltype_downsampled.h5ad  (canonical copy)
    output_root/08_xenium_reference/markers/<type>.csv
    output_root/08_xenium_reference/xenium_celltype_coverage_report.html

MANUSCRIPT PANEL(S)
    Cross-platform bridge; SingleR reference for Xenium panels (Fig 4-6, SF10/SF11).

RUNTIME TIER
    moderate (backed read + 16k-cell DE).
"""

import gc
import os
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

import warnings
warnings.filterwarnings("ignore")
np.random.seed(SEED)

H5AD_PATH  = path("data_root", "2026_final_atlas", "hgsc_atlas_xenium.h5ad")
PANEL_FILE = path("data_root", "2026_final_atlas", "output", "xenium_panel_genes.txt")
OUT_DIR    = path("output_root", "08_xenium_reference")
MARKER_DIR = os.path.join(OUT_DIR, "markers")
OUT_HTML   = os.path.join(OUT_DIR, "xenium_celltype_coverage_report.html")
OUT_H5AD   = obj("xenium_ref")
os.makedirs(MARKER_DIR, exist_ok=True)
os.makedirs(os.path.dirname(OUT_H5AD), exist_ok=True)

MAX_CELLS = 1000

LABEL_ORDER = [
    "Secretory epithelium", "Ciliated epithelium", "Mesothelial", "Fibroblast",
    "Smooth muscle", "Pericyte", "Endothelial", "T cell", "NK cell", "B cell",
    "Plasma cell", "Macrophage", "Conventional dendritic cell",
    "Plasmacytoid dendritic cell", "Neutrophil", "Mast cell",
]

# ── load xenium panel ──────────────────────────────────────────────────
panel_genes = set()
with open(PANEL_FILE) as f:
    for line in f:
        g = line.strip()
        if g:
            panel_genes.add(g)
print(f"Xenium panel: {len(panel_genes)} genes")

# ── load & downsample ─────────────────────────────────────────────────
print(f"\nLoading {H5AD_PATH}...", flush=True)
adata = sc.read_h5ad(H5AD_PATH, backed="r")
print(f"  Shape: {adata.shape[0]:,} × {adata.shape[1]:,}")

obs = adata.obs[["xenium_celltype"]].copy()
obs["xenium_celltype"] = obs["xenium_celltype"].astype(str)

print(f"\nDownsampling to {MAX_CELLS:,} cells per type...", flush=True)
rng = np.random.default_rng(SEED)
keep_idx = []
for ct in LABEL_ORDER:
    ct_idx = obs.index[obs["xenium_celltype"] == ct].tolist()
    n_take = min(MAX_CELLS, len(ct_idx))
    if len(ct_idx) > MAX_CELLS:
        ct_idx = rng.choice(ct_idx, MAX_CELLS, replace=False).tolist()
    keep_idx.extend(ct_idx)
    print(f"  {ct:<35} {n_take:>6,} cells")
print(f"  Total: {len(keep_idx):,} cells")

print("Loading subset into memory...", flush=True)
adata_sub = adata[keep_idx].to_memory()
adata.file.close()
del adata
gc.collect()
print(f"  In-memory shape: {adata_sub.shape}")

print(f"\nSaving downsampled h5ad → {OUT_H5AD}...", flush=True)
adata_sub.write_h5ad(OUT_H5AD)
print(f"  File size: {os.path.getsize(OUT_H5AD) / (1024**3):.2f} GB")

# ── compute markers ───────────────────────────────────────────────────
print("\nComputing Wilcoxon DEG markers (one-vs-rest)...", flush=True)
sc.tl.rank_genes_groups(adata_sub, groupby="xenium_celltype", method="wilcoxon",
                        key_added="rank_genes_xenium")

print("\nExtracting & saving markers...", flush=True)
result = adata_sub.uns["rank_genes_xenium"]
groups = result["names"].dtype.names

celltype_markers = {}
for ct in groups:
    df = pd.DataFrame({
        "names": [str(x) for x in result["names"][ct]],
        "scores": result["scores"][ct].astype(float),
        "logfoldchanges": result["logfoldchanges"][ct].astype(float),
        "pvals": result["pvals"][ct].astype(float),
        "pvals_adj": result["pvals_adj"][ct].astype(float),
    })
    df = df[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"] > 0)].reset_index(drop=True)
    celltype_markers[ct] = df
    safe = ct.replace("/", "-").replace(" ", "_")
    df.to_csv(os.path.join(MARKER_DIR, f"{safe}.csv"), index=False)
    print(f"  {ct}: {len(df)} markers")

del adata_sub
gc.collect()

# ═══════════════════════════════════════════════════════════════════════
#  COVERAGE ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

print("\nAnalysing panel coverage...", flush=True)
records = []
for ct in LABEL_ORDER:
    df = celltype_markers.get(ct, pd.DataFrame())
    if df.empty:
        records.append(dict(celltype=ct, total_markers=0, in_panel_all=0, pct_all=0,
                            in_panel_top20=0, pct_top20=0, in_panel_top50=0, pct_top50=0,
                            top_in=[], top_miss=[]))
        continue
    df["in_panel"] = df["names"].isin(panel_genes)
    total = len(df)
    h20 = int(df.head(20)["in_panel"].sum())
    h50 = int(df.head(50)["in_panel"].sum())
    hall = int(df["in_panel"].sum())
    records.append(dict(
        celltype=ct, total_markers=total,
        in_panel_all=hall, pct_all=round(100 * hall / total, 1) if total else 0,
        in_panel_top20=h20, pct_top20=round(100 * h20 / min(20, total), 1) if total else 0,
        in_panel_top50=h50, pct_top50=round(100 * h50 / min(50, total), 1) if total else 0,
        top_in=df.loc[df["in_panel"], "names"].head(5).tolist(),
        top_miss=df.loc[~df["in_panel"], "names"].head(5).tolist(),
    ))

summary = pd.DataFrame(records)


def risk_from_hits(n):
    if n < 3:  return "HIGH RISK"
    if n <= 5: return "MODERATE"
    return "OK"


summary["risk"] = summary["in_panel_top20"].apply(risk_from_hits)

# ═══════════════════════════════════════════════════════════════════════
#  HTML REPORT
# ═══════════════════════════════════════════════════════════════════════

def rc(risk):
    return {"HIGH RISK": "#e74c3c", "MODERATE": "#f39c12", "OK": "#27ae60"}[risk]


def pct_bar(pct, width=140):
    col = "#e74c3c" if pct < 15 else "#f39c12" if pct < 30 else "#27ae60"
    return (f'<div style="background:#eee;border-radius:4px;width:{width}px;height:18px;'
            f'display:inline-block;vertical-align:middle">'
            f'<div style="background:{col};width:{pct*width/100:.0f}px;height:18px;'
            f'border-radius:4px"></div></div>'
            f' <span style="font-size:12px">{pct:.0f}%</span>')


n_total = len(summary)
n_high  = (summary["risk"] == "HIGH RISK").sum()
n_mod   = (summary["risk"] == "MODERATE").sum()
n_ok    = (summary["risk"] == "OK").sum()

bar_items = ""
for _, row in summary.sort_values("pct_top20").iterrows():
    pct = row["pct_top20"]
    col = rc(row["risk"])
    bar_items += (
        f'<div style="display:flex;align-items:center;margin:4px 0">'
        f'<div style="width:220px;text-align:right;padding-right:10px;font-size:13px">'
        f'{row["celltype"]}</div>'
        f'<div style="background:#eee;width:320px;border-radius:4px;height:22px">'
        f'<div style="background:{col};width:{pct*3.2:.0f}px;height:22px;'
        f'border-radius:4px"></div></div>'
        f' <span style="font-size:13px;margin-left:8px;font-weight:600">'
        f'{row["in_panel_top20"]}/20</span>'
        f' <span style="font-size:12px;color:#888">({pct:.0f}%)</span></div>\n'
    )

detail_cards = []
for _, row in summary.iterrows():
    ct = row["celltype"]
    df_detail = celltype_markers.get(ct, pd.DataFrame())
    if not df_detail.empty:
        df_detail["in_panel"] = df_detail["names"].isin(panel_genes)
    top30 = df_detail.head(30)
    gene_rows = ""
    for _, g in top30.iterrows():
        badge = ('<span style="background:#27ae60;color:#fff;padding:1px 6px;'
                 'border-radius:8px;font-size:11px">IN PANEL</span>'
                 if g["in_panel"] else
                 '<span style="background:#e0e0e0;color:#666;padding:1px 6px;'
                 'border-radius:8px;font-size:11px">NOT IN PANEL</span>')
        gene_rows += (f"<tr><td>{g['names']}</td><td>{g['scores']:.1f}</td>"
                      f"<td>{g['logfoldchanges']:.2f}</td><td>{badge}</td></tr>\n")
    risk_col = rc(row["risk"])
    top_in_str = ", ".join(row["top_in"]) if row["top_in"] else "NONE"
    top_miss_str = ", ".join(row["top_miss"]) if row["top_miss"] else "—"
    detail_cards.append(f"""
    <div class="card">
      <h3 style="margin:0 0 4px 0">{ct}
        <span style="background:{risk_col};color:#fff;padding:2px 10px;border-radius:10px;
               font-size:13px;margin-left:8px">{row['risk']}</span>
      </h3>
      <div style="display:flex;gap:24px;margin-bottom:6px;font-size:13px">
        <div>Top-20 in panel: <b>{row['in_panel_top20']}/20</b> ({row['pct_top20']}%)</div>
        <div>Top-50 in panel: <b>{row['in_panel_top50']}/{min(50, row['total_markers'])}</b>
             ({row['pct_top50']}%)</div>
        <div>All markers: <b>{row['in_panel_all']}/{row['total_markers']}</b>
             ({row['pct_all']}%)</div>
      </div>
      <div style="font-size:12px;margin-bottom:4px">
        <b>Best panel markers:</b> {top_in_str}<br>
        <b>Top missing:</b> <span style="color:#999">{top_miss_str}</span>
      </div>
      <details><summary style="cursor:pointer;font-size:13px;color:#555">
        Show top 30 markers</summary>
      <table class="gene-tbl">
        <tr><th>Gene</th><th>Score</th><th>logFC</th><th>Panel</th></tr>
        {gene_rows}
      </table>
      </details>
    </div>
    """)

table_rows = ""
for _, row in summary.iterrows():
    risk_col = rc(row["risk"])
    top_in_str = ", ".join(row["top_in"]) if row["top_in"] else "NONE"
    top_miss_str = ", ".join(row["top_miss"]) if row["top_miss"] else "—"
    table_rows += f"""<tr>
      <td><span style="background:{risk_col};color:#fff;padding:2px 8px;border-radius:8px;
           font-size:12px">{row['risk']}</span></td>
      <td><b>{row['celltype']}</b></td>
      <td>{row['total_markers']}</td>
      <td>{pct_bar(row['pct_top20'])}</td>
      <td style="text-align:center"><b>{row['in_panel_top20']}</b></td>
      <td>{pct_bar(row['pct_top50'], 90)}</td>
      <td style="font-size:12px">{top_in_str}</td>
      <td style="font-size:12px;color:#999">{top_miss_str}</td>
    </tr>\n"""

html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<title>Xenium Panel Coverage — xenium_celltype (16 types)</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         margin:20px 40px; background:#fafafa; color:#222; }}
  h1 {{ border-bottom:3px solid #2c3e50; padding-bottom:8px; }}
  h2 {{ color:#2c3e50; margin-top:36px; }}
  .card {{ background:#fff; border:1px solid #ddd; border-radius:8px;
           padding:14px 18px; margin:10px 0; }}
  .stat-box {{ display:inline-block; background:#fff; border:1px solid #ddd;
               border-radius:8px; padding:14px 20px; margin:6px; text-align:center; }}
  .stat-box .num {{ font-size:36px; font-weight:700; }}
  .stat-box .lbl {{ font-size:13px; color:#666; }}
  table {{ border-collapse:collapse; font-size:13px; }}
  th, td {{ padding:5px 10px; border-bottom:1px solid #eee; text-align:left; }}
  th {{ background:#f5f5f5; position:sticky; top:0; }}
  .gene-tbl {{ margin-top:6px; }}
  details {{ margin-top:4px; }}
  .filter-btn {{ padding:4px 12px; border-radius:16px; border:1px solid #ccc;
                 background:#fff; cursor:pointer; font-size:13px; margin:2px; }}
  .filter-btn.active {{ background:#2c3e50; color:#fff; border-color:#2c3e50; }}
</style>
</head><body>

<h1>Xenium Panel Coverage — xenium_celltype</h1>
<p style="color:#555">
  Final assessment: <b>16 cell types</b> derived from level-1 with targeted splits
  (Epithelial → Ciliated/Secretory, DC → pDC/cDC, T/NK → T/NK).
  Markers computed via Wilcoxon one-vs-rest on a <b>{MAX_CELLS:,} cells/type</b>
  downsample against the <b>{len(panel_genes)}-gene Xenium panel</b>.
</p>

<h2>Executive Summary</h2>
<div style="display:flex;flex-wrap:wrap">
  <div class="stat-box"><div class="num" style="color:#2c3e50">{n_total}</div>
    <div class="lbl">Cell types</div></div>
  <div class="stat-box"><div class="num" style="color:#e74c3c">{n_high}</div>
    <div class="lbl">HIGH RISK</div></div>
  <div class="stat-box"><div class="num" style="color:#f39c12">{n_mod}</div>
    <div class="lbl">MODERATE</div></div>
  <div class="stat-box"><div class="num" style="color:#27ae60">{n_ok}</div>
    <div class="lbl">OK</div></div>
</div>
<p style="font-size:13px;color:#666">
  <b>Risk criteria (top-20 DEG markers in panel):</b>
  HIGH RISK = &lt;3 &nbsp;|&nbsp; MODERATE = 3–5 &nbsp;|&nbsp; OK = &gt;5
</p>

<h2>Coverage Bar Chart</h2>
{bar_items}

<h2>Summary Table</h2>
<div style="overflow-x:auto">
<table>
  <tr><th>Risk</th><th>Cell Type</th><th># Markers</th><th>Top-20 Coverage</th>
      <th>Hits</th><th>Top-50 Coverage</th><th>Best Panel Markers</th><th>Top Missing</th></tr>
  {table_rows}
</table>
</div>

<h2>Cell-Type Detail Cards</h2>
<div>
  <button class="filter-btn active" onclick="filterCards('all')">All</button>
  <button class="filter-btn" onclick="filterCards('HIGH RISK')" style="color:#e74c3c">HIGH RISK</button>
  <button class="filter-btn" onclick="filterCards('MODERATE')" style="color:#f39c12">MODERATE</button>
  <button class="filter-btn" onclick="filterCards('OK')" style="color:#27ae60">OK</button>
</div>
{"".join(detail_cards)}

<script>
function filterCards(risk) {{
  document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
  event.target.classList.add('active');
  document.querySelectorAll('.card').forEach(c => {{
    if (risk === 'all') {{ c.style.display = ''; return; }}
    const span = c.querySelector('h3 span');
    if (!span) {{ c.style.display = ''; return; }}
    c.style.display = span.textContent.trim() === risk ? '' : 'none';
  }});
}}
</script>

<p style="font-size:11px;color:#aaa;margin-top:40px">
  Generated {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")} ·
  Xenium panel: {len(panel_genes)} genes ·
  Downsample: {MAX_CELLS:,} cells/type · 16 xenium_celltype labels
</p>

</body></html>
"""

with open(OUT_HTML, "w") as f:
    f.write(html)

print(f"\n{'='*60}\nDONE\n{'='*60}")
print(f"  Downsampled h5ad:  {OUT_H5AD}")
print(f"  Marker CSVs:       {MARKER_DIR}/")
print(f"  Coverage report:   {OUT_HTML}")
print(summary["risk"].value_counts().to_string())
