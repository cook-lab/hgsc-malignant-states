#!/usr/bin/env python3
"""
Fig 6 ROI picker — pick an ROI on a whole Xenium tissue (or the WT01 Visium HD bins) and preview candidate ROI rows
in house style. The page's rows "g", "h" and "k" are working labels of an earlier layout; the published 6i tiles
are fixed in 13_fig6_draw_panels.py, which also imports this file as a library (palettes, ROI loader, polygons).

  python roi_server.py            -> http://127.0.0.1:8793  (or $PORT)

Tissues: WT01, WT02 (SingleR labels, raw counts from fig6_viewer/degdata, cell polygons from the raw
Xenium cell_boundaries.parquet). Every cell in the ROI is drawn (no subsampling). Saved ROIs are appended to
$OUTPUT_ROOT/fig6_singler/fig6_rois.json, which 13_fig6_draw_panels.py reads: the 6i ROI by name, and the last saved
Visium window for panels 6a-d.

Row spec syntax (editable in the page): comma-separated tiles; "celltype" = cell-type map; "GENE@types" colours only
those cell types (epi, secb, mac, tnk, stroma, all); "GENE1*GENE2@types" = co-expression (geometric mean of the two
2-98% scaled values). Default types are chosen per gene (see DEFAULT_TYPES).
"""
import io, os, json, gzip, time, threading, datetime
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
import numpy as np, pandas as pd, pyarrow.parquet as pq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm, to_rgba
from scipy.io import mmread; from scipy import sparse; from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get("OUTPUT_ROOT", "output"); DEG = f"{R}/fig6_viewer/degdata"; SG = f"{R}/singler_recompute"; O = f"{R}/fig6_singler"
CACHE = f"{R}/fig6_roi_picker/cache"; os.makedirs(CACHE, exist_ok=True)
RAW = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "202609_final_hgsc_spatial/data/HGSC xenium ")
SAMPLES = {"WT01": f"{RAW}/WT05_WT01/cell_boundaries.parquet", "WT02": f"{RAW}/WT02/cell_boundaries.parquet"}
ROIFILE = f"{O}/fig6_rois.json"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7})
PAL = {"SecA epithelium": "#E6A141", "Intermediate epithelium": "#C08E48", "SecB epithelium": "#9A7D55",
       "Ciliated epithelium": "#E07850", "Mesothelial": "#D4A574", "Fibroblast": "#DDD5CA", "Smooth muscle": "#D14E6C",
       "Pericyte": "#B87A7A", "Endothelial": "#7D4E4E", "Macrophage": "#8FBC8F", "Mast cell": "#8B9B6B",
       "Neutrophil": "#6B8E23", "Conventional dendritic cell": "#2E8B57", "Plasmacytoid dendritic cell": "#A0A0A0",
       "NK cell": "#56AFC4", "T cell": "#87CEFA", "Plasma cell": "#8A5DAF", "B cell": "#5665B6"}
GRAD = LinearSegmentedColormap.from_list("grad", list(zip([0, .28, .55, .80, 1], ["#F6EFE5", "#ECDDD0", "#D89E97", "#A03A4A", "#3A111A"])))
DIV = LinearSegmentedColormap.from_list("div", ["#3F4E9C", "#9FA8D6", "#F4F1EC", "#D89E97", "#A03A4A", "#5A1A26"])
GREY, GREY_EDGE = "#E8E8E8", "#C4C4C4"
EPI = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]
TYPES = {"epi": EPI, "secb": ["SecB epithelium"], "mac": ["Macrophage"], "tnk": ["T cell", "NK cell"],
         "stroma": ["Fibroblast", "Pericyte", "Endothelial", "Smooth muscle"], "all": None}
DEFAULT_TYPES = {"VEGFA": "mac", "CD55": "epi", "IL32": "epi", "LGALS9": "epi", "HAVCR2": "tnk", "LAG3": "tnk", "SLC16A3": "epi",
                 "ICAM1": "epi", "CTSL": "all", "MMP7": "epi", "CD44": "mac", "ITGB6": "epi", "ITGAV": "epi", "BIRC3": "epi",
                 "NFKB2": "epi", "TGM2": "all", "LDHA": "all", "SLC2A1": "all", "NDRG1": "all", "PDCD1": "tnk"}
GLYCO = ["LDHA", "ALDOA", "ENO1", "PGK1", "PDK1", "SLC2A1", "SLC16A3", "NDRG1", "VEGFA"]; RESID = ["C1QC", "MRC1", "A2M", "C7"]
DATA = {}; LOCK = threading.Lock(); RLOCK = threading.Lock()   # matplotlib is not thread-safe

def load(s):
    with LOCK:
        if s in DATA: return DATA[s]
        t0 = time.time()
        d = pd.read_csv(f"{SG}/{s}_singler_cells.csv.gz", usecols=["barcode", "x", "y", "cell_label_singler"])
        d = d[d.cell_label_singler.notna() & (d.cell_label_singler != "")]
        genes = np.array([l.strip() for l in open(f"{DEG}/{s}_genes.txt")])
        cb = pd.read_csv(f"{DEG}/{s}_cells.csv.gz", usecols=["barcode"]).barcode.values
        with gzip.open(f"{DEG}/{s}_counts.mtx.gz", "rb") as f: M = sparse.csr_matrix(mmread(io.BytesIO(f.read()))).T.tocsr()
        tot = np.asarray(M.sum(1)).ravel(); tot[tot == 0] = 1
        L = sparse.diags(1e4 / tot) @ M; L.data = np.log1p(L.data); L = L.tocsc()
        pos = pd.Series(np.arange(len(cb)), index=cb); d = d[d.barcode.isin(pos.index)].reset_index(drop=True)
        row = pos.loc[d.barcode].values
        gi = {g: i for i, g in enumerate(genes)}
        mac = (d.cell_label_singler == "Macrophage").values
        def z(gl):
            M_ = np.column_stack([L[row, gi[g]].toarray().ravel() for g in gl if g in gi])[mac]
            return ((M_ - M_.mean(0)) / (M_.std(0) + 1e-9)).mean(1)
        g_, r_ = z(GLYCO), z(RESID); prog = np.full(len(d), np.nan); prog[mac] = (g_ - g_.mean()) / g_.std() - (r_ - r_.mean()) / r_.std()
        d["mac_program"] = prog
        ni = pd.read_csv(f"{O}/mac_niche_WT.csv.gz"); ni = ni[ni["sample"] == s].set_index("barcode").group
        d["niche"] = d.barcode.map(ni)
        # polygons: cache per sample as npz (cell order = d)
        cp = f"{CACHE}/{s}_poly.npz"
        if not os.path.exists(cp):
            B = pq.read_table(SAMPLES[s], columns=["cell_id", "vertex_x", "vertex_y"]).to_pandas()
            B = B[B.cell_id.isin(set(d.barcode))]
            B["k"] = pd.Categorical(B.cell_id, categories=d.barcode.values).codes
            B = B.sort_values("k", kind="stable"); k = B.k.values
            starts = np.searchsorted(k, np.arange(len(d))); ends = np.searchsorted(k, np.arange(len(d)), side="right")
            np.savez(cp, vx=B.vertex_x.values.astype(np.float32), vy=B.vertex_y.values.astype(np.float32), s=starts, e=ends)
        _z = np.load(cp); P = {k: _z[k] for k in ("vx", "vy", "s", "e")}   # materialise (npz access is lazy)
        cen_x = np.array([P["vx"][a:b].mean() if b > a else np.nan for a, b in zip(P["s"][:2000], P["e"][:2000])])
        off = np.nanmedian(cen_x - d.x.values[:2000])
        DATA[s] = dict(d=d, L=L, row=row, gi=gi, P=P, tree=cKDTree(d[["x", "y"]].values),
                       ext=(float(d.x.min()), float(d.x.max()), float(d.y.min()), float(d.y.max())))
        print(f"[{s}] loaded {len(d):,} cells in {time.time()-t0:.0f}s; polygon-centroid offset {off:+.2f} um", flush=True)
        return DATA[s]

def overview(s, mode):
    fp = f"{CACHE}/{s}_overview_{mode}.png"
    D = load(s)
    if not os.path.exists(fp):
        d = D["d"]; x0, x1, y0, y1 = D["ext"]; W = 9.0; H = W * (y1 - y0) / (x1 - x0)
        fig = plt.figure(figsize=(W, H), dpi=160); ax = fig.add_axes([0, 0, 1, 1])
        lab = d.cell_label_singler.values
        if mode == "celltype":
            ax.scatter(d.x, d.y, s=1.1, c=[PAL.get(l, GREY) for l in lab], lw=0, rasterized=True)
        else:   # niche: SecB brown, SecA orange (faint), niche mac dark green, other mac sage, T/NK blue
            base = np.array([to_rgba(GREY)] * len(d)); base[np.isin(lab, ["SecA epithelium", "Intermediate epithelium"])] = to_rgba("#EFC98A")
            base[lab == "SecB epithelium"] = to_rgba(PAL["SecB epithelium"]); base[lab == "Macrophage"] = to_rgba("#8FBC8F")
            base[(d.niche == "niche").values] = to_rgba("#1F5E36"); base[np.isin(lab, ["T cell", "NK cell"])] = to_rgba("#1E88E5")
            order = np.argsort(np.isin(lab, ["Macrophage", "T cell", "NK cell"]))
            ax.scatter(d.x.values[order], d.y.values[order], s=1.4, c=base[order], lw=0, rasterized=True)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.axis("off"); fig.savefig(fp, dpi=160); plt.close(fig)
    return fp

def roi_idx(D, x, y, size):
    d = D["d"]; h = size / 2
    return np.where((d.x.values >= x - h) & (d.x.values <= x + h) & (d.y.values >= y - h) & (d.y.values <= y + h))[0]

def polys(D, idx):
    P = D["P"]; return [np.column_stack([P["vx"][P["s"][i]:P["e"][i]], P["vy"][P["s"][i]:P["e"][i]]]) for i in idx]

def gene_vals(D, idx, g):
    if g not in D["gi"]: return None
    return D["L"][D["row"][idx], D["gi"][g]].toarray().ravel()

def scaled(v, show):
    lo, hi = np.nanpercentile(v[show], [2, 98]) if show.any() else (0, 1)
    return np.clip((v - lo) / (hi - lo if hi > lo else 1), 0, 1)

def render_row(s, x, y, size, spec, title):
  with RLOCK:
    return _render_row(s, x, y, size, spec, title)

def _render_row(s, x, y, size, spec, title):
    D = load(s); idx = roi_idx(D, x, y, size); d = D["d"].iloc[idx]; pl = polys(D, idx); lab = d.cell_label_singler.values
    tiles = [t.strip() for t in spec.split(",") if t.strip()]
    n = len(tiles); fig = plt.figure(figsize=(1.35 * n, 1.75), dpi=220)
    for j, t in enumerate(tiles):
        ax = fig.add_axes([j / n + 0.01, 0.16, 1 / n - 0.02, 0.68])
        if t.lower() == "celltype":
            faces = [to_rgba(PAL.get(l, GREY)) for l in lab]; ttl = "Cell type"; edge = "#8A8A8A"; order = np.arange(len(pl)); sub = ""
        elif t.lower() == "program":
            m = lab == "Macrophage"; nrm = TwoSlopeNorm(0, -2, 2)
            faces = [DIV(nrm(np.clip(v, -2, 2))) if mm else (to_rgba(PAL["SecB epithelium"], 0.35) if l == "SecB epithelium" else to_rgba(GREY))
                     for v, mm, l in zip(d.mac_program.values, m, lab)]
            ttl = "Mac program"; edge = GREY_EDGE; order = np.argsort(m, kind="stable"); sub = "resident ← → hypoxic"
        else:
            gpart, _, tpart = t.partition("@"); gl = gpart.split("*")
            tk = tpart.strip() or DEFAULT_TYPES.get(gl[0].strip().upper(), "all"); types = TYPES.get(tk, None)
            show = np.ones(len(lab), bool) if types is None else np.isin(lab, types)
            vs = [gene_vals(D, idx, g.strip().upper()) for g in gl]
            if any(v is None for v in vs):
                ax.text(0.5, 0.5, f"{gpart}\nnot on Xenium panel", ha="center", va="center", fontsize=6, color="#999", transform=ax.transAxes)
                ax.axis("off"); ax.set_title(gpart, fontsize=7, fontweight="bold", fontstyle="italic"); continue
            val = scaled(vs[0], show) if len(vs) == 1 else np.sqrt(np.prod([scaled(v, show) for v in vs], axis=0))
            faces = [GRAD(v) if sh else to_rgba(GREY) for v, sh in zip(val, show)]
            order = np.lexsort((np.where(show, val, -1),)); ttl = " × ".join(g.strip().upper() for g in gl); edge = GREY_EDGE
            sub = {"epi": "epithelium", "secb": "SecB", "mac": "macrophages", "tnk": "T / NK cells", "stroma": "stroma", "all": "all cells"}.get(tk, tk)
        ax.add_collection(PolyCollection([pl[i] for i in order], facecolors=[faces[i] for i in order], edgecolors=edge, linewidths=0.06))
        if t.lower() not in ("celltype", "program") and "@tnk" in t.replace(" ", "") or (t.lower() not in ("celltype", "program") and DEFAULT_TYPES.get(t.split("@")[0].strip().upper()) == "tnk" and "@" not in t):
            sel = [i for i in order if show[i]]   # sparse lymphocytes: redraw on top with a dark outline
            ax.add_collection(PolyCollection([pl[i] for i in sel], facecolors=[faces[i] for i in sel], edgecolors="#111", linewidths=0.5))
        ax.set_xlim(x - size / 2, x + size / 2); ax.set_ylim(y - size / 2, y + size / 2); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(ttl, fontsize=7, fontweight="bold", fontstyle="normal" if t.lower() in ("celltype", "program") else "italic", pad=2)
        ax.text(0.5, -0.04, sub, transform=ax.transAxes, ha="center", va="top", fontsize=5.5, color="#555")
        if j == 0:
            L_ = 100 if size >= 300 else 50; xs = x + size / 2 - size * 0.05 - L_; ys = y - size / 2 + size * 0.04
            ax.plot([xs, xs + L_], [ys, ys], color="k", lw=1.2, solid_capstyle="butt"); ax.text(xs + L_ / 2, ys + size * 0.025, f"{L_} µm", ha="center", fontsize=5)
    fig.suptitle(title, fontsize=7, x=0.01, ha="left", y=0.99, color="#444")
    b = io.BytesIO(); fig.savefig(b, format="png", facecolor="white"); plt.close(fig); return b.getvalue()

def roi_stats(s, x, y, size):
    D = load(s); idx = roi_idx(D, x, y, size); d = D["d"].iloc[idx]; lab = d.cell_label_singler
    vc = lab.value_counts()
    tnk = d[lab.isin(["T cell", "NK cell"])]; nm = d[d.niche == "niche"]
    near = 0
    if len(tnk) and len(nm): near = int((cKDTree(nm[["x", "y"]].values).query(tnk[["x", "y"]].values, distance_upper_bound=30)[0] < 30).sum())
    return dict(cells=int(len(d)), SecB=int(vc.get("SecB epithelium", 0)), SecA=int(vc.get("SecA epithelium", 0)),
                Int=int(vc.get("Intermediate epithelium", 0)), mac=int(vc.get("Macrophage", 0)), niche_mac=int(len(nm)),
                T=int(vc.get("T cell", 0)), NK=int(vc.get("NK cell", 0)), TNK_within30_of_niche_mac=near,
                mean_mac_program=float(np.nanmean(d.mac_program)) if d.mac_program.notna().any() else None)

def suggest(s, size, kind):
    D = load(s); d = D["d"]; nm = d[d.niche == "niche"][["x", "y"]].values
    if not len(nm): return []
    tr = cKDTree(nm); h = size / 2
    if kind == "tnk":   # windows centred on niche macs, scored by T/NK cells within 30 um of a niche mac inside the window
        tn = d[d.cell_label_singler.isin(["T cell", "NK cell"])][["x", "y"]].values
        close = tn[tr.query(tn, distance_upper_bound=30)[0] < 30] if len(tn) else np.zeros((0, 2))
        if not len(close): return []
        tc = cKDTree(close); score = np.array([len(v) for v in tc.query_ball_point(nm, r=h, p=np.inf)])
    else:
        score = np.array([len(v) for v in tr.query_ball_point(nm, r=h, p=np.inf)])
    out, taken = [], []
    for i in np.argsort(-score):
        c = nm[i]
        if any(np.max(np.abs(c - t)) < size for t in taken): continue
        taken.append(c); out.append(dict(x=float(c[0]), y=float(c[1]), score=int(score[i])))
        if len(out) == 6: break
    return out


# ============================================================ Visium HD mode (WT01, 16 um bins; the window of panels 6a-d)
VBIN_UM = 16.0; VASPECT = 0.6            # preview window height = 0.6 x width; 13_fig6_draw_panels.py draws the saved window as a square
VCOL = {"other": GREY, "mac (other)": "#CFE3CF", "epithelium": "#EFE6D8", "SecA-hi": PAL["SecA epithelium"], "SecB-hi": PAL["SecB epithelium"],
        "mac (distal)": "#8FBC8F", "mac (SecB niche)": "#2F6B45"}
VORDER = ["other", "epithelium", "mac (other)", "SecA-hi", "SecB-hi", "mac (distal)", "mac (SecB niche)"]
VIS = {}
def vload():
    with LOCK:
        if not VIS:
            B = pd.read_csv(f"{O}/visium_n2_bins_WT01.csv.gz"); VIS["B"] = B
            VIS["ext"] = (float(B.col.min()), float(B.col.max()), float(B.row.min()), float(B.row.max()))
        return VIS
def voverview():
    fp = f"{CACHE}/VISIUM_overview.png"; V = vload()
    if not os.path.exists(fp):
        B = V["B"]; x0, x1, y0, y1 = V["ext"]; Wd = 9.0; Hd = Wd * (y1 - y0) / (x1 - x0)
        fig = plt.figure(figsize=(Wd, Hd), dpi=160); ax = fig.add_axes([0, 0, 1, 1])
        for c in VORDER:
            q = B[B.cls == c]; ax.scatter(q.col, q.row, s=1.6, marker="s", c=VCOL[c], lw=0, rasterized=True)
        ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis("off"); fig.savefig(fp, dpi=160); plt.close(fig)
    return fp
def vwin(x, y, w_um):
    wx = w_um / VBIN_UM; wy = wx * VASPECT; B = vload()["B"]
    return B[(B.col.between(x - wx / 2, x + wx / 2)) & (B.row.between(y - wy / 2, y + wy / 2))].copy(), wx, wy
def vstats(x, y, w_um):
    Bz, wx, wy = vwin(x, y, w_um); vc = Bz.cls.value_counts()
    return dict(cells=int(len(Bz)), SecB=int(vc.get("SecB-hi", 0)), SecA=int(vc.get("SecA-hi", 0)), Int=int(vc.get("epithelium", 0)),
                mac=int(vc.get("mac (SecB niche)", 0) + vc.get("mac (distal)", 0) + vc.get("mac (other)", 0)), niche_mac=int(vc.get("mac (SecB niche)", 0)),
                T=0, NK=0, TNK_within30_of_niche_mac=0, mean_mac_program=None, unit="16 um bins",
                mm=f"{w_um/1000:.2f} x {w_um*VASPECT/1000:.2f} mm")
def vrender(x, y, w_um, title):
  with RLOCK:
    Bz, wx, wy = vwin(x, y, w_um); TWd, THd = 2.5, 2.5 * VASPECT
    fig = plt.figure(figsize=(2 * TWd + 0.6, THd + 0.75), dpi=220)
    s_bin = (TWd * 72 / wx) ** 2 * 0.95
    def tile(i, ttl):
        ax = fig.add_axes([(0.15 + i * (TWd + 0.3)) / (2 * TWd + 0.6), 0.2, TWd / (2 * TWd + 0.6), THd / (THd + 0.75)])
        ax.set_xlim(x - wx / 2, x + wx / 2); ax.set_ylim(y + wy / 2, y - wy / 2); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(ttl, fontsize=7, fontweight="bold", pad=2); return ax
    a = tile(0, "b · SecA / SecB and macrophage bins")
    for c in VORDER:
        q = Bz[Bz.cls == c]; a.scatter(q.col, q.row, s=s_bin, marker="s", c=VCOL[c], lw=0)
    Lb = 250 / VBIN_UM; a.plot([x + wx / 2 - Lb - 3, x + wx / 2 - 3], [y + wy / 2 - 3] * 2, color="k", lw=1.2); a.text(x + wx / 2 - Lb / 2 - 3, y + wy / 2 - 5, "250 µm", fontsize=5, ha="center")
    c_ = tile(1, "c · CA9 (hypoxia), SecB-high outlined")
    if len(Bz):
        nb = cKDTree(Bz[["col", "row"]].values).query_ball_point(Bz[["col", "row"]].values, r=1.0, p=np.inf)
        v = np.array([Bz["CA9"].values[i].mean() for i in nb]); nn = Normalize(0, max(np.quantile(v, 0.98), 1e-6)); o_ = np.argsort(v)
        c_.scatter(Bz.col.values[o_], Bz.row.values[o_], s=s_bin, marker="s", lw=0, c=[GRAD(nn(t)) for t in v[o_]])
        sb = Bz[Bz.cls == "SecB-hi"]; c_.scatter(sb.col, sb.row, s=s_bin, marker="s", facecolors="none", edgecolors=PAL["SecB epithelium"], linewidths=0.45)
    fig.suptitle(title, fontsize=6.5, x=0.01, ha="left", y=0.99, color="#444")
    b = io.BytesIO(); fig.savefig(b, format="png", facecolor="white"); plt.close(fig); return b.getvalue()
def vsuggest(w_um):
    B = vload()["B"]; wx = w_um / VBIN_UM; wy = wx * VASPECT; P = B[B.cls == "mac (SecB niche)"][["col", "row"]].values
    S = B[B.cls == "SecB-hi"][["col", "row"]].values
    sc_ = np.array([((np.abs(P[:, 0] - q[0]) <= wx / 2) & (np.abs(P[:, 1] - q[1]) <= wy / 2)).sum() + 0.5 * ((np.abs(S[:, 0] - q[0]) <= wx / 2) & (np.abs(S[:, 1] - q[1]) <= wy / 2)).sum() for q in P])
    out, taken = [], []
    for i in np.argsort(-sc_):
        q = P[i]
        if any(abs(q[0] - t[0]) < wx and abs(q[1] - t[1]) < wy for t in taken): continue
        taken.append(q); out.append(dict(x=float(q[0]), y=float(q[1]), score=int(sc_[i])))
        if len(out) == 6: break
    return out

class H(BaseHTTPRequestHandler):
    def log_message(self, *a): pass
    def _send(self, code, body, ctype="application/json"):
        self.send_response(code); self.send_header("Content-Type", ctype); self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store"); self.end_headers(); self.wfile.write(body)
    def do_GET(self):
        try:
            if self.path in ("/", "/index.html"): return self._send(200, open(f"{HERE}/index.html", "rb").read(), "text/html; charset=utf-8")
            if self.path.startswith("/overview/VISIUM"): return self._send(200, open(voverview(), "rb").read(), "image/png")
            if self.path.startswith("/meta/VISIUM"): return self._send(200, json.dumps(dict(ext=vload()["ext"], yflip=True, unit_um=VBIN_UM, aspect=VASPECT)).encode())
            if self.path.startswith("/overview/"):
                s, mode = self.path.split("/")[2], self.path.split("/")[3].split("?")[0]
                return self._send(200, open(overview(s, mode), "rb").read(), "image/png")
            if self.path.startswith("/meta/"):
                D = load(self.path.split("/")[2]); return self._send(200, json.dumps(dict(ext=D["ext"])).encode())
            if self.path == "/rois":
                return self._send(200, (open(ROIFILE, "rb").read() if os.path.exists(ROIFILE) else b"[]"))
            self._send(404, b"{}")
        except Exception as e:
            self._send(500, json.dumps(dict(error=str(e))).encode())
    def do_POST(self):
        try:
            q = json.loads(self.rfile.read(int(self.headers["Content-Length"])) or b"{}")
            if q.get("sample") == "VISIUM" and self.path == "/render":
                return self._send(200, vrender(float(q["x"]), float(q["y"]), float(q["size"]), q.get("title", "")), "image/png")
            if q.get("sample") == "VISIUM" and self.path == "/stats": return self._send(200, json.dumps(vstats(float(q["x"]), float(q["y"]), float(q["size"]))).encode())
            if q.get("sample") == "VISIUM" and self.path == "/suggest": return self._send(200, json.dumps(vsuggest(float(q["size"])) if q.get("kind") != "tnk" else []).encode())
            if self.path == "/render":
                png = render_row(q["sample"], float(q["x"]), float(q["y"]), float(q["size"]), q["spec"], q.get("title", ""))
                return self._send(200, png, "image/png")
            if self.path == "/stats": return self._send(200, json.dumps(roi_stats(q["sample"], float(q["x"]), float(q["y"]), float(q["size"]))).encode())
            if self.path == "/suggest": return self._send(200, json.dumps(suggest(q["sample"], float(q["size"]), q.get("kind", "niche"))).encode())
            if self.path == "/save":
                rois = json.load(open(ROIFILE)) if os.path.exists(ROIFILE) else []
                q["saved_at"] = datetime.datetime.now().isoformat(timespec="seconds"); rois.append(q)
                json.dump(rois, open(ROIFILE, "w"), indent=1); return self._send(200, json.dumps(dict(ok=True, n=len(rois))).encode())
            self._send(404, b"{}")
        except Exception as e:
            import traceback; traceback.print_exc(); self._send(500, json.dumps(dict(error=str(e))).encode())

if __name__ == "__main__":
    for s in SAMPLES: load(s); overview(s, "celltype"); overview(s, "niche")
    vload(); voverview()
    port = int(os.environ.get("PORT", "8793"))
    print(f"ROI picker ready: http://127.0.0.1:{port}", flush=True)
    ThreadingHTTPServer(("127.0.0.1", port), H).serve_forever()
