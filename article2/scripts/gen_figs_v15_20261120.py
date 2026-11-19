# -*- coding: utf-8 -*-
"""A2 v15 figures 1-3 (Briefings in Bioinformatics version).

Fig 1  three targets and where a pooled analysis sits (schematic)
Fig 2  simulation audit: (a) S1 calibration, (b) S2 after a common shift,
       (c) S3 split by gene subset, (d) S4 detection among the injected genes
Fig 3  real-data audit: (a) null floor by framework, (b) observed vs null-expected
       counts, (c) the four ZP3 estimates with bootstrap intervals

All numbers are read from frozen JSON; nothing is typed by hand. Pure stdlib SVG,
Arial 10px minimum, single-column-safe palette.
"""
import glob, json, os, sys

sys.stdout.reconfigure(encoding="utf-8")
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

C_POOL = "#b03a2e"      # pooled sample-level
C_CLUST = "#1f4e79"     # cluster permutation of the same statistic
C_PAT = "#4b7f52"       # patient-level test
C_DIRECT = "#7d5ba6"    # direct test of change
C_GREY = "#6b6b6b"
FONT = "Arial, Helvetica, sans-serif"


# ------------------------------------------------------------------ helpers
def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def txt(x, y, s, size=10.5, anchor="start", fill="#000", weight="normal", style=""):
    st = f' font-style="{style}"' if style else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" '
            f'fill="{fill}" font-weight="{weight}" font-family="{FONT}"{st}>{esc(s)}</text>')


def ln(x1, y1, x2, y2, stroke="#000", w=0.8, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{w}"{d}/>')


def rect(x, y, w, h, fill="none", stroke="none", sw=0.8, rx=0, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')


def circ(x, y, r, fill, stroke="none", sw=0.7, op=1.0):
    return (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{sw}" opacity="{op}"/>')


def yaxis(x, ytop, ybot, lo, hi, ticks, label, fmt="{:.0f}", label_x=None):
    out = [ln(x, ytop, x, ybot, "#000", 1.0)]
    for t in ticks:
        yy = ybot - (t - lo) / (hi - lo) * (ybot - ytop)
        out.append(ln(x - 3, yy, x, yy, "#000", 1.0))
        out.append(txt(x - 5, yy + 3.4, fmt.format(t), 9.8, "end", "#333"))
    out.append(txt(label_x if label_x is not None else x - 34, (ytop + ybot) / 2, label,
                   10.2, "middle", "#000").replace(
        "<text ", '<text transform="rotate(-90 %.1f %.1f)" ' % (
            label_x if label_x is not None else x - 34, (ytop + ybot) / 2)))
    return "".join(out)


def vy(x, ybot, ytop, lo, hi, v):
    return ybot - (v - lo) / (hi - lo) * (ybot - ytop)


def jitter(i, n, width):
    if n == 1:
        return 0.0
    return (i - (n - 1) / 2.0) / max(1, n - 1) * width


# ------------------------------------------------------------------ data
def load_cells():
    out = []
    for f in glob.glob(os.path.join(RES, "a2_p2_sim_core_20261118_full*.json")):
        if "smoke" in f:
            continue
        d = json.load(open(f, encoding="utf-8"))
        c = d["cfg"]
        s = {(r["scenario"], r["target"], r["method"]): r for r in d["summary"] if "rate" in r}
        out.append(dict(np=c["n_pat"], rp=c["resp_prop"], icc=c["icc"], N=c["n_outer"], s=s))
    return out


CELLS = load_cells()
GRID = [c for c in CELLS if c["N"] == 500]
MAIN = [c for c in CELLS if c["N"] == 1000][0]


def rate(cell, key):
    return 100.0 * cell["s"][key]["rate"]


def mc_halfwidth(cell, key):
    r = cell["s"][key]
    return 100.0 * (r["ci95"][1] - r["ci95"][0]) / 2.0


# ------------------------------------------------------------------ Figure 1
def fig1():
    W, H = 680, 336
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">',
         rect(0, 0, W, H, "#ffffff")]
    x0, y0, bw, bh, gx, gy = 122, 78, 118, 62, 12, 44
    xR = x0 + 2 * bw + gx
    yB = y0 + bh + gy
    gm = (y0 + bh + yB) / 2.0

    cells = [("Pre R", 0, 0), ("Pre N", 1, 0), ("On R", 0, 1), ("On N", 1, 1)]
    for lab, c, r in cells:
        x, y = x0 + c * (bw + gx), y0 + r * (bh + gy)
        o.append(rect(x, y, bw, bh, "#e8eef7" if lab.endswith("R") else "#f2f2f2", "#333", 1.0))
        o.append(txt(x + bw / 2, y + bh / 2 + 4, lab, 11.5, "middle", "#000", "bold"))
    o.append(txt(x0 - 12, y0 + bh / 2 + 4, "pre", 10.5, "end", "#000", "bold"))
    o.append(txt(x0 - 12, yB + bh / 2 + 4, "on", 10.5, "end", "#000", "bold"))

    # baseline: spanning the two pre-treatment cells, above the grid
    yb = y0 - 16
    o.append(ln(x0 + 4, yb, xR - 4, yb, C_DIRECT, 1.4))
    o.append(ln(x0 + 4, yb, x0 + 4, yb + 8, C_DIRECT, 1.4))
    o.append(ln(xR - 4, yb, xR - 4, yb + 8, C_DIRECT, 1.4))
    o.append(txt((x0 + xR) / 2, yb - 7, "baseline", 10.5, "middle", C_DIRECT, "bold"))

    # on-treatment: spanning the two on-treatment cells, below the grid
    yo = yB + bh + 14
    o.append(ln(x0 + 4, yo, xR - 4, yo, C_DIRECT, 1.4))
    o.append(ln(x0 + 4, yo - 8, x0 + 4, yo, C_DIRECT, 1.4))
    o.append(ln(xR - 4, yo - 8, xR - 4, yo, C_DIRECT, 1.4))
    o.append(txt((x0 + xR) / 2, yo + 15, "on-treatment", 10.5, "middle", C_DIRECT, "bold"))

    # change: paired, inside each column across the gap between the rows
    for c in (0, 1):
        xc = x0 + c * (bw + gx) + bw / 2
        o.append(ln(xc, y0 + bh + 6, xc, yB - 6, C_PAT, 1.4))
        for yy, s in ((y0 + bh + 6, 1), (yB - 6, -1)):
            o.append(f'<path d="M{xc - 3.6:.1f} {yy + 4.4 * s:.1f} L{xc:.1f} {yy:.1f} '
                     f'L{xc + 3.6:.1f} {yy + 4.4 * s:.1f}" fill="none" stroke="{C_PAT}" '
                     f'stroke-width="1.4"/>')
    o.append(txt(x0 - 12, gm + 4, "change", 10.5, "end", C_PAT, "bold"))

    # pooled: one box around all four cells
    o.append(rect(x0 - 8, y0 - 4, xR - x0 + 16, 2 * bh + gy + 8, "none", C_POOL, 1.5, 0, "6 4"))
    o.append(txt(xR + 22, gm - 10, "pooled", 10.8, "start", C_POOL, "bold"))
    o.append(txt(xR + 22, gm + 4, "all four cells;", 10, "start", C_POOL))
    o.append(txt(xR + 22, gm + 17, "weights = pair counts;", 10, "start", C_POOL))
    o.append(txt(xR + 22, gm + 30, "null = none of the three", 10, "start", C_POOL))

    y = H - 40
    o.append(txt(24, y, "Baseline and on-treatment comparisons use one sample per patient; the change question uses the",
                 10.2))
    o.append(txt(24, y + 14, "paired difference within a patient. The pooled analysis uses every sample, so its estimate",
                 10.2))
    o.append(txt(24, y + 28, "is a weighted average of the four cell values and its null is the null of none of the three.",
                 10.2))
    o.append("</svg>")
    p = os.path.join(FIG, "A2v15_Fig1.svg")
    open(p, "w", encoding="utf-8").write("\n".join(o))
    return p


# ------------------------------------------------------------------ Figure 2
def fig2():
    W, H = 680, 560
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">',
         rect(0, 0, W, H, "#ffffff")]

    def ptitle(x, y, letter, title):
        o.append(txt(x, y, letter, 12, "start", "#000", "bold"))
        o.append(txt(x + 14, y, title, 10.4, "start", "#000", "bold"))

    # ---------------- (a) S1 calibration
    ax, ay, aw, ah = 84, 70, 240, 144
    ptitle(24, ay - 10, "a", "S1: complete null")
    lo, hi = 0, 16
    o.append(yaxis(ax, ay, ay + ah, lo, hi, [0, 4, 8, 12, 16], "genes rejected, %"))
    o.append(ln(ax, vy(ax, ay + ah, ay, lo, hi, 5), ax + aw, vy(ax, ay + ah, ay, lo, hi, 5),
                C_GREY, 0.9, "4 3"))
    cats = [("pooled\nsample level", ("S1", "T1", "mw_pooled"), C_POOL),
            ("cluster\npermutation", ("S1", "T1", "pooled_cluster"), C_CLUST),
            ("patient\nlevel", ("S1", "T1", "mw_baseline"), C_PAT)]
    for i, (lab, key, col) in enumerate(cats):
        cx = ax + 32 + i * (aw - 44) / 3.0
        for j, cell in enumerate(GRID):
            o.append(circ(cx + jitter(j, len(GRID), 24), vy(ax, ay + ah, ay, lo, hi, rate(cell, key)),
                          2.2, col, op=0.5))
        mv = rate(MAIN, key)
        o.append(ln(cx - 14, vy(ax, ay + ah, ay, lo, hi, mv), cx + 14,
                    vy(ax, ay + ah, ay, lo, hi, mv), col, 1.8))
        for k, part in enumerate(lab.split("\n")):
            o.append(txt(cx, ay + ah + 15 + 11 * k, part, 10, "middle", "#000"))
    o.append(txt(ax + aw / 2 - 10, ay + ah + 44, "18 cells (dots), main configuration (bar)",
                 9.8, "middle", "#444"))

    # ---------------- (b) S2 after a common shift
    bx, by, bw, bh = 412, 70, 236, 144
    ptitle(352, by - 10, "b", "S2: after a common shift")
    mech = json.load(open(os.path.join(RES, "a2_p2_s2_mechanism_20261118.json"), encoding="utf-8"))
    ratios = {}
    for row in mech["rows"]:
        ratios.setdefault(row["icc"], []).append((row["tau"], row["s2"]))
    lo2, hi2 = 0.4, 2.0
    o.append(yaxis(bx, by, by + bh, lo2, hi2, [0.4, 0.8, 1.2, 1.6, 2.0],
                   "variance ratio of the pooled statistic", fmt="{:.1f}", label_x=bx - 48))
    o.append(ln(bx, vy(bx, by + bh, by, lo2, hi2, 1.0), bx + bw,
                vy(bx, by + bh, by, lo2, hi2, 1.0), C_GREY, 0.9, "4 3"))
    cols = {0.2: "#8fb3d9", 0.5: "#4f81bd", 0.8: "#1f4e79"}
    for icc, pts in sorted(ratios.items()):
        pts = sorted(pts)
        d = " ".join(f'{"M" if i == 0 else "L"}{bx + 14 + (t / 2.0) * (bw - 30):.1f} '
                     f'{vy(bx, by + bh, by, lo2, hi2, v):.1f}' for i, (t, v) in enumerate(pts))
        o.append(f'<path d="{d}" fill="none" stroke="{cols[icc]}" stroke-width="1.6"/>')
        for t, v in pts:
            o.append(circ(bx + 14 + (t / 2.0) * (bw - 30), vy(bx, by + bh, by, lo2, hi2, v), 3.2,
                          cols[icc]))
        o.append(txt(bx + 14 + 8, vy(bx, by + bh, by, lo2, hi2, pts[0][1]) - 5,
                     f"ICC {icc:g}", 10, "start", cols[icc], "bold"))
    for t in (0.0, 0.5, 1.0, 2.0):
        o.append(txt(bx + 14 + (t / 2.0) * (bw - 30), by + bh + 15, f"{t:g}", 10, "middle", "#000"))
    o.append(txt(bx + bw / 2, by + bh + 29, "common shift, log2", 10, "middle", "#000"))
    o.append(txt(bx + bw / 2 - 6, by + bh + 44,
                 "post hoc experiment; grid S2 sits at 1.0", 9.8, "middle", "#444"))

    # ---------------- (c) S3 by gene subset
    cx0, cy, chh = 84, 320, 144
    ptitle(24, cy - 10, "c", "S3: association without differential change")
    part = json.load(open(os.path.join(RES, "a2_p2_partition_audit_20261120.json"), encoding="utf-8"))
    P = part["per_cell"]
    alo, ahi = 0, 20
    o.append(yaxis(cx0, cy, cy + chh, alo, ahi, [0, 5, 10, 15, 20], "genes rejected, %"))
    o.append(ln(cx0, vy(cx0, cy + chh, cy, alo, ahi, 5), cx0 + 152,
                vy(cx0, cy + chh, cy, alo, ahi, 5), C_GREY, 0.9, "4 3"))
    left = [("all 500\ngenes", "all_genes", C_POOL),
            ("450 null\ngenes", "null", C_CLUST)]
    for i, (lab, key, col) in enumerate(left):
        xx = cx0 + 42 + i * 78
        for j, c in enumerate(P):
            o.append(circ(xx + jitter(j, len(P), 22),
                          vy(cx0, cy + chh, cy, alo, ahi, 100 * c["S3"][key]), 2.2, col, op=0.55))
        mv = 100 * part["summary"]["S3_main"][key]
        o.append(ln(xx - 13, vy(cx0, cy + chh, cy, alo, ahi, mv), xx + 13,
                    vy(cx0, cy + chh, cy, alo, ahi, mv), col, 1.8))
        for k, t in enumerate(lab.split("\n")):
            o.append(txt(xx, cy + chh + 15 + 11 * k, t, 10, "middle", "#000"))
    # right sub-scale, association genes only
    rx = 300
    o.append(ln(248, cy - 6, 248, cy + chh, "#c8c8c8", 0.9))
    o.append(yaxis(rx, cy, cy + chh, 0, 100, [0, 25, 50, 75, 100], " ", fmt="{:.0f}"))
    xx = 340
    for j, c in enumerate(P):
        o.append(circ(xx + jitter(j, len(P), 26),
                      vy(rx, cy + chh, cy, 0, 100, 100 * c["S3"]["signal"]), 2.2, C_POOL, op=0.55))
    mv = 100 * part["summary"]["S3_main"]["signal"]
    o.append(ln(xx - 15, vy(rx, cy + chh, cy, 0, 100, mv), xx + 15,
                vy(rx, cy + chh, cy, 0, 100, mv), C_POOL, 1.8))
    o.append(txt(xx, cy + chh + 15, "50 signal", 10, "middle", "#000"))
    o.append(txt(xx, cy + chh + 26, "genes", 10, "middle", "#000"))
    o.append(txt(cx0, cy + chh + 44,
                 "left axis 0-20%, right axis 0-100%; every dot is a pooled test",
                 9.8, "start", "#444"))

    # ---------------- (d) S4 detection among the injected genes
    dx, dy, dw, dh = 412, 320, 236, 144
    ptitle(352, dy - 10, "d", "S4: injected change, 50 genes")
    o.append(yaxis(dx, dy, dy + dh, 0, 100, [0, 25, 50, 75, 100], "genes detected, %"))
    cats4 = [("on-\ntreatment", "mw_on_signal", C_PAT),
             ("direct\nchange", "delta_mwu_signal", C_DIRECT),
             ("pooled", "pooled_nocorr_signal", "#d98880"),
             ("pooled\n+cluster", "signal", C_CLUST)]
    for i, (lab, key, col) in enumerate(cats4):
        xx = dx + 30 + i * 54
        for j, c in enumerate(P):
            o.append(circ(xx + jitter(j, len(P), 22),
                          vy(dx, dy + dh, dy, 0, 100, 100 * c["S4"][key]), 2.2, col, op=0.55))
        mv = 100 * part["summary"]["S4_main"][key]
        o.append(ln(xx - 12, vy(dx, dy + dh, dy, 0, 100, mv), xx + 12,
                    vy(dx, dy + dh, dy, 0, 100, mv), col, 1.8))
        for k, t in enumerate(lab.split("\n")):
            o.append(txt(xx, dy + dh + 15 + 11 * k, t, 10, "middle", "#000"))
    o.append(txt(dx + 46, dy + dh + 44, "denominator: the 50 injected genes", 9.8,
                 "middle", "#444"))

    # footer
    o.append(txt(24, H - 30, "Small dots: 18 grid cells at 500 replicates. Horizontal bars: the main configuration at "
                             "1,000 replicates.", 9.8, "start", "#444"))
    o.append(txt(24, H - 15, "95% Monte Carlo intervals on the cell rates are about +/-0.1 points and are not drawn.",
                 9.8, "start", "#444"))
    o.append("</svg>")
    p = os.path.join(FIG, "A2v15_Fig2.svg")
    open(p, "w", encoding="utf-8").write("\n".join(o))
    return p


# ------------------------------------------------------------------ Figure 3
def fig3():
    W, H = 680, 500
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">',
         rect(0, 0, W, H, "#ffffff")]
    fr = json.load(open(os.path.join(RES, "a2_v14_cluster_null_frozen.json"), encoding="utf-8"))
    perm = fr["perm"]
    f42 = fr["perm_fixed42"]

    def ptitle(x, y, letter, title):
        o.append(txt(x, y, letter, 12, "start", "#000", "bold"))
        o.append(txt(x + 14, y, title, 10.4, "start", "#000", "bold"))

    # ---------------- (a) null floor by framework
    ax, ay, aw, ah = 84, 66, 236, 144
    ptitle(24, ay - 10, "a", "null floor by framework")
    lo, hi = 0, 13
    o.append(yaxis(ax, ay, ay + ah, lo, hi, [0, 4, 8, 12], "genes at P<0.05, %"))
    o.append(ln(ax, vy(ax, ay + ah, ay, lo, hi, 5), ax + aw, vy(ax, ay + ah, ay, lo, hi, 5),
                C_GREY, 0.9, "4 3"))
    bars = [("pre", perm["A_pre"]["null_pct"], "#9dc3a8"),
            ("on", perm["On_only"]["null_pct"], "#6f9e7c"),
            ("pooled", perm["B_naive"]["null_pct"], C_POOL),
            ("patient", perm["C_patient"]["null_pct"], C_PAT)]
    rings = [f42["P_pre42"]["null_pct"], f42["P_on42"]["null_pct"], f42["P_pool84"]["null_pct"]]
    for i, (lab, v, col) in enumerate(bars):
        x = ax + 30 + i * 54
        o.append(rect(x - 12, vy(ax, ay + ah, ay, lo, hi, v), 24,
                      ay + ah - vy(ax, ay + ah, ay, lo, hi, v), col))
        o.append(txt(x, vy(ax, ay + ah, ay, lo, hi, v) - 4, f"{v:.1f}", 9.8, "middle", "#000"))
        o.append(txt(x, ay + ah + 15, lab, 10, "middle", "#000"))
        if i < len(rings):
            o.append(circ(x + 14, vy(ax, ay + ah, ay, lo, hi, rings[i]), 4.2, "#ffffff", C_CLUST, 1.2))
    o.append(txt(ax + aw / 2 - 22, ay + ah + 31,
                 "rings (offset right): same framework on the 42 paired patients", 9.8, "middle", C_CLUST))

    # ---------------- (b) observed vs null-expected
    bx, by, bww, bhh = 420, 66, 228, 144
    ptitle(352, by - 10, "b", "observed vs null-expected hits")
    lo2, hi2 = 0, 5000
    o.append(yaxis(bx, by, by + bhh, lo2, hi2, [0, 1000, 2000, 3000, 4000, 5000],
                   "genes at P<0.05 (count)", fmt="{:.0f}", label_x=bx - 44))
    pairs = [("pre", perm["A_pre"]["observed_scipy"], perm["A_pre"]["null_mean"]),
             ("on", perm["On_only"]["observed_scipy"], perm["On_only"]["null_mean"]),
             ("pooled", perm["B_naive"]["observed_scipy"], perm["B_naive"]["null_mean"]),
             ("patient", perm["C_patient"]["observed_scipy"], perm["C_patient"]["null_mean"])]
    for i, (lab, obs, nul) in enumerate(pairs):
        x = bx + 28 + i * 52
        o.append(rect(x - 5, vy(bx, by + bhh, by, lo2, hi2, obs), 10,
                      by + bhh - vy(bx, by + bhh, by, lo2, hi2, obs), "#3d3d3d"))
        o.append(ln(x - 13, vy(bx, by + bhh, by, lo2, hi2, nul), x + 13,
                    vy(bx, by + bhh, by, lo2, hi2, nul), C_POOL, 1.8))
        o.append(txt(x, by + bhh + 15, lab, 10, "middle", "#000"))
    o.append(txt(bx + bww / 2 - 30, by + bhh + 31, "bars observed, red lines null mean",
                 9.8, "middle", "#444"))

    # ---------------- (c) the four ZP3 estimates on the same 42 patients
    z = json.load(open(os.path.join(RES, "a2_zp3_pooled_interval_diag_20261120.json"), encoding="utf-8"))
    fx = json.load(open(os.path.join(RES, "a2_zp3_fixed42_20261119.json"), encoding="utf-8"))
    cy2, ch2 = 266, 130
    ptitle(24, cy2 - 12, "c", "ZP3, the same 42 patients")
    xl, xr = 150, 620
    lo3, hi3 = -1.0, 0.4

    def X(v):
        return xl + (v - lo3) / (hi3 - lo3) * (xr - xl)

    o.append(ln(X(0), cy2 - 4, X(0), cy2 + ch2, C_GREY, 0.9, "4 3"))
    o.append(ln(xl, cy2 + ch2, xr, cy2 + ch2, "#000", 1.0))
    for t in (-1.0, -0.75, -0.5, -0.25, 0.0, 0.25):
        o.append(ln(X(t), cy2 + ch2, X(t), cy2 + ch2 + 3, "#000", 1.0))
        o.append(txt(X(t), cy2 + ch2 + 15, f"{t:g}", 9.8, "middle", "#333"))
    o.append(txt((xl + xr) / 2, cy2 + ch2 + 30, "rank-biserial r", 10, "middle", "#000"))
    rows = [("baseline", fx["pre42"]["r"], fx["pre42"]["boot_CI"], C_CLUST),
            ("on-treatment", fx["on42"]["r"], fx["on42"]["boot_CI"], C_CLUST),
            ("change", fx["delta42"]["r"], fx["delta42"]["boot_CI"], C_DIRECT),
            ("pooled (84 samples)", z["r"], z["bootstrap"]["CI"], C_POOL)]
    for i, (lab, r, ci, col) in enumerate(rows):
        yy = cy2 + 10 + i * 28
        o.append(txt(xl - 12, yy + 3.6, lab, 10, "end", "#000"))
        o.append(ln(X(ci[0]), yy, X(ci[1]), yy, col, 1.6))
        o.append(ln(X(ci[0]), yy - 4, X(ci[0]), yy + 4, col, 1.6))
        o.append(ln(X(ci[1]), yy - 4, X(ci[1]), yy + 4, col, 1.6))
        o.append(circ(X(r), yy, 3.8, col))
        o.append(txt(X(r), yy - 8, f"{r:.2f}", 9.8, "middle", col))
    o.append(txt(24, H - 40, "Pooled row: 84 samples from 42 patients, r=-0.36, 95% interval -0.64 to -0.07, "
                             "nominal P=0.020,", 9.8, "start", C_POOL))
    o.append(txt(24, H - 25, "cluster-permutation P=0.069. Clustered rows: 42 patients. Intervals are 95% percentile "
                             "patient-cluster bootstrap.", 9.8, "start", C_POOL))
    o.append(txt(24, H - 9, "Dotted lines mark 5% (a) and zero (c).", 9.8, "start", "#444"))
    o.append("</svg>")
    p = os.path.join(FIG, "A2v15_Fig3.svg")
    open(p, "w", encoding="utf-8").write("\n".join(o))
    return p


if __name__ == "__main__":
    for f in (fig1(), fig2(), fig3()):
        print("wrote", f, os.path.getsize(f), "bytes")
