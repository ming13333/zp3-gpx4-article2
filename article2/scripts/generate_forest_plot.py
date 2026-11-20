# -*- coding: utf-8 -*-
"""
Supplementary Figure S1: cross-cohort forest plot of the ZP3-response rank-biserial r.
Reads the frozen result tables rather than hard-coded values.
r_rb = 2U/(n1*n2) - 1; positive means higher ZP3 in responders.
No pooled effect is estimated, because the cohorts differ in tumour type, platform and treatment.
The GSE91061 sensitivity row reports the 105-sample pooled set with both P values, and quotes no
interval for it, since the reference distribution for that row is the patient-cluster permutation.
"""
import os, csv, json, math

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article2", "figures")
os.makedirs(OUT, exist_ok=True)
FROZEN = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article2", "results", "a2_effects_frozen.csv")
SENS = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article2", "results", "a2_effects_sensitivity.csv")

cohorts = []
with open(FROZEN, encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        cohorts.append({
            "lab": (r["Cohort"] + "\n" + r["Tissue"].split(",")[0] + ", " + r["Treatment"]),
            "r": float(r["r_rb"]), "lo": float(r["CI_low"]), "hi": float(r["CI_high"]),
            "nr": int(r["n_resp"]), "nn": int(r["n_non"]),
        })
# fixed cohort order
order = {"GSE78220": 0, "IMvigor210": 1, "GSE91061": 2}
cohorts.sort(key=lambda c: order.get(c["lab"].split("\n")[0], 9))

# GSE91061 sensitivity analysis, including on-treatment samples
sens_note = None
with open(SENS, encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r["Cohort"] == "GSE91061" and "All" in r["Scope"]:
            sens_r, sens_p = float(r["r_rb"]), float(r["P"])

# the pooled 105-sample row is described by both P values; no interval is quoted for it
sens_note = f"GSE91061 pooled, 105 samples: r={sens_r:+.2f}, nominal P={sens_p:.3f}"
CLUST = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                      "article2", "results", "a2_zp3_cluster_perm_v15_20261119.json")
cluster_p = None
with open(CLUST, encoding="utf-8") as _fh:
    _b = json.load(_fh)["B"]
cluster_p = _b["p_addone"] if _b["n_samples"] == 105 else None
if sens_note is not None and cluster_p is not None:
    sens_note += f", cluster P={cluster_p:.3f}"

W, H = 780, 660
x0, x1 = 300, 620          # forest-plot x extent
y0, y1 = 100, 480          # plotting y extent
rmin, rmax = -0.9, 0.9     # r axis, symmetric about zero

def lin(v, vmin, vmax, px0, px1):
    return px0 + (v - vmin) / (vmax - vmin) * (px1 - px0)

def txt(x, y, s, size=11, anchor="start", fill="#000", weight="normal", transform=""):
    s = str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    t = f' transform="{transform}"' if transform else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" fill="{fill}" font-weight="{weight}" font-family="Arial, Helvetica, sans-serif"{t}>{s}</text>'

def line(x1, y1, x2, y2, stroke="#000", width=1, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{width}"{d}/>'

def rect(x, y, w, h, fill, stroke="#000", width=1):
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>'

s = f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Arial, Helvetica, sans-serif">'
s += rect(0, 0, W, H, "#ffffff")
s += txt(40, 34, "Supplementary Figure S1. Cross-cohort effect sizes for the ZP3-response association", size=15, weight="bold")
s += txt(40, 52, "Rank-biserial r (95% CI); positive = higher ZP3 in responders; dashed line = null.", size=10, fill="#555")

# zero line
z = lin(0, rmin, rmax, x0, x1)
s += line(x0, y0, x0, y1, "#ccc", 1)
s += line(z, y0, z, y1, "#333", 1, "4,3")

# axis ticks
for tick, lab in [(-0.75, "-0.75"), (-0.5, "-0.5"), (-0.25, "-0.25"), (0, "0"), (0.25, "0.25"), (0.5, "0.5"), (0.75, "0.75")]:
    xt = lin(tick, rmin, rmax, x0, x1)
    s += line(xt, y0, xt, y0 + 6, "#666", 1)
    s += txt(xt, y0 + 20, lab, 8, "middle", fill="#555")

# one row per cohort
n = len(cohorts)
row_h = (y1 - y0) / n
for i, c in enumerate(cohorts):
    cy = y0 + row_h * (i + 0.5)
    for j, ln in enumerate(c["lab"].split("\n")):
        s += txt(40, cy - 6 + j * 14, ln, 10)
    s += txt(248, cy + 2, f"n={c['nr']+c['nn']}\n({c['nr']}/{c['nn']})", 9, "middle", fill="#555")
    xlo, xhi, xr = lin(c["lo"], rmin, rmax, x0, x1), lin(c["hi"], rmin, rmax, x0, x1), lin(c["r"], rmin, rmax, x0, x1)
    s += line(xlo, cy, xhi, cy, "#666", 2)
    s += line(xlo, cy - 8, xlo, cy + 8, "#666", 2) + line(xhi, cy - 8, xhi, cy + 8, "#666", 2)
    s += rect(xr - 6, cy - 6, 12, 12, "#2980b9", "#000", 1)  # all three primary analyses are non-significant: blue
    s += txt(x1 + 14, cy + 2, f"{c['r']:.3f} ({c['lo']:.2f}, {c['hi']:.2f})", 10, "start")

# legend and notes
s += txt(40, y1 + 40, "Blue dots: not statistically significant (all three primary analyses).", 10, fill="#555")
s += txt(40, y1 + 58, "Effect estimates varied across cohorts; no pooled effect was estimated because the cohorts differ in", 9, fill="#777")
s += txt(40, y1 + 72, "tumor type, treatment, platform and analysis sample (pre-treatment only unless noted).", 9, fill="#777")
if sens_note:
    s += txt(40, y1 + 90, sens_note + " (sensitivity; not a primary analysis).", 9, fill="#c0392b")
s += txt(40, y1 + 108, "Sources: a2_effects_frozen.csv for the three cohort rows; a2_zp3_cluster_perm_v15_20261119.json", size=9, fill="#777")
s += txt(40, y1 + 122, "for the pooled row. Cohort intervals are normal approximations of rank-biserial r.", size=9, fill="#777")
s += "</svg>"

with open(os.path.join(OUT, "A2_SuppS1.svg"), "w", encoding="utf-8") as f:
    f.write(s)
print(f"written A2_SuppS1.svg ({len(s)} bytes); cohorts={[c['lab'].split(chr(10))[0] for c in cohorts]}")
