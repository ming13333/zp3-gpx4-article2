# -*- coding: utf-8 -*-
"""
A2 v14 Figure 6 (replaces the v13 version; same layout, v14 numbers).

Panel a: observed vs cluster-permutation null-expected significant-gene counts
         for four frameworks -> inflated false-positive FLOOR under naive
         pooling (10.1%) versus a valid ~4.5-4.9% elsewhere.
Panel b: excess over null (observed - null mean) by framework. The fall from
         +2,929 (on-treatment only) to +2,083 (naive pooled) is directionally
         consistent but its 95% CI includes zero, so it is labelled exploratory.

Reads results/a2_v14_cluster_null_frozen.json (10,000 permutations).
Pure stdlib SVG.
"""
import os, json

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = os.path.join(ROOT, "article2", "results")
FIG = os.path.join(ROOT, "article2", "figures")
os.makedirs(FIG, exist_ok=True)

D = json.load(open(os.path.join(RES, "a2_v14_cluster_null_frozen.json"), encoding="utf-8"))
PERM = D["perm"]
F42 = D["perm_fixed42"]

ORDER = [
    ("A_pre", "A", ["pre-treatment", "only"], "#1f3b73"),
    ("On_only", "On", ["on-treatment", "only"], "#0f6b4a"),
    ("B_naive", "B", ["naive pooled", "(pre + on)"], "#c0392b"),
    ("C_patient", "C", ["patient-level", "aggregated"], "#5b4b8a"),
]

rows = []
for key, short, lab, col in ORDER:
    p = PERM[key]
    rows.append(dict(key=key, short=short, lab=lab, col=col,
                     n=p["n"], npt=p["n_patients"], nR=p["nR_patients"],
                     obs=p["observed_scipy"], null=p["null_mean"],
                     nullpct=p["null_pct"], permP=p["perm_P"], excess=p["excess"]))

F_ON = F42["P_on42"]["null_pct"]
F_POOL = F42["P_pool84"]["null_pct"]
F_PRE = F42["P_pre42"]["null_pct"]


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def txt(x, y, s, size=11, anchor="start", fill="#000", weight="normal"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" '
            f'fill="{fill}" font-weight="{weight}" '
            f'font-family="Arial, Helvetica, sans-serif">{esc(s)}</text>')


def tspan_text(x, y, spans, size=9.8, fill="#333"):
    """spans = [(text, bold), ...] rendered as one flowing line."""
    inner = "".join(
        f'<tspan font-weight="{"bold" if b else "normal"}" fill="{("#222" if b else fill)}">{esc(t)}</tspan>'
        for t, b in spans)
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" xml:space="preserve" '
            f'font-family="Arial, Helvetica, sans-serif">{inner}</text>')


def line(x1, y1, x2, y2, stroke="#000", width=1, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{width}"{d}/>')


def rect(x, y, w, h, fill, stroke="none", sw=0):
    return (f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


W, H = 940, 742
s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
     f'viewBox="0 0 {W} {H}" font-family="Arial, Helvetica, sans-serif">',
     rect(0, 0, W, H, "#ffffff")]

s.append(txt(40, 32, "Figure 6. Patient-cluster permutation separates an inflated false-positive floor from a genuine on-treatment signal",
             size=14.5, weight="bold"))
s.append(txt(40, 50, "GSE91061, all 22,187 gene rows of the matrix. Null distributions from 10,000 permutations of patient labels, preserving within-patient repeated-measure structure. Mann-Whitney P with tie and continuity correction.",
             size=9.2, fill="#555"))

# ================= PANEL a =================
s.append(txt(40, 84, "a", size=15, weight="bold"))
s.append(txt(62, 84, "Observed vs. null-expected significant genes", size=12, weight="bold", fill="#333"))

AX, AY, AW, AH = 96, 460, 360, 340
YMAX = 5000.0
ay = lambda v: AY - (v / YMAX) * AH

s.append(line(AX, AY, AX + AW, AY, width=1.3))
s.append(line(AX, AY, AX, AY - AH - 6, width=1.3))
for gv in range(0, 5001, 1000):
    y = ay(gv)
    if gv:
        s.append(line(AX, y, AX + AW, y, stroke="#e2e2e2", width=0.7, dash="2 3"))
    s.append(txt(AX - 8, y + 4, f"{gv:,}", size=9, anchor="end", fill="#444"))
s.append(txt(AX - 8, AY - AH - 18, "genes", size=9, anchor="end", fill="#444"))
s.append(txt(AX - 8, AY - AH - 6, "P<0.05", size=9, anchor="end", fill="#444"))

slot = AW / len(rows)
BW = 34
for i, r in enumerate(rows):
    cx = AX + slot * (i + 0.5)
    xo = cx - BW - 4
    s.append(rect(xo, ay(r["obs"]), BW, (r["obs"] / YMAX) * AH, r["col"]))
    s.append(txt(xo + BW / 2, ay(r["obs"]) - 6, f"{r['obs']:,}", size=9.5, anchor="middle", weight="bold"))
    xn = cx + 4
    bad = r["nullpct"] > 7
    s.append(rect(xn, ay(r["null"]), BW, (r["null"] / YMAX) * AH,
                  "#f0b3ab" if bad else "#c9ccd1", stroke=("#c0392b" if bad else "none"), sw=(1.6 if bad else 0)))
    s.append(txt(xn + BW / 2, ay(r["null"]) - 6, f"{r['null']:,.0f}", size=9.5, anchor="middle", fill="#444"))
    s.append(txt(cx, AY + 26, f"{r['nullpct']:.1f}%", size=11.5, anchor="middle",
                 weight="bold", fill=("#c0392b" if bad else "#1a7f5a")))
    for k, part in enumerate(r["lab"]):
        s.append(txt(cx, AY + 48 + k * 13, part, size=10, anchor="middle",
                     weight="bold" if k == 0 else "normal", fill="#000" if k == 0 else "#555"))
    s.append(txt(cx, AY + 48 + len(r["lab"]) * 13, f"n={r['npt']}", size=8.5, anchor="middle", fill="#777"))

s.append(txt(AX + AW / 2, AY + 98, "null share of the nominal 5% threshold (patient-cluster permutation)",
             size=9.2, anchor="middle", fill="#555"))
s.append(txt(AX + AW / 2, AY + 113,
             f"patient set held fixed (42 paired patients): {F_PRE:.2f}% pre-only, {F_ON:.2f}% on-only, {F_POOL:.2f}% pooled",
             size=8.6, anchor="middle", fill="#c0392b"))

lgx, lgy = AX + AW - 168, AY - AH + 16
s.append(rect(lgx, lgy - 9, 12, 11, "#1f3b73"))
s.append(txt(lgx + 17, lgy, "observed", size=9.5, fill="#333"))
s.append(rect(lgx + 84, lgy - 9, 12, 11, "#c9ccd1"))
s.append(txt(lgx + 101, lgy, "null mean", size=9.5, fill="#333"))

# ================= PANEL b =================
BX = 560
s.append(txt(BX, 84, "b", size=15, weight="bold"))
s.append(txt(BX + 22, 84, "Excess over null, by framework", size=12, weight="bold", fill="#333"))

BAX, BAW = BX + 56, 320
BAY0 = 420.0          # y of value 0
BSCALE = 360.0 / 4200.0
byy = lambda v: BAY0 - v * BSCALE

s.append(line(BAX, byy(-900), BAX, byy(3300) - 6, width=1.3))
for gv in (-500, 0, 1000, 2000, 3000):
    y = byy(gv)
    if gv:
        s.append(line(BAX, y, BAX + BAW, y, stroke="#e2e2e2", width=0.7, dash="2 3"))
    s.append(txt(BAX - 8, y + 4, f"{gv:+,}" if gv else "0", size=9, anchor="end", fill="#444"))
s.append(line(BAX, byy(0), BAX + BAW, byy(0), width=1.3))

slotb = BAW / len(rows)
BWb = 40
for i, r in enumerate(rows):
    cx = BAX + slotb * (i + 0.5)
    v = r["excess"]
    col = "#c0392b" if r["key"] == "B_naive" else ("#0f6b4a" if r["key"] == "On_only" else "#9aa0a6")
    top = byy(v) if v > 0 else byy(0)
    hgt = abs(byy(v) - byy(0))
    s.append(rect(cx - BWb / 2, top, BWb, hgt, col))
    s.append(txt(cx, (top - 7) if v > 0 else (byy(0) + hgt + 14), f"{v:+,}",
                 size=10, anchor="middle", weight="bold", fill=col))
    s.append(txt(cx, byy(-900) + 24, r["short"], size=10.5, anchor="middle", weight="bold"))

xc_on = BAX + slotb * 1.5
xc_b = BAX + slotb * 2.5
y_ann = byy(3180)
s.append(line(xc_on + BWb / 2, y_ann, xc_b - BWb / 2, y_ann, stroke="#c0392b", width=1.4, dash="4 3"))
s.append(f'<path d="M{xc_b - BWb/2:.1f} {y_ann:.1f} l-9 -5 l0 10 z" fill="#c0392b"/>')
s.append(txt((xc_on + xc_b) / 2, y_ann - 8, "difference -1,517", size=9.0, anchor="middle",
             weight="bold", fill="#c0392b"))
s.append(txt((xc_on + xc_b) / 2, y_ann - 20, "95% CI -3,423 to +149", size=8.2, anchor="middle",
             fill="#c0392b"))
s.append(txt(BAX + BAW / 2, byy(-900) + 44, "excess = observed - null mean", size=9.2, anchor="middle", fill="#555"))

# ================= TAKE-HOME =================
ty = 596
lines = [
    [("Pseudo-replication is real and quantified: ", True),
     ("under naive pooling the null share of the nominal 5% threshold doubles to 10.1%, so the same test", False)],
    [("", False),
     ("over-reports \"significant\" genes roughly twofold when the null is true.", False)],
    [("The mechanism is the analysis unit, not patient composition: ", True),
     ("holding the 42 paired patients fixed, pooling alone", False)],
    [("", False),
     (f"moves the null share from {F_ON:.2f}% to {F_POOL:.2f}% (pre-only {F_PRE:.2f}%).", False)],
    [("The rise in gene counts is dominated by genuine on-treatment biology: ", True),
     ("on-treatment samples alone give 4,015 genes", False)],
    [("", False), ("at P<0.05 (permutation P=0.0016; excess +2,929) on a valid 4.9% null floor.", False)],
    [("The accompanying fall in excess (to +2,083) is directional but not resolved: ", True),
     ("its 95% CI includes zero, so it is", False)],
    [("", False), ("reported as exploratory rather than as an established second effect.", False)],
]
for j, spans in enumerate(lines):
    s.append(tspan_text(40, ty + j * 14, spans))

s.append("</svg>")

OUT = os.path.join(FIG, "A2_Fig6_v14.svg")
open(OUT, "w", encoding="utf-8").write("\n".join(s))
print("WROTE", OUT)
print("panel a:", [(r["short"], r["obs"], round(r["null"]), r["nullpct"]) for r in rows])
print("panel b:", [(r["short"], r["excess"]) for r in rows])
print("fixed42:", F_PRE, F_ON, F_POOL)
