# -*- coding: utf-8 -*-
"""
A2 v13 Figure 5 generator.

Message (v13 framing): the GSE91061 ZP3 association is confined to analyses
that INCLUDE on-treatment samples. Adding the on-treatment-only tier and a
second patient-level aggregation rule shows that (i) the on-treatment-only
association is itself robust under patient-cluster permutation (P=0.013),
and (ii) the patient-level result is sensitive to the aggregation rule
(pre-else-on P=0.175 vs mean-aggregate P=0.054), so the honest statement is
"absent from the pre-specified framework", not "absent from every framework".

Pure stdlib SVG. Writes figures/A2_Fig5_v13.svg
"""
import os, math

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG = os.path.join(ROOT, "article2", "figures")
os.makedirs(FIG, exist_ok=True)


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def txt(x, y, s, size=12, anchor="start", fill="#000", weight="normal", transform=""):
    t = f' transform="{transform}"' if transform else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" '
            f'fill="{fill}" font-weight="{weight}"{t} '
            f'font-family="Arial, Helvetica, sans-serif">{esc(s)}</text>')


def line(x1, y1, x2, y2, stroke="#000", width=1, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{width}"{d}/>')


def rect(x, y, w, h, fill, stroke="none", sw=0):
    return (f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


W, H = 1080, 668
s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
     f'font-family="Arial, Helvetica, sans-serif">', rect(0, 0, W, H, "#ffffff")]

s.append(txt(40, 32, "Figure 5. The GSE91061 ZP3 association is confined to analyses that include on-treatment samples",
             size=14.5, weight="bold"))
s.append(txt(40, 50, "Left: -log10(P) for the ZP3-response association across eleven analysis frameworks; right: subtype-stratified P values (pre-treatment, n=49).",
             size=9.5, fill="#555"))

# ---------------- LEFT PANEL ----------------
LX0, LY0, LY1 = 118, 470, 120          # axis x, y at zero, y at vmax
VMAX = 2.8
ly = lambda v: LY0 - (v / VMAX) * (LY0 - LY1)

s.append(line(LX0, LY0, LX0, LY1 - 8, width=1.3))
for gv in (0.0, 1.0, 2.0):
    y = ly(gv)
    if gv:
        s.append(line(LX0, y, 700, y, stroke="#e4e4e4", width=0.7, dash="2 3"))
    s.append(txt(LX0 - 8, y + 4, f"{gv:.0f}", size=9, anchor="end", fill="#444"))
s.append(txt(74, (LY0 + LY1) / 2, "-log10 (P)", size=10, anchor="middle",
             transform=f"rotate(-90 74 {(LY0+LY1)/2})"))

# groups: (label lines, -log10 P, colour, perm P or None)
G_INCL = "#c0392b"
G_PRE = "#8d9298"
G_FIRTH = "#b9bec4"
G_AGG = "#5b4b8a"

groups = [
    [("on-treatment", "only", "n=56 (55 pt)"), ("pooled", "pre+on", "n=105 (62 pt)")],
    [("pre-treatment", "only", "n=49"), ("logistic", "univariable", "OR 0.72"),
     ("logistic", "+ IFN-\u03b3", ""), ("logistic", "+ suppression", "")],
    [("Firth", "univariable", "OR 0.82"), ("Firth", "+ IFN-\u03b3", ""), ("Firth", "+ suppression", "")],
    [("patient-level", "pre-else-on", "n=62"), ("patient-level", "mean rule", "n=62")],
]
vals = [0.0110, 0.0081, 0.3522, 0.4545, 0.5154, 0.5772, 0.6060, 0.8300, 0.6660, 0.1753, 0.0539]
cols = [G_INCL, G_INCL, G_PRE, G_PRE, G_PRE, G_PRE, G_FIRTH, G_FIRTH, G_FIRTH, G_AGG, G_AGG]
perm = {0: 0.0131, 1: 0.0279, 2: 0.3573}

BW, SLOT, GAP = 34, 40, 18
x = 150.0
centers, flat = [], []
for gi, g in enumerate(groups):
    for bi, lab in enumerate(g):
        cx = x + BW / 2
        centers.append(cx)
        flat.append((lab, vals[len(flat)], cols[len(flat)], len(flat)))
        x += SLOT
    x += GAP

for lab, v, col, idx in flat:
    cx = centers[idx]
    lg = math.log10(1.0 / v)
    h = (lg / VMAX) * (LY0 - LY1)
    s.append(rect(cx - BW / 2, LY0 - h, BW, h, col))
    s.append(txt(cx, LY0 - h - 6, f"{lg:.2f}", size=8.5, anchor="middle", weight="bold"))
    for j, ln in enumerate(lab):
        s.append(txt(cx, LY0 + 15 + j * 10, ln, size=7.6, anchor="middle",
                     fill="#000" if j == 0 else "#555"))
    if idx in perm:
        s.append(txt(cx, LY0 - h - 17, f"\u2020 perm P={perm[idx]:.3f}", size=7.2, anchor="middle", fill=col))

# P=0.05 threshold
ythr = ly(math.log10(1 / 0.05))
s.append(line(LX0, ythr, 700, ythr, "#333", 1.1, "4,3"))
s.append(txt(704, ythr + 3, "P=0.05", size=8.5, fill="#333"))

# group brackets
gx = 150.0
glabels = ["includes\non-treatment", "pre-treatment\n(nominal)", "pre-treatment\nFirth-penalised",
           "patient-level\naggregation"]
for gi, g in enumerate(groups):
    gstart = gx
    gend = gx + SLOT * len(g) - (SLOT - BW)
    gx += SLOT * len(g) + GAP
    yb = LY0 + 62
    s.append(line(gstart, yb, gend, yb, "#999", 0.9))
    s.append(txt((gstart + gend) / 2, yb + 13, glabels[gi].split("\n")[0], size=7.8,
                 anchor="middle", weight="bold", fill="#444"))
    s.append(txt((gstart + gend) / 2, yb + 23, glabels[gi].split("\n")[1], size=7.8,
                 anchor="middle", fill="#666"))

# ---------------- RIGHT PANEL ----------------
RX0, RY0, RY1 = 800, 470, 120
RVMAX = 1.2
s.append(line(RX0, RY0, RX0, RY1 - 8, width=1.3))
for gv in (0.0, 0.4, 0.8, 1.2):
    y = RY0 - (gv / RVMAX) * (RY0 - RY1)
    s.append(line(RX0, y, RX0 + 5, y, "#333", 1))
    s.append(txt(RX0 - 6, y + 4, f"{gv:.1f}", size=8.5, anchor="end", fill="#444"))
s.append(txt(766, (RY0 + RY1) / 2, "-log10 (P)", size=10, anchor="middle",
             transform=f"rotate(-90 766 {(RY0+RY1)/2})"))
sub = [("Immune", 0.8001, "n=20"), ("Keratin", 0.7692, "n=14")]
for i, (g, pv, nn) in enumerate(sub):
    cx = 860 + i * 95
    lg = math.log10(1 / pv)
    h = (lg / RVMAX) * (RY0 - RY1)
    s.append(rect(cx - 27, RY0 - h, 54, h, "#8d9298"))
    s.append(txt(cx, RY0 - h - 7, f"P={pv:.2f}", size=8.5, anchor="middle", weight="bold"))
    s.append(txt(cx, RY0 + 15, g, size=9.5, anchor="middle"))
    s.append(txt(cx, RY0 + 27, nn, size=8, anchor="middle", fill="#666"))
s.append(txt(907, 108, "Subtype-stratified (pre-treatment)", size=9.5, anchor="middle", weight="bold"))
s.append(txt(907, 124, "No subtype association", size=9, anchor="middle", weight="bold", fill="#c0392b"))
s.append(txt(907, 138, "MITF-high/low: not estimable", size=8, anchor="middle", fill="#888"))

# ---------------- FOOTER ----------------
fy = 588
notes = [
    "Only frameworks that include on-treatment samples reach nominal significance. \u2020Patient-cluster permutation (500 permutations, preserving within-patient",
    "repeated measures) leaves the on-treatment-only association significant (P=0.013) but weakens the pooled association (P=0.028), and returns P=0.357 for the",
    "pre-treatment framework - i.e. the pooled signal is a genuine on-treatment pharmacodynamic association misreported as baseline prediction, not a false positive.",
    "The patient-level result depends on the aggregation rule (pre-else-on P=0.175 vs mean-aggregate P=0.054); no aggregation rule reaches P<0.05.",
]
for j, a in enumerate(notes):
    s.append(txt(40, fy + j * 13, a, size=8.8, fill="#333"))

s.append("</svg>")
OUT = os.path.join(FIG, "A2_Fig5_v13.svg")
open(OUT, "w", encoding="utf-8").write("\n".join(s))
print("WROTE", OUT)
print("bars:", len(flat))
