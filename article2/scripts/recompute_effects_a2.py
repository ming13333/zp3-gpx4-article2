# -*- coding: utf-8 -*-
"""
A2 unified effect-size recomputation (science-audit fix, 2026-08-17)
=====================================================================
Fixes applied:
  1. rank-biserial r direction unified:  r_rb = 2*U/(n1*n2) - 1
     (U from mannwhitneyu(responders, non_responders); positive => higher ZP3 in responders)
  2. GSE91061 gene ID corrected: Entrez 8277 -> 7784 (true ZP3; 8277 is SP5).
  3. GSE91061 primary analysis restricted to pre-treatment (baseline) samples,
     consistent with GSE78220 (baseline) and IMvigor210 (single pre-therapy biopsy).
     All-sample (incl. on-treatment) results reported transparently as sensitivity.
  4. CI by normal approximation of the MWU statistic:
     SE = 2*sqrt((n+1)/(12*n1*n2)); CI = r +/- 1.96*SE.
  5. Power table recomputed: MDE(80%) ~= 2.80*SE (two-sided rank-biserial).

Outputs:
  - output/a2_effects_frozen.csv   (frozen effect table, primary analysis)
  - output/a2_effects_sensitivity.csv (GSE91061 visit-stratified + all-sample)
  - output/power_analysis_results.csv (updated)
"""
import os
import math
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "article2", "results", "a2_effects_frozen.csv")
OUT_S = os.path.join(ROOT, "article2", "results", "a2_effects_sensitivity.csv")
POWER = os.path.join(ROOT, "article2", "results", "power_analysis_results.csv")


def rank_biserial(x_resp, x_non):
    """Unified estimator. Positive r => responders higher. Returns (U, P, r, SE, CI_lo, CI_hi)."""
    x_resp = np.asarray(x_resp, float)
    x_non = np.asarray(x_non, float)
    n1, n2 = len(x_resp), len(x_non)
    u, p = stats.mannwhitneyu(x_resp, x_non, alternative="two-sided")
    r = 2.0 * u / (n1 * n2) - 1.0
    se = 2.0 * math.sqrt((n1 + n2 + 1.0) / (12.0 * n1 * n2))
    lo, hi = r - 1.96 * se, r + 1.96 * se
    return u, p, r, se, lo, hi


def med_iqr(s):
    s = np.asarray(s, float)
    return float(np.median(s)), float(np.percentile(s, 25)), float(np.percentile(s, 75))


def fmt_ci(lo, hi):
    return f"{lo:.2f}, {hi:.2f}"


rows = []

# ============================================================
# 1. GSE78220 (melanoma, anti-PD-1) — baseline samples only
# ============================================================
h2 = os.path.join(ROOT, "article2", "results", "gse78220_zp3_response.csv")
df20 = pd.read_csv(h2)
df20 = df20[df20["timepoint"] == "baseline"].copy()
df20["is_resp"] = df20["response"].isin(["CR", "PR"])
r20 = df20[df20["is_resp"]]["zp3_expr"]
n20 = df20[~df20["is_resp"]]["zp3_expr"]
u20, p20, rb20, se20, lo20, hi20 = rank_biserial(r20, n20)
m_r20, q1_r20, q3_r20 = med_iqr(r20)
m_n20, q1_n20, q3_n20 = med_iqr(n20)
rows.append({
    "Cohort": "GSE78220", "Tissue": "Melanoma", "Treatment": "Anti-PD-1",
    "Scope": "primary", "n_resp": len(r20), "n_non": len(n20),
    "med_resp": round(m_r20, 2), "med_non": round(m_n20, 2),
    "U": round(u20, 2), "P": p20, "r_rb": rb20, "SE": se20,
    "CI_low": lo20, "CI_high": hi20, "CI_str": fmt_ci(lo20, hi20),
})
print(f"[GSE78220] n={len(r20)+len(n20)} (R {len(r20)} / NR {len(n20)}); "
      f"med {m_r20:.2f} vs {m_n20:.2f}; U={u20:.0f} P={p20:.4f} r_rb={rb20:+.3f} "
      f"CI ({lo20:.2f}, {hi20:.2f})")

# ============================================================
# 2. IMvigor210 (urothelial, anti-PD-L1) — single pre-therapy biopsy
# ============================================================
im = os.path.join(ROOT, "output", "immunotherapy_validation")
ex = pd.read_csv(os.path.join(im, "exmat_censored_IMvigor210.csv"), index_col=0)
fd = pd.read_csv(os.path.join(im, "fData_IMvigor210.csv"))
pdt = pd.read_csv(os.path.join(im, "pData_IMvigor210.csv"))
ex.columns = [str(c).strip('"') for c in ex.columns]
fd.columns = [str(c).strip('"') for c in fd.columns]
pdt.columns = [str(c).strip('"') for c in pdt.columns]

# --- explicit row-match validation -----------------------------------
assert len(ex) == len(fd), f"exmat rows {len(ex)} != fData rows {len(fd)}"
fd["row_name"] = ["gene_" + str(i) for i in range(1, len(fd) + 1)]
sym2row = dict(zip(fd["Symbol"].astype(str), fd["row_name"]))
zp3_row = sym2row.get("ZP3")
assert zp3_row is not None, "ZP3 symbol not found in fData"
zp3 = ex.loc[zp3_row].astype(float)
# sanity: ZP3 should not be all-zero and should have a plausible dynamic range
assert (zp3 > 0).mean() > 0.5, f"ZP3 row suspiciously sparse: {(zp3 > 0).mean():.2%} nonzero"
# cross-check with previously frozen results (median 441 / 461)
resp_map = dict(zip(pdt["X"], pdt["binaryResponse"]))
mz_r = np.median(zp3[[c for c in ex.columns if resp_map.get(c) == "CR/PR"]])
mz_n = np.median(zp3[[c for c in ex.columns if resp_map.get(c) == "SD/PD"]])
print(f"[IMvigor210] ZP3 row '{zp3_row}' validated; medians resp={mz_r:.0f} non={mz_n:.0f} "
      f"(previous frozen: 441 / 461)")

# --- primary analysis -------------------------------------------------
zp3_r = zp3[[c for c in ex.columns if resp_map.get(c) == "CR/PR"]]
zp3_n = zp3[[c for c in ex.columns if resp_map.get(c) == "SD/PD"]]
u21, p21, rb21, se21, lo21, hi21 = rank_biserial(zp3_r, zp3_n)
m_r21, q1_r21, q3_r21 = med_iqr(zp3_r)
m_n21, q1_n21, q3_n21 = med_iqr(zp3_n)
rows.append({
    "Cohort": "IMvigor210", "Tissue": "Urothelial carcinoma", "Treatment": "Anti-PD-L1",
    "Scope": "primary", "n_resp": len(zp3_r), "n_non": len(zp3_n),
    "med_resp": round(m_r21, 1), "med_non": round(m_n21, 1),
    "U": round(u21, 2), "P": p21, "r_rb": rb21, "SE": se21,
    "CI_low": lo21, "CI_high": hi21, "CI_str": fmt_ci(lo21, hi21),
})
print(f"[IMvigor210] n={len(zp3_r)+len(zp3_n)} (R {len(zp3_r)} / NR {len(zp3_n)}); "
      f"med {m_r21:.1f} vs {m_n21:.1f}; U={u21:.0f} P={p21:.4f} r_rb={rb21:+.3f} "
      f"CI ({lo21:.2f}, {hi21:.2f})")

# ============================================================
# 3. GSE91061 (melanoma, anti-PD-1) — Entrez 7784 (true ZP3), Pre baseline
# ============================================================
g9 = os.path.join(ROOT, "output", "gse91061_validation")
expr = pd.read_csv(os.path.join(g9, "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"), index_col=0)
ZP3_ENTREZ = 7784   # corrected (was erroneously 8277 = SP5)
assert ZP3_ENTREZ in expr.index, "7784 not in GSE91061 FPKM"
zp3_expr = expr.loc[ZP3_ENTREZ]

# clinical parsing from series matrix (same fields as v2 script)
st = st2 = vs = rs = []
for line in open(os.path.join(g9, "GSE91061_series_matrix.txt"), encoding="utf-8"):
    if line.startswith("!Sample_title") and not st:
        st = [p.strip('"') for p in line.strip().split("\t")[1:]]
    if line.startswith("!Sample_characteristics_ch1") and "visit" in line and not vs:
        vs = [p.split(": ")[-1].strip('"') for p in line.strip().split("\t")[1:] if ": " in p]
    if line.startswith("!Sample_characteristics_ch1") and "response" in line and not rs:
        rs = [p.split(": ")[-1].strip('"') for p in line.strip().split("\t")[1:] if ": " in p]
cl = pd.DataFrame({"sample_title": st, "visit": vs, "response": rs})
rmap = {"PRCR": "Response", "PR": "Response", "CR": "Response",
        "SD": "Non-response", "PD": "Non-response", "UNK": "Unknown"}
cl["binary_response"] = cl["response"].map(rmap)
cl["zp3"] = cl["sample_title"].map(zp3_expr)
cl = cl.dropna(subset=["zp3"])

sens_rows = []
for scope, mask in [("Pre (baseline)", cl["visit"] == "Pre"),
                    ("All samples", pd.Series(True, index=cl.index))]:
    d = cl[mask]
    d = d[d["binary_response"].isin(["Response", "Non-response"])]
    rr = d[d["binary_response"] == "Response"]["zp3"]
    nn = d[d["binary_response"] == "Non-response"]["zp3"]
    if len(rr) >= 3 and len(nn) >= 3:
        u, p, rb, se, lo, hi = rank_biserial(rr, nn)
        m_r, _, _ = med_iqr(rr)
        m_n, _, _ = med_iqr(nn)
        sens_rows.append({
            "Cohort": "GSE91061", "Scope": scope, "n_resp": len(rr), "n_non": len(nn),
            "med_resp": round(m_r, 2), "med_non": round(m_n, 2),
            "U": round(u, 2), "P": p, "r_rb": rb, "SE": se,
            "CI_low": lo, "CI_high": hi, "CI_str": fmt_ci(lo, hi),
        })
        print(f"[GSE91061:{scope}] n={len(rr)+len(nn)} (R {len(rr)} / NR {len(nn)}); "
              f"med {m_r:.2f} vs {m_n:.2f}; U={u:.0f} P={p:.4f} r_rb={rb:+.3f} "
              f"CI ({lo:.2f}, {hi:.2f})")

# primary = Pre baseline
prim9 = sens_rows[0]
rows.append({
    "Cohort": "GSE91061", "Tissue": "Melanoma", "Treatment": "Anti-PD-1",
    "Scope": "primary", "n_resp": prim9["n_resp"], "n_non": prim9["n_non"],
    "med_resp": prim9["med_resp"], "med_non": prim9["med_non"],
    "U": prim9["U"], "P": prim9["P"], "r_rb": prim9["r_rb"], "SE": prim9["SE"],
    "CI_low": prim9["CI_low"], "CI_high": prim9["CI_high"], "CI_str": prim9["CI_str"],
})

# ============================================================
# Save frozen tables
# ============================================================
frz = pd.DataFrame(rows)
frz.to_csv(OUT, index=False)
pd.DataFrame(sens_rows).to_csv(OUT_S, index=False)
print(f"\nsaved -> {OUT}  ({len(frz)} rows)")
print(f"saved -> {OUT_S}  ({len(sens_rows)} rows)")

# ============================================================
# Power table (MDE at 80% power ~ 2.80 * SE; TOST for IMvigor210 only)
# ============================================================
def mde_80(se):
    return 2.80 * se

power_rows = []
for _, r_ in frz.iterrows():
    mde = mde_80(r_["SE"])
    # post-hoc power for observed r (two-sided, alpha 0.05)
    z = abs(r_["r_rb"]) / r_["SE"]
    posthoc = 1.0 - stats.norm.cdf(1.96 - z) + stats.norm.cdf(-1.96 - z)
    row = {"Cohort": r_["Cohort"], "n1": int(r_["n_resp"]), "n2": int(r_["n_non"]),
           "r_obs": round(r_["r_rb"], 4), "se": round(r_["SE"], 4),
           "MDE_80": round(mde, 3), "posthoc_power": round(posthoc, 4)}
    if r_["Cohort"] == "IMvigor210":
        # TOST equivalence vs |r| < 0.15 (pre-specified), using z-test on r/se
        bound = 0.15
        z1 = (bound - abs(r_["r_rb"])) / r_["SE"]
        z2 = (bound + abs(r_["r_rb"])) / r_["SE"]
        p_tost = 1.0 - stats.norm.cdf(z1)  # conservative: lower bound drives
        row.update({"TOST_bound": bound, "TOST_p": round(p_tost, 4),
                    "TOST_equivalent": bool(p_tost < 0.05)})
    power_rows.append(row)

pwr = pd.DataFrame(power_rows)
pwr.to_csv(POWER, index=False)
print("\nupdated power table ->", POWER)
print(pwr.to_string(index=False))
