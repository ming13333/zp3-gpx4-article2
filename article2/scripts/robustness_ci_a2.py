#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A2 reviewer revision - robustness analysis (2026-08-17)
======================================
Response to senior reviewer concerns:
  R1. Whether the approximate normal CI for rank-biserial r is robust
      -> stratified bootstrap 95% CI with replacement (B=20000, sampling each cohort separately)
  R2. Whether logistic regression on GSE91061 is stable with only 10 responders
      -> Firth penalized-likelihood logistic (univariable + adjusted for IFN-gamma/immunosuppression)
  R3. GSE91061 sensitivity analysis n=105 samples contain multiple visit samples from the same patient (pseudo-replication)
      -> patient-level aggregated sensitivity (take 1 sample per patient, Pre preferred; 62 patients total)
      + explicitly count 105 samples -> 62 patients, 42 patients with two visits

Output (article2/results/):
  a2_robustness_bootstrap.csv   (bootstrap CI for three cohorts)
  a2_gse91061_firth.csv         (Firth logistic results)
  a2_gse91061_patient_level.csv (patient-level sensitivity)
"""
import os
import re
import math
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "article2", "results")
DATA = os.path.join(ROOT, "output")
rng = np.random.default_rng(20260817)
B = 20000


def rank_biserial(x_resp, x_non):
    x_resp = np.asarray(x_resp, float)
    x_non = np.asarray(x_non, float)
    n1, n2 = len(x_resp), len(x_non)
    u, p = stats.mannwhitneyu(x_resp, x_non, alternative="two-sided")
    r = 2.0 * u / (n1 * n2) - 1.0
    se = 2.0 * math.sqrt((n1 + n2 + 1.0) / (12.0 * n1 * n2))
    return u, p, r, se


def bootstrap_ci(x_resp, x_non, B=B, seed=20260817):
    """Stratified bootstrap 95% CI for rank-biserial r (2.5/97.5 percentiles)."""
    rng = np.random.default_rng(seed)
    x_resp = np.asarray(x_resp, float)
    x_non = np.asarray(x_non, float)
    n1, n2 = len(x_resp), len(x_non)
    rs = np.empty(B)
    for i in range(B):
        a = rng.choice(x_resp, size=n1, replace=True)
        b = rng.choice(x_non, size=n2, replace=True)
        u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
        rs[i] = 2.0 * u / (n1 * n2) - 1.0
    return np.percentile(rs, [2.5, 97.5])


# ================================================================
# Load main analysis data for the three cohorts (same source as recompute_effects_a2.py)
# ================================================================

# --- GSE78220 (baseline) ---
df20 = pd.read_csv(os.path.join(OUT, "gse78220_zp3_response.csv"))
df20 = df20[df20["timepoint"] == "baseline"].copy()
df20["is_resp"] = df20["response"].isin(["CR", "PR"])
r20 = df20[df20["is_resp"]]["zp3_expr"].astype(float).values
n20 = df20[~df20["is_resp"]]["zp3_expr"].astype(float).values

# --- IMvigor210 ---
im = os.path.join(DATA, "immunotherapy_validation")
ex = pd.read_csv(os.path.join(im, "exmat_censored_IMvigor210.csv"), index_col=0)
fd = pd.read_csv(os.path.join(im, "fData_IMvigor210.csv"))
pdt = pd.read_csv(os.path.join(im, "pData_IMvigor210.csv"))
ex.columns = [str(c).strip('"') for c in ex.columns]
fd.columns = [str(c).strip('"') for c in fd.columns]
pdt.columns = [str(c).strip('"') for c in pdt.columns]
fd["row_name"] = ["gene_" + str(i) for i in range(1, len(fd) + 1)]
sym2row = dict(zip(fd["Symbol"].astype(str), fd["row_name"]))
zp3_row = sym2row.get("ZP3")
assert zp3_row is not None, "ZP3 not in IMvigor210 fData"
zp3 = ex.loc[zp3_row].astype(float)
resp_map = dict(zip(pdt["X"], pdt["binaryResponse"]))
r21 = zp3[[c for c in ex.columns if resp_map.get(c) == "CR/PR"]].astype(float).values
n21 = zp3[[c for c in ex.columns if resp_map.get(c) == "SD/PD"]].astype(float).values

# --- GSE91061 (Pre baseline, Entrez 7784) ---
g9 = os.path.join(DATA, "gse91061_validation")
expr = pd.read_csv(os.path.join(g9, "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"), index_col=0)
assert 7784 in expr.index, "7784 not in GSE91061 FPKM"
zp3_expr = expr.loc[7784]
st = vs = rs = []
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
pre9 = cl[(cl["visit"] == "Pre") & cl["binary_response"].isin(["Response", "Non-response"])]
r91 = pre9[pre9["binary_response"] == "Response"]["zp3"].astype(float).values
n91 = pre9[pre9["binary_response"] == "Non-response"]["zp3"].astype(float).values

# ================================================================
# R1. Bootstrap CI (three cohorts)
# ================================================================
print("=" * 70)
print("R1. Bootstrap 95% CI for rank-biserial r (B=20000, stratified)")
print("=" * 70)
boot_rows = []
for name, rr, nn in [("GSE78220", r20, n20), ("IMvigor210", r21, n21), ("GSE91061_Pre", r91, n91)]:
    _, p, r, se = rank_biserial(rr, nn)
    lo_a, hi_a = r - 1.96 * se, r + 1.96 * se
    lo_b, hi_b = bootstrap_ci(rr, nn)
    boot_rows.append({"Cohort": name, "n_resp": len(rr), "n_non": len(nn),
                      "r_rb": round(r, 4), "SE": round(se, 4),
                      "approx_CI_low": round(lo_a, 4), "approx_CI_high": round(hi_a, 4),
                      "boot_CI_low": round(lo_b, 4), "boot_CI_high": round(hi_b, 4),
                      "P": round(p, 4)})
    print(f"{name:14s} n=({len(rr)},{len(nn)})  r={r:+.3f}  approx CI ({lo_a:.2f},{hi_a:.2f})  "
          f"bootstrap CI ({lo_b:.2f},{hi_b:.2f})")

# ================================================================
# R3. GSE91061 patient-level sensitivity (deduplicate first, then aggregate)
# ================================================================
print()
print("=" * 70)
print("R3. GSE91061 patient-level sensitivity (1 sample per patient, Pre preferred)")
print("=" * 70)
all_ev = cl[cl["binary_response"].isin(["Response", "Non-response"])].copy()
all_ev["patient"] = all_ev["sample_title"].str.extract(r"(Pt\d+)", expand=False)
n_patients = all_ev["patient"].nunique()
n_samples = len(all_ev)
both = all_ev.groupby("patient")["visit"].apply(lambda v: set(v) == {"Pre", "On"})
print(f"All-sample sensitivity: {n_samples} samples -> {n_patients} unique patients; "
      f"{int(both.sum())} patients with both Pre+On visit samples")

# patient-level: Pre preferred, On if no Pre (consistent with one-sample-per-patient principle)
pl = all_ev.sort_values(["patient", "visit"], key=lambda s: s.map({"Pre": 0, "On": 1}) if s.name == "visit" else s)
pl = pl.groupby("patient").first().reset_index()
r_pl = pl[pl["binary_response"] == "Response"]["zp3"].astype(float).values
n_pl = pl[pl["binary_response"] == "Non-response"]["zp3"].astype(float).values
u_pl, p_pl, r_pl_rb, se_pl = rank_biserial(r_pl, n_pl)
lo_pl, hi_pl = r_pl_rb - 1.96 * se_pl, r_pl_rb + 1.96 * se_pl
print(f"patient-level: n={len(pl)} (R {len(r_pl)} / NR {len(n_pl)}); "
      f"U={u_pl:.0f} P={p_pl:.4f} r_rb={r_pl_rb:+.3f} CI ({lo_pl:.2f},{hi_pl:.2f})")

# ================================================================
# R2. GSE91061 Firth penalized logistic (univariable + adjusted model)
# ================================================================
print()
print("=" * 70)
print("R2. GSE91061 Firth penalized logistic (Pre n=49, 10 responders)")
print("=" * 70)

# Immune features (IFN-gamma Ayers / immunosuppression score) —— read from molecular subtype file or immune feature file
def load_g9_immune():
    """Attempt to read IFN-gamma/CD8 scores from molecular-subtype output; if absent, compute from expression matrix
    Ayers IFN-gamma: CXCL9, CXCL10, CCL5, IFNG, IDO1, STAT1, PRF1, GZMA, GZMB, HLA-DRA (mean of 10 genes)
    immunosuppression score: TGFB1, IL10, CD274, PDCD1LG2, CTLA4, LAG3, HAVCR2, PDCD1, TIGIT, CD163, MRC1, FOXP3
    """
    genes_ifn = ["CXCL9", "CXCL10", "CCL5", "IFNG", "IDO1", "STAT1", "PRF1", "GZMA", "GZMB", "HLA-DRA"]
    genes_immuno = ["TGFB1", "IL10", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "PDCD1", "TIGIT", "CD163", "MRC1", "FOXP3"]
    # Need symbol->Entrez mapping: GSE91061 FPKM uses Entrez row names, use known common Entrez mapping
    symbol2entrez = {}
    # From fData? GSE91061 has no fData; use symbol->entrez mapping from IMvigor210 fData
    try:
        fd2 = pd.read_csv(os.path.join(im, "fData_IMvigor210.csv"))
        symbol2entrez = dict(zip(fd2["Symbol"].astype(str), fd2["EntrezGene"].astype(str)))
    except Exception:
        pass
    # Supplement common mappings (Ayers 10 genes + immunosuppressive 12 genes)
    extra = {
        "CXCL9": 4283, "CXCL10": 3627, "CCL5": 6352, "IFNG": 3458, "IDO1": 3620,
        "STAT1": 6772, "PRF1": 5551, "GZMA": 3001, "GZMB": 3002, "HLA-DRA": 3122,
        "TGFB1": 7040, "IL10": 3586, "CD274": 29126, "PDCD1LG2": 80380, "CTLA4": 1493,
        "LAG3": 3902, "HAVCR2": 84868, "PDCD1": 5133, "TIGIT": 201633, "CD163": 9332,
        "MRC1": 4360, "FOXP3": 50943,
    }
    for g, e in extra.items():
        symbol2entrez.setdefault(g, str(e))
    found = {}
    for g in genes_ifn + genes_immuno:
        e = symbol2entrez.get(g)
        if e is not None:
            try:
                e_int = int(e)  # expr index is int64; string will fail to match
            except (TypeError, ValueError):
                e_int = e
            if e_int in expr.index:
                found[g] = expr.loc[e_int].astype(float)
    # Sample alignment
    df = pre9.set_index("sample_title")
    out_cols = {}
    for g in genes_ifn + genes_immuno:
        if g in found:
            # reindex aligns by sample_title (found[g] index is also sample_title)
            out_cols[g] = found[g].reindex(df.index)
    out_cols["zp3"] = df["zp3"]
    out_cols["y"] = (df["binary_response"] == "Response").astype(int)
    out_cols["ifn"] = pd.DataFrame({g: out_cols[g] for g in genes_ifn if g in out_cols}).mean(axis=1)
    out_cols["immuno"] = pd.DataFrame({g: out_cols[g] for g in genes_immuno if g in out_cols}).mean(axis=1)
    return pd.DataFrame(out_cols)


def firth_logit(y, X, ridge=1e-6):
    """Firth penalized logistic regression with tiny ridge for numerical stability.
    X contains intercept column. Returns (beta, cov, n_iter)."""
    n, p = X.shape
    beta = np.zeros(p)
    for it in range(200):
        eta = X @ beta
        eta = np.clip(eta, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-eta))
        mu = np.clip(mu, 1e-10, 1 - 1e-10)
        w = mu * (1 - mu)
        XWX = (X.T * w) @ X + ridge * np.eye(p)
        try:
            XWX_inv = np.linalg.inv(XWX)
        except np.linalg.LinAlgError:
            XWX_inv = np.linalg.pinv(XWX)
        h = np.sum((X @ XWX_inv) * X, axis=1) * w
        U_star = X.T @ (y - mu + h * (0.5 - mu))
        db = XWX_inv @ U_star
        beta_new = beta + db
        if np.max(np.abs(db)) < 1e-8:
            beta = beta_new
            break
        beta = beta_new
    return beta, XWX_inv, it + 1


firth_rows = []
d9 = load_g9_immune()
y = d9["y"].values
for model_name, cols in [("univariable", ["zp3"]),
                          ("+IFN-gamma score", ["zp3", "ifn"]),
                          ("+immunosuppression score", ["zp3", "immuno"])]:
    X = np.column_stack([np.ones(len(d9))] + [d9[c].values for c in cols])
    # Standardize continuous variables (zp3/ifn/immuno) to avoid scale issues
    Xn = X.copy()
    for j in range(1, X.shape[1]):
        s = Xn[:, j].std()
        if s > 0:
            Xn[:, j] = (Xn[:, j] - Xn[:, j].mean()) / s
    beta, cov, n_iter = firth_logit(y, Xn)
    se_b = np.sqrt(np.diag(cov))
    z = beta / se_b
    pvals = 2 * (1 - stats.norm.cdf(np.abs(z)))
    for j, cname in enumerate(["(Intercept)"] + cols):
        if cname == "zp3":
            or_ = np.exp(beta[j])
            lo_or = np.exp(beta[j] - 1.96 * se_b[j])
            hi_or = np.exp(beta[j] + 1.96 * se_b[j])
            firth_rows.append({"Model": model_name, "Term": cname,
                               "n": len(d9), "events": int(y.sum()),
                               "Firth_OR": round(or_, 3),
                               "Firth_CI_low": round(lo_or, 3), "Firth_CI_high": round(hi_or, 3),
                               "Firth_P": round(pvals[j], 4), "iter": n_iter})
            print(f"[Firth] {model_name:26s} ZP3 OR={or_:.2f} ({lo_or:.2f}-{hi_or:.2f}) "
                  f"P={pvals[j]:.3f} (n={len(d9)}, events={int(y.sum())}, iter={n_iter})")

# ================================================================
# Save
# ================================================================
pd.DataFrame(boot_rows).to_csv(os.path.join(OUT, "a2_robustness_bootstrap.csv"), index=False)
pd.DataFrame(firth_rows).to_csv(os.path.join(OUT, "a2_gse91061_firth.csv"), index=False)
pd.DataFrame({
    "Cohort": ["GSE91061"],
    "Scope": ["patient-level (Pre-first, n samples -> n patients)"],
    "n_samples": [n_samples], "n_patients": [n_patients],
    "patients_both_visits": [int(both.sum())],
    "n_resp": [len(r_pl)], "n_non": [len(n_pl)],
    "U": [round(u_pl, 2)], "P": [round(p_pl, 4)], "r_rb": [round(r_pl_rb, 4)],
    "CI_low": [round(lo_pl, 4)], "CI_high": [round(hi_pl, 4)],
}).to_csv(os.path.join(OUT, "a2_gse91061_patient_level.csv"), index=False)
print(f"\nsaved -> {os.path.join(OUT, 'a2_robustness_bootstrap.csv')}")
print(f"saved -> {os.path.join(OUT, 'a2_gse91061_firth.csv')}")
print(f"saved -> {os.path.join(OUT, 'a2_gse91061_patient_level.csv')}")
