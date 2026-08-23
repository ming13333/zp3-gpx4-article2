# -*- coding: utf-8 -*-
"""
A2 P1 阳性对照最终版: 金标准阳性对照 (TMB, PD-L1 IHC IC score) + 表达类标志物
=============================================================================
方向修正: rank-biserial r = 2U/(n1*n2) - 1, U 为 responder 组的 U 统计量。
r>0 表示 responder 组更高。
结果: 两个金标准标志物 (TMB P<0.0001; IC score P=0.0012) 均显著 -> pipeline 灵敏度验证通过。
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

IM = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output", "immunotherapy_validation")
PD = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output", "phase1_knowledge_gap_filling")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article2", "results", "positive_control_results.csv")

def rank_biserial(x_resp, x_non):
    x_resp = np.asarray(x_resp, float); x_non = np.asarray(x_non, float)
    n1, n2 = len(x_resp), len(x_non)
    u, p = stats.mannwhitneyu(x_resp, x_non, alternative="two-sided")
    r = 2 * u / (n1 * n2) - 1
    return u, p, r

rows = []

# ---- IMvigor210: TMB (+direction check) ----
pd_ = pd.read_csv(os.path.join(IM, "pData_IMvigor210.csv"))
pd_.columns = [c.strip('"') for c in pd_.columns]
d = pd_[["X", "binaryResponse", "FMOne.mutation.burden.per.MB", "IC.Level"]].dropna(subset=["binaryResponse"]).copy()
d["tmb"] = pd.to_numeric(d["FMOne.mutation.burden.per.MB"], errors="coerce")
rg = d[d["binaryResponse"] == "CR/PR"]["tmb"].dropna()
ng = d[d["binaryResponse"] != "CR/PR"]["tmb"].dropna()
u, p, r = rank_biserial(rg, ng)
rows.append({"Cohort": "IMvigor210", "Marker": "TMB (per MB)", "Type": "gold-standard (genomic)",
             "n_resp": len(rg), "n_non": len(ng),
             "median_resp": round(rg.median(), 2), "median_non": round(ng.median(), 2),
             "P": p, "rank_biserial_r": round(r, 3), "significant": p < 0.05})

# ---- IMvigor210: PD-L1 IHC (IC score) ----
d2 = d.dropna(subset=["IC.Level"]).copy()
d2["IC_high"] = d2["IC.Level"].isin(["IC2", "IC2+"])
tab = pd.crosstab(d2["IC_high"], d2["binaryResponse"])
from scipy.stats import chi2_contingency
chi2, p_ic, dof, _ = chi2_contingency(tab)
resp_ic2 = tab.loc[True, "CR/PR"]; tot_ic2 = tab.loc[True].sum()
resp_ic01 = tab.loc[False, "CR/PR"]; tot_ic01 = tab.loc[False].sum()
rows.append({"Cohort": "IMvigor210", "Marker": "PD-L1 IHC (IC2+ vs IC0/1)", "Type": "gold-standard (protein IHC)",
             "n_resp": int(resp_ic2), "n_non": int(tot_ic2 - resp_ic2),
             "median_resp": round(100 * resp_ic2 / tot_ic2, 1), "median_non": round(100 * resp_ic01 / tot_ic01, 1),
             "P": p_ic, "rank_biserial_r": np.nan, "significant": p_ic < 0.05})

# ---- IMvigor210: CD274 expression ----
import pandas as pd
ex = pd.read_csv(os.path.join(IM, "exmat_censored_IMvigor210.csv"), index_col=0)
fd = pd.read_csv(os.path.join(IM, "fData_IMvigor210.csv"))
ex.columns = [c.strip('"') for c in ex.columns]
ex.index = [str(i).strip('"') for i in ex.index]
fd["row_name"] = ["gene_" + str(i) for i in range(1, len(fd) + 1)]
sym2row = dict(zip(fd["Symbol"].astype(str), fd["row_name"]))
resp_map = dict(zip(pd_["X"], pd_["binaryResponse"]))
sam_ids = [c for c in ex.columns if c in resp_map]
cd274_r, cd274_n, cxcl9_r, cxcl9_n = [], [], [], []
for s in sam_ids:
    if pd.isna(resp_map[s]): continue
    v_cd = ex.loc[sym2row.get("CD274", ""), s]; v_cx = ex.loc[sym2row.get("CXCL9", ""), s]
    if resp_map[s] == "CR/PR":
        cd274_r.append(v_cd); cxcl9_r.append(v_cx)
    else:
        cd274_n.append(v_cd); cxcl9_n.append(v_cx)
u, p, r = rank_biserial(cd274_r, cd274_n)
rows.append({"Cohort": "IMvigor210", "Marker": "CD274 (PD-L1 mRNA)", "Type": "expression",
             "n_resp": len(cd274_r), "n_non": len(cd274_n),
             "median_resp": round(float(np.median(cd274_r)), 2), "median_non": round(float(np.median(cd274_n)), 2),
             "P": p, "rank_biserial_r": round(r, 3), "significant": p < 0.05})

# ---- GSE91061: IFN-g / CD8 scores (Pre subset; molecular subtype file) ----
sub = pd.read_csv(os.path.join(PD, "gse91061_molecular_subtype_clinical.csv"))
sub = sub[(sub["visit"] == "Pre") & (sub["response"] != "UNK")].copy()
sub["is_resp"] = sub["response"] == "PRCR"
rg2, ng2 = sub[sub["is_resp"]], sub[~sub["is_resp"]]
for mcol, label in [("IFN_gamma_score", "IFN-g score (Ayers)"), ("CD8_score", "CD8 score"),
                    ("expr_CD274", "CD274 (PD-L1 mRNA)"), ("Tcell_score", "T-cell score")]:
    u, p, r = rank_biserial(rg2[mcol], ng2[mcol])
    rows.append({"Cohort": "GSE91061", "Marker": label, "Type": "expression/signature",
                 "n_resp": len(rg2), "n_non": len(ng2),
                 "median_resp": round(float(rg2[mcol].median()), 2), "median_non": round(float(ng2[mcol].median()), 2),
                 "P": p, "rank_biserial_r": round(r, 3), "significant": p < 0.05})

# ---- GSE78220: CD274 + IFN-g signature ----
H2 = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output", "h2_bulk")
IFNG_GENES = ["CXCL9", "CXCL10", "CCL5", "IFNG", "IDO1", "STAT1", "PRF1", "GZMA", "GZMB", "HLA-DRA"]
expr_df = pd.read_excel(os.path.join(H2, "GSE78220_PatientFPKM.xlsx"), sheet_name=0, index_col=0)
resp = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article2", "results", "gse78220_zp3_response.csv"))
resp = resp[resp["timepoint"] == "baseline"].copy()
resp["is_resp"] = resp["response"].isin(["CR", "PR"])
avail = [c for c in expr_df.columns if c in set(resp["sample_id"])]
sub3 = resp[resp["sample_id"].isin(avail)].copy()
sig = expr_df.reindex(IFNG_GENES).loc[:, sub3["sample_id"]].T.fillna(expr_df.reindex(IFNG_GENES).mean(axis=1))
sub3 = sub3.assign(ifng_sig=sig.mean(axis=1).values)
rg3, ng3 = sub3[sub3["is_resp"]], sub3[~sub3["is_resp"]]
u, p, r = rank_biserial(rg3["pdl1_expr"], ng3["pdl1_expr"])
rows.append({"Cohort": "GSE78220", "Marker": "CD274 (PD-L1 mRNA)", "Type": "expression",
             "n_resp": len(rg3), "n_non": len(ng3),
             "median_resp": round(float(rg3["pdl1_expr"].median()), 2), "median_non": round(float(ng3["pdl1_expr"].median()), 2),
             "P": p, "rank_biserial_r": round(r, 3), "significant": p < 0.05})
u, p, r = rank_biserial(rg3["ifng_sig"], ng3["ifng_sig"])
rows.append({"Cohort": "GSE78220", "Marker": "IFN-g signature (Ayers 10g)", "Type": "signature",
             "n_resp": len(rg3), "n_non": len(ng3),
             "median_resp": round(float(rg3["ifng_sig"].median()), 2), "median_non": round(float(ng3["ifng_sig"].median()), 2),
             "P": p, "rank_biserial_r": round(r, 3), "significant": p < 0.05})

res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print("=== A2 positive control results ===")
print(res.to_string(index=False))
print()
n_sig = res["significant"].sum()
print(f"Total: {len(res)}, significant: {n_sig}")
gold = res[res["Type"].str.contains("gold-standard")]
print(f"Gold-standard significant: {gold['significant'].sum()}/{len(gold)}")
print("=> Pipeline sensitivity: ", "PASS" if n_sig >= 1 else "FAIL")