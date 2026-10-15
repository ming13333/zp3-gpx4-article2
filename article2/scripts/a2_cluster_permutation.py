# -*- coding: utf-8 -*-
"""
A2 决定性检验：合并框架 B 的名义 P 值，在"尊重患者聚类"的零分布下是否成立？
用患者层面标签置换（cluster permutation），保留每位患者的样本集合不变。
- 若置换 P 远大于 0.008 ⇒ 名义 P 反保守（伪重复确实使 P 值失真）
- 若置换 P ≈ 0.008     ⇒ 聚类不是问题，合并信号是真实的
同时对框架 A（49 患者 1:1，无聚类）做同样置换作为管线校验。
"""
import csv, json, os, sys
from collections import defaultdict
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
# --- portable paths: everything resolves relative to this file's repository ---
_HERE = os.path.dirname(os.path.abspath(__file__))   # <repo>/article2/scripts
A2    = os.path.dirname(_HERE)                       # <repo>/article2
ROOT  = os.path.dirname(A2)                          # <repo>
FPKM  = os.environ.get("A2_GSE91061_FPKM", os.path.join(
        ROOT, "output", "gse91061_validation",
        "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"))
MASK  = os.path.join(A2, "results", "framework_cohort_masks_v2.json")
OUT   = os.path.join(A2, "results", "a2_cluster_permutation_20261014.json")
if not os.path.exists(FPKM):
    raise SystemExit(
        "Input FPKM matrix not found: %s\n"
        "Obtain GSE91061 from GEO and place "
        "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv at that path, or set A2_GSE91061_FPKM." % FPKM)
ZP3_ENTREZ = "7784"

def load_zp3():
    with open(FPKM, newline="", encoding="utf-8") as f:
        r = csv.reader(f); hdr = next(r); samples = hdr[1:]
        for row in r:
            if row and row[0].strip().strip('"') == ZP3_ENTREZ:
                return dict(zip(samples, [float(x) for x in row[1:]]))
    raise SystemExit("ZP3 not found")

def rr(gr, gn):
    u = stats.mannwhitneyu(gr, gn, alternative="two-sided")
    return 2*u.statistic/(len(gr)*len(gn)) - 1, u.pvalue, u.statistic

zp3 = load_zp3()
masks = json.load(open(MASK, encoding="utf-8"))
labels = masks["labels"]; A = masks["A_pre"]; B = masks["B_naive_pooled"]; C = masks["C_patient_pre_else_on"]

rng = np.random.default_rng(20261014)
res = {}

def cluster_perm(samples, n_perm=20000, tag=""):
    """患者层面标签置换；统计量=样本层 Mann-Whitney U（稿件的统计量）"""
    bypt = defaultdict(list)
    for s in samples: bypt[s.split("_")[0]].append(s)
    pids = list(bypt.keys())
    plab = {p: labels[bypt[p][0]] for p in pids}
    nR = sum(1 for p in pids if plab[p] == "R")
    obs_R = [zp3[s] for s in samples if labels[s] == "R"]
    obs_N = [zp3[s] for s in samples if labels[s] == "N"]
    r_obs, p_obs, u_obs = rr(obs_R, obs_N)
    # 观测统计量用 |r|
    null = np.empty(n_perm)
    idx = np.arange(len(pids))
    for i in range(n_perm):
        perm = rng.permutation(idx)
        Rset = set(idx[perm][:nR])
        gR, gN = [], []
        for k, p in enumerate(pids):
            (gR if k in Rset else gN).extend(zp3[s] for s in bypt[p])
        if len(gR) < 2 or len(gN) < 2:
            null[i] = 0.0; continue
        u = stats.mannwhitneyu(gR, gN, alternative="two-sided").statistic
        null[i] = 2*u/(len(gR)*len(gN)) - 1
    p_perm = float(np.mean(np.abs(null) >= abs(r_obs)))
    print(f"\n  [{tag}] 样本数={len(samples)}  患者数={len(pids)} (R患者={nR})")
    print(f"      样本层 MWU（稿件统计量）: r={r_obs:+.4f}  名义 P={p_obs:.4f}")
    print(f"      患者聚类置换 P（{n_perm}次）: {p_perm:.4f}")
    print(f"      零分布 |r| 95分位 = {np.percentile(np.abs(null),95):.4f}")
    return dict(n_samples=len(samples), n_patients=len(pids), nR_patients=nR,
                r=r_obs, nominal_P=p_obs, cluster_perm_P=p_perm,
                null_absr_p95=float(np.percentile(np.abs(null),95)))

print("="*96)
print("【检验1】框架 B（合并 105 样本 / 62 患者）—— 聚类置换 vs 名义 P")
res["B"] = cluster_perm(B, tag="Framework B pooled")

print("\n"+"="*96)
print("【检验2】框架 A（49 预处理样本 / 49 患者，1:1 无聚类）—— 管线校验")
print("     预期：置换 P ≈ 名义 P（因为无聚类，参数近似成立）")
res["A"] = cluster_perm(A, tag="Framework A pre-only")

print("\n"+"="*96)
print("【检验3】on-treatment-only 层（有聚类）")
on = [s for s in B if "_On_" in s]
res["On"] = cluster_perm(on, tag="On-treatment only")

# 均值聚合（另一种患者层汇总，与稿件 C 的"优先取 Pre"不同）
print("\n"+"="*96)
print("【检验4】患者层均值聚合（Pre+On 取均值；对照稿件 C 的'优先取 Pre'）")
bypt = defaultdict(dict)
for s in B: bypt[s.split("_")[0]]["On" if "_On_" in s else "Pre"] = zp3[s]
lab_pt = {s.split("_")[0]: labels[s] for s in B}
gR = [np.mean(list(v.values())) for k,v in bypt.items() if lab_pt[k]=="R"]
gN = [np.mean(list(v.values())) for k,v in bypt.items() if lab_pt[k]=="N"]
r_m, p_m, _ = rr(gR, gN)
print(f"  均值聚合 (n={len(bypt)}; R={len(gR)}, N={len(gN)}): r={r_m:+.4f} P={p_m:.4f}")
print(f"  （对照：稿件 C 患者层 n=62 P=0.1753 / 本次患者层 n={len(bypt)}）")
res["C_mean_aggregate"] = {"n": len(bypt), "nR": len(gR), "nN": len(gN), "r": r_m, "P": p_m}

json.dump(res, open(OUT,"w",encoding="utf-8"), ensure_ascii=False, indent=1)
print("\nsaved:", OUT)
