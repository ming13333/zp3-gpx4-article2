# -*- coding: utf-8 -*-
"""
决定性检验：稿件的头条证据"832 -> 4334 -> 843"能否成立？
对每个框架建立「患者聚类置换」零分布（保留患者内重复结构 + 样本构成），
得到该设计下"无关联时"应有的显著基因数，与观测数比较。

方法：秩在固定样本集上只算一次；置换只改变患者标签 → 组内秩和。
用带结值校正的正态近似（结值结构在固定样本集下不变，可预计算）。
"""
import csv, json, os, sys
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
OUT   = os.path.join(A2, "results", "a2_transcriptome_cluster_null_20261014.json")
if not os.path.exists(FPKM):
    raise SystemExit(
        "Input FPKM matrix not found: %s\n"
        "Obtain GSE91061 from GEO and place "
        "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv at that path, or set A2_GSE91061_FPKM." % FPKM)

print("loading FPKM ...")
hdr = None; rows = []; ids = []
with open(FPKM, newline="", encoding="utf-8") as f:
    r = csv.reader(f); hdr = next(r); samples_all = [h.strip() for h in hdr[1:]]
    for row in r:
        if not row: continue
        ids.append(row[0].strip().strip('"')); rows.append([float(x) for x in row[1:]])
X = np.array(rows, dtype=np.float64)          # genes x samples
print("matrix:", X.shape, "samples:", len(samples_all))
sid = {s: i for i, s in enumerate(samples_all)}

masks = json.load(open(MASK, encoding="utf-8"))
labels = masks["labels"]
A = masks["A_pre"]; B = masks["B_naive_pooled"]; C = masks["C_patient_pre_else_on"]
ON = [s for s in B if "_On_" in s]

def framework_null(cols, perm_labels=True, n_perm=200, seed=20261014):
    """cols: 样本名列表（决定用哪些列）。返回观测显著数与置换零分布。"""
    idx = [sid[c] for c in cols]
    sub = X[:, idx]                                  # genes x n
    n = sub.shape[1]
    lab = np.array([labels[c] for c in cols])
    # 秩（沿样本轴）
    Rk = stats.rankdata(sub, axis=1)                 # genes x n
    # 结值校正因子（固定）
    # sum(t^3 - t) per gene
    sumt = np.zeros(sub.shape[0])
    for g in range(sub.shape[0]):
        _, cnt = np.unique(sub[g], return_counts=True)
        sumt[g] = np.sum(cnt**3 - cnt)
    corr = (n + 1) - sumt / (n * (n - 1))
    corr = np.maximum(corr, 1e-9)

    def count_sig(maskR):
        nR = int(maskR.sum()); nN = n - nR
        if nR < 2 or nN < 2: return 0
        U = Rk[:, maskR].sum(axis=1) - nR * (nR + 1) / 2
        mu = nR * nN / 2
        sd = np.sqrt(nR * nN / 12 * corr)
        z = (U - mu) / sd
        p = 2 * (1 - stats.norm.cdf(np.abs(z)))
        return int(np.sum(p < 0.05))

    obs = count_sig(lab == "R")

    if not perm_labels:
        return dict(n=n, observed=obs, null_mean=None, null_p95=None, perm_P=None)

    # 患者层置换
    pt_of = np.array([c.split("_")[0] for c in cols])
    pids = list(dict.fromkeys(pt_of.tolist()))
    plab = {p: labels[[c for c in cols if c.split("_")[0] == p][0]] for p in pids}
    nRp = sum(1 for p in pids if plab[p] == "R")
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, dtype=int)
    for i in range(n_perm):
        perm = rng.permutation(len(pids))
        Rset = set(np.array(pids)[perm][:nRp].tolist())
        maskR = np.isin(pt_of, list(Rset))
        null[i] = count_sig(maskR)
    return dict(n=n, observed=obs, null_mean=float(null.mean()),
                null_p95=float(np.percentile(null, 95)), perm_P=float(np.mean(null >= obs)),
                nR_obs=int((lab == "R").sum()), n_patients=len(pids), nR_patients=nRp)

res = {}
print("\n" + "=" * 100)
print("【框架对照】观测显著基因数 vs 患者聚类置换零分布（P<0.05，名义）")
print("  注：本表用同一管线内部自洽计算，用于框架间比较；绝对值与稿件冻结值口径不同，")
print("      仅框架间差值可解释（见 RULES：不引用自算绝对计数）。")
print("=" * 100)
for tag, cols in [("A  pre-only (49)", A), ("On on-treatment-only (56)", ON),
                  ("B  naive pooled (105)", B), ("C  patient-level (62)", C)]:
    r = framework_null(cols, n_perm=200)
    res[tag] = r
    print("\n  [%s]  样本数=%d  (R样本=%s, 患者=%s)" % (tag, r["n"], r.get("nR_obs"), r.get("n_patients")))
    print("      观测显著基因数        = %d" % r["observed"])
    if r["null_mean"] is not None:
        print("      置换零分布均值        = %.0f   (95分位 = %.0f)" % (r["null_mean"], r["null_p95"]))
        print("      置换 P(观测≥零分布)   = %.4f" % r["perm_P"])
        print("      → 超额显著数 = %d" % (r["observed"] - int(r["null_mean"])))

json.dump(res, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\nsaved:", OUT)
