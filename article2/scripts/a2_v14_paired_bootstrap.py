# -*- coding: utf-8 -*-
"""
A2 v14 决定性分析：为「合并是否稀释真信号」这个承重主张提供不确定性量化。

外部评审（2026-10-15）指出：+2,927 -> +2,092 的「稀释」全部来自零分布均值抬升
(1,089 -> 2,243)，而 B 的观测命中反而更多 (4,334 > 4,015)，且该差值没有任何
不确定性量化。本脚本回答两个问题：

【Part 1】固定患者集 + 固定基因集：42 例配对患者（9R/33N），
          pre-only(42) / on-only(42) / pooled(84) 三框架对比。
          ⇒ 只改变「是否加入第二个时点」，患者组成完全不变。
          这是评审要求的隔离实验（分离「多一个时点」与「错计独立样本数」）。

【Part 2】成对患者轨迹自举：对 42 例配对患者重抽样（重复抽中的患者作为不同簇），
          每次重算 on-only 与 pooled 两框架的观测计数 + 内层标签置换零均值，
          给出「两框架超额之差」的 95% CI。

【Part 3】把四框架主表（A / On-only / B / C）的置换次数从 500 提到 10,000，
          报告尾部命中次数，并采用 (b+1)/(B+1)。

口径与冻结管线严格一致：含结值 + 连续性校正的正态近似 Mann-Whitney。
（该向量化实现已复现 scipy 冻结计数，偏差 0.0%。）

输出：results/a2_v14_paired_bootstrap.json
"""
import csv, json, os, sys, time
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")

_HERE = os.path.dirname(os.path.abspath(__file__))
A2    = os.path.dirname(_HERE)
ROOT  = os.path.dirname(A2)
FPKM  = os.environ.get("A2_GSE91061_FPKM", os.path.join(
        ROOT, "output", "gse91061_validation",
        "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"))
MASKS = os.path.join(A2, "results", "framework_cohort_masks_v2.json")
OUT   = os.path.join(A2, "results", "a2_v14_paired_bootstrap.json")
if not os.path.exists(FPKM):
    raise SystemExit("Input FPKM matrix not found: %s" % FPKM)

QUICK = "--quick" in sys.argv
NPERM_MAIN = 200 if QUICK else 10000      # Part 3
NPERM_FIX  = 200 if QUICK else 10000      # Part 1
NBOOT      = 20 if QUICK else 600         # Part 2 outer
NINNER     = 20 if QUICK else 80          # Part 2 inner perms
SEED       = 20261015

# Guard: a smoke test must never overwrite the production artefact. Running with
# --quick writes to a separate _quick file; only a full run touches OUT.
if QUICK:
    OUT = os.path.join(A2, "results", "a2_v14_paired_bootstrap_quick.json")
    print("[--quick] smoke test: writing to %s (the full-run artefact is left untouched)" % OUT)

# ---------------- load ----------------
with open(FPKM, newline="", encoding="utf-8") as f:
    r = csv.reader(f); hdr = next(r)
    samples_all = [h.strip() for h in hdr[1:]]
    rows = []
    for row in r:
        if not row: continue
        rows.append([float(x) for x in row[1:]])
X = np.array(rows, float)
sid = {s: i for i, s in enumerate(samples_all)}
masks = json.load(open(MASKS, encoding="utf-8"))
lab = masks["labels"]
A_pre, B_pool, C_pat = masks["A_pre"], masks["B_naive_pooled"], masks["C_patient_pre_else_on"]
G = X.shape[0]
print("matrix: %d genes x %d samples" % (G, X.shape[1]))

# ---------------- frozen-caliber engine ----------------
def tie_sum(S):
    """sum(t^3 - t) per gene for a sorted (genes x n) matrix, vectorised."""
    g, n = S.shape
    eq = np.empty((g, n), dtype=bool)
    eq[:, 0] = False
    np.equal(S[:, 1:], S[:, :-1], out=eq[:, 1:])
    runidx = np.cumsum(~eq, axis=1) - 1
    nruns = runidx[:, -1] + 1
    off = (np.arange(g, dtype=np.int64) * (n + 1))[:, None]
    flat = (runidx + off).ravel()
    cnt = np.bincount(flat, minlength=g * (n + 1)).reshape(g, n + 1)
    valid = np.arange(n + 1)[None, :] < nruns[:, None]
    t = np.where(valid, cnt, 0).astype(np.float64)
    return (t ** 3 - t).sum(axis=1)

class Frame:
    """Precompute rank matrix + tie correction once; then each label assignment is O(genes)."""
    def __init__(self, cols, labv):
        self.cols = list(cols); self.labv = np.asarray(labv)
        sub = X[:, [sid[c] for c in self.cols]]
        self.n = sub.shape[1]
        self.Rk = stats.rankdata(sub, axis=1)
        st = tie_sum(np.sort(sub, axis=1))
        corr = (self.n + 1) - st / (self.n * (self.n - 1))
        self.corr = np.maximum(corr, 1e-9)
        self.pt = np.array([c.split("_")[0] for c in self.cols])

    def sig_count(self, maskR):
        nR = int(maskR.sum()); nN = self.n - nR
        if nR == 0 or nN == 0: return 0
        U = self.Rk[:, maskR].sum(axis=1) - nR * (nR + 1) / 2
        mu = nR * nN / 2.0
        sd = np.sqrt(nR * nN / 12.0 * self.corr)
        d = np.maximum(np.abs(U - mu) - 0.5, 0.0)
        return int(np.sum(2 * (1 - stats.norm.cdf(d / sd)) < 0.05))

    def patient_perm(self, nperm, rng, obs=None):
        """Permute labels across whole patients (fixed number of R patients)."""
        pids = list(dict.fromkeys(self.pt.tolist()))
        plab = {p: self.labv[self.pt == p][0] for p in pids}
        nRp = sum(1 for p in pids if plab[p] == "R")
        null = np.empty(nperm, dtype=np.int64)
        pids_arr = np.array(pids)
        for i in range(nperm):
            Rset = set(pids_arr[rng.permutation(len(pids))][:nRp].tolist())
            null[i] = self.sig_count(np.isin(self.pt, list(Rset)))
        out = dict(nperm=nperm, n=self.n, n_patients=len(pids), nR_patients=nRp,
                   null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
                   null_pct=round(100 * float(null.mean()) / G, 2))
        if obs is not None:
            b = int(np.sum(null >= obs))
            out["observed"] = int(obs)
            out["perm_P"] = (b + 1) / (nperm + 1)
            out["tail_hits"] = b
            out["excess"] = int(obs - null.mean())
        return out

# ---------------- Part 1: fixed 42-patient set ----------------
print("\n" + "=" * 78)
print("Part 1  固定患者集（42 例配对患者，9R/33N）+ 固定基因集（%d）" % G)
print("=" * 78)
byp = {}
for s in B_pool:
    pid = s.split("_")[0]
    byp.setdefault(pid, {})["on" if "_On_" in s else "pre"] = s
paired = sorted([p for p, v in byp.items() if "pre" in v and "on" in v])
p_pre = [byp[p]["pre"] for p in paired]
p_on  = [byp[p]["on"]  for p in paired]
p_pool = p_pre + p_on
labof = lambda p: lab[byp[p]["pre"]]
print("paired patients = %d  (R=%d, N=%d)" % (len(paired), sum(1 for p in paired if labof(p) == "R"),
                                              sum(1 for p in paired if labof(p) == "N")))

rng = np.random.default_rng(SEED)
fix = {}
for name, cols in [("P_pre42", p_pre), ("P_on42", p_on), ("P_pool84", p_pool)]:
    labv = [lab[c] for c in cols]
    fr = Frame(cols, labv)
    obs = fr.sig_count(np.array(labv) == "R")
    res = fr.patient_perm(NPERM_FIX, rng, obs=obs)
    fix[name] = res
    print("  %-9s n=%-4d 观测=%-5d 零均值=%-8.1f (%.2f%%)  置换P=%.5f (尾部命中 %d/%d)  超额=%+d"
          % (name, res["n"], obs, res["null_mean"], res["null_pct"], res["perm_P"],
             res["tail_hits"], res["nperm"], res["excess"]))

e_on  = fix["P_on42"]["excess"]; e_pool = fix["P_pool84"]["excess"]
print("\n  ⇒ 同一批 42 例患者，仅改变分析单位：")
print("     on-only(42样本) 超额 = %+d ; pooled(84样本) 超额 = %+d ; 差 = %+d" % (e_on, e_pool, e_pool - e_on))
print("     观测计数变化 %d -> %d (%+d) ; 零均值变化 %.1f -> %.1f (%+.1f)"
      % (fix["P_on42"]["observed"], fix["P_pool84"]["observed"],
         fix["P_pool84"]["observed"] - fix["P_on42"]["observed"],
         fix["P_on42"]["null_mean"], fix["P_pool84"]["null_mean"],
         fix["P_pool84"]["null_mean"] - fix["P_on42"]["null_mean"]))

# ---------------- Part 2: paired, response-stratified patient-trajectory bootstrap ----------------
def excess_for(cols, labv, cluster, rng):
    """observed count + inner label-permutation null mean for an arbitrary design."""
    fr = Frame(cols, labv)
    fr.pt = np.asarray(cluster)          # duplicate patient -> its own cluster
    obs = fr.sig_count(np.asarray(labv) == "R")
    r = fr.patient_perm(NINNER, rng, obs=obs)
    return obs, r["null_mean"], r["excess"]

def make_specs(pick, suffix):
    """pick = list of (patient, cluster_tag); returns (on_spec, pool_spec)."""
    on_s, po_s = [], []
    for p, tag in pick:
        on_s.append((byp[p]["on"], tag, labof(p)))
        po_s.append((byp[p]["pre"], tag, labof(p)))
        po_s.append((byp[p]["on"],  tag, labof(p)))
    return on_s, po_s

print("\n" + "=" * 78)
print("Part 2  成对、按应答分层的患者轨迹自举（B=%d，内层置换 %d 次/框架）" % (NBOOT, NINNER))
print("=" * 78)
R_pats = [p for p in paired if labof(p) == "R"]
N_pats = [p for p in paired if labof(p) == "N"]
print("  分层池: R=%d, N=%d" % (len(R_pats), len(N_pats)))
rngb = np.random.default_rng(SEED + 1)
diffs, exc_on, exc_pool, cnt_on, cnt_pool = [], [], [], [], []
t0 = time.time()
for b in range(NBOOT):
    pick = []
    for i in rngb.integers(0, len(R_pats), len(R_pats)): pick.append((R_pats[i], "R%d__b%d" % (i, b)))
    for i in rngb.integers(0, len(N_pats), len(N_pats)): pick.append((N_pats[i], "N%d__b%d" % (i, b)))
    on_s, po_s = make_specs(pick, b)
    sb = {}
    for nm, spec in [("on", on_s), ("pool", po_s)]:
        sb[nm] = excess_for([c for c, _, _ in spec], [l for _, _, l in spec],
                            [t for _, t, _ in spec], rngb)
    cnt_on.append(sb["on"][0]); cnt_pool.append(sb["pool"][0])
    exc_on.append(sb["on"][2]);  exc_pool.append(sb["pool"][2])
    diffs.append(sb["pool"][2] - sb["on"][2])
    if (b + 1) % max(1, NBOOT // 10) == 0:
        el = time.time() - t0
        print("   ... %4d/%d 已用 %.0fs  预计总 %.0fs" % (b + 1, NBOOT, el, el * NBOOT / (b + 1)))
diffs = np.array(diffs, float)
lo, hi = np.percentile(diffs, [2.5, 97.5])
boot = dict(NBOOT=NBOOT, NINNER=NINNER, stratified=True,
            excess_diff_mean=float(diffs.mean()), excess_diff_sd=float(diffs.std(ddof=1)),
            ci95=[float(lo), float(hi)],
            excludes_zero=bool(lo > 0 or hi < 0),
            exc_on_mean=float(np.mean(exc_on)), exc_pool_mean=float(np.mean(exc_pool)),
            cnt_on_mean=float(np.mean(cnt_on)), cnt_pool_mean=float(np.mean(cnt_pool)),
            frac_diff_negative=float(np.mean(diffs < 0)))
print("\n  超额之差 (pooled - on-only): 均值 %+.0f, SD %.0f, 95%%CI [%+.0f, %+.0f]"
      % (boot["excess_diff_mean"], boot["excess_diff_sd"], lo, hi))
print("  CI 是否排除 0: %s ; 自举中差值为负的比例 = %.3f" % (boot["excludes_zero"], boot["frac_diff_negative"]))
print("  自举均值: 观测计数 on=%.0f pool=%.0f ; 超额 on=%+.0f pool=%+.0f"
      % (boot["cnt_on_mean"], boot["cnt_pool_mean"], boot["exc_on_mean"], boot["exc_pool_mean"]))
print("  注意: 自举含重复抽样患者（样本向量完全重复），会使两框架计数同向偏高，"
      "故另做 Part 4 刀切法（无重复抽样）作为对照。")

# ---------------- Part 4: leave-one-patient-out jackknife (duplication-free) ----------------
print("\n" + "=" * 78)
print("Part 4  留一患者刀切法（42 次，无重复抽样）")
print("=" * 78)
rngj = np.random.default_rng(SEED + 3)
jk = []
for i, drop in enumerate(paired):
    keep = [p for p in paired if p != drop]
    pick = [(p, p) for p in keep]
    on_s, po_s = make_specs(pick, 0)
    a = excess_for([c for c, _, _ in on_s],  [l for _, _, l in on_s],  [t for _, t, _ in on_s],  rngj)
    b2 = excess_for([c for c, _, _ in po_s], [l for _, _, l in po_s], [t for _, t, _ in po_s], rngj)
    jk.append((b2[2] - a[2], a[2], b2[2], a[1], b2[1], a[0], b2[0]))
    if (i + 1) % 10 == 0:
        print("   ... %d/42" % (i + 1))
jk = np.array(jk, float)
th = jk[:, 0]; n = len(th)
se_jk = float(np.sqrt((n - 1) / n * np.sum((th - th.mean()) ** 2)))
lo_j, hi_j = th.mean() - 1.96 * se_jk, th.mean() + 1.96 * se_jk
jack = dict(n=n, theta_mean=float(th.mean()), se=se_jk, ci95=[float(lo_j), float(hi_j)],
            excludes_zero=bool(lo_j > 0 or hi_j < 0),
            excess_on_mean=float(jk[:, 1].mean()), excess_pool_mean=float(jk[:, 2].mean()),
            null_on_mean=float(jk[:, 3].mean()), null_pool_mean=float(jk[:, 4].mean()),
            cnt_on_mean=float(jk[:, 5].mean()), cnt_pool_mean=float(jk[:, 6].mean()))
print("\n  超额之差 均值 %+.0f, 刀切 SE %.0f, 95%%CI [%+.0f, %+.0f], 排除 0: %s"
      % (jack["theta_mean"], jack["se"], lo_j, hi_j, jack["excludes_zero"]))
print("  留一均值: 观测 on=%.0f pool=%.0f ; 零均值 on=%.0f pool=%.0f ; 超额 on=%+.0f pool=%+.0f"
      % (jack["cnt_on_mean"], jack["cnt_pool_mean"], jack["null_on_mean"],
         jack["null_pool_mean"], jack["excess_on_mean"], jack["excess_pool_mean"]))

# ---------------- Part 3: 4-framework main table at 10,000 permutations ----------------
print("\n" + "=" * 78)
print("Part 3  四框架主表，置换次数 %d（(b+1)/(B+1)）" % NPERM_MAIN)
print("=" * 78)
ON = [s for s in B_pool if "_On_" in s]
rngm = np.random.default_rng(SEED + 2)
main = {}
for name, cols in [("A_pre", A_pre), ("On_only", ON), ("B_naive", B_pool), ("C_patient", C_pat)]:
    labv = [lab[c] for c in cols]
    fr = Frame(cols, labv)
    obs = fr.sig_count(np.array(labv) == "R")
    res = fr.patient_perm(NPERM_MAIN, rngm, obs=obs)
    main[name] = res
    print("  %-10s n=%-4d 观测=%-5d 零均值=%-8.1f (%.2f%%)  置换P=%.5f  尾部命中=%d/%d  超额=%+d"
          % (name, res["n"], obs, res["null_mean"], res["null_pct"], res["perm_P"],
             res["tail_hits"], res["nperm"], res["excess"]))

json.dump(dict(seed=SEED, n_genes=G, part1_fixed42=fix, part2_bootstrap=boot,
               part3_main_10k=main, part4_jackknife=jack),
          open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\nsaved:", OUT)
