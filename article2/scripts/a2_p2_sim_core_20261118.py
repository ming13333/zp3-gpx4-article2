# -*- coding: utf-8 -*-
"""
A2 P2 round 1: the 2x2 core simulation (S1-S4)
Protocol: reports/A2_P2_simulation_protocol_20261118.md (frozen)

DGP: X_git = mu_g + b_gi + a_g*R_i + tau_g*t + delta_g*R_i*t + eps_git

Data layout: columns are patient-major, so patient i has baseline at 2i and on-treatment at 2i+1.
  => baseline columns = X[:, 0::2], on-treatment columns = X[:, 1::2], Delta = X[:,1::2]-X[:,0::2]
  => patient-level labels are expanded to sample level with np.repeat(lab, 2)

Analyst view (protocol section 3): T1 uses the baseline columns only; T2 the on-treatment columns;
T3 the Delta only. The cross-time-point (pooled) procedure is the target-mismatch control.

Usage:
  python a2_p2_sim_core_20261118.py            # smoke test (G=200, n_outer=30)
  A2_FULL=1 python a2_p2_sim_core_20261118.py  # full run (G=500, n_outer=1000)
"""
import json, os, sys, time
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
_HERE = os.path.dirname(os.path.abspath(__file__))
A2 = os.path.dirname(_HERE)
OUTDIR = os.path.join(A2, "results")
FULL = os.environ.get("A2_FULL", "0") == "1"

CFG = dict(
    G          = int(os.environ.get("A2_G", "500" if FULL else "200")),
    n_pat      = int(os.environ.get("A2_NPAT", "80")),
    resp_prop  = float(os.environ.get("A2_RESP", "0.25")),
    icc        = float(os.environ.get("A2_ICC", "0.5")),
    a0         = 0.5,      # baseline between-group difference (log2)
    tau0       = 1.0,      # common time effect
    delta0     = 0.6,      # response-specific change
    sig_frac   = 0.10,     # fraction of signal genes
    n_outer    = int(os.environ.get("A2_OUTER", "1000" if FULL else "30")),
    n_perm     = int(os.environ.get("A2_PERM", "200" if FULL else "60")),
    alpha      = 0.05,
    seed       = 20261118,
)
SUF = os.environ.get("A2_SUFFIX", "")
TAG = ("full" if FULL else "smoke") + (("_"+SUF) if SUF else "")
OUT = os.path.join(OUTDIR, "a2_p2_sim_core_20261118_%s.json" % TAG)

# ---------------------------------------------------------------- DGP
def simulate(rng, cfg, scen):
    """Return X (G x 2P), lab (P,) bool responder indicator, truth dict"""
    G, P = cfg["G"], cfg["n_pat"]
    nR = int(round(P * cfg["resp_prop"]))
    lab = np.zeros(P, dtype=bool); lab[rng.choice(P, nR, replace=False)] = True
    k  = max(1, int(round(G * cfg["sig_frac"])))

    a = np.zeros(G); tau = np.zeros(G); d = np.zeros(G)
    if scen == "S2":
        tau[:] = cfg["tau0"]                       # the same common time effect in both groups
    elif scen == "S3":
        a[:k] = cfg["a0"]
    elif scen == "S4":
        d[:k] = cfg["delta0"]
    # S1: all zero

    mu = np.zeros(G)   # gene baseline fixed at 0: every test is within-gene, so mu cannot affect
                       # any statistic, and 0 is what makes the total variance 1 (otherwise
    sb = np.sqrt(cfg["icc"]); se = np.sqrt(1.0 - cfg["icc"])
    b  = rng.normal(0, sb, (G, P))
    e0 = rng.normal(0, se, (G, P)); e1 = rng.normal(0, se, (G, P))
    R  = lab.astype(float)

    pre = mu[:, None] + b + a[:, None] * R[None, :] + e0
    on  = mu[:, None] + b + a[:, None] * R[None, :] + tau[:, None] + d[:, None] * R[None, :] + e1
    X = np.empty((G, 2 * P)); X[:, 0::2] = pre; X[:, 1::2] = on
    truth = dict(a=a, tau=tau, delta=d, sig_idx=np.arange(k))
    return X, lab, truth

# ---------------------------------------------------------------- methods
def _tie_term(Xsub):
    n = Xsub.shape[1]
    Xs = np.sort(Xsub, axis=1)
    hastie = np.any(np.diff(Xs, axis=1) == 0, axis=1)
    tie = np.zeros(Xsub.shape[0])
    for i in np.where(hastie)[0]:
        _, c = np.unique(Xsub[i], return_counts=True)
        tie[i] = np.sum(c.astype(float) ** 3 - c) / (n * (n - 1.0))
    return tie

def mwu_p(Xsub, isR, tie=None):
    """Per-gene Mann-Whitney U with tie and continuity correction (normal approximation)"""
    n = Xsub.shape[1]; nR = int(isR.sum()); nN = n - nR
    if nR < 2 or nN < 2:
        return np.full(Xsub.shape[0], np.nan)
    Rk = stats.rankdata(Xsub, axis=1)
    if tie is None:
        tie = _tie_term(Xsub)
    U  = Rk[:, isR].sum(axis=1) - nR * (nR + 1.0) / 2.0
    mu = nR * nN / 2.0
    var = nR * nN / 12.0 * ((n + 1.0) - tie)
    dd = U - mu
    dd = np.where(dd > 0, dd - 0.5, np.where(dd < 0, dd + 0.5, dd))
    z = dd / np.sqrt(np.where(var > 0, var, np.nan))
    return 2.0 * (1.0 - stats.norm.cdf(np.abs(z)))

def delta_welch(D, isR):
    """Welch t on the paired Delta: under a balanced design this equals the interaction test of a random-intercept model"""
    dR, dN = D[:, isR], D[:, ~isR]
    nR, nN = dR.shape[1], dN.shape[1]
    mR, mN = dR.mean(1), dN.mean(1)
    vR, vN = dR.var(1, ddof=1), dN.var(1, ddof=1)
    se = np.sqrt(vR / nR + vN / nN)
    est = mR - mN
    with np.errstate(divide="ignore", invalid="ignore"):
        t = est / se
        df = (vR / nR + vN / nN) ** 2 / ((vR / nR) ** 2 / (nR - 1) + (vN / nN) ** 2 / (nN - 1))
        p = 2.0 * stats.t.sf(np.abs(t), df)
        tc = stats.t.ppf(0.975, df)
    return p, est, est - tc * se, est + tc * se

def pooled_cluster_perm(Xpool, lab, tie, nperm, rng, nconst):
    """Pool every sample, then permute patient labels (2 samples per patient -> np.repeat)"""
    Rk = stats.rankdata(Xpool, axis=1); n = Xpool.shape[1]
    def zstat(isR_s):
        nR = int(isR_s.sum()); nN = n - nR
        if nR < 2 or nN < 2:
            return np.full(Xpool.shape[0], np.nan)
        U  = Rk[:, isR_s].sum(axis=1) - nR * (nR + 1.0) / 2.0
        mu = nR * nN / 2.0
        var = nR * nN / 12.0 * nconst
        dd = U - mu
        dd = np.where(dd > 0, dd - 0.5, np.where(dd < 0, dd + 0.5, dd))
        return dd / np.sqrt(np.where(var > 0, var, np.nan))
    zobs = zstat(np.repeat(lab, 2))
    cnt = np.zeros(Xpool.shape[0])
    for _ in range(nperm):
        zp = zstat(np.repeat(rng.permutation(lab), 2))
        cnt += (np.abs(zp) >= np.abs(zobs) - 1e-12).astype(float)
    return (cnt + 1.0) / (nperm + 1.0)

# ---------------------------------------------------------------- generator checks
def generator_checks(cfg, rng, n_draw=20):
    """Averaging over independent draws reduces single-draw noise: with 200 genes the SE is about
    0.02, so a 0.02 tolerance cannot be judged from one draw. The variance check uses the
    baseline columns only, because a common time shift raises the pooled variance, which is
    not a defect."""

    rep = []
    for scen in ["S1", "S2", "S3", "S4"]:
        bd, cd, icc, vv, rf = [], [], [], [], []
        for _ in range(n_draw):
            X, lab, tr = simulate(rng, cfg, scen)
            pre, on = X[:, 0::2], X[:, 1::2]
            D = on - pre
            bd.append((pre[:, lab].mean(1) - pre[:, ~lab].mean(1)).mean())
            cd.append((D[:, lab].mean(1) - D[:, ~lab].mean(1)).mean())
            icc.append(np.mean([np.corrcoef(pre[i], on[i])[0, 1] for i in range(min(150, X.shape[0]))]))
            vv.append(pre.var()); rf.append(lab.mean())
        rep.append(dict(scen=scen, mean_baseline_diff=float(np.mean(bd)),
                        mean_change_diff=float(np.mean(cd)), emp_icc=float(np.mean(icc)),
                        var_baseline=float(np.mean(vv)), resp_frac=float(np.mean(rf)),
                        se_baseline_diff=float(np.std(bd) / np.sqrt(n_draw)),
                        se_change_diff=float(np.std(cd) / np.sqrt(n_draw))))
    checks, tol = {}, 0.02
    checks["S1_T1_zero"] = abs(rep[0]["mean_baseline_diff"]) < tol
    checks["S1_T3_zero"] = abs(rep[0]["mean_change_diff"]) < tol
    checks["S2_T1_zero"] = abs(rep[1]["mean_baseline_diff"]) < tol
    checks["S2_T3_zero"] = abs(rep[1]["mean_change_diff"]) < tol
    checks["S3_T3_zero"] = abs(rep[2]["mean_change_diff"]) < tol
    checks["S4_T1_zero"] = abs(rep[3]["mean_baseline_diff"]) < tol
    checks["ICC_within_0.05"]  = all(abs(r["emp_icc"] - cfg["icc"]) < 0.05 for r in rep)
    checks["resp_frac_ok"]     = all(abs(r["resp_frac"] - cfg["resp_prop"]) < 0.02 for r in rep)
    checks["var_base_within5"] = all(abs(r["var_baseline"] - 1.0) < 0.05 for r in rep)
    return rep, checks

# ---------------------------------------------------------------- main loop
def main():
    t0 = time.time(); cfg = CFG
    rng = np.random.default_rng(cfg["seed"])
    print("=" * 78)
    print("P2 round 1, 2x2 core simulation  [%s]  G=%d  P=%d  resp=%.2f  ICC=%.2f  n_outer=%d  n_perm=%d"
          % (TAG, cfg["G"], cfg["n_pat"], cfg["resp_prop"], cfg["icc"], cfg["n_outer"], cfg["n_perm"]))
    print("=" * 78)

    reps, checks = generator_checks(cfg, rng)
    print("\n[generator checks]")
    for r in reps:
        print("  %s  mean baseline diff=%+.4f (SE %.4f)  mean change diff=%+.4f (SE %.4f)  empirical ICC=%.3f  baseline var=%.3f  responder fraction=%.3f"
              % (r["scen"], r["mean_baseline_diff"], r["se_baseline_diff"],
                 r["mean_change_diff"], r["se_change_diff"],
                 r["emp_icc"], r["var_baseline"], r["resp_frac"]))
    print("  verdict:", "ALL PASS" if all(checks.values()) else "FAIL")
    for k, v in checks.items():
        if not v: print("    ✗", k)
    if not all(checks.values()):
        print("  generator checks failed, stopping.")
        json.dump(dict(cfg=cfg, checks=checks, reps=reps), open(OUT, "w", encoding="utf-8"),
                  ensure_ascii=False, indent=2)
        return

    # accumulators: scen -> method -> arrays (rejmat holds rep x gene booleans for rep-level MC intervals)
    METHODS = ["mw_baseline", "mw_on", "mw_pooled", "pooled_cluster", "delta_mwu", "delta_welch"]
    acc = {}
    for scen in ["S1", "S2", "S3", "S4"]:
        acc[scen] = {mth: dict(rej=np.zeros(cfg["G"]), est=np.zeros(cfg["G"]),
                               cover=np.zeros(cfg["G"]), nan=np.zeros(cfg["G"]),
                               rejmat=np.zeros((cfg["n_outer"], cfg["G"]), dtype=bool))
                     for mth in METHODS}
        acc[scen]["_truth"] = None

    print("\n[simulating]")
    for scen in ["S1", "S2", "S3", "S4"]:
        tS = time.time()
        for rep_i in range(cfg["n_outer"]):
            X, lab, tr = simulate(rng, cfg, scen)
            if acc[scen]["_truth"] is None: acc[scen]["_truth"] = tr
            pre, on, D = X[:, 0::2], X[:, 1::2], X[:, 1::2] - X[:, 0::2]

            # --- T1: baseline columns only
            p1 = mwu_p(pre, lab)
            # --- T2: on-treatment columns only
            p2 = mwu_p(on, lab)
            # --- cross-time-point (target-mismatch control)
            pP = mwu_p(X, np.repeat(lab, 2))
            tieX = _tie_term(X)
            pC = pooled_cluster_perm(X, lab, tieX, cfg["n_perm"], rng,
                                     (X.shape[1] + 1.0) - tieX)
            # --- T3: paired Delta only
            p3r = mwu_p(D, lab)
            p3m, est, lo, hi = delta_welch(D, lab)

            for mth, p in [("mw_baseline", p1), ("mw_on", p2), ("mw_pooled", pP),
                           ("pooled_cluster", pC), ("delta_mwu", p3r), ("delta_welch", p3m)]:
                a = acc[scen][mth]
                ok = ~np.isnan(p)
                rej_ok = (p[ok] < cfg["alpha"])
                a["rej"][ok] += rej_ok.astype(float)
                a["nan"][~ok] += 1.0
                a["rejmat"][rep_i, ok] = rej_ok
            a = acc[scen]["delta_welch"]
            truth_d = tr["delta"]
            a["est"] += est
            a["cover"] += ((lo <= truth_d) & (truth_d <= hi)).astype(float)
            if (rep_i + 1) % max(1, cfg["n_outer"] // 5) == 0:
                print("  %s  %d/%d  (%.0fs)" % (scen, rep_i + 1, cfg["n_outer"], time.time() - tS))

    # ---------------- summary
    def summarize(scen, mth, idx, N):
        """Monte Carlo intervals are estimated at the level of the outer replicate: form the
        proportion within each replicate, then take the mean and SE across replicates. This
        holds even when genes are correlated."""

        a = acc[scen][mth]
        per_rep = a["rejmat"][:, idx].mean(axis=1)     # length N
        pbar = float(per_rep.mean())
        se = float(per_rep.std(ddof=1) / np.sqrt(N)) if N > 1 else 0.0
        n = len(idx)
        return dict(rate=pbar, ci95=[pbar - 1.96 * se, pbar + 1.96 * se],
                    se_rep=se, n_genes=int(n),
                    fail_rate=float(a["nan"][idx].sum() / (n * N)))

    N = cfg["n_outer"]; G = cfg["G"]; k = max(1, int(round(G * cfg["sig_frac"])))
    sig = np.arange(k); nul = np.arange(k, G)

    # true state of every (scenario, target)
    PLAN = {
        "S1": dict(T1=("null", None), T2=("null", None), T3=("null", None)),
        "S2": dict(T1=("null", None), T2=("null", None), T3=("null", None)),
        "S3": dict(T1=("alt", "sig"), T2=("alt", "sig"), T3=("null", None)),
        "S4": dict(T1=("null", None), T2=("alt", "sig"), T3=("alt", "sig")),
    }
    MT = {"T1": ["mw_baseline", "mw_pooled", "pooled_cluster"],
          "T2": ["mw_on", "mw_pooled", "pooled_cluster"],
          "T3": ["delta_mwu", "delta_welch", "mw_pooled", "pooled_cluster"]}
    CORRECT = {"T1": ["mw_baseline"], "T2": ["mw_on"], "T3": ["delta_mwu", "delta_welch"]}

    summary = []
    print("\n" + "=" * 78); print("Summary (Type I error / power, alpha=0.05)"); print("=" * 78)
    print("%-4s %-4s %-16s %-8s %-22s %s" % ("scen", "tgt", "method", "kind", "rate [95% MC]", "n_genes"))
    for scen in ["S1", "S2", "S3", "S4"]:
        for tgt in ["T1", "T2", "T3"]:
            kind, which = PLAN[scen][tgt]
            for mth in MT[tgt]:
                if kind == "null":
                    s = summarize(scen, mth, np.arange(G), N); lab_kind = "TypeI"
                else:
                    s = summarize(scen, mth, sig, N); lab_kind = "power"
                tagm = mth + ("" if mth in CORRECT[tgt] else " *")
                print("%-4s %-4s %-16s %-8s %.4f [%.4f, %.4f]      %d"
                      % (scen, tgt, tagm, lab_kind, s["rate"], s["ci95"][0], s["ci95"][1], s["n_genes"]))
                summary.append(dict(scenario=scen, target=tgt, method=mth, kind=lab_kind,
                                    correct=(mth in CORRECT[tgt]), **s))
        # bias and coverage of delta_welch
        a = acc[scen]["delta_welch"]; truth_d = acc[scen]["_truth"]["delta"]
        bias_sig = float((a["est"][sig] / N - truth_d[sig]).mean())
        bias_nul = float((a["est"][nul] / N - truth_d[nul]).mean())
        cov = float(a["cover"].mean() / N)
        print("     [delta_welch] bias(signal)=%+.4f  bias(null)=%+.4f  CI coverage=%.4f" % (bias_sig, bias_nul, cov))
        summary.append(dict(scenario=scen, target="T3", method="delta_welch", kind="diag",
                            bias_signal=bias_sig, bias_null=bias_nul, coverage=cov))

    out = dict(cfg=cfg, tag=TAG, checks=checks, generator_reps=reps, summary=summary,
               elapsed_s=round(time.time() - t0, 1))
    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\nwritten: %s  (%.0fs)" % (OUT, time.time() - t0))

if __name__ == "__main__":
    main()
