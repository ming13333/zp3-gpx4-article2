# -*- coding: utf-8 -*-
"""
A2 v15: recompute ZP3 on one analysis population (the 42 strictly paired patients) for all three targets.

External review (2026-11-19, point 3): the text treated the baseline (49 samples), on-treatment
(56 samples) and change (42 paired patients) comparisons as if they came from one population, and
inferred the change effect from the two marginal r values. This script recomputes all three on the same population, as the parallel column of Table 2A.

Conventions match the frozen pipeline (scripts/a2_cluster_permutation.py and
scripts/a2_v14_paired_bootstrap.py:
  - Mann-Whitney U with tie and continuity correction (normal approximation)
  - r = 2U/(n1*n2) - 1 (rank-biserial; a linear transform of P(R>N) + 0.5*P(tie))
  - SE(r) = sqrt((n1+n2+1) / (3*n1*n2)); 95% CI = r ± 1.96 SE
  - patient-level permutation: labels reassigned by patient, P = (b+1)/(B+1), two-sided on abs(r)

Output: results/a2_zp3_fixed42_20261119.json
"""
import csv, json, math, os, sys
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")

_HERE = os.path.dirname(os.path.abspath(__file__))
A2 = os.path.dirname(_HERE)
ROOT = os.path.dirname(A2)
FPKM = os.environ.get("A2_GSE91061_FPKM", os.path.join(
    ROOT, "output", "gse91061_validation",
    "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"))
MASKS = os.path.join(A2, "results", "framework_cohort_masks_v2.json")
OUT = os.path.join(A2, "results", "a2_zp3_fixed42_20261119.json")

NPERM = 10000
NBOOT = 10000
SEED = 20261119


def r_of(U, n1, n2):
    return 2.0 * U / (n1 * n2) - 1.0


def se_of(n1, n2):
    return math.sqrt((n1 + n2 + 1.0) / (3.0 * n1 * n2))


def mwu(x, y):
    """Two independent samples; returns U (for x), two-sided P (normal approximation with tie and continuity correction) and r."""
    res = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic",
                             use_continuity=True)
    U = float(res.statistic)
    n1, n2 = len(x), len(y)
    return U, float(res.pvalue), r_of(U, n1, n2), n1, n2


def perm_p(x, y, nperm, rng):
    """Patient-level label permutation (here equal to sample level, one sample per patient), two-sided."""
    obs = r_of(stats.mannwhitneyu(x, y, alternative="two-sided",
                                  method="asymptotic").statistic, len(x), len(y))
    pool = np.concatenate([x, y])
    n1 = len(x)
    b = 0
    for _ in range(nperm):
        idx = rng.permutation(pool.size)
        xx, yy = pool[idx[:n1]], pool[idx[n1:]]
        U = stats.mannwhitneyu(xx, yy, alternative="two-sided",
                               method="asymptotic").statistic
        if abs(r_of(U, n1, pool.size - n1)) >= abs(obs) - 1e-12:
            b += 1
    return obs, (b + 1.0) / (nperm + 1.0), int(b)


def boot_r(x, y, nboot, rng):
    """Patient-level stratified bootstrap: resample patients within group, duplicates kept as separate draws. Returns the 2.5/97.5 percentiles of r."""
    bs = np.empty(nboot)
    for i in range(nboot):
        bx = rng.choice(x, x.size, replace=True)
        by = rng.choice(y, y.size, replace=True)
        u = stats.mannwhitneyu(bx, by, alternative="two-sided", method="asymptotic").statistic
        bs[i] = r_of(u, bx.size, by.size)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return float(lo), float(hi)


def main():
    # two independent streams: permutation and bootstrap are each reproducible, and editing one does not move the other
    rng = np.random.default_rng(SEED)          # permutation
    rng_b = np.random.default_rng(SEED + 1)    # bootstrap

    # ---- load expression
    with open(FPKM, newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        hdr = next(r)
        samples = [h.strip() for h in hdr[1:]]
        row = None
        for rec in r:
            if rec and rec[0].strip() == "7784":
                row = rec[1:]
                break
        if row is None:
            raise SystemExit("ZP3 row (Entrez 7784) not found")
    expr = {s: float(v) for s, v in zip(samples, row)}

    masks = json.load(open(MASKS, encoding="utf-8"))
    lab = masks["labels"]

    # ---- build the strictly paired 42
    # the paired 42 are drawn from the 105 annotated samples (B_naive_pooled),
    # not from every column of the matrix
    pre_by_pt, on_by_pt = {}, {}
    for s in masks["B_naive_pooled"]:
        parts = s.split("_")
        pt, tp = parts[0], parts[1]
        (pre_by_pt if tp == "Pre" else on_by_pt)[pt] = s
    paired = sorted(set(pre_by_pt) & set(on_by_pt))
    pre42 = [pre_by_pt[p] for p in paired]
    on42 = [on_by_pt[p] for p in paired]

    def xy(sample_list, subset=None):
        ss = [s for s in sample_list if subset is None or s in subset]
        x = np.array([expr[s] for s in ss if lab[s] == "R"])
        y = np.array([expr[s] for s in ss if lab[s] == "N"])
        return x, y, ss

    out = dict(provenance=dict(
        fpkm=os.path.basename(FPKM), entrez="7784", nperm=NPERM,
        nboot=NBOOT, seed=SEED,
        note="fixed 42-patient subset; one sample per time point per patient"))

    # ---- reproduction checks against frozen outputs
    a_pre = masks["A_pre"]
    o56 = sorted({s for s in samples if s.split("_")[1] == "On"} & set(lab))
    checks = {}
    for name, sl, target_r, target_p in [
        ("pre49", a_pre, -0.19487179487179485, 0.3522362191881916),
        ("on56", o56, -0.47048300536672627, 0.011014075190481976),
    ]:
        x, y, _ = xy(sl)
        U, P, rv, n1, n2 = mwu(x, y)
        checks[name] = dict(n1=n1, n2=n2, r=rv, P=P,
                            dr=abs(rv - target_r), dP=abs(P - target_p),
                            pass_=bool(abs(rv - target_r) < 1e-9 and abs(P - target_p) < 1e-9))
    out["reproduction"] = checks

    # ---- three targets on the same 42 patients
    x, y, _ = xy(pre42); U, P, rv, n1, n2 = mwu(x, y)
    _, Pp, b = perm_p(x, y, NPERM, rng)
    blo, bhi = boot_r(x, y, NBOOT, rng_b)
    out["pre42"] = dict(n1=n1, n2=n2, U=U, r=rv, P=P, SE=se_of(n1, n2),
                        CI=[rv - 1.96 * se_of(n1, n2), rv + 1.96 * se_of(n1, n2)],
                        boot_CI=[blo, bhi], nboot=NBOOT,
                        perm_P=Pp, tail_hits=b, nperm=NPERM,
                        median_R=float(np.median(x)), median_N=float(np.median(y)))

    x, y, _ = xy(on42); U, P, rv, n1, n2 = mwu(x, y)
    _, Pp, b = perm_p(x, y, NPERM, rng)
    blo, bhi = boot_r(x, y, NBOOT, rng_b)
    out["on42"] = dict(n1=n1, n2=n2, U=U, r=rv, P=P, SE=se_of(n1, n2),
                       CI=[rv - 1.96 * se_of(n1, n2), rv + 1.96 * se_of(n1, n2)],
                       boot_CI=[blo, bhi], nboot=NBOOT,
                       perm_P=Pp, tail_hits=b, nperm=NPERM,
                       median_R=float(np.median(x)), median_N=float(np.median(y)))

    # within-patient change
    dR = np.array([expr[on_by_pt[p]] - expr[pre_by_pt[p]] for p in paired if lab[on_by_pt[p]] == "R"])
    dN = np.array([expr[on_by_pt[p]] - expr[pre_by_pt[p]] for p in paired if lab[on_by_pt[p]] == "N"])
    U, P, rv, n1, n2 = mwu(dR, dN)
    se = se_of(n1, n2)

    # patient-level stratified bootstrap of r (resample patients within group)
    lo, hi = boot_r(dR, dN, NBOOT, rng_b)

    out["delta42"] = dict(n1=n1, n2=n2, U=U, r=rv, P=P, SE=se,
                          CI=[rv - 1.96 * se, rv + 1.96 * se],
                          boot_CI=[float(lo), float(hi)], nboot=NBOOT,
                          median_R=float(np.median(dR)), median_N=float(np.median(dN)),
                          n_paired=len(paired))

    # ---- pooled 84 samples in the same 42 patients, and its four-cell partition
    # (each patient contributes one pre and one on sample, so the four cells are
    #  equal-sized: 9*33 = 297 pairs each, 1,188 in total)
    samp84 = [s for s in pre42] + [s for s in on42]
    rmask84 = np.array([lab[s] == "R" for s in samp84])
    pid84 = np.array([s.split("_")[0] for s in samp84])
    tp84 = np.array(["Pre" if "_Pre_" in s else "On" for s in samp84])
    v84 = np.array([expr[s] for s in samp84])
    R84 = stats.rankdata(v84)
    nR84 = int(rmask84.sum()); nN84 = v84.size - nR84
    U84 = R84[rmask84].sum() - nR84 * (nR84 + 1.0) / 2.0
    r84 = r_of(U84, nR84, nN84)

    # patient-cluster permutation of the same statistic (patients carry 2 samples)
    pids84 = np.array(paired)
    nRp = len({p for p in pids84 if lab[pre_by_pt[p]] == "R"})
    b = 0
    for _ in range(NPERM):
        chosen = rng.permutation(pids84)[:nRp]
        m = np.isin(pid84, chosen)
        a = int(m.sum()); c = v84.size - a
        u = R84[m].sum() - a * (a + 1.0) / 2.0
        if abs(r_of(u, a, c)) >= abs(r84):
            b += 1
    pooled_perm_P = (b + 1.0) / (NPERM + 1.0)

    # stratified patient-level bootstrap (resample patients, keep both samples)
    bs = np.empty(NBOOT)
    for i in range(NBOOT):
        pr = rng_b.choice([p for p in pids84 if lab[pre_by_pt[p]] == "R"], nRp, replace=True)
        pn = rng_b.choice([p for p in pids84 if lab[pre_by_pt[p]] == "N"], len(pids84) - nRp, replace=True)
        take = set(pr.tolist()) | set(pn.tolist())
        # a resampled patient may appear twice; keep multiplicity by index
        sel = np.concatenate([np.where(pid84 == p)[0] for p in list(pr) + list(pn)])
        mm = rmask84[sel]
        vv = v84[sel]
        rk = stats.rankdata(vv)
        a = int(mm.sum()); c = vv.size - a
        u = rk[mm].sum() - a * (a + 1.0) / 2.0
        bs[i] = r_of(u, a, c)
    lo84, hi84 = np.percentile(bs, [2.5, 97.5])
    del take

    nominal_P84 = float(stats.mannwhitneyu(v84[rmask84], v84[~rmask84],
                                           alternative="two-sided",
                                           method="asymptotic",
                                           use_continuity=True).pvalue)
    out["pooled84"] = dict(n_samples=int(v84.size), n_patients=len(paired),
                           nR_patients=nRp, n1=nR84, n2=nN84, U=float(U84),
                           r=float(r84), nominal_P=nominal_P84,
                           perm_P=float(pooled_perm_P),
                           tail_hits=int(b), nperm=NPERM,
                           boot_CI=[float(lo84), float(hi84)], nboot=NBOOT)

    # four cells of the balanced partition
    cells = []
    for rt in ("Pre", "On"):
        for nt in ("Pre", "On"):
            xR = np.array([expr[s] for s in samp84 if lab[s] == "R" and tp84[list(samp84).index(s)] == rt])
            xN = np.array([expr[s] for s in samp84 if lab[s] == "N" and tp84[list(samp84).index(s)] == nt])
            u = stats.mannwhitneyu(xR, xN, alternative="two-sided",
                                   method="asymptotic").statistic
            npr = xR.size * xN.size
            cells.append(dict(responder=rt, non_responder=nt, pairs=int(npr),
                              win_rate=float(u / npr), r=float(r_of(u, xR.size, xN.size))))
    mean_cell_r = float(np.mean([c["r"] for c in cells]))
    out["cells84"] = dict(cells=cells, total_pairs=int(sum(c["pairs"] for c in cells)),
                          mean_cell_r=mean_cell_r, pooled_r=float(r84),
                          identity_gap=abs(mean_cell_r - r84))

    # ---- cross-check against frozen paired_delta
    fr = json.load(open(os.path.join(A2, "results", "a2_mechanism_decomp_20261014.json"),
                        encoding="utf-8"))["paired_delta"]
    out["reproduction"]["delta42"] = dict(
        r=rv, P=P, dr=abs(rv - fr["r"]), dP=abs(P - fr["P"]),
        pass_=bool(abs(rv - fr["r"]) < 1e-9 and abs(P - fr["P"]) < 1e-9))

    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    print("== reproduction ==")
    for k, v in out["reproduction"].items():
        print("  %-8s r=%+.4f  P=%.4g  dr=%.1e  dP=%.1e  %s"
              % (k, v["r"], v["P"], v["dr"], v["dP"], "OK" if v["pass_"] else "FAIL"))
    print("\n== same population (42 paired patients) ==")
    for k in ("pre42", "on42", "delta42"):
        v = out[k]
        print("  %-8s n=%d/%d  r=%+.4f (95%% CI %+.2f to %+.2f)  P=%.4f%s"
              % (k, v["n1"], v["n2"], v["r"], v["CI"][0], v["CI"][1], v["P"],
                 ("  paCI %+.2f to %+.2f  perm P=%.4f"
                  % (v["boot_CI"][0], v["boot_CI"][1], v["perm_P"]))
                 if "perm_P" in v else ("  paCI %+.2f to %+.2f" % (v["boot_CI"][0], v["boot_CI"][1]))))
    pv = out["pooled84"]
    print("  %-8s n=%d samples/%d patients  r=%+.4f  nominal P=%.4f  cluster P=%.4f  paCI %+.2f to %+.2f"
          % ("pooled84", pv["n_samples"], pv["n_patients"], pv["r"], pv["nominal_P"],
             pv["perm_P"],
             pv["boot_CI"][0], pv["boot_CI"][1]))
    print("\n== four cells in the balanced 42-patient set ==")
    for c in out["cells84"]["cells"]:
        print("  %s R x %s N  pairs=%d  win=%.3f  r=%+.4f"
              % (c["responder"], c["non_responder"], c["pairs"], c["win_rate"], c["r"]))
    print("  mean of four cell r = %+.4f ; pooled r = %+.4f ; gap=%.2e"
          % (out["cells84"]["mean_cell_r"], out["cells84"]["pooled_r"],
             out["cells84"]["identity_gap"]))
    print("\nwritten:", OUT)


if __name__ == "__main__":
    main()
