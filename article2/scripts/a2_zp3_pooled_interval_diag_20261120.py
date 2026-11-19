# -*- coding: utf-8 -*-
"""
A2 v16: why the pooled bootstrap interval and the patient-cluster permutation P
need not agree, and what each one is.

Motivation (external review, 2026-11-20, item 1): in the 42 paired patients the
pooled 95% percentile interval excludes zero (-0.640 to -0.067) while the
patient-cluster permutation returns P = 0.069. A reader will read that as an
arithmetic error. It is not: the two quantify different reference distributions,
and this script measures both on the same 84 samples.

Reported here, all on the same 84 samples / 42 patients:
  (a) the point estimate and the percentile bootstrap interval, replayed with the
      same random stream and in the same call order as the frozen script, so that
      the interval reproduces exactly;
  (b) the bootstrap distribution's SD and its mass at r >= 0, i.e. what a test
      inverted from the same resampling would have returned;
  (c) the permutation null distribution of r (10,000 draws, same stream), its SD,
      its 2.5/97.5 percentiles and the tail count behind P = 0.069;
  (d) two analytic reference SDs for r: the independent-sample approximation
      sqrt((n+1)/(3*n1*n2)) and the patient-level approximation
      sqrt((n_pat+1)/(3*nR_pat*nN_pat));
  (e) the implied variance inflation of the permutation null against the
      independent-sample approximation.

Design effect aside: the inflation is what a calibrated test should show, because
the pooled statistic is computed on 84 samples but the labels are assigned to 42
patients.

Output: results/a2_zp3_pooled_interval_diag_20261120.json
"""
import importlib.util, json, os, sys, math
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")

_HERE = os.path.dirname(os.path.abspath(__file__))
A2 = os.path.dirname(_HERE)
ROOT = os.path.dirname(A2)
FROZEN = os.path.join(_HERE, "a2_zp3_fixed42_20261119.py")
FPKM = os.path.join(ROOT, "output", "gse91061_validation",
                    "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv")
MASKS = os.path.join(A2, "results", "framework_cohort_masks_v2.json")
FROZEN_JSON = os.path.join(A2, "results", "a2_zp3_fixed42_20261119.json")
OUT = os.path.join(A2, "results", "a2_zp3_pooled_interval_diag_20261120.json")

NPERM = 10000
NBOOT = 10000
SEED = 20261119


def load_module(path):
    """Import the frozen script so that every statistic is computed by the same
    code that produced the published values (no re-implementation)."""
    spec = importlib.util.spec_from_file_location("frozen_fixed42", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    fr = load_module(FROZEN)
    frozen = json.load(open(FROZEN_JSON, encoding="utf-8"))
    rng = np.random.default_rng(SEED)        # permutation, same seed and order
    rng_b = np.random.default_rng(SEED + 1)  # bootstrap, same seed and order

    # ---- rebuild the 42 paired patients exactly as the frozen script does
    with open(FPKM, newline="", encoding="utf-8") as f:
        import csv
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
    pre_by_pt, on_by_pt = {}, {}
    for s in masks["B_naive_pooled"]:
        parts = s.split("_")
        pt, tp = parts[0], parts[1]
        (pre_by_pt if tp == "Pre" else on_by_pt)[pt] = s
    paired = sorted(set(pre_by_pt) & set(on_by_pt))
    pre42 = [pre_by_pt[p] for p in paired]
    on42 = [on_by_pt[p] for p in paired]

    # ---- (a) replay the bootstrap stream in the frozen call order
    # pre42 -> on42 -> delta42 -> pooled84; the first three are replays whose only
    # purpose is to advance rng_b exactly as the frozen script did.
    def boot_replay(x, y):
        return fr.boot_r(x, y, NBOOT, rng_b)

    xp = np.array([expr[s] for s in pre42 if lab[s] == "R"])
    yp = np.array([expr[s] for s in pre42 if lab[s] == "N"])
    xo = np.array([expr[s] for s in on42 if lab[s] == "R"])
    yo = np.array([expr[s] for s in on42 if lab[s] == "N"])
    dR = np.array([expr[on_by_pt[p]] - expr[pre_by_pt[p]] for p in paired
                   if lab[on_by_pt[p]] == "R"])
    dN = np.array([expr[on_by_pt[p]] - expr[pre_by_pt[p]] for p in paired
                   if lab[on_by_pt[p]] == "N"])
    replays = {"pre42": boot_replay(xp, yp), "on42": boot_replay(xo, yo),
               "delta42": boot_replay(dR, dN)}

    # ---- the 84 pooled samples
    samp84 = list(pre42) + list(on42)
    rmask = np.array([lab[s] == "R" for s in samp84])
    pid = np.array([s.split("_")[0] for s in samp84])
    v = np.array([expr[s] for s in samp84])
    Rk = stats.rankdata(v)
    n1 = int(rmask.sum()); n2 = v.size - n1
    U = Rk[rmask].sum() - n1 * (n1 + 1.0) / 2.0
    r_obs = fr.r_of(U, n1, n2)

    pids = np.array(paired)
    nRp = len({p for p in pids if lab[pre_by_pt[p]] == "R"})
    nNp = len(pids) - nRp

    # stratified patient-cluster percentile bootstrap (patients resampled within
    # group, every sample a resampled patient contributes is kept, multiplicity
    # preserved because selection is by index)
    bs = np.empty(NBOOT)
    rs = [p for p in pids if lab[pre_by_pt[p]] == "R"]
    ns = [p for p in pids if lab[pre_by_pt[p]] == "N"]
    for i in range(NBOOT):
        pr = rng_b.choice(rs, nRp, replace=True)
        pn = rng_b.choice(ns, nNp, replace=True)
        sel = np.concatenate([np.where(pid == p)[0] for p in list(pr) + list(pn)])
        mm = rmask[sel]; vv = v[sel]
        rk = stats.rankdata(vv)
        a = int(mm.sum()); c = vv.size - a
        u = rk[mm].sum() - a * (a + 1.0) / 2.0
        bs[i] = fr.r_of(u, a, c)

    lo, hi = np.percentile(bs, [2.5, 97.5])
    boot = dict(CI=[float(lo), float(hi)], sd=float(bs.std(ddof=1)),
                median=float(np.median(bs)),
                mass_r_ge_0=float((bs >= 0).mean()),
                # a one-sided test inverted from the same resampling
                inv_P_two_sided=float(2 * min((bs >= 0).mean(), (bs <= 0).mean())),
                nboot=NBOOT)

    # ---- (c) permutation null of r, same stream
    # the frozen script draws the pre42 and on42 permutations before the pooled
    # one; replay them so that the pooled draw sequence is identical
    for x_, y_ in ((xp, yp), (xo, yo)):
        fr.perm_p(x_, y_, NPERM, rng)

    null = np.empty(NPERM)
    for k in range(NPERM):
        chosen = rng.permutation(pids)[:nRp]
        m = np.isin(pid, chosen)
        a = int(m.sum()); c = v.size - a
        u = Rk[m].sum() - a * (a + 1.0) / 2.0
        null[k] = fr.r_of(u, a, c)
    b_hits = int((np.abs(null) >= abs(r_obs)).sum())
    perm = dict(perm_P=float((b_hits + 1.0) / (NPERM + 1.0)), tail_hits=b_hits,
                nperm=NPERM, sd=float(null.std(ddof=1)),
                mean=float(null.mean()),
                q=[float(np.percentile(null, 2.5)), float(np.percentile(null, 97.5))],
                nR_patients=int(nRp), nN_patients=int(nNp))

    # ---- (d) analytic reference SDs
    se_iid = fr.se_of(n1, n2)                      # 84 samples as independent
    se_pt = fr.se_of(nRp, nNp)                     # 42 patients as the unit
    ref = dict(se_independent_samples=float(se_iid),
               se_patient_level=float(se_pt),
               var_inflation_vs_independent=float((perm["sd"] / se_iid) ** 2),
               var_inflation_vs_patient=float((perm["sd"] / se_pt) ** 2),
               boot_sd_over_se_independent=float(boot["sd"] / se_iid))

    # ---- reproduction guard against the frozen output
    fz = frozen["pooled84"]
    checks = {
        "r_match": abs(r_obs - fz["r"]) < 1e-12,
        "CI_match": abs(lo - fz["boot_CI"][0]) < 1e-12 and abs(hi - fz["boot_CI"][1]) < 1e-12,
        "perm_P_match": abs(perm["perm_P"] - fz["perm_P"]) < 1e-12,
        "tail_hits_match": b_hits == fz["tail_hits"],
        "replay_pre42_CI": abs(replays["pre42"][0] - frozen["pre42"]["boot_CI"][0]) < 1e-12,
        "replay_on42_CI": abs(replays["on42"][0] - frozen["on42"]["boot_CI"][0]) < 1e-12,
        "replay_delta42_CI": abs(replays["delta42"][0] - frozen["delta42"]["boot_CI"][0]) < 1e-12,
    }
    checks = {k: bool(v) for k, v in checks.items()}
    assert all(checks.values()), checks

    out = dict(provenance=dict(
        frozen_script=os.path.basename(FROZEN), seed=SEED, nperm=NPERM, nboot=NBOOT,
        note=("percentile bootstrap, patients resampled within response group, both "
              "samples per patient retained, multiplicity preserved"),
        framing=("the bootstrap interval describes the precision of the estimate; the "
                 "permutation P tests exchangeability of response labels across whole "
                 "patient records; they are not inverses of each other")),
        r=float(r_obs), n_samples=int(v.size), n_patients=int(len(pids)),
        n1=n1, n2=n2, bootstrap=boot, permutation=perm, reference=ref,
        reproduction=checks)
    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    print("pooled r = %+.4f  (84 samples / 42 patients; %d R vs %d N samples, %d vs %d patients)"
          % (r_obs, n1, n2, nRp, nNp))
    print("percentile bootstrap (B=%d): %+.4f to %+.4f   SD=%.4f   P(r>=0)=%.4f"
          % (NBOOT, lo, hi, boot["sd"], boot["mass_r_ge_0"]))
    print("permutation (B=%d): P=%.4f (b=%d)   null SD=%.4f   null 2.5-97.5%%: %+.3f to %+.3f"
          % (NPERM, perm["perm_P"], b_hits, perm["sd"], perm["q"][0], perm["q"][1]))
    print("reference SDs: independent samples %.4f ; patient level %.4f" % (se_iid, se_pt))
    print("inflation of the permutation null: x%.2f vs independent, x%.2f vs patient level"
          % (ref["var_inflation_vs_independent"], ref["var_inflation_vs_patient"]))
    print("reproduction:", "ALL PASS" if all(checks.values()) else checks)
    print("written:", OUT)


if __name__ == "__main__":
    main()
