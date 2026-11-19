# -*- coding: utf-8 -*-
"""
A2 v16: the rank-scale mechanism behind the S1 inflation, and a check of the
approximation quoted in the text.

Claim under test (Methods/Results, S1): in the balanced two-time-point design with
no time shift, the null variance of the pooled statistic is inflated by a design
effect 1 + rho_S, where rho_S is the within-patient correlation of the rank scale,
not the product-moment ICC rho of the expression values. For bivariate normal pairs
with correlation rho, the rank (Spearman) correlation is

    rho_S = (6/pi) * arcsin(rho/2),

so the predicted variance ratio is 1 + (6/pi)*arcsin(rho/2). At rho = 0.5 this is
1.4826 against the simulated 1.4816.

This script does not run a new simulation grid. It checks:
  (a) the identity rho_S = (6/pi)*arcsin(rho/2) by Monte Carlo on 200,000 normal
      pairs per ICC level;
  (b) the predicted ratio against the frozen S2 mechanism output
      (results/a2_p2_s2_mechanism_20261118.json, tau = 0 rows), which measured the
      permutation variance ratio directly on the generator at 500 genes;
  (c) the observed within-patient rank correlation in a generator draw at each ICC,
      to show that the design effect tracks it rather than rho.

Conditions: balanced (one sample per time point per patient), continuous Gaussian
errors, no time effect, no association. The approximation is not carried to S2
(a shared shift) or to any real cohort.

Output: results/a2_p2_icc_rank_approx_20261120.json
"""
import json, os, sys
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
_HERE = os.path.dirname(os.path.abspath(__file__))
A2 = os.path.dirname(_HERE)
MECH = os.path.join(A2, "results", "a2_p2_s2_mechanism_20261118.json")
OUT = os.path.join(A2, "results", "a2_p2_icc_rank_approx_20261120.json")

SEED = 20261120
NPAIR = 200000
RHOS = (0.2, 0.5, 0.8)


def rho_s_closed(rho):
    return (6.0 / np.pi) * np.arcsin(rho / 2.0)


def main():
    rng = np.random.default_rng(SEED)
    mech = json.load(open(MECH, encoding="utf-8"))
    frozen = {(r["tau"], r["icc"]): r["s2"] for r in mech["rows"]}
    frozen_sd = {(r["tau"], r["icc"]): r["sd"] for r in mech["rows"]}

    rows = []
    for rho in RHOS:
        # (a) Monte Carlo rank correlation of paired normal draws
        z1 = rng.normal(size=NPAIR)
        z2 = rho * z1 + np.sqrt(1.0 - rho ** 2) * rng.normal(size=NPAIR)
        rho_s_mc = float(stats.spearmanr(z1, z2).statistic)
        rho_s_cf = float(rho_s_closed(rho))
        pred = 1.0 + rho_s_cf
        s2 = float(frozen[(0.0, rho)])

        # (b) implied Type I error of the pooled test if the null SD is sqrt(1+rho_S)
        implied = float(2 * (1 - stats.norm.cdf(1.96 / np.sqrt(pred))))

        rows.append(dict(icc=rho, rho_s_mc=rho_s_mc,
                         rho_s_closed=rho_s_cf,
                         predicted_variance_ratio=pred,
                         simulated_variance_ratio=s2,
                         simulated_sd=frozen_sd[(0.0, rho)],
                         rel_error=float(pred / s2 - 1.0),
                         naive_one_plus_icc=1.0 + rho,
                         rel_error_naive=float((1.0 + rho) / s2 - 1.0),
                         implied_type1_from_prediction=implied))

    # empirical within-patient rank correlation on a single generator draw, to show
    # that the design effect follows the rank scale rather than the raw ICC
    def gen_pairs(icc, P=80, G=4000):
        nR = int(round(P * 0.25))
        lab = np.zeros(P, dtype=bool); lab[rng.choice(P, nR, replace=False)] = True
        sb, se = np.sqrt(icc), np.sqrt(1.0 - icc)
        b = rng.normal(0, sb, (G, P))
        pre = b + rng.normal(0, se, (G, P))
        on = b + rng.normal(0, se, (G, P))
        return pre, on

    for row in rows:
        pre, on = gen_pairs(row["icc"])
        rs = [stats.spearmanr(pre[g], on[g]).statistic for g in range(pre.shape[0])]
        row["observed_within_patient_rank_corr"] = float(np.mean(rs))
        row["observed_within_patient_pearson"] = float(
            np.mean([np.corrcoef(pre[g], on[g])[0, 1] for g in range(pre.shape[0])]))

    out = dict(provenance=dict(
        frozen_mechanism=os.path.basename(MECH), n_pairs=NPAIR, seed=SEED,
        note=("no new grid cell; (a) is a Monte Carlo check of the normal-scores "
              "identity, (b) compares with the frozen tau=0 mechanism rows")),
        formula="variance_ratio = 1 + (6/pi)*arcsin(rho/2)",
        rows=rows)
    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    print("%-5s %-8s %-8s %-9s %-9s %-8s %-9s" %
          ("ICC", "rhoS_mc", "rhoS_cf", "1+rhoS", "s2(frozen)", "rel.err", "rel.err_1+ICC"))
    for r in rows:
        print("%-5.1f %-8.4f %-8.4f %-9.4f %-9.4f %+8.2f%% %+9.2f%%"
              % (r["icc"], r["rho_s_mc"], r["rho_s_closed"],
                 r["predicted_variance_ratio"], r["simulated_variance_ratio"],
                 100 * r["rel_error"], 100 * r["rel_error_naive"]))
    print()
    for r in rows:
        print("ICC %.1f: rank corr in generator %.4f (Pearson %.4f)"
              % (r["icc"], r["observed_within_patient_rank_corr"],
                 r["observed_within_patient_pearson"]))
    print("written:", OUT)


if __name__ == "__main__":
    main()
