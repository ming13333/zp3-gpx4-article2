# -*- coding: utf-8 -*-
"""
Four-way decomposition of the pooled Mann-Whitney U statistic (GSE91061, ZP3).

Purpose
-------
The pooled framework (B) compares 105 samples as if independent. Those 105
samples split into four cells by (response label R/N) x (time point Pre/On).
The Mann-Whitney U statistic is a sum over all R-vs-N sample pairs, so it can
be partitioned exactly into four contributions:

    U = U(Pre_R vs Pre_N) + U(Pre_R vs On_N)
      + U(On_R  vs Pre_N) + U(On_R  vs On_N)

Two of the four cells are same-time-point comparisons (the comparisons a
baseline-predictive question would want); the other two are cross-time-point
comparisons, which compare a pre-treatment responder with an on-treatment
non-responder (or vice versa). The decomposition reports how much of the
pooled statistic each cell contributes.

This is a descriptive decomposition of the observed statistic. It is not a
test of time-point x response interaction.

Outputs: results/a2_pooled_fourway_decomp_20261119.json
"""
import csv, json, os, sys
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")

A2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.dirname(A2)
FPKM = os.environ.get("A2_GSE91061_FPKM", os.path.join(
    PROJ, "output", "gse91061_validation",
    "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"))
MASK = os.path.join(A2, "results", "framework_cohort_masks_v2.json")
OUT = os.path.join(A2, "results", "a2_pooled_fourway_decomp_20261119.json")

ZP3_ENTREZ = "7784"


def load_zp3():
    with open(FPKM, newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        hdr = next(r)
        samples = hdr[1:]
        for row in r:
            if row and row[0].strip() == ZP3_ENTREZ:
                vals = [float(x) if x not in ("", "NA") else np.nan for x in row[1:]]
                return dict(zip(samples, vals))
    raise SystemExit("ZP3 row not found")


def main():
    expr = load_zp3()
    masks = json.load(open(MASK, encoding="utf-8"))
    labels = masks["labels"]
    pooled = masks["B_naive_pooled"]

    cells = {"Pre_R": [], "Pre_N": [], "On_R": [], "On_N": []}
    for s in pooled:
        if s not in expr:
            raise SystemExit(f"sample missing from matrix: {s}")
        tp = "Pre" if "_Pre_" in s else ("On" if "_On_" in s else None)
        lab = labels.get(s)
        if tp is None or lab not in ("R", "N"):
            raise SystemExit(f"cannot classify {s} tp={tp} lab={lab}")
        cells[f"{tp}_{lab}"].append(expr[s])

    for k, v in cells.items():
        print(f"  {k:6s} n={len(v):3d}  median={np.median(v):.4f}")
    nR = len(cells["Pre_R"]) + len(cells["On_R"])
    nN = len(cells["Pre_N"]) + len(cells["On_N"])
    print(f"  total  nR={nR} nN={nN} n={nR + nN}")
    assert nR + nN == 105, "pooled set must be 105 samples"

    # --- four-way partition of U: count pairs where R > N (ties = 0.5) ---
    def wins(a, b):
        a = np.asarray(a)[:, None]
        b = np.asarray(b)[None, :]
        return float(np.sum(a > b) + 0.5 * np.sum(a == b))

    pairs = [
        ("Pre_R_vs_Pre_N", cells["Pre_R"], cells["Pre_N"]),
        ("Pre_R_vs_On_N", cells["Pre_R"], cells["On_N"]),
        ("On_R_vs_Pre_N", cells["On_R"], cells["Pre_N"]),
        ("On_R_vs_On_N", cells["On_R"], cells["On_N"]),
    ]
    parts = []
    U_total = 0.0
    denom = 0
    for name, a, b in pairs:
        w = wins(a, b)
        npairs = len(a) * len(b)
        U_total += w
        denom += npairs
        parts.append({
            "cell": name,
            "n_R": len(a),
            "n_N": len(b),
            "pairs": npairs,
            "pair_weight_pct": round(100.0 * npairs / (nR * nN), 3),
            "wins_R_gt_N": round(w, 2),
            "win_rate_pct": round(100.0 * w / npairs, 4) if npairs else None,
            "contribution_to_U_pct": None,  # filled below
        })

    for p in parts:
        p["contribution_to_U_pct"] = round(100.0 * p["wins_R_gt_N"] / U_total, 3)

    r_pooled = 2.0 * U_total / denom - 1.0

    # --- independent cross-checks ---
    allR = cells["Pre_R"] + cells["On_R"]
    allN = cells["Pre_N"] + cells["On_N"]
    u_scipy = stats.mannwhitneyu(allR, allN, alternative="two-sided")
    P_scipy = float(u_scipy.pvalue)

    # same-time-point only (a basis for comparison, not a test of interaction)
    U_same = parts[0]["wins_R_gt_N"] + parts[3]["wins_R_gt_N"]
    pairs_same = parts[0]["pairs"] + parts[3]["pairs"]
    r_same = 2.0 * U_same / pairs_same - 1.0
    U_cross = parts[1]["wins_R_gt_N"] + parts[2]["wins_R_gt_N"]
    pairs_cross = parts[1]["pairs"] + parts[2]["pairs"]
    r_cross = 2.0 * U_cross / pairs_cross - 1.0

    res = {
        "source": {
            "expr": os.path.basename(FPKM),
            "masks": os.path.basename(MASK),
            "gene": f"ZP3 (Entrez {ZP3_ENTREZ})",
            "framework": "B_naive_pooled",
            "n_samples": nR + nN,
            "n_R": nR,
            "n_N": nN,
        },
        "cells": {k: {"n": len(v), "median": round(float(np.median(v)), 6)} for k, v in cells.items()},
        "four_way": parts,
        "totals": {
            "U": round(U_total, 2),
            "n_pairs": denom,
            "rank_biserial_r": round(r_pooled, 6),
            "scipy_P_two_sided": round(P_scipy, 6),
            "scipy_U": float(u_scipy.statistic),
        },
        "same_timepoint_pooled": {
            "U": round(U_same, 2),
            "n_pairs": pairs_same,
            "rank_biserial_r": round(r_same, 6),
        },
        "cross_timepoint_pooled": {
            "U": round(U_cross, 2),
            "n_pairs": pairs_cross,
            "rank_biserial_r": round(r_cross, 6),
        },
        "note": (
            "Descriptive decomposition of the observed pooled U statistic into four "
            "response x time-point cells. The two cross-time-point cells compare "
            "samples drawn at different treatment states and therefore address a "
            "different question from a baseline association. This partition is not a "
            "test of time-point x response interaction and does not estimate the "
            "proportion of true positives among pooled hits."
        ),
    }

    print("\n=== four-way partition ===")
    for p in parts:
        print(f"  {p['cell']:16s} pairs={p['pairs']:5d} ({p['pair_weight_pct']:5.1f}% of all)  "
              f"win={p['win_rate_pct']:6.2f}%  contributes {p['contribution_to_U_pct']:5.1f}% of U")
    print(f"\n  U_total={U_total:.2f}  n_pairs={denom}  r={r_pooled:.6f}  scipy P={P_scipy:.6f}")
    print(f"  same-timepoint  r={r_same:.6f} (pairs={pairs_same})")
    print(f"  cross-timepoint r={r_cross:.6f} (pairs={pairs_cross})")

    # consistency with frozen v14 value r = -0.362672
    frozen = -0.3626723223753977
    print(f"\n  frozen r (a2_mechanism_decomp B_naive) = {frozen:.6f}")
    print(f"  |difference| = {abs(r_pooled - frozen):.2e}")
    res["consistency_with_frozen"] = {
        "frozen_r": frozen,
        "abs_diff": float(abs(r_pooled - frozen)),
        "matches": bool(abs(r_pooled - frozen) < 1e-6),
    }

    json.dump(res, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"\nwritten: {OUT}")


if __name__ == "__main__":
    main()
