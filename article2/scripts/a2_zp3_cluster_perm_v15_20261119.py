# -*- coding: utf-8 -*-
"""
A2 v15: bring the ZP3-level patient-cluster permutation P onto (b+1)/(B+1).

Background: scripts/a2_cluster_permutation.py (frozen 2026-10-14) reports b/B, whereas the
transcriptome-level P values use (b+1)/(B+1). Two estimators in one manuscript cannot both be
described accurately in Methods. This script replays the permutation distribution with exactly
the same random stream and call order as the original (seed 20261014, order B -> A -> On) and
reports b, B, b/B and (b+1)/(B+1), asserting that b matches the frozen values
(558 / 7145 / 262): only the estimator changes, not the permutation.
Output: results/a2_zp3_cluster_perm_v15_20261119.json
"""
import csv, json, os, sys
from collections import defaultdict

import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")

_HERE = os.path.dirname(os.path.abspath(__file__))
A2 = os.path.dirname(_HERE)
ROOT = os.path.dirname(A2)
FPKM = os.environ.get("A2_GSE91061_FPKM", os.path.join(
    ROOT, "output", "gse91061_validation",
    "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"))
MASK = os.path.join(A2, "results", "framework_cohort_masks_v2.json")
OUT = os.path.join(A2, "results", "a2_zp3_cluster_perm_v15_20261119.json")

NPERM = 20000
SEED = 20261014
FROZEN_B = {"B": 558, "A": 7145, "On": 262}


def load_zp3():
    with open(FPKM, newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        hdr = next(r)
        samples = [h.strip() for h in hdr[1:]]
        for row in r:
            if row and row[0].strip().strip('"') == "7784":
                return dict(zip(samples, [float(x) for x in row[1:]]))
    raise SystemExit("ZP3 (Entrez 7784) not found")


def main():
    zp3 = load_zp3()
    masks = json.load(open(MASK, encoding="utf-8"))
    labels = masks["labels"]
    sets = {"B": masks["B_naive_pooled"], "A": masks["A_pre"]}
    sets["On"] = [s for s in sets["B"] if "_On_" in s]

    rng = np.random.default_rng(SEED)
    out = dict(seed=SEED, nperm=NPERM,
               note="replay of a2_cluster_permutation.py rng stream; estimator changed to (b+1)/(B+1)")

    for tag in ("B", "A", "On"):                     # identical call order to the frozen script
        samples = sets[tag]
        bypt = defaultdict(list)
        for s in samples:
            bypt[s.split("_")[0]].append(s)
        pids = list(bypt.keys())
        plab = {p: labels[bypt[p][0]] for p in pids}
        nR = sum(1 for p in pids if plab[p] == "R")

        vals = np.array([zp3[s] for s in samples])
        Rk = stats.rankdata(vals)
        n = vals.size
        # patient index per sample, in the same sample order as `samples`
        pid_of = np.array([pids.index(s.split("_")[0]) for s in samples])

        obs_mask = np.array([labels[s] == "R" for s in samples])
        nR_s = int(obs_mask.sum()); nN_s = n - nR_s
        U = Rk[obs_mask].sum() - nR_s * (nR_s + 1.0) / 2.0
        r_obs = 2.0 * U / (nR_s * nN_s) - 1.0

        idx = np.arange(len(pids))
        b = 0
        null = np.empty(NPERM)
        for i in range(NPERM):
            perm = rng.permutation(idx)
            Rset = set(idx[perm][:nR].tolist())
            m = np.array([p in Rset for p in pid_of])
            nR_i = int(m.sum()); nN_i = n - nR_i
            Ui = Rk[m].sum() - nR_i * (nR_i + 1.0) / 2.0
            ri = 2.0 * Ui / (nR_i * nN_i) - 1.0
            null[i] = ri
            if abs(ri) >= abs(r_obs):   # same comparison as the frozen implementation (no tolerance), else the one exact tie is over-counted
                b += 1

        out[tag] = dict(n_samples=n, n_patients=len(pids), nR_patients=nR,
                        r=float(r_obs), b=int(b), B=NPERM,
                        p_plain=b / NPERM, p_addone=(b + 1.0) / (NPERM + 1.0),
                        frozen_b=FROZEN_B[tag], b_matches=bool(b == FROZEN_B[tag]))
        print("%-3s n=%3d  r=%+.4f  b=%5d/%d  b/B=%.5f  (b+1)/(B+1)=%.5f  frozen_b=%d %s"
              % (tag, n, r_obs, b, NPERM, b / NPERM, (b + 1.0) / (NPERM + 1.0),
                 FROZEN_B[tag], "OK" if b == FROZEN_B[tag] else "MISMATCH"))

    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("written:", OUT)


if __name__ == "__main__":
    main()
