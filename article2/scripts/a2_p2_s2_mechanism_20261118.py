# -*- coding: utf-8 -*-
"""
P2 mechanism experiment: is the S2 variance mis-calibration really a mismatch between patient-level permutation and the iid variance formula?

Reviewer point 6 (2026-11-18) required demoting any claim that S2 'can only be caused by ...',
and checking the mechanism with a tau x within-patient-correlation factorial experiment, comparing the empirical variance of pooled U against the iid approximation.

Design (a=0, delta=0, i.e. data under the complete null with only a common time effect tau injected):
  for each (tau, ICC):
    1. generate one dataset (G genes x 2P samples, patient-major)
    2. under patient-level permutation, repeat B times and compute the standardised statistic z_g on the pooled samples
    3. empirical variance s2 = mean_g Var_perm(z_g)
       -- the iid variance formula assumes s2 = 1; s2 > 1 means the realised null variance exceeds the assumed one
    4. implied two-sided rejection rate from s2: 2*(1-Phi(1.96/s))
    5. compare with the S2 Type I error measured in the main simulation

Prediction: s2 close to 1 at tau=0; s2 rises monotonically with tau and ICC.

Output: results/a2_p2_s2_mechanism_20261118.json
"""
import json, os, sys, time
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
_HERE = os.path.dirname(os.path.abspath(__file__))
A2 = os.path.dirname(_HERE)
OUT = os.path.join(A2, "results", "a2_p2_s2_mechanism_20261118.json")

G, P, RESP = 500, 80, 0.25
NPERM = 2000
SEED = 20261118
rng = np.random.default_rng(SEED)

def gen(tau, icc):
    nR = int(round(P * RESP))
    lab = np.zeros(P, dtype=bool); lab[rng.choice(P, nR, replace=False)] = True
    sb, se = np.sqrt(icc), np.sqrt(1.0 - icc)
    b  = rng.normal(0, sb, (G, P))
    e0 = rng.normal(0, se, (G, P)); e1 = rng.normal(0, se, (G, P))
    pre = b + e0
    on  = b + tau + e1
    X = np.empty((G, 2 * P)); X[:, 0::2] = pre; X[:, 1::2] = on
    return X, lab

def tie_term(X):
    n = X.shape[1]
    Xs = np.sort(X, axis=1)
    hastie = np.any(np.diff(Xs, axis=1) == 0, axis=1)
    t = np.zeros(X.shape[0])
    for i in np.where(hastie)[0]:
        _, c = np.unique(X[i], return_counts=True)
        t[i] = np.sum(c.astype(float) ** 3 - c) / (n * (n - 1.0))
    return t

def pooled_z_stats(X, lab, nperm):
    """Return (z_obs, empirical variance s2). z is standardised by the iid variance formula."""
    n = X.shape[1]
    Rk = stats.rankdata(X, axis=1)
    tie = tie_term(X)
    nconst = (n + 1.0) - tie
    def z(isR_s):
        nR = int(isR_s.sum()); nN = n - nR
        U = Rk[:, isR_s].sum(axis=1) - nR * (nR + 1.0) / 2.0
        mu = nR * nN / 2.0
        var = nR * nN / 12.0 * nconst
        dd = U - mu
        dd = np.where(dd > 0, dd - 0.5, np.where(dd < 0, dd + 0.5, dd))
        return dd / np.sqrt(var)
    zobs = z(np.repeat(lab, 2))
    Z = np.empty((nperm, X.shape[0]))
    for k in range(nperm):
        Z[k] = z(np.repeat(rng.permutation(lab), 2))
    return zobs, Z

rows = []
print("=" * 84)
print("S2 mechanism experiment: effect of a common time effect tau and within-patient correlation ICC on the pooled null variance")
print("  (a=0, delta=0; patient-level permutation B=%d; s2=1 means agreement with the iid variance formula)" % NPERM)
print("=" * 84)
print("%-6s %-6s %-10s %-10s %-12s %-14s" % ("tau", "ICC", "s2", "sqrt(s2)", "implied TypeI", "note"))
for tau in [0.0, 0.5, 1.0, 2.0]:
    for icc in [0.2, 0.5, 0.8]:
        t0 = time.time()
        X, lab = gen(tau, icc)
        zobs, Z = pooled_z_stats(X, lab, NPERM)
        s2 = float(np.mean(Z.var(axis=0, ddof=1)))
        s = np.sqrt(s2)
        implied = float(2 * (1 - stats.norm.cdf(1.96 / s)))
        rows.append(dict(tau=tau, icc=icc, s2=s2, sd=s, implied_type1=implied,
                         var_z_obs=float(zobs.var())))
        print("%-6.1f %-6.1f %-10.4f %-10.4f %-12.4f %.0fs" % (tau, icc, s2, s, implied, time.time()-t0))

json.dump(dict(G=G, P=P, resp=RESP, nperm=NPERM, seed=SEED, rows=rows),
          open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
print("\nwritten:", OUT)
