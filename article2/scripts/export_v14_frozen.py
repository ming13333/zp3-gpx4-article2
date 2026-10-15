# -*- coding: utf-8 -*-
"""Export results/a2_v14_cluster_null_frozen.json from results/a2_v14_paired_bootstrap.json.

Why this exists
---------------
The v13 frozen file (a2_v13_cluster_null_frozen.json) and the downstream consumers
(figure 6 generator, supplementary-table builder, number audit) all share one schema.
The v14 permutation run lives in a2_v14_paired_bootstrap.json, whose Part 3 holds the
same four frameworks at 10,000 permutations, and whose Part 1 holds the fixed-42-patient
isolation experiment.

Rather than fork the schema, this exporter restates the v14 numbers in the frozen schema
so that the figure/table/audit tooling keeps working unchanged. The v14 file is therefore
a *derived view*; the authoritative artefact is a2_v14_paired_bootstrap.json, which is what
re-running `a2_v14_paired_bootstrap.py` reproduces.

Sections:
  frozen_caliber_observed : observation counts (identical to v13; the implementation is unchanged)
  vectorized_check        : continuity-on/off counts (identical to v13; copied forward)
  perm                    : the four frameworks at 10,000 permutations  (Part 3)
  perm_fixed42            : the fixed-42-patient isolation experiment    (Part 1)
"""
import io, json, os, sys
sys.stdout.reconfigure(encoding="utf-8")

A2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BOOT = os.path.join(A2, "results", "a2_v14_paired_bootstrap.json")
V13 = os.path.join(A2, "results", "a2_v13_cluster_null_frozen.json")
OUT = os.path.join(A2, "results", "a2_v14_cluster_null_frozen.json")

b = json.load(io.open(BOOT, encoding="utf-8"))
v13 = json.load(io.open(V13, encoding="utf-8"))

FRAMES = ["A_pre", "On_only", "B_naive", "C_patient"]
p3 = b["part3_main_10k"]


def conv(d):
    return {
        "n": d["n"],
        "n_patients": d["n_patients"],
        "nR_patients": d["nR_patients"],
        "observed_vectorized": d["observed"],
        "observed_scipy": d["observed"],
        "null_mean": d["null_mean"],
        "null_sd": d.get("null_sd"),
        "null_pct": d["null_pct"],
        "perm_P": d["perm_P"],
        "tail_hits": d.get("tail_hits"),
        "nperm": d.get("nperm", 10000),
        "excess": d["excess"],
    }


perm = {}
for f in FRAMES:
    perm[f] = conv(p3[f])
    # observed_scipy is the scipy-caliber observation count, which is independent of the
    # permutation count and therefore carried over from the v13 frozen file verbatim.
    perm[f]["observed_scipy"] = v13["frozen_caliber_observed"][f]["sig"]

FIXMAP = {"P_pre42": "P_pre42", "P_on42": "P_on42", "P_pool84": "P_pool84"}
fixed42 = {}
for k, v in b["part1_fixed42"].items():
    fixed42[FIXMAP.get(k, k)] = conv(v)
    fixed42[FIXMAP.get(k, k)]["n_patients"] = v["n_patients"]
    fixed42[FIXMAP.get(k, k)]["nR_patients"] = v["nR_patients"]

out = {
    "_provenance": {
        "derived_from": "results/a2_v14_paired_bootstrap.json",
        "generator": "scripts/export_v14_frozen.py",
        "seed": b["seed"],
        "n_genes": b["n_genes"],
        "note": ("Derived view restating the v14 numbers in the v13 frozen schema. "
                 "Authoritative artefact = a2_v14_paired_bootstrap.json."),
    },
    "frozen_caliber_observed": v13["frozen_caliber_observed"],
    "vectorized_check": v13["vectorized_check"],
    "perm": perm,
    "perm_fixed42": fixed42,
}

json.dump(out, io.open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("WROTE", OUT)
print()
print("%-10s %6s %8s %10s %8s %9s %8s" % ("frame", "n", "observed", "null_mean", "null%", "permP", "excess"))
for f in FRAMES:
    d = perm[f]
    print("%-10s %6d %8d %10.1f %8.2f %9.4f %8d" % (f, d["n"], d["observed_vectorized"],
                                                    d["null_mean"], d["null_pct"], d["perm_P"], d["excess"]))
print()
print("%-10s %6s %8s %10s %8s %9s %8s" % ("fixed42", "n", "observed", "null_mean", "null%", "permP", "excess"))
for f in ["P_pre42", "P_on42", "P_pool84"]:
    d = fixed42[f]
    print("%-10s %6d %8d %10.1f %8.2f %9.4f %8d" % (f, d["n"], d["observed_vectorized"],
                                                    d["null_mean"], d["null_pct"], d["perm_P"], d["excess"]))
