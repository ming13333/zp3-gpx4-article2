# -*- coding: utf-8 -*-
"""
A2 v16: the denominators behind the S3 and S4 percentages.

The published summaries were computed over different gene subsets depending on
which target was being scored, and the text did not say so. In the frozen
summaries:

  * "TypeI" rows are averaged over all 500 genes per replicate;
  * "power" rows are averaged over the 50 genes per replicate that carry the
    injected non-zero term (10% of 500);
  * the pooled procedure appears in both roles because it is scored against each
    target's truth: all-gene rows come from the target that is null in that
    scenario (S3 change, S4 baseline), signal rows from the target that is
    non-null.

This script recovers the quantity that was never reported - the rejection share
among the 450 genes without an injected term - from the identity

    all = (50/500) * signal + (450/500) * null

which is exact, because the per-replicate mean over genes is a mean of two
disjoint blocks and k/G is exactly 1/10. No new simulation is run; the script
asserts the block sizes from the frozen JSON (n_genes = 50 and 500) and the
replicate count, and writes the per-cell decomposition.

Output: results/a2_p2_partition_audit_20261120.json
"""
import glob, json, os, sys

sys.stdout.reconfigure(encoding="utf-8")
_HERE = os.path.dirname(os.path.abspath(__file__))
A2 = os.path.dirname(_HERE)
RES = os.path.join(A2, "results")
OUT = os.path.join(RES, "a2_p2_partition_audit_20261120.json")

G, K = 500, 50
FRAC = K / float(G)


def load(p):
    d = json.load(open(p, encoding="utf-8"))
    c = d["cfg"]
    rows = {}
    for s in d["summary"]:
        if "rate" in s:
            rows[(s["scenario"], s["target"], s["method"])] = s
    return c, rows


def cell(c, rows):
    """signal / all / null shares for the pooled cluster-corrected test."""
    out = {}
    for scen in ("S3", "S4"):
        sig = rows[(scen, "T1" if scen == "S3" else "T2", "pooled_cluster")]
        allg = rows[(scen, "T3" if scen == "S3" else "T1", "pooled_cluster")]
        assert sig["n_genes"] == K and allg["n_genes"] == G, (scen, sig["n_genes"], allg["n_genes"])
        out[scen] = dict(signal=sig["rate"], all_genes=allg["rate"],
                         null=(allg["rate"] - FRAC * sig["rate"]) / (1 - FRAC),
                         n_genes_signal=K, n_genes_null=G - K)
    d3 = rows[("S3", "T3", "delta_mwu")]
    assert d3["n_genes"] == G
    out["S3"]["delta_mwu_all"] = d3["rate"]
    d4 = rows[("S4", "T3", "delta_mwu")]
    assert d4["n_genes"] == K
    out["S4"]["delta_mwu_signal"] = d4["rate"]
    out["S4"]["mw_on_signal"] = rows[("S4", "T2", "mw_on")]["rate"]
    out["S4"]["pooled_nocorr_signal"] = rows[("S4", "T2", "mw_pooled")]["rate"]
    out["S1"] = dict(pooled=rows[("S1", "T1", "mw_pooled")]["rate"],
                     cluster=rows[("S1", "T1", "pooled_cluster")]["rate"])
    return out


def main():
    files = [f for f in glob.glob(os.path.join(RES, "a2_p2_sim_core_20261118_full*.json"))
             if "smoke" not in f]
    cells, mains = [], []
    for f in files:
        c, rows = load(f)
        rec = dict(np=c["n_pat"], rp=c["resp_prop"], icc=c["icc"],
                   n_outer=c["n_outer"], source=os.path.basename(f), **cell(c, rows))
        (mains if c["n_outer"] == 1000 else cells).append(rec)
    assert len(cells) == 18, len(cells)
    assert len(mains) == 1
    main_cell = mains[0]

    def rng(recs, scen, key):
        v = [r[scen][key] for r in recs]
        return [min(v), max(v)]

    summary = {}
    for scen in ("S3", "S4", "S1"):
        summary[scen] = {k: rng(cells, scen, k) for k in cells[0][scen]}
        summary[scen + "_main"] = dict(main_cell[scen])

    out = dict(provenance=dict(
        inputs="results/a2_p2_sim_core_20261118_full*.json (18 grid cells + main)",
        identity="all = (K/G)*signal + (1-K/G)*null, K=50, G=500",
        note="derived, not simulated; block sizes asserted from n_genes"),
        grid_cells=len(cells), summary=summary,
        per_cell=[{k: v for k, v in r.items() if k != "source"} for r in cells])
    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    def pct(x):
        return "%.2f%%" % (100 * x)

    print("S3  pooled cluster-corrected, per replicate")
    for k, lab in (("all_genes", "all 500 genes"), ("signal", "50 signal genes"),
                   ("null", "450 null genes"), ("delta_mwu_all", "direct change test")):
        print("   %-20s %s  (main %s)" % (lab, " - ".join(pct(x) for x in summary["S3"][k]),
                                          pct(summary["S3_main"][k])))
    print("S4  per replicate")
    for k in ("mw_on_signal", "delta_mwu_signal", "pooled_nocorr_signal", "signal", "null"):
        print("   %-20s %s  (main %s)" % (k, " - ".join(pct(x) for x in summary["S4"][k]),
                                          pct(summary["S4_main"][k])))
    print("S1  %s / %s" % (pct(summary["S1_main"]["pooled"]), pct(summary["S1_main"]["cluster"])))
    print("written:", OUT)


if __name__ == "__main__":
    main()
