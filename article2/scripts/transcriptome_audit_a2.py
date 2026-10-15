# -*- coding: utf-8 -*-
"""
A2 supplemental Step-4: transcriptome-wide audit of temporal sampling bias
and pseudo-replication in GSE91061.

For every gene in the frozen expression matrix (~22k Entrez rows), compute a
two-sided Mann-Whitney responder-vs-non-responder association under:
  Framework A  pre-treatment only            n=49  (10 R / 39 N)
  Framework B  naive pooled pre+on samples   n=105 (23 R / 82 N)  [pseudo-replicated]
  Framework C  patient-level, pre-else-on    n=62  (14 R / 48 N)

Quantify:
  - # genes P<0.05 in A vs B vs C (vs nominal 5% expectation)
  - "B-gained" = P<0.05 in B but NOT in A  (temporal sampling + pseudo-replication)
  - "B-gained & C-lost" = the subset where the naive signal also disappears under
    patient-level aggregation (structurally spurious for a baseline question)
  - BH-FDR q-values per framework; genes q<0.05 in B vs A
  - IFN-gamma / inflammation gene-set overlap among B-gained (Fisher exact),
    with every gene-set ID verified live against NCBI esummary (red line:
    no memory-only Entrez IDs)

Outputs (new *_v2 files only; frozen CSVs untouched):
  article2/results/genome_wide_pvalues_v2.csv
  article2/results/Supp_Table_S6_v2.csv      (top 50 B-gained, annotated)
  article2/results/inflation_summary_a2.json
"""
import os, json, sys, time
import numpy as np
import pandas as pd
from scipy import stats
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
G9DIR = os.path.join(ROOT, "output", "gse91061_validation")
RES = os.path.join(ROOT, "article2", "results")
MASKS = os.path.join(RES, "framework_cohort_masks_v2.json")
RESP_MAP = {"PRCR": "R", "PR": "R", "CR": "R", "SD": "N", "PD": "N"}

IFN_SET = {  # (entrez, expected_symbol) -- verified at runtime against NCBI
    "IFNG": (3458, "IFNG"), "CXCL9": (4283, "CXCL9"), "CXCL10": (3627, "CXCL10"),
    "CXCL11": (6373, "CXCL11"), "GBP1": (2633, "GBP1"), "GBP2": (2634, "GBP2"),
    "STAT1": (6772, "STAT1"), "IDO1": (3620, "IDO1"), "IRF1": (3659, "IRF1"),
    # 2026-08-30 NCBI [Symbol]+human[ORGN] verified: TAP1=6890 (6863=TAC1),
    # PSMB9=5698 (5683=PSMA2), GZMA=3001 (3099=HK2)
    "TAP1": (6890, "TAP1"), "PSMB9": (5698, "PSMB9"), "CCL5": (6352, "CCL5"),
    "GZMA": (3001, "GZMA"), "GZMB": (3002, "GZMB"), "PRF1": (5551, "PRF1"),
    "CD8A": (925, "CD8A"), "CD274": (29126, "CD274"), "PDCD1": (5133, "PDCD1"),
    "LAG3": (3902, "LAG3"), "TIGIT": (201633, "TIGIT"),
}

def esummary_symbols(ids):
    """Batch Entrez -> (symbol, status). One request, <=20 ids."""
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
           "?db=gene&retmode=json&id=" + ",".join(str(i) for i in ids))
    with urllib.request.urlopen(url, timeout=60) as r:
        d = json.load(r)
    out = {}
    for i in ids:
        rec = d["result"].get(str(i))
        if rec:
            out[i] = (rec.get("nomenclaturesymbol") or rec.get("name"), rec.get("description", ""))
        else:
            out[i] = (None, None)
    return out

def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(len(p), np.nan)
    idx = np.where(ok)[0]
    pv = p[idx]
    order = np.argsort(pv)
    n = len(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    q[idx] = out
    return q

def main():
    masks = json.load(open(MASKS, encoding="utf-8"))
    A = masks["A_pre"]; B = masks["B_naive_pooled"]; C = masks["C_patient_pre_else_on"]
    lab = masks["labels"]
    assert (len(A), len(B), len(C)) == (49, 105, 62)

    expr = pd.read_csv(os.path.join(G9DIR, "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"),
                       index_col=0)
    X = expr.to_numpy(float)
    yA = np.array([lab[s] for s in A])
    yB = np.array([lab[s] for s in B])
    yC = np.array([lab[s] for s in C])
    cols = {s: i for i, s in enumerate(expr.columns)}
    iA = [cols[s] for s in A]; iB = [cols[s] for s in B]; iC = [cols[s] for s in C]

    def fw_test(idx, y):
        ia = [i for i, g in zip(idx, y) if g == "R"]
        ib = [i for i, g in zip(idx, y) if g == "N"]
        sub = X[:, ia], X[:, ib]
        ps = np.full(X.shape[0], np.nan)
        rs = np.full(X.shape[0], np.nan)
        n1, n2 = len(ia), len(ib)
        for g in range(X.shape[0]):
            x1, x2 = sub[0][g], sub[1][g]
            if np.all(x1 == x1[0]) and np.all(x2 == x2[0]) and x1[0] == x2[0]:
                ps[g] = 1.0; rs[g] = 0.0
                continue
            try:
                u, p = stats.mannwhitneyu(x1, x2, alternative="two-sided")
                ps[g] = p; rs[g] = 2 * u / (n1 * n2) - 1
            except ValueError:
                pass
        return ps, rs

    print("testing frameworks A/B/C across", X.shape[0], "genes ...")
    pA, rA = fw_test(iA, yA)
    pB, rB = fw_test(iB, yB)
    pC, rC = fw_test(iC, yC)
    qA, qB, qC = bh(pA), bh(pB), bh(pC)

    sigA, sigB, sigC = pA < 0.05, pB < 0.05, pC < 0.05
    gained = sigB & ~sigA                     # apparent signal from B only
    gained_lost = gained & ~sigC              # also vanishes at patient level
    n_tested = int(np.sum(~np.isnan(pA) & ~np.isnan(pB) & ~np.isnan(pC)))

    tab = pd.DataFrame({"Entrez": expr.index.astype(int),
                        "P_A": pA, "r_A": rA, "Q_A": qA,
                        "P_B": pB, "r_B": rB, "Q_B": qB,
                        "P_C": pC, "r_C": rC, "Q_C": qC,
                        "sig_A": sigA, "sig_B": sigB, "sig_C": sigC,
                        "B_gained": gained, "B_gained_C_lost": gained_lost})
    tab.to_csv(os.path.join(RES, "genome_wide_pvalues_v2.csv"), index=False)
    print("wrote genome_wide_pvalues_v2.csv")

    # Top 50 B-gained by P_B, annotate via NCBI esummary (live verification)
    top = tab[tab.B_gained].nsmallest(50, "P_B")
    ids = top.Entrez.tolist()
    sym = {}
    for k in range(0, len(ids), 20):
        sym.update(esummary_symbols(ids[k:k+20]))
        time.sleep(0.35)
    s6 = top.copy()
    s6["Symbol"] = [sym.get(i, (None, None))[0] for i in s6.Entrez]
    s6["Description"] = [sym.get(i, (None, None))[1] for i in s6.Entrez]
    s6.to_csv(os.path.join(RES, "Supp_Table_S6_v2.csv"), index=False)
    print("wrote Supp_Table_S6_v2.csv (top 50 B-gained with symbols)")

    # IFN-gamma gene-set check: verify IDs live, then Fisher overlap among gained
    ids_chk = [v[0] for v in IFN_SET.values()]
    chk = esummary_symbols(ids_chk)
    bad = {g: (v, chk.get(v[0], (None,))[0]) for g, v in IFN_SET.items()
           if chk.get(v[0], (None,))[0] != v[1]}
    if bad:
        print("!! IFN gene-set ID verification FAILED:", bad); sys.exit(1)
    print("IFN gene-set: 20/20 Entrez IDs verified against NCBI")
    set_rows = np.isin(tab.Entrez.values, ids_chk)
    a = int(np.sum(gained & set_rows)); b = int(np.sum(gained & ~set_rows))
    c = int(np.sum(~gained & set_rows)); d0 = int(np.sum(~gained & ~set_rows))
    orr, pf = stats.fisher_exact([[a, b], [c, d0]], alternative="greater")
    print(f"B-gained IFN-set overlap: {a}/{int(set_rows.sum())} genes (Fisher p={pf:.3g}, OR={orr:.1f})")

    summary = {
        "n_genes_in_matrix": int(X.shape[0]),
        "n_tested_all_frameworks": n_tested,
        "nominal_expectation_p05": round(0.05 * n_tested),
        "sig_P05": {"A_pre": int(sigA.sum()), "B_naive": int(sigB.sum()),
                    "C_patient": int(sigC.sum())},
        "sig_q05": {"A_pre": int(np.nansum(qA < 0.05)), "B_naive": int(np.nansum(qB < 0.05)),
                    "C_patient": int(np.nansum(qC < 0.05))},
        "B_gained_vs_A": int(gained.sum()),
        "B_gained_and_lost_in_C": int(gained_lost.sum()),
        "inflation_fold_B_over_A": round(float(sigB.sum()) / max(sigA.sum(), 1), 2),
        "IFN_set_overlap": {"in_gained": a, "set_size": int(set_rows.sum()),
                            "fisher_p": float(pf), "odds_ratio": float(orr)},
        "cohorts": {"A": "n=49 (10R/39N)", "B": "n=105 samples (23R/82N)",
                    "C": "n=62 patients (14R/48N)"},
    }
    with open(os.path.join(RES, "inflation_summary_a2.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
