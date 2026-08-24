#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
freeze_subtype_response_a2.py
Freeze GSE91061 Pre baseline (n=49) molecular subtype-stratified ZP3-response Mann-Whitney P values,
closing the Fig5 right-panel "Immune / Keratin" two-bar (P=0.8001 / 0.7692)
"claimed on figure but no frozen source" gap (same-source audit as freeze_standard_logistic_a2.py).

Replicates the "stratify by molecular subtype" convention of gse91061_molecular_subtype.py::analyze_response:
  - Expression data: GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv (Entrez ID rows)
  - ZP3 Entrez = 7784 (2026-08-17 audit fix)
  - Pre baseline samples only (visit == 'Pre')
  - binary_response: PRCR/PR/CR = Response; SD/PD = Non-response
  - Molecular subtypes: Immune / Keratin / MITF-high / MITF-low (thresholded by score median, assign_subtypes)
  - Per subtype: Mann-Whitney U (ZP3 | Response vs Non-response, two-sided)
    - Either group < 2 samples -> not estimable

Output: article2/results/a2_gse91061_subtype_response.csv
Columns: subtype, n_total, n_response, n_nonresponse, zp3_mean_response,
    zp3_mean_nonresponse, mwu_p, estimable
"""
import os, sys
import pandas as pd
import numpy as np
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
DATA_DIR = os.path.join(ROOT, "output", "gse91061_validation")
RESULTS_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

sys.path.insert(0, SCRIPT_DIR)
import gse91061_molecular_subtype as sub  # reuse helpers


def load_clinical_aligned():
    """Replicate main() data preparation and assign molecular subtypes."""
    expr = pd.read_csv(os.path.join(DATA_DIR, "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"), index_col=0)

    sample_titles, gsm_ids, visits, responses = [], [], [], []
    with open(os.path.join(DATA_DIR, "GSE91061_series_matrix.txt"), "r", encoding="utf-8") as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("!Sample_title"):
            parts = line.strip().split("\t")
            if len(parts) > 1:
                sample_titles = [p.strip('"') for p in parts[1:]]
        elif line.startswith("!Sample_geo_accession"):
            parts = line.strip().split("\t")
            if len(parts) > 1:
                gsm_ids = [p.strip('"') for p in parts[1:]]
    for line in lines:
        if line.startswith("!Sample_characteristics_ch1"):
            parts = line.strip().split("\t")
            if len(parts) > 1 and "visit" in parts[1].lower():
                visits = [p.split(": ")[-1].strip('"') for p in parts[1:] if ": " in p]
            elif len(parts) > 1 and "response" in parts[1].lower():
                responses = [p.split(": ")[-1].strip('"') for p in parts[1:] if ": " in p]

    clinical = pd.DataFrame({
        "sample_title": sample_titles, "gsm_id": gsm_ids,
        "visit": visits, "response": responses,
    })
    common_samples = [t for t in clinical["sample_title"] if t in expr.columns]
    clinical_aligned = clinical[clinical["sample_title"].isin(common_samples)].set_index("sample_title")
    expr_aligned = expr[common_samples]

    entrez_to_symbol = sub.get_gene_mapping()
    expr_symbol = expr_aligned.rename(index=entrez_to_symbol)
    if expr_symbol.index.duplicated().any():
        expr_symbol = expr_symbol.groupby(expr_symbol.index).mean()

    subtype_scores = sub.compute_subtype_scores(expr_symbol)
    clinical_aligned = clinical_aligned.join(subtype_scores)
    clinical_aligned = sub.assign_subtypes(clinical_aligned)

    immune_features = sub.compute_immune_features(expr_symbol)
    clinical_aligned = clinical_aligned.join(immune_features)

    zp3 = expr_aligned.loc[7784]
    clinical_aligned["ZP3"] = zp3[clinical_aligned.index]
    return clinical_aligned


def subtype_mwu():
    df = load_clinical_aligned()
    pre_df = df[df["visit"] == "Pre"].copy()
    response_map = {"PRCR": "Response", "PR": "Response", "CR": "Response",
                    "SD": "Non-response", "PD": "Non-response"}
    pre_df["binary_response"] = pre_df["response"].map(response_map)
    valid = pre_df[pre_df["binary_response"].notna()].copy()
    n_pre = len(valid)

    rows = []
    for st in sorted(valid["Molecular_subtype"].unique()):
        sub_df = valid[valid["Molecular_subtype"] == st]
        sr = sub_df[sub_df["binary_response"] == "Response"]
        sn = sub_df[sub_df["binary_response"] == "Non-response"]
        rec = {
            "subtype": st,
            "n_total": len(sub_df),
            "n_response": len(sr),
            "n_nonresponse": len(sn),
            "zp3_mean_response": round(float(sr["ZP3"].mean()), 4) if len(sr) else None,
            "zp3_mean_nonresponse": round(float(sn["ZP3"].mean()), 4) if len(sn) else None,
        }
        if len(sr) >= 2 and len(sn) >= 2:
            stat, p = stats.mannwhitneyu(sr["ZP3"], sn["ZP3"], alternative="two-sided")
            rec["mwu_p"] = round(float(p), 4)
            rec["estimable"] = True
        else:
            rec["mwu_p"] = None
            rec["estimable"] = False
        rows.append(rec)

    out = pd.DataFrame(rows, columns=["subtype", "n_total", "n_response", "n_nonresponse",
                                      "zp3_mean_response", "zp3_mean_nonresponse", "mwu_p", "estimable"])
    return out, n_pre


def main():
    out, n_pre = subtype_mwu()
    print(f"GSE91061 Pre baseline n={n_pre}")
    print(out.to_string(index=False))

    declared = {"Immune": 0.8001, "Keratin": 0.7692}
    ok = True
    for st, dp in declared.items():
        row = out[out["subtype"] == st]
        if row.empty or not bool(row["estimable"].iloc[0]):
            print(f"  {st}: MISSING / not estimable! declared={dp}")
            ok = False
            continue
        got = float(row["mwu_p"].iloc[0])
        match = abs(got - dp) < 0.001
        ok = ok and match
        print(f"  {st}: frozen p={got:.4f}  Fig5 declared={dp}  {'OK' if match else 'MISMATCH!'}")

    path = os.path.join(RESULTS_DIR, "a2_gse91061_subtype_response.csv")
    out.to_csv(path, index=False)
    print(f"\nWROTE {path}")
    print("ALL MATCH" if ok else "MISMATCH DETECTED — DO NOT SHIP")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
