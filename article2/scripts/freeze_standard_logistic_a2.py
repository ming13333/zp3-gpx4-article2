#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
freeze_standard_logistic_a2.py
冻结 GSE91061 Pre 基线 (n=49) 标准 logistic 回归结果，封堵 Fig5 左面板
"logistic univariable / +IFN-g / +Supp" 三柱 (P=0.4545/0.5154/0.5772) 的
"图上声明但无冻结源" 缺口。

复刻 gse91061_molecular_subtype.py::multivariate_analysis 的精确口径:
  - 表达数据: GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv (Entrez ID 行)
  - ZP3 Entrez = 7784 (2026-08-17 审计修复)
  - 仅 Pre 基线样本 (visit == 'Pre')
  - binary_response: PRCR/PR/CR = 1 (Response); SD/PD = 0 (Non-response)
  - 3 个模型: (1) ZP3 单独; (2) ZP3 + IFN-gamma score; (3) ZP3 + immunosuppression score
  - 协变量均按列标准化 (z-score)

输出: article2/results/a2_gse91061_logistic_standard.csv
列: model, term, n, OR, ci_low, ci_high, p
"""
import os, sys
import pandas as pd
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
DATA_DIR = os.path.join(ROOT, "output", "gse91061_validation")
RESULTS_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

sys.path.insert(0, SCRIPT_DIR)
import gse91061_molecular_subtype as sub  # reuse get_gene_mapping / compute_immune_features


def load_clinical_aligned():
    """复刻 main() 的数据准备 (lines 64-137)，返回带 ZP3 + 免疫评分的临床表。"""
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

    immune_features = sub.compute_immune_features(expr_symbol)
    clinical_aligned = clinical_aligned.join(immune_features)

    zp3_entrez = 7784
    zp3 = expr_aligned.loc[zp3_entrez]
    clinical_aligned["ZP3"] = zp3[clinical_aligned.index]
    return clinical_aligned


def standard_logistic():
    from statsmodels.formula.api import logit

    df = load_clinical_aligned()
    pre_df = df[df["visit"] == "Pre"].copy()
    response_map = {"PRCR": 1, "PR": 1, "CR": 1, "SD": 0, "PD": 0}
    pre_df["binary_response"] = pre_df["response"].map(response_map)
    valid = pre_df[pre_df["binary_response"].notna()].copy()
    n = len(valid)

    valid["ZP3_scaled"] = (valid["ZP3"] - valid["ZP3"].mean()) / valid["ZP3"].std()
    if "IFN_gamma_score" in valid.columns:
        valid["IFN_scaled"] = (valid["IFN_gamma_score"] - valid["IFN_gamma_score"].mean()) / valid["IFN_gamma_score"].std()
    if "Immunosuppression_score" in valid.columns:
        valid["Supp_scaled"] = (valid["Immunosuppression_score"] - valid["Immunosuppression_score"].mean()) / valid["Immunosuppression_score"].std()

    rows = []
    models = [
        ("Model1_ZP3_only", "binary_response ~ ZP3_scaled", ["ZP3_scaled"]),
        ("Model2_ZP3_IFNg", "binary_response ~ ZP3_scaled + IFN_scaled", ["ZP3_scaled", "IFN_scaled"]),
        ("Model3_ZP3_Supp", "binary_response ~ ZP3_scaled + Supp_scaled", ["ZP3_scaled", "Supp_scaled"]),
    ]
    for mname, formula, terms in models:
        try:
            m = logit(formula, data=valid).fit(disp=False)
            for term in terms:
                orr = float(np.exp(m.params[term]))
                ci = m.conf_int().loc[term]
                rows.append({
                    "model": mname, "term": term, "n": n,
                    "OR": round(orr, 4),
                    "ci_low": round(float(np.exp(ci[0])), 4),
                    "ci_high": round(float(np.exp(ci[1])), 4),
                    "p": round(float(m.pvalues[term]), 4),
                })
        except Exception as e:
            print(f"  {mname} failed: {e}")

    out = pd.DataFrame(rows, columns=["model", "term", "n", "OR", "ci_low", "ci_high", "p"])
    return out, n


def main():
    out, n = standard_logistic()
    print(f"GSE91061 Pre baseline n={n}")
    print(out.to_string(index=False))

    # 自检: Fig5 左面板声明值
    declared = {"Model1_ZP3_only": 0.4545, "Model2_ZP3_IFNg": 0.5154, "Model3_ZP3_Supp": 0.5772}
    ok = True
    for mname, dp in declared.items():
        got = float(out.loc[out["model"] == mname, "p"].iloc[0])
        match = abs(got - dp) < 0.001
        ok = ok and match
        print(f"  {mname}: frozen p={got:.4f}  Fig5 declared={dp}  {'OK' if match else 'MISMATCH!'}")

    path = os.path.join(RESULTS_DIR, "a2_gse91061_logistic_standard.csv")
    out.to_csv(path, index=False)
    print(f"\nWROTE {path}")
    print("ALL MATCH" if ok else "MISMATCH DETECTED — DO NOT SHIP")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
