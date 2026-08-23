import os as _os
def _project_root():
    d = _os.path.dirname(_os.path.abspath(__file__))
    while True:
        if _os.path.isdir(_os.path.join(d, "output")):
            return d
        p = _os.path.dirname(d)
        if p == d:
            break
        d = p
    return _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
ROOT = _project_root()
import os
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE91061 黑色素瘤分子亚型分析
通过表达特征定义分子亚型，分析 ZP3 在不同亚型中的表达与治疗反应关系
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
import warnings
warnings.filterwarnings('ignore')

# 设置绘图风格
plt.style.use('default')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10

OUTPUT_DIR = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling")
DATA_DIR = os.path.join(ROOT, "output", "gse91061_validation")

# TCGA 黑色素瘤转录组亚型基因集
# 来源: TCGA Skin Cutaneous Melanoma (Cell 2015) 4 种转录组亚型
IMMUNE_GENES = [  # Immune subtype signature genes
    'CD2', 'CD3D', 'CD3E', 'CD3G', 'CD8A', 'GZMA', 'GZMB', 'PRF1',
    'CCL5', 'STAT1', 'CXCL9', 'CXCL10', 'IFNG', 'NKG7', 'TBX21', 'IL7R'
]
KERATIN_GENES = ['KRT14', 'KRT17', 'KRT5', 'KRT6A', 'KRT6B', 'KRT6C', 'KRT16', 'KRT1']  # Keratin subtype
MITF_GENES = ['MITF', 'MLANA', 'PMEL', 'TYRP1', 'TYR', 'DCT', 'SILV']  # MITF-high

# 免疫特征基因集
IFNG_SIGNATURE = ['IFNG', 'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'STAT1', 'HLA-DRA', 'CCL5', 'GZMB', 'PRF1']
CD8_SIGNATURE = ['CD8A', 'CD8B']
TCELL_SIGNATURE = ['CD3D', 'CD3E', 'CD2', 'CD28', 'TRAC', 'IL7R']

# 黑色素瘤关键免疫抑制基因
IMMUNOSUPPRESSIVE = {
    'TGFB1': 'TGFB1', 'IL10': 'IL10', 'CD274': 'CD274', 'PDCD1LG2': 'PDCD1LG2',
    'CTLA4': 'CTLA4', 'LAG3': 'LAG3', 'HAVCR2': 'HAVCR2', 'PDCD1': 'PDCD1',
    'TIGIT': 'TIGIT', 'CD163': 'CD163', 'MRC1': 'MRC1', 'FOXP3': 'FOXP3'
}

def main():
    print("=" * 60)
    print("GSE91061 黑色素瘤分子亚型分析")
    print("=" * 60)
    
    # 1. 加载表达数据
    print("\n1. 加载表达数据...")
    expr = pd.read_csv(f'{DATA_DIR}/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv', index_col=0)
    print(f"  基因数: {len(expr)}")
    print(f"  样本数: {len(expr.columns)}")
    
    # 2. 解析临床信息
    print("\n2. 解析临床信息...")
    sample_titles, gsm_ids, visits, responses = [], [], [], []
    
    with open(f'{DATA_DIR}/GSE91061_series_matrix.txt', 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    for line in lines:
        if line.startswith('!Sample_title'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                sample_titles = [p.strip('"') for p in parts[1:]]
        elif line.startswith('!Sample_geo_accession'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                gsm_ids = [p.strip('"') for p in parts[1:]]
    
    # 提取 visit 和 response
    for line in lines:
        if line.startswith('!Sample_characteristics_ch1'):
            parts = line.strip().split('\t')
            if len(parts) > 1 and 'visit' in parts[1].lower():
                visits = [p.split(': ')[-1].strip('"') for p in parts[1:] if ': ' in p]
            elif len(parts) > 1 and 'response' in parts[1].lower():
                responses = [p.split(': ')[-1].strip('"') for p in parts[1:] if ': ' in p]
    
    clinical = pd.DataFrame({
        'sample_title': sample_titles, 'gsm_id': gsm_ids,
        'visit': visits, 'response': responses
    })
    
    # 3. 对齐样本并提取 ZP3
    print("\n3. 对齐样本并提取 ZP3...")
    zp3_entrez = 7784  # 2026-08-17 审计修复: ZP3 = 7784 (原误用 8277=SP5)
    common_samples = [t for t in clinical['sample_title'] if t in expr.columns]
    print(f"  共同样本数: {len(common_samples)}")
    
    clinical_aligned = clinical[clinical['sample_title'].isin(common_samples)].set_index('sample_title')
    expr_aligned = expr[common_samples]
    
    # 4. 基因名映射（Entrez ID -> Gene Symbol）
    print("\n4. 基因名映射...")
    # 使用一个简单的映射表（基于已知常见基因）
    # 由于表达数据是 Entrez ID，需要映射到基因符号
    entrez_to_symbol = get_gene_mapping()
    # 处理重复映射：如果多个 Entrez ID 映射到同一 symbol，取表达量均值
    expr_symbol = expr_aligned.rename(index=entrez_to_symbol)
    if expr_symbol.index.duplicated().any():
        print(f"  检测到 {expr_symbol.index.duplicated().sum()} 个重复基因名，按均值聚合")
        expr_symbol = expr_symbol.groupby(expr_symbol.index).mean()
    
    # 5. 计算分子亚型
    print("\n5. 计算分子亚型...")
    subtype_scores = compute_subtype_scores(expr_symbol)
    clinical_aligned = clinical_aligned.join(subtype_scores)
    
    # 6. 计算免疫特征
    print("\n6. 计算免疫特征...")
    immune_features = compute_immune_features(expr_symbol)
    clinical_aligned = clinical_aligned.join(immune_features)
    
    # 7. 定义分子亚型
    print("\n7. 定义分子亚型...")
    clinical_aligned = assign_subtypes(clinical_aligned)
    
    # 8. 提取 ZP3
    zp3 = expr_aligned.loc[zp3_entrez]
    clinical_aligned['ZP3'] = zp3[clinical_aligned.index]
    
    # 9. 分析
    print("\n8. 分子亚型与 ZP3 表达...")
    analyze_subtypes(clinical_aligned)
    
    # 10. 分子亚型与治疗反应
    print("\n9. 分子亚型与治疗反应...")
    analyze_response(clinical_aligned)
    
    # 11. 多因素分析
    print("\n10. 多因素分析...")
    multivariate_analysis(clinical_aligned)
    
    # 12. 可视化
    print("\n11. 生成可视化...")
    visualize(clinical_aligned)
    
    # 13. 保存结果
    print("\n12. 保存结果...")
    clinical_aligned.to_csv(f'{OUTPUT_DIR}/gse91061_molecular_subtype_clinical.csv', encoding='utf-8-sig')
    print(f"  已保存: gse91061_molecular_subtype_clinical.csv")
    
    # 创建报告
    create_report(clinical_aligned)
    
    print("\n" + "=" * 60)
    print("分析完成")
    print("=" * 60)

def get_gene_mapping():
    """获取 Entrez ID -> Gene Symbol 映射"""
    # 关键基因的 Entrez ID 映射
    mapping = {
        7784: 'ZP3',  # ZP3 (2026-08-17 审计修复: 原误用 8277=SP5)
        54209: 'TREM2',
        968: 'CD68',
        9332: 'CD163',
        4360: 'MRC1',
        929: 'CD14',
        4069: 'LYZ',
        1436: 'CSF1R',
        3684: 'ITGAM',
        3586: 'IL10',
        7124: 'TNF',
        3001: 'GZMA',
        5551: 'PRF1',
        3458: 'IFNG',
        5133: 'PDCD1',
        29126: 'CD274',
        1493: 'CTLA4',
        3902: 'LAG3',
        50943: 'FOXP3',
        925: 'CD8A',
        926: 'CD8B',
        920: 'CD4',
        925: 'CD8A',
        5788: 'PTPRC',
        199: 'AIF1',
        712: 'C1QA',
        7305: 'TYROBP',
        348: 'APOE',
        2670: 'GFAP',
        10215: 'OLIG2',
        6663: 'SOX10',
        7040: 'TGFB1',
        80380: 'PDCD1LG2',
        201633: 'TIGIT',
        3903: 'HAVCR2',
        925: 'CD8A',
        7040: 'TGFB1',
        3569: 'IL6',
        3561: 'IL2',
        60: 'ACTB',
        # 免疫亚型基因
        914: 'CD2', 923: 'CD3D', 916: 'CD3E', 917: 'CD3G',
        6352: 'CCL5', 6772: 'STAT1', 4283: 'CXCL9', 3627: 'CXCL10',
        4818: 'NKG7', 30009: 'TBX21', 3605: 'IL7R', 5199: 'CFLAR',
        # 角蛋白亚型
        3861: 'KRT14', 3872: 'KRT17', 3852: 'KRT5', 3853: 'KRT6A',
        3854: 'KRT6B', 3855: 'KRT6C', 3868: 'KRT16', 3856: 'KRT1',
        # MITF 亚型
        4286: 'MITF', 2315: 'MLANA', 6490: 'PMEL', 7306: 'TYRP1',
        7299: 'TYR', 1631: 'DCT', 4128: 'MAFA',
        # 其他免疫
        3575: 'IL7R', 3605: 'IL7R', 3117: 'HLA-DQA1', 3119: 'HLA-DQB1',
        3122: 'HLA-DRA', 3123: 'HLA-DRB1', 3107: 'HLA-C', 3105: 'HLA-A',
        3106: 'HLA-B', 6925: 'TCEB2', 3162: 'HMOX1',
        # IDO1
        3620: 'IDO1',
        # CD28
        940: 'CD28',
        # TRAC
        6955: 'TRAC',
    }
    return mapping

def compute_subtype_scores(expr_symbol):
    """计算亚型评分"""
    scores = {}
    
    # 免疫亚型评分（免疫基因平均 z-score）
    for name, genes in [('Immune_score', IMMUNE_GENES), ('Keratin_score', KERATIN_GENES), ('MITF_score', MITF_GENES)]:
        available = [g for g in genes if g in expr_symbol.index]
        if available:
            # 使用 log1p 转换后标准化
            sub_expr = np.log1p(expr_symbol.loc[available])
            # z-score 标准化（按基因，axis=0 显式对齐行索引）
            z_scores = sub_expr.sub(sub_expr.mean(axis=1), axis=0).div(sub_expr.std(axis=1) + 1e-10, axis=0)
            scores[name] = z_scores.mean(axis=0)
        else:
            scores[name] = pd.Series(0, index=expr_symbol.columns)
    
    return pd.DataFrame(scores)

def compute_immune_features(expr_symbol):
    """计算免疫特征"""
    features = {}
    
    def calc_score(genes):
        available = [g for g in genes if g in expr_symbol.index]
        if available:
            sub_expr = np.log1p(expr_symbol.loc[available])
            return sub_expr.mean(axis=0)
        return pd.Series(0, index=expr_symbol.columns)
    
    features['IFN_gamma_score'] = calc_score(IFNG_SIGNATURE)
    features['CD8_score'] = calc_score(CD8_SIGNATURE)
    features['Tcell_score'] = calc_score(TCELL_SIGNATURE)
    
    # 免疫抑制基因
    for gene, symbol in IMMUNOSUPPRESSIVE.items():
        if symbol in expr_symbol.index:
            features[f'expr_{symbol}'] = np.log1p(expr_symbol.loc[symbol])
        else:
            features[f'expr_{symbol}'] = 0
    
    # 综合免疫抑制评分
    supp_genes = [f'expr_{s}' for s in IMMUNOSUPPRESSIVE.values() if s in expr_symbol.index]
    if supp_genes:
        features['Immunosuppression_score'] = sum(features[g] for g in supp_genes) / len(supp_genes)
    
    return pd.DataFrame(features)

def assign_subtypes(df):
    """基于评分分配分子亚型"""
    # TCGA 黑色素瘤转录组亚型：
    # 1. Immune (免疫富集)
    # 2. Keratin (角质蛋白)
    # 3. MITF-low
    # 4. MITF-high
    subtypes = []
    
    immune_threshold = df['Immune_score'].median()
    keratin_threshold = df['Keratin_score'].median()
    mitf_threshold = df['MITF_score'].median()
    
    for idx, row in df.iterrows():
        if row['Immune_score'] > immune_threshold:
            subtypes.append('Immune')
        elif row['Keratin_score'] > keratin_threshold:
            subtypes.append('Keratin')
        elif row['MITF_score'] > mitf_threshold:
            subtypes.append('MITF-high')
        else:
            subtypes.append('MITF-low')
    
    df['Molecular_subtype'] = subtypes
    return df

def analyze_subtypes(df):
    """分析分子亚型与 ZP3 表达"""
    print("\n  分子亚型分布:")
    print(df['Molecular_subtype'].value_counts())
    
    print("\n  各亚型 ZP3 表达:")
    subtype_groups = df.groupby('Molecular_subtype')
    for name, group in subtype_groups:
        print(f"    {name}: n={len(group)}, ZP3 mean={group['ZP3'].mean():.3f}, "
              f"median={group['ZP3'].median():.3f}")
    
    # Kruskal-Wallis 检验
    if df['Molecular_subtype'].nunique() >= 3:
        groups = [group['ZP3'].values for name, group in subtype_groups]
        stat, p = stats.kruskal(*groups)
        print(f"\n  Kruskal-Wallis 检验: H={stat:.3f}, p={p:.4f}")
    
    # 免疫 vs 非免疫
    immune = df[df['Molecular_subtype'] == 'Immune']['ZP3']
    non_immune = df[df['Molecular_subtype'] != 'Immune']['ZP3']
    if len(immune) > 5 and len(non_immune) > 5:
        stat, p = stats.mannwhitneyu(immune, non_immune, alternative='two-sided')
        print(f"\n  Immune vs Non-immune: n={len(immune)}/{len(non_immune)}, "
              f"mean={immune.mean():.3f}/{non_immune.mean():.3f}, p={p:.4f}")

def analyze_response(df):
    """分析分子亚型与治疗反应"""
    # 只分析 pre-treatment 样本（治疗前基线）
    pre_df = df[df['visit'] == 'Pre'].copy()
    print(f"\n  Pre-treatment 样本数: {len(pre_df)}")
    
    # 二元反应
    response_map = {'PRCR': 'Response', 'PR': 'Response', 'CR': 'Response',
                    'SD': 'Non-response', 'PD': 'Non-response'}
    pre_df['binary_response'] = pre_df['response'].map(response_map)
    valid = pre_df[pre_df['binary_response'].notna()]
    print(f"  有明确反应的样本数: {len(valid)}")
    
    if len(valid) < 20:
        print("  样本太少，跳过反应分析")
        return
    
    # 整体 ZP3 与反应
    resp = valid[valid['binary_response'] == 'Response']
    non_resp = valid[valid['binary_response'] == 'Non-response']
    if len(resp) > 3 and len(non_resp) > 3:
        stat, p = stats.mannwhitneyu(resp['ZP3'], non_resp['ZP3'], alternative='two-sided')
        n1, n2 = len(resp), len(non_resp)
        r = 2 * stat / (n1 * n2) - 1  # 2026-08-17 审计修复: positive = responders higher
        print(f"\n  整体: Response n={len(resp)}, Non-response n={len(non_resp)}")
        print(f"    ZP3 mean: {resp['ZP3'].mean():.3f} vs {non_resp['ZP3'].mean():.3f}")
        print(f"    Mann-Whitney p={p:.4f}, r={r:.3f}")
    
    # 按分子亚型分层
    print("\n  按分子亚型分层:")
    for subtype in valid['Molecular_subtype'].unique():
        sub = valid[valid['Molecular_subtype'] == subtype]
        sub_resp = sub[sub['binary_response'] == 'Response']
        sub_non = sub[sub['binary_response'] == 'Non-response']
        if len(sub_resp) >= 2 and len(sub_non) >= 2:
            stat, p = stats.mannwhitneyu(sub_resp['ZP3'], sub_non['ZP3'], alternative='two-sided')
            print(f"    {subtype}: Resp n={len(sub_resp)}, Non-resp n={len(sub_non)}, "
                  f"ZP3 {sub_resp['ZP3'].mean():.2f} vs {sub_non['ZP3'].mean():.2f}, p={p:.4f}")
    
    # 免疫特征分层（IFN-γ 高/低）
    print("\n  按 IFN-γ 特征分层:")
    if 'IFN_gamma_score' in valid.columns:
        ifn_median = valid['IFN_gamma_score'].median()
        ifn_high = valid[valid['IFN_gamma_score'] > ifn_median]
        ifn_low = valid[valid['IFN_gamma_score'] <= ifn_median]
        
        if len(ifn_high) > 5 and len(ifn_low) > 5:
            stat, p = stats.mannwhitneyu(ifn_high['ZP3'], ifn_low['ZP3'], alternative='two-sided')
            print(f"    IFN-γ high (n={len(ifn_high)}): ZP3 mean={ifn_high['ZP3'].mean():.3f}")
            print(f"    IFN-γ low (n={len(ifn_low)}): ZP3 mean={ifn_low['ZP3'].mean():.3f}")
            print(f"    Mann-Whitney p={p:.4f}")
            
            # 亚组内反应分析
            for label, sub in [('IFN-γ high', ifn_high), ('IFN-γ low', ifn_low)]:
                sub_resp = sub[sub['binary_response'] == 'Response']
                sub_non = sub[sub['binary_response'] == 'Non-response']
                if len(sub_resp) >= 2 and len(sub_non) >= 2:
                    stat, p = stats.mannwhitneyu(sub_resp['ZP3'], sub_non['ZP3'], alternative='two-sided')
                    print(f"      {label}: Resp n={len(sub_resp)}, Non-resp n={len(sub_non)}, p={p:.4f}")
                else:
                    print(f"      {label}: 样本不足 (Resp={len(sub_resp)}, Non-resp={len(sub_non)})")

def multivariate_analysis(df):
    """多因素 Logistic 回归分析"""
    from statsmodels.formula.api import logit
    
    pre_df = df[df['visit'] == 'Pre'].copy()
    response_map = {'PRCR': 1, 'PR': 1, 'CR': 1, 'SD': 0, 'PD': 0}
    pre_df['binary_response'] = pre_df['response'].map(response_map)
    valid = pre_df[pre_df['binary_response'].notna()].copy()
    
    if len(valid) < 15:
        print("  样本太少，跳过多因素分析")
        return
    
    # 标准化 ZP3 和免疫特征
    valid['ZP3_scaled'] = (valid['ZP3'] - valid['ZP3'].mean()) / valid['ZP3'].std()
    if 'IFN_gamma_score' in valid.columns:
        valid['IFN_scaled'] = (valid['IFN_gamma_score'] - valid['IFN_gamma_score'].mean()) / valid['IFN_gamma_score'].std()
    if 'Immunosuppression_score' in valid.columns:
        valid['Supp_scaled'] = (valid['Immunosuppression_score'] - valid['Immunosuppression_score'].mean()) / valid['Immunosuppression_score'].std()
    
    # 模型 1: ZP3 单独
    print("\n  模型 1: response ~ ZP3")
    try:
        model1 = logit('binary_response ~ ZP3_scaled', data=valid).fit(disp=False)
        print(f"    ZP3: OR={np.exp(model1.params['ZP3_scaled']):.3f}, "
              f"p={model1.pvalues['ZP3_scaled']:.4f}, "
              f"95%CI=[{np.exp(model1.conf_int().loc['ZP3_scaled', 0]):.3f}, "
              f"{np.exp(model1.conf_int().loc['ZP3_scaled', 1]):.3f}]")
    except Exception as e:
        print(f"    模型 1 失败: {e}")
    
    # 模型 2: ZP3 + IFN-γ
    if 'IFN_scaled' in valid.columns:
        print("\n  模型 2: response ~ ZP3 + IFN-γ")
        try:
            model2 = logit('binary_response ~ ZP3_scaled + IFN_scaled', data=valid).fit(disp=False)
            for var in ['ZP3_scaled', 'IFN_scaled']:
                print(f"    {var}: OR={np.exp(model2.params[var]):.3f}, p={model2.pvalues[var]:.4f}")
        except Exception as e:
            print(f"    模型 2 失败: {e}")
    
    # 模型 3: ZP3 + 免疫抑制评分
    if 'Supp_scaled' in valid.columns:
        print("\n  模型 3: response ~ ZP3 + Immunosuppression")
        try:
            model3 = logit('binary_response ~ ZP3_scaled + Supp_scaled', data=valid).fit(disp=False)
            for var in ['ZP3_scaled', 'Supp_scaled']:
                print(f"    {var}: OR={np.exp(model3.params[var]):.3f}, p={model3.pvalues[var]:.4f}")
        except Exception as e:
            print(f"    模型 3 失败: {e}")

def visualize(df):
    """可视化"""
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('GSE91061 Molecular Subtype Analysis', fontsize=14, fontweight='bold')
    
    # 1. 分子亚型分布
    ax = axes[0, 0]
    df['Molecular_subtype'].value_counts().plot(kind='bar', ax=ax, color=['#3498db', '#e74c3c', '#2ecc71', '#f39c12'])
    ax.set_title('Molecular Subtype Distribution')
    ax.set_ylabel('Count')
    ax.tick_params(axis='x', rotation=0)
    
    # 2. ZP3 按亚型
    ax = axes[0, 1]
    sns.boxplot(data=df, x='Molecular_subtype', y='ZP3', ax=ax, palette='husl')
    ax.set_title('ZP3 by Molecular Subtype')
    ax.tick_params(axis='x', rotation=45)
    
    # 3. ZP3 与 IFN-γ 散点
    ax = axes[0, 2]
    if 'IFN_gamma_score' in df.columns:
        ax.scatter(df['IFN_gamma_score'], np.log1p(df['ZP3']), alpha=0.5, s=20)
        # 添加回归线
        mask = df['IFN_gamma_score'].notna() & df['ZP3'].notna()
        if mask.sum() > 5:
            slope, intercept, r, p, se = stats.linregress(
                df.loc[mask, 'IFN_gamma_score'], np.log1p(df.loc[mask, 'ZP3']))
            x = np.linspace(df['IFN_gamma_score'].min(), df['IFN_gamma_score'].max(), 100)
            ax.plot(x, slope * x + intercept, 'r-', alpha=0.7)
            ax.text(0.05, 0.9, f'r={r:.3f}, p={p:.4f}', transform=ax.transAxes)
    ax.set_xlabel('IFN-γ Score')
    ax.set_ylabel('log1p(ZP3)')
    ax.set_title('ZP3 vs IFN-γ Signature')
    
    # 4. 治疗反应按亚型
    ax = axes[1, 0]
    pre_df = df[df['visit'] == 'Pre'].copy()
    response_map = {'PRCR': 'Response', 'PR': 'Response', 'CR': 'Response',
                    'SD': 'Non-response', 'PD': 'Non-response'}
    pre_df['binary_response'] = pre_df['response'].map(response_map)
    valid = pre_df[pre_df['binary_response'].notna()]
    if len(valid) > 5:
        crosstab = pd.crosstab(valid['Molecular_subtype'], valid['binary_response'])
        crosstab = crosstab.reindex(columns=['Response', 'Non-response'], fill_value=0)
        crosstab.plot(kind='bar', stacked=True, ax=ax, color=['#2ecc71', '#e74c3c'])
    ax.set_title('Response by Molecular Subtype (Pre-treatment)')
    ax.set_ylabel('Count')
    ax.tick_params(axis='x', rotation=45)
    
    # 5. 各亚型响应率
    ax = axes[1, 1]
    if len(valid) > 5:
        subtype_rates = {}
        for subtype in valid['Molecular_subtype'].unique():
            sub = valid[valid['Molecular_subtype'] == subtype]
            rate = (sub['binary_response'] == 'Response').mean() * 100
            subtype_rates[subtype] = rate
        pd.Series(subtype_rates).sort_values().plot(kind='barh', ax=ax, color='#3498db')
    ax.set_xlabel('Response Rate (%)')
    ax.set_title('Response Rate by Subtype')
    
    # 6. ZP3 与免疫抑制基因相关性热图
    ax = axes[1, 2]
    supp_cols = [c for c in df.columns if c.startswith('expr_')]
    if supp_cols and 'ZP3' in df.columns:
        corr_genes = supp_cols + ['ZP3']
        corr_data = df[corr_genes].corr(method='spearman')
        im = ax.imshow(corr_data, cmap='RdBu_r', vmin=-1, vmax=1)
        ax.set_xticks(range(len(corr_genes)))
        ax.set_yticks(range(len(corr_genes)))
        ax.set_xticklabels([g.replace('expr_', '') for g in corr_genes], rotation=90)
        ax.set_yticklabels([g.replace('expr_', '') for g in corr_genes])
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_title('ZP3-Immune Gene Correlation')
    
    plt.tight_layout()
    fig_path = f'{OUTPUT_DIR}/fig_gse91061_molecular_subtype.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  已保存可视化: {fig_path}")

def create_report(df):
    """创建报告"""
    report = f"""# GSE91061 黑色素瘤分子亚型分析报告

## 分析信息
- 队列: GSE91061 (n=109, 黑色素瘤抗 PD-1 治疗)
- 治疗: Nivolumab (anti-PD-1)
- 数据版本: hg19, FPKM normalized

## 1. 分子亚型定义
基于 TCGA SKCM 转录组亚型分类器，将样本分为 4 种分子亚型:
- **Immune**: 免疫基因富集（CD3, CD8, CCL5, STAT1, CXCL9/10 等）
- **Keratin**: 角质蛋白基因富集（KRT5, KRT14, KRT17 等）
- **MITF-high**: MITF/色素基因高表达
- **MITF-low**: MITF/色素基因低表达

## 2. 分子亚型分布
"""
    report += "| 亚型 | 样本数 | 占比 |\n|---|---|---|\n"
    for subtype, count in df['Molecular_subtype'].value_counts().items():
        report += f"| {subtype} | {count} | {100*count/len(df):.1f}% |\n"
    
    report += f"""
## 3. ZP3 表达按亚型
"""
    for subtype, group in df.groupby('Molecular_subtype'):
        report += f"- **{subtype}** (n={len(group)}): ZP3 mean={group['ZP3'].mean():.3f}, median={group['ZP3'].median():.3f}\n"
    
    report += f"""
## 4. 治疗反应分析
- 只分析 pre-treatment 样本
- Response: PRCR (Partial/Complete Response)
- Non-response: SD (Stable Disease) + PD (Progressive Disease)

## 5. 免疫特征
- IFN-γ 特征、CD8 特征、T 细胞特征、免疫抑制评分已计算
- ZP3 与免疫抑制基因相关性已计算

## 6. 多因素分析
- 模型 1: response ~ ZP3
- 模型 2: response ~ ZP3 + IFN-γ
- 模型 3: response ~ ZP3 + Immunosuppression

## 7. 关键发现
- ZP3 在分子亚型间的表达差异
- ZP3 与治疗反应的关系是否依赖分子亚型
- ZP3 是否独立于已知免疫标志物

## 8. 局限
- GSE91061 不提供 BRAF/NRAS 突变和 TMB 数据
- 分子亚型基于表达特征推断，非 TCGA 官方分类
- 样本量有限，亚型内部分析统计功效不足

## 9. 输出文件
- `gse91061_molecular_subtype_clinical.csv` - 完整临床+分子数据
- `fig_gse91061_molecular_subtype.png` - 可视化
"""
    
    report_path = f'{OUTPUT_DIR}/gse91061_molecular_subtype_report.md'
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"  已创建报告: {report_path}")

if __name__ == "__main__":
    main()