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
GSE91061 melanoma molecular subtype analysis
Define molecular subtypes via expression signatures, analyze ZP3 expression and treatment response across subtypes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
import warnings
warnings.filterwarnings('ignore')

# Set plot style
plt.style.use('default')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10

OUTPUT_DIR = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling")
DATA_DIR = os.path.join(ROOT, "output", "gse91061_validation")

# TCGA melanoma transcriptome subtype gene sets
# Source: TCGA Skin Cutaneous Melanoma (Cell 2015) 4 transcriptome subtypes
IMMUNE_GENES = [  # Immune subtype signature genes
    'CD2', 'CD3D', 'CD3E', 'CD3G', 'CD8A', 'GZMA', 'GZMB', 'PRF1',
    'CCL5', 'STAT1', 'CXCL9', 'CXCL10', 'IFNG', 'NKG7', 'TBX21', 'IL7R'
]
KERATIN_GENES = ['KRT14', 'KRT17', 'KRT5', 'KRT6A', 'KRT6B', 'KRT6C', 'KRT16', 'KRT1']  # Keratin subtype
MITF_GENES = ['MITF', 'MLANA', 'PMEL', 'TYRP1', 'TYR', 'DCT', 'SILV']  # MITF-high

# Immune signature gene sets
IFNG_SIGNATURE = ['IFNG', 'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'STAT1', 'HLA-DRA', 'CCL5', 'GZMB', 'PRF1']
CD8_SIGNATURE = ['CD8A', 'CD8B']
TCELL_SIGNATURE = ['CD3D', 'CD3E', 'CD2', 'CD28', 'TRAC', 'IL7R']

# Key immunosuppressive genes in melanoma
IMMUNOSUPPRESSIVE = {
    'TGFB1': 'TGFB1', 'IL10': 'IL10', 'CD274': 'CD274', 'PDCD1LG2': 'PDCD1LG2',
    'CTLA4': 'CTLA4', 'LAG3': 'LAG3', 'HAVCR2': 'HAVCR2', 'PDCD1': 'PDCD1',
    'TIGIT': 'TIGIT', 'CD163': 'CD163', 'MRC1': 'MRC1', 'FOXP3': 'FOXP3'
}

def main():
    print("=" * 60)
    print("GSE91061 melanoma molecular subtype analysis")
    print("=" * 60)
    
    # 1. Load expression data
    print("\n1. Loading expression data...")
    expr = pd.read_csv(f'{DATA_DIR}/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv', index_col=0)
    print(f"  Number of genes: {len(expr)}")
    print(f"  Number of samples: {len(expr.columns)}")
    
    # 2. Parse clinical information
    print("\n2. Parsing clinical information...")
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
    
    # Extract visit and response
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
    
    # 3. Align samples and extract ZP3
    print("\n3. Aligning samples and extracting ZP3...")
    zp3_entrez = 7784  # 2026-08-17 audit fix: ZP3 = 7784 (previously incorrectly used 8277=SP5)
    common_samples = [t for t in clinical['sample_title'] if t in expr.columns]
    print(f"  Common sample count: {len(common_samples)}")
    
    clinical_aligned = clinical[clinical['sample_title'].isin(common_samples)].set_index('sample_title')
    expr_aligned = expr[common_samples]
    
    # 4. Gene name mapping (Entrez ID -> Gene Symbol)
    print("\n4. Mapping gene names...")
    # Use a simple mapping table (based on known common genes)
    # Since expression data uses Entrez IDs, map them to gene symbols
    entrez_to_symbol = get_gene_mapping()
    # Handle duplicate mappings: if multiple Entrez IDs map to same symbol, take mean expression
    expr_symbol = expr_aligned.rename(index=entrez_to_symbol)
    if expr_symbol.index.duplicated().any():
        print(f"  Detected {expr_symbol.index.duplicated().sum()} duplicated gene names, aggregating by mean")
        expr_symbol = expr_symbol.groupby(expr_symbol.index).mean()
    
    # 5. Compute molecular subtypes
    print("\n5. Computing molecular subtypes...")
    subtype_scores = compute_subtype_scores(expr_symbol)
    clinical_aligned = clinical_aligned.join(subtype_scores)
    
    # 6. Compute immune features
    print("\n6. Computing immune features...")
    immune_features = compute_immune_features(expr_symbol)
    clinical_aligned = clinical_aligned.join(immune_features)
    
    # 7. Define molecular subtypes
    print("\n7. Defining molecular subtypes...")
    clinical_aligned = assign_subtypes(clinical_aligned)
    
    # 8. Extract ZP3
    zp3 = expr_aligned.loc[zp3_entrez]
    clinical_aligned['ZP3'] = zp3[clinical_aligned.index]
    
    # 9. Analysis
    print("\n8. Molecular subtypes and ZP3 expression...")
    analyze_subtypes(clinical_aligned)
    
    # 10. Molecular subtypes and treatment response
    print("\n9. Molecular subtypes and treatment response...")
    analyze_response(clinical_aligned)
    
    # 11. Multivariate analysis
    print("\n10. Multivariate analysis...")
    multivariate_analysis(clinical_aligned)
    
    # 12. Visualization
    print("\n11. Generating visualization...")
    visualize(clinical_aligned)
    
    # 13. Save results
    print("\n12. Saving results...")
    clinical_aligned.to_csv(f'{OUTPUT_DIR}/gse91061_molecular_subtype_clinical.csv', encoding='utf-8-sig')
    print(f"  Saved: gse91061_molecular_subtype_clinical.csv")
    
    # Create report
    create_report(clinical_aligned)
    
    print("\n" + "=" * 60)
    print("Analysis complete")
    print("=" * 60)

def get_gene_mapping():
    """Get Entrez ID -> Gene Symbol mapping"""
    # Entrez ID mapping for key genes
    mapping = {
        7784: 'ZP3',  # ZP3 (2026-08-17 audit fix: previously incorrectly used 8277=SP5)
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
        # Immune subtype genes
        914: 'CD2', 923: 'CD3D', 916: 'CD3E', 917: 'CD3G',
        6352: 'CCL5', 6772: 'STAT1', 4283: 'CXCL9', 3627: 'CXCL10',
        4818: 'NKG7', 30009: 'TBX21', 3605: 'IL7R', 5199: 'CFLAR',
        # Keratin subtype
        3861: 'KRT14', 3872: 'KRT17', 3852: 'KRT5', 3853: 'KRT6A',
        3854: 'KRT6B', 3855: 'KRT6C', 3868: 'KRT16', 3856: 'KRT1',
        # MITF subtype
        4286: 'MITF', 2315: 'MLANA', 6490: 'PMEL', 7306: 'TYRP1',
        7299: 'TYR', 1631: 'DCT', 4128: 'MAFA',
        # Other immune
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
    """Calculate subtype scores"""
    scores = {}
    
    # Immune subtype score (mean z-score of immune genes)
    for name, genes in [('Immune_score', IMMUNE_GENES), ('Keratin_score', KERATIN_GENES), ('MITF_score', MITF_GENES)]:
        available = [g for g in genes if g in expr_symbol.index]
        if available:
            # Standardize after log1p transformation
            sub_expr = np.log1p(expr_symbol.loc[available])
            # z-score standardization (by gene, axis=0 explicitly aligns row index)
            z_scores = sub_expr.sub(sub_expr.mean(axis=1), axis=0).div(sub_expr.std(axis=1) + 1e-10, axis=0)
            scores[name] = z_scores.mean(axis=0)
        else:
            scores[name] = pd.Series(0, index=expr_symbol.columns)
    
    return pd.DataFrame(scores)

def compute_immune_features(expr_symbol):
    """Calculate immune features"""
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
    
    # Immunosuppressive genes
    for gene, symbol in IMMUNOSUPPRESSIVE.items():
        if symbol in expr_symbol.index:
            features[f'expr_{symbol}'] = np.log1p(expr_symbol.loc[symbol])
        else:
            features[f'expr_{symbol}'] = 0
    
    # Combined immunosuppression score
    supp_genes = [f'expr_{s}' for s in IMMUNOSUPPRESSIVE.values() if s in expr_symbol.index]
    if supp_genes:
        features['Immunosuppression_score'] = sum(features[g] for g in supp_genes) / len(supp_genes)
    
    return pd.DataFrame(features)

def assign_subtypes(df):
    """Assign molecular subtypes based on scores"""
    # TCGA melanoma transcriptomic subtypes:
    # 1. Immune (immune enriched)
    # 2. Keratin (keratin)
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
    """Analyze molecular subtypes and ZP3 expression"""
    print("\n  Molecular subtype distribution:")
    print(df['Molecular_subtype'].value_counts())
    
    print("\n  ZP3 expression by subtype:")
    subtype_groups = df.groupby('Molecular_subtype')
    for name, group in subtype_groups:
        print(f"    {name}: n={len(group)}, ZP3 mean={group['ZP3'].mean():.3f}, "
              f"median={group['ZP3'].median():.3f}")
    
    # Kruskal-Wallis test
    if df['Molecular_subtype'].nunique() >= 3:
        groups = [group['ZP3'].values for name, group in subtype_groups]
        stat, p = stats.kruskal(*groups)
        print(f"\n  Kruskal-Wallis test: H={stat:.3f}, p={p:.4f}")
    
    # Immune vs non-immune
    immune = df[df['Molecular_subtype'] == 'Immune']['ZP3']
    non_immune = df[df['Molecular_subtype'] != 'Immune']['ZP3']
    if len(immune) > 5 and len(non_immune) > 5:
        stat, p = stats.mannwhitneyu(immune, non_immune, alternative='two-sided')
        print(f"\n  Immune vs Non-immune: n={len(immune)}/{len(non_immune)}, "
              f"mean={immune.mean():.3f}/{non_immune.mean():.3f}, p={p:.4f}")

def analyze_response(df):
    """Analyze molecular subtypes and treatment response"""
    # Analyze only pre-treatment samples (baseline before treatment)
    pre_df = df[df['visit'] == 'Pre'].copy()
    print(f"\n  Pre-treatment sample count: {len(pre_df)}")
    
    # Binary response
    response_map = {'PRCR': 'Response', 'PR': 'Response', 'CR': 'Response',
                    'SD': 'Non-response', 'PD': 'Non-response'}
    pre_df['binary_response'] = pre_df['response'].map(response_map)
    valid = pre_df[pre_df['binary_response'].notna()]
    print(f"  Number of samples with definite response: {len(valid)}")
    
    if len(valid) < 20:
        print("  Too few samples, skipping response analysis")
        return
    
    # Overall ZP3 vs response
    resp = valid[valid['binary_response'] == 'Response']
    non_resp = valid[valid['binary_response'] == 'Non-response']
    if len(resp) > 3 and len(non_resp) > 3:
        stat, p = stats.mannwhitneyu(resp['ZP3'], non_resp['ZP3'], alternative='two-sided')
        n1, n2 = len(resp), len(non_resp)
        r = 2 * stat / (n1 * n2) - 1  # 2026-08-17 audit fix: positive = responders higher
        print(f"\n  Overall: Response n={len(resp)}, Non-response n={len(non_resp)}")
        print(f"    ZP3 mean: {resp['ZP3'].mean():.3f} vs {non_resp['ZP3'].mean():.3f}")
        print(f"    Mann-Whitney p={p:.4f}, r={r:.3f}")
    
    # Stratify by molecular subtype
    print("\n  Stratified by molecular subtype:")
    for subtype in valid['Molecular_subtype'].unique():
        sub = valid[valid['Molecular_subtype'] == subtype]
        sub_resp = sub[sub['binary_response'] == 'Response']
        sub_non = sub[sub['binary_response'] == 'Non-response']
        if len(sub_resp) >= 2 and len(sub_non) >= 2:
            stat, p = stats.mannwhitneyu(sub_resp['ZP3'], sub_non['ZP3'], alternative='two-sided')
            print(f"    {subtype}: Resp n={len(sub_resp)}, Non-resp n={len(sub_non)}, "
                  f"ZP3 {sub_resp['ZP3'].mean():.2f} vs {sub_non['ZP3'].mean():.2f}, p={p:.4f}")
    
    # Immune signature stratification (IFN-γ high/low)
    print("\n  Stratified by IFN-γ signature:")
    if 'IFN_gamma_score' in valid.columns:
        ifn_median = valid['IFN_gamma_score'].median()
        ifn_high = valid[valid['IFN_gamma_score'] > ifn_median]
        ifn_low = valid[valid['IFN_gamma_score'] <= ifn_median]
        
        if len(ifn_high) > 5 and len(ifn_low) > 5:
            stat, p = stats.mannwhitneyu(ifn_high['ZP3'], ifn_low['ZP3'], alternative='two-sided')
            print(f"    IFN-γ high (n={len(ifn_high)}): ZP3 mean={ifn_high['ZP3'].mean():.3f}")
            print(f"    IFN-γ low (n={len(ifn_low)}): ZP3 mean={ifn_low['ZP3'].mean():.3f}")
            print(f"    Mann-Whitney p={p:.4f}")
            
            # Response analysis within subgroup
            for label, sub in [('IFN-γ high', ifn_high), ('IFN-γ low', ifn_low)]:
                sub_resp = sub[sub['binary_response'] == 'Response']
                sub_non = sub[sub['binary_response'] == 'Non-response']
                if len(sub_resp) >= 2 and len(sub_non) >= 2:
                    stat, p = stats.mannwhitneyu(sub_resp['ZP3'], sub_non['ZP3'], alternative='two-sided')
                    print(f"      {label}: Resp n={len(sub_resp)}, Non-resp n={len(sub_non)}, p={p:.4f}")
                else:
                    print(f"      {label}: insufficient samples (Resp={len(sub_resp)}, Non-resp={len(sub_non)})")

def multivariate_analysis(df):
    """Multivariate logistic regression analysis"""
    from statsmodels.formula.api import logit
    
    pre_df = df[df['visit'] == 'Pre'].copy()
    response_map = {'PRCR': 1, 'PR': 1, 'CR': 1, 'SD': 0, 'PD': 0}
    pre_df['binary_response'] = pre_df['response'].map(response_map)
    valid = pre_df[pre_df['binary_response'].notna()].copy()
    
    if len(valid) < 15:
        print("  Too few samples, skipping multivariate analysis")
        return
    
    # Standardize ZP3 and immune features
    valid['ZP3_scaled'] = (valid['ZP3'] - valid['ZP3'].mean()) / valid['ZP3'].std()
    if 'IFN_gamma_score' in valid.columns:
        valid['IFN_scaled'] = (valid['IFN_gamma_score'] - valid['IFN_gamma_score'].mean()) / valid['IFN_gamma_score'].std()
    if 'Immunosuppression_score' in valid.columns:
        valid['Supp_scaled'] = (valid['Immunosuppression_score'] - valid['Immunosuppression_score'].mean()) / valid['Immunosuppression_score'].std()
    
    # Model 1: ZP3 alone
    print("\n  Model 1: response ~ ZP3")
    try:
        model1 = logit('binary_response ~ ZP3_scaled', data=valid).fit(disp=False)
        print(f"    ZP3: OR={np.exp(model1.params['ZP3_scaled']):.3f}, "
              f"p={model1.pvalues['ZP3_scaled']:.4f}, "
              f"95%CI=[{np.exp(model1.conf_int().loc['ZP3_scaled', 0]):.3f}, "
              f"{np.exp(model1.conf_int().loc['ZP3_scaled', 1]):.3f}]")
    except Exception as e:
        print(f"    Model 1 failed: {e}")
    
    # Model 2: ZP3 + IFN-γ
    if 'IFN_scaled' in valid.columns:
        print("\n  Model 2: response ~ ZP3 + IFN-γ")
        try:
            model2 = logit('binary_response ~ ZP3_scaled + IFN_scaled', data=valid).fit(disp=False)
            for var in ['ZP3_scaled', 'IFN_scaled']:
                print(f"    {var}: OR={np.exp(model2.params[var]):.3f}, p={model2.pvalues[var]:.4f}")
        except Exception as e:
            print(f"    Model 2 failed: {e}")
    
    # Model 3: ZP3 + immunosuppression score
    if 'Supp_scaled' in valid.columns:
        print("\n  Model 3: response ~ ZP3 + Immunosuppression")
        try:
            model3 = logit('binary_response ~ ZP3_scaled + Supp_scaled', data=valid).fit(disp=False)
            for var in ['ZP3_scaled', 'Supp_scaled']:
                print(f"    {var}: OR={np.exp(model3.params[var]):.3f}, p={model3.pvalues[var]:.4f}")
        except Exception as e:
            print(f"    Model 3 failed: {e}")

def visualize(df):
    """Visualization"""
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('GSE91061 Molecular Subtype Analysis', fontsize=14, fontweight='bold')
    
    # 1. Molecular subtype distribution
    ax = axes[0, 0]
    df['Molecular_subtype'].value_counts().plot(kind='bar', ax=ax, color=['#3498db', '#e74c3c', '#2ecc71', '#f39c12'])
    ax.set_title('Molecular Subtype Distribution')
    ax.set_ylabel('Count')
    ax.tick_params(axis='x', rotation=0)
    
    # 2. ZP3 by subtype
    ax = axes[0, 1]
    sns.boxplot(data=df, x='Molecular_subtype', y='ZP3', ax=ax, palette='husl')
    ax.set_title('ZP3 by Molecular Subtype')
    ax.tick_params(axis='x', rotation=45)
    
    # 3. ZP3 vs IFN-γ scatter
    ax = axes[0, 2]
    if 'IFN_gamma_score' in df.columns:
        ax.scatter(df['IFN_gamma_score'], np.log1p(df['ZP3']), alpha=0.5, s=20)
        # Add regression line
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
    
    # 4. Treatment response by subtype
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
    
    # 5. Response rate by subtype
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
    
    # 6. ZP3 and immunosuppressive gene correlation heatmap
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
    print(f"  Saved visualization: {fig_path}")

def create_report(df):
    """Create report"""
    report = f"""# GSE91061 Melanoma Molecular Subtype Analysis Report

## Analysis Information
- Cohort: GSE91061 (n=109, melanoma anti-PD-1 treatment)
- Treatment: Nivolumab (anti-PD-1)
- Data version: hg19, FPKM normalized

## 1. Molecular subtype definition
Based on TCGA SKCM transcriptome subtype classifier, samples are divided into 4 molecular subtypes:
- **Immune**: immune gene enrichment (CD3, CD8, CCL5, STAT1, CXCL9/10, etc.)
- **Keratin**: keratin gene enrichment (KRT5, KRT14, KRT17, etc.)
- **MITF-high**: high expression of MITF/pigmentation genes
- **MITF-low**: low expression of MITF/pigmentation genes

## 2. Molecular subtype distribution
"""
    report += "| Subtype | N | Percentage |\n|---|---|---|\n"
    for subtype, count in df['Molecular_subtype'].value_counts().items():
        report += f"| {subtype} | {count} | {100*count/len(df):.1f}% |\n"
    
    report += f"""
## 3. ZP3 expression by subtype
"""
    for subtype, group in df.groupby('Molecular_subtype'):
        report += f"- **{subtype}** (n={len(group)}): ZP3 mean={group['ZP3'].mean():.3f}, median={group['ZP3'].median():.3f}\n"
    
    report += f"""
## 4. Treatment response analysis
- Only analyze pre-treatment samples
- Response: PRCR (Partial/Complete Response)
- Non-response: SD (Stable Disease) + PD (Progressive Disease)

## 5. Immune features
- IFN-γ signature, CD8 signature, T-cell signature, and immunosuppression score have been calculated
- Correlation between ZP3 and immunosuppressive genes has been calculated

## 6. Multivariate analysis
- Model 1: response ~ ZP3
- Model 2: response ~ ZP3 + IFN-γ
- Model 3: response ~ ZP3 + Immunosuppression

## 7. Key findings
- ZP3 expression differences among molecular subtypes
- Whether the relationship between ZP3 and treatment response depends on molecular subtype
- Whether ZP3 is independent of known immune markers

## 8. Limitations
- GSE91061 does not provide BRAF/NRAS mutations and TMB data
- Molecular subtypes are inferred from expression signatures, not official TCGA classification
- Limited sample size, insufficient statistical power for within-subtype analyses

## 9. Output files
- `gse91061_molecular_subtype_clinical.csv` - Complete clinical + molecular data
- `fig_gse91061_molecular_subtype.png` - Visualization
"""
    
    report_path = f'{OUTPUT_DIR}/gse91061_molecular_subtype_report.md'
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"  Report created: {report_path}")

if __name__ == "__main__":
    main()
