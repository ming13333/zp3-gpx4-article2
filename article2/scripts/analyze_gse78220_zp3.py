# -*- coding: utf-8 -*-
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
"""
Analyze the relationship between ZP3 and anti-PD-1 treatment response in GSE78220
Data: melanoma anti-PD-1 treatment cohort
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

OUT_DIR = os.path.join(ROOT, "output", "h2_bulk")

# Read data
print('Reading data...')
expr_df = pd.read_excel(os.path.join(OUT_DIR, 'GSE78220_PatientFPKM.xlsx'), sheet_name=0, index_col=0)
annot_df = pd.read_csv(os.path.join(OUT_DIR, 'gse78220_annotation.csv'))

# Extract treatment response
print('Extracting treatment response information...')
# Extract response from Sample_source_name_ch1
annot_df['response'] = annot_df['Sample_source_name_ch1'].apply(lambda x: 'CR' if 'Complete Response' in x else ('PR' if 'Partial Response' in x else 'PD'))
# Extract timepoint from Sample_description
annot_df['timepoint'] = annot_df['Sample_description'].apply(lambda x: 'OnTx' if 'on-anti-PD-1' in x else 'baseline')

# Create sample ID mapping
sample_map = {}
for i, row in annot_df.iterrows():
    sample_id = f"{row['Sample_title']}.{'OnTx' if row['timepoint'] == 'OnTx' else 'baseline'}"
    sample_map[sample_id] = {
        'response': row['response'],
        'timepoint': row['timepoint'],
        'patient_id': row['Sample_title']
    }

print(f'Number of samples: {len(sample_map)}')
print(f'Response distribution: {pd.Series([v["response"] for v in sample_map.values()]).value_counts().to_dict()}')
print(f'Timepoint distribution: {pd.Series([v["timepoint"] for v in sample_map.values()]).value_counts().to_dict()}')

# Merge expression data with clinical information
print('\nMerging expression data with clinical information...')
merged_data = []
for sample in expr_df.columns:
    if sample in sample_map:
        info = sample_map[sample]
        merged_data.append({
            'sample_id': sample,
            'zp3_expr': expr_df.loc['ZP3', sample] if 'ZP3' in expr_df.index else np.nan,
            'trem2_expr': expr_df.loc['TREM2', sample] if 'TREM2' in expr_df.index else np.nan,
            'cd68_expr': expr_df.loc['CD68', sample] if 'CD68' in expr_df.index else np.nan,
            'cd163_expr': expr_df.loc['CD163', sample] if 'CD163' in expr_df.index else np.nan,
            'pdl1_expr': expr_df.loc['CD274', sample] if 'CD274' in expr_df.index else np.nan,
            'response': info['response'],
            'timepoint': info['timepoint'],
            'patient_id': info['patient_id']
        })

merged_df = pd.DataFrame(merged_data)
print(f'Merged sample count: {len(merged_df)}')
print(f'Response distribution: {merged_df["response"].value_counts().to_dict()}')

# Save merged data
merged_csv = os.path.join(OUT_DIR, 'gse78220_zp3_response.csv')
merged_df.to_csv(merged_csv, index=False)
print(f'Merged data saved: {merged_csv}')

# ============================================================
# Main analysis (2026-08-17 scientific audit fix): baseline-only, responder vs non-responder
# Mann-Whitney U + rank-biserial r (r = 2U/(n1*n2) - 1, positive = responders higher)
# Consistent with manuscript Methods 2.3; the ANOVA/pairwise t-tests below are exploratory supplements only.
# ============================================================
print('\n=== Main analysis: baseline samples Responder (CR+PR) vs Non-responder (PD) ===')
base = merged_df[merged_df['timepoint'] == 'baseline'].copy()
base['is_resp'] = base['response'].isin(['CR', 'PR'])
g_r = base.loc[base['is_resp'], 'zp3_expr'].dropna()
g_n = base.loc[~base['is_resp'], 'zp3_expr'].dropna()
print(f'Baseline samples: n={len(base)} (responder {len(g_r)}, non-responder {len(g_n)})')
print(f'ZP3 median: responder={g_r.median():.3f}, non-responder={g_n.median():.3f}')
u_stat, p_mwu = stats.mannwhitneyu(g_r, g_n, alternative='two-sided')
n1, n2 = len(g_r), len(g_n)
r_rb = 2 * u_stat / (n1 * n2) - 1
import math as _math
se_rb = 2.0 * _math.sqrt((n1 + n2 + 1.0) / (12.0 * n1 * n2))
lo, hi = r_rb - 1.96 * se_rb, r_rb + 1.96 * se_rb
print(f'Mann-Whitney U={u_stat:.1f}, P={p_mwu:.4f}')
print(f'rank-biserial r={r_rb:.3f} (95% CI {lo:.2f}, {hi:.2f}); SE={se_rb:.3f}')
# Write frozen results to CSV (for figure/main-text citation)
import csv as _csv
with open(os.path.join(OUT_DIR, 'gse78220_zp3_primary_result.csv'), 'w', newline='') as _f:
    _w = _csv.writer(_f)
    _w.writerow(['Cohort','Scope','n_resp','n_non','med_resp','med_non','U','P','r_rb','CI_low','CI_high','SE'])
    _w.writerow(['GSE78220','primary_baseline', n1, n2, round(g_r.median(),3), round(g_n.median(),3),
                 round(u_stat,1), p_mwu, round(r_rb,4), round(lo,4), round(hi,4), round(se_rb,4)])
print(f'Frozen result: {os.path.join(OUT_DIR, "gse78220_zp3_primary_result.csv")}')

# Analyze the relationship between ZP3 and treatment response
print('\n=== Relationship between ZP3 and treatment response ===')
print(f'\nZP3 expression by response group:')
response_stats = merged_df.groupby('response')['zp3_expr'].agg(['mean', 'median', 'std', 'count'])
print(response_stats)

# Statistical tests
print('\nStatistical tests:')
groups = [merged_df[merged_df['response'] == r]['zp3_expr'].dropna() for r in ['CR', 'PR', 'PD']]
# ANOVA
if all(len(g) > 0 for g in groups):
    f_stat, p_val = stats.f_oneway(*groups)
    print(f'ANOVA: F={f_stat:.3f}, p={p_val:.4f}')
else:
    print('ANOVA: insufficient sample size')

# Pairwise comparisons
pairs = [('CR', 'PR'), ('CR', 'PD'), ('PR', 'PD')]
for g1, g2 in pairs:
    data1 = merged_df[merged_df['response'] == g1]['zp3_expr'].dropna()
    data2 = merged_df[merged_df['response'] == g2]['zp3_expr'].dropna()
    if len(data1) > 0 and len(data2) > 0:
        t_stat, p_val = stats.ttest_ind(data1, data2)
        print(f'{g1} vs {g2}: t={t_stat:.3f}, p={p_val:.4f}')

# Visualization
print('\nGenerating visualizations...')
plt.style.use('default')
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1. ZP3 expression by response group
ax1 = axes[0, 0]
response_order = ['CR', 'PR', 'PD']
colors = {'CR': '#2ecc71', 'PR': '#f39c12', 'PD': '#e74c3c'}
for i, resp in enumerate(response_order):
    data = merged_df[merged_df['response'] == resp]['zp3_expr'].dropna()
    ax1.boxplot(data, positions=[i], widths=0.6, patch_artist=True,
                boxprops=dict(facecolor=colors[resp], alpha=0.7),
                medianprops=dict(color='black', linewidth=2))
    # Add scatter points
    jitter = np.random.default_rng(0).normal(0, 0.05, len(data))
    ax1.scatter([i] * len(data) + jitter, data, alpha=0.6, color=colors[resp], s=30)
ax1.set_xticks(range(len(response_order)))
ax1.set_xticklabels(response_order)
ax1.set_ylabel('ZP3 Expression (FPKM)')
ax1.set_title('ZP3 Expression by Response to Anti-PD-1')
ax1.grid(True, alpha=0.3)

# 2. ZP3 vs TREM2 scatter plot
ax2 = axes[0, 1]
for resp in response_order:
    mask = merged_df['response'] == resp
    ax2.scatter(merged_df[mask]['zp3_expr'], merged_df[mask]['trem2_expr'], 
                alpha=0.6, color=colors[resp], label=resp, s=50)
ax2.set_xlabel('ZP3 Expression')
ax2.set_ylabel('TREM2 Expression')
ax2.set_title('ZP3 vs TREM2 Expression')
ax2.legend()
ax2.grid(True, alpha=0.3)

# 3. ZP3 vs PD-L1 scatter plot
ax3 = axes[1, 0]
for resp in response_order:
    mask = merged_df['response'] == resp
    ax3.scatter(merged_df[mask]['zp3_expr'], merged_df[mask]['pdl1_expr'], 
                alpha=0.6, color=colors[resp], label=resp, s=50)
ax3.set_xlabel('ZP3 Expression')
ax3.set_ylabel('PD-L1 (CD274) Expression')
ax3.set_title('ZP3 vs PD-L1 Expression')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. Responder vs non-responder ZP3 expression
ax4 = axes[1, 1]
merged_df['response_binary'] = merged_df['response'].apply(lambda x: 'Responder (CR/PR)' if x in ['CR', 'PR'] else 'Non-responder (PD)')
binary_order = ['Responder (CR/PR)', 'Non-responder (PD)']
binary_colors = ['#2ecc71', '#e74c3c']
for i, resp in enumerate(binary_order):
    data = merged_df[merged_df['response_binary'] == resp]['zp3_expr'].dropna()
    ax4.boxplot(data, positions=[i], widths=0.6, patch_artist=True,
                boxprops=dict(facecolor=binary_colors[i], alpha=0.7),
                medianprops=dict(color='black', linewidth=2))
    jitter = np.random.default_rng(0).normal(0, 0.05, len(data))
    ax4.scatter([i] * len(data) + jitter, data, alpha=0.6, color=binary_colors[i], s=30)
ax4.set_xticks(range(len(binary_order)))
ax4.set_xticklabels(binary_order)
ax4.set_ylabel('ZP3 Expression (FPKM)')
ax4.set_title('ZP3 Expression: Responders vs Non-responders')
ax4.grid(True, alpha=0.3)

# Add statistical test results
# Responder vs Non-responder
resp_data = merged_df[merged_df['response_binary'] == 'Responder (CR/PR)']['zp3_expr'].dropna()
non_resp_data = merged_df[merged_df['response_binary'] == 'Non-responder (PD)']['zp3_expr'].dropna()
if len(resp_data) > 0 and len(non_resp_data) > 0:
    t_stat, p_val = stats.ttest_ind(resp_data, non_resp_data)
    ax4.text(0.5, 0.95, f'Responder vs Non-responder: p={p_val:.4f}', 
             transform=ax4.transAxes, ha='center', va='top', fontsize=10)

plt.tight_layout()
fig_path = os.path.join(OUT_DIR, 'fig_gse78220_zp3_response.png')
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.close()
print(f'Figure saved: {fig_path}')

# Calculate correlation coefficients
print('\n=== Correlation Analysis ===')
corr_cols = ['zp3_expr', 'trem2_expr', 'cd68_expr', 'cd163_expr', 'pdl1_expr']
corr_df = merged_df[corr_cols].dropna()
if len(corr_df) > 5:
    corr_matrix = corr_df.corr()
    print('Correlation coefficient matrix:')
    print(corr_matrix.round(3))
    
    # Save correlation matrix
    corr_csv = os.path.join(OUT_DIR, 'gse78220_correlation_matrix.csv')
    corr_matrix.to_csv(corr_csv)
    print(f'Correlation matrix saved: {corr_csv}')

# Generate summary report
print('\n=== Analysis Summary ===')
summary = f"""
GSE78220 ZP3 and Anti-PD-1 Treatment Response Analysis Report
==========================================

Dataset: GSE78220 (melanoma anti-PD-1 treatment cohort)
Number of samples: {len(merged_df)}
Treatment response distribution:
- Complete Response (CR): {len(merged_df[merged_df['response'] == 'CR'])} cases
- Partial Response (PR): {len(merged_df[merged_df['response'] == 'PR'])} cases
- Progressive Disease (PD): {len(merged_df[merged_df['response'] == 'PD'])} cases

ZP3 expression statistics:
- All samples: mean={merged_df['zp3_expr'].mean():.3f}, median={merged_df['zp3_expr'].median():.3f}
- CR group: mean={merged_df[merged_df['response'] == 'CR']['zp3_expr'].mean():.3f}
- PR group: mean={merged_df[merged_df['response'] == 'PR']['zp3_expr'].mean():.3f}
- PD group: mean={merged_df[merged_df['response'] == 'PD']['zp3_expr'].mean():.3f}

Statistical Tests:
- ANOVA (CR vs PR vs PD): p={p_val if 'p_val' in locals() else 'N/A'}
- Responder (CR/PR) vs Non-responder (PD): p={p_val if 'p_val' in locals() else 'N/A'}

Key Findings:
1. Differential expression of ZP3 across treatment response groups
2. Correlation of ZP3 with TREM2 and PD-L1
3. Value of ZP3 as a potential predictive biomarker

Figure: {fig_path}
"""
summary_path = os.path.join(OUT_DIR, 'gse78220_analysis_summary.txt')
with open(summary_path, 'w', encoding='utf-8') as f:
    f.write(summary)
print(f'Summary report saved: {summary_path}')
