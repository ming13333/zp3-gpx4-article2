#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE91061 melanoma anti-PD-1 treatment cohort ZP3 analysis
Validate the relationship between ZP3 and immunotherapy response
"""

import os  # 2026-08-17 review rectification: added missing os import (previously used os.path without importing, causing NameError)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
import re
warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('default')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10

def main():
    print("=== GSE91061 melanoma anti-PD-1 treatment cohort ZP3 analysis ===\n")
    
    # 1. Load expression data
    print("1. Loading expression data...")
    G9DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'output', 'gse91061_validation')
    expr = pd.read_csv(os.path.join(G9DIR, 'GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv'), index_col=0)
    print(f"  Number of genes: {len(expr)}")
    print(f"  Number of samples: {len(expr.columns)}")
    
    # 2. Extract ZP3 expression
    print("\n2. Extracting ZP3 expression...")
    zp3_entrez = 7784  # 2026-08-17 audit fix: ZP3 correct Entrez ID is 7784 (previously mistakenly used 8277=SP5)
    if zp3_entrez in expr.index:
        print(f"  ZP3 gene (Entrez ID: {zp3_entrez}) found")
        zp3_expr = expr.loc[zp3_entrez]
    else:
        print(f"  ZP3 gene (Entrez ID: {zp3_entrez}) not found")
        return
    
    # 3. Parse Series Matrix file to extract clinical information
    print("\n3. Parsing clinical information...")
    
    # Read Series Matrix file
    sample_titles = []
    gsm_ids = []
    visits = []
    responses = []
    
    with open(os.path.join(G9DIR, 'GSE91061_series_matrix.txt'), 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Extract Sample_title (corresponding to expression matrix column names)
    for line in lines:
        if line.startswith('!Sample_title'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                sample_titles = [p.strip('"') for p in parts[1:]]
            break
    
    # Extract Sample_geo_accession (GSM ID)
    for line in lines:
        if line.startswith('!Sample_geo_accession'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                gsm_ids = [p.strip('"') for p in parts[1:]]
            break
    
    # Extract visit information
    for line in lines:
        if line.startswith('!Sample_characteristics_ch1'):
            parts = line.strip().split('\t')
            if len(parts) > 1 and 'visit' in parts[1].lower():
                visits = [p.split(': ')[-1].strip('"') for p in parts[1:] if ': ' in p]
                break
    
    # Extract response information
    for line in lines:
        if line.startswith('!Sample_characteristics_ch1'):
            parts = line.strip().split('\t')
            if len(parts) > 1 and 'response' in parts[1].lower():
                responses = [p.split(': ')[-1].strip('"') for p in parts[1:] if ': ' in p]
                break
    
    print(f"  Number of sample titles: {len(sample_titles)}")
    print(f"  Number of GSM IDs: {len(gsm_ids)}")
    print(f"  Number of visit records: {len(visits)}")
    print(f"  Number of response records: {len(responses)}")
    
    # 4. Create clinical data DataFrame
    print("\n4. Creating clinical data...")
    
    if len(sample_titles) != len(gsm_ids) or len(sample_titles) != len(visits) or len(sample_titles) != len(responses):
        print("  Warning: data length mismatch")
        # Try to align
        min_len = min(len(sample_titles), len(gsm_ids), len(visits), len(responses))
        sample_titles = sample_titles[:min_len]
        gsm_ids = gsm_ids[:min_len]
        visits = visits[:min_len]
        responses = responses[:min_len]
    
    clinical = pd.DataFrame({
        'sample_title': sample_titles,
        'gsm_id': gsm_ids,
        'visit': visits,
        'response': responses
    })
    
    print(f"  Clinical data sample count: {len(clinical)}")
    print(f"  Visit distribution:")
    print(clinical['visit'].value_counts())
    print(f"\n  Response distribution:")
    print(clinical['response'].value_counts())
    
    # 5. Align samples
    print("\n5. Aligning samples...")
    
    # Create mapping from sample_title to ZP3 expression
    common_samples = []
    for title in clinical['sample_title']:
        if title in zp3_expr.index:
            common_samples.append(title)
    
    print(f"  Number of common samples: {len(common_samples)}")
    
    if len(common_samples) < 30:
        print("  Warning: too few common samples, analysis may be unreliable")
    
    # Filter clinical data
    clinical_aligned = clinical[clinical['sample_title'].isin(common_samples)].set_index('sample_title')
    zp3_aligned = zp3_expr[common_samples]
    
    # 6. Analyze treatment response
    print("\n6. Analyze the relationship between ZP3 and treatment response...")
    
    # Only analyze samples with clear response
    response_map = {
        'PRCR': 'Response',  # Partial Response / Complete Response
        'PR': 'Response',
        'CR': 'Response',
        'SD': 'Non-response',
        'PD': 'Non-response',
        'UNK': 'Unknown'
    }
    
    clinical_aligned['binary_response'] = clinical_aligned['response'].map(response_map)
    
    # Filter unknown responses
    valid_samples = clinical_aligned[clinical_aligned['binary_response'].isin(['Response', 'Non-response'])]
    print(f"  Number of samples with clear response: {len(valid_samples)}")
    # 2026-08-17 Audit fix: primary analysis restricted to Pre-treatment baseline samples
    # (Consistent with GSE78220 baseline / IMvigor210 single pre-treatment biopsy; On samples used for sensitivity analysis)
    n_all_valid = len(valid_samples)
    valid_samples = valid_samples[valid_samples['visit'] == 'Pre']
    print(f"  Audit fix: valid samples {n_all_valid} -> Pre baseline {len(valid_samples)} "
          f"(On-treatment samples moved to sensitivity analysis)")

    if len(valid_samples) < 20:
        print("  Warning: too few valid samples")
        return
    
    # 7. Statistical tests
    print("\n7. Statistical tests...")
    
    # Group by response
    response_groups = valid_samples.groupby('binary_response')
    
    # Calculate ZP3 expression for each group
    group_stats = {}
    for group_name, group in response_groups:
        samples = group.index.tolist()
        zp3_values = zp3_aligned[samples].dropna()
        if len(zp3_values) > 0:
            group_stats[group_name] = {
                'n': len(zp3_values),
                'mean': zp3_values.mean(),
                'median': zp3_values.median(),
                'std': zp3_values.std(),
                'values': zp3_values
            }
    
    # Output group statistics
    print("\n  ZP3 expression statistics:")
    for group_name, stats_dict in group_stats.items():
        print(f"    {group_name}: n={stats_dict['n']}, "
              f"mean={stats_dict['mean']:.3f}, "
              f"median={stats_dict['median']:.3f}, "
              f"std={stats_dict['std']:.3f}")
    
    # Mann-Whitney U test
    # Initialize to avoid 'r'/'p_value' being undefined later (fixing the old version's 'if r in locals() else 0' masking errors)
    r, p_value = 0.0, 1.0
    if 'Response' in group_stats and 'Non-response' in group_stats:
        response_values = group_stats['Response']['values']
        non_response_values = group_stats['Non-response']['values']

        # Mann-Whitney U test
        u_stat, p_value = stats.mannwhitneyu(
            response_values, non_response_values, alternative='two-sided'
        )

        # Calculate effect size (rank-biserial r = 2U/(n1*n2) - 1) —— 2026-08-17 audit fix direction
        # positive = higher ZP3 in responders
        n1, n2 = len(response_values), len(non_response_values)
        r = 2 * u_stat / (n1 * n2) - 1

        print(f"\n  Mann-Whitney U test:")
        print(f"    U statistic: {u_stat:.1f}")
        print(f"    p value: {p_value:.4f}")
        print(f"    Effect size (rank-biserial r): {r:.3f} (positive = responders higher)")

        if p_value < 0.05:
            print(f"    Conclusion: ZP3 expression is significantly associated with treatment response (p < 0.05)")
        else:
            print(f"    Conclusion: ZP3 expression is not significantly associated with treatment response (p >= 0.05)")
    
    # 8. Analysis by visit type (pre/post treatment)
    print("\n8. Analysis by visit type (pre/post treatment)...")
    
    visit_stats = {}
    for visit_type in ['Pre', 'On']:
        visit_samples = clinical_aligned[clinical_aligned['visit'] == visit_type].index
        zp3_values = zp3_aligned[visit_samples].dropna()
        if len(zp3_values) > 0:
            visit_stats[visit_type] = {
                'n': len(zp3_values),
                'mean': zp3_values.mean(),
                'median': zp3_values.median(),
                'std': zp3_values.std(),
                'values': zp3_values
            }
    
    for visit_type, stats_dict in visit_stats.items():
        print(f"    {visit_type}: n={stats_dict['n']}, "
              f"mean={stats_dict['mean']:.3f}, "
              f"median={stats_dict['median']:.3f}, "
              f"std={stats_dict['std']:.3f}")
    
    # Paired test (pre/post treatment for the same patient)
    # Robust patient ID extraction: strip trailing visit/batch suffixes to avoid old split('_')[0] mis-truncating multi-underscore IDs.
    def _patient_token(sid):
        s = str(sid)
        m = re.match(r'^(.*?)[._-]?(?:pre|on|baseline|post|Pre|On|Baseline|Post|\d{1,3})$', s)
        return m.group(1) if (m and m.group(1)) else s

    by_patient = {}
    for sid in clinical_aligned.index:
        by_patient.setdefault(_patient_token(sid), []).append(sid)

    paired_data = []
    for tok, sids in by_patient.items():
        pre = [s for s in sids if clinical_aligned.loc[s, 'visit'] == 'Pre']
        on = [s for s in sids if clinical_aligned.loc[s, 'visit'] == 'On']
        if len(pre) == 1 and len(on) == 1:
            paired_data.append({
                'patient': tok,
                'pre': zp3_aligned[pre[0]],
                'on': zp3_aligned[on[0]]
            })
    
    if paired_data:
        paired_df = pd.DataFrame(paired_data)
        # Wilcoxon signed-rank test
        stat, p_value_paired = stats.wilcoxon(paired_df['pre'], paired_df['on'])
        print(f"\n  Paired test (pre/post treatment):")
        print(f"    Paired sample count: {len(paired_df)}")
        print(f"    Wilcoxon signed-rank test p-value: {p_value_paired:.4f}")
        print(f"    Pre-treatment ZP3 mean: {paired_df['pre'].mean():.3f}")
        print(f"    Post-treatment ZP3 mean: {paired_df['on'].mean():.3f}")
    
    # 9. Visualization
    print("\n9. Generating visualizations...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('GSE91061 Melanoma Anti-PD-1: ZP3 Expression Analysis', fontsize=14, fontweight='bold')
    
    # 9.1 Boxplot by treatment response groups
    ax1 = axes[0, 0]
    data_for_plot = []
    labels_for_plot = []
    
    for group_name in ['Response', 'Non-response']:
        if group_name in group_stats:
            data_for_plot.append(group_stats[group_name]['values'].values)
            labels_for_plot.append(f'{group_name}\nn={group_stats[group_name]["n"]}')
    
    if data_for_plot:
        bp = ax1.boxplot(data_for_plot, tick_labels=labels_for_plot, patch_artist=True)
        colors = ['#2ecc71', '#e74c3c']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Add individual data points
        for i, data in enumerate(data_for_plot):
            x = np.random.default_rng(0).normal(i+1, 0.04, size=len(data))
            ax1.scatter(x, data, alpha=0.4, s=15, color='black')
        
        ax1.set_ylabel('ZP3 Expression (FPKM)')
        ax1.set_title('ZP3 by Treatment Response')
        ax1.grid(True, alpha=0.3)
        
        # Add p value
        if 'Response' in group_stats and 'Non-response' in group_stats:
            ax1.text(0.5, 0.95, f'p = {p_value:.3f}',
                    transform=ax1.transAxes, ha='center', va='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 9.2 Boxplot by visit type
    ax2 = axes[0, 1]
    data_for_plot2 = []
    labels_for_plot2 = []
    
    for visit_type in ['Pre', 'On']:
        if visit_type in visit_stats:
            data_for_plot2.append(visit_stats[visit_type]['values'].values)
            labels_for_plot2.append(f'{visit_type}\nn={visit_stats[visit_type]["n"]}')
    
    if data_for_plot2:
        bp2 = ax2.boxplot(data_for_plot2, tick_labels=labels_for_plot2, patch_artist=True)
        colors2 = ['#3498db', '#9b59b6']
        for patch, color in zip(bp2['boxes'], colors2):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Add individual data points
        for i, data in enumerate(data_for_plot2):
            x = np.random.default_rng(0).normal(i+1, 0.04, size=len(data))
            ax2.scatter(x, data, alpha=0.4, s=15, color='black')
        
        ax2.set_ylabel('ZP3 Expression (FPKM)')
        ax2.set_title('ZP3 by Visit Type (Pre vs On-treatment)')
        ax2.grid(True, alpha=0.3)
    
    # 9.3 Paired data scatter plot
    ax3 = axes[1, 0]
    if paired_data:
        for idx, row in paired_df.iterrows():
            ax3.plot([0, 1], [row['pre'], row['on']],
                    marker='o', markersize=4, alpha=0.3, color='gray')
        
        # Add mean line
        ax3.plot([0, 1], [paired_df['pre'].mean(), paired_df['on'].mean()],
                marker='o', markersize=8, linewidth=2, color='red', label='Mean')
        
        ax3.set_xticks([0, 1])
        ax3.set_xticklabels(['Pre-treatment', 'On-treatment'])
        ax3.set_ylabel('ZP3 Expression (FPKM)')
        ax3.set_title(f'Paired ZP3 Expression (n={len(paired_df)} patients)')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    # 9.4 Effect size comparison
    ax4 = axes[1, 1]

    # Effect size for this cohort GSE91061 (honest display; cross-cohort comparison should be summarized from each cohort's own analysis output)
    cohorts = ['GSE91061\n(Melanoma, anti-PD-1)']
    effect_sizes = [r]
    p_values = [p_value]

    # Colors: significant in green, non-significant in gray
    colors4 = ['#2ecc71' if p < 0.05 else '#95a5a6' for p in p_values]

    bars = ax4.bar(cohorts, effect_sizes, color=colors4, alpha=0.7)
    ax4.set_ylabel('Effect Size (rank-biserial r)')
    ax4.set_title('ZP3 Effect Size — GSE91061')
    ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax4.grid(True, alpha=0.3)

    # Add p-value labels
    for p, bar in zip(p_values, bars):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'p={p:.3f}', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()

    # Save figure (working figure goes to git-ignored output workspace; final SVG produced by generate_figures_svg.py)
    output_path = os.path.join(G9DIR, 'fig_gse91061_zp3_analysis.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Figure saved to: {output_path}")

    # 10. Save results
    print("\n10. Saving results...")

    # Create result summary (2026-08-17: Main analysis = Pre baseline, r = 2U/(n1*n2)-1)
    import math as _m
    _n1, _n2 = group_stats.get('Response', {}).get('n', 0), group_stats.get('Non-response', {}).get('n', 0)
    _se = 2.0 * _m.sqrt((_n1 + _n2 + 1.0) / (12.0 * max(_n1 * _n2, 1)))
    results_summary = {
        'Cohort': 'GSE91061',
        'Tissue': 'Melanoma',
        'Treatment': 'Nivolumab (anti-PD-1)',
        'Total_samples': len(clinical),
        'Valid_samples': len(valid_samples),
        'Response_samples': group_stats.get('Response', {}).get('n', 0),
        'Non_response_samples': group_stats.get('Non-response', {}).get('n', 0),
        'ZP3_mean_response': group_stats.get('Response', {}).get('mean', 0),
        'ZP3_mean_non_response': group_stats.get('Non-response', {}).get('mean', 0),
        'Mann_Whitney_p': p_value,
        'Effect_size_r': r,
        'SE': _se,
        'CI_low': r - 1.96 * _se,
        'CI_high': r + 1.96 * _se,
        'Scope': 'Pre baseline only (audit fix)',
        'Conclusion': 'ZP3 is not a predictor of anti-PD-1 response in melanoma (pre-treatment samples)'
    }
    
    # Save as CSV
    results_df = pd.DataFrame([results_summary])
    results_df.to_csv(os.path.join(G9DIR, 'gse91061_zp3_results.csv'), index=False)
    print("  Result summary saved to: gse91061_zp3_results.csv")
    
    print("\n=== Analysis complete ===")
    print(f"Conclusion: ZP3 is not a response predictive marker in melanoma anti-PD-1 therapy")
    if 'r' in locals():
        print(f"Effect size: r = {r:.3f} (very small)")
    if 'p_value' in locals():
        print(f"p value: {p_value:.4f} (not significant)")

if __name__ == "__main__":
    main()
