#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE91061 Melanoma anti-PD-1 treatment cohort ZP3 analysis
Validate the relationship between ZP3 and immunotherapy response
"""

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

def parse_series_matrix(matrix_file):
    """Parse Series Matrix file to extract clinical information"""
    sample_ids = []
    visits = []
    responses = []
    tissues = []
    
    with open(matrix_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('!Series_sample_id'):
                # Extract sample IDs
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    sample_ids = parts[1].strip('"').split()
            elif line.startswith('!Sample_characteristics_ch1'):
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    # Extract visit information (first ch1 row)
                    if not visits:
                        visits = [p.split(': ')[-1].strip('"') for p in parts[1:] if 'visit' in p.lower()]
                    # Extract response information (second ch1 row)
                    elif not responses:
                        responses = [p.split(': ')[-1].strip('"') for p in parts[1:] if 'response' in p.lower()]
                    # Extract tissue information (third ch1 row)
                    elif not tissues:
                        tissues = [p.split(': ')[-1].strip('"') for p in parts[1:] if 'tissue' in p.lower()]
    
    # Create clinical data DataFrame
    if sample_ids and responses:
        clinical = pd.DataFrame({
            'sample_id': sample_ids,
            'visit': visits if len(visits) == len(sample_ids) else ['UNK'] * len(sample_ids),
            'response': responses,
            'tissue': tissues if len(tissues) == len(sample_ids) else ['melanoma'] * len(sample_ids)
        })
        return clinical
    else:
        return None

def create_sample_mapping(expr_columns):
    """Create mapping from expression matrix column names to GSM sample names"""
    # Expression matrix column name format: Pt1_Pre_AD101148-6
    # Clinical data sample name format: GSM2420259
    # Need to extract the mapping from the Series Matrix file
    
    # Read the Series Matrix file to obtain the mapping from GSM to patient information
    mapping = {}
    with open(os.path.join(G9DIR, 'GSE91061_series_matrix.txt'), 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    # Find the Sample_title line (contains GSM ID and sample description)
    sample_titles = []
    for line in lines:
        if line.startswith('!Sample_title'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                sample_titles = parts[1:]
            break
    
    # Find the Sample_geo_accession line (GSM ID)
    gsm_ids = []
    for line in lines:
        if line.startswith('!Sample_geo_accession'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                gsm_ids = [g.strip('"') for g in parts[1:]]  # Remove quotes from each GSM
            break
    
    # Create mapping from GSM to sample description
    gsm_to_desc = {}
    for i, gsm in enumerate(gsm_ids):
        if i < len(sample_titles):
            desc = sample_titles[i].strip('"')
            gsm_to_desc[gsm] = desc
    
    # Create mapping from sample description to GSM
    desc_to_gsm = {desc: gsm for gsm, desc in gsm_to_desc.items()}
    
    # Extract patient ID and visit type from sample description
    # Format example: "Pt1_Pre_AD101148-6"
    expr_to_gsm = {}
    for expr_sample in expr_columns:
        # Directly use expression matrix column names as sample description
        if expr_sample in desc_to_gsm:
            expr_to_gsm[expr_sample] = desc_to_gsm[expr_sample]
    
    return expr_to_gsm

def main():
    print("=== GSE91061 melanoma anti-PD-1 treatment cohort ZP3 analysis ===\n")
    
    # 1. Load clinical data
    print("1. Loading clinical data...")
    clinical = parse_series_matrix('GSE91061_series_matrix.txt')
    if clinical is None:
        print("Error: Unable to parse clinical data")
        return
    
    print(f"  Total samples: {len(clinical)}")
    print(f"  Visit type distribution:")
    print(clinical['visit'].value_counts())
    print(f"\n  Treatment response distribution:")
    print(clinical['response'].value_counts())
    
    # 2. Load expression data
    print("\n2. Loading expression data...")
    G9DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'output', 'gse91061_validation')
    expr = pd.read_csv(os.path.join(G9DIR, 'GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv'), index_col=0)
    print(f"  Gene count: {len(expr)}")
    print(f"  Sample count: {len(expr.columns)}")
    
    # 3. Check ZP3 gene
    # Note: expression matrix uses Entrez Gene IDs; need to convert ZP3 to numeric ID
    # ZP3's Entrez Gene ID is 7784 (2026-08-17 audit fix; this script has been superseded by v2, kept only for historical purposes)
    zp3_entrez = 7784
    if zp3_entrez in expr.index:
        print(f"\n3. ZP3 gene (Entrez ID: {zp3_entrez}) found")
        zp3_expr = expr.loc[zp3_entrez]
    else:
        print(f"\n3. ZP3 gene (Entrez ID: {zp3_entrez}) not found")
        print("  Trying to find other ZP3-related genes...")
        # Find possible ZP3-related genes
        zp3_genes = [idx for idx in expr.index if 'ZP3' in str(idx)]
        if zp3_genes:
            print(f"  Found {len(zp3_genes)} ZP3-related genes: {zp3_genes[:5]}")
        else:
            print("  No ZP3-related genes found")
            return
    
    # 4. Align samples
    print("\n4. Aligning samples...")
    
    # Create mapping
    sample_mapping = create_sample_mapping(expr.columns)
    print(f"  Number of successfully mapped samples: {len(sample_mapping)}")
    
    if len(sample_mapping) < 30:
        print("  Warning: too few successfully mapped samples, analysis may be unreliable")
    
    # Reorganize clinical data to match expression matrix
    mapped_clinical = []
    for expr_sample, gsm in sample_mapping.items():
        if gsm in clinical['sample_id'].values:
            row = clinical[clinical['sample_id'] == gsm].iloc[0]
            mapped_clinical.append({
                'sample_id': expr_sample,
                'gsm': gsm,
                'visit': row['visit'],
                'response': row['response'],
                'tissue': row['tissue']
            })
    
    clinical_aligned = pd.DataFrame(mapped_clinical).set_index('sample_id')
    zp3_aligned = zp3_expr[clinical_aligned.index]
    
    common_samples = clinical_aligned.index.tolist()
    print(f"  Final common sample count: {len(common_samples)}")
    
    # 5. Analyze treatment response
    print("\n5. Analyze the relationship between ZP3 and treatment response...")
    
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
    
    # Filter out unknown responses
    valid_samples = clinical_aligned[clinical_aligned['binary_response'].isin(['Response', 'Non-response'])]
    print(f"  Number of samples with clear response: {len(valid_samples)}")
    
    if len(valid_samples) < 20:
        print("  Warning: too few valid samples")
        return
    
    # 6. Statistical test
    print("\n6. Statistical test...")
    
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
    # Initialize to avoid subsequent 'r'/'p_value' being undefined (fixes old version where 'if r in locals() else 0' masked errors)
    r, p_value = 0.0, 1.0
    if 'Response' in group_stats and 'Non-response' in group_stats:
        response_values = group_stats['Response']['values']
        non_response_values = group_stats['Non-response']['values']

        # Mann-Whitney U test
        u_stat, p_value = stats.mannwhitneyu(
            response_values, non_response_values, alternative='two-sided'
        )

        # Calculate effect size (rank-biserial r = 1 - 2U/(n1*n2)) — standard Mann-Whitney effect size
        n1, n2 = len(response_values), len(non_response_values)
        r = 1 - (2 * u_stat) / (n1 * n2)

        print(f"\n  Mann-Whitney U test:")
        print(f"    U statistic: {u_stat:.1f}")
        print(f"    p value: {p_value:.4f}")
        print(f"    Effect size (rank-biserial r): {r:.3f}")

        if p_value < 0.05:
            print(f"    Conclusion: ZP3 expression is significantly associated with treatment response (p < 0.05)")
        else:
            print(f"    Conclusion: ZP3 expression is not significantly associated with treatment response (p >= 0.05)")
    
    # 7. Analysis by visit type (before/after treatment)
    print("\n7. Analysis by visit type (before/after treatment)...")
    
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
    
    # Paired test (before/after treatment for the same patient)
    # Robustly extract patient ID: remove trailing visit/batch suffix (pre/on/baseline/post or digits),
    # Avoid the old split('_')[0] incorrectly truncating IDs with multiple underscores (e.g., TCGA-AB-1234-01).
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
        stat, p_value = stats.wilcoxon(paired_df['pre'], paired_df['on'])
        print(f"\n  Paired test (before/after treatment):")
        print(f"    Number of paired samples: {len(paired_df)}")
        print(f"    Wilcoxon signed-rank test p value: {p_value:.4f}")
        print(f"    Pre-treatment ZP3 mean: {paired_df['pre'].mean():.3f}")
        print(f"    Post-treatment ZP3 mean: {paired_df['on'].mean():.3f}")
    
    # 8. Visualization
    print("\n8. Generating visualizations...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('GSE91061 Melanoma Anti-PD-1: ZP3 Expression Analysis', fontsize=14, fontweight='bold')
    
    # 8.1 Boxplot by treatment response group
    ax1 = axes[0, 0]
    data_for_plot = []
    labels_for_plot = []
    
    for group_name in ['Response', 'Non-response']:
        if group_name in group_stats:
            data_for_plot.append(group_stats[group_name]['values'].values)
            labels_for_plot.append(f'{group_name}\nn={group_stats[group_name]["n"]}')
    
    if data_for_plot:
        bp = ax1.boxplot(data_for_plot, labels=labels_for_plot, patch_artist=True)
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
        
        # Add p-value
        if 'Response' in group_stats and 'Non-response' in group_stats:
            ax1.text(0.5, 0.95, f'p = {p_value:.3f}', 
                    transform=ax1.transAxes, ha='center', va='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 8.2 Boxplot by visit type group
    ax2 = axes[0, 1]
    data_for_plot2 = []
    labels_for_plot2 = []
    
    for visit_type in ['Pre', 'On']:
        if visit_type in visit_stats:
            data_for_plot2.append(visit_stats[visit_type]['values'].values)
            labels_for_plot2.append(f'{visit_type}\nn={visit_stats[visit_type]["n"]}')
    
    if data_for_plot2:
        bp2 = ax2.boxplot(data_for_plot2, labels=labels_for_plot2, patch_artist=True)
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
    
    # 8.3 Paired data scatter plot
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
    
    # 8.4 Effect size (this cohort GSE91061, presented honestly; cross-cohort comparison should be summarized from each cohort's own analysis output)
    ax4 = axes[1, 1]
    cohorts = ['GSE91061\n(Melanoma, anti-PD-1)']
    effect_sizes = [r]
    p_values = [p_value]
    colors4 = ['#2ecc71' if p < 0.05 else '#95a5a6' for p in p_values]
    bars = ax4.bar(cohorts, effect_sizes, color=colors4, alpha=0.7)
    ax4.set_ylabel('Effect Size (rank-biserial r)')
    ax4.set_title('ZP3 Effect Size — GSE91061')
    ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax4.grid(True, alpha=0.3)
    for p, bar in zip(p_values, bars):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'p={p:.3f}', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()

    # Save figure
    output_path = 'fig_gse91061_zp3_analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Figure saved to: {output_path}")

    # 9. Save results
    print("\n9. Save results...")

    # Create results summary
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
        'Conclusion': ('ZP3 expression is significantly associated with anti-PD-1 response in melanoma '
                      '(higher ZP3 linked to poorer response)' if p_value < 0.05
                      else 'ZP3 expression shows no significant association with anti-PD-1 response in melanoma')
    }
    
    # Save as CSV
    results_df = pd.DataFrame([results_summary])
    results_df.to_csv(os.path.join(G9DIR, 'gse91061_zp3_results.csv'), index=False)
    print("  Results summary saved to: gse91061_zp3_results.csv")
    
    print("\n=== Analysis complete ===")
    _nr_mean = group_stats.get('Non-response', {}).get('mean', 0)
    _r_mean = group_stats.get('Response', {}).get('mean', 0)
    _direction = ("High ZP3 expression is associated with poorer treatment response (non-response)"
                  if _nr_mean > _r_mean
                  else "High ZP3 expression is associated with better treatment response (response)")
    _effect_label = "medium" if abs(r) >= 0.3 else ("small" if abs(r) >= 0.1 else "very small")
    _sig_label = "significant association" if p_value < 0.05 else "no significant association"
    print(f"Conclusion: ZP3 shows {_sig_label} with treatment response to anti-PD-1 in melanoma (Mann-Whitney U, p={p_value:.4f})")
    print(f"Direction: {_direction} (responder mean {_r_mean:.3f} vs non-responder {_nr_mean:.3f})")
    print(f"Effect size: rank-biserial r = {r:.3f} ({_effect_label} effect)")
    print(f"p value: {p_value:.4f} ({'significant' if p_value < 0.05 else 'not significant'})")

if __name__ == "__main__":
    main()
