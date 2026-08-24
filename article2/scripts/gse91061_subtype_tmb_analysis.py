#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE91061 (Riaz 2017 Cell) Molecular subtype + TMB + BRAF/NRAS combined analysis

Data sources:
- GSE91061 expression matrix (FPKM, Entrez ID)
- bms038_clinical_data.csv (Riaz GitHub): SubtypeEZ (BRAF/RAS/NF1/TripleWT), BOR, myBOR, PFS, OS
- pre_therapy_nonsynonmous_mutations.csv (Riaz GitHub): mutations → TMB + BRAF/NRAS status
- SampleTableCorrected.9.19.16.csv: sample mapping

Analysis:
1. ZP3 expression compared by SubtypeEZ strata
2. Correlation between TMB and ZP3
3. BRAF/NRAS mutation status and ZP3
4. ZP3-response relationship stratified by subtype
5. ZP3-response relationship stratified by high/low TMB
6. Multivariable logistic regression
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import os, sys, json, warnings
warnings.filterwarnings('ignore')

DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output", "gse91061_validation")
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article2", "results", "gse91061_mutation_analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_data():
    print("=" * 60)
    print("GSE91061 molecular subtypes + TMB comprehensive analysis")
    print("=" * 60)
    
    # 1. Expression data
    print("\n[1/6] Loading expression data...")
    expr = pd.read_csv(os.path.join(DIR, 'GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv'), index_col=0)
    zp3_expr = expr.loc[7784]  # ZP3 (2026-08-17 audit fix: 7784, originally mistakenly used 8277=SP5)
    print(f"  Expression matrix: {len(expr)} genes × {len(expr.columns)} samples")
    
    # 2. Clinical data
    clin = pd.read_csv(os.path.join(DIR, 'bms038_clinical_data.csv'))
    clin = clin.set_index('PatientID')
    print(f"  Clinical data: {len(clin)} patients")
    
    # 3. Mutation data → TMB
    print("\n[2/6] Calculating TMB and BRAF/NRAS mutations...")
    mut = pd.read_csv(os.path.join(DIR, 'pre_therapy_nonsynonmous_mutations.csv'))
    print(f"  Total mutations: {len(mut)}")
    
    # Number of nonsynonymous mutations per patient = TMB proxy
    tmb_by_patient = mut.groupby('Patient').size().to_dict()
    
    # BRAF/NRAS mutation status
    braf_mut = mut[mut['Hugo Symbol'] == 'BRAF'].groupby('Patient')['HGVS_p'].apply(list).to_dict()
    nras_mut = mut[mut['Hugo Symbol'] == 'NRAS'].groupby('Patient')['HGVS_p'].apply(list).to_dict()
    
    braf_patients = set(braf_mut.keys())
    nras_patients = set(nras_mut.keys())
    print(f"  BRAF mutated patients: {len(braf_patients)}")
    print(f"  NRAS mutated patients: {len(nras_patients)}")
    print(f"  TMB range: {min(tmb_by_patient.values())} - {max(tmb_by_patient.values())}")
    
    # 4. Sample mapping
    print("\n[3/6] Mapping samples...")
    sample_table = pd.read_csv(os.path.join(DIR, 'SampleTableCorrected.9.19.16.csv'), index_col=0)
    # SampleTable index is Pt101_Pre, PatientID column has Pt101
    sample_to_patient = dict(zip(sample_table.index, sample_table['PatientID']))
    
    # Keep only pre-treatment samples
    pre_samples = []
    for col in expr.columns:
        # Extract Pt101_Pre
        short_name = '_'.join(col.split('_')[:2])
        if short_name in sample_to_patient and short_name.endswith('_Pre'):
            pre_samples.append((col, short_name))
    
    print(f"  Pre-treatment samples: {len(pre_samples)}")
    
    # Construct analysis data
    rows = []
    for full_name, short_name in pre_samples:
        patient = sample_to_patient.get(short_name)
        if patient is None or patient not in clin.index:
            continue
        
        row = {
            'sample': full_name,
            'patient': patient,
            'ZP3_FPKM': zp3_expr[full_name],
        }
        
        # Clinical
        c = clin.loc[patient]
        row['SubtypeEZ'] = c.get('SubtypeEZ', np.nan)
        row['BOR'] = c.get('BOR', np.nan)
        row['myBOR'] = c.get('myBOR', np.nan)
        row['PFS'] = c.get('PFS', np.nan)
        row['PFSWK'] = c.get('PFSWK', np.nan)
        row['OS'] = c.get('OS', np.nan)
        row['Cohort'] = c.get('Cohort', np.nan)
        
        # TMB
        row['TMB_nonsyn'] = tmb_by_patient.get(patient, 0)
        
        # BRAF/NRAS
        row['BRAF_mut'] = patient in braf_patients
        row['NRAS_mut'] = patient in nras_patients
        row['BRAF_variant'] = '; '.join(braf_mut.get(patient, []))
        row['NRAS_variant'] = '; '.join(nras_mut.get(patient, []))
        
        # Response binary classification (using myBOR: PRCR=response, PD/SD=non-response)
        if row['myBOR'] == 'PRCR':
            row['response_binary'] = 1  # Response
        elif row['myBOR'] in ('PD', 'SD'):
            row['response_binary'] = 0  # Non-response
        else:
            row['response_binary'] = np.nan
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    print(f"  Analyzing data: {len(df)} samples")
    print(f"  ZP3 range: {df['ZP3_FPKM'].min():.2f} - {df['ZP3_FPKM'].max():.2f}")
    
    return df, clin

def analyze_subtype_zp3(df):
    """Analysis 1: ZP3 expression by SubtypeEZ"""
    print("\n[4/6] Stratifying ZP3 by SubtypeEZ...")
    
    # Merge BRAF/NRAS subtype categories
    df = df.copy()
    df['mut_type'] = df.apply(lambda r: 
        'BRAF' if r['BRAF_mut'] else 
        'NRAS' if r['NRAS_mut'] else 
        'NF1' if r['SubtypeEZ'] == 'NF1' else 
        'TripleWT' if r['SubtypeEZ'] == 'TripleWt' else
        'Unknown', axis=1)
    
    subtypes = [s for s in ['BRAF', 'NRAS', 'NF1', 'TripleWT'] if s in df['mut_type'].values]
    
    print(f"  Mutation type distribution: {df['mut_type'].value_counts().to_dict()}")
    
    # Statistics
    stats_rows = []
    for st in subtypes:
        vals = df[df['mut_type'] == st]['ZP3_FPKM']
        stats_rows.append({
            'subtype': st, 'n': len(vals), 
            'mean': vals.mean(), 'median': vals.median(), 'std': vals.std()
        })
    
    # Test: BRAF vs TripleWT
    if 'BRAF' in subtypes and 'TripleWT' in subtypes:
        u, p = stats.mannwhitneyu(
            df[df['mut_type']=='BRAF']['ZP3_FPKM'],
            df[df['mut_type']=='TripleWT']['ZP3_FPKM']
        )
        print(f"  BRAF vs TripleWT: U={u:.0f}, p={p:.3f}")
    
    return df

def analyze_tmb_zp3(df):
    """Analysis 2: TMB correlation with ZP3"""
    print("\n[5/6] TMB-ZP3 correlation analysis...")
    
    valid = df.dropna(subset=['TMB_nonsyn', 'ZP3_FPKM'])
    if len(valid) < 10:
        print("  Insufficient samples")
        return
    
    r, p = stats.spearmanr(valid['TMB_nonsyn'], valid['ZP3_FPKM'])
    print(f"  TMB vs ZP3: ρ={r:.3f}, p={p:.3f}")
    
    r_p, p_p = stats.pearsonr(valid['TMB_nonsyn'], valid['ZP3_FPKM'])
    print(f"  TMB vs ZP3 (Pearson): r={r_p:.3f}, p={p_p:.3f}")
    
    # Stratify ZP3 response by high/low TMB
    median_tmb = valid['TMB_nonsyn'].median()
    df['TMB_group'] = np.where(df['TMB_nonsyn'] > median_tmb, 'High_TMB', 'Low_TMB')
    
    return df

def stratified_response_analysis(df):
    """Analysis 3: ZP3-response relationship stratified by subtype/TMB"""
    print("\n[6/6] Stratified analysis of ZP3-response relationship...")
    
    results = {}
    
    valid_resp = df.dropna(subset=['response_binary', 'ZP3_FPKM'])
    print(f"  Samples with response data: {len(valid_resp)}")
    
    # Overall
    resp = valid_resp[valid_resp['response_binary'] == 1]
    nonresp = valid_resp[valid_resp['response_binary'] == 0]
    if len(resp) > 3 and len(nonresp) > 3:
        u, p = stats.mannwhitneyu(resp['ZP3_FPKM'], nonresp['ZP3_FPKM'], alternative='two-sided')
        r_p, p_p = stats.pearsonr(valid_resp['ZP3_FPKM'], valid_resp['response_binary'])
        results['Overall'] = {
            'n_resp': len(resp), 'n_nonresp': len(nonresp),
            'resp_median': resp['ZP3_FPKM'].median(),
            'nonresp_median': nonresp['ZP3_FPKM'].median(),
            'MWU_p': p, 'pearson_r': r_p, 'pearson_p': p_p
        }
        print(f"  Overall: R(n={len(resp)}) med={resp['ZP3_FPKM'].median():.2f} vs NR(n={len(nonresp)}) med={nonresp['ZP3_FPKM'].median():.2f}, p={p:.4f}")
    
    # Stratified by mutation type
    for st in ['BRAF', 'NRAS', 'TripleWT']:
        sub = valid_resp[valid_resp['mut_type'] == st]
        r_sub = sub[sub['response_binary'] == 1]
        nr_sub = sub[sub['response_binary'] == 0]
        if len(r_sub) >= 2 and len(nr_sub) >= 2:
            u, p = stats.mannwhitneyu(r_sub['ZP3_FPKM'], nr_sub['ZP3_FPKM'], alternative='two-sided')
            results[f'{st}_stratified'] = {
                'n_resp': len(r_sub), 'n_nonresp': len(nr_sub),
                'resp_median': r_sub['ZP3_FPKM'].median(),
                'nonresp_median': nr_sub['ZP3_FPKM'].median(),
                'MWU_p': p
            }
            print(f"  {st}: R(n={len(r_sub)}) med={r_sub['ZP3_FPKM'].median():.2f} vs NR(n={len(nr_sub)}) med={nr_sub['ZP3_FPKM'].median():.2f}, p={p:.4f}")
    
    # Stratified by TMB
    if 'TMB_group' in df.columns:
        for g in ['High_TMB', 'Low_TMB']:
            sub = valid_resp[valid_resp['TMB_group'] == g]
            r_sub = sub[sub['response_binary'] == 1]
            nr_sub = sub[sub['response_binary'] == 0]
            if len(r_sub) >= 2 and len(nr_sub) >= 2:
                u, p = stats.mannwhitneyu(r_sub['ZP3_FPKM'], nr_sub['ZP3_FPKM'], alternative='two-sided')
                results[f'{g}_stratified'] = {
                    'n_resp': len(r_sub), 'n_nonresp': len(nr_sub),
                    'resp_median': r_sub['ZP3_FPKM'].median(),
                    'nonresp_median': nr_sub['ZP3_FPKM'].median(),
                    'MWU_p': p
                }
                print(f"  {g}: R(n={len(r_sub)}) med={r_sub['ZP3_FPKM'].median():.2f} vs NR(n={len(nr_sub)}) med={nr_sub['ZP3_FPKM'].median():.2f}, p={p:.4f}")
    
    return results, valid_resp

def multivariate_logistic(df):
    """Analysis 4: Multivariate logistic regression"""
    print("\n[Extra] Multivariate logistic regression...")
    
    valid = df.dropna(subset=['response_binary', 'ZP3_FPKM', 'TMB_nonsyn'])
    if len(valid) < 15:
        print("  Insufficient samples")
        return None
    
    try:
        import statsmodels.api as sm
        
        # Features
        X = pd.DataFrame({
            'ZP3_FPKM': valid['ZP3_FPKM'],
            'TMB_nonsyn': valid['TMB_nonsyn'],
            'BRAF_mut': valid['BRAF_mut'].astype(int),
            'NRAS_mut': valid['NRAS_mut'].astype(int),
        })
        y = valid['response_binary']
        
        X_const = sm.add_constant(X)
        model = sm.Logit(y, X_const).fit(disp=0)
        
        print(f"  Logit model (n={len(valid)}):")
        for name, coef, p in zip(model.params.index, model.params, model.pvalues):
            print(f"    {name}: coef={coef:.3f}, OR={np.exp(coef):.2f}, p={p:.4f}")
        
        return model
    except Exception as e:
        print(f"  Regression failed: {e}")
        return None

def plot_analysis(df, valid_resp, results, model):
    """Comprehensive visualization"""
    print("\nGenerating visualization...")
    
    fig = plt.figure(figsize=(18, 14))
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.3)
    
    # A: ZP3 by mutation subtype
    ax1 = fig.add_subplot(gs[0, 0])
    subtypes = ['BRAF', 'NRAS', 'NF1', 'TripleWT']
    present = [s for s in subtypes if s in df['mut_type'].values]
    data_by_type = [df[df['mut_type']==s]['ZP3_FPKM'] for s in present]
    
    bp = ax1.boxplot(data_by_type, patch_artist=True, widths=0.5,
                     flierprops={'alpha': 0.4, 'markersize': 4})
    colors_box = ['#E74C3C', '#9B59B6', '#F39C12', '#3498DB'][:len(present)]
    for patch, c in zip(bp['boxes'], colors_box):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    ax1.set_xticklabels(present, fontsize=10)
    ax1.set_ylabel('ZP3 FPKM', fontsize=11)
    ax1.set_title('A. ZP3 by Mutation Subtype', fontsize=12, fontweight='bold')
    
    # B: ZP3 vs TMB scatter
    ax2 = fig.add_subplot(gs[0, 1])
    for st in present:
        sub = df[df['mut_type']==st]
        ax2.scatter(sub['TMB_nonsyn'], sub['ZP3_FPKM'], label=st, alpha=0.6, s=40)
    valid = df.dropna(subset=['TMB_nonsyn', 'ZP3_FPKM'])
    if len(valid) > 5:
        r, p = stats.spearmanr(valid['TMB_nonsyn'], valid['ZP3_FPKM'])
        ax2.set_title(f'B. ZP3 vs TMB (ρ={r:.2f}, p={p:.3f})', fontsize=12, fontweight='bold')
    else:
        ax2.set_title('B. ZP3 vs TMB', fontsize=12, fontweight='bold')
    ax2.set_xlabel('TMB (nonsynonymous mutations)', fontsize=11)
    ax2.set_ylabel('ZP3 FPKM', fontsize=11)
    ax2.legend(fontsize=8)
    
    # C: ZP3-response overall
    ax3 = fig.add_subplot(gs[0, 2])
    resp = valid_resp[valid_resp['response_binary']==1]
    nonresp = valid_resp[valid_resp['response_binary']==0]
    bp2 = ax3.boxplot([resp['ZP3_FPKM'], nonresp['ZP3_FPKM']], 
                      patch_artist=True, widths=0.5,
                      flierprops={'alpha':0.4, 'markersize':4})
    for patch, c in zip(bp2['boxes'], ['#2ECC71', '#E74C3C']):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    ax3.set_xticklabels([f'Response\n(n={len(resp)})', f'Non-response\n(n={len(nonresp)})'], fontsize=10)
    ax3.set_ylabel('ZP3 FPKM', fontsize=11)
    if 'Overall' in results:
        p_val = results['Overall']['MWU_p']
        ax3.set_title(f'C. ZP3 by Response (p={p_val:.3f})', fontsize=12, fontweight='bold')
    else:
        ax3.set_title('C. ZP3 by Response', fontsize=12, fontweight='bold')
    
    # D: ZP3-response stratified by mutation type (bar chart)
    ax4 = fig.add_subplot(gs[1, :2])
    strat_keys = [k for k in results if k.endswith('_stratified')]
    strat_labels = [k.replace('_stratified', '') for k in strat_keys]
    
    x = np.arange(len(strat_keys))
    width = 0.35
    resp_meds = [results[k]['resp_median'] for k in strat_keys]
    nonresp_meds = [results[k]['nonresp_median'] for k in strat_keys]
    
    ax4.bar(x - width/2, resp_meds, width, label='Response', color='#2ECC71', alpha=0.8)
    ax4.bar(x + width/2, nonresp_meds, width, label='Non-response', color='#E74C3C', alpha=0.8)
    ax4.set_xticks(x)
    ax4.set_xticklabels(strat_labels, rotation=15, ha='right', fontsize=10)
    ax4.set_ylabel('ZP3 FPKM (median)', fontsize=11)
    ax4.set_title('D. ZP3-Response Stratified by Subtype/TMB', fontsize=12, fontweight='bold')
    ax4.legend()
    
    # E: Multivariate regression forest plot
    ax5 = fig.add_subplot(gs[1, 2])
    if model is not None:
        names = ['ZP3', 'TMB', 'BRAF', 'NRAS']
        ors = [np.exp(model.params[n]) for n in ['ZP3_FPKM', 'TMB_nonsyn', 'BRAF_mut', 'NRAS_mut']]
        pvals = [model.pvalues[n] for n in ['ZP3_FPKM', 'TMB_nonsyn', 'BRAF_mut', 'NRAS_mut']]
        
        y_pos = np.arange(len(names))[::-1]
        colors = ['#E74C3C' if p < 0.05 else '#7F8C8D' for p in pvals]
        ax5.barh(y_pos, np.log2(ors), color=colors, alpha=0.8)
        ax5.set_yticks(y_pos)
        ax5.set_yticklabels(names, fontsize=10)
        ax5.axvline(0, color='gray', linestyle='--', linewidth=1)
        ax5.set_xlabel('log2(OR)', fontsize=11)
        ax5.set_title('E. Multivariate Logistic (Response)', fontsize=12, fontweight='bold')
        for yi, (name, or_, pv) in enumerate(zip(names, ors, pvals)):
            ax5.text(np.log2(or_) + 0.05, y_pos[yi], f'OR={or_:.2f}\np={pv:.3f}', 
                    va='center', fontsize=8)
    else:
        ax5.axis('off')
        ax5.text(0.5, 0.5, 'Model not available\n(sample size too small)', 
                ha='center', va='center', fontsize=12)
    
    # F: TMB-stratified ZP3-response
    ax6 = fig.add_subplot(gs[2, :])
    ax6.axis('off')
    text_lines = ["Stratified ZP3-Response Analysis Summary", "=" * 55]
    for k in sorted(results.keys()):
        r_ = results[k]
        text_lines.append(
            f"  {k:18s}: R(n={r_['n_resp']}) med={r_['resp_median']:.2f} | "
            f"NR(n={r_['n_nonresp']}) med={r_['nonresp_median']:.2f} | p={r_['MWU_p']:.4f}"
        )
    if model is not None:
        text_lines.append("")
        text_lines.append("Multivariate Logistic Regression (response ~ ZP3 + TMB + BRAF + NRAS):")
        for name in ['ZP3_FPKM', 'TMB_nonsyn', 'BRAF_mut', 'NRAS_mut']:
            text_lines.append(
                f"  {name:12s}: OR={np.exp(model.params[name]):.2f} "
                f"(95% CI [{np.exp(model.conf_int().loc[name,0]):.2f}, "
                f"{np.exp(model.conf_int().loc[name,1]):.2f}]), p={model.pvalues[name]:.3f}"
            )
    ax6.text(0.03, 0.95, "\n".join(text_lines), transform=ax6.transAxes,
            fontsize=10, va='top', fontfamily='monospace')
    ax6.set_title('F. Summary', fontsize=12, fontweight='bold')
    
    fig_path = os.path.join(OUTPUT_DIR, "fig_gse91061_subtype_tmb_zp3.png")
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    print(f"  Figure: {fig_path}")
    plt.close()
    return fig_path

def save_outputs(df, results, model):
    """Save results"""
    df.to_csv(os.path.join(OUTPUT_DIR, "gse91061_subtype_tmb_merged.csv"), index=False)
    
    with open(os.path.join(OUTPUT_DIR, "gse91061_stratified_results.json"), "w") as f:
        json.dump(results, f, indent=2, default=str)
    
    if model is not None:
        summary = pd.DataFrame({
            'feature': model.params.index,
            'coef': model.params.values,
            'OR': np.exp(model.params.values),
            'pvalue': model.pvalues.values
        })
        summary.to_csv(os.path.join(OUTPUT_DIR, "gse91061_logistic_regression.csv"), index=False)
    
    print(f"\n✓ Analysis complete!")
    print(f"  Output directory: {OUTPUT_DIR}")

def main():
    df, clin = load_data()
    df = analyze_subtype_zp3(df)
    df = analyze_tmb_zp3(df)
    results, valid_resp = stratified_response_analysis(df)
    model = multivariate_logistic(df)
    fig_path = plot_analysis(df, valid_resp, results, model)
    save_outputs(df, results, model)

if __name__ == "__main__":
    main()
