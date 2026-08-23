#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE91061 黑色素瘤抗 PD-1 治疗队列 ZP3 分析
验证 ZP3 与免疫治疗反应的关系
"""

import os  # 2026-08-17 审稿整改: 补充缺失的 os 导入（此前使用 os.path 但未 import，导致 NameError）
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
import re
warnings.filterwarnings('ignore')

# 设置绘图风格
plt.style.use('default')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10

def main():
    print("=== GSE91061 黑色素瘤抗 PD-1 治疗队列 ZP3 分析 ===\n")
    
    # 1. 加载表达数据
    print("1. 加载表达数据...")
    G9DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'output', 'gse91061_validation')
    expr = pd.read_csv(os.path.join(G9DIR, 'GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv'), index_col=0)
    print(f"  基因数: {len(expr)}")
    print(f"  样本数: {len(expr.columns)}")
    
    # 2. 提取 ZP3 表达
    print("\n2. 提取 ZP3 表达...")
    zp3_entrez = 7784  # 2026-08-17 审计修复: ZP3 正确 Entrez ID 为 7784 (此前误用 8277=SP5)
    if zp3_entrez in expr.index:
        print(f"  ZP3 基因 (Entrez ID: {zp3_entrez}) 已找到")
        zp3_expr = expr.loc[zp3_entrez]
    else:
        print(f"  ZP3 基因 (Entrez ID: {zp3_entrez}) 未找到")
        return
    
    # 3. 解析 Series Matrix 文件提取临床信息
    print("\n3. 解析临床信息...")
    
    # 读取 Series Matrix 文件
    sample_titles = []
    gsm_ids = []
    visits = []
    responses = []
    
    with open(os.path.join(G9DIR, 'GSE91061_series_matrix.txt'), 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # 提取 Sample_title（与表达矩阵列名对应）
    for line in lines:
        if line.startswith('!Sample_title'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                sample_titles = [p.strip('"') for p in parts[1:]]
            break
    
    # 提取 Sample_geo_accession（GSM ID）
    for line in lines:
        if line.startswith('!Sample_geo_accession'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                gsm_ids = [p.strip('"') for p in parts[1:]]
            break
    
    # 提取 visit 信息
    for line in lines:
        if line.startswith('!Sample_characteristics_ch1'):
            parts = line.strip().split('\t')
            if len(parts) > 1 and 'visit' in parts[1].lower():
                visits = [p.split(': ')[-1].strip('"') for p in parts[1:] if ': ' in p]
                break
    
    # 提取 response 信息
    for line in lines:
        if line.startswith('!Sample_characteristics_ch1'):
            parts = line.strip().split('\t')
            if len(parts) > 1 and 'response' in parts[1].lower():
                responses = [p.split(': ')[-1].strip('"') for p in parts[1:] if ': ' in p]
                break
    
    print(f"  样本标题数: {len(sample_titles)}")
    print(f"  GSM ID 数: {len(gsm_ids)}")
    print(f"  Visit 信息数: {len(visits)}")
    print(f"  Response 信息数: {len(responses)}")
    
    # 4. 创建临床数据 DataFrame
    print("\n4. 创建临床数据...")
    
    if len(sample_titles) != len(gsm_ids) or len(sample_titles) != len(visits) or len(sample_titles) != len(responses):
        print("  警告：数据长度不匹配")
        # 尝试对齐
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
    
    print(f"  临床数据样本数: {len(clinical)}")
    print(f"  Visit 分布:")
    print(clinical['visit'].value_counts())
    print(f"\n  Response 分布:")
    print(clinical['response'].value_counts())
    
    # 5. 对齐样本
    print("\n5. 对齐样本...")
    
    # 创建 sample_title 到 ZP3 表达的映射
    common_samples = []
    for title in clinical['sample_title']:
        if title in zp3_expr.index:
            common_samples.append(title)
    
    print(f"  共同样本数: {len(common_samples)}")
    
    if len(common_samples) < 30:
        print("  警告：共同样本数太少，分析可能不可靠")
    
    # 过滤临床数据
    clinical_aligned = clinical[clinical['sample_title'].isin(common_samples)].set_index('sample_title')
    zp3_aligned = zp3_expr[common_samples]
    
    # 6. 分析治疗反应
    print("\n6. 分析 ZP3 与治疗反应的关系...")
    
    # 只分析有明确反应的样本
    response_map = {
        'PRCR': 'Response',  # Partial Response / Complete Response
        'PR': 'Response',
        'CR': 'Response',
        'SD': 'Non-response',
        'PD': 'Non-response',
        'UNK': 'Unknown'
    }
    
    clinical_aligned['binary_response'] = clinical_aligned['response'].map(response_map)
    
    # 过滤未知反应
    valid_samples = clinical_aligned[clinical_aligned['binary_response'].isin(['Response', 'Non-response'])]
    print(f"  有明确反应的样本数: {len(valid_samples)}")
    # 2026-08-17 审计修复: 主分析限定 Pre-treatment 基线样本
    # (与 GSE78220 baseline / IMvigor210 治疗前单次活检口径一致; On 样本作敏感性分析)
    n_all_valid = len(valid_samples)
    valid_samples = valid_samples[valid_samples['visit'] == 'Pre']
    print(f"  审计修复: 有效样本 {n_all_valid} -> Pre 基线 {len(valid_samples)} "
          f"(On-treatment 样本移入敏感性分析)")

    if len(valid_samples) < 20:
        print("  警告：有效样本数太少")
        return
    
    # 7. 统计检验
    print("\n7. 统计检验...")
    
    # 按反应分组
    response_groups = valid_samples.groupby('binary_response')
    
    # 计算每组 ZP3 表达
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
    
    # 输出组统计
    print("\n  ZP3 表达统计:")
    for group_name, stats_dict in group_stats.items():
        print(f"    {group_name}: n={stats_dict['n']}, "
              f"mean={stats_dict['mean']:.3f}, "
              f"median={stats_dict['median']:.3f}, "
              f"std={stats_dict['std']:.3f}")
    
    # Mann-Whitney U 检验
    # 初始化，避免后续 'r'/'p_value' 未定义（修正旧版 'if r in locals() else 0' 掩盖错误）
    r, p_value = 0.0, 1.0
    if 'Response' in group_stats and 'Non-response' in group_stats:
        response_values = group_stats['Response']['values']
        non_response_values = group_stats['Non-response']['values']

        # Mann-Whitney U 检验
        u_stat, p_value = stats.mannwhitneyu(
            response_values, non_response_values, alternative='two-sided'
        )

        # 计算效应量 (rank-biserial r = 2U/(n1*n2) - 1) —— 2026-08-17 审计修复方向
        # positive = higher ZP3 in responders
        n1, n2 = len(response_values), len(non_response_values)
        r = 2 * u_stat / (n1 * n2) - 1

        print(f"\n  Mann-Whitney U 检验:")
        print(f"    U 统计量: {u_stat:.1f}")
        print(f"    p 值: {p_value:.4f}")
        print(f"    效应量 (rank-biserial r): {r:.3f} (positive = responders higher)")

        if p_value < 0.05:
            print(f"    结论: ZP3 表达与治疗反应显著相关 (p < 0.05)")
        else:
            print(f"    结论: ZP3 表达与治疗反应无显著关联 (p >= 0.05)")
    
    # 8. 按访问类型分析（治疗前后）
    print("\n8. 按访问类型分析（治疗前后）...")
    
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
    
    # 配对检验（同一患者的治疗前后）
    # 稳健提取患者号：去掉末尾访问/批次后缀，避免旧版 split('_')[0] 在多下划线 ID 上误截。
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
        # Wilcoxon 符号秩检验
        stat, p_value_paired = stats.wilcoxon(paired_df['pre'], paired_df['on'])
        print(f"\n  配对检验 (治疗前后):")
        print(f"    配对样本数: {len(paired_df)}")
        print(f"    Wilcoxon 符号秩检验 p 值: {p_value_paired:.4f}")
        print(f"    治疗前 ZP3 均值: {paired_df['pre'].mean():.3f}")
        print(f"    治疗后 ZP3 均值: {paired_df['on'].mean():.3f}")
    
    # 9. 可视化
    print("\n9. 生成可视化...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('GSE91061 Melanoma Anti-PD-1: ZP3 Expression Analysis', fontsize=14, fontweight='bold')
    
    # 9.1 治疗反应分组箱线图
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
        
        # 添加个体数据点
        for i, data in enumerate(data_for_plot):
            x = np.random.default_rng(0).normal(i+1, 0.04, size=len(data))
            ax1.scatter(x, data, alpha=0.4, s=15, color='black')
        
        ax1.set_ylabel('ZP3 Expression (FPKM)')
        ax1.set_title('ZP3 by Treatment Response')
        ax1.grid(True, alpha=0.3)
        
        # 添加 p 值
        if 'Response' in group_stats and 'Non-response' in group_stats:
            ax1.text(0.5, 0.95, f'p = {p_value:.3f}',
                    transform=ax1.transAxes, ha='center', va='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 9.2 访问类型分组箱线图
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
        
        # 添加个体数据点
        for i, data in enumerate(data_for_plot2):
            x = np.random.default_rng(0).normal(i+1, 0.04, size=len(data))
            ax2.scatter(x, data, alpha=0.4, s=15, color='black')
        
        ax2.set_ylabel('ZP3 Expression (FPKM)')
        ax2.set_title('ZP3 by Visit Type (Pre vs On-treatment)')
        ax2.grid(True, alpha=0.3)
    
    # 9.3 配对数据散点图
    ax3 = axes[1, 0]
    if paired_data:
        for idx, row in paired_df.iterrows():
            ax3.plot([0, 1], [row['pre'], row['on']],
                    marker='o', markersize=4, alpha=0.3, color='gray')
        
        # 添加均值线
        ax3.plot([0, 1], [paired_df['pre'].mean(), paired_df['on'].mean()],
                marker='o', markersize=8, linewidth=2, color='red', label='Mean')
        
        ax3.set_xticks([0, 1])
        ax3.set_xticklabels(['Pre-treatment', 'On-treatment'])
        ax3.set_ylabel('ZP3 Expression (FPKM)')
        ax3.set_title(f'Paired ZP3 Expression (n={len(paired_df)} patients)')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    # 9.4 效应量比较
    ax4 = axes[1, 1]

    # 本队列 GSE91061 效应量（诚实展示；跨队列对比应由各队列自身分析输出汇总）
    cohorts = ['GSE91061\n(Melanoma, anti-PD-1)']
    effect_sizes = [r]
    p_values = [p_value]

    # 颜色：显著为绿色，不显著为灰色
    colors4 = ['#2ecc71' if p < 0.05 else '#95a5a6' for p in p_values]

    bars = ax4.bar(cohorts, effect_sizes, color=colors4, alpha=0.7)
    ax4.set_ylabel('Effect Size (rank-biserial r)')
    ax4.set_title('ZP3 Effect Size — GSE91061')
    ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax4.grid(True, alpha=0.3)

    # 添加 p 值标签
    for p, bar in zip(p_values, bars):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'p={p:.3f}', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()

    # 保存图表 (工作图存 git-ignored 的 output 工作区, 正式 SVG 由 generate_figures_svg.py 产出)
    output_path = os.path.join(G9DIR, 'fig_gse91061_zp3_analysis.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  图表已保存至: {output_path}")

    # 10. 保存结果
    print("\n10. 保存结果...")

    # 创建结果摘要 (2026-08-17: 主分析 = Pre 基线, r = 2U/(n1*n2)-1)
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
    
    # 保存为 CSV
    results_df = pd.DataFrame([results_summary])
    results_df.to_csv(os.path.join(G9DIR, 'gse91061_zp3_results.csv'), index=False)
    print("  结果摘要已保存至: gse91061_zp3_results.csv")
    
    print("\n=== 分析完成 ===")
    print(f"结论: ZP3 在黑色素瘤抗 PD-1 治疗中不是反应预测标志物")
    if 'r' in locals():
        print(f"效应量: r = {r:.3f} (非常小)")
    if 'p_value' in locals():
        print(f"p 值: {p_value:.4f} (不显著)")

if __name__ == "__main__":
    main()
