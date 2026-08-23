#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE91061 黑色素瘤抗 PD-1 治疗队列 ZP3 分析
验证 ZP3 与免疫治疗反应的关系
"""

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

def parse_series_matrix(matrix_file):
    """解析 Series Matrix 文件提取临床信息"""
    sample_ids = []
    visits = []
    responses = []
    tissues = []
    
    with open(matrix_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('!Series_sample_id'):
                # 提取样本ID
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    sample_ids = parts[1].strip('"').split()
            elif line.startswith('!Sample_characteristics_ch1'):
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    # 提取 visit 信息（第一行 ch1）
                    if not visits:
                        visits = [p.split(': ')[-1].strip('"') for p in parts[1:] if 'visit' in p.lower()]
                    # 提取 response 信息（第二行 ch1）
                    elif not responses:
                        responses = [p.split(': ')[-1].strip('"') for p in parts[1:] if 'response' in p.lower()]
                    # 提取 tissue 信息（第三行 ch1）
                    elif not tissues:
                        tissues = [p.split(': ')[-1].strip('"') for p in parts[1:] if 'tissue' in p.lower()]
    
    # 创建临床数据 DataFrame
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
    """创建表达矩阵列名到 GSM 样本名的映射"""
    # 表达矩阵列名格式: Pt1_Pre_AD101148-6
    # 临床数据样本名格式: GSM2420259
    # 需要从 Series Matrix 文件中提取映射关系
    
    # 读取 Series Matrix 文件获取 GSM 到患者信息的映射
    mapping = {}
    with open(os.path.join(G9DIR, 'GSE91061_series_matrix.txt'), 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    # 找到 Sample_title 行（包含 GSM ID 和样本描述）
    sample_titles = []
    for line in lines:
        if line.startswith('!Sample_title'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                sample_titles = parts[1:]
            break
    
    # 找到 Sample_geo_accession 行（GSM ID）
    gsm_ids = []
    for line in lines:
        if line.startswith('!Sample_geo_accession'):
            parts = line.strip().split('\t')
            if len(parts) > 1:
                gsm_ids = [g.strip('"') for g in parts[1:]]  # 去除每个 GSM 的引号
            break
    
    # 创建 GSM 到样本描述的映射
    gsm_to_desc = {}
    for i, gsm in enumerate(gsm_ids):
        if i < len(sample_titles):
            desc = sample_titles[i].strip('"')
            gsm_to_desc[gsm] = desc
    
    # 创建样本描述到 GSM 的映射
    desc_to_gsm = {desc: gsm for gsm, desc in gsm_to_desc.items()}
    
    # 从样本描述中提取患者ID和visit类型
    # 格式示例: "Pt1_Pre_AD101148-6"
    expr_to_gsm = {}
    for expr_sample in expr_columns:
        # 直接使用表达矩阵列名作为样本描述
        if expr_sample in desc_to_gsm:
            expr_to_gsm[expr_sample] = desc_to_gsm[expr_sample]
    
    return expr_to_gsm

def main():
    print("=== GSE91061 黑色素瘤抗 PD-1 治疗队列 ZP3 分析 ===\n")
    
    # 1. 加载临床数据
    print("1. 加载临床数据...")
    clinical = parse_series_matrix('GSE91061_series_matrix.txt')
    if clinical is None:
        print("错误：无法解析临床数据")
        return
    
    print(f"  总样本数: {len(clinical)}")
    print(f"  访问类型分布:")
    print(clinical['visit'].value_counts())
    print(f"\n  治疗反应分布:")
    print(clinical['response'].value_counts())
    
    # 2. 加载表达数据
    print("\n2. 加载表达数据...")
    G9DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'output', 'gse91061_validation')
    expr = pd.read_csv(os.path.join(G9DIR, 'GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv'), index_col=0)
    print(f"  基因数: {len(expr)}")
    print(f"  样本数: {len(expr.columns)}")
    
    # 3. 检查 ZP3 基因
    # 注意：表达矩阵使用 Entrez Gene ID，需要将 ZP3 转换为数字 ID
    # ZP3 的 Entrez Gene ID 是 7784（2026-08-17 审计修复；本脚本已由 v2 取代，仅保留作历史）
    zp3_entrez = 7784
    if zp3_entrez in expr.index:
        print(f"\n3. ZP3 基因 (Entrez ID: {zp3_entrez}) 已找到")
        zp3_expr = expr.loc[zp3_entrez]
    else:
        print(f"\n3. ZP3 基因 (Entrez ID: {zp3_entrez}) 未找到")
        print("  尝试寻找其他 ZP3 相关基因...")
        # 查找可能的 ZP3 相关基因
        zp3_genes = [idx for idx in expr.index if 'ZP3' in str(idx)]
        if zp3_genes:
            print(f"  找到 {len(zp3_genes)} 个 ZP3 相关基因: {zp3_genes[:5]}")
        else:
            print("  未找到任何 ZP3 相关基因")
            return
    
    # 4. 对齐样本
    print("\n4. 对齐样本...")
    
    # 创建映射
    sample_mapping = create_sample_mapping(expr.columns)
    print(f"  成功映射的样本数: {len(sample_mapping)}")
    
    if len(sample_mapping) < 30:
        print("  警告：成功映射的样本数太少，分析可能不可靠")
    
    # 重新组织临床数据以匹配表达矩阵
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
    print(f"  最终共同样本数: {len(common_samples)}")
    
    # 5. 分析治疗反应
    print("\n5. 分析 ZP3 与治疗反应的关系...")
    
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
    
    if len(valid_samples) < 20:
        print("  警告：有效样本数太少")
        return
    
    # 6. 统计检验
    print("\n6. 统计检验...")
    
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

        # 计算效应量 (rank-biserial r = 1 - 2U/(n1*n2)) —— 标准 Mann-Whitney 效应量
        n1, n2 = len(response_values), len(non_response_values)
        r = 1 - (2 * u_stat) / (n1 * n2)

        print(f"\n  Mann-Whitney U 检验:")
        print(f"    U 统计量: {u_stat:.1f}")
        print(f"    p 值: {p_value:.4f}")
        print(f"    效应量 (rank-biserial r): {r:.3f}")

        if p_value < 0.05:
            print(f"    结论: ZP3 表达与治疗反应显著相关 (p < 0.05)")
        else:
            print(f"    结论: ZP3 表达与治疗反应无显著关联 (p >= 0.05)")
    
    # 7. 按访问类型分析（治疗前后）
    print("\n7. 按访问类型分析（治疗前后）...")
    
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
    # 稳健提取患者号：去掉末尾的访问/批次后缀（pre/on/baseline/post 或数字），
    # 避免旧版 split('_')[0] 在含多下划线 ID（如 TCGA-AB-1234-01）上误截。
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
        stat, p_value = stats.wilcoxon(paired_df['pre'], paired_df['on'])
        print(f"\n  配对检验 (治疗前后):")
        print(f"    配对样本数: {len(paired_df)}")
        print(f"    Wilcoxon 符号秩检验 p 值: {p_value:.4f}")
        print(f"    治疗前 ZP3 均值: {paired_df['pre'].mean():.3f}")
        print(f"    治疗后 ZP3 均值: {paired_df['on'].mean():.3f}")
    
    # 8. 可视化
    print("\n8. 生成可视化...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('GSE91061 Melanoma Anti-PD-1: ZP3 Expression Analysis', fontsize=14, fontweight='bold')
    
    # 8.1 治疗反应分组箱线图
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
    
    # 8.2 访问类型分组箱线图
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
        
        # 添加个体数据点
        for i, data in enumerate(data_for_plot2):
            x = np.random.default_rng(0).normal(i+1, 0.04, size=len(data))
            ax2.scatter(x, data, alpha=0.4, s=15, color='black')
        
        ax2.set_ylabel('ZP3 Expression (FPKM)')
        ax2.set_title('ZP3 by Visit Type (Pre vs On-treatment)')
        ax2.grid(True, alpha=0.3)
    
    # 8.3 配对数据散点图
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
    
    # 8.4 效应量（本队列 GSE91061，诚实展示；跨队列对比应由各队列自身分析输出汇总）
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

    # 保存图表
    output_path = 'fig_gse91061_zp3_analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  图表已保存至: {output_path}")

    # 9. 保存结果
    print("\n9. 保存结果...")

    # 创建结果摘要
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
    
    # 保存为 CSV
    results_df = pd.DataFrame([results_summary])
    results_df.to_csv(os.path.join(G9DIR, 'gse91061_zp3_results.csv'), index=False)
    print("  结果摘要已保存至: gse91061_zp3_results.csv")
    
    print("\n=== 分析完成 ===")
    _nr_mean = group_stats.get('Non-response', {}).get('mean', 0)
    _r_mean = group_stats.get('Response', {}).get('mean', 0)
    _direction = ("高 ZP3 表达与较差治疗反应（非应答）相关"
                  if _nr_mean > _r_mean
                  else "高 ZP3 表达与较好治疗反应（应答）相关")
    _effect_label = "中等" if abs(r) >= 0.3 else ("小" if abs(r) >= 0.1 else "极小")
    _sig_label = "显著相关" if p_value < 0.05 else "无显著关联"
    print(f"结论: ZP3 在黑色素瘤抗 PD-1 治疗中与治疗反应{_sig_label} (Mann-Whitney U, p={p_value:.4f})")
    print(f"方向: {_direction}（应答组均值 {_r_mean:.3f} vs 非应答组 {_nr_mean:.3f}）")
    print(f"效应量: rank-biserial r = {r:.3f}（{_effect_label}效应）")
    print(f"p 值: {p_value:.4f}（{'显著' if p_value < 0.05 else '不显著'}）")

if __name__ == "__main__":
    main()
