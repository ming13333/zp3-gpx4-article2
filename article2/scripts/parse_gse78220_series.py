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
解析 GSE78220 series matrix 文件，提取临床信息
"""
import os
import gzip
import re
import pandas as pd

OUT_DIR = os.path.join(ROOT, "output", "h2_bulk")
series_file = os.path.join(OUT_DIR, 'GSE78220_series_matrix.txt.gz')

print('解析 series matrix 文件...')
lines = []
with gzip.open(series_file, 'rt') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            break
        if line.startswith('!'):
            lines.append(line.rstrip())

print(f'读取 {len(lines)} 行元数据')

# 提取字段
fields = {}
for line in lines:
    m = re.match(r'!(\w+)\s+(.*)', line, re.S)
    if not m:
        continue
    key = m.group(1)
    # 按制表符分割，但注意引号内的制表符
    # 简单处理：按制表符分割，然后去掉引号
    cols = m.group(2).split('\t')
    cols = [c.strip('"') for c in cols]
    fields[key] = cols

print(f'提取到 {len(fields)} 个字段')

# 查看关键字段
key_fields = ['Sample_description', 'Sample_characteristics_ch1', 'Sample_characteristics_ch2', 
              'Sample_characteristics_ch3', 'Sample_characteristics_ch4', 'Sample_characteristics_ch5']
for field in key_fields:
    if field in fields:
        values = fields[field]
        print(f'\n{field} (长度 {len(values)}):')
        # 显示前几个值
        for i, v in enumerate(values[:5]):
            print(f'  {i}: {v}')
        if len(values) > 5:
            print(f'  ... (共 {len(values)} 个)')

# 提取样本描述和细胞类型
if 'Sample_description' in fields:
    descriptions = fields['Sample_description']
    print(f'\n样本描述 (前10个):')
    for i, desc in enumerate(descriptions[:10]):
        print(f'  {i}: {desc}')

# 提取特征
print('\n所有特征字段:')
for key, values in fields.items():
    if key.startswith('Sample_characteristics'):
        print(f'{key}: {len(values)} 个值')
        # 显示唯一值
        unique_vals = list(set(values))
        if len(unique_vals) <= 10:
            print(f'  唯一值: {unique_vals}')
        else:
            print(f'  唯一值数量: {len(unique_vals)}')
            print(f'  示例: {unique_vals[:5]}')

# 尝试提取生存信息
print('\n尝试提取生存信息...')
survival_keywords = ['survival', 'os', 'overall', 'time', 'month', 'day', 'year', 'status', 'event', 'alive', 'dead', 'death']
for key, values in fields.items():
    if any(keyword in key.lower() for keyword in survival_keywords):
        print(f'{key}: {values[:5]}...')

# 保存解析结果
parsed_csv = os.path.join(OUT_DIR, 'gse78220_series_matrix_parsed.csv')
parsed_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in fields.items()]))
parsed_df.to_csv(parsed_csv, index=False)
print(f'\n解析结果已保存: {parsed_csv}')

# 提取样本注释表
print('\n构建样本注释表...')
if 'Sample_description' in fields and 'Sample_characteristics_ch1' in fields:
    # 假设每个字段的长度相同
    n_samples = len(fields['Sample_description'])
    print(f'样本数: {n_samples}')
    
    # 构建注释表
    annotation_data = {}
    for key, values in fields.items():
        if len(values) == n_samples:
            annotation_data[key] = values
    
    annotation_df = pd.DataFrame(annotation_data)
    print(f'注释表形状: {annotation_df.shape}')
    print(f'注释表列: {list(annotation_df.columns)}')
    
    # 保存注释表
    annot_csv = os.path.join(OUT_DIR, 'gse78220_annotation.csv')
    annotation_df.to_csv(annot_csv, index=False)
    print(f'注释表已保存: {annot_csv}')
    
    # 显示前几行
    print('\n注释表预览:')
    print(annotation_df.head().to_string())
else:
    print('无法构建注释表：字段长度不匹配')