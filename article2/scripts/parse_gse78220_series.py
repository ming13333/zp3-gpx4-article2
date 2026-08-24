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
Parse GSE78220 series matrix file, extract clinical information
"""
import os
import gzip
import re
import pandas as pd

OUT_DIR = os.path.join(ROOT, "output", "h2_bulk")
series_file = os.path.join(OUT_DIR, 'GSE78220_series_matrix.txt.gz')

print('Parsing series matrix file...')
lines = []
with gzip.open(series_file, 'rt') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            break
        if line.startswith('!'):
            lines.append(line.rstrip())

print(f'Read {len(lines)} metadata lines')

# Extract fields
fields = {}
for line in lines:
    m = re.match(r'!(\w+)\s+(.*)', line, re.S)
    if not m:
        continue
    key = m.group(1)
    # Split by tab, but note tabs inside quotes
    # Simple handling: split by tab, then remove quotes
    cols = m.group(2).split('\t')
    cols = [c.strip('"') for c in cols]
    fields[key] = cols

print(f'Extracted {len(fields)} fields')

# Inspect key fields
key_fields = ['Sample_description', 'Sample_characteristics_ch1', 'Sample_characteristics_ch2', 
              'Sample_characteristics_ch3', 'Sample_characteristics_ch4', 'Sample_characteristics_ch5']
for field in key_fields:
    if field in fields:
        values = fields[field]
        print(f'\n{field} (length {len(values)}):')
        # Show first few values
        for i, v in enumerate(values[:5]):
            print(f'  {i}: {v}')
        if len(values) > 5:
            print(f'  ... ({len(values)} total)')

# Extract sample descriptions and cell types
if 'Sample_description' in fields:
    descriptions = fields['Sample_description']
    print(f'\nSample descriptions (first 10):')
    for i, desc in enumerate(descriptions[:10]):
        print(f'  {i}: {desc}')

# Extract features
print('\nAll feature fields:')
for key, values in fields.items():
    if key.startswith('Sample_characteristics'):
        print(f'{key}: {len(values)} values')
        # Show unique values
        unique_vals = list(set(values))
        if len(unique_vals) <= 10:
            print(f'  Unique values: {unique_vals}')
        else:
            print(f'  Number of unique values: {len(unique_vals)}')
            print(f'  Examples: {unique_vals[:5]}')

# Try to extract survival information
print('\nTrying to extract survival information...')
survival_keywords = ['survival', 'os', 'overall', 'time', 'month', 'day', 'year', 'status', 'event', 'alive', 'dead', 'death']
for key, values in fields.items():
    if any(keyword in key.lower() for keyword in survival_keywords):
        print(f'{key}: {values[:5]}...')

# Save parsed results
parsed_csv = os.path.join(OUT_DIR, 'gse78220_series_matrix_parsed.csv')
parsed_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in fields.items()]))
parsed_df.to_csv(parsed_csv, index=False)
print(f'\nParsed results saved: {parsed_csv}')

# Extract sample annotation table
print('\nBuilding sample annotation table...')
if 'Sample_description' in fields and 'Sample_characteristics_ch1' in fields:
    # Assume all fields have the same length
    n_samples = len(fields['Sample_description'])
    print(f'Number of samples: {n_samples}')
    
    # Build annotation table
    annotation_data = {}
    for key, values in fields.items():
        if len(values) == n_samples:
            annotation_data[key] = values
    
    annotation_df = pd.DataFrame(annotation_data)
    print(f'Annotation table shape: {annotation_df.shape}')
    print(f'Annotation table columns: {list(annotation_df.columns)}')
    
    # Save annotation table
    annot_csv = os.path.join(OUT_DIR, 'gse78220_annotation.csv')
    annotation_df.to_csv(annot_csv, index=False)
    print(f'Annotation table saved: {annot_csv}')
    
    # Show first few rows
    print('\nAnnotation table preview:')
    print(annotation_df.head().to_string())
else:
    print('Unable to build annotation table: field lengths do not match')
