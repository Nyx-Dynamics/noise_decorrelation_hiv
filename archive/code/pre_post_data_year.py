from pathlib import Path

import numpy as np
import pandas as pd

# Path to your ratio comparison input (adjust if needed)
p = Path('data/extracted_expanded/data_ratios_comparison/bayesian_inputs_3_1_1.csv')
df = pd.read_csv(p)

# Ensure publication_year exists and is numeric when possible
def to_num(s):
    return pd.to_numeric(s, errors='coerce')

if 'publication_year' not in df.columns:
    # Create a column filled with NaN (float) so we can coalesce from other cols
    df['publication_year'] = np.nan
else:
    df['publication_year'] = to_num(df['publication_year'])

# Try to coalesce from any plausible year columns that might exist
candidates = [
    'publication_year', 'PublicationYear', 'pub_year',
    'measurement_year', 'scan_year', 'Year', 'YearPublished'
]

# Start with NaN (float), fill from candidates where present
year_numeric = pd.Series(np.nan, index=df.index, dtype='float64')
for c in candidates:
    if c in df.columns:
        year_numeric = year_numeric.fillna(to_num(df[c]))

# If still missing, map from
year_map = {
    'Chang': 2002,
    'Chang_2002': 2002,
    'Sailasuta_2012': 2012,
    'Sailasuta et al. 2016': 2016,
    'Young_2014': 2014,
    'Mohamed et al. 2010': 2010,
    'Valcour_2015': 2015,
}
if 'Study' in df.columns:
    year_numeric = year_numeric.fillna(pd.Series(df['Study']).map(year_map))

# Commit merged year back to publication_year (numeric)
df['publication_year'] = to_num(year_numeric)

# Compute art_era using 2006/2007 split (NaN → 'unknown')
cutoff = 2006
art = pd.Series('unknown', index=df.index, dtype='object')
art[(df['publication_year'].notna()) & (df['publication_year'] <= cutoff)] = 'pre_modern'
art[(df['publication_year'].notna()) & (df['publication_year'] >= cutoff + 1)] = 'post_modern'

df['art_era'] = art
df['art_era_idx'] = df['art_era'].map({'pre_modern': 0, 'post_modern': 1, 'unknown': -1})

print('Counts by art_era:\n', df['art_era'].value_counts(dropna=False))
print('Rows in pre_modern:', (df['art_era'] == 'pre_modern').sum())
print('Rows in post_modern:', (df['art_era'] == 'post_modern').sum())
print('Rows unknown:', (df['art_era'] == 'unknown').sum())

# Optional: backup the CSV once
backup = p.with_suffix('.backup.csv')
if not backup.exists():
    df.to_csv(backup, index=False)

# Write back
df.to_csv(p, index=False)
print('Updated CSV written to', p)