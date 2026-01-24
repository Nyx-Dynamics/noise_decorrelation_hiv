from quantum.data_loaders import load_bayesian_inputs

ratio = '3_1_1'
df, _ = load_bayesian_inputs(ratio)
print('Counts by art_era:\n', df['art_era'].value_counts(dropna=False))
print('\nUnique studies in pre_modern:', df[df['art_era']=='pre_modern']['study'].nunique() if 'study' in df.columns else 'n/a')
print('Rows in pre_modern:', (df['art_era']=='pre_modern').sum())