from pathlib import Path

import pandas as pd

root = Path('data/extracted_expanded/data_ratios_comparison')
src = root / 'bayesian_inputs_3_1_1.csv'
out = root / 'bayesian_inputs_3_1_1_NAAcr_only.csv'

df = pd.read_csv(src)
mask = df['Metabolite'].astype(str).str.replace(' ', '').str.upper().eq('NAA/CR')
naa = df[mask].copy()
naa = naa.rename(columns={'Mean': 'NAA_mean', 'SE': 'NAA_SE', 'Study': 'study'})
naa['Phase'] = naa['Phase'].fillna('Control')
naa.to_csv(out, index=False)
print('Wrote:', out)
