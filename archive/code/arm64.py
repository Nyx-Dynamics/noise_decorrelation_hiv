python - <<'PY'
from pathlib import Path

import pandas as pd

root = Path('data/extracted_expanded/data_ratios_comparison')
src = root / 'bayesian_inputs_3_1_1.csv'
out = root / 'bayesian_inputs_3_1_1_NAAcr_only.csv'

# Load and filter
(df := pd.read_csv(src))
naa = df[df['Metabolite'].str.upper().str.replace(' ', '') == 'NAA/CR'].copy()
# Standardize required columns
naa = naa.rename(columns={'Mean': 'NAA_mean', 'SE': 'NAA_SE', 'Study':'study'})
# Optional: ensure Phase is present and consistent
naa['Phase'] = naa['Phase'].fillna('Control')
naa.to_csv(out, index=False)
print('Wrote:', out)

