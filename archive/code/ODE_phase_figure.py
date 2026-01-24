import pandas as pd, matplotlib.pyplot as plt
from pathlib import Path
v4dir = sorted((Path('results/enzyme_v4/3_1_1/both')).iterdir())[-1]
pred = pd.read_csv(v4dir / 'predictions.csv')
fig, ax = plt.subplots(1,2, figsize=(10,4))
ax[0].bar(pred['Phase'], pred['NAA_Cr']); ax[0].set_title('ODE NAA/Cr (3:1:1)')
ax[1].bar(pred['Phase'], pred['Cho_Cr']); ax[1].set_title('ODE Cho/Cr (3:1:1)')
Path('figures/figures').mkdir(parents=True, exist_ok=True)
plt.tight_layout(); plt.savefig('figures/figures/ode_phase_ratios_3_1_1_both.png', dpi=200, bbox_inches='tight')
print('✓ Wrote figures/figures/ode_phase_ratios_3_1_1_both.png')
