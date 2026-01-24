from pathlib import Path

import arviz as az

run_dir = Path('results/bayesian_v3_6/3_1_1/both/20251123T223826Z_a8c794cc')
idata = az.from_netcdf(run_dir / 'trace_v3_6.nc')
print('max R-hat:', az.rhat(idata).to_array().max().values)
print('min ESS bulk:', az.ess(idata, method='bulk').to_array().min().values)
print('min ESS tail:', az.ess(idata, method='tail').to_array().min().values)
# Optional: energy diagnostic plot saved next to the run for review
import matplotlib.pyplot as plt
az.plot_energy(idata)
plt.savefig(run_dir / 'energy_diagnostic.png', dpi=200, bbox_inches='tight'); plt.close()
