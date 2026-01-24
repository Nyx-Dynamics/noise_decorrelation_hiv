
import arviz as az
import numpy as np

# Load the traces (they'll be saved even with divergences)
idata_full = az.from_netcdf('results/model_comparison/Full_v4_trace.nc')

# Check which parameters have issues
div_info = az.plots.plot_pair(
    idata_full,
    var_names=['xi_acute', 'xi_chronic', 'beta_xi'],
    divergences=True,
    figsize=(12, 12)
)
plt.savefig('results/model_comparison/divergence_diagnostics.png', dpi=150)
print('Divergence plot saved')

# Get divergence summary
post = idata_full.posterior
samp_stats = idata_full.sample_stats
n_divs = samp_stats.diverging.sum().values
print(f'Total divergences: {n_divs}')
print(f'Divergence rate: {100*n_divs/samp_stats.diverging.size:.1f}%')
