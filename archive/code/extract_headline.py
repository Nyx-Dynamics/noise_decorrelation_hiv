import arviz as az
from pathlib import Path
import numpy as np

# Point this at your chosen run (replace if you’ve picked a newer one)
run = Path('results/bayesian_v3_6/3_1_1/both/20251125T155428Z_f39492cb')
idata = az.from_netcdf(run / 'trace_v3_6.nc')

# 1) Summary of key parameters (use the actual variable names)
print(az.summary(
    idata,
    var_names=['beta_xi','xi_nm_acute','xi_nm_chronic','xi_nm_control'],
    kind='stats'
))

# 2) Probability xi_acute < xi_chronic
xac = idata.posterior['xi_nm_acute'].values.ravel()
xch = idata.posterior['xi_nm_chronic'].values.ravel()
print('P(xi_acute < xi_chronic) =', float((xac < xch).mean()))

# Optional: confirm chains/ESS/R-hat for these variables
print(az.summary(idata, var_names=['beta_xi','xi_nm_acute','xi_nm_chronic'], kind='diagnostics'))
