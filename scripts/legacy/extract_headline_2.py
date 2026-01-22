import arviz as az
from pathlib import Path
import numpy as np
run = Path('results/bayesian_v3_6/3_1_1/both')
run = sorted([p for p in run.iterdir() if p.is_dir()])[-1]
idata = az.from_netcdf(run / 'trace_v3_6.nc')
print(az.summary(idata, var_names=['beta_xi','xi_nm_acute','xi_nm_chronic','xi_nm_control'], kind='stats'))
xac = idata.posterior['xi_nm_acute'].values.ravel()
xch = idata.posterior['xi_nm_chronic'].values.ravel()
print('P(xi_acute < xi_chronic) =', float((xac < xch).mean()))
