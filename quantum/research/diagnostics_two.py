import arviz as az
from pathlib import Path
import numpy as np
base = Path('results/bayesian_v3_6/3_1_1/both')
run = sorted([p for p in base.iterdir() if p.is_dir()])[-1]
idata = az.from_netcdf(run / 'trace_v3_6.nc')
print('Run:', run.name)
vars = ['beta_xi','xi_nm_control','xi_nm_acute','xi_nm_chronic']
print(az.summary(idata, var_names=vars, kind='stats'))
xictrl = idata.posterior['xi_nm_control'].values.ravel()
xia = idata.posterior['xi_nm_acute'].values.ravel()
xic = idata.posterior['xi_nm_chronic'].values.ravel()
print('P(xi_nm_acute < xi_nm_control) =', float((xia < xictrl).mean()))
print('P(xi_nm_acute < xi_nm_chronic) =', float((xia < xic).mean()))