import arviz as az
from pathlib import Path
base = Path('results/bayesian_v3_6/3_1_1/both')
run = sorted([p for p in base.iterdir() if p.is_dir()])[-1]
idata = az.from_netcdf(run / 'trace_v3_6.nc')
print('Run:', run.name)
print(az.summary(idata, var_names=['beta_xi','xi_acute','xi_chronic'], kind='stats'))
# Example: compute P(xi_acute < xi_chronic)
xia = idata.posterior['xi_acute'].values.ravel()
xic = idata.posterior['xi_chronic'].values.ravel()
import numpy as np
p = float(np.mean(xia < xic))
print('P(xi_acute < xi_chronic) =', p)