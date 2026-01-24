import arviz as az
from pathlib import Path
import numpy as np
run = Path('results/bayesian_v3_6/3_1_1/both')
run = sorted([p for p in run.iterdir() if p.is_dir()])[-1]
idata = az.from_netcdf(run / 'trace_v3_6.nc')
print('Run:', run.name)
print('Total divergences:', int(idata.sample_stats['diverging'].sum()))
print('BFMI per chain:', az.bfmi(idata))
print(az.summary(idata, var_names=['beta_xi','xi_nm_acute','xi_nm_chronic','xi_nm_control','study_sd','se_scale','era_effect'], kind='diagnostics'))
