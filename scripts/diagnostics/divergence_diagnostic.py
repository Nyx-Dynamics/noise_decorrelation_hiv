import arviz as az
from pathlib import Path
import numpy as np
run = Path('results/bayesian_v3_6/3_1_1/both/20251125T155428Z_f39492cb')
idata = az.from_netcdf(run / 'trace_v3_6.nc')
print('BFMI per chain:', az.bfmi(idata))
print(az.summary(idata, var_names=['beta_xi','xi_nm_acute','xi_nm_chronic','xi_nm_control','study_sd','se_scale','era_effect'], kind='diagnostics'))
