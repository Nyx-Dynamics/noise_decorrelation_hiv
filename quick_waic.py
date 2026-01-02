# Save as: quick_waic_loo.py and run: python quick_waic_loo.py
import arviz as az
import xarray as xr
from pathlib import Path
import pandas as pd

with_path = Path('quantum/quantum/results_v3_6/runs/with_valcour_waic/trace_with_valcour.nc')
no_path   = Path('quantum/quantum/results_v3_6/runs/no_valcour_waic/trace_no_valcour.nc')

# Load traces
idata_with = az.from_netcdf(str(with_path))
idata_no   = az.from_netcdf(str(no_path))

# Select core NAA likelihoods only (keeps memory low)
core_vars = ['NAA_ratio_acute_obs','NAA_ratio_chronic_obs','NAA_ratio_control_obs']

def make_combined_llk(idata):
    if not hasattr(idata, 'log_likelihood') or idata.log_likelihood is None:
        raise RuntimeError('No log_likelihood present; re-run model with idata_kwargs={"log_likelihood": True}.')
    llk = idata.log_likelihood
    use = [v for v in core_vars if v in llk.data_vars]
    if not use:
        # Fallback: sum everything (may use more memory)
        use = list(llk.data_vars)
    total = None
    for v in use:
        arr = llk[v]
        total = arr if total is None else (total + arr)
    new_llk = xr.Dataset({'log_likelihood': total})
    # Rebuild idata with a single combined log_likelihood
    out = az.InferenceData(**{g: getattr(idata, g) for g in idata.groups() if g != 'log_likelihood'})
    out.add_groups({'log_likelihood': new_llk})
    return out
idata_with2 = make_combined_llk(idata_with)
idata_no2   = make_combined_llk(idata_no)

# Compute WAIC and PSIS-LOO
waic_with = az.waic(idata_with2)
waic_no   = az.waic(idata_no2)
loo_with  = az.loo(idata_with2)
loo_no    = az.loo(idata_no2)

rows = []
rows.append({'model':'with_valcour', 'waic':waic_with.waic, 'p_waic':waic_with.p_waic, 'elpd_waic':waic_with.elpd_waic,
             'loo':loo_with.loo, 'p_loo':loo_with.p_loo, 'elpd_loo':loo_with.elpd_loo})
rows.append({'model':'no_valcour', 'waic':waic_no.waic, 'p_waic':waic_no.p_waic, 'elpd_waic':waic_no.elpd_waic,
             'loo':loo_no.loo, 'p_loo':loo_no.p_loo, 'elpd_loo':loo_no.elpd_loo})

df = pd.DataFrame(rows).set_index('model')
print('\nWAIC/LOO per model (core NAA terms only):')
print(df)

# Simple pairwise comparison (lower WAIC/LOO is better)
print('\nΔWAIC = WAIC(no_valcour) - WAIC(with_valcour):', float(df.loc['no_valcour','waic'] - df.loc['with_valcour','waic']))
print('ΔLOO  = LOO(no_valcour)  - LOO(with_valcour):', float(df.loc['no_valcour','loo']  - df.loc['with_valcour','loo']))


# Faarr = llk[v]
#         total = arr if total is None else (total + arr)
#     new_llk = xr.Dataset({'log_likelihood': total})
#     # Rebuild idata with a single combined log_likelihood
#     out = az.InferenceData(**{g: getattr(idata, g) for g in idata.groups() if g != 'log_likelihood'})
#     out.add_groups({'log_likelihood': new_llk})
#     return out
#
# idata_with2 = make_combined_llk(idata_with)
# idata_no2   = make_combined_llk(idata_no)
#
# # Compute WAIC and PSIS-LOO
# waic_with = az.waic(idata_with2)llback: sum ever