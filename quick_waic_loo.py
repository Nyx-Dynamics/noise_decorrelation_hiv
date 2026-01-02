import arviz as az
import xarray as xr
from pathlib import Path
import pandas as pd

with_path = Path('quantum/quantum/results_v3_6/runs/with_valcour_waic/trace_with_valcour.nc')
no_path = Path('quantum/quantum/results_v3_6/runs/no_valcour_waic/trace_no_valcour.nc')

# Load traces
idata_with = az.from_netcdf(str(with_path))
idata_no = az.from_netcdf(str(no_path))

# Use core NAA likelihoods (keeps memory low)
CORE_VARS = ['NAA_ratio_acute_obs', 'NAA_ratio_chronic_obs', 'NAA_ratio_control_obs']


def combine_llk_core(idata):
    if not hasattr(idata, 'log_likelihood') or idata.log_likelihood is None:
        raise RuntimeError('No log_likelihood present; re-run model with idata_kwargs={"log_likelihood": True}.')
    llk = idata.log_likelihood
    use = [v for v in CORE_VARS if v in llk.data_vars]
    if not use:
        # Fallback: sum everything (may be heavier)
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


idata_with2 = combine_llk_core(idata_with)
idata_no2 = combine_llk_core(idata_no)

# Compute WAIC and PSIS-LOO (ELPDData objects)
waic_with = az.waic(idata_with2)
waic_no = az.waic(idata_no2)
loo_with = az.loo(idata_with2)
loo_no = az.loo(idata_no2)
# Convert to conventional WAIC/LOO numbers
row_with = {
    'model': 'with_valcour',
    'elpd_waic': float(waic_with.elpd_waic),
    'p_waic': float(waic_with.p_waic),
    'waic': float(-2.0 * waic_with.elpd_waic),
    'elpd_loo': float(loo_with.elpd_loo),
    'p_loo': float(loo_with.p_loo),
    'loo': float(-2.0 * loo_with.elpd_loo),
}
row_no = {
    'model': 'no_valcour',
    'elpd_waic': float(waic_no.elpd_waic),
    'p_waic': float(waic_no.p_waic),
    'waic': float(-2.0 * waic_no.elpd_waic),
    'elpd_loo': float(loo_no.elpd_loo),
    'p_loo': float(loo_no.p_loo),
    'loo': float(-2.0 * loo_no.elpd_loo),
}

df = pd.DataFrame([row_with, row_no]).set_index('model')
print('\nWAIC/LOO per model (core NAA terms only):')
print(df)

# Pairwise deltas (lower is better)
print('\nΔWAIC = WAIC(no_valcour) - WAIC(with_valcour):', df.loc['no_valcour','waic'] - df.loc['with_valcour','waic'])
print('ΔLOO  = LOO(no_valcour)  - LOO(with_valcour):', df.loc['no_valcour','loo']  - df.loc['with_valcour','loo'])
