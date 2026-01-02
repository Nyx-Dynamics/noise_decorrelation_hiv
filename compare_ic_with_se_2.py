#!/usr/bin/env python3
import argparse
from pathlib import Path
import arviz as az
import xarray as xr
import pandas as pd

# Strict common-set comparison: prefer Young BG (Acute/Chronic/Control) group means
YOUNG_SET = [
    'NAA_Young_BG_Acute_mean_obs',
    'NAA_Young_BG_Chronic_mean_obs',
    'NAA_Young_BG_Control_mean_obs',
]
SAIL_SET = [
    'NAA_Sailasuta_BG_Acute_mean_obs',
    'NAA_Sailasuta_BG_Chronic_mean_obs',
    'NAA_Sailasuta_BG_Control_mean_obs',
]

REQUIRED_SETS = [YOUNG_SET, SAIL_SET]
REQUIRED_SETS = [YOUNG_SET, SAIL_SET]

def pick_common_vars(idata_with, idata_no):
    llw = set(idata_with.log_likelihood.data_vars)
    lln = set(idata_no.log_likelihood.data_vars)
    for candidate in REQUIRED_SETS:
        if all(v in llw for v in candidate) and all(v in lln for v in candidate):
            return candidate
    raise RuntimeError("Could not find a strict common set of acute/chronic/control group-mean variables in both traces.\n"
                       "Expected Young BG or Sailasuta BG group means to be present in BOTH traces.")

def build_common_idata(idata, use_vars):
    llk = idata.log_likelihood
    total = None
    for v in use_vars:
        arr = llk[v]
        total = arr if total is None else (total + arr)
    new_llk = xr.Dataset({'log_likelihood': total})
    out = az.InferenceData(**{g: getattr(idata, g) for g in idata.groups() if g != 'log_likelihood'})
    out.add_groups({'log_likelihood': new_llk})
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Compare IC with SE (strict common set)')
    ap.add_argument('--with-trace', required=True, help='Path to with‑Valcour .nc')
    ap.add_argument('--no-trace',   required=True, help='Path to no‑Valcour .nc')
    args = ap.parse_args()

    p_with = Path(args.with_trace)
    p_no   = Path(args.no_trace)

    id_with = az.from_netcdf(str(p_with))
    id_no   = az.from_netcdf(str(p_no))
    if not hasattr(id_with, 'log_likelihood') or id_with.log_likelihood is None:
        raise RuntimeError(f'{p_with} has no log_likelihood; re-run sampling with idata_kwargs={"log_likelihood": True}.')
    if not hasattr(id_no, 'log_likelihood') or id_no.log_likelihood is None:
        raise RuntimeError(f'{p_no} has no log_likelihood; re-run sampling with idata_kwargs={"log_likelihood": True}.')

    use_vars = pick_common_vars(id_with, id_no)
    print('Using common LL vars:', use_vars)
id_with_c = build_common_idata(id_with, use_vars)
id_no_c   = build_common_idata(id_no, use_vars)

# Per-model IC
loo_with, loo_no   = az.loo(id_with_c), az.loo(id_no_c)
waic_with, waic_no = az.waic(id_with_c), az.waic(id_no_c)
df = pd.DataFrame([
        {'model':'with_valcour', 'elpd_loo': float(loo_with.elpd_loo), 'p_loo': float(loo_with.p_loo), 'LOO': float(-2*loo_with.elpd_loo),
                                      'elpd_waic': float(waic_with.elpd_waic), 'p_waic': float(waic_with.p_waic), 'WAIC': float(-2*waic_with.elpd_waic)},
        {'model':'no_valcour',   'elpd_loo': float(loo_no.elpd_loo),   'p_loo': float(loo_no.p_loo),   'LOO': float(-2*loo_no.elpd_loo),
                                      'elpd_waic': float(waic_no.elpd_waic),   'p_waic': float(waic_no.p_waic),   'WAIC': float(-2*waic_no.elpd_waic)}
    ]).set_index('model')
print('\n=== Per-model IC (common Young/Sailasuta BG means) ===')
print(df)

# Pairwise Δ with SE (now n_points is identical → compare works)
comp_loo  = az.compare({'with_valcour': id_with_c, 'no_valcour': id_no_c}, ic='loo')
comp_waic = az.compare({'with_valcour': id_with_c, 'no_valcour': id_no_c}, ic='waic')
print('\n=== Pairwise comparison with SE (PSIS‑LOO) ===')
print(comp_loo)
print('\n=== Pairwise comparison with SE (WAIC) ===')
print(comp_waic)