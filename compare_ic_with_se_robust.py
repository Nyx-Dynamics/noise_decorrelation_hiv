#!/usr/bin/env python3
import argparse
from pathlib import Path
import arviz as az
import xarray as xr
import pandas as pd

PREF_CHRONIC = [
    'NAA_ratio_chronic_obs',                  # individual chronic
    'NAA_Young_BG_Chronic_mean_obs',          # young BG chronic
    'NAA_Sailasuta_BG_Chronic_mean_obs',      # sailasuta BG chronic
]
PREF_CONTROL = [
    'NAA_ratio_control_obs',                  # individual control
    'NAA_Young_BG_Control_mean_obs',          # young BG control
    'NAA_Sailasuta_BG_Control_mean_obs',      # sailasuta BG control
]

# Utility: pick ONE name present in both sets from a preference list
def pick_one_common(pref_list, set_a, set_b):
    for name in pref_list:
        if name in set_a and name in set_b:
            return name
    return None

# Fallback: pick first common var whose name contains the phase token
def pick_by_substring(token, set_a, set_b):
    cand = sorted(set_a.intersection(set_b))
    for v in cand:
        if token.lower() in v.lower():
            return v
    return None

# Build a combined InferenceData with a single log_likelihood from a list of var names
def combine_llk(idata, use_vars):
    if not use_vars:
        raise RuntimeError('No variables specified to combine.')
    llk = idata.log_likelihood
    total = None
    for v in use_vars:
        if v not in llk:
            raise RuntimeError(f'Requested LL var not found in trace: {v}')
        arr = llk[v]
        total = arr if total is None else (total + arr)
    new_llk = xr.Dataset({'log_likelihood': total})
    out = az.InferenceData(**{g: getattr(idata, g) for g in idata.groups() if g != 'log_likelihood'})
    out.add_groups({'log_likelihood': new_llk})
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Compare IC with SE (robust common set, chronic+control fallback)')
    ap.add_argument('--with-trace', required=True, help='Path to with‑Valcour .nc')
    ap.add_argument('--no-trace',   required=True, help='Path to no‑Valcour .nc')
    args = ap.parse_args()

    p_with = Path(args.with_trace)
    p_no   = Path(args.no_trace)

    id_with = az.from_netcdf(str(p_with))
    id_no   = az.from_netcdf(str(p_no))
if not hasattr(id_with, 'log_likelihood') or id_with.log_likelihood is None:
        raise RuntimeError(f'{p_with} has no log_likelihood; re-run with idata_kwargs={"log_likelihood": True}.')
if not hasattr(id_no, 'log_likelihood') or id_no.log_likelihood is None:
        raise RuntimeError(f'{p_no} has no log_likelihood; re-run with idata_kwargs={"log_likelihood": True}.')

set_with = set(id_with.log_likelihood.data_vars)
set_no   = set(id_no.log_likelihood.data_vars)

print('Available LL vars (WITH Valcour):')
print(sorted(set_with))
print('\nAvailable LL vars (NO Valcour):')
print(sorted(set_no))

# Try to select ONE chronic and ONE control var present in BOTH traces
chronic = pick_one_common(PREF_CHRONIC, set_with, set_no) or pick_by_substring('Chronic', set_with, set_no)
control = pick_one_common(PREF_CONTROL, set_with, set_no) or pick_by_substring('Control', set_with, set_no)

if not chronic or not control:
        raise RuntimeError('Could not find a usable common set (need at least one chronic and one control LL var present in BOTH traces).')

print(f"\nChosen common chronic: {chronic}")
print(f"Chosen common control: {control}")

id_with_c = combine_llk(id_with, [chronic, control])
id_no_c   = combine_llk(id_no,   [chronic, control])

# Per‑model IC (chronic+control)
loo_with, loo_no   = az.loo(id_with_c), az.loo(id_no_c)
waic_with, waic_no = az.waic(id_with_c), az.waic(id_no_c)

per_model = pd.DataFrame([
        {'model':'with_valcour', 'elpd_loo': float(loo_with.elpd_loo), 'p_loo': float(loo_with.p_loo), 'LOO': float(-2*loo_with.elpd_loo),
                                 'elpd_waic': float(waic_with.elpd_waic), 'p_waic': float(waic_with.p_waic), 'WAIC': float(-2*waic_with.elpd_waic)},
        {'model':'no_valcour',   'elpd_loo': float(loo_no.elpd_loo),   'p_loo': float(loo_no.p_loo),   'LOO': float(-2*loo_no.elpd_loo),
                                 'elpd_waic': float(waic_no.elpd_waic),   'p_waic': float(waic_no.p_waic),   'WAIC': float(-2*waic_no.elpd_waic)}
    ]).set_index('model')

print('\n=== Per‑model IC (common chronic+control) ===')
print(per_model)

# Δ with SE (now n_points identical → compare works)
comp_loo  = az.compare({'with_valcour': id_with_c, 'no_valcour': id_no_c}, ic='loo')
comp_waic = az.compare({'with_valcour': id_with_c, 'no_valcour': id_no_c}, ic='waic')

print('\n=== Pairwise comparison with SE (PSIS‑LOO) ===')
print(comp_loo)
print('\n=== Pairwise comparison with SE (WAIC) ===')
print(comp_waic)
