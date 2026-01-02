#!/usr/bin/env python3
import argparse
from pathlib import Path
import arviz as az
import xarray as xr
import pandas as pd

# Harmonized core logic: acute via individuals; if absent, fall back to acute group means
CORE_FALLBACKS = {
    'acute':   {'primary': ['NAA_ratio_acute_obs'],
                'fallback': lambda vars_available: [v for v in vars_available if v.startswith('NAA_') and v.endswith('_Acute_mean_obs')]},
    'chronic': {'primary': ['NAA_ratio_chronic_obs'],
                'fallback': lambda vars_available: [v for v in vars_available if v.startswith('NAA_') and v.endswith('_Chronic_mean_obs')]},
    'control': {'primary': ['NAA_ratio_control_obs'],
                'fallback': lambda vars_available: [v for v in vars_available if v.startswith('NAA_') and v.endswith('_Control_mean_obs')]},
}
def pick_core_vars(llk):
    vars_available = set(llk.data_vars)
    chosen = []
    for phase in ['acute','chronic','control']:
        prim = [v for v in CORE_FALLBACKS[phase]['primary'] if v in vars_available]
        if prim:
            chosen.extend(prim)
        else:
            fb = CORE_FALLBACKS[phase]['fallback'](vars_available)
            chosen.extend(sorted(fb))
    seen, ordered = set(), []
    for v in chosen:
        if v not in seen:
            ordered.append(v); seen.add(v)
    return ordered
def to_core_idata(nc_path: Path) -> az.InferenceData:
    idata = az.from_netcdf(str(nc_path))
    if not hasattr(idata, 'log_likelihood') or idata.log_likelihood is None:
        raise RuntimeError(f'{nc_path} has no log_likelihood; re‑run sampling with idata_kwargs={"log_likelihood": True}.')
    llk = idata.log_likelihood
    use = pick_core_vars(llk)
    if not use:
        raise RuntimeError(f'Could not determine core vars in {nc_path.name}.')
    total = None
    for v in use:
        arr = llk[v]
        total = arr if total is None else (total + arr)
    new_llk = xr.Dataset({'log_likelihood': total})
    id2 = az.InferenceData(**{g: getattr(idata, g) for g in idata.groups() if g != 'log_likelihood'})
    id2.add_groups({'log_likelihood': new_llk})
    return id2

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Compare IC with SE (core, harmonized)')
    ap.add_argument('--with-trace', required=True, help='Path to with‑Valcour .nc')
    ap.add_argument('--no-trace',   required=True, help='Path to no‑Valcour .nc')
args = ap.parse_args()
p_with = Path(args.with_trace)
p_no   = Path(args.no_trace)

id_with = to_core_idata(p_with)
id_no   = to_core_idata(p_no)

loo_with, loo_no = az.loo(id_with), az.loo(id_no)
waic_with, waic_no = az.waic(id_with), az.waic(id_no)

comp_loo  = az.compare({'with_valcour': id_with, 'no_valcour': id_no}, ic='loo')
comp_waic = az.compare({'with_valcour': id_with, 'no_valcour': id_no}, ic='waic')

print('\n=== Per‑model IC (core) ===')
df = pd.DataFrame([
        {'model':'with_valcour', 'elpd_loo': float(loo_with.elpd_loo), 'p_loo': float(loo_with.p_loo), 'LOO': float(-2*loo_with.elpd_loo),
                                      'elpd_waic': float(waic_with.elpd_waic), 'p_waic': float(waic_with.p_waic), 'WAIC': float(-2*waic_with.elpd_waic)},
        {'model':'no_valcour',   'elpd_loo': float(loo_no.elpd_loo),   'p_loo': float(loo_no.p_loo),   'LOO': float(-2*loo_no.elpd_loo),
                                      'elpd_waic': float(waic_no.elpd_waic),   'p_waic': float(waic_no.p_waic),   'WAIC': float(-2*waic_no.elpd_waic)}
    ]).set_index('model')
print(df)
print('\n=== Pairwise comparison with SE (PSIS‑LOO) ===')
print(comp_loo)
print('\n=== Pairwise comparison with SE (WAIC) ===')
print(comp_waic)

dL = float(loo_no.elpd_loo - loo_with.elpd_loo)
dW = float(waic_no.elpd_waic - waic_with.elpd_waic)
try:
    se_loo = float(comp_loo.loc['no_valcour','se']) if 'se' in comp_loo.columns else float('nan')
except Exception:
    se_loo = float('nan')
    try:
        se_waic = float(comp_waic.loc['no_valcour','se']) if 'se' in comp_waic.columns else float('nan')
    except Exception:
        se_waic = float('nan')

    print(f"\nΔ elpd_loo (no - with) = {dL:.3f} ± {se_loo:.3f}  (higher is better)")
    print(f"Δ elpd_waic(no - with) = {dW:.3f}  (higher is better; WAIC = -2*eelpd_waic)")