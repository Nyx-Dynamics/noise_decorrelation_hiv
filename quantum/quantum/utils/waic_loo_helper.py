#!/usr/bin/env python3
"""
WAIC/PSIS-LOO helper for Bayesian v3.6 runs.

Usage (from repo root):

  # Auto-discover up to the N most recent traces under results_v3_6/runs/
  python quantum/quantum/utils/waic_loo_helper.py --auto 5

  # Or provide explicit comma-separated list of .nc files (absolute or relative paths)
  python quantum/quantum/utils/waic_loo_helper.py --traces path/to/trace1.nc,path/to/trace2.nc

The script prints a comparison table (ArviZ az.compare) for WAIC and LOO,
and saves a CSV `waic_loo_compare.csv` into the most recent run folder if found,
otherwise in the current working directory.
"""

import argparse
from pathlib import Path
import sys
import arviz as az
import pandas as pd
import xarray as xr


def find_recent_traces(base: Path, limit: int = 5):
    traces = []
    if not base.exists():
        return traces
    # Collect all .nc files under runs/* (one level down)
    for run_dir in sorted(base.glob('*/'), key=lambda p: p.stat().st_mtime, reverse=True):
        for nc in sorted(run_dir.glob('*.nc')):
            traces.append(nc)
        for nc in sorted((run_dir / 'original').glob('*.nc')):
            traces.append(nc)
        if len(traces) >= limit:
            break
    return traces[:limit]


def load_idatas(paths):
    models = {}
    for p in paths:
        try:
            idata = az.from_netcdf(str(p))
            key = p.stem
            # Ensure unique keys
            k = key
            k_i = 1
            while k in models:
                k = f"{key}_{k_i}"
                k_i += 1
            models[k] = idata
        except Exception as e:
            print(f"⚠ Failed to load {p}: {e}")
    return models


def combine_log_likelihood(idata: az.InferenceData) -> az.InferenceData:
    """Return a copy of idata with a single combined log_likelihood variable.

    Some ArviZ versions require `var_name` when multiple observed log-likelihood
    arrays are present. This routine sums all available log_likelihood data_vars
    into a single `log_likelihood` array so WAIC/LOO can be computed without
    additional arguments.
    """
    try:
        llk = idata.log_likelihood
    except Exception:
        return idata
    if llk is None or not hasattr(llk, 'data_vars'):
        return idata
    vars_list = list(llk.data_vars)
    if len(vars_list) <= 1 and ('log_likelihood' in vars_list or len(vars_list) == 1):
        # Already single-var or already named `log_likelihood` — nothing to do
        return idata
    # Sum all log-likelihood arrays pointwise across observed variables
    try:
        total = None
        for v in vars_list:
            arr = llk[v]
            # Sum over all dimensions except 'chain' and 'draw'
            # to get total log-likelihood for each sample
            dims_to_sum = [d for d in arr.dims if d not in ('chain', 'draw')]
            if dims_to_sum:
                arr = arr.sum(dim=dims_to_sum)
            
            if total is None:
                total = arr
            else:
                total = total + arr
        
        # New log_likelihood has dims (chain, draw)
        new_llk = xr.Dataset({'log_likelihood': total})
        # Create a shallow copy of idata and replace log_likelihood dataset
        groups = {g: getattr(idata, g) for g in idata.groups() if g != 'log_likelihood'}
        idata2 = az.InferenceData(**groups)
        idata2.add_groups({'log_likelihood': new_llk})
        return idata2
    except Exception as e:
        print(f"⚠ Failed to combine log_likelihood: {e}")
        # Fallback to original if combination fails
        return idata


def most_recent_run_dir(base_runs: Path):
    if not base_runs.exists():
        return None
    runs = [p for p in base_runs.glob('*/') if p.is_dir()]
    if not runs:
        return None
    runs.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return runs[0]


def main():
    parser = argparse.ArgumentParser(description='Compute WAIC/LOO for Bayesian v3.6 traces')
    parser.add_argument('--auto', type=int, default=0, help='Auto-discover up to N most recent traces')
    parser.add_argument('--traces', type=str, default='', help='Comma-separated list of .nc trace files')
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[3]
    runs_base = repo_root / 'quantum' / 'quantum' / 'results_v3_6' / 'runs'

    paths = []
    if args.traces.strip():
        for tok in args.traces.split(','):
            p = Path(tok.strip())
            if not p.is_absolute():
                p = (repo_root / p).resolve()
            if p.exists():
                paths.append(p)
            else:
                print(f"⚠ Trace not found: {p}")
    elif args.auto and args.auto > 0:
        paths = find_recent_traces(runs_base, limit=args.auto)
        if not paths:
            print('⚠ No recent traces found under results_v3_6/runs/')
    else:
        print('Provide --traces list or --auto N.')
        sys.exit(1)

    if not paths:
        print('❌ No valid traces to compare.')
        sys.exit(1)

    print('Comparing traces:')
    for p in paths:
        print('  -', p)

    models = load_idatas(paths)
    if len(models) < 2:
        print('⚠ Provide at least two traces for a meaningful comparison.')

    # Filter models that actually contain log_likelihood
    def _has_llk(idata):
        try:
            return hasattr(idata, 'log_likelihood') and idata.log_likelihood is not None and len(idata.log_likelihood.data_vars) > 0
        except Exception:
            return False

    valid_models = {name: idata for name, idata in models.items() if _has_llk(idata)}
    invalid = [name for name in models.keys() if name not in valid_models]
    if invalid:
        print('\n⚠ Skipping models without log_likelihood (re-run those fits to enable WAIC/LOO):')
        for name in invalid:
            print('  -', name)
        print('  Hint: re-run via Makefile targets (e.g., `make bayes-v36`) to regenerate traces with log_likelihood.')

    # Compute WAIC and LOO per model (use all for per-model table, but comparisons only on valid subset)
    rows = []
    for name, idata in models.items():
        # Ensure single combined log_likelihood for ArviZ versions that require it
        idata_for_ic = combine_log_likelihood(idata)
        try:
            waic = az.waic(idata_for_ic)
        except Exception as e:
            print(f"⚠ WAIC failed for {name}: {e}")
            waic = None
        try:
            loo = az.loo(idata_for_ic)
        except Exception as e:
            print(f"⚠ LOO failed for {name}: {e}")
            loo = None
        rows.append({
            'model': name,
            'waic': getattr(waic, 'waic', float('nan')) if waic is not None else float('nan'),
            'p_waic': getattr(waic, 'p_waic', float('nan')) if waic is not None else float('nan'),
            'elpd_waic': getattr(waic, 'elpd_waic', float('nan')) if waic is not None else float('nan'),
            'loo': getattr(loo, 'loo', float('nan')) if loo is not None else float('nan'),
            'p_loo': getattr(loo, 'p_loo', float('nan')) if loo is not None else float('nan'),
            'elpd_loo': getattr(loo, 'elpd_loo', float('nan')) if loo is not None else float('nan')
        })

    df = pd.DataFrame(rows).set_index('model')
    print('\nWAIC/LOO per model:')
    print(df)

    if len(valid_models) >= 2:
        try:
            comp_waic = az.compare({k: combine_log_likelihood(v) for k, v in valid_models.items()}, ic='waic')
            print('\nWAIC comparison (lower is better):')
            print(comp_waic)
        except Exception as e:
            print(f"⚠ az.compare WAIC failed: {e}")
            comp_waic = None
        try:
            comp_loo = az.compare({k: combine_log_likelihood(v) for k, v in valid_models.items()}, ic='loo')
            print('\nLOO comparison (lower is better):')
            print(comp_loo)
        except Exception as e:
            print(f"⚠ az.compare LOO failed: {e}")
            comp_loo = None
    else:
        comp_waic = comp_loo = None
        print('\nℹ Not enough valid models with log_likelihood for a pairwise comparison. Re-run at least two fits to regenerate traces with log_likelihood enabled.')

    # Save to the most recent run folder if available
    out_dir = most_recent_run_dir(runs_base) or Path.cwd()
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        out_df = df.copy()
        out_df.to_csv(out_dir / 'waic_loo_per_model.csv')
        if comp_waic is not None:
            comp_waic.to_csv(out_dir / 'waic_compare.csv')
        if comp_loo is not None:
            comp_loo.to_csv(out_dir / 'loo_compare.csv')
        print(f"\n✅ Saved WAIC/LOO tables under: {out_dir}")
    except Exception as e:
        print(f"⚠ Could not save WAIC/LOO CSVs: {e}")


if __name__ == '__main__':
    main()
