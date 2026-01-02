#!/usr/bin/env python3
"""
Inspect WAIC/PSIS-LOO for one or two traces and optionally compare them.

Usage examples (from repo root):
  # Single trace (v3.6 or v4) — core NAA likelihood terms only
  python inspect_waic.py --trace quantum/quantum/results_v3_6/runs/with_valcour_nopseudo/trace_with_valcour_nopseudo.nc --mode core

  # You can also pass a run directory; the tool will pick the trace*.nc inside
  python inspect_waic.py --trace /Users/acdstudpro/Library/Mobile\ Documents/com~apple~CloudDocs/PycharmProjects\ -\ Noise\ Decorrelation\ copy/results_v3_6/runs/with_valcour_waic --mode core

  # Two traces (with vs. without Valcour) — prints per-trace metrics and deltas
  python inspect_waic.py \
    --trace quantum/quantum/results_v3_6/runs/with_valcour_nopseudo/trace_with_valcour_nopseudo.nc \
            quantum/quantum/results_v3_6/runs/no_valcour_nopseudo/trace_no_valcour_nopseudo.nc \
    --mode core

Notes:
- mode=core will sum only the three core NAA observed terms (acute/chronic/control),
  which is memory-light and robust. mode=all sums all observed terms (heavier).
- Traces must include log_likelihood or be generated with idata_kwargs={"log_likelihood": True}.
"""

import argparse
from pathlib import Path
import arviz as az
import xarray as xr
import pandas as pd
import numpy as np

CORE_VARS = ['NAA_ratio_acute_obs','NAA_ratio_chronic_obs','NAA_ratio_control_obs']

def _pick_core_vars(llk) -> list:
    """Return a list of log_likelihood data_vars to represent core NAA terms
    (acute, chronic, control). Prefer individual acute if present; otherwise
    fall back to acute group-mean variables (e.g., NAA_Young_BG_Acute_mean_obs,
    NAA_Sailasuta_BG_Acute_mean_obs). Also try group-mean fallbacks for
    chronic/control if individual terms are absent.
    """
    vars_available = set(llk.data_vars)

    chosen = []
    # Acute
    if 'NAA_ratio_acute_obs' in vars_available:
        chosen.append('NAA_ratio_acute_obs')
    else:
        # Fallback to any acute group-mean NAA obs present
        acute_gm = [v for v in vars_available if v.startswith('NAA_') and v.endswith('_Acute_mean_obs')]
        # Use all acute group means available (e.g., Young + Sailasuta)
        chosen.extend(sorted(acute_gm))

    # Chronic
    if 'NAA_ratio_chronic_obs' in vars_available:
        chosen.append('NAA_ratio_chronic_obs')
    else:
        chronic_gm = [v for v in vars_available if v.startswith('NAA_') and v.endswith('_Chronic_mean_obs')]
        # If multiple, include all
        if chronic_gm:
            chosen.extend(sorted(chronic_gm))

    # Control
    if 'NAA_ratio_control_obs' in vars_available:
        chosen.append('NAA_ratio_control_obs')
    else:
        control_gm = [v for v in vars_available if v.startswith('NAA_') and v.endswith('_Control_mean_obs')]
        if control_gm:
            chosen.extend(sorted(control_gm))

    # Ensure uniqueness & order preserved
    seen = set(); ordered = []
    for v in chosen:
        if v not in seen:
            ordered.append(v); seen.add(v)
    return ordered

def combine_llk(idata: az.InferenceData, mode: str = 'core') -> az.InferenceData:
    """Return a copy of idata with a single combined log_likelihood variable.
    mode='core' sums only the three NAA likelihoods; mode='all' sums all observed terms.
    """
    if not hasattr(idata, 'log_likelihood') or idata.log_likelihood is None:
        raise RuntimeError('No log_likelihood present; re-run sampling with idata_kwargs={"log_likelihood": True}.')
    llk = idata.log_likelihood
    if mode == 'core':
        # Build a harmonized three-term core set with sensible fallbacks
        use = _pick_core_vars(llk)
        if not use:
            # Fallback: sum everything if core cannot be determined
            use = list(llk.data_vars)
    else:
        use = list(llk.data_vars)
    total = None
    for v in use:
        arr = llk[v]
        total = arr if total is None else (total + arr)
    new_llk = xr.Dataset({'log_likelihood': total})
    out = az.InferenceData(**{g: getattr(idata, g) for g in idata.groups() if g != 'log_likelihood'})
    out.add_groups({'log_likelihood': new_llk})
    return out

def _summarize_trace(idata: az.InferenceData) -> dict:
    """Return basic sampler stats for additional context (divergences, chains, draws)."""
    out = {'divergences': None, 'n_chains': None, 'n_draws': None}
    try:
        ss = idata.sample_stats
        # divergences may be (chain, draw)
        if 'diverging' in ss:
            div = ss['diverging']
            # xarray: sum over all dims
            out['divergences'] = int(np.asarray(div).sum())
        # infer chains & draws
        for group in ['posterior', 'sample_stats']:
            if hasattr(idata, group) and getattr(idata, group) is not None:
                sel = getattr(idata, group)
                if 'chain' in sel.dims and 'draw' in sel.dims:
                    out['n_chains'] = int(sel.sizes['chain'])
                    out['n_draws'] = int(sel.sizes['draw'])
                    break
    except Exception:
        pass
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Inspect WAIC/LOO for one or two traces (with optional comparison)')
    ap.add_argument('--trace', nargs='+', required=True, help='Path(s) to .nc trace file(s). Provide one or two paths.')
    ap.add_argument('--mode', choices=['core','all'], default='core', help='Combine core NAA terms only or all observed terms')
    args = ap.parse_args()

    paths = [Path(p) for p in args.trace]
    if len(paths) > 2:
        raise SystemExit('Provide at most two traces for comparison.')

    def _resolve_to_trace(p: Path) -> Path:
        """Return a .nc trace path. If p is a directory, pick a suitable trace*.nc in it."""
        if p.is_file():
            return p
        if p.is_dir():
            cands = sorted(p.glob('trace*.nc'))
            if not cands:
                raise SystemExit(f"No trace .nc files found under directory: {p}")
            # Preference order: exact conventional names, then those without '_r', else first
            prefer = ['trace_with_valcour.nc', 'trace_no_valcour.nc']
            for name in prefer:
                for c in cands:
                    if c.name == name:
                        return c
            no_r = [c for c in cands if '_r.' not in c.name and not c.name.endswith('_r.nc')]
            return no_r[0] if no_r else cands[0]
        # Try to expand tilde etc.
        q = Path(str(p))
        if q.exists():
            return q
        raise SystemExit(f"Path not found: {p}")

    results = []
    chosen_terms = []
    for p in paths:
        trace_path = _resolve_to_trace(p)
        idata = az.from_netcdf(str(trace_path))
        # Summarize sampler
        sam_summary = _summarize_trace(idata)
        idata2 = combine_llk(idata, mode=args.mode)
        waic = az.waic(idata2)
        loo  = az.loo(idata2)
        # Count pointwise terms
        try:
            combined_ll = idata2.log_likelihood['log_likelihood']
            n_points = combined_ll.shape[-1]
        except Exception:
            n_points = 'unknown'
        # Store
        results.append({
            'name': trace_path.name,
            'path': str(trace_path),
            'elpd_waic': float(waic.elpd_waic),
            'p_waic': float(waic.p_waic),
            'waic': float(-2.0 * waic.elpd_waic),
            'elpd_loo': float(loo.elpd_loo),
            'p_loo': float(loo.p_loo),
            'loo': float(-2.0 * loo.elpd_loo),
            'n_points': n_points,
            **sam_summary,
        })
        try:
            if args.mode == 'core':
                original_llk = idata.log_likelihood
                chosen = _pick_core_vars(original_llk)
            else:
                chosen = list(idata.log_likelihood.data_vars)
        except Exception:
            chosen = []
        chosen_terms.append((trace_path.name, chosen))

    # Render
    print(f"\nWAIC/LOO (mode={args.mode})")
    df = pd.DataFrame(results)
    display_cols = ['name','waic','p_waic','elpd_waic','loo','p_loo','elpd_loo','divergences','n_chains','n_draws','n_points']
    print(df[display_cols].to_string(index=False))

    # If two traces: compute deltas (trace2 - trace1)
    if len(results) == 2:
        a, b = results[0], results[1]
        d_waic = b['waic'] - a['waic']
        d_loo  = b['loo']  - a['loo']
        print("\nPairwise deltas (trace2 - trace1):")
        print(f"  ΔWAIC = {d_waic:.3f}  (negative favors trace2)")
        print(f"  ΔLOO  = {d_loo:.3f}  (negative favors trace2)")
        # Special labeling if filenames hint at with/no_valcour
        names = (a['name'].lower(), b['name'].lower())
        if ('with_valcour' in names[0] and 'no_valcour' in names[1]) or ('with_valcour' in names[1] and 'no_valcour' in names[0]):
            if 'with_valcour' in names[0]:
                d_w = a['waic']
                d_n = b['waic']
                print(f"  ΔWAIC (no - with) = {d_n - d_w:.3f}")
                d_w = a['loo']
                d_n = b['loo']
                print(f"  ΔLOO  (no - with) = {d_n - d_w:.3f}")
            else:
                d_w = b['waic']
                d_n = a['waic']
                print(f"  ΔWAIC (no - with) = {d_n - d_w:.3f}")
                d_w = b['loo']
                d_n = a['loo']
                print(f"  ΔLOO  (no - with) = {d_n - d_w:.3f}")

    # Reveal which LL terms were included
    try:
        for name, chosen in chosen_terms:
            if chosen:
                print(f"\nIncluded LL terms for {name}: {chosen}")
    except Exception:
        pass