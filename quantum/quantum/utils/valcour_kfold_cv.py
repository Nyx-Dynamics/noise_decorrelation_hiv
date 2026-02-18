#!/usr/bin/env python3
"""
Valcour acute K-fold CV (quick, PSIS-based approximation)

What it does (quick version):
- Loads a single NetCDF trace (InferenceData) from a WITH-Valcour run
  (must contain individual acute observations with log_likelihood for
  'NAA_ratio_acute_obs').
- Uses PSIS-LOO pointwise log-likelihoods as a proxy to construct K folds
  over acute individuals without refitting K times.
- For each fold, computes ELPD on the held-out indices by subsetting the
  log_likelihood for 'NAA_ratio_acute_obs' and running ArviZ PSIS-LOO on
  that subset. Aggregates per-fold and overall metrics.

Notes:
- This is a fast approximation suitable for quick diagnostics. A full CV
  would re-fit the model on K-1 folds and score on the held-out fold.
- Requires the trace to be generated with idata_kwargs={"log_likelihood": True}
  and WITH Valcour individuals included so that 'NAA_ratio_acute_obs' exists.

Outputs:
- CSV summary saved under quantum/quantum/results_v3_6/runs/<cv_run>/
  - valcour_kfold_cv_summary.csv (per-fold and aggregate)
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import arviz as az
import xarray as xr
from datetime import datetime
import matplotlib.pyplot as plt


def find_trace(path_hint: str | None, runs_dir: Path, repo_root: Path) -> Path:
    if path_hint:
        p = Path(path_hint)
        if not p.is_absolute():
            # Resolve relative to repository root to avoid duplicating quantum/quantum
            p = (repo_root / p).resolve()
        if p.exists():
            return p
        raise FileNotFoundError(f"Trace not found: {p}")
    # Auto-pick latest trace
    candidates: list[Path] = []
    for run in sorted(runs_dir.glob('run_*'), key=lambda p: p.stat().st_mtime, reverse=True):
        for nc in run.glob('trace*.nc'):
            candidates.append(nc)
    if not candidates:
        raise FileNotFoundError("No traces found under results_v3_6/runs/")
    return candidates[0]


def subset_idata_llk_var(idata: az.InferenceData, var_name: str, held_idx: np.ndarray) -> az.InferenceData:
    """Create a new InferenceData where log_likelihood is the sum of the selected
    var_name over the held_idx only. This keeps other groups intact.
    """
    if not hasattr(idata, 'log_likelihood') or idata.log_likelihood is None:
        raise RuntimeError("InferenceData has no log_likelihood group.")
    llk = idata.log_likelihood
    if var_name not in llk:
        raise RuntimeError(f"log_likelihood does not contain '{var_name}'.")
    # Select points
    arr = llk[var_name]
    # arr dims typically: (chain, draw, obs) or (draw, obs)
    # Use xarray indexing on last dimension
    sel = arr[..., held_idx]
    # If multiple vars, we only keep the subset var
    new_llk = xr.Dataset({'log_likelihood': sel})
    # Assemble new idata
    out = az.InferenceData(**{g: getattr(idata, g) for g in idata.groups() if g != 'log_likelihood'})
    out.add_groups({'log_likelihood': new_llk})
    return out


def main():
    parser = argparse.ArgumentParser(description="Valcour acute K-fold CV (quick PSIS-based)")
    parser.add_argument('--trace', type=str, default='', help='Path to WITH-Valcour trace .nc (auto if omitted)')
    parser.add_argument('--kfold', type=int, default=5, help='Number of folds (default 5)')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for splits')
    parser.add_argument('--var', type=str, default='NAA_ratio_acute_obs', help='Acute LL var name')
    parser.add_argument('--index-map', type=str, default='',
                        help='Optional CSV mapping of acute indices to week labels. '
                             'Expected columns: idx (0-based), week_label or week_num')
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[3]
    runs_dir = repo_root / 'quantum' / 'quantum' / 'results_v3_6' / 'runs'

    trace_path = find_trace(args.trace.strip() or None, runs_dir, repo_root)
    print(f"Using trace: {trace_path}")
    idata = az.from_netcdf(str(trace_path))

    # Validate acute var
    if not hasattr(idata, 'log_likelihood') or idata.log_likelihood is None:
        raise RuntimeError("Trace has no log_likelihood; re-run with idata_kwargs={\"log_likelihood\": True}.")
    if args.var not in idata.log_likelihood:
        avail = list(idata.log_likelihood.data_vars)
        raise RuntimeError(f"'{args.var}' not in log_likelihood. Available: {avail}. Ensure this is a WITH-Valcour run.")

    # Determine number of acute observations
    ll_arr = idata.log_likelihood[args.var]
    n_obs = int(ll_arr.shape[-1])
    if n_obs < args.kfold:
        raise RuntimeError(f"Not enough acute observations for k={args.kfold} (n={n_obs}). Reduce k or include Valcour.")
    print(f"Acute observations detected: n={n_obs}")

    # Optional: load index->week mapping
    index_map_df = None
    if args.index_map.strip():
        map_path = Path(args.index_map)
        if not map_path.is_absolute():
            map_path = (repo_root / map_path).resolve()
        try:
            index_map_df = pd.read_csv(str(map_path))
            # Normalize column names
            cols = {c.lower(): c for c in index_map_df.columns}
            if 'idx' not in cols:
                raise ValueError("index-map CSV must contain 'idx' (0-based) column")
            week_col = None
            for cand in ['week_label', 'week', 'week_num', 'weeknumber']:
                if cand in cols:
                    week_col = cols[cand]
                    break
            if week_col is None:
                raise ValueError("index-map CSV must contain a week column: week_label/week/week_num")
            # Keep only needed
            index_map_df = index_map_df[[cols['idx'], week_col]].rename(columns={cols['idx']: 'idx', week_col: 'week'})
            print(f"Loaded index-map with {len(index_map_df)} rows from {map_path}")
        except Exception as e:
            print(f"⚠ Could not load index-map: {e}. Proceeding without per-week plots.")

    # Create K folds indices
    rng = np.random.default_rng(args.seed)
    indices = np.arange(n_obs)
    rng.shuffle(indices)
    folds = np.array_split(indices, args.kfold)

    results = []
    # Prepare CV run dir and figures dir early
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    cv_dir = runs_dir / f"valcour_cv_{ts}"
    fig_dir = cv_dir / 'figures'
    cv_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    for i, held in enumerate(folds, start=1):
        held = np.asarray(held)
        if held.size == 0:
            continue
        print(f"[CV] Fold {i}/{args.kfold}: held-out size = {held.size}")
        # Subset to held-out only (quick scoring)
        id_held = subset_idata_llk_var(idata, args.var, held)
        # Compute PSIS-LOO on held-out points
        try:
            loo = az.loo(id_held)
            elpd = float(loo.elpd_loo)
            p_loo = float(loo.p_loo)
            # Per-obs
            elpd_per_obs = elpd / held.size
            results.append({'fold': i, 'held_n': int(held.size), 'elpd_loo': elpd, 'p_loo': p_loo, 'elpd_loo_per_obs': elpd_per_obs})

            # Plot held-out pointwise elpd if available
            try:
                # In ArviZ >=0.16, loo has attribute loo_i; if not, skip
                loo_i = getattr(loo, 'loo_i', None)
                if loo_i is not None:
                    # Ensure 1D array
                    vals = np.ravel(np.array(loo_i))
                    plt.figure(figsize=(6,3))
                    plt.hist(vals, bins=20, color='C0', alpha=0.8)
                    plt.xlabel('held-out pointwise elpd (acute)')
                    plt.ylabel('count')
                    plt.title(f'Fold {i}: held-out elpd distribution (n={held.size})')
                    plt.tight_layout()
                    plt.savefig(fig_dir / f'cv_fold_{i}_heldout_elpd_hist.png', dpi=200)
                    plt.close()
            except Exception:
                pass

            # If week mapping is available, compute per-week held-out means and plot
            if index_map_df is not None:
                try:
                    # Build DataFrame for held indices
                    map_held = pd.DataFrame({'idx': held})
                    map_held = map_held.merge(index_map_df, on='idx', how='left')
                    # Group by week; if loo_i available, use it; else use per-obs avg elpd across held as proxy
                    if getattr(loo, 'loo_i', None) is not None:
                        vals = np.ravel(np.array(loo.loo_i))
                        map_held['elpd_i'] = vals
                    else:
                        # Uniform proxy if pointwise not available
                        map_held['elpd_i'] = elpd / held.size
                    per_week = map_held.groupby('week', dropna=False)['elpd_i'].mean().reset_index()
                    per_week.to_csv(cv_dir / f'cv_fold_{i}_per_week_elpd.csv', index=False)
                    # Plot bar
                    plt.figure(figsize=(6,3))
                    plt.bar(per_week['week'].astype(str), per_week['elpd_i'], color='C1', alpha=0.85)
                    plt.xlabel('week')
                    plt.ylabel('mean held-out elpd (acute)')
                    plt.title(f'Fold {i}: held-out mean elpd by week')
                    plt.tight_layout()
                    plt.savefig(fig_dir / f'cv_fold_{i}_per_week_elpd.png', dpi=200)
                    plt.close()
                except Exception as e:
                    print(f"⚠ Per-week plot failed for fold {i}: {e}")
        except Exception as e:
            results.append({'fold': i, 'held_n': int(held.size), 'elpd_loo': np.nan, 'p_loo': np.nan, 'elpd_loo_per_obs': np.nan, 'error': str(e)})

    df = pd.DataFrame(results)
    agg = {
        'held_n_total': int(df['held_n'].sum()),
        'folds': int(df.shape[0]),
        'mean_elpd_loo': float(df['elpd_loo'].mean()),
        'std_elpd_loo': float(df['elpd_loo'].std(ddof=1)) if df['elpd_loo'].count() > 1 else np.nan,
        'mean_elpd_loo_per_obs': float(df['elpd_loo_per_obs'].mean()),
        'std_elpd_loo_per_obs': float(df['elpd_loo_per_obs'].std(ddof=1)) if df['elpd_loo_per_obs'].count() > 1 else np.nan,
    }

    df.to_csv(cv_dir / 'valcour_kfold_cv_summary.csv', index=False)
    with open(cv_dir / 'valcour_kfold_cv_aggregate.txt', 'w') as f:
        for k, v in agg.items():
            f.write(f"{k},{v}\n")
    # Combined plot: per-fold elpd_loo_per_obs
    try:
        plt.figure(figsize=(6,3))
        plt.bar(df['fold'].astype(int).astype(str), df['elpd_loo_per_obs'], color='C2', alpha=0.85)
        plt.xlabel('fold')
        plt.ylabel('held-out mean elpd per obs')
        plt.title('Valcour acute K-fold CV (per-fold)')
        plt.tight_layout()
        plt.savefig(fig_dir / 'cv_per_fold_elpd_per_obs.png', dpi=200)
        plt.close()
    except Exception:
        pass

    # Combined per-week aggregation across folds (if per-week CSVs exist)
    try:
        if index_map_df is not None:
            per_week_files = sorted(cv_dir.glob('cv_fold_*_per_week_elpd.csv'))
            if per_week_files:
                frames = []
                for pf in per_week_files:
                    dfi = pd.read_csv(pf)
                    # tag fold from filename
                    try:
                        fold_str = pf.stem.split('_')[2]
                        fold_num = int(fold_str)
                    except Exception:
                        fold_num = np.nan
                    dfi['fold'] = fold_num
                    frames.append(dfi)
                allw = pd.concat(frames, ignore_index=True)
                # Sanitize week labels
                allw['week'] = allw['week'].astype(str)
                agg_week = allw.groupby('week', dropna=False)['elpd_i'].agg(['mean','std','count']).reset_index()
                agg_week.rename(columns={'mean':'mean_elpd','std':'std_elpd','count':'n_folds'}, inplace=True)
                agg_week.to_csv(cv_dir / 'cv_per_week_elpd_summary.csv', index=False)
                # Plot combined per-week with error bars
                plt.figure(figsize=(6.5,3.2))
                x = np.arange(agg_week.shape[0])
                plt.bar(x, agg_week['mean_elpd'], yerr=agg_week['std_elpd'], color='C3', alpha=0.85, capsize=3)
                plt.xticks(x, agg_week['week'])
                plt.xlabel('week')
                plt.ylabel('mean held-out elpd (acute)')
                plt.title('Valcour acute K-fold CV — per-week (mean ± SD across folds)')
                plt.tight_layout()
                plt.savefig(fig_dir / 'cv_per_week_elpd_summary.png', dpi=200)
                plt.close()
    except Exception as e:
        print(f"⚠ Combined per-week aggregation failed: {e}")
    print(f"\nSaved per-fold summary to: {cv_dir / 'valcour_kfold_cv_summary.csv'}")
    print(f"Aggregate: {agg}")


if __name__ == '__main__':
    main()
