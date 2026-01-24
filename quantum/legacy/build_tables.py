#!/usr/bin/env python3
"""
Build concise ablation and ratio‑configuration tables for the manuscript.

Features
- Ablation table (Main): compares Full vs Constant‑xi vs No‑protection (beta_xi=0)
  * Computes WAIC/ELPD via ArviZ when possible from trace NetCDFs; otherwise reads a fallback metrics CSV.
  * Optionally computes a simple PPC error (%) by aggregating posterior predictive CSVs
    to the three canonical conditions [healthy, acute_HIV, chronic_HIV].
- Ratio‑configuration table (Main): collates key stats per configuration/run
  * beta_xi mean [94% HDI], P(xi_acute < xi_chronic), and optional PPC error (%)

Inputs
- Trace NetCDFs (e.g., trace_v3_6.nc or *_plot.nc) for each model/config, when available
- Posterior predictive CSVs (posterior_predictive*.csv) for PPC error (optional)
- Optional metrics CSV fallback with columns:
    model, config, waic, elpd_diff, elpd_se

Outputs
- Markdown tables saved under --outdir (default: tables/)
    ablation_table.md, ratio_config_table.md

Examples
1) Ablation from traces and predictions (preferred):
   python -m quantum.build_tables \
     --ablation \
     --full-trace results/bayesian_v3_6/3_1_1/both/FINAL/trace_v3_6.nc \
     --constxi-trace results/.../const_xi/trace_v3_6.nc \
     --noprot-trace results/.../beta0/trace_v3_6.nc \
     --full-pred results/bayesian_v3_6/3_1_1/both/FINAL/posterior_predictive.csv \
     --constxi-pred results/.../const_xi/posterior_predictive.csv \
     --noprot-pred results/.../beta0/posterior_predictive.csv \
     --outdir tables

2) Ratio config summary (multiple runs):
   python -m quantum.build_tables \
     --ratio-config \
     --entry 3_1_1:both:results/bayesian_v3_6/3_1_1/both/FINAL/trace_v3_6.nc:results/.../posterior_predictive.csv \
     --entry 3_1_1:pre_modern:results/.../trace_v3_6.nc:results/.../posterior_predictive.csv \
     --entry 3_1_1:post_modern:results/.../trace_v3_6.nc:results/.../posterior_predictive.csv \
     --outdir tables
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Tuple, Dict, List

import arviz as az
import numpy as np
import pandas as pd

# -----------------------------
# Utilities
# -----------------------------

def hdi94(x: np.ndarray) -> Tuple[float, float]:
    try:
        h = az.hdi(x, hdi_prob=0.94)
        return float(h[0]), float(h[1])
    except Exception:
        lo, hi = np.quantile(x, [0.03, 0.97])
        return float(lo), float(hi)


def ppc_error_pct(pred_csv: Path) -> Optional[float]:
    """Compute a simple 3-condition PPC % error (mean absolute percent error across conditions).
    Accepts either the raw runner CSV (Phase/NAA_mean_observed/NAA_mean_predicted)
    or the plotting-friendly CSV (condition/NAA_obs/NAA_pred). Aggregates to 3 conditions.
    """
    if not pred_csv or not Path(pred_csv).exists():
        return None
    df = pd.read_csv(pred_csv)
    if {'condition', 'NAA_obs', 'NAA_pred'}.issubset(df.columns):
        data = df[['condition', 'NAA_obs', 'NAA_pred']].copy()
    else:
        # Map runner columns
        if not {'Phase', 'NAA_mean_observed', 'NAA_mean_predicted'}.issubset(df.columns):
            return None
        phase_to_cond = {'Control': 'healthy', 'Acute': 'acute_HIV', 'Chronic': 'chronic_HIV'}
        data = pd.DataFrame({
            'condition': df['Phase'].map(phase_to_cond).fillna('unknown'),
            'NAA_obs': df['NAA_mean_observed'],
            'NAA_pred': df['NAA_mean_predicted'],
        })
    # Aggregate to 3 conditions
    order = ['healthy', 'acute_HIV', 'chronic_HIV']
    agg = data.groupby('condition', as_index=False).mean(numeric_only=True)
    agg = agg[agg['condition'].isin(order)].copy()
    if agg.empty:
        return None
    err = (agg['NAA_pred'] - agg['NAA_obs']).abs() / agg['NAA_obs'].replace(0, np.nan)
    return float(err.mean() * 100.0)


# -----------------------------
# Ablation table
# -----------------------------

def waic_from_idata(trace_nc: Optional[Path]) -> Optional[Tuple[float, float]]:
    if trace_nc is None or not Path(trace_nc).exists():
        return None
    try:
        idata = az.from_netcdf(trace_nc)
        w = az.waic(idata)
        # Return WAIC and ELPD_diff (relative to itself = 0)
        return float(w.waic), 0.0
    except Exception:
        return None


def extract_ablation(ppc_paths: Dict[str, Optional[Path]], trace_paths: Dict[str, Optional[Path]], metrics_fallback: Optional[Path]) -> pd.DataFrame:
    rows = []
    # Attempt WAIC from traces
    waics: Dict[str, Optional[float]] = {}
    for key, p in trace_paths.items():
        val = waic_from_idata(p)
        waics[key] = val[0] if val is not None else None
    # If any WAIC missing and fallback provided, read it
    if any(v is None for v in waics.values()) and metrics_fallback and Path(metrics_fallback).exists():
        m = pd.read_csv(metrics_fallback)
        # Expect rows for model variants: full, constxi, noprot
        for key in ('full', 'constxi', 'noprot'):
            if waics.get(key) is None:
                sub = m[m['variant'] == key]
                if not sub.empty and 'waic' in sub.columns:
                    waics[key] = float(sub['waic'].iloc[0])
    # Compute ΔWAIC vs full when possible
    d_waic = {k: None for k in ('full', 'constxi', 'noprot')}
    if waics.get('full') is not None:
        full = waics['full']
        for k, v in waics.items():
            if v is not None:
                d_waic[k] = float(v - full)
    # PPC errors
    ppc_err = {k: ppc_error_pct(ppc_paths.get(k)) if ppc_paths.get(k) else None for k in ('full', 'constxi', 'noprot')}
    # Build table
    def row(name: str, label: str) -> Dict[str, Optional[float]]:
        return {
            'Variant': label,
            'ΔWAIC vs Full': d_waic[name],
            'PPC error (%)': ppc_err[name],
        }
    rows.append(row('full', 'Full (phase‑ξ; β_ξ free)'))
    rows.append(row('constxi', 'Constant‑ξ (β_ξ free)'))
    rows.append(row('noprot', 'No protection (β_ξ = 0)'))
    return pd.DataFrame(rows)


# -----------------------------
# Ratio configuration table
# -----------------------------

def extract_ratio_entry(label: str, trace_nc: Path, pred_csv: Optional[Path]) -> Dict[str, object]:
    idata = az.from_netcdf(trace_nc)
    post = idata.posterior
    beta = np.asarray(post['beta_xi']).ravel()
    # Variable names per our final runner
    xi_a = np.asarray(post.get('xi_nm_acute'))
    xi_c = np.asarray(post.get('xi_nm_chronic'))
    if xi_a is None or xi_c is None:
        # Legacy aliases
        xi_a = np.asarray(post.get('xi_acute_nm'))
        xi_c = np.asarray(post.get('xi_chronic_nm'))
    xi_a = xi_a.values.ravel() if hasattr(xi_a, 'values') else xi_a.ravel()
    xi_c = xi_c.values.ravel() if hasattr(xi_c, 'values') else xi_c.ravel()

    lo, hi = hdi94(beta)
    p = float((xi_a < xi_c).mean())
    ppc = ppc_error_pct(pred_csv) if pred_csv else None
    return {
        'Configuration': label,
        'β_ξ mean [94% HDI]': f"{np.mean(beta):.2f} [{lo:.2f}, {hi:.2f}]",
        'P(ξ_acute < ξ_chronic)': f"{p:.2f}",
        'PPC error (%)': None if ppc is None else f"{ppc:.1f}",
    }


# -----------------------------
# CLI
# -----------------------------

def parse_args():
    p = argparse.ArgumentParser(description='Build concise ablation and ratio‑config tables')
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--ablation', action='store_true', help='Build ablation table')
    g.add_argument('--ratio-config', action='store_true', help='Build ratio‑config table')

    # Ablation inputs
    p.add_argument('--full-trace')
    p.add_argument('--constxi-trace')
    p.add_argument('--noprot-trace')
    p.add_argument('--full-pred')
    p.add_argument('--constxi-pred')
    p.add_argument('--noprot-pred')
    p.add_argument('--metrics-fallback', help='CSV with columns: variant, waic, elpd_diff, elpd_se')

    # Ratio config inputs: label:trace[:pred]
    p.add_argument('--entry', action='append', default=[], help='Format label:trace_path[:pred_path], e.g. 3_1_1:both:path/to/trace.nc:path/to/pred.csv')

    p.add_argument('--outdir', default='tables', help='Output directory for markdown tables')
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.ablation:
        trace_paths = {
            'full': Path(args.full_trace) if args.full_trace else None,
            'constxi': Path(args.constxi_trace) if args.constxi_trace else None,
            'noprot': Path(args.noprot_trace) if args.noprot_trace else None,
        }
        ppc_paths = {
            'full': Path(args.full_pred) if args.full_pred else None,
            'constxi': Path(args.constxi_pred) if args.constxi_pred else None,
            'noprot': Path(args.noprot_pred) if args.noprot_pred else None,
        }
        df = extract_ablation(ppc_paths, trace_paths, Path(args.metrics_fallback) if args.metrics_fallback else None)
        md = df.to_markdown(index=False)
        out = outdir / 'ablation_table.md'
        out.write_text(md)
        print(f"✓ Wrote {out}")

    if args.ratio_config:
        rows = []
        for e in args.entry:
            parts = e.split(':')
            if len(parts) < 3:
                print(f"WARN: --entry '{e}' malformed; expected label:era:trace[:pred]")
                continue
            label = f"{parts[0]} ({parts[1]})"
            trace = Path(parts[2])
            pred = Path(parts[3]) if len(parts) > 3 else None
            try:
                rows.append(extract_ratio_entry(label, trace, pred))
            except Exception as ex:
                print(f"WARN: failed to process entry '{e}': {ex}")
        if rows:
            df = pd.DataFrame(rows)
            md = df.to_markdown(index=False)
            out = outdir / 'ratio_config_table.md'
            out.write_text(md)
            print(f"✓ Wrote {out}")


if __name__ == '__main__':
    main()
