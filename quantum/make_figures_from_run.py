#!/usr/bin/env python3
"""
Generate manuscript figures for a specific v3.6 run, robustly.

This script:
- Accepts explicit paths (--trace, --pred) or a run directory (--run-dir)
- Builds a plotting-friendly NetCDF that:
  * adds legacy variable aliases expected by older plotting code
    (xi_acute_nm, xi_chronic_nm, xi_healthy_nm)
  * computes protection factors Pi_healthy, Pi_acute, Pi_chronic
    from beta_xi and xi per-phase using the canonical mapping
      Pi(ξ) = (ξ_ref / ξ) ** beta_xi, with ξ_ref = 0.66 by default
- Calls quantum.generate_figures.{generate_figure1,2,3,4}
  passing explicit paths so no simulation branch is triggered.
- Writes outputs to --outdir (default figures/figures).

Usage examples:

# Using a run directory (auto-detects trace and pred files)
python -m quantum.make_figures_from_run \
  --run-dir results/bayesian_v3_6/3_1_1/both/20251125T160846Z_ce3b0657 \
  --outdir figures/figures

# Using explicit files
python -m quantum.make_figures_from_run \
  --trace results/.../trace_v3_6.nc \
  --pred  results/.../posterior_predictive_v3_6.csv \
  --outdir figures/figures
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Optional

import arviz as az
import numpy as np
import pandas as pd

# Lazy import to avoid circulars at module import time

def _build_plot_ready_trace(trace_in: Path, plot_out: Path, xi_ref: float = 0.66) -> Path:
    """Create a NetCDF trace with legacy aliases and Pi_* if missing.

    - Adds aliases:
        xi_acute_nm   <- xi_nm_acute
        xi_chronic_nm <- xi_nm_chronic
        xi_healthy_nm <- xi_nm_control
    - Adds protection factors:
        Pi_healthy, Pi_acute, Pi_chronic
    """
    idata = az.from_netcdf(trace_in)
    post = idata.posterior

    # Add legacy aliases if available
    assigns = {}
    if 'xi_nm_acute' in post and 'xi_acute_nm' not in post:
        assigns['xi_acute_nm'] = post['xi_nm_acute']
    if 'xi_nm_chronic' in post and 'xi_chronic_nm' not in post:
        assigns['xi_chronic_nm'] = post['xi_nm_chronic']
    if 'xi_nm_control' in post and 'xi_healthy_nm' not in post:
        assigns['xi_healthy_nm'] = post['xi_nm_control']

    # Compute Pi_* if not present but ingredients are available
    beta = post.get('beta_xi')
    xi_h = post.get('xi_nm_control')
    xi_a = post.get('xi_nm_acute')
    xi_c = post.get('xi_nm_chronic')
    if beta is not None:
        if xi_h is not None and 'Pi_healthy' not in post:
            assigns['Pi_healthy'] = (xi_ref / xi_h) ** beta
        if xi_a is not None and 'Pi_acute' not in post:
            assigns['Pi_acute'] = (xi_ref / xi_a) ** beta
        if xi_c is not None and 'Pi_chronic' not in post:
            assigns['Pi_chronic'] = (xi_ref / xi_c) ** beta
        # Provide a placeholder astrocyte compensation parameter for plotting if absent
        if 'astrocyte_comp' not in post:
            assigns['astrocyte_comp'] = beta*0 + 1.0  # unit baseline, uncorrelated with beta_xi

    if assigns:
        idata2 = idata.copy()
        idata2.posterior = post.assign(**assigns)
    else:
        idata2 = idata

    # Write to desired output file
    plot_out.parent.mkdir(parents=True, exist_ok=True)
    # If plot_out already exists, give it a new unique name
    out_path = plot_out
    if out_path.exists():
        stem = out_path.stem
        suffix = out_path.suffix
        for k in range(1, 100):
            candidate = out_path.with_name(f"{stem}_{k}{suffix}")
            if not candidate.exists():
                out_path = candidate
                break
    az.to_netcdf(idata2, out_path)
    return out_path


def _prepare_predictions_csv(pred_in: Path, out_csv: Path) -> Path:
    """Create a predictions CSV with the column names expected by generate_figure3.

    Input columns from runner:
      - 'Phase' in {'Control','Acute','Chronic'}
      - 'NAA_mean_observed'
      - 'NAA_mean_predicted'
    Output columns expected by plotting code:
      - 'condition' in {'healthy','acute_HIV','chronic_HIV'}
      - 'NAA_obs'
      - 'NAA_pred'
    """
    df = pd.read_csv(pred_in)
    # Map condition labels
    phase_to_cond = {
        'Control': 'healthy',
        'Acute': 'acute_HIV',
        'Chronic': 'chronic_HIV',
    }
    if 'condition' not in df.columns:
        if 'Phase' in df.columns:
            df['condition'] = df['Phase'].map(phase_to_cond).fillna('unknown')
        else:
            # Fallback: assume order healthy/acute/chronic if three rows
            default = ['healthy', 'acute_HIV', 'chronic_HIV']
            df['condition'] = (default + ['unknown'] * len(df))[:len(df)]
    # Create expected NAA columns
    if 'NAA_obs' not in df.columns and 'NAA_mean_observed' in df.columns:
        df['NAA_obs'] = df['NAA_mean_observed']
    if 'NAA_pred' not in df.columns and 'NAA_mean_predicted' in df.columns:
        df['NAA_pred'] = df['NAA_mean_predicted']
    # Keep only required columns
    keep = ['condition', 'NAA_obs', 'NAA_pred']
    missing = [c for c in keep if c not in df.columns]
    if missing:
        raise ValueError(f"Predictions CSV missing required columns after mapping: {missing}")

    # Aggregate to the three canonical conditions expected by the plotting code
    # If multiple rows per condition exist, take the mean values per condition.
    order = ['healthy', 'acute_HIV', 'chronic_HIV']
    agg = df[keep].groupby('condition', as_index=False).mean(numeric_only=True)
    agg = agg[agg['condition'].isin(order)].copy()
    agg['condition'] = pd.Categorical(agg['condition'], order, ordered=True)
    agg = agg.sort_values('condition')

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    agg.to_csv(out_csv, index=False)
    return out_csv


def _resolve_paths(args) -> tuple[Path, Path]:
    # Resolve trace
    if args.trace:
        trace = Path(args.trace).expanduser().resolve()
    elif args.run_dir:
        rd = Path(args.run_dir).expanduser().resolve()
        # Prefer already-prepared plot/alias files if present
        for name in ("trace_v3_6_plot.nc", "trace_v3_6_alias.nc", "trace_v3_6.nc"):
            cand = rd / name
            if cand.exists():
                trace = cand
                break
        else:
            raise FileNotFoundError(f"No trace file found in {rd} (looked for trace_v3_6*.nc)")
    else:
        raise ValueError("Provide either --trace or --run-dir")

    # Resolve predictions CSV
    if args.pred:
        pred = Path(args.pred).expanduser().resolve()
    elif args.run_dir:
        rd = Path(args.run_dir).expanduser().resolve()
        for name in ("posterior_predictive_v3_6.csv", "posterior_predictive.csv"):
            cand = rd / name
            if cand.exists():
                pred = cand
                break
        else:
            raise FileNotFoundError(f"No posterior predictive CSV found in {rd} (looked for posterior_predictive*_v3_6.csv)")
    else:
        raise ValueError("Provide either --pred or --run-dir")

    return trace, pred


def main():
    p = argparse.ArgumentParser(description="Generate manuscript figures for a specific v3.6 run")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--run-dir", help="Path to a run directory containing trace_v3_6.nc and posterior_predictive.csv")
    g.add_argument("--trace", help="Path to trace_v3_6.nc (or already-aliased) file")
    p.add_argument("--pred", help="Path to posterior_predictive_v3_6.csv (or posterior_predictive.csv)")
    p.add_argument("--xi-ref", type=float, default=0.66, help="Reference xi (nm) for Pi computation; default 0.66")
    p.add_argument("--outdir", default="figures/figures", help="Output directory for figures")
    p.add_argument("--force-rebuild", action="store_true", help="Force rebuilding the plot-ready NetCDF with aliases and Pi_*")
    args = p.parse_args()

    trace_in, pred_in = _resolve_paths(args)

    # Prepare a plot-ready trace with aliases + Pi_* if needed
    if args.run_dir:
        rd_abs = Path(args.run_dir).expanduser().resolve()
        plot_nc = rd_abs / "trace_v3_6_plot.nc"
    else:
        # If only a trace was provided, create the plot file next to it
        plot_nc = Path(trace_in).with_name("trace_v3_6_plot.nc").resolve()

    # Only (re)build if we don't already have a valid plot file
    need_build = True
    if plot_nc.exists() and not args.force_rebuild:
        try:
            _ = az.from_netcdf(plot_nc)
            need_build = False
            print(f"Using existing plot-ready trace: {plot_nc}")
        except Exception:
            need_build = True
    if need_build:
        print(f"Building plot-ready trace at: {plot_nc}")
        plot_nc = _build_plot_ready_trace(trace_in, plot_nc, xi_ref=float(args.xi_ref))
        print(f"✓ Wrote plot-ready trace: {plot_nc}")

    # Ensure output dir exists and chdir so helper functions save there
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    # Build a plotting-friendly predictions CSV with expected columns
    pred_plot = _prepare_predictions_csv(Path(pred_in), outdir / "posterior_predictive_for_fig3.csv")

    cwd = Path.cwd()
    try:
        os.chdir(outdir)
        # Import plotting helpers only now (ensures project package import works from any cwd)
        from quantum.generate_figures import (
            generate_figure1, generate_figure2, generate_figure3, generate_figure4
        )
        print("Generating Figure 1...")
        generate_figure1()
        print("Generating Figure 2...")
        print(f"Using trace file: {plot_nc}")
        generate_figure2(str(plot_nc))
        print("Generating Figure 3...")
        print(f"Using predictions file: {pred_plot}")
        generate_figure3(str(plot_nc), str(pred_plot))
        print("Generating Figure 4...")
        generate_figure4(str(plot_nc))
        print(f"\n✓ All figures written to {outdir}")
    finally:
        os.chdir(cwd)


if __name__ == "__main__":
    main()
