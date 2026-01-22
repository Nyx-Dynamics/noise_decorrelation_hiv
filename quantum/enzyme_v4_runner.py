#!/usr/bin/env python3
"""
CLI runner for Bayesian Enzyme v4.0 validation on specified cohort configuration.

Runs the mechanistic enzyme kinetics model (NAT8L/ASPA/Kennedy) with
canonical protection mapping:
    Pi_xi(ξ) = (ξ_ref / ξ) ** beta_xi,  beta_xi > 0

Outputs
- NetCDF trace: trace_enzyme_v4.nc
- Simple predictions CSV: predictions.csv  (Phase, NAA_Cr, Cho_Cr)
- Summary CSV: summary_v4.csv
- Optional figures (PDF) if matplotlib is available
- Results directory: results/enzyme_v4/<ratio>/<era>/<TIMESTAMP_RUNID>/

Example
python -m quantum.enzyme_v4_runner \
  --ratio 3_1_1 --era both --xi-ref 0.66 \
  --draws 3000 --tune 2000 --chains 4 --target-accept 0.95 --seed 42
"""
from __future__ import annotations

import argparse
import json
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import pymc as pm
import arviz as az

try:
    import matplotlib.pyplot as plt  # optional
except Exception:
    plt = None

from .bayesian_enzyme_v4 import build_enzyme_model, forward_model_enzyme, CONDITIONS


def build_outdir(ratio: str, era: str) -> Path:
    base = Path('results') / 'enzyme_v4' / ratio / era
    ts = datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
    runid = f"{ts}_{np.base_repr(int(time.time()*1e6), 36).lower()}"
    out = base / runid
    out.mkdir(parents=True, exist_ok=True)
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description='Run Bayesian Enzyme v4.0 validation')
    p.add_argument('--ratio', default='3_1_1', help='Cohort ratio key (e.g., 3_1_1)')
    p.add_argument('--era', default='both', choices=['both', 'pre_modern', 'post_modern'], help='ART era subset')
    p.add_argument('--xi-ref', type=float, default=0.66, help='Reference ξ in nm')
    p.add_argument('--draws', type=int, default=3000)
    p.add_argument('--tune', type=int, default=2000)
    p.add_argument('--chains', type=int, default=4)
    p.add_argument('--target-accept', type=float, default=0.95)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--outdir', default=None, help='Override output directory')
    return p.parse_args()


def main():
    args = parse_args()

    # Prepare output directory
    outdir = Path(args.outdir) if args.outdir else build_outdir(args.ratio, args.era)
    outdir.mkdir(parents=True, exist_ok=True)

    # Build model and sample
    print('Building enzyme v4.0 model...')
    with build_enzyme_model() as model:
        print('Sampling...')
        idata = pm.sample(
            draws=args.draws,
            tune=args.tune,
            chains=args.chains,
            random_seed=args.seed,
            target_accept=args.target_accept,
            return_inferencedata=True,
            progressbar=True,
        )

    # Save trace
    trace_path = outdir / 'trace_enzyme_v4.nc'
    az.to_netcdf(idata, trace_path)

    # Posterior means for prediction
    post = idata.posterior
    xi_acute = float(post['xi_acute'].mean()) if 'xi_acute' in post else 0.6e-9
    xi_chronic = float(post['xi_chronic'].mean()) if 'xi_chronic' in post else 0.8e-9
    beta_xi = float(post['beta_xi'].mean()) if 'beta_xi' in post else 2.0
    gamma_coh = float(post['gamma_coh'].mean()) if 'gamma_coh' in post else 1.5
    vdam_a = float(post['viral_damage_acute'].mean()) if 'viral_damage_acute' in post else 0.95
    vdam_c = float(post['viral_damage_chronic'].mean()) if 'viral_damage_chronic' in post else 0.90
    mem_a = float(post['membrane_acute'].mean()) if 'membrane_acute' in post else 2.0
    mem_c = float(post['membrane_chronic'].mean()) if 'membrane_chronic' in post else 1.2

    NAA_pred, Cho_pred = forward_model_enzyme(
        xi_acute=xi_acute,
        xi_chronic=xi_chronic,
        beta_xi=beta_xi,
        gamma_coh=gamma_coh,
        viral_damage_acute=vdam_a,
        viral_damage_chronic=vdam_c,
        membrane_acute=mem_a,
        membrane_chronic=mem_c,
    )

    # Save predictions CSV in manuscript schema
    pred_df = pd.DataFrame({
        'Phase': ['Control', 'Acute', 'Chronic'],
        'NAA_Cr': np.asarray(NAA_pred).astype(float),
        'Cho_Cr': np.asarray(Cho_pred).astype(float),
    })
    pred_df.to_csv(outdir / 'predictions.csv', index=False)

    # Save simple summary
    summ = az.summary(idata, var_names=['beta_xi','xi_acute','xi_chronic','gamma_coh','viral_damage_acute','viral_damage_chronic','membrane_acute','membrane_chronic'], kind='stats')
    summ.to_csv(outdir / 'summary_v4.csv')

    # Optional figures
    if plt is not None:
        try:
            # Posterior violin / KDE for beta_xi
            b = post['beta_xi'].values.ravel()
            plt.figure(figsize=(4,3))
            az.plot_kde(b); plt.title('β_ξ (enzyme v4)'); plt.tight_layout()
            plt.savefig(outdir / 'v4_beta_xi_kde.pdf'); plt.close()

            # Phase bar plots
            plt.figure(figsize=(6,3))
            plt.subplot(1,2,1); plt.bar(['Ctrl','Acute','Chronic'], pred_df['NAA_Cr']); plt.title('NAA/Cr')
            plt.subplot(1,2,2); plt.bar(['Ctrl','Acute','Chronic'], pred_df['Cho_Cr']); plt.title('Cho/Cr')
            plt.tight_layout(); plt.savefig(outdir / 'v4_phase_bars.pdf'); plt.close()
        except Exception:
            pass

    # Write a minimal run_info.json
    (outdir / 'run_info.json').write_text(json.dumps({
        'ratio': args.ratio,
        'era': args.era,
        'xi_ref': args.xi_ref,
        'draws': args.draws,
        'tune': args.tune,
        'chains': args.chains,
        'target_accept': args.target_accept,
        'seed': args.seed,
        'trace': str(trace_path),
        'predictions': str(outdir / 'predictions.csv'),
        'summary': str(outdir / 'summary_v4.csv'),
    }, indent=2))

    print('\n✅ Enzyme v4.0 run complete')
    print('  Output:', outdir)
    print('  Trace :', trace_path)
    print('  Preds :', outdir / 'predictions.csv')


if __name__ == '__main__':
    main()
