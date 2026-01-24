"""
Run Enzyme v3 forward model (ODE) on ratio datasets with ART-era partitioning.

Implements Option A (xi → Pi_xi) as the sole modulation path for enzyme activity.
Retains enzyme_activity_fold only for validation in the annotated inputs.

Outputs per run (non-overwriting):
  results/enzyme_v3/<ratio>/<era>/<run_id>/
    - predictions.csv (per-Phase NAA/Cr and Cho/Cr)
    - inputs_annotated.csv (original + derived columns and validation deltas)
    - run_manifest.json (provenance, parameters, checksums)
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

from quantum.data_loaders import load_enzyme_inputs
from quantum.enzyme_kinetics import EnzymeKinetics
from quantum.final_calibrated_model_v3 import (
    compute_viral_damage,
    compute_membrane_turnover,
    concentration_to_ratio,
)
from quantum.utils.run_manifest import (
    make_run_id,
    base_environment,
    write_manifest,
    Manifest,
)


def parse_args():
    p = argparse.ArgumentParser(
        description='ODE forward model v3 on ratio datasets (Option A: xi→Pi_xi)'
    )
    p.add_argument('--ratio', required=True, choices=['3_1_1', '1_1_1', '1_2_1', '2_1'],
                   help='Dataset ratio under data/.../data_ratios_comparison/')
    p.add_argument('--era', default='both', choices=['pre_modern', 'post_modern', 'both'],
                   help='ART era filter; both keeps all (including unknown)')
    p.add_argument('--output-root', default='results/enzyme_v3', help='Root for outputs')
    p.add_argument('--run-notes', default='', help='Freeform notes recorded in manifest')
    p.add_argument('--beta-xi', type=float, default=1.89, help='Protection exponent β_ξ')
    p.add_argument('--xi-baseline-nm', type=float, default=0.66, help='Baseline ξ (nm) for Π_ξ computation')
    p.add_argument('--duration-days', type=float, default=60.0, help='ODE integration duration in days')
    return p.parse_args()


def aggregate_predictions(df: pd.DataFrame, duration_days: float) -> pd.DataFrame:
    """
    For each Phase, run the ODE with average Pi_xi and report steady-state
    NAA/Cr and Cho/Cr.
    """
    rows = []
    for phase, sub in df.groupby('Phase'):
        if 'Pi_xi' not in sub.columns:
            continue
        Pi_xi = float(np.nanmean(sub['Pi_xi']))
        # Coherence term is already embedded via ξ→Π_ξ under Option A; set η_coh=1.0
        enzymes = EnzymeKinetics(Pi_xi=Pi_xi, eta_coh=1.0, viral_damage_factor=compute_viral_damage_phase(phase))
        # Membrane turnover based on clinical phase mapping used in v3
        mem_turn = compute_membrane_turnover_phase(phase)
        NAA_molar, Cho_molar = enzymes.integrate(duration_days=duration_days, membrane_turnover=mem_turn)
        # Convert to MRS ratios using calibrated Cr
        NAA_Cr, Cho_Cr = concentration_to_ratio(NAA_molar, Cho_molar)
        rows.append({
            'Phase': phase,
            'Pi_xi_mean': Pi_xi,
            'NAA_mM': NAA_molar * 1e3,
            'Cho_mM': Cho_molar * 1e3,
            'NAA_Cr': NAA_Cr,
            'Cho_Cr': Cho_Cr,
        })
    return pd.DataFrame(rows)


# Map dataset Phase labels to v3 condition names
_PHASE_TO_CONDITION: Dict[str, str] = {
    'Control': 'healthy',
    'control': 'healthy',
    'Healthy': 'healthy',
    'Acute': 'acute_HIV',
    'acute': 'acute_HIV',
    'Chronic': 'chronic_HIV',
    'chronic': 'chronic_HIV',
}


def compute_viral_damage_phase(phase: str) -> float:
    cond = _PHASE_TO_CONDITION.get(phase, 'healthy')
    return compute_viral_damage(cond)


def compute_membrane_turnover_phase(phase: str) -> float:
    cond = _PHASE_TO_CONDITION.get(phase, 'healthy')
    return compute_membrane_turnover(cond)


def main():
    args = parse_args()
    df, input_path = load_enzyme_inputs(args.ratio, beta_xi=args.beta_xi, xi_baseline_nm=args.xi_baseline_nm)

    # Filter by era if requested (keep unknown by default in 'both')
    if args.era != 'both':
        df = df[df['art_era'] == args.era]

    # Prepare output directory
    run_id = make_run_id()
    out_dir = Path(args.output_root) / args.ratio / args.era / run_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Annotated inputs (retain enzyme_activity_fold for validation only)
    inputs_annot_path = out_dir / 'inputs_annotated.csv'
    df.to_csv(inputs_annot_path, index=False)

    # Aggregate ODE predictions per Phase
    preds = aggregate_predictions(df, duration_days=args.duration_days)
    pred_path = out_dir / 'predictions.csv'
    preds.to_csv(pred_path, index=False)

    # Compose manifest
    env = base_environment()
    manifest = Manifest(
        run_id=run_id,
        timestamp_utc=run_id.split('_')[0],
        model_variant='enzyme_v3_ode',
        ratio=args.ratio,
        era=args.era,
        input_files={'enzyme_inputs': str(input_path)},
        cli_args={
            'beta_xi': str(args.beta_xi),
            'xi_baseline_nm': str(args.xi_baseline_nm),
            'duration_days': str(args.duration_days),
            'run_notes': args.run_notes,
        },
        mechanism={'beta_xi': args.beta_xi, 'xi_baseline_nm': args.xi_baseline_nm},
        data_summary={
            'n_rows': int(len(df)),
            'n_unknown_era': int((df['art_era'] == 'unknown').sum()) if 'art_era' in df.columns else 0,
            'n_phases': int(df['Phase'].nunique()) if 'Phase' in df.columns else 0,
        },
        environment=env,
        outputs={
            'predictions_csv': str(pred_path),
            'inputs_annotated_csv': str(inputs_annot_path),
        },
        validation_notes={
            'primary_parameter': 'xi_estimate_nm → Pi_xi',
            'derived_parameter': 'Pi_xi computed via β_ξ (Option A)',
            'validation_column': 'enzyme_activity_fold retained for cross-check only',
            'era_handling': 'categorical covariate; does not modify xi→Pi_xi',
        },
    )
    write_manifest(out_dir / 'run_manifest.json', manifest)

    print("\n✨ Enzyme v3 ODE run complete")
    print(f"  Output: {out_dir}")
    if not preds.empty:
        print("  Predictions (head):")
        print(preds.head().to_string(index=False))


if __name__ == '__main__':
    main()
