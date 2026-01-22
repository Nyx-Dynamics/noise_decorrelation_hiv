"""
Bayesian v3.6 runner with ART-era hierarchical effect while preserving
core mechanism Option A (xi → Pi_xi).

- Loads bayesian_inputs_<ratio>.csv via data_loaders.load_bayesian_inputs
- Treats ART era as a categorical covariate (additive effect), not altering xi→Pi_xi
- Uses a simple forward mapping for NAA/Cr: baseline * Pi_xi_phase with study and era effects
- Retains validation-only fields in the annotated inputs written to outputs
- Writes per-run manifest with git/env details and checksums

Assumptions about input CSV (flexible but expected columns):
  - 'Phase' with values among {Acute, Chronic, Control} (case-insensitive tolerated)
  - 'NAA_mean' and 'NAA_SE' (mean and standard error of NAA/Cr)
  - Optional 'study' for hierarchical random effects; if missing, a single study is assumed
  - Optional 'publication_year' or mapping inferred by loader for ART era

Outputs per run (non-overwriting):
  results/bayesian_v3_6/<ratio>/<era>/<run_id>/
    - trace_v3_6.nc (InferenceData NetCDF)
    - posterior_predictive.csv
    - inputs_annotated.csv
    - run_manifest.json
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import arviz as az
import pandas as pd
import pymc as pm
import numpy as np

from quantum.data_loaders import load_bayesian_inputs
from quantum.utils.run_manifest import make_run_id, base_environment, write_manifest, Manifest

# Phase label normalization
_PHASE_NORM: Dict[str, str] = {
    'control': 'Control', 'healthy': 'Control', 'hc': 'Control', 'null': 'Control',
    'acute': 'Acute', 'acute_hiv': 'Acute',
    'chronic': 'Chronic', 'chronic_hiv': 'Chronic'
}


def norm_phase(val: str) -> str:
    if pd.isna(val):
        return 'Control'
    key = str(val).strip().lower()
    return _PHASE_NORM.get(key, val)


def parse_args():
    p = argparse.ArgumentParser(description='Bayesian v3.6 with ART-era hierarchical effect (Option A mechanism)')
    p.add_argument('--ratio', required=True, choices=['3_1_1', '1_1_1', '1_2_1', '2_1'])
    p.add_argument('--era', default='both', choices=['pre_modern', 'post_modern', 'both'])
    p.add_argument('--draws', type=int, default=1500)
    p.add_argument('--tune', type=int, default=1000)
    p.add_argument('--chains', type=int, default=4)
    p.add_argument('--target-accept', type=float, default=0.92)
    p.add_argument('--seed', type=int, default=2025)
    p.add_argument('--output-root', default='results/bayesian_v3_6')
    p.add_argument('--run-notes', default='')
    # Mechanism parameters
    p.add_argument('--beta-xi-mean', type=float, default=1.89)
    p.add_argument('--beta-xi-sd', type=float, default=0.25)
    p.add_argument('--xi-baseline-nm', type=float, default=0.66)
    p.add_argument('--baseline-naa-cr', type=float, default=2.2, help='Healthy baseline NAA/Cr for scaling')
    # Legacy mirroring
    p.add_argument('--legacy-output-root', default=None, help='Optional path to mirror this run (copy) into an external folder (e.g., iCloud runs/run_YYYYMMDD_HHMMSS)')
    p.add_argument('--legacy-mirror-mode', default='copy', choices=['copy', 'symlink'], help='Mirror mode for legacy-output-root: copy files or create symlinks')
    return p.parse_args()


def prepare(df: pd.DataFrame, era: str) -> pd.DataFrame:
    df = df.copy()
    # Filter era if requested (keep unknown when era == 'both')
    if era != 'both' and 'art_era' in df.columns:
        df = df[df['art_era'] == era]

    # Normalize phases
    if 'Phase' in df.columns:
        df['Phase'] = df['Phase'].apply(norm_phase)
    else:
        raise ValueError("bayesian_inputs CSV must include 'Phase' column")

    # Auto-detect schema: filter to NAA/Cr metabolite if present and map Mean/SE → NAA_mean/NAA_SE
    if 'Metabolite' in df.columns:
        meta_norm = df['Metabolite'].astype(str).str.replace(' ', '', regex=False).str.upper()
        naa_mask = meta_norm.eq('NAA/CR') | meta_norm.eq('NAA/CRR') | meta_norm.eq('NAA/CRATIO')
        # Default to rows that look like NAA/Cr; if none match, keep all rows
        if naa_mask.any():
            df = df[naa_mask].copy()
    # If expected columns missing but generic Mean/SE exist, rename for NAA/Cr rows
    if ('NAA_mean' not in df.columns or 'NAA_SE' not in df.columns) and {'Mean', 'SE'}.issubset(df.columns):
        df = df.rename(columns={'Mean': 'NAA_mean', 'SE': 'NAA_SE'})

    # Require NAA mean/SE; drop rows where missing (documented in manifest size)
    needed = ['NAA_mean', 'NAA_SE']
    missing_cols = [c for c in needed if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in bayesian inputs: {missing_cols}")
    df = df[~df['NAA_mean'].isna() & ~df['NAA_SE'].isna()].copy()

    # Study id for random effects (optional)
    if 'study' not in df.columns:
        # try to derive a study name from Study/SourceFile if present
        if 'Study' in df.columns:
            df['study'] = df['Study']
        elif 'SourceFile' in df.columns:
            df['study'] = df['SourceFile'].astype(str).str.extract(r'([^/\\]+)')[0].fillna('pooled')
        else:
            df['study'] = 'pooled'
    df['study_id'] = df['study'].astype('category').cat.codes

    # Era index mapping for hierarchical covariate
    if 'art_era_idx' not in df.columns:
        # Backfill via art_era if present; else unknown = -1
        if 'art_era' in df.columns:
            df['art_era_idx'] = df['art_era'].map({'pre_modern': 0, 'post_modern': 1}).fillna(-1).astype(int)
        else:
            df['art_era_idx'] = -1

    return df


def build_model(df: pd.DataFrame, beta_xi_mean: float, beta_xi_sd: float, xi_baseline_nm: float, baseline_naacr: float):
    # Indexing helpers
    phases = ['Control', 'Acute', 'Chronic']
    phase_idx = df['Phase'].astype('category').cat.set_categories(phases).cat.codes.values
    study_ids = df['study_id'].values
    era_idx = df['art_era_idx'].values

    N = len(df)
    n_studies = int(df['study_id'].nunique())

    with pm.Model() as model:
        # Mechanism: beta_xi prior
        beta_xi = pm.TruncatedNormal('beta_xi', mu=beta_xi_mean, sigma=beta_xi_sd, lower=0.0)

        # Latent xi per phase (log-parameterization for improved geometry)
        log_xi_ctrl = pm.Normal('log_xi_control', mu=np.log(0.66), sigma=0.15)
        log_xi_acute = pm.Normal('log_xi_acute', mu=np.log(0.61), sigma=0.15)
        log_xi_chron = pm.Normal('log_xi_chronic', mu=np.log(0.80), sigma=0.15)
        xi_nm_ctrl = pm.Deterministic('xi_nm_control', pm.math.exp(log_xi_ctrl))
        xi_nm_acute = pm.Deterministic('xi_nm_acute', pm.math.exp(log_xi_acute))
        xi_nm_chron = pm.Deterministic('xi_nm_chronic', pm.math.exp(log_xi_chron))
        log_xi = pm.math.stack([log_xi_ctrl, log_xi_acute, log_xi_chron])

        # Protection factor per phase (canonical mapping: (xi_ref/xi)^{beta_xi} with beta_xi>0)
        Pi_xi_phase = pm.math.exp((-beta_xi) * (log_xi - np.log(xi_baseline_nm)))

        # Era additive effect (pre, post); unknown era gets 0 via masking
        era_effect = pm.Normal('era_effect', mu=0.0, sigma=0.10, shape=2)

        # Study random effects
        study_sd = pm.HalfNormal('study_sd', sigma=0.10)
        study_offset = pm.Normal('study_offset', mu=0.0, sigma=1.0, shape=n_studies)
        study_effect = pm.Deterministic('study_effect', study_offset * study_sd)

        # Observation noise (per-observation scaling around reported SE)
        se_scale = pm.HalfNormal('se_scale', sigma=0.20)

        # Expected mean per observation
        mu_base = baseline_naacr * Pi_xi_phase[phase_idx]
        # Era term (mask unknown=-1 to 0)
        era_term = pm.math.switch(pm.math.eq(era_idx, -1), 0.0, era_effect[era_idx])
        # Additive effects applied multiplicatively to baseline via (1 + effects)
        mu = mu_base * (1.0 + era_term + study_effect[study_ids])

        # Likelihood using reported SE (scaled)
        sigma = pm.math.maximum(1e-6, se_scale * df['NAA_SE'].values)
        pm.Normal('naa_obs', mu=mu, sigma=sigma, observed=df['NAA_mean'].values)

    return model


def main():
    args = parse_args()

    # Load and prepare data
    df, input_path = load_bayesian_inputs(args.ratio)
    df = prepare(df, args.era)

    # Output dir
    run_id = make_run_id()
    out_dir = Path(args.output_root) / args.ratio / args.era / run_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save annotated inputs for transparency
    inputs_annot_path = out_dir / 'inputs_annotated.csv'
    df.to_csv(inputs_annot_path, index=False)

    # Build and sample
    model = build_model(
        df=df,
        beta_xi_mean=args.beta_xi_mean,
        beta_xi_sd=args.beta_xi_sd,
        xi_baseline_nm=args.xi_baseline_nm,
        baseline_naacr=args.baseline_naa_cr,
    )

    with model:
        idata = pm.sample(
            draws=args.draws,
            tune=args.tune,
            chains=args.chains,
            random_seed=args.seed,
            target_accept=args.target_accept,
            return_inferencedata=True,
            progressbar=True,
        )
        ppc = pm.sample_posterior_predictive(idata, return_inferencedata=True)
        idata.extend(ppc)

    # Save trace
    trace_path = out_dir / 'trace_v3_6.nc'
    az.to_netcdf(idata, trace_path)

    # Posterior predictive summary for observed rows
    pp_means = idata.posterior_predictive['naa_obs'].stack(draws=("chain", "draw")).mean(dim='draws').values
    pp_df = pd.DataFrame({
        'Phase': df['Phase'].values,
        'study': df['study'].values,
        'art_era': df.get('art_era', pd.Series(['unknown'] * len(df))).values,
        'NAA_mean_observed': df['NAA_mean'].values,
        'NAA_SE_observed': df['NAA_SE'].values,
        'NAA_mean_predicted': pp_means,
    })
    pp_path = out_dir / 'posterior_predictive.csv'
    pp_df.to_csv(pp_path, index=False)

    # Manifest
    env = base_environment()
    manifest = Manifest(
        run_id=run_id,
        timestamp_utc=run_id.split('_')[0],
        model_variant='bayesian_v3_6',
        ratio=args.ratio,
        era=args.era,
        input_files={'bayesian_inputs': str(input_path)},
        cli_args={
            'draws': str(args.draws), 'tune': str(args.tune), 'chains': str(args.chains), 'seed': str(args.seed),
            'beta_xi_mean': str(args.beta_xi_mean), 'beta_xi_sd': str(args.beta_xi_sd),
            'xi_baseline_nm': str(args.xi_baseline_nm), 'baseline_naacr': str(args.baseline_naa_cr),
            'run_notes': args.run_notes,
        },
        mechanism={'beta_xi_mean': args.beta_xi_mean, 'beta_xi_sd': args.beta_xi_sd, 'xi_baseline_nm': args.xi_baseline_nm},
        data_summary={
            'n_rows': int(len(df)),
            'n_studies': int(df['study'].nunique()),
            'n_phases': int(df['Phase'].nunique()),
            'n_unknown_era': int((df.get('art_era', pd.Series([])) == 'unknown').sum()) if 'art_era' in df.columns else 0,
        },
        environment=env,
        outputs={
            'trace_netcdf': str(trace_path),
            'posterior_predictive_csv': str(pp_path),
            'inputs_annotated_csv': str(inputs_annot_path),
        },
        validation_notes={
            'primary_parameter': 'xi_estimate_nm (mechanistic) — not overridden by era',
            'derived_parameter': 'Pi_xi from xi via beta_xi',
            'era_handling': 'categorical covariate (additive), unknown era retained as 0-effect',
            'enzyme_activity_fold': 'retained in data (if present) for validation only; not used in model',
        }
    )
    write_manifest(out_dir / 'run_manifest.json', manifest)

    print("\n✅ Bayesian v3.6 run complete")
    print(f"  Output: {out_dir}")
    print("  Trace:", trace_path)
    print("  PPD head:\n", pp_df.head().to_string(index=False))


if __name__ == '__main__':
    main()
