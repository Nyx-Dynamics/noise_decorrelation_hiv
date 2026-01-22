"""
Bayesian v3 runner with ART-era hierarchical covariate and Option A (xi → Pi_xi).

- Loads bayesian_inputs_<ratio>.csv via quantum.data_loaders
- Attaches ART era (pre ≤2006; post ≥2007; unknown retained unless filtered out)
- Mechanistic core preserved: xi_estimate_nm → Pi_xi using beta_xi prior N(1.89, 0.25)
- Adds additive ART-era effect and study random effects; does NOT alter xi→Pi_xi mapping
- Writes outputs per run (non-overwriting):
    results/bayesian_v3/<ratio>/<era>/<run_id>/
      - trace.nc
      - posterior_predictive.csv
      - inputs_annotated.csv
      - run_manifest.json
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

try:
    import pymc as pm  # type: ignore
    import arviz as az  # type: ignore
except Exception as e:  # pragma: no cover - handled at runtime
    pm = None
    az = None

from quantum.data_loaders import load_bayesian_inputs
from quantum.utils.run_manifest import (
    make_run_id,
    base_environment,
    write_manifest,
    Manifest,
)


PHASE_NORMALIZE: Dict[str, str] = {
    'control': 'control', 'Control': 'control', 'Healthy': 'control', 'healthy': 'control',
    'Acute': 'acute', 'acute': 'acute',
    'Chronic': 'chronic', 'chronic': 'chronic'
}


def parse_args():
    p = argparse.ArgumentParser(description='Bayesian v3 runner (Option A; ART-era as covariate)')
    p.add_argument('--ratio', required=True, choices=['3_1_1', '1_1_1', '1_2_1', '2_1'])
    p.add_argument('--era', default='both', choices=['pre_modern', 'post_modern', 'both'])
    p.add_argument('--output-root', default='results/bayesian_v3')
    p.add_argument('--draws', type=int, default=1500)
    p.add_argument('--tune', type=int, default=1000)
    p.add_argument('--chains', type=int, default=4)
    p.add_argument('--target-accept', type=float, default=0.9)
    p.add_argument('--seed', type=int, default=2025)
    p.add_argument('--beta-xi-mean', type=float, default=1.89)
    p.add_argument('--beta-xi-sd', type=float, default=0.25)
    p.add_argument('--xi-baseline-nm', type=float, default=0.66)
    p.add_argument('--run-notes', default='')
    # Legacy mirroring (parity with v3.6)
    p.add_argument('--legacy-output-root', default=None, help='Optional path to mirror this run (copy/symlink) into an external folder (e.g., iCloud runs/run_YYYYMMDD_HHMMSS)')
    p.add_argument('--legacy-mirror-mode', default='copy', choices=['copy', 'symlink'], help='Mirror mode for legacy-output-root: copy files or create symlinks')
    return p.parse_args()


def _prepare_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    # Normalize phase/condition fields
    if 'phase' in d.columns:
        d['phase_norm'] = d['phase'].map(PHASE_NORMALIZE).fillna(d['phase'])
    elif 'Phase' in d.columns:
        d['phase_norm'] = d['Phase'].map(PHASE_NORMALIZE).fillna(d['Phase'])
    else:
        # Try Condition
        d['phase_norm'] = d.get('Condition', 'control')
        d['phase_norm'] = d['phase_norm'].map(PHASE_NORMALIZE).fillna(d['phase_norm'])

    # Auto-detect schema and filter to NAA/Cr metabolite when present
    if 'Metabolite' in d.columns:
        meta_norm = d['Metabolite'].astype(str).str.replace(' ', '', regex=False).str.upper()
        naa_mask = meta_norm.eq('NAA/CR') | meta_norm.eq('NAA/CRR') | meta_norm.eq('NAA/CRATIO')
        if naa_mask.any():
            d = d[naa_mask].copy()

    # Study identifiers
    if 'study' not in d.columns:
        # attempt to infer from provenance or file columns
        prov_col = None
        for c in ['study', 'Study', 'provenance', 'source_file', 'SourceFile']:
            if c in d.columns:
                prov_col = c
                break
        if prov_col is not None:
            d['study'] = d[prov_col].astype(str).str.extract(r'([^/\\]+)')[0].fillna('UnknownStudy')
        else:
            d['study'] = 'UnknownStudy'

    d['study_id'] = d['study'].astype('category').cat.codes

    # Observed means/SE: try common column names
    # Preferred NAA/Cr and Cho/Cr columns
    naa_cols = ['NAA_Cr', 'NAA_over_Cr', 'NAA_over_Cr_mean', 'naa_cr_mean']
    cho_cols = ['Cho_Cr', 'Cho_over_Cr', 'Cho_over_Cr_mean', 'cho_cr_mean']

    # Initialize targets
    d['naa_cr_mean'] = np.nan
    d['cho_cr_mean'] = np.nan
    d['naa_cr_se'] = np.nan

    # If generic Mean/SE present (and we filtered to NAA/Cr), map them
    if {'Mean', 'SE'}.issubset(d.columns):
        d['naa_cr_mean'] = pd.to_numeric(d['Mean'], errors='coerce')
        d['naa_cr_se'] = pd.to_numeric(d['SE'], errors='coerce')
    else:
        # Try specific columns
        for c in naa_cols:
            if c in d.columns:
                d['naa_cr_mean'] = pd.to_numeric(d[c], errors='coerce')
                break
        for c in cho_cols:
            if c in d.columns:
                d['cho_cr_mean'] = pd.to_numeric(d[c], errors='coerce')
                break
        # SE columns variants
        for c in ['SE', 'SE_NAA_Cr', 'naa_cr_se']:
            if c in d.columns:
                d['naa_cr_se'] = pd.to_numeric(d[c], errors='coerce')
                break

    # If standard error missing, set a weakly-informative default to avoid zero-variance
    if not np.isfinite(d['naa_cr_se']).any():
        d['naa_cr_se'] = 0.1

    # Encode era
    # art_era and art_era_idx already attached by loader; ensure availability
    if 'art_era' not in d.columns:
        d['art_era'] = 'unknown'
    if 'art_era_idx' not in d.columns:
        d['art_era_idx'] = d['art_era'].map({'pre_modern': 0, 'post_modern': 1, 'unknown': -1})

    return d


def _filter_by_era(df: pd.DataFrame, era: str) -> pd.DataFrame:
    if era == 'both':
        return df
    return df[df['art_era'] == era]


def build_model(data: pd.DataFrame, beta_xi_mean: float, beta_xi_sd: float, xi_baseline_nm: float):
    """
    Bayesian v3 model:
      - Mechanistic core: xi→Pi_xi via beta_xi prior
      - Per-phase baseline xi (latent), mapped to NAA/Cr via Pi_xi
      - ART era additive effect (2 levels) on the mean outcome (does not modify Pi_xi)
      - Study random effects for heterogeneity
    """
    if pm is None:
        raise RuntimeError("pymc is not available. Please install pymc to run Bayesian inference.")

    # Indexing helpers
    phase_cat = pd.Categorical(data['phase_norm'], categories=['control', 'acute', 'chronic'])
    phase_idx = phase_cat.codes

    study_idx = data['study_id'].values.astype(int)
    n_studies = int(data['study_id'].nunique())

    era_idx = data['art_era_idx'].values.astype(int)
    # Replace unknown (-1) with 0 (pre) for effect indexing but keep a mask
    unknown_mask = era_idx < 0
    era_idx_safe = era_idx.copy()
    era_idx_safe[unknown_mask] = 0

    observed = data['naa_cr_mean'].values.astype(float)
    se = data['naa_cr_se'].values.astype(float)

    with pm.Model() as model:
        # Mechanism parameter
        beta_xi = pm.Normal('beta_xi', mu=beta_xi_mean, sigma=beta_xi_sd)

        # Latent per-phase xi (in nm) with weak priors around plausible values
        xi_control_nm = pm.Normal('xi_control_nm', mu=0.66, sigma=0.08)
        xi_acute_nm   = pm.Normal('xi_acute_nm',   mu=0.61, sigma=0.08)
        xi_chronic_nm = pm.Normal('xi_chronic_nm', mu=0.80, sigma=0.08)
        xi_by_phase_nm = pm.math.stack([xi_control_nm, xi_acute_nm, xi_chronic_nm])

        # Map xi to Pi_xi
        Pi_xi_by_phase = (xi_by_phase_nm / xi_baseline_nm) ** (-beta_xi)

        # Baseline NAA/Cr scale parameter (healthy reference), weak prior
        baseline_ratio = pm.Normal('baseline_ratio', mu=2.2, sigma=0.3)

        # Era additive effect (does not change Pi_xi)
        era_effect = pm.Normal('era_effect', mu=0.0, sigma=0.15, shape=2)  # [pre, post]

        # Study random effects
        study_sd = pm.HalfNormal('study_sd', sigma=0.20)
        study_raw = pm.Normal('study_raw', mu=0.0, sigma=1.0, shape=n_studies)
        study_effect = pm.Deterministic('study_effect', study_raw * study_sd)

        # Expected mean per observation
        Pi_obs = Pi_xi_by_phase[phase_idx]
        mu = baseline_ratio * Pi_obs * (1.0 + era_effect[era_idx_safe] + study_effect[study_idx])

        # Observation noise combines reported SE (when present) with additional residual
        sigma_resid = pm.HalfNormal('sigma_resid', sigma=0.15)
        sigma = pm.math.sqrt(sigma_resid**2 + se**2)

        pm.Normal('naa_cr_obs', mu=mu, sigma=sigma, observed=observed)

        # Posterior predictive for fitted points
        pm.Deterministic('mu_fitted', mu)

    return model


def run_inference(model, draws: int, tune: int, chains: int, target_accept: float, seed: int):
    with model:
        trace = pm.sample(draws=draws, tune=tune, chains=chains, target_accept=target_accept,
                          random_seed=seed, return_inferencedata=True, progressbar=True)
        pm.sample_posterior_predictive(trace, extend_inferencedata=True)
    return trace


def save_outputs(trace, data: pd.DataFrame, out_dir: Path) -> Tuple[Path, Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    # Trace
    trace_path = out_dir / 'trace.nc'
    az.to_netcdf(trace, trace_path)

    # Posterior predictive summary (means per observation)
    mu_post = trace.posterior['mu_fitted']  # dims: chain, draw, obs
    mu_mean = mu_post.mean(dim=['chain', 'draw']).values
    pp_df = data.copy()
    pp_df['mu_fitted_mean'] = mu_mean
    pp_path = out_dir / 'posterior_predictive.csv'
    pp_df.to_csv(pp_path, index=False)

    # Save annotated inputs
    inputs_path = out_dir / 'inputs_annotated.csv'
    data.to_csv(inputs_path, index=False)

    return trace_path, pp_path, inputs_path


def main():
    args = parse_args()
    # Load and filter inputs
    df, input_path = load_bayesian_inputs(args.ratio)
    df = _prepare_dataframe(df)
    if args.era != 'both':
        df = _filter_by_era(df, args.era)

    # Prepare output directory
    run_id = make_run_id()
    out_dir = Path(args.output_root) / args.ratio / args.era / run_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build & run model
    model = build_model(df, beta_xi_mean=args.beta_xi_mean, beta_xi_sd=args.beta_xi_sd,
                        xi_baseline_nm=args.xi_baseline_nm)
    trace = run_inference(model, draws=args.draws, tune=args.tune, chains=args.chains,
                          target_accept=args.target_accept, seed=args.seed)

    # Save outputs
    trace_path, pp_path, inputs_path = save_outputs(trace, df, out_dir)

    # Manifest
    env = base_environment()
    manifest = Manifest(
        run_id=run_id,
        timestamp_utc=run_id.split('_')[0],
        model_variant='bayesian_v3',
        ratio=args.ratio,
        era=args.era,
        input_files={'bayesian_inputs': str(input_path)},
        cli_args={
            'draws': str(args.draws), 'tune': str(args.tune), 'chains': str(args.chains),
            'target_accept': str(args.target_accept), 'seed': str(args.seed),
            'beta_xi_mean': str(args.beta_xi_mean), 'beta_xi_sd': str(args.beta_xi_sd),
            'xi_baseline_nm': str(args.xi_baseline_nm), 'run_notes': args.run_notes,
        },
        mechanism={'beta_xi_mean': args.beta_xi_mean, 'beta_xi_sd': args.beta_xi_sd,
                   'xi_baseline_nm': args.xi_baseline_nm},
        data_summary={
            'n_rows': int(len(df)),
            'n_unknown_era': int((df.get('art_era', pd.Series([])) == 'unknown').sum()),
            'n_studies': int(df['study_id'].nunique()),
        },
        environment=env,
        outputs={
            'trace_nc': str(trace_path),
            'posterior_predictive_csv': str(pp_path),
            'inputs_annotated_csv': str(inputs_path),
        },
        validation_notes={
            'primary_parameter': 'xi_estimate_nm',
            'derived_parameter': 'Pi_xi computed from xi via β_ξ (Option A)',
            'era_handling': 'categorical covariate; does not modify xi→Pi_xi',
        }
    )
    write_manifest(out_dir / 'run_manifest.json', manifest)

    print("\n✨ Bayesian v3 run complete")
    print(f"  Output: {out_dir}")


if __name__ == '__main__':
    main()
