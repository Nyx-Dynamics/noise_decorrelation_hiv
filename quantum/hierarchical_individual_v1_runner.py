"""
Hierarchical Individual-Level Model Runner
==========================================

This script implements a hierarchical Bayesian model that integrates:
1. Individual-level data from Valcour 2015 (n=62 patients, multiple regions)
2. Aggregate study-level data from other landmark cohorts (Sailasuta, Young, Chang)

It formally tests the Noise-Mediated Neuroprotection hypothesis by:
- Modeling individual variability within the Valcour study.
- Modeling study-level differences across cohorts.
- Estimating the latent protection exponent beta_xi across all evidence.
"""

import numpy as np
import pandas as pd
import pymc as pm
import sys
from pathlib import Path

# Check for required packages before importing
try:
    import numpy as np
except ImportError:
    print("ERROR: numpy is not installed or not accessible.")
    print("\nPlease run one of the following:")
    print("  make install              # Install all dependencies")
    print("  pip install numpy         # Install numpy only")
    print("\nOr activate your virtual environment:")
    print("  source .venv/bin/activate")
    print("  # or")
    print("  source noiseenv/bin/activate")
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas is not installed.")
    print("Run: pip install pandas")
    sys.exit(1)

try:
    import pymc as pm
except ImportError:
    print("ERROR: pymc is not installed.")
    print("Run: pip install pymc")
    sys.exit(1)

try:
    import arviz as az
except ImportError:
    print("ERROR: arviz is not installed.")
    print("Run: pip install arviz")
    sys.exit(1)

import argparse
import os

# --- 1. Load Data ---

def load_data():
    # Get project root relative to this script
    base_path = Path(__file__).resolve().parent.parent
    
    # A. Aggregate Data (derived from model_comparison_hierarchical_v1.py)
    aggregate_data = [
        {'study': 'Sailasuta_2012', 'phase': 'acute', 'NAA_mean': 1.134, 'NAA_sd': 0.14, 'n': 31},
        {'study': 'Sailasuta_2012', 'phase': 'chronic', 'NAA_mean': 1.000, 'NAA_sd': 0.14, 'n': 26},
        {'study': 'Sailasuta_2012', 'phase': 'healthy', 'NAA_mean': 1.077, 'NAA_sd': 0.13, 'n': 10},
        {'study': 'Young_2014', 'phase': 'acute', 'NAA_mean': 1.125, 'NAA_sd': 0.20, 'n': 9},
        {'study': 'Young_2014', 'phase': 'chronic', 'NAA_mean': 1.05, 'NAA_sd': 0.15, 'n': 18},
        {'study': 'Young_2014', 'phase': 'healthy', 'NAA_mean': 1.15, 'NAA_sd': 0.15, 'n': 19},
        {'study': 'Chang_2002', 'phase': 'chronic', 'NAA_mean': 7.96 / 8.76, 'NAA_sd': 0.15, 'n': 15},
        {'study': 'Chang_2002', 'phase': 'healthy', 'NAA_mean': 1.0, 'NAA_sd': 0.09, 'n': 15},
    ]
    df_agg = pd.DataFrame(aggregate_data)
    
    # B. Individual Data (Valcour 2015)
    # We use BGNAA as it is the most responsive region mentioned in DATA_INVENTORY
    ind_path = base_path / 'data/individual/VALCOUR_2015_INDIVIDUAL_PATIENTS.csv'
    try:
        if not ind_path.exists():
            print(f"Warning: Individual data file not found at {ind_path}")
            df_ind = pd.DataFrame(columns=['BGNAA', 'NAA_ratio', 'phase', 'study'])
        else:
            df_ind = pd.read_csv(ind_path)
            # Filter for rows that have BGNAA data
            df_ind = df_ind.dropna(subset=['BGNAA'])
            # Valcour 2015 is Acute HIV (n=~62)
            # We need to normalize Valcour data. DATA_INVENTORY says healthy BG NAA is ~8.76 mM
            # Valcour BGNAA is absolute. Let's use 8.76 as reference to match ratios.
            df_ind['NAA_ratio'] = df_ind['BGNAA'] / 8.76
            df_ind['phase'] = 'acute'
            df_ind['study'] = 'Valcour_2015'
    except Exception as e:
        print(f"Error loading individual data: {e}")
        df_ind = pd.DataFrame(columns=['BGNAA', 'NAA_ratio', 'phase', 'study'])

    return df_agg, df_ind

# --- 2. Build Model ---

def build_model(df_agg, df_ind):
    # Combine data indices
    agg_studies = set(df_agg['study'].unique())
    ind_studies = set(df_ind['study'].unique()) if not df_ind.empty else set()
    all_studies = sorted(list(agg_studies | ind_studies))
    study_map = {s: i for i, s in enumerate(all_studies)}
    
    phase_map = {'healthy': 0, 'acute': 1, 'chronic': 2}
    
    # Aggregated obs
    y_agg = df_agg['NAA_mean'].values
    sd_agg = df_agg['NAA_sd'].values
    study_idx_agg = df_agg['study'].map(study_map).values
    phase_idx_agg = df_agg['phase'].map(phase_map).values
    
    # Individual obs
    has_ind = not df_ind.empty
    if has_ind:
        y_ind = df_ind['NAA_ratio'].values
        study_idx_ind = df_ind['study'].map(study_map).values
        phase_idx_ind = df_ind['phase'].map(phase_map).values
    
    with pm.Model() as model:
        # Global mechanism
        beta_xi = pm.Normal('beta_xi', mu=1.9, sigma=0.2)
        
        # Latent xi per phase (nm)
        xi_healthy = pm.TruncatedNormal('xi_healthy', mu=0.66, sigma=0.1, lower=0.4)
        xi_acute = pm.TruncatedNormal('xi_acute', mu=0.55, sigma=0.1, lower=0.3)
        xi_chronic = pm.TruncatedNormal('xi_chronic', mu=0.80, sigma=0.1, lower=0.6)
        
        xi_phases = pm.math.stack([xi_healthy, xi_acute, xi_chronic])
        xi_ref = 0.80
        
        # Protection factor: (xi_ref / xi)^beta_xi
        # We model the observed NAA ratio as: baseline * ProtectionFactor * CoherenceLoss
        # For simplicity, we merge baseline and coherence into a study-level effect
        
        study_baseline = pm.Normal('study_baseline', mu=1.0, sigma=0.1, shape=len(all_studies))
        
        # Observation sigma
        sigma_model = pm.HalfNormal('sigma_model', sigma=0.1)
        
        # --- Likelihood A: Aggregated Data ---
        # Pred: baseline[study] * (xi_ref / xi[phase])^beta_xi
        pi_agg = (xi_ref / xi_phases[phase_idx_agg])**beta_xi
        mu_agg = study_baseline[study_idx_agg] * pi_agg
        
        # Total uncertainty = model error + measurement error (SE = SD/sqrt(n))
        se_agg = sd_agg / np.sqrt(df_agg['n'].values)
        sigma_total_agg = pm.math.sqrt(sigma_model**2 + se_agg**2)
        pm.Normal('obs_agg', mu=mu_agg, sigma=sigma_total_agg, observed=y_agg)
        
        # --- Likelihood B: Individual Data ---
        if has_ind:
            pi_ind = (xi_ref / xi_phases[phase_idx_ind])**beta_xi
            mu_ind = study_baseline[study_idx_ind] * pi_ind
            pm.Normal('obs_ind', mu=mu_ind, sigma=sigma_model, observed=y_ind)
        
        # Deterministics
        pm.Deterministic('xi_diff', xi_chronic - xi_acute)
        pm.Deterministic('p_acute_protected', pm.math.switch(xi_acute < xi_chronic, 1, 0))
        
    return model

# --- 3. Run ---

def parse_args():
    p = argparse.ArgumentParser(description='Hierarchical Individual-Level Model Runner')
    p.add_argument('--draws', type=int, default=2000)
    p.add_argument('--tune', type=int, default=1000)
    p.add_argument('--chains', type=int, default=4)
    p.add_argument('--target-accept', type=float, default=0.95)
    p.add_argument('--seed', type=int, default=42)
    return p.parse_args()

def main():
    args = parse_args()
    print("="*60)
    print("Hierarchical Individual-Level Model Runner (v1)")
    print("="*60)
    print(f"Python: {sys.version}")
    print(f"NumPy: {np.__version__}")
    print(f"Pandas: {pd.__version__}")
    print(f"PyMC: {pm.__version__}")
    print(f"ArviZ: {az.__version__}")
    print("="*60)
    
    print("\nLoading data...")
    df_agg, df_ind = load_data()
    print(f"Loaded {len(df_agg)} aggregate points and {len(df_ind)} individual points.")
    
    print("\nBuilding hierarchical individual model...")
    model = build_model(df_agg, df_ind)
    
    print("\nSampling (this may take several minutes)...")
    with model:
        idata = pm.sample(
            draws=args.draws,
            tune=args.tune,
            chains=args.chains,
            target_accept=args.target_accept,
            random_seed=args.seed,
            return_inferencedata=True
        )
        
    # Analysis
    print("\n--- RESULTS ---")
    summary = az.summary(idata, var_names=['beta_xi', 'xi_healthy', 'xi_acute', 'xi_chronic', 'xi_diff'])
    print(summary)
    
    p_protected = np.mean(idata.posterior['xi_diff'].values > 0)
    print(f"\nProbability ξ_acute < ξ_chronic: {p_protected:.4f}")
    
    # Save results
    # Save relative to script location
    base_path = Path(__file__).resolve().parent.parent
    out_dir = base_path / 'results/hierarchical_individual_v1'
    out_dir.mkdir(parents=True, exist_ok=True)
    az.to_netcdf(idata, str(out_dir / 'trace.nc'))
    summary.to_csv(out_dir / 'summary.csv')
    print(f"\nResults saved to {out_dir}")
    print("="*60)

if __name__ == '__main__':
    main()
