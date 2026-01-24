"""
Regional Hierarchical Bayesian Model
====================================

This script implements a hierarchical Bayesian model to test the 
Regional Sensitivity Mapping of the Noise-Mediated Neuroprotection hypothesis.

It uses multi-regional data from Valcour 2015 (BG, FWM, FGM, PGM) 
to determine if the protection factor (beta_xi) or the noise structure (xi) 
varies across brain regions.
"""

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
from pathlib import Path
import argparse
import os

def load_data():
    base_path = Path(__file__).resolve().parent.parent
    ind_path = base_path / 'data/individual/VALCOUR_2015_INDIVIDUAL_PATIENTS.csv'
    
    if not ind_path.exists():
        print(f"Warning: Individual data file not found at {ind_path}")
        return pd.DataFrame(columns=['patient_id', 'region', 'NAA_ratio', 'phase'])
        
    df = pd.read_csv(ind_path)
    
    # Brain regions mentioned in Valcour 2015
    regions = {
        'BG': 'BGNAA',
        'FWM': 'FWMNAA',
        'FGM': 'FGMNAA',
        'PGM': 'PGMNAA'
    }
    
    # References for normalization (from DATA_INVENTORY / Valcour 2015 contexts)
    # Healthy references: BG ~ 8.76, others assumed similar or from general MRS literature
    # We use 8.76 as a global reference to match ratios if others aren't specified.
    ref_naa = 8.76
    
    rows = []
    for _, patient in df.iterrows():
        p_id = _
        for reg_code, col in regions.items():
            if col in df.columns and pd.notna(patient[col]):
                rows.append({
                    'patient_id': p_id,
                    'region': reg_code,
                    'NAA_ratio': patient[col] / ref_naa,
                    'phase': 'acute' # Valcour 2015 is an acute cohort
                })
    
    df_long = pd.DataFrame(rows)
    return df_long

def build_regional_model(df):
    if df.empty:
        # Return a dummy model if data is missing for build phase
        with pm.Model() as model:
            pm.Normal('beta_xi_global', mu=1.8, sigma=0.2)
        return model

    regions = sorted(df['region'].unique())
    region_map = {r: i for i, r in enumerate(regions)}
    region_idx = df['region'].map(region_map).values
    
    unique_patients = sorted(df['patient_id'].unique())
    patient_map = {p: i for i, p in enumerate(unique_patients)}
    patient_idx = df['patient_id'].map(patient_map).values
    n_patients = len(unique_patients)
    
    with pm.Model() as model:
        # Global mechanism prior
        beta_xi_global = pm.Normal('beta_xi_global', mu=1.8, sigma=0.2)
        
        # Region-specific effects on beta_xi (is protection stronger in some regions?)
        beta_xi_offset = pm.Normal('beta_xi_offset', mu=0, sigma=0.2, shape=len(regions))
        beta_xi_region = pm.Deterministic('beta_xi_region', beta_xi_global + beta_xi_offset)
        
        # Noise structure per region (latent xi)
        # Acute phase xi is generally ~0.55-0.70 nm
        xi_acute_region = pm.TruncatedNormal('xi_acute_region', mu=0.6, sigma=0.1, lower=0.3, upper=0.9, shape=len(regions))
        
        # Reference xi
        xi_ref = 0.80
        
        # Hierarchical patient baseline (accounting for individual differences in NAA)
        patient_baseline = pm.Normal('patient_baseline', mu=1.0, sigma=0.1, shape=n_patients)
        
        # Model error
        sigma = pm.HalfNormal('sigma', sigma=0.1)
        
        # Expected NAA ratio
        # mu = patient_baseline * (xi_ref / xi_acute_region[region])^beta_xi_region[region]
        pi_region = (xi_ref / xi_acute_region[region_idx])**beta_xi_region[region_idx]
        mu = patient_baseline[patient_idx] * pi_region
        
        # Likelihood
        pm.Normal('obs', mu=mu, sigma=sigma, observed=df['NAA_ratio'].values)
        
    return model

def parse_args():
    p = argparse.ArgumentParser(description='Regional Hierarchical Bayesian Model Runner')
    p.add_argument('--draws', type=int, default=1000)
    p.add_argument('--tune', type=int, default=1000)
    p.add_argument('--chains', type=int, default=4)
    p.add_argument('--target-accept', type=float, default=0.95)
    p.add_argument('--seed', type=int, default=42)
    return p.parse_args()

def main():
    args = parse_args()
    print("="*60)
    print("Regional Hierarchical Bayesian Model Runner (v1)")
    print("="*60)

    print("Loading regional data from Valcour 2015...")
    df = load_data()
    if df.empty:
        print("ERROR: No data loaded. Check data/individual/VALCOUR_2015_INDIVIDUAL_PATIENTS.csv")
        return

    print(f"Loaded {len(df)} observations across {len(df['region'].unique())} regions.")
    
    print("Building regional hierarchical model...")
    model = build_regional_model(df)
    
    print("Sampling...")
    with model:
        idata = pm.sample(
            draws=args.draws,
            tune=args.tune,
            chains=args.chains,
            target_accept=args.target_accept,
            random_seed=args.seed,
            return_inferencedata=True
        )
    
    print("\n--- REGIONAL RESULTS ---")
    summary = az.summary(idata, var_names=['beta_xi_region', 'xi_acute_region'])
    print(summary)
    
    # Identify region names for clarity
    regions = sorted(df['region'].unique())
    for i, region in enumerate(regions):
        print(f"Index {i}: {region}")
    
    # Save results
    base_path = Path(__file__).resolve().parent.parent
    out_dir = base_path / 'results/regional_hierarchical_v1'
    out_dir.mkdir(parents=True, exist_ok=True)
    az.to_netcdf(idata, str(out_dir / 'trace.nc'))
    summary.to_csv(out_dir / 'summary.csv')
    
    # Generate Visuals
    print("\nGenerating regional visuals...")
    try:
        import matplotlib.pyplot as plt
        
        # Fig 1: Regional beta_xi
        plt.figure(figsize=(10, 6))
        az.plot_forest(idata, var_names=['beta_xi_region'], combined=True)
        plt.axvline(1.0, color='r', linestyle='--', label='Linear')
        plt.title('Regional Protection Exponent (beta_xi)')
        plt.tight_layout()
        plt.savefig(out_dir / 'fig1_regional_beta.png', dpi=200)
        plt.close()
        
        # Fig 2: Regional xi_acute
        plt.figure(figsize=(10, 6))
        az.plot_forest(idata, var_names=['xi_acute_region'], combined=True)
        plt.title('Regional Noise Correlation Length (xi_acute)')
        plt.tight_layout()
        plt.savefig(out_dir / 'fig2_regional_xi.png', dpi=200)
        plt.close()

        # Touch a dummy visuals.png for legacy compatibility
        (out_dir / 'visuals.png').touch()
        
        print(f"Regional visuals saved to {out_dir}")
    except Exception as e:
        print(f"Error generating regional visuals: {e}")

    print(f"\nResults saved to {out_dir}")
    print("="*60)

if __name__ == '__main__':
    main()
