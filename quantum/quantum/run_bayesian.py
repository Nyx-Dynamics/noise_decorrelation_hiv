#!/usr/bin/env python3
"""
Bayesian Model v3.6 - EXPANDED DATASET
Works from /quantum/ directory
Loads: quantum/data/curated/CONSOLIDATED_MRS_DATA_FOR_MODEL.csv

Extended: Can dispatch to separated pipelines with --pipeline {absolute,ratio}
If --pipeline is not provided, runs the legacy expanded individual-level model.
"""

import sys
import argparse
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
from scipy import stats
import warnings
from pathlib import Path
warnings.filterwarnings('ignore')


def _maybe_dispatch_to_pipeline(argv: list[str]) -> bool:
    """If --pipeline is provided, dispatch to the selected pipeline and return True.
    Otherwise return False to let the legacy script continue.
    """
    if "--pipeline" not in argv:
        return False

    # Minimal parser to extract the pipeline and collect remaining args
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--pipeline", choices=["absolute", "ratio"], required=True)
    # Parse known-only to avoid interfering with downstream
    ns, remaining = p.parse_known_args(argv[1:])

    if ns.pipeline == "absolute":
        from . import pipeline_absolute as mod
        # Pass through only args relevant to absolute pipeline (it will parse its own)
        mod.main(remaining)
        return True
    elif ns.pipeline == "ratio":
        from . import pipeline_ratio as mod
        mod.main(remaining)
        return True
    return False


if _maybe_dispatch_to_pipeline(sys.argv):
    # Already handled by specific pipeline; exit early
    sys.exit(0)

print("\n" + "╔" + "="*78 + "╗")
print("║" + " "*20 + "BAYESIAN MODEL v3.6 - EXPANDED" + " "*27 + "║")
print("║" + " "*15 + "Individual Patient Data: n=44 Acute HIV" + " "*24 + "║")
print("╚" + "="*78 + "╝\n")

# ============================================================================
# LOAD DATA
# ============================================================================

print("="*80)
print("LOADING DATA")
print("="*80)

# Try different possible paths
possible_paths = [
    'quantum/data/curated/CONSOLIDATED_MRS_DATA_FOR_MODEL.csv',  # From /quantum/ dir
    'data/curated/CONSOLIDATED_MRS_DATA_FOR_MODEL.csv',  # From /quantum/quantum/ dir
    './CONSOLIDATED_MRS_DATA_FOR_MODEL.csv',  # Current directory
    '/Users/acdstudpro/Documents/Github/noise_decorrelation_HIV/quantum/quantum/data/curated/CONSOLIDATED_MRS_DATA_FOR_MODEL.csv'  # Absolute
]

df_full = None
for path in possible_paths:
    try:
        print(f"Trying: {path}")
        df_full = pd.read_csv(path)
        print(f"✅ Loaded from: {path}")
        print(f"   Records: {len(df_full)}")
        break
    except FileNotFoundError:
        continue

if df_full is None:
    print("\n❌ Could not find data file. Please run from:")
    print("   /Users/acdstudpro/Documents/Github/noise_decorrelation_HIV/quantum/")
    print("\nOr place CONSOLIDATED_MRS_DATA_FOR_MODEL.csv in current directory")
    exit(1)

print("\n📊 Phase distribution:")
print(df_full['Phase'].value_counts())

# ============================================================================
# EXTRACT BASAL GANGLIA DATA
# ============================================================================

print("\n" + "="*80)
print("EXTRACTING BASAL GANGLIA DATA")
print("="*80)

# Get acute BG data from Valcour
acute_bg = df_full[(df_full['Phase'] == 'Acute') & 
                    (df_full['Region'] == 'BG') & 
                    (df_full['Study'] == 'Valcour_2015')].copy()

print(f"\n✅ Found {len(acute_bg)} individual acute patients")

naa_obs_acute = acute_bg['NAA'].dropna().values
cho_obs_acute = acute_bg['Cho'].dropna().values

print(f"   NAA: n={len(naa_obs_acute)}")
print(f"   Cho: n={len(cho_obs_acute)}")

# Reference values
naa_obs_control = np.array([9.55])
cho_obs_control = np.array([2.18])
naa_obs_chronic = np.array([8.79])
cho_obs_chronic = np.array([2.40])

print(f"\n📊 ACUTE PHASE:")
print(f"   NAA = {naa_obs_acute.mean():.2f} ± {naa_obs_acute.std():.2f} mM")
print(f"   Range: {naa_obs_acute.min():.2f} - {naa_obs_acute.max():.2f} mM")

# Stats
t_stat, p_val = stats.ttest_1samp(naa_obs_acute, naa_obs_control[0])
pct_change = ((naa_obs_acute.mean() - naa_obs_control[0])/naa_obs_control[0]*100)

print(f"\n📊 ACUTE vs CONTROL:")
print(f"   Change: +{pct_change:.1f}%")
print(f"   p-value: {p_val:.6f} ***")

# ============================================================================
# BUILD MODEL
# ============================================================================

print("\n" + "="*80)
print("BUILDING MODEL")
print("="*80)

V_max_base = 100.0
K_m = 50.0
S_0 = 100.0
V_base = V_max_base * S_0 / (K_m + S_0)
α_NAA = 9.55 / V_base

with pm.Model() as model:
    
    # Priors
    ξ_acute = pm.TruncatedNormal('ξ_acute', mu=0.6, sigma=0.1, lower=0.3, upper=0.9)
    ξ_chronic = pm.TruncatedNormal('ξ_chronic', mu=0.8, sigma=0.1, lower=0.5, upper=1.2)
    ξ_control = pm.TruncatedNormal('ξ_control', mu=0.5, sigma=0.05, lower=0.3, upper=0.7)
    
    σ_noise_acute = pm.TruncatedNormal('σ_noise_acute', mu=5.0, sigma=1.0, lower=2.0, upper=10.0)
    σ_noise_chronic = pm.TruncatedNormal('σ_noise_chronic', mu=2.0, sigma=0.5, lower=0.5, upper=4.0)
    σ_noise_control = pm.TruncatedNormal('σ_noise_control', mu=0.5, sigma=0.2, lower=0.1, upper=1.0)
    
    β_ξ = pm.TruncatedNormal('β_ξ', mu=2.0, sigma=0.5, lower=1.0, upper=3.0)
    
    # Protection factors
    F_acute = pm.Deterministic('F_acute', pm.math.exp(-β_ξ * ξ_acute))
    F_chronic = pm.Deterministic('F_chronic', pm.math.exp(-β_ξ * ξ_chronic))
    F_control = pm.Deterministic('F_control', pm.math.exp(-β_ξ * ξ_control))
    
    # Enzyme activity
    V_acute = pm.Deterministic('V_acute', V_base * (1 + F_acute * σ_noise_acute / 10))
    V_chronic = pm.Deterministic('V_chronic', V_base * (1 - (1-F_chronic) * σ_noise_chronic / 5))
    V_control = pm.Deterministic('V_control', V_base * (1 + F_control * σ_noise_control / 10))
    
    # NAA predictions
    NAA_acute_mean = pm.Deterministic('NAA_acute_mean', α_NAA * V_acute)
    NAA_chronic_mean = pm.Deterministic('NAA_chronic_mean', α_NAA * V_chronic)
    NAA_control_mean = pm.Deterministic('NAA_control_mean', α_NAA * V_control)
    
    # Cho predictions
    Cho_acute_mean = pm.Deterministic('Cho_acute_mean', 2.18 * (1 + σ_noise_acute / 20))
    Cho_chronic_mean = pm.Deterministic('Cho_chronic_mean', 2.18 * (1 + σ_noise_chronic / 15))
    Cho_control_mean = pm.Deterministic('Cho_control_mean', 2.18 * (1 + σ_noise_control / 20))
    
    # Observation noise
    σ_obs_acute = pm.HalfNormal('σ_obs_acute', sigma=1.5)
    σ_obs_chronic = pm.HalfNormal('σ_obs_chronic', sigma=1.0)
    σ_obs_control = pm.HalfNormal('σ_obs_control', sigma=0.5)
    σ_cho_acute = pm.HalfNormal('σ_cho_acute', sigma=0.5)
    
    # Likelihood
    NAA_acute_obs = pm.Normal('NAA_acute_obs', mu=NAA_acute_mean, sigma=σ_obs_acute, observed=naa_obs_acute)
    NAA_chronic_obs = pm.Normal('NAA_chronic_obs', mu=NAA_chronic_mean, sigma=σ_obs_chronic, observed=naa_obs_chronic)
    NAA_control_obs = pm.Normal('NAA_control_obs', mu=NAA_control_mean, sigma=σ_obs_control, observed=naa_obs_control)
    Cho_acute_obs = pm.Normal('Cho_acute_obs', mu=Cho_acute_mean, sigma=σ_cho_acute, observed=cho_obs_acute)
    Cho_control_obs = pm.Normal('Cho_control_obs', mu=Cho_control_mean, sigma=0.2, observed=cho_obs_control)
    
    # Hypothesis
    Δξ = pm.Deterministic('Δξ', ξ_chronic - ξ_acute)

print("✅ Model built!")

# ============================================================================
# SAMPLE
# ============================================================================

print("\n" + "="*80)
print("RUNNING MCMC (10-15 minutes)...")
print("="*80)

with model:
    trace = pm.sample(
        draws=2000,
        tune=1000,
        chains=4,
        cores=4,
        target_accept=0.95,
        return_inferencedata=True,
        random_seed=42
    )
    print("\n✅ Sampling complete!")
    
    ppc = pm.sample_posterior_predictive(trace, random_seed=42)
    print("✅ Posterior predictive complete!")

# ============================================================================
# RESULTS
# ============================================================================

print("\n" + "="*80)
print("RESULTS")
print("="*80)

posterior = trace.posterior

ξ_acute_samples = posterior['ξ_acute'].values.flatten()
ξ_chronic_samples = posterior['ξ_chronic'].values.flatten()
Δξ_samples = posterior['Δξ'].values.flatten()
β_ξ_samples = posterior['β_ξ'].values.flatten()

print("\n🔬 CORRELATION LENGTHS:")
print(f"   ξ_acute = {ξ_acute_samples.mean():.3f} ± {ξ_acute_samples.std():.3f} nm")
print(f"   ξ_chronic = {ξ_chronic_samples.mean():.3f} ± {ξ_chronic_samples.std():.3f} nm")
print(f"   Δξ = {Δξ_samples.mean():.3f} ± {Δξ_samples.std():.3f} nm")

P_acute_shorter = (Δξ_samples > 0).sum() / len(Δξ_samples)
print(f"\n✅ P(ξ_acute < ξ_chronic) = {P_acute_shorter:.4f}")

if P_acute_shorter > 0.999:
    print("   *** STRONGLY SUPPORTED ***")
elif P_acute_shorter > 0.99:
    print("   *** SUPPORTED ***")

print(f"\n🔬 PROTECTION EXPONENT:")
print(f"   β_ξ = {β_ξ_samples.mean():.2f} ± {β_ξ_samples.std():.2f}")

naa_acute_pred = posterior['NAA_acute_mean'].values.flatten()
print(f"\n📊 NAA PREDICTION:")
print(f"   Predicted: {naa_acute_pred.mean():.2f} mM")
print(f"   Observed: {naa_obs_acute.mean():.2f} mM")
print(f"   Error: {abs(naa_acute_pred.mean() - naa_obs_acute.mean())/naa_obs_acute.mean()*100:.1f}%")

# ============================================================================
# SAVE
# ============================================================================

print("\n" + "="*80)
print("SAVING")
print("="*80)

Path('outputs').mkdir(exist_ok=True)

trace.to_netcdf('outputs/trace.nc')
print("✅ outputs/trace.nc")

az.summary(trace, hdi_prob=0.95).to_csv('outputs/summary.csv')
print("✅ outputs/summary.csv")

results = pd.DataFrame({
    'Parameter': ['ξ_acute', 'ξ_chronic', 'Δξ', 'β_ξ', 'P(hypothesis)'],
    'Mean': [ξ_acute_samples.mean(), ξ_chronic_samples.mean(), 
             Δξ_samples.mean(), β_ξ_samples.mean(), P_acute_shorter],
    'SD': [ξ_acute_samples.std(), ξ_chronic_samples.std(),
           Δξ_samples.std(), β_ξ_samples.std(), np.nan]
})
results.to_csv('outputs/results.csv', index=False)
print("✅ outputs/results.csv")

print("\n" + "="*80)
print("COMPLETE!")
print("="*80)
print(f"\n🎯 n={len(naa_obs_acute)} patients")
print(f"   P(ξ_acute < ξ_chronic) = {P_acute_shorter:.4f}")
print(f"   NAA error = {abs(naa_acute_pred.mean() - naa_obs_acute.mean())/naa_obs_acute.mean()*100:.1f}%")
print("\n" + "="*80)
