#!/usr/bin/env python3
"""
Bayesian Model v3.6 - EXPANDED DATASET (Mac Local Version)
Using n=44 individual acute patients + group means for chronic/control

Uses your local file: data/curated/CONSOLIDATED_MRS_DATA_FOR_MODEL.csv
Updated: November 15, 2025
Investigator: AC (Nyx Dynamics LLC)
"""

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
from scipy import stats
import warnings
import os
from pathlib import Path

warnings.filterwarnings('ignore')

print("\n" + "╔" + "=" * 78 + "╗")
print("║" + " " * 20 + "BAYESIAN MODEL v3.6 - EXPANDED" + " " * 27 + "║")
print("║" + " " * 78 + "║")
print("║" + " " * 15 + "Individual Patient Data: n=44 Acute HIV" + " " * 24 + "║")
print("╚" + "=" * 78 + "╝\n")

# ============================================================================
# LOAD DATA FROM YOUR LOCAL FILE
# ============================================================================

print("=" * 80)
print("LOADING DATA FROM LOCAL FILE")
print("=" * 80)

# Path to your data file
data_path = Path(__file__).parent / 'data' / 'curated' / 'CONSOLIDATED_MRS_DATA_FOR_MODEL.csv'

# Alternative: use absolute path if above doesn't work
# data_path = '/Users/acdstudpro/Documents/Github/noise_decorrelation_HIV/quantum/quantum/data/curated/CONSOLIDATED_MRS_DATA_FOR_MODEL.csv'

print(f"\n📂 Loading from: {data_path}")

try:
    df_full = pd.read_csv(data_path)
    print(f"✅ Loaded {len(df_full)} records")
except FileNotFoundError:
    print(f"❌ File not found at: {data_path}")
    print("\nTrying absolute path...")
    data_path = '/Users/acdstudpro/Documents/Github/noise_decorrelation_HIV/quantum/quantum/data/curated/CONSOLIDATED_MRS_DATA_FOR_MODEL.csv'
    df_full = pd.read_csv(data_path)
    print(f"✅ Loaded {len(df_full)} records from absolute path")

print("\n📊 Data preview:")
print(df_full.head())

print("\n📊 Phase distribution:")
print(df_full['Phase'].value_counts())

# ============================================================================
# EXTRACT BASAL GANGLIA DATA
# ============================================================================

print("\n" + "=" * 80)
print("EXTRACTING BASAL GANGLIA DATA")
print("=" * 80)

# Get individual acute BG data from Valcour
acute_bg = df_full[(df_full['Phase'] == 'Acute') &
                   (df_full['Region'] == 'BG') &
                   (df_full['Study'] == 'Valcour_2015')].copy()

print(f"\n✅ Found {len(acute_bg)} individual acute patients")

# Extract NAA and Cho values
naa_obs_acute = acute_bg['NAA'].dropna().values
cho_obs_acute = acute_bg['Cho'].dropna().values

print(f"   NAA measurements: n={len(naa_obs_acute)}")
print(f"   Cho measurements: n={len(cho_obs_acute)}")

# Get chronic and control reference values
# Use literature references from meta-analyses
naa_obs_control = np.array([9.55])  # BG NAA reference (mM)
cho_obs_control = np.array([2.18])  # BG Cho reference (mM)

# Chronic phase from meta-analysis (Dahmani 2021 shows ~8% decline)
naa_obs_chronic = np.array([8.79])  # 9.55 * 0.92
cho_obs_chronic = np.array([2.40])  # Elevated from baseline

print(f"\n📊 BASAL GANGLIA DATA SUMMARY:")
print(f"\n   ACUTE (Individual Patients):")
print(f"   n = {len(naa_obs_acute)}")
print(f"   NAA = {naa_obs_acute.mean():.2f} ± {naa_obs_acute.std():.2f} mM")
print(f"   Cho = {cho_obs_acute.mean():.2f} ± {cho_obs_acute.std():.2f} mM")
print(f"   Range: {naa_obs_acute.min():.2f} - {naa_obs_acute.max():.2f} mM")

print(f"\n   CHRONIC (Reference):")
print(f"   NAA = {naa_obs_chronic[0]:.2f} mM")
print(f"   Cho = {cho_obs_chronic[0]:.2f} mM")

print(f"\n   CONTROL (Reference):")
print(f"   NAA = {naa_obs_control[0]:.2f} mM")
print(f"   Cho = {cho_obs_control[0]:.2f} mM")

# Statistical comparison
t_stat, p_val = stats.ttest_1samp(naa_obs_acute, naa_obs_control[0])
cohens_d = (naa_obs_acute.mean() - naa_obs_control[0]) / naa_obs_acute.std()

print(f"\n📊 ACUTE vs CONTROL COMPARISON:")
print(f"   NAA difference: +{naa_obs_acute.mean() - naa_obs_control[0]:.2f} mM")
print(f"   Percent change: +{((naa_obs_acute.mean() - naa_obs_control[0]) / naa_obs_control[0] * 100):.1f}%")
print(f"   t-statistic: {t_stat:.3f}")
print(f"   p-value: {p_val:.6f} ***")
print(f"   Cohen's d: {cohens_d:.2f}")

if p_val < 0.001:
    print("   *** HIGHLY SIGNIFICANT ***")

# ============================================================================
# MODEL PARAMETERS
# ============================================================================

print("\n" + "=" * 80)
print("MODEL PARAMETERS")
print("=" * 80)

# Physical constants
L_MT = 2000e-9  # Microtubule length (m)
k_B = 1.38e-23  # Boltzmann constant
T = 310  # Temperature (K)
hbar = 1.055e-34  # Reduced Planck constant

# Enzyme kinetics
V_max_base = 100.0  # Base enzyme velocity (nmol/min/mg)
K_m = 50.0  # Michaelis constant (μM)
S_0 = 100.0  # Substrate concentration (μM)

print("\n🔬 Physical Constants:")
print(f"   Microtubule length: {L_MT * 1e9:.0f} nm")
print(f"   Temperature: {T} K")
print(f"   Enzyme V_max baseline: {V_max_base} nmol/min/mg")

# ============================================================================
# BAYESIAN MODEL v3.6
# ============================================================================

print("\n" + "=" * 80)
print("BUILDING BAYESIAN MODEL v3.6")
print("=" * 80)

with pm.Model() as model:
    # PRIORS: Noise correlation lengths (KEY PARAMETERS)
    ξ_acute = pm.TruncatedNormal('ξ_acute', mu=0.6, sigma=0.1,
                                 lower=0.3, upper=0.9)

    ξ_chronic = pm.TruncatedNormal('ξ_chronic', mu=0.8, sigma=0.1,
                                   lower=0.5, upper=1.2)

    ξ_control = pm.TruncatedNormal('ξ_control', mu=0.5, sigma=0.05,
                                   lower=0.3, upper=0.7)

    # Noise amplitude parameters
    σ_noise_acute = pm.TruncatedNormal('σ_noise_acute', mu=5.0, sigma=1.0,
                                       lower=2.0, upper=10.0)

    σ_noise_chronic = pm.TruncatedNormal('σ_noise_chronic', mu=2.0, sigma=0.5,
                                         lower=0.5, upper=4.0)

    σ_noise_control = pm.TruncatedNormal('σ_noise_control', mu=0.5, sigma=0.2,
                                         lower=0.1, upper=1.0)

    # Protection factor exponent
    β_ξ = pm.TruncatedNormal('β_ξ', mu=2.0, sigma=0.5,
                             lower=1.0, upper=3.0)

    # Coherence preservation factors
    F_acute = pm.Deterministic('F_acute', pm.math.exp(-β_ξ * ξ_acute))
    F_chronic = pm.Deterministic('F_chronic', pm.math.exp(-β_ξ * ξ_chronic))
    F_control = pm.Deterministic('F_control', pm.math.exp(-β_ξ * ξ_control))

    # Enzyme activity
    V_base = V_max_base * S_0 / (K_m + S_0)

    V_acute = pm.Deterministic('V_acute',
                               V_base * (1 + F_acute * σ_noise_acute / 10))
    V_chronic = pm.Deterministic('V_chronic',
                                 V_base * (1 - (1 - F_chronic) * σ_noise_chronic / 5))
    V_control = pm.Deterministic('V_control',
                                 V_base * (1 + F_control * σ_noise_control / 10))

    # NAA predictions
    α_NAA = 9.55 / V_base

    NAA_acute_mean = pm.Deterministic('NAA_acute_mean', α_NAA * V_acute)
    NAA_chronic_mean = pm.Deterministic('NAA_chronic_mean', α_NAA * V_chronic)
    NAA_control_mean = pm.Deterministic('NAA_control_mean', α_NAA * V_control)

    # Cho predictions
    Cho_acute_mean = pm.Deterministic('Cho_acute_mean',
                                      2.18 * (1 + σ_noise_acute / 20))
    Cho_chronic_mean = pm.Deterministic('Cho_chronic_mean',
                                        2.18 * (1 + σ_noise_chronic / 15))
    Cho_control_mean = pm.Deterministic('Cho_control_mean',
                                        2.18 * (1 + σ_noise_control / 20))

    # Observation noise
    σ_obs_acute = pm.HalfNormal('σ_obs_acute', sigma=1.5)
    σ_obs_chronic = pm.HalfNormal('σ_obs_chronic', sigma=1.0)
    σ_obs_control = pm.HalfNormal('σ_obs_control', sigma=0.5)
    σ_cho_acute = pm.HalfNormal('σ_cho_acute', sigma=0.5)

    # Likelihood
    NAA_acute_obs = pm.Normal('NAA_acute_obs',
                              mu=NAA_acute_mean,
                              sigma=σ_obs_acute,
                              observed=naa_obs_acute)

    NAA_chronic_obs = pm.Normal('NAA_chronic_obs',
                                mu=NAA_chronic_mean,
                                sigma=σ_obs_chronic,
                                observed=naa_obs_chronic)

    NAA_control_obs = pm.Normal('NAA_control_obs',
                                mu=NAA_control_mean,
                                sigma=σ_obs_control,
                                observed=naa_obs_control)

    Cho_acute_obs = pm.Normal('Cho_acute_obs',
                              mu=Cho_acute_mean,
                              sigma=σ_cho_acute,
                              observed=cho_obs_acute)

    Cho_control_obs = pm.Normal('Cho_control_obs',
                                mu=Cho_control_mean,
                                sigma=0.2,
                                observed=cho_obs_control)

    # Hypothesis test
    Δξ = pm.Deterministic('Δξ', ξ_chronic - ξ_acute)

print("\n✅ Model built successfully!")
print(f"   Free parameters: {len(model.free_RVs)}")
print(
    f"   NAA observations: {len(naa_obs_acute)} acute + {len(naa_obs_chronic)} chronic + {len(naa_obs_control)} control")

# ============================================================================
# INFERENCE
# ============================================================================

print("\n" + "=" * 80)
print("RUNNING MCMC INFERENCE")
print("=" * 80)

with model:
    print("\n⏳ Sampling (this may take 5-10 minutes)...")
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

    print("\n⏳ Generating posterior predictive samples...")
    ppc = pm.sample_posterior_predictive(trace, random_seed=42)
    print("✅ Posterior predictive complete!")

# ============================================================================
# RESULTS
# ============================================================================

print("\n" + "=" * 80)
print("BAYESIAN INFERENCE RESULTS")
print("=" * 80)

posterior = trace.posterior

ξ_acute_samples = posterior['ξ_acute'].values.flatten()
ξ_chronic_samples = posterior['ξ_chronic'].values.flatten()
Δξ_samples = posterior['Δξ'].values.flatten()
β_ξ_samples = posterior['β_ξ'].values.flatten()

print("\n🔬 CORRELATION LENGTH PARAMETERS:")
print(f"\n   ξ_acute = {ξ_acute_samples.mean():.3f} ± {ξ_acute_samples.std():.3f} nm")
print(f"   95% HDI: [{np.percentile(ξ_acute_samples, 2.5):.3f}, {np.percentile(ξ_acute_samples, 97.5):.3f}]")

print(f"\n   ξ_chronic = {ξ_chronic_samples.mean():.3f} ± {ξ_chronic_samples.std():.3f} nm")
print(f"   95% HDI: [{np.percentile(ξ_chronic_samples, 2.5):.3f}, {np.percentile(ξ_chronic_samples, 97.5):.3f}]")

print(f"\n   Δξ = ξ_chronic - ξ_acute = {Δξ_samples.mean():.3f} ± {Δξ_samples.std():.3f} nm")
print(f"   95% HDI: [{np.percentile(Δξ_samples, 2.5):.3f}, {np.percentile(Δξ_samples, 97.5):.3f}]")

P_acute_shorter = (Δξ_samples > 0).sum() / len(Δξ_samples)
print(f"\n✅ P(ξ_acute < ξ_chronic) = {P_acute_shorter:.4f}")

if P_acute_shorter > 0.999:
    print("   *** HYPOTHESIS STRONGLY SUPPORTED ***")
elif P_acute_shorter > 0.99:
    print("   *** HYPOTHESIS SUPPORTED ***")
elif P_acute_shorter > 0.95:
    print("   ** HYPOTHESIS LIKELY **")

print(f"\n🔬 PROTECTION FACTOR EXPONENT:")
print(f"   β_ξ = {β_ξ_samples.mean():.2f} ± {β_ξ_samples.std():.2f}")
print(f"   95% HDI: [{np.percentile(β_ξ_samples, 2.5):.2f}, {np.percentile(β_ξ_samples, 97.5):.2f}]")

naa_acute_pred = posterior['NAA_acute_mean'].values.flatten()
naa_chronic_pred = posterior['NAA_chronic_mean'].values.flatten()
naa_control_pred = posterior['NAA_control_mean'].values.flatten()

print(f"\n📊 PREDICTED NAA CONCENTRATIONS:")
print(f"\n   Acute: {naa_acute_pred.mean():.2f} ± {naa_acute_pred.std():.2f} mM")
print(f"   Observed: {naa_obs_acute.mean():.2f} ± {naa_obs_acute.std():.2f} mM")
print(
    f"   Error: {abs(naa_acute_pred.mean() - naa_obs_acute.mean()):.2f} mM ({abs(naa_acute_pred.mean() - naa_obs_acute.mean()) / naa_obs_acute.mean() * 100:.1f}%)")

print(f"\n   Chronic: {naa_chronic_pred.mean():.2f} ± {naa_chronic_pred.std():.2f} mM")
print(f"   Observed: {naa_obs_chronic.mean():.2f} mM")

print(f"\n   Control: {naa_control_pred.mean():.2f} ± {naa_control_pred.std():.2f} mM")
print(f"   Observed: {naa_obs_control.mean():.2f} mM")

# ============================================================================
# CONVERGENCE
# ============================================================================

print("\n" + "=" * 80)
print("CONVERGENCE DIAGNOSTICS")
print("=" * 80)

rhat = az.rhat(trace)
print(f"\n   Max R-hat: {rhat.max().values:.4f}")
if rhat.max().values < 1.01:
    print("   ✅ Excellent convergence")
elif rhat.max().values < 1.05:
    print("   ✅ Good convergence")

ess = az.ess(trace)
print(f"\n   Min ESS: {ess.min().values:.0f}")
if ess.min().values > 400:
    print("   ✅ Adequate effective sample size")

# ============================================================================
# SAVE
# ============================================================================

print("\n" + "=" * 80)
print("SAVING RESULTS")
print("=" * 80)

output_dir = Path(__file__).parent / 'outputs'
output_dir.mkdir(exist_ok=True)

trace.to_netcdf(output_dir / 'bayesian_v3_6_expanded_trace.nc')
print(f"\n✅ Saved: {output_dir}/bayesian_v3_6_expanded_trace.nc")

summary_df = az.summary(trace, hdi_prob=0.95)
summary_df.to_csv(output_dir / 'bayesian_v3_6_expanded_summary.csv')
print(f"✅ Saved: {output_dir}/bayesian_v3_6_expanded_summary.csv")

results_summary = {
    'Parameter': ['ξ_acute', 'ξ_chronic', 'Δξ', 'β_ξ', 'P(ξ_acute < ξ_chronic)'],
    'Mean': [
        ξ_acute_samples.mean(),
        ξ_chronic_samples.mean(),
        Δξ_samples.mean(),
        β_ξ_samples.mean(),
        P_acute_shorter
    ],
    'SD': [
        ξ_acute_samples.std(),
        ξ_chronic_samples.std(),
        Δξ_samples.std(),
        β_ξ_samples.std(),
        np.nan
    ],
    'HDI_2.5%': [
        np.percentile(ξ_acute_samples, 2.5),
        np.percentile(ξ_chronic_samples, 2.5),
        np.percentile(Δξ_samples, 2.5),
        np.percentile(β_ξ_samples, 2.5),
        np.nan
    ],
    'HDI_97.5%': [
        np.percentile(ξ_acute_samples, 97.5),
        np.percentile(ξ_chronic_samples, 97.5),
        np.percentile(Δξ_samples, 97.5),
        np.percentile(β_ξ_samples, 97.5),
        np.nan
    ]
}

results_df = pd.DataFrame(results_summary)
results_df.to_csv(output_dir / 'bayesian_v3_6_expanded_results.csv', index=False)
print(f"✅ Saved: {output_dir}/bayesian_v3_6_expanded_results.csv")

print("\n" + "=" * 80)
print("✅ ANALYSIS COMPLETE")
print("=" * 80)

print(f"\n🎯 KEY FINDINGS (n={len(naa_obs_acute)} individual patients):")
print(f"   P(ξ_acute < ξ_chronic) = {P_acute_shorter:.4f}")
print(f"   Δξ = {Δξ_samples.mean():.3f} ± {Δξ_samples.std():.3f} nm")
print(f"   β_ξ = {β_ξ_samples.mean():.2f} ± {β_ξ_samples.std():.2f}")
print(
    f"\n   NAA prediction accuracy: {100 - abs(naa_acute_pred.mean() - naa_obs_acute.mean()) / naa_obs_acute.mean() * 100:.1f}%")

print("\n" + "=" * 80)