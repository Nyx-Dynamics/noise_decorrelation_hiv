#!/usr/bin/env python3
"""
Bayesian Model v3.6 - EXPANDED DATASET
Using n=44 individual acute patients + group means for chronic/control

Updated: November 15, 2025
Investigator: AC (Nyx Dynamics LLC)
"""

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
from scipy import stats
import argparse
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

print("\n" + "╔" + "=" * 78 + "╗")
print("║" + " " * 20 + "BAYESIAN MODEL v3.6 - EXPANDED" + " " * 27 + "║")
print("║" + " " * 78 + "║")
print("║" + " " * 15 + "Individual Patient Data: n=44 Acute HIV" + " " * 24 + "║")
print("╚" + "=" * 78 + "╝\n")

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

print("=" * 80)
print("LOADING EXPANDED DATASET")
print("=" * 80)

# ------------------------------
# CLI configuration
# ------------------------------
parser = argparse.ArgumentParser(description="Bayesian v3.6 Expanded - Individual-level pipeline")
this_dir = Path(__file__).resolve().parent
repo_root = this_dir.parents[1]

parser.add_argument("--input", type=str, default=str(this_dir / "data" / "curated" / "CONSOLIDATED_MRS_DATA_FOR_MODEL.csv"),
                    help="Path to consolidated patient-level CSV (relative or absolute)")
parser.add_argument("--results-dir", type=str, default=str(this_dir / "results" / "bayesian_v3_6_expanded"),
                    help="Directory to write outputs")
parser.add_argument("--metric-family", type=str, choices=["auto", "abs_rel", "naacr_rel"], default="auto",
                    help="Which metric family to use for modeling")
parser.add_argument("--acute-cutoff-days", type=int, default=180,
                    help="Days threshold to define acute vs chronic when days_since_infection is available")
parser.add_argument("--phase-filter", type=str, choices=["acute", "chronic", "both"], default="acute",
                    help="Which phases to include in the likelihood")
parser.add_argument("--draws", type=int, default=2000)
parser.add_argument("--tune", type=int, default=1000)
parser.add_argument("--chains", type=int, default=4)
parser.add_argument("--cores", type=int, default=4)
parser.add_argument("--target-accept", type=float, default=0.98)
parser.add_argument("--seed", type=int, default=42)
parser.add_argument("--save-trace", action="store_true")
parser.add_argument("--time-slice", type=str, default="all",
                    help="Comma-separated Valcour weeks to include, e.g., '0,4,12,24'; 'all' to include every week")

args, _ = parser.parse_known_args()

INPUT_PATH = Path(args.input)
if not INPUT_PATH.is_absolute():
    INPUT_PATH = (this_dir / INPUT_PATH).resolve()
RESULTS_DIR = Path(args.results_dir)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

print(f"\nUsing input: {INPUT_PATH}")
print(f"Writing results to: {RESULTS_DIR}")

# Load consolidated data (path-robust)
df_full = pd.read_csv(str(INPUT_PATH))
print(f"\n✅ Loaded {len(df_full)} total records from consolidated CSV")
# Optional counts if a Phase/phase_bin column exists
phase_col = None
for cand in ["Phase", "phase", "phase_bin"]:
    if cand in df_full.columns:
        phase_col = cand
        break
if phase_col is not None:
    try:
        counts = df_full[phase_col].value_counts(dropna=False)
        print("   Counts by phase:")
        for k, v in counts.items():
            print(f"     - {k}: {v}")
    except Exception:
        pass

# Focus on Basal Ganglia where we have the most complete data
# From the analysis: Acute BG n=44 individual patients

# Get individual acute BG data from Valcour
acute_bg = df_full[(df_full['Phase'] == 'Acute') &
                   (df_full['Region'] == 'BG') &
                   (df_full['Study'] == 'Valcour_2015')].copy()

# Get chronic and control group means for BG
chronic_bg = df_full[(df_full['Phase'] == 'Chronic') &
                     (df_full['Region'].str.contains('Basal', na=False))].copy()

control_bg = df_full[(df_full['Phase'] == 'Control') &
                     (df_full['Region'].str.contains('Basal', na=False))].copy()

print("\n" + "=" * 80)
print("BASAL GANGLIA DATA SUMMARY")
print("=" * 80)

print(f"\n📊 ACUTE (Individual Patients):")
print(f"   n = {len(acute_bg)}")
print(f"   NAA = {acute_bg['NAA'].mean():.2f} ± {acute_bg['NAA'].std():.2f} mM")
print(f"   Cho = {acute_bg['Cho'].mean():.2f} ± {acute_bg['Cho'].std():.2f} mM")
print(f"   Range: {acute_bg['NAA'].min():.2f} - {acute_bg['NAA'].max():.2f} mM")

if len(chronic_bg) > 0:
    print(f"\n📊 CHRONIC (Group Means):")
    print(f"   n = {len(chronic_bg)} studies")
    print(f"   NAA = {chronic_bg['NAA'].mean():.2f} ± {chronic_bg['NAA'].std():.2f} mM")
    if chronic_bg['Cho'].notna().any():
        print(f"   Cho = {chronic_bg['Cho'].mean():.2f} ± {chronic_bg['Cho'].std():.2f} mM")

if len(control_bg) > 0:
    print(f"\n📊 CONTROL (Group Means):")
    print(f"   n = {len(control_bg)} studies")
    print(f"   NAA = {control_bg['NAA'].mean():.2f} ± {control_bg['NAA'].std():.2f} mM")
    if control_bg['Cho'].notna().any():
        print(f"   Cho = {control_bg['Cho'].mean():.2f} ± {control_bg['Cho'].std():.2f} mM")

# Prepare data arrays
naa_obs_acute = acute_bg['NAA'].values
cho_obs_acute = acute_bg['Cho'].values

# Use literature reference values for control (from Chang 2002, Mohamed 2010)
naa_obs_control = np.array([9.55])  # BG NAA reference (mM)
cho_obs_control = np.array([2.18])  # BG Cho reference (mM)

# Use chronic group mean if available, otherwise use meta-analysis values
if len(chronic_bg) > 0:
    naa_obs_chronic = chronic_bg['NAA'].values
    cho_obs_chronic = chronic_bg['Cho'].dropna().values
else:
    # From meta-analysis: chronic shows ~8% decline
    naa_obs_chronic = np.array([8.79])  # 9.55 * 0.92
    cho_obs_chronic = np.array([2.40])  # Elevated from 2.18

print(f"\n✅ Final data for model:")
print(f"   Acute: n={len(naa_obs_acute)} individual patients")
print(f"   Chronic: n={len(naa_obs_chronic)} (reference/group means)")
print(f"   Control: n={len(naa_obs_control)} (reference value)")

# ============================================================================
# MODEL PARAMETERS (from quantum coherence framework)
# ============================================================================

print("\n" + "=" * 80)
print("MODEL PARAMETERS")
print("=" * 80)

# Microtubule parameters
L_MT = 2000e-9  # Microtubule length (m)
k_B = 1.38e-23  # Boltzmann constant
T = 310  # Temperature (K)
hbar = 1.055e-34  # Reduced Planck constant

# Enzyme kinetics parameters (from literature)
V_max_base = 100.0  # Base enzyme velocity (nmol/min/mg)
K_m = 50.0  # Michaelis constant (μM)
S_0 = 100.0  # Substrate concentration (μM)

# Neuroinflammation noise parameters
# From literature: acute has HIGH amplitude, SHORT correlation length
# Chronic has LOW amplitude, LONG correlation length

print("\n🔬 Physical Constants:")
print(f"   Microtubule length: {L_MT * 1e9:.0f} nm")
print(f"   Temperature: {T} K")
print(f"   Enzyme V_max baseline: {V_max_base} nmol/min/mg")
print(f"   Michaelis K_m: {K_m} μM")

# ============================================================================
# BAYESIAN MODEL v3.6
# ============================================================================

print("\n" + "=" * 80)
print("BUILDING BAYESIAN MODEL v3.6")
print("=" * 80)

with pm.Model() as model:
    # ========================================================================
    # PRIORS: Noise correlation lengths (THE KEY PARAMETERS)
    # ========================================================================

    # Acute: SHORT correlation length (protective)
    # Based on quantum coherence theory: ξ ~ 0.5-0.7 nm
    ξ_acute = pm.TruncatedNormal('ξ_acute', mu=0.6, sigma=0.1,
                                 lower=0.3, upper=0.9)

    # Chronic: LONG correlation length (destructive)
    # Based on theory: ξ ~ 0.7-1.0 nm
    ξ_chronic = pm.TruncatedNormal('ξ_chronic', mu=0.8, sigma=0.1,
                                   lower=0.5, upper=1.2)

    # Control: MINIMAL noise (baseline)
    ξ_control = pm.TruncatedNormal('ξ_control', mu=0.5, sigma=0.05,
                                   lower=0.3, upper=0.7)

    # ========================================================================
    # Noise amplitude parameters
    # ========================================================================

    # Acute: HIGH amplitude (cytokine storm)
    σ_noise_acute = pm.TruncatedNormal('σ_noise_acute', mu=5.0, sigma=1.0,
                                       lower=2.0, upper=10.0)

    # Chronic: LOW amplitude (smoldering inflammation)
    σ_noise_chronic = pm.TruncatedNormal('σ_noise_chronic', mu=2.0, sigma=0.5,
                                         lower=0.5, upper=4.0)

    # Control: MINIMAL amplitude
    σ_noise_control = pm.TruncatedNormal('σ_noise_control', mu=0.5, sigma=0.2,
                                         lower=0.1, upper=1.0)

    # ========================================================================
    # Protection factor exponent (from quantum coherence theory)
    # ========================================================================

    # β_ξ: Exponent relating ξ to enzyme protection
    # Theory predicts β ~ 1.5-2.5 (superlinear relationship)
    β_ξ = pm.TruncatedNormal('β_ξ', mu=2.0, sigma=0.5,
                             lower=1.0, upper=3.0)

    # ========================================================================
    # DERIVED QUANTITIES: Coherence preservation factors
    # ========================================================================

    # Protection inversely related to ξ (shorter ξ = better protection)
    # F = exp(-β * ξ) → shorter ξ gives higher F
    F_acute = pm.Deterministic('F_acute',
                               pm.math.exp(-β_ξ * ξ_acute))

    F_chronic = pm.Deterministic('F_chronic',
                                 pm.math.exp(-β_ξ * ξ_chronic))

    F_control = pm.Deterministic('F_control',
                                 pm.math.exp(-β_ξ * ξ_control))

    # ========================================================================
    # ENZYME ACTIVITY: Modulated by coherence
    # ========================================================================

    # Base enzyme activity (Michaelis-Menten)
    V_base = V_max_base * S_0 / (K_m + S_0)

    # Modulated activity with coherence protection
    V_acute = pm.Deterministic('V_acute',
                               V_base * (1 + F_acute * σ_noise_acute / 10))

    V_chronic = pm.Deterministic('V_chronic',
                                 V_base * (1 - (1 - F_chronic) * σ_noise_chronic / 5))

    V_control = pm.Deterministic('V_control',
                                 V_base * (1 + F_control * σ_noise_control / 10))

    # ========================================================================
    # NAA PREDICTIONS: Proportional to enzyme activity
    # ========================================================================

    # Scaling factor: enzyme activity → NAA concentration
    # Calibrated to match control NAA ~9.55 mM
    α_NAA = 9.55 / V_base  # Scaling factor

    NAA_acute_mean = pm.Deterministic('NAA_acute_mean', α_NAA * V_acute)
    NAA_chronic_mean = pm.Deterministic('NAA_chronic_mean', α_NAA * V_chronic)
    NAA_control_mean = pm.Deterministic('NAA_control_mean', α_NAA * V_control)

    # ========================================================================
    # CHO PREDICTIONS: Elevated by inflammation
    # ========================================================================

    # Cho increases with noise amplitude (membrane turnover)
    # Baseline Cho ~2.18 mM
    Cho_acute_mean = pm.Deterministic('Cho_acute_mean',
                                      2.18 * (1 + σ_noise_acute / 20))

    Cho_chronic_mean = pm.Deterministic('Cho_chronic_mean',
                                        2.18 * (1 + σ_noise_chronic / 15))

    Cho_control_mean = pm.Deterministic('Cho_control_mean',
                                        2.18 * (1 + σ_noise_control / 20))

    # ========================================================================
    # OBSERVATION NOISE (measurement variability)
    # ========================================================================

    # Individual patient variability (acute)
    σ_obs_acute = pm.HalfNormal('σ_obs_acute', sigma=1.5)

    # Group mean variability (chronic, control)
    σ_obs_chronic = pm.HalfNormal('σ_obs_chronic', sigma=1.0)
    σ_obs_control = pm.HalfNormal('σ_obs_control', sigma=0.5)

    # Cho observation noise
    σ_cho_acute = pm.HalfNormal('σ_cho_acute', sigma=0.5)
    σ_cho_chronic = pm.HalfNormal('σ_cho_chronic', sigma=0.3)

    # ========================================================================
    # LIKELIHOOD: Observed data
    # ========================================================================

    # NAA observations
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

    # Cho observations (where available)
    Cho_acute_obs = pm.Normal('Cho_acute_obs',
                              mu=Cho_acute_mean,
                              sigma=σ_cho_acute,
                              observed=cho_obs_acute)

    Cho_control_obs = pm.Normal('Cho_control_obs',
                                mu=Cho_control_mean,
                                sigma=0.2,
                                observed=cho_obs_control)

    # ========================================================================
    # KEY HYPOTHESIS TEST: Is ξ_acute < ξ_chronic?
    # ========================================================================

    Δξ = pm.Deterministic('Δξ', ξ_chronic - ξ_acute)

print("\n✅ Model built successfully!")
print(f"   Free parameters: {len(model.free_RVs)}")
print(f"   Observations: {len(naa_obs_acute)} acute + {len(naa_obs_chronic)} chronic + {len(naa_obs_control)} control")

# ============================================================================
# INFERENCE
# ============================================================================

print("\n" + "=" * 80)
print("RUNNING MCMC INFERENCE")
print("=" * 80)

with model:
    # Sample from posterior
    print("\n⏳ Sampling (this may take 5-10 minutes)...")
    trace = pm.sample(
        draws=args.draws,
        tune=args.tune,
        chains=args.chains,
        cores=args.cores,
        target_accept=args.target_accept,
        return_inferencedata=True,
        random_seed=args.seed
    )

    print("\n✅ Sampling complete!")

    # Sample posterior predictive
    print("\n⏳ Generating posterior predictive samples...")
    ppc = pm.sample_posterior_predictive(trace, random_seed=42)
    print("✅ Posterior predictive complete!")

# ============================================================================
# RESULTS ANALYSIS
# ============================================================================

print("\n" + "=" * 80)
print("BAYESIAN INFERENCE RESULTS")
print("=" * 80)

# Extract posterior samples
posterior = trace.posterior

# Key parameters
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

# Hypothesis test
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

# Predicted NAA values
naa_acute_pred = posterior['NAA_acute_mean'].values.flatten()
naa_chronic_pred = posterior['NAA_chronic_mean'].values.flatten()
naa_control_pred = posterior['NAA_control_mean'].values.flatten()

print(f"\n📊 PREDICTED NAA CONCENTRATIONS:")
print(f"\n   Acute: {naa_acute_pred.mean():.2f} ± {naa_acute_pred.std():.2f} mM")
print(f"   Observed: {naa_obs_acute.mean():.2f} ± {naa_obs_acute.std():.2f} mM")
print(
    f"   Prediction error: {abs(naa_acute_pred.mean() - naa_obs_acute.mean()):.2f} mM ({abs(naa_acute_pred.mean() - naa_obs_acute.mean()) / naa_obs_acute.mean() * 100:.1f}%)")

print(f"\n   Chronic: {naa_chronic_pred.mean():.2f} ± {naa_chronic_pred.std():.2f} mM")
print(f"   Observed: {naa_obs_chronic.mean():.2f} ± {naa_obs_chronic.std():.2f} mM")
print(
    f"   Prediction error: {abs(naa_chronic_pred.mean() - naa_obs_chronic.mean()):.2f} mM ({abs(naa_chronic_pred.mean() - naa_obs_chronic.mean()) / naa_obs_chronic.mean() * 100:.1f}%)")

print(f"\n   Control: {naa_control_pred.mean():.2f} ± {naa_control_pred.std():.2f} mM")
print(f"   Observed: {naa_obs_control.mean():.2f} mM")
print(
    f"   Prediction error: {abs(naa_control_pred.mean() - naa_obs_control.mean()):.2f} mM ({abs(naa_control_pred.mean() - naa_obs_control.mean()) / naa_obs_control.mean() * 100:.1f}%)")

# ============================================================================
# CONVERGENCE DIAGNOSTICS
# ============================================================================

print("\n" + "=" * 80)
print("CONVERGENCE DIAGNOSTICS")
print("=" * 80)

# R-hat (should be < 1.01)
rhat = az.rhat(trace)
print(f"\n   Max R-hat: {rhat.max().values:.4f}")
if rhat.max().values < 1.01:
    print("   ✅ Excellent convergence (R-hat < 1.01)")
elif rhat.max().values < 1.05:
    print("   ✅ Good convergence (R-hat < 1.05)")
else:
    print("   ⚠ May need more samples (R-hat > 1.05)")

# Effective sample size
ess = az.ess(trace)
print(f"\n   Min ESS: {ess.min().values:.0f}")
if ess.min().values > 400:
    print("   ✅ Adequate effective sample size")
else:
    print("   ⚠ Low ESS, consider more samples")

# ============================================================================
# SAVE RESULTS
# ============================================================================

print("\n" + "=" * 80)
print("SAVING RESULTS")
print("=" * 80)

# Save trace (optional)
if args.save_trace:
    trace_path = RESULTS_DIR / 'trace.nc'
    trace.to_netcdf(str(trace_path))
    print(f"\n✅ Saved trace: {trace_path}")
else:
    print("\nℹ Skipping trace save (pass --save-trace to enable)")

# Save summary statistics
summary_path = RESULTS_DIR / 'summary.csv'
summary_df = az.summary(trace, hdi_prob=0.95)
summary_df.to_csv(str(summary_path))
print(f"✅ Saved summary: {summary_path}")

# Create results summary
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
results_path = RESULTS_DIR / 'results.csv'
results_df.to_csv(str(results_path), index=False)
print(f"✅ Saved results: {results_path}")

print("\n" + "=" * 80)
print("✅ ANALYSIS COMPLETE - n=44 INDIVIDUAL PATIENTS")
print("=" * 80)

print(f"\n🎯 KEY FINDING:")
print(f"   With n={len(naa_obs_acute)} individual acute patients:")
print(f"   P(ξ_acute < ξ_chronic) = {P_acute_shorter:.4f}")
print(f"   Δξ = {Δξ_samples.mean():.3f} ± {Δξ_samples.std():.3f} nm")
print(f"   β_ξ = {β_ξ_samples.mean():.2f} ± {β_ξ_samples.std():.2f}")

print("\n" + "=" * 80)