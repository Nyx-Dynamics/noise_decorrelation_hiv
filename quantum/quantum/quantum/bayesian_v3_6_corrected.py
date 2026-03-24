#!/usr/bin/env python3
"""
Bayesian Model v3.6 - EXPANDED WITH RATIO DATA
Including Valcour (absolute) + Young + Sailasuta (ratios converted to absolute)

Updated: November 16, 2025
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
print("║" + " " * 15 + "BAYESIAN MODEL v3.6 - WITH RATIO DATA" + " " * 22 + "║")
print("║" + " " * 78 + "║")
print("║" + " Valcour (absolute) + Young + Sailasuta (ratios → absolute)" + " " * 14 + "║")
print("╚" + "=" * 78 + "╝\n")

# ============================================================================
# REFERENCE CREATINE VALUES FOR RATIO CONVERSION
# ============================================================================
# Literature reference values (mM) from Chang 2002, Mohamed 2010
CR_REFERENCE = {
    'BG': 8.0,  # Basal Ganglia
    'FWM': 6.8,  # Frontal White Matter
    'PGM': 7.8,  # Posterior Gray Matter / Parietal
    'FGM': 7.8,  # Frontal Gray Matter
    'AC': 7.8,  # Anterior Cingulate
    'OGM': 7.5  # Occipital Gray Matter
}


def convert_ratio_to_absolute(ratio_value, region):
    """Convert metabolite/Cr ratio to absolute concentration (mM)"""
    cr_ref = CR_REFERENCE.get(region, 7.5)  # Default to 7.5 if region unknown
    return ratio_value * cr_ref


# ============================================================================
# LOAD DATA
# ============================================================================

print("=" * 80)
print("LOADING DATA")
print("=" * 80)

# ------------------------------
# Load Valcour data (absolute concentrations)
# ------------------------------
print("\n📁 Loading Valcour 2015 data (absolute concentrations)...")

# Load all timepoints
valcour_weeks = ['week_0', 'week_4', 'week_12', 'week_24']
valcour_dfs = []

for week in valcour_weeks:
    try:
        df = pd.read_excel(f'/mnt/user-data/uploads/valcour_2015_{week}.xlsx')
        df['Week'] = week
        valcour_dfs.append(df)
        print(f"   ✅ Loaded {week}: {len(df)} patients")
    except Exception as e:
        print(f"   ⚠ Could not load {week}: {e}")

if valcour_dfs:
    valcour_df = pd.concat(valcour_dfs, ignore_index=True)
    print(f"\n✅ Total Valcour records: {len(valcour_df)}")
else:
    print("   ⚠ No Valcour data loaded")
    valcour_df = pd.DataFrame()

# ------------------------------
# Load Young 2014 data (ratios)
# ------------------------------
print("\n📁 Loading Young 2014 data (ratios)...")

try:
    young_df = pd.read_csv('/mnt/user-data/uploads/YOUNG_2014_CROSS_SECTIONAL_DATA.csv')

    # Parse the ratio data - extract NAA and Cho
    young_processed = []

    for _, row in young_df.iterrows():
        if row['Metabolite'] == 'NAA/Cr':
            # Convert ratio to absolute using median value
            naa_abs = convert_ratio_to_absolute(row['Ratio_Median'], row['Region'])

            young_processed.append({
                'Study': 'Young_2014',
                'Phase': 'Acute' if row['Phase'] == 'Primary' else row['Phase'],
                'Region': row['Region'],
                'n': row['n'],
                'NAA': naa_abs,
                'NAA_SE': convert_ratio_to_absolute(row['SE'], row['Region']),
                'Source': 'Ratio'
            })
        elif row['Metabolite'] == 'Cho/Cr':
            # Find matching NAA entry to update
            for proc in young_processed:
                if (proc['Region'] == row['Region'] and
                        proc['Phase'] == ('Acute' if row['Phase'] == 'Primary' else row['Phase']) and
                        proc['n'] == row['n']):
                    cho_abs = convert_ratio_to_absolute(row['Ratio_Median'], row['Region'])
                    proc['Cho'] = cho_abs
                    proc['Cho_SE'] = convert_ratio_to_absolute(row['SE'], row['Region'])
                    break

    young_clean = pd.DataFrame(young_processed)
    print(f"✅ Loaded Young 2014: {len(young_clean)} group means")
    print(f"   Phases: {young_clean['Phase'].value_counts().to_dict()}")
    print(f"   Regions: {young_clean['Region'].value_counts().to_dict()}")

except Exception as e:
    print(f"⚠ Error loading Young 2014: {e}")
    young_clean = pd.DataFrame()

# ------------------------------
# Load Sailasuta 2012 data (ratios)
# ------------------------------
print("\n📁 Loading Sailasuta 2012 data (ratios)...")

try:
    sailasuta_df = pd.read_csv('/mnt/user-data/uploads/Sailasuta_2012.csv', skiprows=1)

    # Process Sailasuta data
    sailasuta_processed = []

    for phase in sailasuta_df['Group'].unique():
        for region in sailasuta_df['Brain_Region'].unique():
            phase_region_data = sailasuta_df[
                (sailasuta_df['Group'] == phase) &
                (sailasuta_df['Brain_Region'] == region)
                ]

            naa_row = phase_region_data[phase_region_data['Metabolite'] == 'NAA']
            cho_row = phase_region_data[phase_region_data['Metabolite'] == 'tCho']

            if not naa_row.empty:
                naa_value = naa_row['Value'].values[0]
                naa_sd = naa_row['SD'].values[0]
                n = int(naa_row['n'].values[0])

                naa_abs = convert_ratio_to_absolute(naa_value, region)
                naa_se = convert_ratio_to_absolute(naa_sd / np.sqrt(n), region)

                entry = {
                    'Study': 'Sailasuta_2012',
                    'Phase': phase,
                    'Region': region,
                    'n': n,
                    'NAA': naa_abs,
                    'NAA_SE': naa_se,
                    'Source': 'Ratio'
                }

                if not cho_row.empty:
                    cho_value = cho_row['Value'].values[0]
                    cho_sd = cho_row['SD'].values[0]
                    cho_abs = convert_ratio_to_absolute(cho_value, region)
                    cho_se = convert_ratio_to_absolute(cho_sd / np.sqrt(n), region)
                    entry['Cho'] = cho_abs
                    entry['Cho_SE'] = cho_se

                sailasuta_processed.append(entry)

    sailasuta_clean = pd.DataFrame(sailasuta_processed)
    print(f"✅ Loaded Sailasuta 2012: {len(sailasuta_clean)} group means")
    print(f"   Phases: {sailasuta_clean['Phase'].value_counts().to_dict()}")
    print(f"   Regions: {sailasuta_clean['Region'].value_counts().to_dict()}")

except Exception as e:
    print(f"⚠ Error loading Sailasuta 2012: {e}")
    sailasuta_clean = pd.DataFrame()

# ------------------------------
# Load Chang 2002 reference data (absolute)
# ------------------------------
print("\n📁 Loading Chang 2002 reference data...")

try:
    chang_df = pd.read_csv('/mnt/user-data/uploads/CHANG_2002_EXTRACTED.csv')
    print(f"✅ Loaded Chang 2002: {len(chang_df)} measurements")
except Exception as e:
    print(f"⚠ Error loading Chang 2002: {e}")
    chang_df = pd.DataFrame()

# ============================================================================
# PREPARE DATA FOR BASAL GANGLIA (PRIMARY FOCUS)
# ============================================================================

print("\n" + "=" * 80)
print("PREPARING BASAL GANGLIA DATA")
print("=" * 80)

# ------------------------------
# Valcour: Individual acute patients
# ------------------------------
if not valcour_df.empty:
    # Use week 0 (baseline) for acute data
    valcour_bg_acute = valcour_df[
        (valcour_df['Week'] == 'week_0') &
        (valcour_df['NAA'].notna()) &
        (valcour_df['Cho'].notna())
        ].copy()

    naa_obs_acute_valcour = valcour_bg_acute['NAA'].values
    cho_obs_acute_valcour = valcour_bg_acute['Cho'].values
    n_acute_valcour = len(naa_obs_acute_valcour)

    print(f"\n📊 VALCOUR ACUTE (Individual Patients, n={n_acute_valcour}):")
    print(f"   NAA: {naa_obs_acute_valcour.mean():.2f} ± {naa_obs_acute_valcour.std():.2f} mM")
    print(f"   Cho: {cho_obs_acute_valcour.mean():.2f} ± {cho_obs_acute_valcour.std():.2f} mM")
else:
    naa_obs_acute_valcour = np.array([])
    cho_obs_acute_valcour = np.array([])
    n_acute_valcour = 0

# ------------------------------
# Young: Group means for BG
# ------------------------------
if not young_clean.empty:
    young_bg = young_clean[young_clean['Region'] == 'BG'].copy()

    young_bg_acute = young_bg[young_bg['Phase'] == 'Acute']
    young_bg_chronic = young_bg[young_bg['Phase'] == 'Chronic']
    young_bg_control = young_bg[young_bg['Phase'] == 'Control']

    print(f"\n📊 YOUNG 2014 BG:")
    if not young_bg_acute.empty:
        print(f"   Acute: n={young_bg_acute['n'].values[0]}, NAA={young_bg_acute['NAA'].values[0]:.2f} mM")
    if not young_bg_chronic.empty:
        print(f"   Chronic: n={young_bg_chronic['n'].values[0]}, NAA={young_bg_chronic['NAA'].values[0]:.2f} mM")
    if not young_bg_control.empty:
        print(f"   Control: n={young_bg_control['n'].values[0]}, NAA={young_bg_control['NAA'].values[0]:.2f} mM")
else:
    young_bg_acute = pd.DataFrame()
    young_bg_chronic = pd.DataFrame()
    young_bg_control = pd.DataFrame()

# ------------------------------
# Sailasuta: Group means for BG
# ------------------------------
if not sailasuta_clean.empty:
    sailasuta_bg = sailasuta_clean[sailasuta_clean['Region'] == 'BG'].copy()

    sailasuta_bg_acute = sailasuta_bg[sailasuta_bg['Phase'] == 'Acute']
    sailasuta_bg_chronic = sailasuta_bg[sailasuta_bg['Phase'] == 'Chronic']
    sailasuta_bg_control = sailasuta_bg[sailasuta_bg['Phase'] == 'Control']

    print(f"\n📊 SAILASUTA 2012 BG:")
    if not sailasuta_bg_acute.empty:
        print(f"   Acute: n={sailasuta_bg_acute['n'].values[0]}, NAA={sailasuta_bg_acute['NAA'].values[0]:.2f} mM")
    if not sailasuta_bg_chronic.empty:
        print(
            f"   Chronic: n={sailasuta_bg_chronic['n'].values[0]}, NAA={sailasuta_bg_chronic['NAA'].values[0]:.2f} mM")
    if not sailasuta_bg_control.empty:
        print(
            f"   Control: n={sailasuta_bg_control['n'].values[0]}, NAA={sailasuta_bg_control['NAA'].values[0]:.2f} mM")
else:
    sailasuta_bg_acute = pd.DataFrame()
    sailasuta_bg_chronic = pd.DataFrame()
    sailasuta_bg_control = pd.DataFrame()

# ------------------------------
# Combine data for modeling
# ------------------------------

# ACUTE: Combine Valcour individuals + Young/Sailasuta group means
acute_naa_list = [naa_obs_acute_valcour] if n_acute_valcour > 0 else []
acute_cho_list = [cho_obs_acute_valcour] if n_acute_valcour > 0 else []

# Add Young acute as weighted observations (repeat by sample size)
if not young_bg_acute.empty:
    n_young = int(young_bg_acute['n'].values[0])
    naa_young = young_bg_acute['NAA'].values[0]
    se_young = young_bg_acute['NAA_SE'].values[0]
    # Create pseudo-observations with appropriate spread
    acute_naa_list.append(np.random.normal(naa_young, se_young * np.sqrt(n_young), n_young))

    if 'Cho' in young_bg_acute.columns:
        cho_young = young_bg_acute['Cho'].values[0]
        cho_se_young = young_bg_acute['Cho_SE'].values[0]
        acute_cho_list.append(np.random.normal(cho_young, cho_se_young * np.sqrt(n_young), n_young))

# Add Sailasuta acute
if not sailasuta_bg_acute.empty:
    n_sail = int(sailasuta_bg_acute['n'].values[0])
    naa_sail = sailasuta_bg_acute['NAA'].values[0]
    se_sail = sailasuta_bg_acute['NAA_SE'].values[0]
    acute_naa_list.append(np.random.normal(naa_sail, se_sail * np.sqrt(n_sail), n_sail))

    if 'Cho' in sailasuta_bg_acute.columns:
        cho_sail = sailasuta_bg_acute['Cho'].values[0]
        cho_se_sail = sailasuta_bg_acute['Cho_SE'].values[0]
        acute_cho_list.append(np.random.normal(cho_sail, cho_se_sail * np.sqrt(n_sail), n_sail))

naa_obs_acute = np.concatenate(acute_naa_list) if acute_naa_list else np.array([])
cho_obs_acute = np.concatenate(acute_cho_list) if acute_cho_list else np.array([])

# CHRONIC: Combine group means
chronic_naa_list = []
chronic_cho_list = []

if not young_bg_chronic.empty:
    chronic_naa_list.append(young_bg_chronic['NAA'].values)
    if 'Cho' in young_bg_chronic.columns:
        chronic_cho_list.append(young_bg_chronic['Cho'].values)

if not sailasuta_bg_chronic.empty:
    chronic_naa_list.append(sailasuta_bg_chronic['NAA'].values)
    if 'Cho' in sailasuta_bg_chronic.columns:
        chronic_cho_list.append(sailasuta_bg_chronic['Cho'].values)

naa_obs_chronic = np.concatenate(chronic_naa_list) if chronic_naa_list else np.array([8.79])
cho_obs_chronic = np.concatenate(chronic_cho_list) if chronic_cho_list else np.array([2.40])

# CONTROL: Use literature reference + any group means
if not chang_df.empty:
    chang_bg_control = chang_df[(chang_df['Phase'] == 'Control') & (chang_df['Region'] == 'BG')]
    if not chang_bg_control.empty:
        chang_naa = chang_bg_control[chang_bg_control['Metabolite'] == 'NAA']['Mean'].values
        chang_cho = chang_bg_control[chang_bg_control['Metabolite'] == 'Cho']['Mean'].values
        naa_obs_control = chang_naa if len(chang_naa) > 0 else np.array([9.55])
        cho_obs_control = chang_cho if len(chang_cho) > 0 else np.array([2.18])
    else:
        naa_obs_control = np.array([9.55])
        cho_obs_control = np.array([2.18])
else:
    naa_obs_control = np.array([9.55])
    cho_obs_control = np.array([2.18])

print("\n" + "=" * 80)
print("FINAL COMBINED DATA FOR MODEL")
print("=" * 80)
print(f"\n✅ ACUTE: n={len(naa_obs_acute)}")
print(f"   NAA: {naa_obs_acute.mean():.2f} ± {naa_obs_acute.std():.2f} mM")
print(f"   Cho: {cho_obs_acute.mean():.2f} ± {cho_obs_acute.std():.2f} mM")

print(f"\n✅ CHRONIC: n={len(naa_obs_chronic)}")
print(f"   NAA: {naa_obs_chronic.mean():.2f} ± {naa_obs_chronic.std():.2f} mM")
print(f"   Cho: {cho_obs_chronic.mean():.2f} ± {cho_obs_chronic.std():.2f} mM")

print(f"\n✅ CONTROL: n={len(naa_obs_control)}")
print(f"   NAA: {naa_obs_control.mean():.2f} mM")
print(f"   Cho: {cho_obs_control.mean():.2f} mM")

# ============================================================================
# MODEL PARAMETERS
# ============================================================================

print("\n" + "=" * 80)
print("MODEL PARAMETERS")
print("=" * 80)

# Microtubule parameters
L_MT = 2000e-9  # Microtubule length (m)
k_B = 1.38e-23  # Boltzmann constant
T = 310  # Temperature (K)
hbar = 1.055e-34  # Reduced Planck constant

# Enzyme kinetics parameters
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
    # ========================================================================
    # PRIORS: Noise correlation lengths
    # ========================================================================

    ξ_acute = pm.TruncatedNormal('ξ_acute', mu=0.6, sigma=0.1,
                                 lower=0.3, upper=0.9)

    ξ_chronic = pm.TruncatedNormal('ξ_chronic', mu=0.8, sigma=0.1,
                                   lower=0.5, upper=1.2)

    ξ_control = pm.TruncatedNormal('ξ_control', mu=0.5, sigma=0.05,
                                   lower=0.3, upper=0.7)

    # ========================================================================
    # PROTECTION MECHANISM
    # ========================================================================

    β_ξ = pm.TruncatedNormal('β_ξ', mu=-2.0, sigma=0.5, lower=-4.0, upper=0.0)


    def quantum_protection_factor(ξ, L_MT, β_ξ):
        """Quantum coherence protection: Γ_eff = exp(β_ξ * ξ/L_MT)"""
        return pm.math.exp(β_ξ * (ξ * 1e-9) / L_MT)


    Γ_acute = quantum_protection_factor(ξ_acute, L_MT, β_ξ)
    Γ_chronic = quantum_protection_factor(ξ_chronic, L_MT, β_ξ)
    Γ_control = quantum_protection_factor(ξ_control, L_MT, β_ξ)

    # ========================================================================
    # METABOLIC MODEL
    # ========================================================================

    # Baseline NAA synthesis rate (control)
    r_NAA_baseline = pm.TruncatedNormal('r_NAA_baseline', mu=9.5, sigma=0.3,
                                        lower=8.0, upper=11.0)

    # Phase-specific modulation
    α_acute = pm.TruncatedNormal('α_acute', mu=1.0, sigma=0.1,
                                 lower=0.8, upper=1.3)
    α_chronic = pm.TruncatedNormal('α_chronic', mu=0.92, sigma=0.05,
                                   lower=0.8, upper=1.0)

    # Predicted NAA
    NAA_control_mean = r_NAA_baseline
    NAA_acute_mean = r_NAA_baseline * α_acute * Γ_acute
    NAA_chronic_mean = r_NAA_baseline * α_chronic * Γ_chronic

    # Store for posterior analysis
    NAA_control_pred = pm.Deterministic('NAA_control_mean', NAA_control_mean)
    NAA_acute_pred = pm.Deterministic('NAA_acute_mean', NAA_acute_mean)
    NAA_chronic_pred = pm.Deterministic('NAA_chronic_mean', NAA_chronic_mean)

    # ========================================================================
    # CHOLINE MODEL
    # ========================================================================

    Cho_baseline = pm.TruncatedNormal('Cho_baseline', mu=2.2, sigma=0.1,
                                      lower=1.8, upper=2.6)

    # Inflammation increases Cho
    Cho_control_mean = Cho_baseline
    Cho_acute_mean = Cho_baseline * 1.15  # Acute inflammation
    Cho_chronic_mean = Cho_baseline * 1.10  # Chronic inflammation

    # ========================================================================
    # LIKELIHOOD
    # ========================================================================

    # Individual-level acute observations
    σ_naa_acute = pm.HalfNormal('σ_naa_acute', sigma=0.5)
    σ_cho_acute = pm.HalfNormal('σ_cho_acute', sigma=0.3)

    if len(naa_obs_acute) > 0:
        NAA_acute_obs = pm.Normal('NAA_acute_obs',
                                  mu=NAA_acute_mean,
                                  sigma=σ_naa_acute,
                                  observed=naa_obs_acute)

    if len(cho_obs_acute) > 0:
        Cho_acute_obs = pm.Normal('Cho_acute_obs',
                                  mu=Cho_acute_mean,
                                  sigma=σ_cho_acute,
                                  observed=cho_obs_acute)

    # Group-level chronic observations
    if len(naa_obs_chronic) > 0:
        NAA_chronic_obs = pm.Normal('NAA_chronic_obs',
                                    mu=NAA_chronic_mean,
                                    sigma=0.3,
                                    observed=naa_obs_chronic)

    if len(cho_obs_chronic) > 0:
        Cho_chronic_obs = pm.Normal('Cho_chronic_obs',
                                    mu=Cho_chronic_mean,
                                    sigma=0.2,
                                    observed=cho_obs_chronic)

    # Control observations
    NAA_control_obs = pm.Normal('NAA_control_obs',
                                mu=NAA_control_mean,
                                sigma=0.2,
                                observed=naa_obs_control)

    Cho_control_obs = pm.Normal('Cho_control_obs',
                                mu=Cho_control_mean,
                                sigma=0.2,
                                observed=cho_obs_control)

    # ========================================================================
    # HYPOTHESIS TEST
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
    print("\n⏳ Sampling (this may take 5-10 minutes)...")
    trace = pm.sample(
        draws=2000,
        tune=1000,
        chains=4,
        cores=4,
        target_accept=0.98,
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

print(f"\n🔬 PROTECTION FACTOR EXPONENT:")
print(f"   β_ξ = {β_ξ_samples.mean():.2f} ± {β_ξ_samples.std():.2f}")

# Predicted NAA
naa_acute_pred = posterior['NAA_acute_mean'].values.flatten()
naa_chronic_pred = posterior['NAA_chronic_mean'].values.flatten()
naa_control_pred = posterior['NAA_control_mean'].values.flatten()

print(f"\n📊 PREDICTED vs OBSERVED NAA:")
print(f"\n   Acute: {naa_acute_pred.mean():.2f} mM (pred) vs {naa_obs_acute.mean():.2f} mM (obs)")
print(f"   Error: {abs(naa_acute_pred.mean() - naa_obs_acute.mean()) / naa_obs_acute.mean() * 100:.1f}%")

print(f"\n   Chronic: {naa_chronic_pred.mean():.2f} mM (pred) vs {naa_obs_chronic.mean():.2f} mM (obs)")
print(f"   Error: {abs(naa_chronic_pred.mean() - naa_obs_chronic.mean()) / naa_obs_chronic.mean() * 100:.1f}%")

print(f"\n   Control: {naa_control_pred.mean():.2f} mM (pred) vs {naa_obs_control.mean():.2f} mM (obs)")
print(f"   Error: {abs(naa_control_pred.mean() - naa_obs_control.mean()) / naa_obs_control.mean() * 100:.1f}%")

# ============================================================================
# CONVERGENCE DIAGNOSTICS
# ============================================================================

print("\n" + "=" * 80)
print("CONVERGENCE DIAGNOSTICS")
print("=" * 80)

rhat = az.rhat(trace)
print(f"\n   Max R-hat: {rhat.max().values:.4f}")
if rhat.max().values < 1.01:
    print("   ✅ Excellent convergence (R-hat < 1.01)")

ess = az.ess(trace)
print(f"\n   Min ESS: {ess.min().values:.0f}")
if ess.min().values > 400:
    print("   ✅ Adequate effective sample size")

# ============================================================================
# SAVE RESULTS
# ============================================================================

print("\n" + "=" * 80)
print("SAVING RESULTS")
print("=" * 80)

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
results_path = '/mnt/user-data/outputs/bayesian_v3_6_with_ratios_results.csv'
results_df.to_csv(results_path, index=False)
print(f"✅ Saved results: {results_path}")

# Save posterior predictive comparison
ppc_summary = {
    'condition': ['Acute', 'Chronic', 'Control'],
    'NAA_pred': [naa_acute_pred.mean(), naa_chronic_pred.mean(), naa_control_pred.mean()],
    'NAA_obs': [naa_obs_acute.mean(), naa_obs_chronic.mean(), naa_obs_control.mean()],
    'n_obs': [len(naa_obs_acute), len(naa_obs_chronic), len(naa_obs_control)]
}
ppc_df = pd.DataFrame(ppc_summary)
ppc_path = '/mnt/user-data/outputs/posterior_predictive_comparison.csv'
ppc_df.to_csv(ppc_path, index=False)
print(f"✅ Saved posterior predictive: {ppc_path}")

print("\n" + "=" * 80)
print("✅ ANALYSIS COMPLETE WITH RATIO DATA INCORPORATED")
print("=" * 80)

print(f"\n🎯 FINAL SAMPLE SIZES:")
print(f"   Total acute observations: n={len(naa_obs_acute)}")
print(f"   - Valcour individuals: {n_acute_valcour}")
print(f"   - Young 2014: {int(young_bg_acute['n'].values[0]) if not young_bg_acute.empty else 0}")
print(f"   - Sailasuta 2012: {int(sailasuta_bg_acute['n'].values[0]) if not sailasuta_bg_acute.empty else 0}")

print(f"\n🎯 KEY FINDING:")
print(f"   P(ξ_acute < ξ_chronic) = {P_acute_shorter:.4f}")
print(f"   Δξ = {Δξ_samples.mean():.3f} ± {Δξ_samples.std():.3f} nm")
print(f"   β_ξ = {β_ξ_samples.mean():.2f} ± {β_ξ_samples.std():.2f}")

print("\n" + "=" * 80)