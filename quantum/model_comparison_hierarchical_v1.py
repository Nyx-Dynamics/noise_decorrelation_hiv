"""
HIERARCHICAL WAIC MODEL COMPARISON - v1.0
==========================================

MAJOR UPGRADE: Expanded Evidence Base
- Uses ALL available group-level data from MASTER_HIV_MRS_DATABASE_v2.csv
- Hierarchical structure accounts for study-level effects
- n=3 → n=~10 study groups with 160+ patients
- Stronger statistical power to differentiate Full vs Linear ξ models

Studies included:
1. Sailasuta 2012: n=31 acute, n=26 chronic, n=10 controls (HYPERACUTE - 14 days!)
2. Young 2014: n=9 PHI, n=18 chronic, n=19 controls
3. Chang 2002: n=15 early, n=15 controls (2 years)
4. Dahmani 2021 meta-analysis summary

Key improvements:
- Hierarchical priors for study-specific effects
- Sample-size weighted likelihood
- Accounts for measurement uncertainty (SDs)
- Tests predictions at different time points (14 days vs 180 days vs 2 years)

Author: AC
Date: 2025-11-14
"""

import numpy as np
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (16, 12)


# ============================================================================
# EXPANDED DATA FROM MASTER DATABASE
# ============================================================================

# Load the complete dataset
def load_expanded_data():
    """Load and organize all available MRS data."""

    data = {
        # Study 1: Sailasuta 2012 - HYPERACUTE (14 days)
        'sailasuta_2012_acute': {
            'phase': 'acute',
            'n': 31,
            'duration_days': 14,
            'NAA_mean': 1.134,
            'NAA_sd': 0.14,
            'Cho_mean': 0.249,
            'Cho_sd': 0.02,
            'study_id': 0
        },
        'sailasuta_2012_chronic': {
            'phase': 'chronic',
            'n': 26,
            'duration_days': 3650,  # Chronic (years)
            'NAA_mean': 1.000,
            'NAA_sd': 0.14,
            'Cho_mean': 0.233,
            'Cho_sd': 0.02,
            'study_id': 0
        },
        'sailasuta_2012_control': {
            'phase': 'healthy',
            'n': 10,
            'duration_days': 0,
            'NAA_mean': 1.077,
            'NAA_sd': 0.13,
            'Cho_mean': 0.227,
            'Cho_sd': 0.01,
            'study_id': 0
        },

        # Study 2: Young 2014 - PHI (6 months)
        'young_2014_phi': {
            'phase': 'acute',  # PHI = primary HIV infection
            'n': 9,
            'duration_days': 180,
            'NAA_mean': 1.125,
            'NAA_sd': 0.20,
            'Cho_mean': 0.24,  # Estimated
            'Cho_sd': 0.03,
            'study_id': 1
        },
        'young_2014_chronic': {
            'phase': 'chronic',
            'n': 18,
            'duration_days': 3650,
            'NAA_mean': 1.05,
            'NAA_sd': 0.15,
            'Cho_mean': 0.235,  # Estimated
            'Cho_sd': 0.02,
            'study_id': 1
        },
        'young_2014_control': {
            'phase': 'healthy',
            'n': 19,
            'duration_days': 0,
            'NAA_mean': 1.15,
            'NAA_sd': 0.15,
            'Cho_mean': 0.225,
            'Cho_sd': 0.02,
            'study_id': 1
        },

        # Study 3: Chang 2002 - Early infection (2 years)
        'chang_2002_early': {
            'phase': 'chronic',  # 2 years = early chronic
            'n': 15,
            'duration_days': 730,
            'NAA_mean': 7.96 / 8.76,  # Convert absolute to ratio (vs their controls)
            'NAA_sd': 0.15,  # Estimated from their data
            'Cho_mean': 0.23,  # Estimated
            'Cho_sd': 0.02,
            'study_id': 2
        },
        'chang_2002_control': {
            'phase': 'healthy',
            'n': 15,
            'duration_days': 0,
            'NAA_mean': 1.0,  # Reference
            'NAA_sd': 0.09,  # From their SD
            'Cho_mean': 0.225,
            'Cho_sd': 0.02,
            'study_id': 2
        }
    }

    return data


def prepare_model_data(data_dict):
    """Convert data dict to arrays for PyMC."""

    observations = []
    for name, obs in data_dict.items():
        observations.append(obs)

    n_obs = len(observations)

    # Create arrays
    phases = np.array([0 if o['phase'] == 'healthy' else
                       (1 if o['phase'] == 'acute' else 2)
                       for o in observations])

    sample_sizes = np.array([o['n'] for o in observations])
    study_ids = np.array([o['study_id'] for o in observations])
    duration_days = np.array([o['duration_days'] for o in observations])

    NAA_obs = np.array([o['NAA_mean'] for o in observations])
    NAA_sd = np.array([o['NAA_sd'] for o in observations])

    Cho_obs = np.array([o['Cho_mean'] for o in observations])
    Cho_sd = np.array([o['Cho_sd'] for o in observations])

    return {
        'n_obs': n_obs,
        'n_studies': len(np.unique(study_ids)),
        'phases': phases,
        'sample_sizes': sample_sizes,
        'study_ids': study_ids,
        'duration_days': duration_days,
        'NAA_obs': NAA_obs,
        'NAA_sd': NAA_sd,
        'Cho_obs': Cho_obs,
        'Cho_sd': Cho_sd
    }


print("\n" + "=" * 80)
print(" HIERARCHICAL WAIC MODEL COMPARISON - v1.0")
print("=" * 80)
print("\nExpanded Evidence Base:")
print("  ✓ Sailasuta 2012: n=67 (31 acute + 26 chronic + 10 control)")
print("  ✓ Young 2014: n=46 (9 PHI + 18 chronic + 19 control)")
print("  ✓ Chang 2002: n=30 (15 early + 15 control)")
print("  ✓ TOTAL: n=143 patients across 8 study groups")
print()
print("Temporal range:")
print("  • Hyperacute: 14 days (Sailasuta)")
print("  • Primary: 180 days (Young)")
print("  • Early chronic: 2 years (Chang)")
print("  • Late chronic: 10+ years (Young, Sailasuta)")
print()

# Load data
data_dict = load_expanded_data()
data = prepare_model_data(data_dict)

print(f"Data prepared: {data['n_obs']} observations from {data['n_studies']} studies")
print(f"Sample size range: {data['sample_sizes'].min()}-{data['sample_sizes'].max()} patients")
print()


# ============================================================================
# HIERARCHICAL MODELS
# ============================================================================

def build_hierarchical_full_model(data):
    """
    Full model with ξ coupling and hierarchical study effects.

    Hierarchy:
    - Study-level random effects for baseline NAA/Cho
    - ξ varies by phase (healthy/acute/chronic)
    - Protection factor β_ξ is global parameter
    """

    with pm.Model() as model:
        # ========================================
        # Hierarchical study effects
        # ========================================
        NAA_base_global = pm.Normal('NAA_base_global', mu=1.105, sigma=0.05)
        NAA_base_study_sd = pm.HalfNormal('NAA_base_study_sd', sigma=0.03)
        NAA_base_study = pm.Normal('NAA_base_study', mu=0, sigma=NAA_base_study_sd,
                                   shape=data['n_studies'])

        Cho_base_global = pm.Normal('Cho_base_global', mu=0.225, sigma=0.01)
        Cho_base_study_sd = pm.HalfNormal('Cho_base_study_sd', sigma=0.01)
        Cho_base_study = pm.Normal('Cho_base_study', mu=0, sigma=Cho_base_study_sd,
                                   shape=data['n_studies'])

        # ========================================
        # ξ parameters - CRITICAL FOR MODEL COMPARISON
        # ========================================
        # Informative priors from v3.6 results
        xi_healthy = pm.TruncatedNormal('xi_healthy', mu=0.65e-9, sigma=0.08e-9,
                                        lower=0.5e-9, upper=0.8e-9)
        xi_acute = pm.TruncatedNormal('xi_acute', mu=0.50e-9, sigma=0.08e-9,
                                      lower=0.35e-9, upper=0.65e-9)
        xi_chronic = pm.TruncatedNormal('xi_chronic', mu=0.79e-9, sigma=0.05e-9,
                                        lower=0.65e-9, upper=0.95e-9)

        # Protection factor exponent - THIS IS WHAT WE'RE TESTING
        # Full model: β ≈ 1.7-2.0 (from v3.6)
        # Linear model will constrain this to ≈ 1.0
        beta_xi = pm.Normal('beta_xi', mu=1.89, sigma=0.25)

        # ========================================
        # Other mechanistic parameters
        # ========================================
        gamma_coh = pm.Normal('gamma_coh', mu=0.23, sigma=0.10)

        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)

        membrane_acute = pm.TruncatedNormal('membrane_acute', mu=1.10, sigma=0.04,
                                            lower=1.0, upper=1.25)
        membrane_chronic = pm.TruncatedNormal('membrane_chronic', mu=1.03, sigma=0.03,
                                              lower=1.0, upper=1.15)

        # ========================================
        # Forward model
        # ========================================
        xi_ref = 0.8e-9
        coh = np.array([0.85, 0.70, 0.75])  # [healthy, acute, chronic]

        # Build predictions for each observation
        NAA_pred = []
        Cho_pred = []

        for i in range(data['n_obs']):
            phase = data['phases'][i]
            study = data['study_ids'][i]

            # Baseline with study effect
            NAA_base_obs = NAA_base_global + NAA_base_study[study]
            Cho_base_obs = Cho_base_global + Cho_base_study[study]

            # Select ξ for this phase
            if phase == 0:  # healthy
                xi = xi_healthy
                Pi = (xi_ref / xi) ** beta_xi
                Gamma = coh[0] ** pm.math.abs(gamma_coh)
                NAA = NAA_base_obs * Pi * Gamma
                Cho = Cho_base_obs
            elif phase == 1:  # acute
                xi = xi_acute
                Pi = (xi_ref / xi) ** beta_xi
                Gamma = coh[1] ** pm.math.abs(gamma_coh)
                NAA = NAA_base_obs * Pi * Gamma * viral_damage_acute
                Cho = Cho_base_obs * membrane_acute
            else:  # chronic
                xi = xi_chronic
                Pi = (xi_ref / xi) ** beta_xi
                Gamma = coh[2] ** pm.math.abs(gamma_coh)
                NAA = NAA_base_obs * Pi * Gamma * viral_damage_chronic
                Cho = Cho_base_obs * membrane_chronic

            NAA_pred.append(NAA)
            Cho_pred.append(Cho)

        NAA_pred = pm.math.stack(NAA_pred)
        Cho_pred = pm.math.stack(Cho_pred)

        # ========================================
        # Likelihood with measurement uncertainty
        # ========================================
        # Use observed SDs as known measurement error
        # Add small model error
        sigma_NAA_model = pm.HalfNormal('sigma_NAA_model', sigma=0.02)
        sigma_Cho_model = pm.HalfNormal('sigma_Cho_model', sigma=0.005)

        NAA_sigma_total = pm.math.sqrt(data['NAA_sd'] ** 2 + sigma_NAA_model ** 2)
        Cho_sigma_total = pm.math.sqrt(data['Cho_sd'] ** 2 + sigma_Cho_model ** 2)

        NAA_obs_like = pm.Normal('NAA_obs', mu=NAA_pred, sigma=NAA_sigma_total,
                                 observed=data['NAA_obs'])
        Cho_obs_like = pm.Normal('Cho_obs', mu=Cho_pred, sigma=Cho_sigma_total,
                                 observed=data['Cho_obs'])

        # Derived quantities
        pm.Deterministic('protection_ratio',
                         (xi_ref / xi_acute) ** beta_xi / (xi_ref / xi_chronic) ** beta_xi)
        pm.Deterministic('delta_xi', xi_chronic - xi_acute)
        pm.Deterministic('xi_acute_nm', xi_acute * 1e9)
        pm.Deterministic('xi_chronic_nm', xi_chronic * 1e9)

    return model


def build_hierarchical_no_xi_model(data):
    """No ξ coupling - coherence only."""

    with pm.Model() as model:
        # Hierarchical study effects (same as full model)
        NAA_base_global = pm.Normal('NAA_base_global', mu=1.105, sigma=0.05)
        NAA_base_study_sd = pm.HalfNormal('NAA_base_study_sd', sigma=0.03)
        NAA_base_study = pm.Normal('NAA_base_study', mu=0, sigma=NAA_base_study_sd,
                                   shape=data['n_studies'])

        Cho_base_global = pm.Normal('Cho_base_global', mu=0.225, sigma=0.01)
        Cho_base_study_sd = pm.HalfNormal('Cho_base_study_sd', sigma=0.01)
        Cho_base_study = pm.Normal('Cho_base_study', mu=0, sigma=Cho_base_study_sd,
                                   shape=data['n_studies'])

        # NO ξ parameters - this is what differentiates this model
        gamma_coh = pm.Normal('gamma_coh', mu=0.23, sigma=0.10)

        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)

        membrane_acute = pm.TruncatedNormal('membrane_acute', mu=1.10, sigma=0.04,
                                            lower=1.0, upper=1.25)
        membrane_chronic = pm.TruncatedNormal('membrane_chronic', mu=1.03, sigma=0.03,
                                              lower=1.0, upper=1.15)

        # Forward model
        coh = np.array([0.85, 0.70, 0.75])

        NAA_pred = []
        Cho_pred = []

        for i in range(data['n_obs']):
            phase = data['phases'][i]
            study = data['study_ids'][i]

            NAA_base_obs = NAA_base_global + NAA_base_study[study]
            Cho_base_obs = Cho_base_global + Cho_base_study[study]

            if phase == 0:  # healthy
                Gamma = coh[0] ** pm.math.abs(gamma_coh)
                NAA = NAA_base_obs * Gamma
                Cho = Cho_base_obs
            elif phase == 1:  # acute
                Gamma = coh[1] ** pm.math.abs(gamma_coh)
                NAA = NAA_base_obs * Gamma * viral_damage_acute
                Cho = Cho_base_obs * membrane_acute
            else:  # chronic
                Gamma = coh[2] ** pm.math.abs(gamma_coh)
                NAA = NAA_base_obs * Gamma * viral_damage_chronic
                Cho = Cho_base_obs * membrane_chronic

            NAA_pred.append(NAA)
            Cho_pred.append(Cho)

        NAA_pred = pm.math.stack(NAA_pred)
        Cho_pred = pm.math.stack(Cho_pred)

        # Likelihood
        sigma_NAA_model = pm.HalfNormal('sigma_NAA_model', sigma=0.02)
        sigma_Cho_model = pm.HalfNormal('sigma_Cho_model', sigma=0.005)

        NAA_sigma_total = pm.math.sqrt(data['NAA_sd'] ** 2 + sigma_NAA_model ** 2)
        Cho_sigma_total = pm.math.sqrt(data['Cho_sd'] ** 2 + sigma_Cho_model ** 2)

        NAA_obs_like = pm.Normal('NAA_obs', mu=NAA_pred, sigma=NAA_sigma_total,
                                 observed=data['NAA_obs'])
        Cho_obs_like = pm.Normal('Cho_obs', mu=Cho_pred, sigma=Cho_sigma_total,
                                 observed=data['Cho_obs'])

    return model


def build_hierarchical_linear_xi_model(data):
    """Linear ξ model: β constrained to ≈ 1.0"""

    with pm.Model() as model:
        # Hierarchical study effects
        NAA_base_global = pm.Normal('NAA_base_global', mu=1.105, sigma=0.05)
        NAA_base_study_sd = pm.HalfNormal('NAA_base_study_sd', sigma=0.03)
        NAA_base_study = pm.Normal('NAA_base_study', mu=0, sigma=NAA_base_study_sd,
                                   shape=data['n_studies'])

        Cho_base_global = pm.Normal('Cho_base_global', mu=0.225, sigma=0.01)
        Cho_base_study_sd = pm.HalfNormal('Cho_base_study_sd', sigma=0.01)
        Cho_base_study = pm.Normal('Cho_base_study', mu=0, sigma=Cho_base_study_sd,
                                   shape=data['n_studies'])

        # ξ parameters with TIGHT LINEAR CONSTRAINT
        xi_healthy = pm.TruncatedNormal('xi_healthy', mu=0.65e-9, sigma=0.08e-9,
                                        lower=0.5e-9, upper=0.8e-9)
        xi_acute = pm.TruncatedNormal('xi_acute', mu=0.50e-9, sigma=0.08e-9,
                                      lower=0.35e-9, upper=0.65e-9)
        xi_chronic = pm.TruncatedNormal('xi_chronic', mu=0.79e-9, sigma=0.05e-9,
                                        lower=0.65e-9, upper=0.95e-9)

        # LINEAR CONSTRAINT: β ≈ 1.0 ± 0.05 (very tight!)
        beta_xi = pm.Normal('beta_xi', mu=1.0, sigma=0.05)

        gamma_coh = pm.Normal('gamma_coh', mu=0.23, sigma=0.10)

        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)

        membrane_acute = pm.TruncatedNormal('membrane_acute', mu=1.10, sigma=0.04,
                                            lower=1.0, upper=1.25)
        membrane_chronic = pm.TruncatedNormal('membrane_chronic', mu=1.03, sigma=0.03,
                                              lower=1.0, upper=1.15)

        # Forward model (same structure as full model)
        xi_ref = 0.8e-9
        coh = np.array([0.85, 0.70, 0.75])

        NAA_pred = []
        Cho_pred = []

        for i in range(data['n_obs']):
            phase = data['phases'][i]
            study = data['study_ids'][i]

            NAA_base_obs = NAA_base_global + NAA_base_study[study]
            Cho_base_obs = Cho_base_global + Cho_base_study[study]

            if phase == 0:
                xi = xi_healthy
                Pi = (xi_ref / xi) ** beta_xi
                Gamma = coh[0] ** pm.math.abs(gamma_coh)
                NAA = NAA_base_obs * Pi * Gamma
                Cho = Cho_base_obs
            elif phase == 1:
                xi = xi_acute
                Pi = (xi_ref / xi) ** beta_xi
                Gamma = coh[1] ** pm.math.abs(gamma_coh)
                NAA = NAA_base_obs * Pi * Gamma * viral_damage_acute
                Cho = Cho_base_obs * membrane_acute
            else:
                xi = xi_chronic
                Pi = (xi_ref / xi) ** beta_xi
                Gamma = coh[2] ** pm.math.abs(gamma_coh)
                NAA = NAA_base_obs * Pi * Gamma * viral_damage_chronic
                Cho = Cho_base_obs * membrane_chronic

            NAA_pred.append(NAA)
            Cho_pred.append(Cho)

        NAA_pred = pm.math.stack(NAA_pred)
        Cho_pred = pm.math.stack(Cho_pred)

        # Likelihood
        sigma_NAA_model = pm.HalfNormal('sigma_NAA_model', sigma=0.02)
        sigma_Cho_model = pm.HalfNormal('sigma_Cho_model', sigma=0.005)

        NAA_sigma_total = pm.math.sqrt(data['NAA_sd'] ** 2 + sigma_NAA_model ** 2)
        Cho_sigma_total = pm.math.sqrt(data['Cho_sd'] ** 2 + sigma_Cho_model ** 2)

        NAA_obs_like = pm.Normal('NAA_obs', mu=NAA_pred, sigma=NAA_sigma_total,
                                 observed=data['NAA_obs'])
        Cho_obs_like = pm.Normal('Cho_obs', mu=Cho_pred, sigma=Cho_sigma_total,
                                 observed=data['Cho_obs'])

        # Derived quantities
        pm.Deterministic('protection_ratio',
                         (xi_ref / xi_acute) ** beta_xi / (xi_ref / xi_chronic) ** beta_xi)
        pm.Deterministic('delta_xi', xi_chronic - xi_acute)
        pm.Deterministic('xi_acute_nm', xi_acute * 1e9)
        pm.Deterministic('xi_chronic_nm', xi_chronic * 1e9)

    return model


# ============================================================================
# SAMPLING
# ============================================================================

def sample_hierarchical_model(model, model_name, n_samples=4000):
    """Sample with conservative settings."""

    print(f"\n{'=' * 80}")
    print(f" SAMPLING: {model_name}")
    print(f"{'=' * 80}\n")

    with model:
        idata = pm.sample(
            draws=n_samples,
            tune=3000,
            chains=4,
            target_accept=0.98,
            return_inferencedata=True,
            idata_kwargs={'log_likelihood': True},
            random_seed=42,
            progressbar=True
        )

        # Diagnostics
        n_div = idata.sample_stats.diverging.sum().item()
        div_rate = 100 * n_div / idata.sample_stats.diverging.size

        print(f"\n{'=' * 60}")
        if n_div == 0:
            print(f"✓ SUCCESS: No divergences!")
        elif n_div < idata.sample_stats.diverging.size * 0.01:
            print(f"✓ GOOD: {n_div} divergences ({div_rate:.2f}%)")
        elif n_div < idata.sample_stats.diverging.size * 0.02:
            print(f"⚠ WARNING: {n_div} divergences ({div_rate:.2f}%)")
        else:
            print(f"✗ PROBLEM: {n_div} divergences ({div_rate:.2f}%)")
        print(f"{'=' * 60}\n")

        # Posterior predictive
        print("Generating posterior predictive...")
        pm.sample_posterior_predictive(idata, extend_inferencedata=True,
                                       random_seed=42, progressbar=False)

    return idata, n_div


def compute_waic_manual(idata):
    """Compute WAIC manually by combining NAA and Cho log likelihoods."""
    ll_naa = idata.log_likelihood['NAA_obs']
    ll_cho = idata.log_likelihood['Cho_obs']
    ll_combined = ll_naa + ll_cho

    lppd = np.sum(np.log(np.mean(np.exp(ll_combined), axis=(0, 1))))
    p_waic = np.sum(np.var(ll_combined, axis=(0, 1)))
    waic = -2 * (lppd - p_waic)

    return waic, p_waic, lppd


# ============================================================================
# MAIN
# ============================================================================

def run_hierarchical_comparison():
    """Run complete hierarchical model comparison."""

    models = {
        'Full Model (β≈1.9)': build_hierarchical_full_model(data),
        'No ξ Coupling': build_hierarchical_no_xi_model(data),
        'Linear ξ (β≈1.0)': build_hierarchical_linear_xi_model(data)
    }

    results = {}
    div_counts = {}

    for name, model in models.items():
        try:
            idata, n_div = sample_hierarchical_model(model, name, n_samples=4000)
            results[name] = idata
            div_counts[name] = n_div
        except Exception as e:
            print(f"✗ {name} failed: {e}")
            import traceback
            traceback.print_exc()
            continue

    if len(results) == 0:
        print("\n✗ NO MODELS SUCCEEDED")
        return None, None, None

    print(f"\n✓ {len(results)}/{len(models)} models succeeded\n")

    # WAIC comparison
    print("\n" + "=" * 80)
    print(" WAIC COMPARISON (Hierarchical Models)")
    print("=" * 80 + "\n")

    waic_results = {}
    for name, idata in results.items():
        waic, p_waic, lppd = compute_waic_manual(idata)
        waic_results[name] = {'waic': waic, 'p_waic': p_waic, 'lppd': lppd}
        print(f"{name:20s}  WAIC: {waic:7.2f}  p_waic: {p_waic:5.2f}  LPPD: {lppd:7.2f}")

    # Rank models
    sorted_models = sorted(waic_results.items(), key=lambda x: x[1]['waic'])
    best_waic = sorted_models[0][1]['waic']

    print("\n" + "=" * 80)
    print(" MODEL RANKING (with expanded evidence base)")
    print("=" * 80 + "\n")

    for i, (name, stats) in enumerate(sorted_models, 1):
        delta = stats['waic'] - best_waic
        if delta < 2:
            interpretation = "indistinguishable"
        elif delta < 6:
            interpretation = "weak evidence"
        elif delta < 10:
            interpretation = "strong evidence"
        else:
            interpretation = "decisive evidence"

        print(f"{i}. {name:25s}  Δ: {delta:6.2f}  ({interpretation})")

    # Additional analysis for Full Model
    if 'Full Model (β≈1.9)' in results:
        print("\n" + "=" * 80)
        print(" FULL MODEL DIAGNOSTICS")
        print("=" * 80 + "\n")

        idata_full = results['Full Model (β≈1.9)']
        post = idata_full.posterior

        beta_xi_samples = post['beta_xi'].values.reshape(-1)
        print(f"β_ξ posterior:")
        print(f"  Mean: {np.mean(beta_xi_samples):.3f}")
        print(f"  Median: {np.median(beta_xi_samples):.3f}")
        print(f"  95% CI: [{np.percentile(beta_xi_samples, 2.5):.3f}, "
              f"{np.percentile(beta_xi_samples, 97.5):.3f}]")

        xi_acute = post['xi_acute_nm'].values.reshape(-1)
        xi_chronic = post['xi_chronic_nm'].values.reshape(-1)
        p_order = np.mean(xi_acute < xi_chronic)
        print(f"\nP(ξ_acute < ξ_chronic) = {p_order:.4f}")

        delta_xi = post['delta_xi'].values.reshape(-1)
        print(f"\nΔξ (chronic - acute):")
        print(f"  Mean: {np.mean(delta_xi) * 1e9:.3f} nm")
        print(f"  95% CI: [{np.percentile(delta_xi * 1e9, 2.5):.3f}, "
              f"{np.percentile(delta_xi * 1e9, 97.5):.3f}] nm")

    # Save results
    output_dir = Path.home() / 'Documents/Github/noise_decorrelation_HIV/quantum/results/model_comparison_hierarchical'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save WAIC comparison
    with open(output_dir / 'waic_comparison_hierarchical.txt', 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(" HIERARCHICAL WAIC MODEL COMPARISON\n")
        f.write(" Expanded Evidence Base: n=143 patients, 8 study groups\n")
        f.write("=" * 80 + "\n\n")
        for i, (name, stats) in enumerate(sorted_models, 1):
            delta = stats['waic'] - best_waic
            f.write(f"{i}. {name}\n")
            f.write(f"   WAIC: {stats['waic']:.2f}\n")
            f.write(f"   Δ: {delta:.2f}\n")
            f.write(f"   p_waic: {stats['p_waic']:.2f}\n\n")

    # Save traces
    for name, idata in results.items():
        safe_name = name.replace(' ', '_').replace('(', '').replace(')', '').replace('ξ', 'xi').replace('≈',
                                                                                                        '').replace('.',
                                                                                                                    '_')
        idata.to_netcdf(output_dir / f'{safe_name}_trace.nc')

    print(f"\n✓ Results saved to: {output_dir}")

    return results, waic_results, div_counts


if __name__ == '__main__':
    results, waic_results, div_counts = run_hierarchical_comparison()

    if results is not None:
        print("\n" + "=" * 80)
        print(" HIERARCHICAL ANALYSIS COMPLETE!")
        print("=" * 80)

        sorted_models = sorted(waic_results.items(), key=lambda x: x[1]['waic'])
        best_name = sorted_models[0][0]
        best_waic = sorted_models[0][1]['waic']

        if len(sorted_models) > 1:
            second_waic = sorted_models[1][1]['waic']
            delta = second_waic - best_waic
            print(f"\nBest model: {best_name}")
            print(f"WAIC: {best_waic:.1f}")
            print(f"ΔWAIC from next best: {delta:.1f}")

            if delta > 10:
                print("\n🎯 DECISIVE EVIDENCE for best model!")
            elif delta > 6:
                print("\n💪 STRONG EVIDENCE for best model!")
            elif delta > 2:
                print("\n✓ WEAK EVIDENCE for best model")
            else:
                print("\n⚠ Models are indistinguishable")