"""
WAIC Model Comparison - FIXED VERSION 2.0
==========================================

FIXES:
1. Proper WAIC computation with manual combination of log_likelihoods
2. Increased target_accept to 0.98
3. Longer sampling (4000 draws, 3000 tune)
4. Informative priors from bayesian_enzyme_v4 results
5. Better parameterization for Linear ξ model

Author: AC
Date: 2025-11-13
"""

import numpy as np
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)

# ============================================================================
# DATA
# ============================================================================

CONDITIONS = ['healthy', 'acute', 'chronic']
NAA_OBS = np.array([1.105, 1.135, 1.005])
CHO_OBS = np.array([0.225, 0.245, 0.235])

print("\n" + "=" * 80)
print(" WAIC MODEL COMPARISON - FIXED v2.0")
print("=" * 80)
print("\nImprovements:")
print("  ✓ Manual WAIC computation (combines NAA + Cho)")
print("  ✓ target_accept = 0.98 (was 0.92)")
print("  ✓ draws = 4000, tune = 3000 (was 2000/1500)")
print("  ✓ Informative priors from v3.6 results")
print("  ✓ Fixed Linear ξ parameterization")
print()


# ============================================================================
# MODELS
# ============================================================================

def build_full_model_improved():
    """Full model with improved priors."""

    with pm.Model() as model:
        # Tighter ξ priors based on v3.6 results
        xi_mean = pm.TruncatedNormal('xi_mean', mu=0.6e-9, sigma=0.08e-9,
                                     lower=0.5e-9, upper=0.7e-9)
        delta_acute = pm.HalfNormal('delta_acute', sigma=0.05e-9)
        delta_chronic = pm.HalfNormal('delta_chronic', sigma=0.05e-9)

        xi_healthy = pm.Deterministic('xi_healthy', xi_mean)
        xi_acute = pm.Deterministic('xi_acute', xi_mean - delta_acute)
        xi_chronic = pm.Deterministic('xi_chronic', xi_mean + delta_chronic)

        # Informative priors from v3.6
        beta_xi = pm.Normal('beta_xi', mu=1.89, sigma=0.20)
        gamma_coh = pm.Normal('gamma_coh', mu=0.23, sigma=0.10)

        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)

        membrane_acute = pm.TruncatedNormal('membrane_acute', mu=1.1, sigma=0.04,
                                            lower=1.0, upper=1.25)
        membrane_chronic = pm.TruncatedNormal('membrane_chronic', mu=1.03, sigma=0.03,
                                              lower=1.0, upper=1.15)

        # Forward model
        xi_ref = 0.8e-9
        coh_healthy = 0.85
        coh_acute = 0.70
        coh_chronic = 0.75

        Pi_healthy = pm.Deterministic('Pi_healthy', (xi_ref / xi_healthy) ** beta_xi)
        Pi_acute = pm.Deterministic('Pi_acute', (xi_ref / xi_acute) ** beta_xi)
        Pi_chronic = pm.Deterministic('Pi_chronic', (xi_ref / xi_chronic) ** beta_xi)

        Gamma_healthy = coh_healthy ** pm.math.abs(gamma_coh)
        Gamma_acute = coh_acute ** pm.math.abs(gamma_coh)
        Gamma_chronic = coh_chronic ** pm.math.abs(gamma_coh)

        NAA_base = 1.105
        NAA_pred = pm.math.stack([
            NAA_base * Pi_healthy * Gamma_healthy,
            NAA_base * Pi_acute * Gamma_acute * viral_damage_acute,
            NAA_base * Pi_chronic * Gamma_chronic * viral_damage_chronic
        ])

        Cho_base = 0.225
        Cho_pred = pm.math.stack([
            Cho_base,
            Cho_base * membrane_acute,
            Cho_base * membrane_chronic
        ])

        # Likelihood
        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.06)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.03)

        NAA_obs = pm.Normal('NAA_obs', mu=NAA_pred, sigma=sigma_NAA, observed=NAA_OBS)
        Cho_obs = pm.Normal('Cho_obs', mu=Cho_pred, sigma=sigma_Cho, observed=CHO_OBS)

        protection_ratio = pm.Deterministic('protection_ratio', Pi_acute / Pi_chronic)

    return model


def build_no_xi_model_improved():
    """No ξ coupling."""

    with pm.Model() as model:
        gamma_coh = pm.Normal('gamma_coh', mu=0.23, sigma=0.10)

        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)

        membrane_acute = pm.TruncatedNormal('membrane_acute', mu=1.1, sigma=0.04,
                                            lower=1.0, upper=1.25)
        membrane_chronic = pm.TruncatedNormal('membrane_chronic', mu=1.03, sigma=0.03,
                                              lower=1.0, upper=1.15)

        coh_healthy = 0.85
        coh_acute = 0.70
        coh_chronic = 0.75

        Gamma_healthy = coh_healthy ** pm.math.abs(gamma_coh)
        Gamma_acute = coh_acute ** pm.math.abs(gamma_coh)
        Gamma_chronic = coh_chronic ** pm.math.abs(gamma_coh)

        NAA_base = 1.105
        NAA_pred = pm.math.stack([
            NAA_base * Gamma_healthy,
            NAA_base * Gamma_acute * viral_damage_acute,
            NAA_base * Gamma_chronic * viral_damage_chronic
        ])

        Cho_base = 0.225
        Cho_pred = pm.math.stack([
            Cho_base,
            Cho_base * membrane_acute,
            Cho_base * membrane_chronic
        ])

        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.06)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.03)

        NAA_obs = pm.Normal('NAA_obs', mu=NAA_pred, sigma=sigma_NAA, observed=NAA_OBS)
        Cho_obs = pm.Normal('Cho_obs', mu=Cho_pred, sigma=sigma_Cho, observed=CHO_OBS)

    return model


def build_linear_xi_model_improved():
    """Linear ξ with β ≈ 1 via tight prior."""

    with pm.Model() as model:
        xi_mean = pm.TruncatedNormal('xi_mean', mu=0.6e-9, sigma=0.08e-9,
                                     lower=0.5e-9, upper=0.7e-9)
        delta_acute = pm.HalfNormal('delta_acute', sigma=0.05e-9)
        delta_chronic = pm.HalfNormal('delta_chronic', sigma=0.05e-9)

        xi_healthy = pm.Deterministic('xi_healthy', xi_mean)
        xi_acute = pm.Deterministic('xi_acute', xi_mean - delta_acute)
        xi_chronic = pm.Deterministic('xi_chronic', xi_mean + delta_chronic)

        # Use tight prior to enforce β ≈ 1
        beta_xi = pm.Normal('beta_xi', mu=1.0, sigma=0.02)

        gamma_coh = pm.Normal('gamma_coh', mu=0.23, sigma=0.10)

        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)

        membrane_acute = pm.TruncatedNormal('membrane_acute', mu=1.1, sigma=0.04,
                                            lower=1.0, upper=1.25)
        membrane_chronic = pm.TruncatedNormal('membrane_chronic', mu=1.03, sigma=0.03,
                                              lower=1.0, upper=1.15)

        # Forward model
        xi_ref = 0.8e-9
        coh_healthy = 0.85
        coh_acute = 0.70
        coh_chronic = 0.75

        Pi_healthy = pm.Deterministic('Pi_healthy', (xi_ref / xi_healthy) ** beta_xi)
        Pi_acute = pm.Deterministic('Pi_acute', (xi_ref / xi_acute) ** beta_xi)
        Pi_chronic = pm.Deterministic('Pi_chronic', (xi_ref / xi_chronic) ** beta_xi)

        Gamma_healthy = coh_healthy ** pm.math.abs(gamma_coh)
        Gamma_acute = coh_acute ** pm.math.abs(gamma_coh)
        Gamma_chronic = coh_chronic ** pm.math.abs(gamma_coh)

        NAA_base = 1.105
        NAA_pred = pm.math.stack([
            NAA_base * Pi_healthy * Gamma_healthy,
            NAA_base * Pi_acute * Gamma_acute * viral_damage_acute,
            NAA_base * Pi_chronic * Gamma_chronic * viral_damage_chronic
        ])

        Cho_base = 0.225
        Cho_pred = pm.math.stack([
            Cho_base,
            Cho_base * membrane_acute,
            Cho_base * membrane_chronic
        ])

        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.06)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.03)

        NAA_obs = pm.Normal('NAA_obs', mu=NAA_pred, sigma=sigma_NAA, observed=NAA_OBS)
        Cho_obs = pm.Normal('Cho_obs', mu=Cho_pred, sigma=sigma_Cho, observed=CHO_OBS)

        protection_ratio = pm.Deterministic('protection_ratio', Pi_acute / Pi_chronic)

    return model


# ============================================================================
# SAMPLING
# ============================================================================

def sample_model_improved(model, model_name, n_samples=4000):
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

def run_comparison_improved():
    """Run complete improved comparison."""

    models = {
        'Full Model': build_full_model_improved(),
        'No ξ Coupling': build_no_xi_model_improved(),
        'Linear ξ (β≈1)': build_linear_xi_model_improved()
    }

    results = {}
    div_counts = {}

    for name, model in models.items():
        try:
            idata, n_div = sample_model_improved(model, name, n_samples=4000)
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

    # WAIC comparison - manual computation
    print("\n" + "=" * 80)
    print(" WAIC COMPARISON (Manual Computation)")
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
    print(" MODEL RANKING")
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

        print(f"{i}. {name:20s}  Δ: {delta:6.2f}  ({interpretation})")

    # Save results
    output_dir = Path.home() / 'Documents/Github/noise_decorrelation_HIV/quantum/results/model_comparison_v2'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save WAIC comparison
    with open(output_dir / 'waic_comparison.txt', 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(" WAIC MODEL COMPARISON\n")
        f.write("=" * 80 + "\n\n")
        for i, (name, stats) in enumerate(sorted_models, 1):
            delta = stats['waic'] - best_waic
            f.write(f"{i}. {name}\n")
            f.write(f"   WAIC: {stats['waic']:.2f}\n")
            f.write(f"   Δ: {delta:.2f}\n")
            f.write(f"   p_waic: {stats['p_waic']:.2f}\n\n")

    # Save traces
    for name, idata in results.items():
        safe_name = name.replace(' ', '_').replace('(', '').replace(')', '').replace('ξ', 'xi').replace('≈', '')
        idata.to_netcdf(output_dir / f'{safe_name}_trace.nc')

    print(f"\n✓ Results saved to: {output_dir}")

    return results, waic_results, div_counts


if __name__ == '__main__':
    results, waic_results, div_counts = run_comparison_improved()

    if results is not None:
        print("\n" + "=" * 80)
        print(" ANALYSIS COMPLETE!")
        print("=" * 80)

        # Get best model
        sorted_models = sorted(waic_results.items(), key=lambda x: x[1]['waic'])
        best_name = sorted_models[0][0]
        best_waic = sorted_models[0][1]['waic']

        if len(sorted_models) > 1:
            second_waic = sorted_models[1][1]['waic']
            delta = second_waic - best_waic
            print(f"\nBest model: {best_name}")
            print(f"WAIC: {best_waic:.1f}")
            print(f"ΔWAIC from next best: {delta:.1f}")