"""
WAIC Model Comparison - CLEAN VERSION
======================================

Fixes for Mac + PyMC compatibility:
1. Non-centered parameterization for ξ
2. Proper log_likelihood computation  
3. Mac file paths
4. PyMC 5.x compatible sampling (no step_kwargs issues)

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

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

# ============================================================================
# OBSERVED DATA
# ============================================================================

CONDITIONS = ['healthy', 'acute', 'chronic']
NAA_OBS = np.array([1.105, 1.135, 1.005])
CHO_OBS = np.array([0.225, 0.245, 0.235])

print("\n" + "="*80)
print(" WAIC MODEL COMPARISON - CLEAN VERSION")
print("="*80)
print("\nKey improvements:")
print("  ✓ Non-centered ξ parameterization")
print("  ✓ No step_kwargs issues")
print("  ✓ Mac-compatible paths")
print("  ✓ Proper log_likelihood")
print()

# ============================================================================
# MODELS
# ============================================================================

def build_full_model():
    """Full model with fixed parameterization."""
    
    with pm.Model() as model:
        # Non-centered ξ
        xi_mean = pm.TruncatedNormal('xi_mean', mu=0.6e-9, sigma=0.15e-9,
                                      lower=0.4e-9, upper=0.9e-9)
        delta_acute = pm.HalfNormal('delta_acute', sigma=0.1e-9)
        delta_chronic = pm.HalfNormal('delta_chronic', sigma=0.1e-9)
        
        xi_healthy = pm.Deterministic('xi_healthy', xi_mean)
        xi_acute = pm.Deterministic('xi_acute', xi_mean - delta_acute)
        xi_chronic = pm.Deterministic('xi_chronic', xi_mean + delta_chronic)
        
        # Other parameters
        beta_xi = pm.TruncatedNormal('beta_xi', mu=1.8, sigma=0.5,
                                      lower=0.5, upper=3.5)
        gamma_coh = pm.TruncatedNormal('gamma_coh', mu=1.0, sigma=0.3,
                                        lower=0.3, upper=2.0)
        
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
        
        Gamma_healthy = coh_healthy ** gamma_coh
        Gamma_acute = coh_acute ** gamma_coh
        Gamma_chronic = coh_chronic ** gamma_coh
        
        NAA_base = 1.105
        NAA_healthy = pm.Deterministic('NAA_healthy',
            NAA_base * Pi_healthy * Gamma_healthy
        )
        NAA_acute = pm.Deterministic('NAA_acute',
            NAA_base * Pi_acute * Gamma_acute * viral_damage_acute
        )
        NAA_chronic = pm.Deterministic('NAA_chronic',
            NAA_base * Pi_chronic * Gamma_chronic * viral_damage_chronic
        )
        
        NAA_pred = pm.math.stack([NAA_healthy, NAA_acute, NAA_chronic])
        
        Cho_base = 0.225
        Cho_pred = pm.math.stack([Cho_base, Cho_base * membrane_acute,
                                  Cho_base * membrane_chronic])
        
        # Likelihood
        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.06)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.03)
        
        NAA_obs = pm.Normal('NAA_obs', mu=NAA_pred, sigma=sigma_NAA, observed=NAA_OBS)
        Cho_obs = pm.Normal('Cho_obs', mu=Cho_pred, sigma=sigma_Cho, observed=CHO_OBS)
        
        protection_ratio = pm.Deterministic('protection_ratio', Pi_acute / Pi_chronic)
        
    return model


def build_no_xi_model():
    """No ξ coupling."""
    
    with pm.Model() as model:
        gamma_coh = pm.TruncatedNormal('gamma_coh', mu=1.0, sigma=0.3,
                                        lower=0.3, upper=2.0)
        
        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)
        
        membrane_acute = pm.TruncatedNormal('membrane_acute', mu=1.1, sigma=0.04,
                                            lower=1.0, upper=1.25)
        membrane_chronic = pm.TruncatedNormal('membrane_chronic', mu=1.03, sigma=0.03,
                                               lower=1.0, upper=1.15)
        
        coh_healthy = 0.85
        coh_acute = 0.70
        coh_chronic = 0.75
        
        Gamma_healthy = coh_healthy ** gamma_coh
        Gamma_acute = coh_acute ** gamma_coh
        Gamma_chronic = coh_chronic ** gamma_coh
        
        NAA_base = 1.105
        NAA_healthy = NAA_base * Gamma_healthy
        NAA_acute = NAA_base * Gamma_acute * viral_damage_acute
        NAA_chronic = NAA_base * Gamma_chronic * viral_damage_chronic
        
        NAA_pred = pm.math.stack([NAA_healthy, NAA_acute, NAA_chronic])
        
        Cho_base = 0.225
        Cho_pred = pm.math.stack([Cho_base, Cho_base * membrane_acute,
                                  Cho_base * membrane_chronic])
        
        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.06)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.03)
        
        NAA_obs = pm.Normal('NAA_obs', mu=NAA_pred, sigma=sigma_NAA, observed=NAA_OBS)
        Cho_obs = pm.Normal('Cho_obs', mu=Cho_pred, sigma=sigma_Cho, observed=CHO_OBS)
        
    return model


def build_linear_xi_model():
    """Linear ξ coupling (β=1)."""
    
    with pm.Model() as model:
        xi_mean = pm.TruncatedNormal('xi_mean', mu=0.6e-9, sigma=0.15e-9,
                                      lower=0.4e-9, upper=0.9e-9)
        delta_acute = pm.HalfNormal('delta_acute', sigma=0.1e-9)
        delta_chronic = pm.HalfNormal('delta_chronic', sigma=0.1e-9)
        
        xi_healthy = pm.Deterministic('xi_healthy', xi_mean)
        xi_acute = pm.Deterministic('xi_acute', xi_mean - delta_acute)
        xi_chronic = pm.Deterministic('xi_chronic', xi_mean + delta_chronic)
        
        beta_xi = 1.0  # FIXED
        
        gamma_coh = pm.TruncatedNormal('gamma_coh', mu=1.0, sigma=0.3,
                                        lower=0.3, upper=2.0)
        
        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)
        
        membrane_acute = pm.TruncatedNormal('membrane_acute', mu=1.1, sigma=0.04,
                                            lower=1.0, upper=1.25)
        membrane_chronic = pm.TruncatedNormal('membrane_chronic', mu=1.03, sigma=0.03,
                                               lower=1.0, upper=1.15)
        
        xi_ref = 0.8e-9
        coh_healthy = 0.85
        coh_acute = 0.70
        coh_chronic = 0.75
        
        Pi_healthy = pm.Deterministic('Pi_healthy', (xi_ref / xi_healthy) ** beta_xi)
        Pi_acute = pm.Deterministic('Pi_acute', (xi_ref / xi_acute) ** beta_xi)
        Pi_chronic = pm.Deterministic('Pi_chronic', (xi_ref / xi_chronic) ** beta_xi)
        
        Gamma_healthy = coh_healthy ** gamma_coh
        Gamma_acute = coh_acute ** gamma_coh
        Gamma_chronic = coh_chronic ** gamma_coh
        
        NAA_base = 1.105
        NAA_healthy = NAA_base * Pi_healthy * Gamma_healthy
        NAA_acute = NAA_base * Pi_acute * Gamma_acute * viral_damage_acute
        NAA_chronic = NAA_base * Pi_chronic * Gamma_chronic * viral_damage_chronic
        
        NAA_pred = pm.math.stack([NAA_healthy, NAA_acute, NAA_chronic])
        
        Cho_base = 0.225
        Cho_pred = pm.math.stack([Cho_base, Cho_base * membrane_acute,
                                  Cho_base * membrane_chronic])
        
        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.06)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.03)
        
        NAA_obs = pm.Normal('NAA_obs', mu=NAA_pred, sigma=sigma_NAA, observed=NAA_OBS)
        Cho_obs = pm.Normal('Cho_obs', mu=Cho_pred, sigma=sigma_Cho, observed=CHO_OBS)
        
        protection_ratio = pm.Deterministic('protection_ratio', Pi_acute / Pi_chronic)
        
    return model


# ============================================================================
# SAMPLING
# ============================================================================

def sample_model(model, model_name, n_samples=2000):
    """Sample with clean API (no step_kwargs issues)."""
    
    print(f"\n{'='*80}")
    print(f" SAMPLING: {model_name}")
    print(f"{'='*80}\n")
    
    with model:
        # Clean sampling - no step_kwargs!
        idata = pm.sample(
            draws=n_samples,
            tune=1500,
            chains=4,
            target_accept=0.92,
            return_inferencedata=True,
            idata_kwargs={'log_likelihood': True},
            random_seed=42,
            progressbar=True
        )
        
        # Check divergences
        n_div = idata.sample_stats.diverging.sum().values
        div_rate = 100 * n_div / idata.sample_stats.diverging.size
        
        print(f"\n{'='*60}")
        if n_div == 0:
            print(f"✓ SUCCESS: No divergences!")
        elif n_div < 50:
            print(f"⚠ {n_div} divergences ({div_rate:.1f}%) - acceptable")
        else:
            print(f"❌ {n_div} divergences ({div_rate:.1f}%) - concerning!")
        print(f"{'='*60}\n")
        
        # Posterior predictive
        if n_div < 500:
            print("Generating posterior predictive...")
            pm.sample_posterior_predictive(idata, extend_inferencedata=True,
                                           random_seed=42, progressbar=False)
    
    return idata, n_div


# ============================================================================
# MAIN
# ============================================================================

def run_comparison():
    """Run complete comparison."""
    
    models = {
        'Full Model': build_full_model(),
        'No ξ Coupling': build_no_xi_model(),
        'Linear ξ (β=1)': build_linear_xi_model()
    }
    
    results = {}
    div_counts = {}
    
    for name, model in models.items():
        try:
            idata, n_div = sample_model(model, name, n_samples=2000)
            results[name] = idata
            div_counts[name] = n_div
        except Exception as e:
            print(f"❌ {name} failed: {e}")
            continue
    
    if len(results) == 0:
        print("\n❌ NO MODELS SUCCEEDED")
        return None, None, None, None
    
    print(f"\n✓ {len(results)}/{len(models)} models succeeded\n")
    
    # WAIC comparison
    print("\n" + "="*80)
    print(" WAIC COMPARISON")
    print("="*80 + "\n")
    
    try:
        comparison = az.compare(results, ic='waic', scale='deviance')
        print(comparison)
    except Exception as e:
        print(f"❌ WAIC failed: {e}")
        comparison = None
    
    # Save results
    output_dir = Path.home() / 'Documents/Github/noise_decorrelation_HIV/quantum/quantum/results/model_comparison_clean'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if comparison is not None:
        comparison.to_csv(output_dir / 'waic_comparison.csv')
    
    for name, idata in results.items():
        safe_name = name.replace(' ', '_').replace('(', '').replace(')', '').replace('ξ', 'xi')
        idata.to_netcdf(output_dir / f'{safe_name}_trace.nc')
    
    print(f"\n✓ Results saved to: {output_dir}")
    
    return results, comparison, div_counts, output_dir


if __name__ == '__main__':
    results, comparison, div_counts, output_dir = run_comparison()
    
    if results is not None:
        print("\n" + "="*80)
        print(" ANALYSIS COMPLETE!")
        print("="*80)
        if comparison is not None:
            print(f"\nBest model: {comparison.index[0]}")
            print(f"WAIC: {comparison['waic'].iloc[0]:.1f}")
