"""
Working Model Comparison with WAIC/LOO
======================================

Fixes from failed run:
1. Non-centered parameterization for ALL correlated parameters
2. Log-space sampling for positive parameters
3. Proper log_likelihood computation
4. Mac-compatible file paths
5. Conservative sampling settings
6. Informative priors that actually constrain the space

Author: AC
Date: 2025-11-12
"""

import numpy as np
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# OBSERVED DATA
# ============================================================================

CONDITIONS = ['healthy', 'acute', 'chronic']
NAA_OBS = np.array([1.105, 1.135, 1.005])
CHO_OBS = np.array([0.225, 0.245, 0.235])

print("\n" + "="*80)
print(" FIXED MODEL COMPARISON - WORKING VERSION")
print("="*80)
print("\nKey Improvements:")
print("  ✓ Non-centered parameterization (hierarchical)")
print("  ✓ Log-space sampling for bounded parameters")
print("  ✓ Proper log_likelihood computation")
print("  ✓ Mac-compatible file paths")
print("  ✓ Conservative sampling (fewer divergences)")
print()

# ============================================================================
# MODEL 1: FULL MODEL (FIXED)
# ============================================================================

def build_full_model_fixed():
    """
    Full model with FIXED parameterization.
    
    Key changes:
    - Hierarchical ξ (mean + offsets)
    - Informative priors from v3.6/v4 posteriors
    - No truncated normals (use HalfNormal + transforms)
    """
    
    with pm.Model() as model:
        
        # ===== HIERARCHICAL ξ PARAMETERIZATION =====
        # Sample mean ξ, then offsets (reduces correlation)
        
        # Mean ξ (around 0.6 nm based on v3.6)
        xi_mean_raw = pm.Normal('xi_mean_raw', mu=0, sigma=1)
        xi_mean = pm.Deterministic('xi_mean', 0.6e-9 + 0.15e-9 * xi_mean_raw)
        
        # Acute offset (negative, so acute < mean)
        xi_acute_offset = pm.HalfNormal('xi_acute_offset', sigma=0.1e-9)
        xi_acute = pm.Deterministic('xi_acute', xi_mean - xi_acute_offset)
        
        # Chronic offset (positive, so chronic > mean)
        xi_chronic_offset = pm.HalfNormal('xi_chronic_offset', sigma=0.1e-9)
        xi_chronic = pm.Deterministic('xi_chronic', xi_mean + xi_chronic_offset)
        
        # Healthy ξ (use mean as reference)
        xi_healthy = pm.Deterministic('xi_healthy', xi_mean)
        
        # ===== COUPLING EXPONENTS =====
        # Based on v3.6: β_ξ ≈ 1.7, γ_coh ≈ 1.1
        
        beta_xi = pm.Normal('beta_xi', mu=1.7, sigma=0.5)  # Allow negative? No, add constraint
        gamma_coh = pm.Normal('gamma_coh', mu=1.0, sigma=0.3)
        
        # ===== VIRAL DAMAGE =====
        # These were fine, keep as-is
        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)
        
        # ===== MEMBRANE EFFECTS =====
        # Use log-normal (always positive, less constrained)
        membrane_acute = pm.LogNormal('membrane_acute', mu=np.log(1.1), sigma=0.05)
        membrane_chronic = pm.LogNormal('membrane_chronic', mu=np.log(1.03), sigma=0.04)
        
        # ===== FORWARD MODEL =====
        xi_ref = 0.8e-9
        coh_healthy = 0.85
        coh_acute = 0.70
        coh_chronic = 0.75
        
        # Protection factors (add small constant to avoid division by zero)
        Pi_healthy = pm.Deterministic('Pi_healthy', 
            (xi_ref / (xi_healthy + 1e-12)) ** pm.math.abs(beta_xi)
        )
        Pi_acute = pm.Deterministic('Pi_acute',
            (xi_ref / (xi_acute + 1e-12)) ** pm.math.abs(beta_xi)
        )
        Pi_chronic = pm.Deterministic('Pi_chronic',
            (xi_ref / (xi_chronic + 1e-12)) ** pm.math.abs(beta_xi)
        )
        
        # Coherence enhancement (use abs to avoid negative powers)
        Gamma_healthy = coh_healthy ** pm.math.abs(gamma_coh)
        Gamma_acute = coh_acute ** pm.math.abs(gamma_coh)
        Gamma_chronic = coh_chronic ** pm.math.abs(gamma_coh)
        
        # NAA predictions
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
        
        # Cho predictions
        Cho_base = 0.225
        Cho_healthy = Cho_base
        Cho_acute = Cho_base * membrane_acute
        Cho_chronic = Cho_base * membrane_chronic
        
        Cho_pred = pm.math.stack([Cho_healthy, Cho_acute, Cho_chronic])
        
        # ===== LIKELIHOOD =====
        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.08)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.04)
        
        # CRITICAL: observed= must be set for log_likelihood computation
        NAA_likelihood = pm.Normal('NAA_obs', mu=NAA_pred, sigma=sigma_NAA, 
                                    observed=NAA_OBS)
        Cho_likelihood = pm.Normal('Cho_obs', mu=Cho_pred, sigma=sigma_Cho,
                                    observed=CHO_OBS)
        
        # Derived quantities
        protection_ratio = pm.Deterministic('protection_ratio', 
                                           Pi_acute / (Pi_chronic + 1e-12))
        
    return model


# ============================================================================
# MODEL 2: NO ξ (SIMPLIFIED)
# ============================================================================

def build_no_xi_model():
    """Model without ξ dependence - coherence only."""
    
    with pm.Model() as model:
        gamma_coh = pm.Normal('gamma_coh', mu=1.0, sigma=0.3)
        
        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)
        
        membrane_acute = pm.LogNormal('membrane_acute', mu=np.log(1.1), sigma=0.05)
        membrane_chronic = pm.LogNormal('membrane_chronic', mu=np.log(1.03), sigma=0.04)
        
        # Forward model - NO ξ protection
        coh_healthy = 0.85
        coh_acute = 0.70
        coh_chronic = 0.75
        
        Gamma_healthy = coh_healthy ** pm.math.abs(gamma_coh)
        Gamma_acute = coh_acute ** pm.math.abs(gamma_coh)
        Gamma_chronic = coh_chronic ** pm.math.abs(gamma_coh)
        
        NAA_base = 1.105
        
        NAA_healthy = NAA_base * Gamma_healthy
        NAA_acute = NAA_base * Gamma_acute * viral_damage_acute
        NAA_chronic = NAA_base * Gamma_chronic * viral_damage_chronic
        
        NAA_pred = pm.math.stack([NAA_healthy, NAA_acute, NAA_chronic])
        
        Cho_base = 0.225
        Cho_pred = pm.math.stack([Cho_base, Cho_base * membrane_acute, 
                                  Cho_base * membrane_chronic])
        
        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.08)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.04)
        
        NAA_obs = pm.Normal('NAA_obs', mu=NAA_pred, sigma=sigma_NAA, observed=NAA_OBS)
        Cho_obs = pm.Normal('Cho_obs', mu=Cho_pred, sigma=sigma_Cho, observed=CHO_OBS)
        
    return model


# ============================================================================
# MODEL 3: LINEAR ξ (β = 1)
# ============================================================================

def build_linear_xi_model():
    """Model with linear ξ dependence (β fixed to 1)."""
    
    with pm.Model() as model:
        # Same hierarchical ξ
        xi_mean_raw = pm.Normal('xi_mean_raw', mu=0, sigma=1)
        xi_mean = pm.Deterministic('xi_mean', 0.6e-9 + 0.15e-9 * xi_mean_raw)
        
        xi_acute_offset = pm.HalfNormal('xi_acute_offset', sigma=0.1e-9)
        xi_acute = pm.Deterministic('xi_acute', xi_mean - xi_acute_offset)
        
        xi_chronic_offset = pm.HalfNormal('xi_chronic_offset', sigma=0.1e-9)
        xi_chronic = pm.Deterministic('xi_chronic', xi_mean + xi_chronic_offset)
        
        xi_healthy = pm.Deterministic('xi_healthy', xi_mean)
        
        # FIXED β = 1
        beta_xi = 1.0
        
        gamma_coh = pm.Normal('gamma_coh', mu=1.0, sigma=0.3)
        
        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)
        
        membrane_acute = pm.LogNormal('membrane_acute', mu=np.log(1.1), sigma=0.05)
        membrane_chronic = pm.LogNormal('membrane_chronic', mu=np.log(1.03), sigma=0.04)
        
        # Forward model with LINEAR ξ
        xi_ref = 0.8e-9
        coh_healthy = 0.85
        coh_acute = 0.70
        coh_chronic = 0.75
        
        Pi_healthy = pm.Deterministic('Pi_healthy', xi_ref / (xi_healthy + 1e-12))
        Pi_acute = pm.Deterministic('Pi_acute', xi_ref / (xi_acute + 1e-12))
        Pi_chronic = pm.Deterministic('Pi_chronic', xi_ref / (xi_chronic + 1e-12))
        
        Gamma_healthy = coh_healthy ** pm.math.abs(gamma_coh)
        Gamma_acute = coh_acute ** pm.math.abs(gamma_coh)
        Gamma_chronic = coh_chronic ** pm.math.abs(gamma_coh)
        
        NAA_base = 1.105
        
        NAA_healthy = NAA_base * Pi_healthy * Gamma_healthy
        NAA_acute = NAA_base * Pi_acute * Gamma_acute * viral_damage_acute
        NAA_chronic = NAA_base * Pi_chronic * Gamma_chronic * viral_damage_chronic
        
        NAA_pred = pm.math.stack([NAA_healthy, NAA_acute, NAA_chronic])
        
        Cho_base = 0.225
        Cho_pred = pm.math.stack([Cho_base, Cho_base * membrane_acute,
                                  Cho_base * membrane_chronic])
        
        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.08)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.04)
        
        NAA_obs = pm.Normal('NAA_obs', mu=NAA_pred, sigma=sigma_NAA, observed=NAA_OBS)
        Cho_obs = pm.Normal('Cho_obs', mu=Cho_pred, sigma=sigma_Cho, observed=CHO_OBS)
        
        protection_ratio = pm.Deterministic('protection_ratio', 
                                           Pi_acute / (Pi_chronic + 1e-12))
        
    return model


# ============================================================================
# SAMPLING WITH PROPER LOG_LIKELIHOOD
# ============================================================================

def sample_model_fixed(model, model_name, n_samples=1000, n_tune=1500):
    """
    Sample with settings optimized to avoid divergences.
    CRITICAL: Computes log_likelihood properly for WAIC/LOO.
    """
    
    print(f"\n{'='*80}")
    print(f" Sampling: {model_name}")
    print(f"{'='*80}\n")
    
    with model:
        # Sample with conservative settings
        idata = pm.sample(
            draws=n_samples,
            tune=n_tune,
            chains=4,
            target_accept=0.85,  # Lower = more exploration, fewer divergences
            max_treedepth=12,    # Allow deeper trees
            return_inferencedata=True,
            random_seed=42,
            progressbar=True,
            idata_kwargs={'log_likelihood': True}  # CRITICAL!
        )
        
        # Check divergences
        n_div = idata.sample_stats.diverging.sum().values
        div_rate = 100 * n_div / idata.sample_stats.diverging.size
        
        print(f"\n{'='*60}")
        print(f"Sampling Complete: {model_name}")
        print(f"  Divergences: {n_div} ({div_rate:.1f}%)")
        
        if n_div == 0:
            print(f"  ✓ EXCELLENT: No divergences")
        elif n_div < 50:
            print(f"  ✓ GOOD: Acceptable divergence rate")
        elif n_div < 200:
            print(f"  ⚠ WARNING: Many divergences but may be usable")
        else:
            print(f"  ❌ PROBLEM: Too many divergences, results unreliable")
        
        print(f"{'='*60}\n")
        
        # Generate posterior predictive
        print(f"Generating posterior predictive for {model_name}...")
        pm.sample_posterior_predictive(idata, extend_inferencedata=True, random_seed=42)
        
    return idata


# ============================================================================
# MAIN COMPARISON
# ============================================================================

def main():
    """Run complete comparison with fixed models."""
    
    # Output directory (Mac-compatible)
    output_dir = Path.cwd() / 'results' / 'model_comparison_fixed'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Output directory: {output_dir}")
    
    # Build and sample models
    models = {
        'Full_Fixed': build_full_model_fixed(),
        'No_xi': build_no_xi_model(),
        'Linear_xi': build_linear_xi_model()
    }
    
    idata_dict = {}
    
    for name, model in models.items():
        try:
            idata = sample_model_fixed(model, name, n_samples=1000, n_tune=1500)
            idata_dict[name] = idata
            
            # Save immediately
            save_path = output_dir / f'{name}_trace.nc'
            idata.to_netcdf(save_path)
            print(f"✓ Saved: {save_path}\n")
            
        except Exception as e:
            print(f"❌ {name} failed: {e}\n")
            continue
    
    # Compare models
    if len(idata_dict) >= 2:
        print("\n" + "="*80)
        print(" MODEL COMPARISON: WAIC")
        print("="*80 + "\n")
        
        try:
            comp_waic = az.compare(idata_dict, ic='waic', scale='deviance')
            print(comp_waic)
            print()
            
            comp_waic.to_csv(output_dir / 'waic_comparison.csv')
            print(f"✓ Saved WAIC comparison")
            
        except Exception as e:
            print(f"WAIC comparison failed: {e}")
        
        print("\n" + "="*80)
        print(" MODEL COMPARISON: LOO")
        print("="*80 + "\n")
        
        try:
            comp_loo = az.compare(idata_dict, ic='loo', scale='deviance')
            print(comp_loo)
            print()
            
            comp_loo.to_csv(output_dir / 'loo_comparison.csv')
            print(f"✓ Saved LOO comparison")
            
        except Exception as e:
            print(f"LOO comparison failed: {e}")
        
        # Create comparison plot
        try:
            fig, axes = plt.subplots(2, 1, figsize=(10, 8))
            
            # WAIC plot
            az.plot_compare(comp_waic, insample_dev=False, ax=axes[0])
            axes[0].set_title('WAIC Comparison', fontweight='bold')
            
            # LOO plot
            az.plot_compare(comp_loo, insample_dev=False, ax=axes[1])
            axes[1].set_title('LOO Comparison', fontweight='bold')
            
            plt.tight_layout()
            plt.savefig(output_dir / 'model_comparison.png', dpi=300, bbox_inches='tight')
            print(f"✓ Saved comparison plot\n")
            
        except Exception as e:
            print(f"Plotting failed: {e}")
    
    print("\n" + "="*80)
    print(" ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {output_dir}")
    print("\nNext steps:")
    print("  1. Check divergence rates (should be < 50)")
    print("  2. Review WAIC comparison table")
    print("  3. Examine parameter summaries")
    print()


if __name__ == '__main__':
    main()
