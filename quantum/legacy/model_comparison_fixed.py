"""
WAIC Model Comparison - FIXED VERSION
======================================

Fixes all issues from previous run:
1. Non-centered parameterization for ξ (eliminates divergences)
2. Proper log_likelihood computation
3. Cross-platform path handling
4. Robust error handling
5. Informative progress messages

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

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

# ============================================================================
# OBSERVED DATA (Sailasuta et al. 2012)
# ============================================================================

CONDITIONS = ['healthy', 'acute', 'chronic']
NAA_OBS = np.array([1.105, 1.135, 1.005])
CHO_OBS = np.array([0.225, 0.245, 0.235])

print("\n" + "="*80)
print(" FIXED WAIC MODEL COMPARISON")
print("="*80)
print("\nKey improvements:")
print("  ✓ Non-centered ξ parameterization")
print("  ✓ Proper log_likelihood computation")
print("  ✓ Cross-platform path handling")
print("  ✓ Target_accept = 0.92 (balanced)")
print()

# ============================================================================
# ROBUST MODEL BUILDING
# ============================================================================

def build_full_model_robust():
    """
    Full model with FIXED parameterization.
    
    Key changes:
    1. Non-centered ξ (mean + offsets)
    2. Wider bounds
    3. Better prior choices
    """
    
    with pm.Model() as model:
        # ===== NON-CENTERED ξ PARAMETERIZATION =====
        # Sample mean and offsets (eliminates correlation!)
        xi_mean = pm.TruncatedNormal('xi_mean', mu=0.6e-9, sigma=0.15e-9,
                                      lower=0.4e-9, upper=0.9e-9)
        
        # Offsets from mean (always positive)
        delta_acute = pm.HalfNormal('delta_acute', sigma=0.1e-9)
        delta_chronic = pm.HalfNormal('delta_chronic', sigma=0.1e-9)
        
        # Derived values (guaranteed ordering)
        xi_healthy = pm.Deterministic('xi_healthy', xi_mean)
        xi_acute = pm.Deterministic('xi_acute', xi_mean - delta_acute)
        xi_chronic = pm.Deterministic('xi_chronic', xi_mean + delta_chronic)
        
        # ===== OTHER PARAMETERS =====
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
        
        # ===== FORWARD MODEL =====
        xi_ref = 0.8e-9
        coh_healthy = 0.85
        coh_acute = 0.70
        coh_chronic = 0.75
        
        # Protection factors
        Pi_healthy = pm.Deterministic('Pi_healthy', (xi_ref / xi_healthy) ** beta_xi)
        Pi_acute = pm.Deterministic('Pi_acute', (xi_ref / xi_acute) ** beta_xi)
        Pi_chronic = pm.Deterministic('Pi_chronic', (xi_ref / xi_chronic) ** beta_xi)
        
        # Coherence enhancement
        Gamma_healthy = coh_healthy ** gamma_coh
        Gamma_acute = coh_acute ** gamma_coh
        Gamma_chronic = coh_chronic ** gamma_coh
        
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
        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.06)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.03)
        
        NAA_obs = pm.Normal('NAA_obs', mu=NAA_pred, sigma=sigma_NAA, observed=NAA_OBS)
        Cho_obs = pm.Normal('Cho_obs', mu=Cho_pred, sigma=sigma_Cho, observed=CHO_OBS)
        
        # Derived
        protection_ratio = pm.Deterministic('protection_ratio', Pi_acute / Pi_chronic)
        
    return model


def build_no_xi_model():
    """No ξ coupling - quantum effects removed."""
    
    with pm.Model() as model:
        gamma_coh = pm.TruncatedNormal('gamma_coh', mu=1.0, sigma=0.3,
                                        lower=0.3, upper=2.0)
        
        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)
        
        membrane_acute = pm.TruncatedNormal('membrane_acute', mu=1.1, sigma=0.04,
                                            lower=1.0, upper=1.25)
        membrane_chronic = pm.TruncatedNormal('membrane_chronic', mu=1.03, sigma=0.03,
                                               lower=1.0, upper=1.15)
        
        # Forward model - NO ξ PROTECTION
        coh_healthy = 0.85
        coh_acute = 0.70
        coh_chronic = 0.75
        
        Gamma_healthy = coh_healthy ** gamma_coh
        Gamma_acute = coh_acute ** gamma_coh
        Gamma_chronic = coh_chronic ** gamma_coh
        
        NAA_base = 1.105
        
        # No Pi factors!
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
    """Linear ξ coupling (β = 1 fixed)."""
    
    with pm.Model() as model:
        # Non-centered ξ
        xi_mean = pm.TruncatedNormal('xi_mean', mu=0.6e-9, sigma=0.15e-9,
                                      lower=0.4e-9, upper=0.9e-9)
        delta_acute = pm.HalfNormal('delta_acute', sigma=0.1e-9)
        delta_chronic = pm.HalfNormal('delta_chronic', sigma=0.1e-9)
        
        xi_healthy = pm.Deterministic('xi_healthy', xi_mean)
        xi_acute = pm.Deterministic('xi_acute', xi_mean - delta_acute)
        xi_chronic = pm.Deterministic('xi_chronic', xi_mean + delta_chronic)
        
        # FIXED β = 1
        beta_xi = 1.0
        
        gamma_coh = pm.TruncatedNormal('gamma_coh', mu=1.0, sigma=0.3,
                                        lower=0.3, upper=2.0)
        
        viral_damage_acute = pm.Beta('viral_damage_acute', alpha=9, beta=1)
        viral_damage_chronic = pm.Beta('viral_damage_chronic', alpha=19, beta=1)
        
        membrane_acute = pm.TruncatedNormal('membrane_acute', mu=1.1, sigma=0.04,
                                            lower=1.0, upper=1.25)
        membrane_chronic = pm.TruncatedNormal('membrane_chronic', mu=1.03, sigma=0.03,
                                               lower=1.0, upper=1.15)
        
        # Forward model with LINEAR protection
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
# ROBUST SAMPLING
# ============================================================================

def sample_model_robust(model, model_name, n_samples=2000):
    """
    Sample with proper error handling and log_likelihood computation.
    """
    
    print(f"\n{'='*80}")
    print(f" SAMPLING: {model_name}")
    print(f"{'='*80}\n")
    
    with model:
        # Sample
        idata = pm.sample(
            draws=n_samples,
            tune=1500,
            chains=4,
            target_accept=0.92,  # Balanced
            return_inferencedata=True,
            idata_kwargs={'log_likelihood': True},  # KEY FIX!
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
        if n_div < 500:  # Only if not catastrophic
            print("Generating posterior predictive...")
            pm.sample_posterior_predictive(idata, extend_inferencedata=True, 
                                           random_seed=42, progressbar=False)
    
    return idata, n_div


# ============================================================================
# COMPARISON ANALYSIS
# ============================================================================

def run_comparison():
    """Run complete WAIC comparison."""
    
    # Build models
    models = {
        'Full Model': build_full_model_robust(),
        'No ξ Coupling': build_no_xi_model(),
        'Linear ξ (β=1)': build_linear_xi_model()
    }
    
    # Sample all models
    results = {}
    div_counts = {}
    
    for name, model in models.items():
        try:
            idata, n_div = sample_model_robust(model, name, n_samples=2000)
            results[name] = idata
            div_counts[name] = n_div
        except Exception as e:
            print(f"❌ {name} failed: {e}")
            continue
    
    # Check if we have valid results
    if len(results) == 0:
        print("\n❌ NO MODELS SUCCEEDED - Check error messages above")
        return None, None
    
    print(f"\n✓ {len(results)}/{len(models)} models succeeded\n")
    
    # Compute WAIC
    print("\n" + "="*80)
    print(" WAIC COMPARISON")
    print("="*80 + "\n")
    
    try:
        comparison = az.compare(results, ic='waic', scale='deviance')
        print(comparison)
        print()
    except Exception as e:
        print(f"❌ WAIC comparison failed: {e}")
        comparison = None
    
    # Compute LOO
    print("\n" + "="*80)
    print(" LOO COMPARISON")
    print("="*80 + "\n")
    
    try:
        comparison_loo = az.compare(results, ic='loo', scale='deviance')
        print(comparison_loo)
        print()
    except Exception as e:
        print(f"❌ LOO comparison failed: {e}")
        comparison_loo = None
    
    return results, comparison, comparison_loo, div_counts


# ============================================================================
# VISUALIZATION
# ============================================================================

def plot_results(results, comparison, div_counts):
    """Create comparison plots."""
    
    if comparison is None:
        print("Cannot plot - no comparison data")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: WAIC comparison
    ax = axes[0, 0]
    models = comparison.index.tolist()
    waic = comparison['waic'].values
    se = comparison['se'].values
    
    colors = ['#2ecc71' if i == 0 else '#3498db' for i in range(len(models))]
    y_pos = np.arange(len(models))
    
    ax.barh(y_pos, waic, xerr=se, color=colors, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(models)
    ax.set_xlabel('WAIC (lower is better)', fontweight='bold')
    ax.set_title('A. Model Comparison (WAIC)', fontweight='bold')
    ax.grid(alpha=0.3, axis='x')
    
    # Panel 2: Divergences
    ax = axes[0, 1]
    model_names = list(div_counts.keys())
    divs = list(div_counts.values())
    
    colors_div = ['#2ecc71' if d < 50 else '#e74c3c' for d in divs]
    ax.bar(range(len(model_names)), divs, color=colors_div, alpha=0.7)
    ax.set_xticks(range(len(model_names)))
    ax.set_xticklabels(model_names, rotation=45, ha='right')
    ax.set_ylabel('Number of Divergences', fontweight='bold')
    ax.set_title('B. Sampling Diagnostics', fontweight='bold')
    ax.axhline(50, color='orange', linestyle='--', label='Acceptable threshold')
    ax.legend()
    ax.grid(alpha=0.3, axis='y')
    
    # Panel 3: Prediction errors
    ax = axes[1, 0]
    for name, idata in results.items():
        ppc = idata.posterior_predictive
        NAA_pred = ppc['NAA_obs'].mean(dim=['chain', 'draw']).values
        errors = 100 * (NAA_pred - NAA_OBS) / NAA_OBS
        ax.plot(CONDITIONS, errors, marker='o', label=name, linewidth=2)
    
    ax.axhline(0, color='black', linestyle='--', alpha=0.3)
    ax.set_ylabel('Prediction Error (%)', fontweight='bold')
    ax.set_xlabel('Condition', fontweight='bold')
    ax.set_title('C. NAA Prediction Errors', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Panel 4: Summary table
    ax = axes[1, 1]
    ax.axis('off')
    
    summary_text = f"""
    MODEL COMPARISON SUMMARY
    
    Best Model: {comparison.index[0]}
    WAIC: {comparison['waic'].iloc[0]:.1f} ± {comparison['se'].iloc[0]:.1f}
    
    Δ WAIC from next best:
    {comparison['d_waic'].iloc[1]:.1f}
    
    Interpretation:
    """
    
    delta = comparison['d_waic'].iloc[1] if len(comparison) > 1 else 0
    if delta > 10:
        summary_text += "DECISIVE evidence for best model"
    elif delta > 6:
        summary_text += "STRONG evidence for best model"
    elif delta > 2:
        summary_text += "WEAK evidence for best model"
    else:
        summary_text += "Models indistinguishable"
    
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    
    # Save
    output_dir = Path('results/model_comparison_fixed')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    plt.savefig(output_dir / 'comparison_results.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/comparison_results.png")
    
    return fig


# ============================================================================
# MAIN
# ============================================================================

if __name__ == '__main__':
    results, comparison, comparison_loo, div_counts = run_comparison()
    
    if results is not None and comparison is not None:
        fig = plot_results(results, comparison, div_counts)
        
        # Save results
        output_dir = Path('results/model_comparison_fixed')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        comparison.to_csv(output_dir / 'waic_comparison.csv')
        if comparison_loo is not None:
            comparison_loo.to_csv(output_dir / 'loo_comparison.csv')
        
        for name, idata in results.items():
            safe_name = name.replace(' ', '_').replace('(', '').replace(')', '')
            idata.to_netcdf(output_dir / f'{safe_name}_trace.nc')
        
        print("\n" + "="*80)
        print(" ANALYSIS COMPLETE!")
        print("="*80)
        print(f"\nResults in: {output_dir}")
        print("\nKey finding:")
        print(f"  Best model: {comparison.index[0]}")
        print(f"  WAIC: {comparison['waic'].iloc[0]:.1f}")
        print(f"  Δ WAIC: {comparison['d_waic'].iloc[1]:.1f}")
        
        plt.show()
    
    else:
        print("\n❌ Analysis failed - check error messages above")
