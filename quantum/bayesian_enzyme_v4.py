"""
BAYESIAN INFERENCE v4.0: ENZYME KINETICS MODEL
===============================================

Replaces phenomenological compensation with mechanistic enzyme modulation.

KEY CHANGES FROM v3.6:
- Uses enzyme_kinetics.py for forward model
- Infers enzyme parameters instead of abstract compensation
- More testable predictions
- Better mechanistic understanding

WORKFLOW:
1. Import enzyme_kinetics module
2. Define PyMC model with enzyme parameters as priors
3. Run MCMC sampling
4. Compare to v3.6 results
"""

import numpy as np
import pymc as pm
import arviz as az
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Import your enzyme kinetics module
from enzyme_kinetics import (
    EnzymeKinetics,
    compute_protection_factor,
    coherence_modulation,
    ENZYME
)


# =============================================================================
# CLINICAL DATA (Sailasuta et al. 2012)
# =============================================================================

CLINICAL_DATA = {
    'healthy': {'NAA': 1.105, 'Cho': 0.225},
    'acute': {'NAA': 1.135, 'Cho': 0.245},
    'chronic': {'NAA': 1.005, 'Cho': 0.235}
}

# Convert to arrays for easier handling
CONDITIONS = ['healthy', 'acute', 'chronic']
NAA_OBS = np.array([CLINICAL_DATA[c]['NAA'] for c in CONDITIONS])
CHO_OBS = np.array([CLINICAL_DATA[c]['Cho'] for c in CONDITIONS])


import pytensor.tensor as pt
from pytensor.compile.ops import as_op

# =============================================================================
# FORWARD MODEL: ENZYME KINETICS
# =============================================================================

@as_op(itypes=[pt.dscalar] * 8, otypes=[pt.dvector, pt.dvector])
def forward_model_enzyme_op(xi_acute, xi_chronic, beta_xi, gamma_coh,
                           viral_damage_acute, viral_damage_chronic,
                           membrane_acute, membrane_chronic):
    """
    Wrapper for forward_model_enzyme to be used in PyMC.
    """
    return forward_model_enzyme(
        xi_acute=xi_acute,
        xi_chronic=xi_chronic,
        beta_xi=beta_xi,
        gamma_coh=gamma_coh,
        viral_damage_acute=viral_damage_acute,
        viral_damage_chronic=viral_damage_chronic,
        membrane_acute=membrane_acute,
        membrane_chronic=membrane_chronic
    )

def forward_model_enzyme(xi_acute, xi_chronic, beta_xi, gamma_coh,
                        viral_damage_acute, viral_damage_chronic,
                        membrane_acute, membrane_chronic,
                        coh_acute=0.95, coh_chronic=0.80):
    """
    Forward model using enzyme kinetics.
    
    Parameters
    ----------
    xi_acute : float
        Correlation length in acute HIV (m)
    xi_chronic : float
        Correlation length in chronic HIV (m)
    beta_xi : float
        Protection factor exponent
    gamma_coh : float
        Coherence coupling exponent
    viral_damage_acute : float
        Viral damage factor in acute (0-1)
    viral_damage_chronic : float
        Viral damage factor in chronic (0-1)
    membrane_acute : float
        Membrane turnover in acute (>1 = elevated)
    membrane_chronic : float
        Membrane turnover in chronic (>1 = elevated)
    coh_acute : float
        Coherence in acute phase
    coh_chronic : float
        Coherence in chronic phase
        
    Returns
    -------
    NAA_pred : array
        [NAA_healthy, NAA_acute, NAA_chronic] in MRS units
    Cho_pred : array
        [Cho_healthy, Cho_acute, Cho_chronic] in MRS units
    """
    
    # Healthy baseline
    xi_healthy = 0.75e-9  # Reference
    Pi_healthy = compute_protection_factor(xi_healthy, beta_xi=beta_xi)
    eta_healthy = coherence_modulation(0.85, gamma=gamma_coh)
    
    enzymes_healthy = EnzymeKinetics(
        Pi_xi=Pi_healthy,
        eta_coh=eta_healthy,
        viral_damage_factor=1.0
    )
    NAA_h, Cho_h = enzymes_healthy.integrate(
        duration_days=60,
        membrane_turnover=1.0
    )
    
    # Acute HIV
    Pi_acute = compute_protection_factor(xi_acute, beta_xi=beta_xi)
    eta_acute = coherence_modulation(coh_acute, gamma=gamma_coh)
    
    enzymes_acute = EnzymeKinetics(
        Pi_xi=Pi_acute,
        eta_coh=eta_acute,
        viral_damage_factor=viral_damage_acute
    )
    NAA_a, Cho_a = enzymes_acute.integrate(
        duration_days=60,
        membrane_turnover=membrane_acute
    )
    
    # Chronic HIV
    Pi_chronic = compute_protection_factor(xi_chronic, beta_xi=beta_xi)
    eta_chronic = coherence_modulation(coh_chronic, gamma=gamma_coh)
    
    enzymes_chronic = EnzymeKinetics(
        Pi_xi=Pi_chronic,
        eta_coh=eta_chronic,
        viral_damage_factor=viral_damage_chronic
    )
    NAA_c, Cho_c = enzymes_chronic.integrate(
        duration_days=60,
        membrane_turnover=membrane_chronic
    )
    
    # Convert from molar to MRS units (relative to creatine)
    # Assume creatine = 8 mM constant
    creatine = 8.0e-3
    
    NAA_pred = np.array([NAA_h, NAA_a, NAA_c]) / creatine
    Cho_pred = np.array([Cho_h, Cho_a, Cho_c]) / creatine
    
    return NAA_pred, Cho_pred


# =============================================================================
# BAYESIAN MODEL
# =============================================================================

def build_enzyme_model():
    """
    Build PyMC model with enzyme kinetics.
    
    KEY PARAMETERS:
    - xi_acute, xi_chronic: Noise correlation lengths
    - beta_xi: Protection factor exponent (expect ~2)
    - gamma_coh: Coherence coupling exponent
    - viral_damage_*: Direct damage to enzymes
    - membrane_*: Membrane turnover rates
    """
    
    with pm.Model() as model:
        
        # =====================================================================
        # PRIORS
        # =====================================================================
        
        # Noise correlation lengths (informed by v3.6)
        xi_acute = pm.TruncatedNormal(
            'xi_acute',
            mu=0.50e-9,
            sigma=0.15e-9,
            lower=0.35e-9,
            upper=0.70e-9
        )
        
        xi_chronic = pm.TruncatedNormal(
            'xi_chronic',
            mu=0.78e-9,
            sigma=0.10e-9,
            lower=0.70e-9,
            upper=0.90e-9
        )
        
        # Protection factor exponent (informed by v3.6: β = 1.731)
        beta_xi = pm.TruncatedNormal(
            'beta_xi',
            mu=1.75,
            sigma=0.50,
            lower=0.5,
            upper=3.5
        )
        
        # Coherence coupling exponent
        gamma_coh = pm.TruncatedNormal(
            'gamma_coh',
            mu=1.5,
            sigma=0.50,
            lower=0.5,
            upper=3.0
        )
        
        # Viral damage factors (0-1, 1 = no damage)
        viral_damage_acute = pm.Beta(
            'viral_damage_acute',
            alpha=19,  # Mean ~0.95
            beta=1
        )
        
        viral_damage_chronic = pm.Beta(
            'viral_damage_chronic',
            alpha=9,  # Mean ~0.90
            beta=1
        )
        
        # Membrane turnover (>1 = elevated)
        membrane_acute = pm.TruncatedNormal(
            'membrane_acute',
            mu=2.0,
            sigma=0.5,
            lower=1.0,
            upper=4.0
        )
        
        membrane_chronic = pm.TruncatedNormal(
            'membrane_chronic',
            mu=1.2,
            sigma=0.3,
            lower=1.0,
            upper=2.0
        )
        
        # =====================================================================
        # FORWARD MODEL
        # =====================================================================
        
        NAA_pred, Cho_pred = forward_model_enzyme_op(
            xi_acute,
            xi_chronic,
            beta_xi,
            gamma_coh,
            viral_damage_acute,
            viral_damage_chronic,
            membrane_acute,
            membrane_chronic
        )
        
        # =====================================================================
        # LIKELIHOOD
        # =====================================================================
        
        # Observation noise (similar to v3.6)
        sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.06)
        sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.03)
        
        # Likelihood
        NAA_likelihood = pm.Normal(
            'NAA_obs',
            mu=NAA_pred,
            sigma=sigma_NAA,
            observed=NAA_OBS
        )
        
        Cho_likelihood = pm.Normal(
            'Cho_obs',
            mu=Cho_pred,
            sigma=sigma_Cho,
            observed=CHO_OBS
        )
        
        # =====================================================================
        # DERIVED QUANTITIES
        # =====================================================================
        
        # Protection factors for reporting
        Pi_acute = pm.Deterministic(
            'Pi_acute',
            (0.8e-9 / xi_acute) ** beta_xi
        )
        
        Pi_chronic = pm.Deterministic(
            'Pi_chronic',
            (0.8e-9 / xi_chronic) ** beta_xi
        )
        
        # Protection ratio
        protection_ratio = pm.Deterministic(
            'protection_ratio',
            Pi_acute / Pi_chronic
        )
        
    return model


# =============================================================================
# SAMPLING
# =============================================================================

def run_inference(n_samples=2000, n_chains=4, target_accept=0.99):
    """
    Run Bayesian inference with enzyme model.
    
    Parameters
    ----------
    n_samples : int
        Number of samples per chain
    n_chains : int
        Number of MCMC chains
    target_accept : float
        Target acceptance rate
        
    Returns
    -------
    idata : InferenceData
        ArviZ InferenceData object with results
    """
    
    print("=" * 80)
    print(" BAYESIAN INFERENCE v4.0: ENZYME KINETICS MODEL")
    print("=" * 80)
    print()
    
    # Build model
    print("Building model...")
    model = build_enzyme_model()
    
    # Sample
    print(f"Sampling: {n_samples} × {n_chains} chains...")
    print(f"Target accept: {target_accept}")
    print()
    
    with model:
        idata = pm.sample(
            draws=n_samples,
            tune=1000,
            chains=n_chains,
            target_accept=target_accept,
            return_inferencedata=True,
            random_seed=42
        )
        
        # Posterior predictive
        print("\nGenerating posterior predictive...")
        pm.sample_posterior_predictive(
            idata,
            extend_inferencedata=True,
            random_seed=42
        )
    
    return idata


# =============================================================================
# ANALYSIS
# =============================================================================

def analyze_results(idata):
    """
    Analyze and compare enzyme model results to v3.6.
    
    Parameters
    ----------
    idata : InferenceData
        Inference results
    """
    
    print("\n" + "=" * 80)
    print(" ENZYME MODEL RESULTS")
    print("=" * 80)
    print()
    
    # Summary statistics
    summary = az.summary(idata, var_names=[
        'xi_acute', 'xi_chronic', 'beta_xi', 'gamma_coh',
        'viral_damage_acute', 'viral_damage_chronic',
        'membrane_acute', 'membrane_chronic',
        'Pi_acute', 'Pi_chronic', 'protection_ratio'
    ])
    
    print(summary)
    print()
    
    # Key findings
    print("\n" + "=" * 80)
    print(" KEY FINDINGS")
    print("=" * 80)
    print()
    
    # Protection factor exponent
    beta_samples = idata.posterior['beta_xi'].values.flatten()
    beta_median = np.median(beta_samples)
    beta_hdi = az.hdi(idata, var_names=['beta_xi'], hdi_prob=0.94)['beta_xi'].values
    
    print(f"Protection Factor Exponent β_ξ:")
    print(f"  Median: {beta_median:.3f}")
    print(f"  94% HDI: [{beta_hdi[0]:.3f}, {beta_hdi[1]:.3f}]")
    print()
    
    # Noise correlation lengths
    xi_acute_samples = idata.posterior['xi_acute'].values.flatten() * 1e9
    xi_chronic_samples = idata.posterior['xi_chronic'].values.flatten() * 1e9
    
    print(f"Noise Correlation Lengths:")
    print(f"  ξ_acute:   {np.median(xi_acute_samples):.3f} nm")
    print(f"  ξ_chronic: {np.median(xi_chronic_samples):.3f} nm")
    print()
    
    # P(ξ_acute < ξ_chronic)
    p_acute_less = np.mean(xi_acute_samples < xi_chronic_samples)
    print(f"P(ξ_acute < ξ_chronic) = {p_acute_less:.4f}")
    print()
    
    # Prediction errors
    ppc = idata.posterior_predictive
    NAA_pred_mean = ppc['NAA_obs'].mean(dim=['chain', 'draw']).values
    Cho_pred_mean = ppc['Cho_obs'].mean(dim=['chain', 'draw']).values
    
    NAA_errors = 100 * (NAA_pred_mean - NAA_OBS) / NAA_OBS
    Cho_errors = 100 * (Cho_pred_mean - CHO_OBS) / CHO_OBS
    
    print("Prediction Errors:")
    print(f"  {'Condition':12s}  {'NAA Error':>12s}  {'Cho Error':>12s}")
    print("  " + "-" * 40)
    for i, cond in enumerate(CONDITIONS):
        print(f"  {cond:12s}  {NAA_errors[i]:>11.2f}%  {Cho_errors[i]:>11.2f}%")
    print()
    
    # Comparison to v3.6
    print("\n" + "=" * 80)
    print(" COMPARISON TO v3.6")
    print("=" * 80)
    print()
    
    print("v3.6 Results:")
    print("  β_ξ = 1.731 (94% HDI: [0.846, 2.792])")
    print("  Chronic NAA error: +0.4%")
    print()
    
    print("v4.0 Results (Enzyme Model):")
    print(f"  β_ξ = {beta_median:.3f} (94% HDI: [{beta_hdi[0]:.3f}, {beta_hdi[1]:.3f}])")
    print(f"  Chronic NAA error: {NAA_errors[2]:+.1f}%")
    print()
    
    # Model comparison
    rms_error = np.sqrt(np.mean(NAA_errors**2 + Cho_errors**2))
    print(f"RMS Error: {rms_error:.2f}%")
    print()
    
    return summary


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_results(idata, output_dir='results/enzyme_v4'):
    """
    Create publication-quality figures.
    
    Parameters
    ----------
    idata : InferenceData
        Inference results
    output_dir : str
        Output directory for figures
    """
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Figure 1: Posterior distributions
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    az.plot_posterior(
        idata,
        var_names=['beta_xi', 'xi_acute', 'xi_chronic',
                   'Pi_acute', 'Pi_chronic', 'protection_ratio'],
        hdi_prob=0.94,
        ax=axes.flatten()
    )
    
    plt.tight_layout()
    plt.savefig(output_path / 'v4_posteriors.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 2: Predicted vs Observed
    ppc = idata.posterior_predictive
    NAA_pred_mean = ppc['NAA_obs'].mean(dim=['chain', 'draw']).values
    Cho_pred_mean = ppc['Cho_obs'].mean(dim=['chain', 'draw']).values
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # NAA
    ax1.scatter(NAA_OBS, NAA_pred_mean, s=100)
    ax1.plot([0.9, 1.2], [0.9, 1.2], 'k--', alpha=0.5)
    ax1.set_xlabel('Observed NAA')
    ax1.set_ylabel('Predicted NAA')
    ax1.set_title('NAA: Enzyme Model v4.0')
    
    for i, cond in enumerate(CONDITIONS):
        ax1.annotate(cond, (NAA_OBS[i], NAA_pred_mean[i]),
                    xytext=(5, 5), textcoords='offset points')
    
    # Cho
    ax2.scatter(CHO_OBS, Cho_pred_mean, s=100)
    ax2.plot([0.2, 0.26], [0.2, 0.26], 'k--', alpha=0.5)
    ax2.set_xlabel('Observed Cho')
    ax2.set_ylabel('Predicted Cho')
    ax2.set_title('Cho: Enzyme Model v4.0')
    
    for i, cond in enumerate(CONDITIONS):
        ax2.annotate(cond, (CHO_OBS[i], Cho_pred_mean[i]),
                    xytext=(5, 5), textcoords='offset points')
    
    plt.tight_layout()
    plt.savefig(output_path / 'v4_pred_vs_obs.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Figures saved to {output_path}/")


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Run complete enzyme model inference."""
    
    # Run inference
    idata = run_inference(n_samples=2000, n_chains=4)
    
    # Analyze results
    summary = analyze_results(idata)
    
    # Save results
    output_dir = Path('results/enzyme_v4')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    idata.to_netcdf(output_dir / 'trace_v4.nc')
    summary.to_csv(output_dir / 'summary_v4.csv')
    
    # Plot results
    plot_results(idata, output_dir=str(output_dir))
    
    print("\n" + "=" * 80)
    print(" INFERENCE COMPLETE")
    print("=" * 80)
    print()
    print(f"Results saved to {output_dir}/")
    print()
    
    return idata


if __name__ == "__main__":
    idata = main()
