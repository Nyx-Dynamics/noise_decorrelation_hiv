"""
Model Comparison - INDIVIDUAL LEVEL DATA
=========================================

Compares three models using individual patient data (n=62 acute HIV patients):
1. Full Model: ξ coupling with β_ξ estimated (nonlinear protection)
2. Linear Model: ξ coupling with β_ξ = 1 (linear protection)
3. No Coupling Model: No ξ effect (null model)

Uses WAIC for model comparison.

Data: Individual patient measurements from VALCOUR_2015_INDIVIDUAL_PATIENTS.csv
      Basal ganglia NAA (BGNAA column)

Author: AC
Date: 2024-11-15
"""

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Configuration
sns.set_style("whitegrid")
RESULTS_DIR = Path("quantum/results/model_comparison_individual")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

print("\n" + "="*80)
print(" MODEL COMPARISON - INDIVIDUAL LEVEL DATA")
print("="*80)

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Load individual patient data
data = pd.read_csv("../data/individual/VALCOUR_2015_INDIVIDUAL_PATIENTS.csv")

# Extract basal ganglia NAA and viral load
NAA_obs = data['BGNAA'].values  # Individual measurements
logVL = data['logpVL'].values   # Log viral load (proxy for inflammation)
n_patients = len(NAA_obs)

# Remove any NaN values
mask = ~np.isnan(NAA_obs) & ~np.isnan(logVL)
NAA_obs = NAA_obs[mask]
logVL = logVL[mask]
n_patients = len(NAA_obs)

print(f"\nIndividual Patient Data:")
print(f"  N = {n_patients} acute HIV patients")
print(f"  NAA mean: {NAA_obs.mean():.3f} ± {NAA_obs.std():.3f}")
print(f"  NAA range: [{NAA_obs.min():.3f}, {NAA_obs.max():.3f}]")
print(f"  log(VL) mean: {logVL.mean():.3f} ± {logVL.std():.3f}")

# Standardize viral load for modeling
logVL_std = (logVL - logVL.mean()) / logVL.std()

# ============================================================================
# MODEL DEFINITIONS
# ============================================================================

def build_full_model():
    """Full model: ξ depends on viral load with estimated β_ξ (nonlinear)"""
    with pm.Model() as model:
        # Priors
        NAA_baseline = pm.Normal("NAA_baseline", mu=10.0, sigma=1.0)
        
        # ξ as function of viral load (higher VL → potentially shorter ξ)
        xi_intercept = pm.Normal("xi_intercept", mu=0.7, sigma=0.2)
        xi_slope = pm.Normal("xi_slope", mu=-0.1, sigma=0.05)  # Negative: higher VL → shorter ξ
        
        # Individual ξ values
        xi_i = pm.Deterministic("xi_i", 
            pm.math.exp(xi_intercept + xi_slope * logVL_std))
        
        # Protection exponent (key parameter)
        beta_xi = pm.Normal("beta_xi", mu=2.0, sigma=0.5)
        
        # NAA prediction for each patient
        NAA_pred = pm.Deterministic("NAA_pred",
            NAA_baseline * (1 / xi_i)**beta_xi)
        
        # Likelihood (individual observations)
        sigma_obs = pm.HalfNormal("sigma_obs", sigma=2.0)
        pm.Normal("NAA_likelihood", mu=NAA_pred, sigma=sigma_obs, observed=NAA_obs)
        
        # Compute pointwise log-likelihood for WAIC
        pm.Deterministic("log_lik", pm.logp(pm.Normal.dist(mu=NAA_pred, sigma=sigma_obs), NAA_obs))
        
    return model

def build_linear_model():
    """Linear model: ξ depends on viral load with β_ξ = 1 (linear)"""
    with pm.Model() as model:
        # Priors (same as full model)
        NAA_baseline = pm.Normal("NAA_baseline", mu=10.0, sigma=1.0)
        
        xi_intercept = pm.Normal("xi_intercept", mu=0.7, sigma=0.2)
        xi_slope = pm.Normal("xi_slope", mu=-0.1, sigma=0.05)
        
        # Individual ξ values
        xi_i = pm.Deterministic("xi_i",
            pm.math.exp(xi_intercept + xi_slope * logVL_std))
        
        # FIXED: β_ξ = 1 (linear)
        beta_xi = pm.Deterministic("beta_xi", pm.math.constant(1.0))
        
        # NAA prediction (linear coupling)
        NAA_pred = pm.Deterministic("NAA_pred",
            NAA_baseline * (1 / xi_i))
        
        # Likelihood
        sigma_obs = pm.HalfNormal("sigma_obs", sigma=2.0)
        pm.Normal("NAA_likelihood", mu=NAA_pred, sigma=sigma_obs, observed=NAA_obs)
        
        # Compute pointwise log-likelihood for WAIC
        pm.Deterministic("log_lik", pm.logp(pm.Normal.dist(mu=NAA_pred, sigma=sigma_obs), NAA_obs))
        
    return model

def build_no_coupling_model():
    """No coupling model: NAA doesn't depend on ξ (null model)"""
    with pm.Model() as model:
        # Direct relationship: NAA ~ intercept + slope * VL
        NAA_intercept = pm.Normal("NAA_intercept", mu=10.0, sigma=1.0)
        NAA_slope = pm.Normal("NAA_slope", mu=0.0, sigma=0.5)
        
        NAA_pred = pm.Deterministic("NAA_pred",
            NAA_intercept + NAA_slope * logVL_std)
        
        # Likelihood
        sigma_obs = pm.HalfNormal("sigma_obs", sigma=2.0)
        pm.Normal("NAA_likelihood", mu=NAA_pred, sigma=sigma_obs, observed=NAA_obs)
        
        # Compute pointwise log-likelihood for WAIC
        pm.Deterministic("log_lik", pm.logp(pm.Normal.dist(mu=NAA_pred, sigma=sigma_obs), NAA_obs))
        
    return model

# ============================================================================
# FIT MODELS
# ============================================================================

models = {
    'Full_Model': build_full_model(),
    'Linear_xi_β=1': build_linear_model(),
    'No_xi_Coupling': build_no_coupling_model()
}

traces = {}

print("\n" + "="*80)
print(" FITTING MODELS")
print("="*80)

for name, model in models.items():
    print(f"\n{name}:")
    print("-" * 40)
    
    with model:
        trace = pm.sample(
            draws=2000,
            tune=1000,
            chains=4,
            cores=4,
            target_accept=0.95,
            return_inferencedata=True
        )
        traces[name] = trace
        
        # Save trace
        trace.to_netcdf(RESULTS_DIR / f"{name}_trace.nc")

        # Convergence check (robust extraction)
        rhat_da = az.rhat(trace).to_array()
        # xarray DataArray -> numpy scalar via .values then float
        max_rhat = float(rhat_da.max().values)
        print(f"  Max R̂: {max_rhat:.4f}")
        
        if max_rhat > 1.01:
            print("  ⚠️  WARNING: Poor convergence!")
        else:
            print("  ✓ Converged")

# ============================================================================
# WAIC COMPARISON
# ============================================================================

print("\n" + "="*80)
print(" WAIC MODEL COMPARISON")
print("="*80)

# Compute WAIC for each model
waic_results = {}
for name, trace in traces.items():
    waic = az.waic(trace, var_name="log_lik")
    waic_results[name] = waic
    print(f"\n{name}:")
    print(f"  WAIC: {waic.waic:.2f}")
    print(f"  SE: {waic.waic_se:.2f}")
    print(f"  pWAIC: {waic.p_waic:.2f}")

# Compare models
comparison = az.compare(waic_results, ic='waic')
print("\n" + "-"*80)
print("Model Ranking (lower WAIC = better):")
print(comparison)

# Save comparison
comparison.to_csv(RESULTS_DIR / "waic_comparison.txt", sep='\t')

# ============================================================================
# PARAMETER ESTIMATES
# ============================================================================

print("\n" + "="*80)
print(" PARAMETER ESTIMATES (Full Model)")
print("="*80)

summary = az.summary(traces['Full_Model'], hdi_prob=0.94)
print(summary[['mean', 'sd', 'hdi_3%', 'hdi_97%', 'r_hat']])
summary.to_csv(RESULTS_DIR / "full_model_summary.csv")

# ============================================================================
# VISUALIZATION
# ============================================================================

print("\n" + "="*80)
print(" GENERATING FIGURES")
print("="*80)

# 1. Posterior distributions for Full Model
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle("Full Model - Key Parameters (Individual Data)", fontsize=14, fontweight='bold')

az.plot_posterior(traces['Full_Model'], var_names=["beta_xi"], 
                  hdi_prob=0.94, ax=axes[0,0])
axes[0,0].set_title("Protection Exponent β_ξ")

az.plot_posterior(traces['Full_Model'], var_names=["xi_intercept", "xi_slope"], 
                  hdi_prob=0.94, ax=axes[0,1])
axes[0,1].set_title("ξ ~ Viral Load Relationship")

az.plot_posterior(traces['Full_Model'], var_names=["NAA_baseline"], 
                  hdi_prob=0.94, ax=axes[1,0])
axes[1,0].set_title("Baseline NAA")

# WAIC comparison bar plot
ax = axes[1,1]
waic_vals = [waic_results[name].waic for name in models.keys()]
waic_ses = [waic_results[name].waic_se for name in models.keys()]
x_pos = np.arange(len(models))
ax.bar(x_pos, waic_vals, yerr=waic_ses, capsize=5)
ax.set_xticks(x_pos)
ax.set_xticklabels(models.keys(), rotation=45, ha='right')
ax.set_ylabel('WAIC')
ax.set_title('Model Comparison (lower is better)')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(RESULTS_DIR / "individual_level_comparison.png", dpi=300, bbox_inches='tight')
print(f"✓ Saved: {RESULTS_DIR / 'individual_level_comparison.png'}")

# 2. Predictions vs Observations
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle("Predictions vs Observations (Individual Patients)", fontsize=14, fontweight='bold')

for i, (name, trace) in enumerate(traces.items()):
    ax = axes[i]
    post = trace.posterior
    
    # Get posterior predictions
    NAA_pred_post = post['NAA_pred'].values.reshape(-1, n_patients)
    NAA_pred_mean = NAA_pred_post.mean(axis=0)
    
    # Scatter plot
    ax.scatter(NAA_obs, NAA_pred_mean, alpha=0.5, s=30)
    
    # Perfect prediction line
    lim = [NAA_obs.min() - 1, NAA_obs.max() + 1]
    ax.plot(lim, lim, 'k--', alpha=0.5, label='Perfect prediction')
    
    # Calculate R²
    residuals = NAA_obs - NAA_pred_mean
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((NAA_obs - NAA_obs.mean())**2)
    r2 = 1 - (ss_res / ss_tot)
    
    ax.set_xlabel('Observed NAA')
    ax.set_ylabel('Predicted NAA')
    ax.set_title(f"{name}\nR² = {r2:.3f}")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_aspect('equal')

plt.tight_layout()
plt.savefig(RESULTS_DIR / "predictions_scatter.png", dpi=300, bbox_inches='tight')
print(f"✓ Saved: {RESULTS_DIR / 'predictions_scatter.png'}")

# 3. ξ vs Viral Load (Full Model only)
fig, ax = plt.subplots(figsize=(10, 6))

post = traces['Full_Model'].posterior
xi_i_post = post['xi_i'].values.reshape(-1, n_patients)
xi_i_mean = xi_i_post.mean(axis=0)
xi_i_std = xi_i_post.std(axis=0)

# Sort by viral load for plotting
sort_idx = np.argsort(logVL)
ax.plot(logVL[sort_idx], xi_i_mean[sort_idx], 'b-', linewidth=2, label='Mean ξ')
ax.fill_between(logVL[sort_idx], 
                 (xi_i_mean - xi_i_std)[sort_idx],
                 (xi_i_mean + xi_i_std)[sort_idx],
                 alpha=0.3, label='±1 SD')

ax.set_xlabel('log(Viral Load)')
ax.set_ylabel('Noise Correlation Length ξ (nm)')
ax.set_title('ξ vs Viral Load (Individual Patients)')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(RESULTS_DIR / "xi_vs_viral_load.png", dpi=300, bbox_inches='tight')
print(f"✓ Saved: {RESULTS_DIR / 'xi_vs_viral_load.png'}")

print("\n" + "="*80)
print(" ANALYSIS COMPLETE")
print("="*80)
print(f"\nResults saved to: {RESULTS_DIR}/")
print("\nKey files:")
print("  - waic_comparison.txt")
print("  - full_model_summary.csv")
print("  - individual_level_comparison.png")
print("  - predictions_scatter.png")
print("  - xi_vs_viral_load.png")
print("  - *_trace.nc (MCMC chains)")
