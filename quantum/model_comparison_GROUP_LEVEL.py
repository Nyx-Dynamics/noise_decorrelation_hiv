"""
Model Comparison - GROUP LEVEL DATA
====================================

Compares three models using group-level summary statistics (n=3: control, acute, chronic):
1. Full Model: ξ coupling with β_ξ estimated (nonlinear protection)
2. Linear Model: ξ coupling with β_ξ = 1 (linear protection)  
3. No Coupling Model: No ξ effect (null model)

Uses WAIC for model comparison.

Data: Group-level means from CRITICAL_STUDIES_COMPLETE_DATA.csv
      Basal ganglia region (most consistent across studies)

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
RESULTS_DIR = Path("quantum/results/model_comparison_group")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

print("\n" + "="*80)
print(" MODEL COMPARISON - GROUP LEVEL DATA")
print("="*80)

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Load group-level data
data = pd.read_csv("../data/extracted/CRITICAL_STUDIES_COMPLETE_DATA.csv")

# Filter for basal ganglia (most reliable region)
bg_data = data[data['Region'] == 'Basal_Ganglia'].copy()

# Debug: Show what phases we have
print("\nPhases found in Basal Ganglia data:")
print(bg_data['Phase'].value_counts())

# CRITICAL: Use ONLY ratio-based studies (Young2014, Sailasuta2012)
# Chang2002 uses mmol/kg (absolute), can't mix with ratios!
ratio_studies = ['Young2014_Spudich', 'Sailasuta2012']
bg_data = bg_data[bg_data['Study'].isin(ratio_studies)].copy()

print(f"\nUsing ratio-based studies: {ratio_studies}")
print(f"Studies excluded (different units): Chang2002, Sailasuta2016")

# CRITICAL: Filter for only Control, Acute, Chronic phases
bg_data = bg_data[bg_data['Phase'].isin(['Control', 'Acute', 'Chronic'])].copy()

print(f"\nAfter filtering, kept {len(bg_data)} rows with Control/Acute/Chronic")

# Get summary statistics for each phase
summary = bg_data.groupby('Phase').agg({
    'NAA_mean': 'mean',
    'NAA_SE': lambda x: np.sqrt(np.mean(x**2)),  # Pool SEs
    'N': 'sum'
}).reset_index()

# Ensure we have exactly 3 phases in correct order
phase_order = ['Control', 'Acute', 'Chronic']
summary = summary.set_index('Phase').reindex(phase_order).reset_index()

# Extract values for the three conditions
NAA_obs = summary['NAA_mean'].values
NAA_se = summary['NAA_SE'].values
N_obs = summary['N'].values.astype(int)

print("\nObserved Data (Basal Ganglia):")
print(f"  Control: NAA = {NAA_obs[0]:.3f} ± {NAA_se[0]:.3f} (n={N_obs[0]})")
print(f"  Acute:   NAA = {NAA_obs[1]:.3f} ± {NAA_se[1]:.3f} (n={N_obs[1]})")
print(f"  Chronic: NAA = {NAA_obs[2]:.3f} ± {NAA_se[2]:.3f} (n={N_obs[2]})")

# ============================================================================
# MODEL DEFINITIONS
# ============================================================================

def build_full_model():
    """Full model: ξ coupling with estimated β_ξ (nonlinear)"""
    with pm.Model() as model:
        # Priors
        NAA_baseline = pm.Normal("NAA_baseline", mu=1.1, sigma=0.1)
        
        # ξ parameters (non-centered)
        xi_control_raw = pm.Normal("xi_control_raw", mu=0, sigma=1)
        xi_acute_raw = pm.Normal("xi_acute_raw", mu=0, sigma=1)
        xi_chronic_raw = pm.Normal("xi_chronic_raw", mu=0, sigma=1)
        
        xi_scale = pm.HalfNormal("xi_scale", sigma=0.2)
        xi_control = pm.Deterministic("xi_control", 0.7 + xi_control_raw * xi_scale)
        xi_acute = pm.Deterministic("xi_acute", 0.6 + xi_acute_raw * xi_scale)
        xi_chronic = pm.Deterministic("xi_chronic", 0.8 + xi_chronic_raw * xi_scale)
        
        # Protection exponent (key parameter)
        beta_xi = pm.Normal("beta_xi", mu=2.0, sigma=0.5)
        
        # NAA predictions with ξ coupling
        NAA_control = pm.Deterministic("NAA_control", 
            NAA_baseline * (1 / xi_control)**beta_xi)
        NAA_acute = pm.Deterministic("NAA_acute",
            NAA_baseline * (1 / xi_acute)**beta_xi)
        NAA_chronic = pm.Deterministic("NAA_chronic",
            NAA_baseline * (1 / xi_chronic)**beta_xi)
        
        NAA_pred = pm.math.stack([NAA_control, NAA_acute, NAA_chronic])
        
        # Likelihood (group-level: observed means with known SEs)
        sigma_obs = pm.HalfNormal("sigma_obs", sigma=0.05)
        obs_error = pm.math.sqrt(NAA_se**2 + sigma_obs**2)
        
        # Store log-likelihood for WAIC
        pm.Normal("obs", mu=NAA_pred, sigma=obs_error, observed=NAA_obs)
        
    return model

def build_linear_model():
    """Linear model: ξ coupling with β_ξ = 1 (linear protection)"""
    with pm.Model() as model:
        # Priors (same as full model)
        NAA_baseline = pm.Normal("NAA_baseline", mu=1.1, sigma=0.1)
        
        xi_control_raw = pm.Normal("xi_control_raw", mu=0, sigma=1)
        xi_acute_raw = pm.Normal("xi_acute_raw", mu=0, sigma=1)
        xi_chronic_raw = pm.Normal("xi_chronic_raw", mu=0, sigma=1)
        
        xi_scale = pm.HalfNormal("xi_scale", sigma=0.2)
        xi_control = pm.Deterministic("xi_control", 0.7 + xi_control_raw * xi_scale)
        xi_acute = pm.Deterministic("xi_acute", 0.6 + xi_acute_raw * xi_scale)
        xi_chronic = pm.Deterministic("xi_chronic", 0.8 + xi_chronic_raw * xi_scale)
        
        # FIXED: β_ξ = 1 (linear)
        beta_xi = pm.Deterministic("beta_xi", pm.math.constant(1.0))
        
        # NAA predictions (linear coupling)
        NAA_control = pm.Deterministic("NAA_control", 
            NAA_baseline * (1 / xi_control))
        NAA_acute = pm.Deterministic("NAA_acute",
            NAA_baseline * (1 / xi_acute))
        NAA_chronic = pm.Deterministic("NAA_chronic",
            NAA_baseline * (1 / xi_chronic))
        
        NAA_pred = pm.math.stack([NAA_control, NAA_acute, NAA_chronic])
        
        # Likelihood (same as full model)
        sigma_obs = pm.HalfNormal("sigma_obs", sigma=0.05)
        obs_error = pm.math.sqrt(NAA_se**2 + sigma_obs**2)
        
        pm.Normal("obs", mu=NAA_pred, sigma=obs_error, observed=NAA_obs)
        
    return model

def build_no_coupling_model():
    """No coupling model: No ξ effect (null model)"""
    with pm.Model() as model:
        # NAA levels directly estimated (no ξ coupling)
        NAA_control = pm.Normal("NAA_control", mu=1.1, sigma=0.1)
        NAA_acute = pm.Normal("NAA_acute", mu=1.1, sigma=0.1)
        NAA_chronic = pm.Normal("NAA_chronic", mu=1.0, sigma=0.1)
        
        NAA_pred = pm.math.stack([NAA_control, NAA_acute, NAA_chronic])
        
        # Likelihood
        sigma_obs = pm.HalfNormal("sigma_obs", sigma=0.05)
        obs_error = pm.math.sqrt(NAA_se**2 + sigma_obs**2)
        
        pm.Normal("obs", mu=NAA_pred, sigma=obs_error, observed=NAA_obs)
        
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
        
        # Compute log-likelihood for WAIC
        pm.compute_log_likelihood(trace)
        
        traces[name] = trace
        
        # Save trace
        trace.to_netcdf(RESULTS_DIR / f"{name}_trace.nc")
        
        # Convergence check
        rhat = az.rhat(trace)
        max_rhat = float(rhat.to_array().max().values)
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
    waic = az.waic(trace)  # Will automatically find log-likelihood from observed variable
    waic_results[name] = waic
    print(f"\n{name}:")
    print(f"  WAIC: {waic.elpd_waic:.2f}")
    print(f"  SE: {waic.se:.2f}")
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
fig.suptitle("Full Model - Key Parameters", fontsize=14, fontweight='bold')

az.plot_posterior(traces['Full_Model'], var_names=["beta_xi"], 
                  hdi_prob=0.94, ax=axes[0,0])
axes[0,0].set_title("Protection Exponent β_ξ")

az.plot_posterior(traces['Full_Model'], var_names=["xi_acute", "xi_chronic"], 
                  hdi_prob=0.94, ax=axes[0,1])
axes[0,1].set_title("Correlation Lengths ξ")

az.plot_posterior(traces['Full_Model'], var_names=["NAA_control", "NAA_acute", "NAA_chronic"], 
                  hdi_prob=0.94, ax=axes[1,0])
axes[1,0].set_title("NAA Predictions")

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
plt.savefig(RESULTS_DIR / "group_level_comparison.png", dpi=300, bbox_inches='tight')
print(f"✓ Saved: {RESULTS_DIR / 'group_level_comparison.png'}")

# 2. Predictions vs Observations
fig, ax = plt.subplots(figsize=(10, 6))

for i, (name, trace) in enumerate(traces.items()):
    post = trace.posterior
    
    # Get predictions
    NAA_control_post = post['NAA_control'].values.flatten()
    NAA_acute_post = post['NAA_acute'].values.flatten()
    NAA_chronic_post = post['NAA_chronic'].values.flatten()
    
    means = [NAA_control_post.mean(), NAA_acute_post.mean(), NAA_chronic_post.mean()]
    
    offset = (i - 1) * 0.1
    ax.scatter([0 + offset, 1 + offset, 2 + offset], means, 
               label=name, s=100, alpha=0.7)

# Observed data
ax.errorbar([0, 1, 2], NAA_obs, yerr=NAA_se, fmt='ko', 
            markersize=10, capsize=5, label='Observed', linewidth=2)

ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['Control', 'Acute', 'Chronic'])
ax.set_ylabel('NAA (relative units)')
ax.set_title('Model Predictions vs Observed Data')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(RESULTS_DIR / "predictions_vs_obs.png", dpi=300, bbox_inches='tight')
print(f"✓ Saved: {RESULTS_DIR / 'predictions_vs_obs.png'}")

print("\n" + "="*80)
print(" ANALYSIS COMPLETE")
print("="*80)
print(f"\nResults saved to: {RESULTS_DIR}/")
print("\nKey files:")
print("  - waic_comparison.txt")
print("  - full_model_summary.csv")
print("  - group_level_comparison.png")
print("  - predictions_vs_obs.png")
print("  - *_trace.nc (MCMC chains)")
