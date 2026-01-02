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
import pandas as pd
from pathlib import Path
import argparse
from typing import Optional, List, Union, Dict

# Optional heavy dependencies (allow importing this module just to use the loader)
try:
    import pymc as pm  # type: ignore
except Exception:
    pm = None  # type: ignore
try:
    import arviz as az  # type: ignore
except Exception:
    az = None  # type: ignore
try:
    import matplotlib.pyplot as plt  # type: ignore
except Exception:
    plt = None  # type: ignore

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

# Convert to arrays for easier handling (defaults; may be overridden by loader)
CONDITIONS = ['healthy', 'acute', 'chronic']
NAA_OBS = np.array([CLINICAL_DATA[c]['NAA'] for c in CONDITIONS])
CHO_OBS = np.array([CLINICAL_DATA[c]['Cho'] for c in CONDITIONS])


# =============================================================================
# OPTIONAL: GROUP-LEVEL DATA LOADER (larger dataset)
# =============================================================================

def load_group_summary(
    csv_path: Union[str, Path] = Path('data/extracted/CRITICAL_STUDIES_COMPLETE_DATA.csv'),
    region: str = 'Basal_Ganglia',
    studies: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
):
    """
    Load and summarize larger dataset to 3 phase means for NAA and Cho.

    - Filters to Basal Ganglia by default (most reliable across studies)
    - Uses ratio-based studies by default (Young2014_Spudich, Sailasuta2012)
    - Returns per-phase means and pooled SEs along with counts.

    Returns dict with keys: NAA_obs, Cho_obs, NAA_se, Cho_se, N
    (each np.array of length 3 in order [Control, Acute, Chronic])
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Group data CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    df = df[df['Region'] == region].copy()

    if studies is None:
        # Restrict to ratio-based studies by default
        studies = ['Young2014_Spudich', 'Sailasuta2012']
    df = df[df['Study'].isin(studies)].copy()

    if phases is None:
        phases = ['Control', 'Acute', 'Chronic']
    df = df[df['Phase'].isin(phases)].copy()

    # Aggregate means and pooled SE for each phase
    def _pool_se(x: pd.Series) -> float:
        # Conservative: root-mean of squared SEs
        x = pd.to_numeric(x, errors='coerce').dropna()
        if len(x) == 0:
            return np.nan
        return float(np.sqrt(np.mean(np.square(x))))

    agg_spec = {
        'NAA_mean': 'mean',
        'NAA_SE': _pool_se,
        'N': 'sum'
    }
    has_cho_cols = 'Cho_mean' in df.columns and 'Cho_SE' in df.columns
    if has_cho_cols:
        agg_spec.update({'Cho_mean': 'mean', 'Cho_SE': _pool_se})
    grp = df.groupby('Phase').agg(agg_spec)

    # Ensure order
    grp = grp.reindex(['Control', 'Acute', 'Chronic'])

    NAA_obs = grp['NAA_mean'].to_numpy(dtype=float)
    Cho_obs = (grp['Cho_mean'].to_numpy(dtype=float)
               if 'Cho_mean' in grp.columns else None)
    NAA_se = grp['NAA_SE'].to_numpy(dtype=float)
    Cho_se = (grp['Cho_SE'].to_numpy(dtype=float)
              if 'Cho_SE' in grp.columns else None)
    N = grp['N'].fillna(0).to_numpy(dtype=int)

    # Basic sanity: all three phases present
    if len(NAA_obs) != 3 or np.any(~np.isfinite(NAA_obs)):
        raise RuntimeError("Failed to build 3-phase NAA summary from group data.")

    return {
        'NAA_obs': NAA_obs,
        'Cho_obs': Cho_obs,
        'NAA_se': NAA_se,
        'Cho_se': Cho_se,
        'N': N,
        'phases': ['Control', 'Acute', 'Chronic']
    }


# =============================================================================
# FORWARD MODEL: ENZYME KINETICS
# =============================================================================

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

def build_enzyme_model(observed_data: Optional[Dict] = None):
    """
    Build PyMC model with enzyme kinetics.
    
    KEY PARAMETERS:
    - xi_acute, xi_chronic: Noise correlation lengths
    - beta_xi: Protection factor exponent (expect ~2)
    - gamma_coh: Coherence coupling exponent
    - viral_damage_*: Direct damage to enzymes
    - membrane_*: Membrane turnover rates
    """
    if pm is None:
        raise ImportError("pymc is required to build the model; install pymc>=5.")
    
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
        
        NAA_pred, Cho_pred = forward_model_enzyme(
            xi_acute=xi_acute,
            xi_chronic=xi_chronic,
            beta_xi=beta_xi,
            gamma_coh=gamma_coh,
            viral_damage_acute=viral_damage_acute,
            viral_damage_chronic=viral_damage_chronic,
            membrane_acute=membrane_acute,
            membrane_chronic=membrane_chronic
        )
        
        # =====================================================================
        # LIKELIHOOD
        # =====================================================================
        
        if observed_data is None:
            # Default: use clinical triplet with free noise scales
            sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=0.06)
            sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=0.03)
            NAA_obs_vec = NAA_OBS
            Cho_obs_vec = CHO_OBS
        else:
            # Use provided group-level summary; set weakly-informative noise around pooled SE
            NAA_obs_vec = np.asarray(observed_data.get('NAA_obs'), dtype=float)
            Cho_obs_vec = np.asarray(observed_data.get('Cho_obs'), dtype=float)
            # Use mean(SE) as scale of HalfNormal prior; fallback to defaults if NaN
            naa_se = observed_data.get('NAA_se')
            cho_se = observed_data.get('Cho_se')
            def _safe_scale(arr, default):
                try:
                    arr = np.asarray(arr, dtype=float)
                    m = float(np.nanmean(arr))
                    if np.isfinite(m) and m > 0:
                        return m
                except Exception:
                    pass
                return default
            sigma_NAA = pm.HalfNormal('sigma_NAA', sigma=_safe_scale(naa_se, 0.06))
            sigma_Cho = pm.HalfNormal('sigma_Cho', sigma=_safe_scale(cho_se, 0.03))

        # Likelihood terms (vector of length 3)
        NAA_likelihood = pm.Normal(
            'NAA_obs',
            mu=NAA_pred,
            sigma=sigma_NAA,
            observed=NAA_obs_vec
        )
        
        Cho_likelihood = pm.Normal(
            'Cho_obs',
            mu=Cho_pred,
            sigma=sigma_Cho,
            observed=Cho_obs_vec
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

def run_inference(n_samples=2000, n_chains=4, target_accept=0.99, observed_data: Optional[Dict] = None):
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
    if pm is None or az is None:
        raise ImportError("pymc and arviz are required to run inference; please install them in your environment.")
    print("Building model...")
    model = build_enzyme_model(observed_data=observed_data)
    
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
            random_seed=42,
            idata_kwargs={"log_likelihood": True}
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
    
    # Determine observed vectors from idata if available (supports larger dataset)
    obs = getattr(idata, 'observed_data', None)
    if obs is not None and 'NAA_obs' in obs:
        naa_obs_vec = obs['NAA_obs'].values
        # Flatten to 1D (length 3 expected)
        naa_obs_vec = np.array(naa_obs_vec).ravel()
    else:
        naa_obs_vec = NAA_OBS
    if obs is not None and 'Cho_obs' in obs:
        cho_obs_vec = np.array(obs['Cho_obs'].values).ravel()
    else:
        cho_obs_vec = None

    # Prediction errors
    ppc = idata.posterior_predictive
    NAA_pred_mean = ppc['NAA_obs'].mean(dim=['chain', 'draw']).values
    NAA_pred_mean = np.array(NAA_pred_mean).ravel()
    NAA_errors = 100 * (NAA_pred_mean - naa_obs_vec) / naa_obs_vec

    if 'Cho_obs' in ppc:
        Cho_pred_mean = np.array(ppc['Cho_obs'].mean(dim=['chain', 'draw']).values).ravel()
    else:
        Cho_pred_mean = None

    print("Prediction Errors:")
    if cho_obs_vec is not None and Cho_pred_mean is not None:
        print(f"  {'Phase':12s}  {'NAA Error':>12s}  {'Cho Error':>12s}")
        print("  " + "-" * 40)
        phases = ['Control', 'Acute', 'Chronic']
        Cho_errors = 100 * (Cho_pred_mean - cho_obs_vec) / cho_obs_vec
        for i, ph in enumerate(phases):
            print(f"  {ph:12s}  {NAA_errors[i]:>11.2f}%  {Cho_errors[i]:>11.2f}%")
        print()
    else:
        print(f"  {'Phase':12s}  {'NAA Error':>12s}")
        print("  " + "-" * 26)
        phases = ['Control', 'Acute', 'Chronic']
        for i, ph in enumerate(phases):
            print(f"  {ph:12s}  {NAA_errors[i]:>11.2f}%")
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
    if cho_obs_vec is not None and Cho_pred_mean is not None:
        Cho_errors = 100 * (Cho_pred_mean - cho_obs_vec) / cho_obs_vec
        rms_error = np.sqrt(np.mean(NAA_errors**2 + Cho_errors**2))
    else:
        rms_error = np.sqrt(np.mean(NAA_errors**2))
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
    
    if plt is None:
        raise ImportError("matplotlib is required for plotting; please install matplotlib.")
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

def main(argv: Optional[List[str]] = None):
    """Run complete enzyme model inference.

    CLI options:
      --use-group-data       Use larger dataset (group summary) instead of fixed clinical triplet
      --data-csv PATH        Path to CRITICAL_STUDIES_COMPLETE_DATA.csv
      --region NAME          Brain region to filter (default: Basal_Ganglia)
      --studies LIST         Comma-separated study names to include (default: Young2014_Spudich,Sailasuta2012)
      --samples N            Samples per chain (default: 2000)
      --chains N             Number of chains (default: 4)
      --target-accept X      NUTS target accept (default: 0.99)
    """

    ap = argparse.ArgumentParser(description='v4 Enzyme model with optional larger dataset ingestion')
    ap.add_argument('--use-group-data', action='store_true', help='Use group-level larger dataset (BG ratios)')
    ap.add_argument('--data-csv', type=str, default=str(Path('data/extracted/CRITICAL_STUDIES_COMPLETE_DATA.csv')))
    ap.add_argument('--region', type=str, default='Basal_Ganglia')
    ap.add_argument('--studies', type=str, default='Young2014_Spudich,Sailasuta2012')
    ap.add_argument('--samples', type=int, default=2000)
    ap.add_argument('--chains', type=int, default=4)
    ap.add_argument('--target-accept', type=float, default=0.99)
    args = ap.parse_args(argv)

    observed_data = None
    if args.use_group_data:
        studies = [s.strip() for s in args.studies.split(',') if s.strip()]
        print("Loading group-level data:")
        print(f"  CSV:     {args.data_csv}")
        print(f"  Region:  {args.region}")
        print(f"  Studies: {studies}")
        observed_data = load_group_summary(args.data_csv, region=args.region, studies=studies)
        # Report summary
        print("Group summary (BG ratios):")
        for i, ph in enumerate(['Control', 'Acute', 'Chronic']):
            print(f"  {ph:<7}: NAA={observed_data['NAA_obs'][i]:.3f} (SE≈{observed_data['NAA_se'][i]:.3f}), "
                  f"Cho={observed_data['Cho_obs'][i]:.3f} (SE≈{observed_data['Cho_se'][i]:.3f}), N={observed_data['N'][i]}")

    # Run inference
    idata = run_inference(n_samples=args.samples, n_chains=args.chains, target_accept=args.target_accept, observed_data=observed_data)
    
    # Analyze results
    summary = analyze_results(idata)
    
    # Save results
    output_dir = Path('results/enzyme_v4')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Tag outputs if group data used
    suffix = '_group' if observed_data is not None else ''
    idata.to_netcdf(output_dir / f'trace_v4{suffix}.nc')
    summary.to_csv(output_dir / f'summary_v4{suffix}.csv')
    
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
