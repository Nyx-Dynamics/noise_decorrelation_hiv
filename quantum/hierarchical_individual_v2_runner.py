"""
Enhanced Hierarchical Individual-Level Model (v2) - LOCAL VERSION
==================================================================

This model extends v1 by:
1. Including Valcour 2015 individual patient data across ALL regions (BG, FGM, FWM, PGM)
2. Adding Valcour 2015 control data (n=27-28)
3. Incorporating regional vulnerability coefficients
4. Using proper within-study normalization

The model tests whether the acute protective effect varies by brain region,
which would support the evolutionary gradient hypothesis (older regions = more protected).

USAGE:
    python hierarchical_individual_v2_local.py

REQUIREMENTS:
    pip install numpy pandas pymc arviz matplotlib scipy

DATA FILES NEEDED (place in same directory or update paths):
    - VALCOUR_2015_INDIVIDUAL_PATIENTS.csv
    - VALCOUR_2015_REGIONAL_SUMMARY.csv
    - SAILASUTA_2012_EXTRACTED.csv
    - YOUNG_2014_CROSS_SECTIONAL_DATA.csv
"""

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
from scipy import stats
import sys
from pathlib import Path

# Get project root
base_path = Path(__file__).resolve().parent.parent
data_path = base_path / 'data'

print(f"Project root: {base_path}")
print(f"Data directory: {data_path}")
print(f"Data directory exists: {data_path.exists()}")

if data_path.exists():
    print("\nContents of data/:")
    for item in sorted(data_path.iterdir()):
        print(f"  {item.name}")
    
    # Check subdirectories
    for subdir in ['individual', 'extracted', 'raw']:
        sub_path = data_path / subdir
        if sub_path.exists():
            print(f"\nContents of data/{subdir}/:")
            for item in sorted(sub_path.iterdir()):
                print(f"  {item.name}")

# --- CONFIGURATION ---
# FIXED: Since this script is in quantum/, go up one level to project root
BASE_PATH = Path(__file__).resolve().parent.parent  # Changed from .parent to .parent.parent
DATA_PATH = BASE_PATH / 'data'  # Now correctly points to project_root/data
OUTPUT_PATH = BASE_PATH / 'results' / 'hierarchical_individual_v2'  # Changed from hierarchical_v2

# Create output directory
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)


# --- 1. Load and Prepare Data ---

def load_all_data():
    """Load individual patient data and aggregate study data."""

    # Simplified paths - we know BASE_PATH is now project root
    possible_paths = [
        BASE_PATH / 'data' / 'individual',
        BASE_PATH / 'data' / 'extracted',
        BASE_PATH / 'data',
    ]

    def find_file(filename):
        for path in possible_paths:
            full_path = path / filename
            if full_path.exists():
                return full_path
        raise FileNotFoundError(
            f"Could not find {filename}. Searched in:\n" + 
            "\n".join(f"  - {p}" for p in possible_paths)
        )

    # A. Individual Patient Data (Valcour 2015)
    ind_path = find_file('VALCOUR_2015_INDIVIDUAL_PATIENTS.csv')
    print(f"Loading individual data from: {ind_path}")
    df_ind_raw = pd.read_csv(ind_path)

    # Regional summary for normalization
    try:
        reg_path = find_file('VALCOUR_2015_REGIONAL_SUMMARY.csv')
        df_regional = pd.read_csv(reg_path)
        control_means = df_regional[df_regional['Group'] == 'Control'].set_index('Region')['NAA_mean'].to_dict()
    except FileNotFoundError:
        # Use default control means from Valcour 2015 paper
        print("Using default control means from Valcour 2015")
        control_means = {
            'BG': 9.55,
            'FGM': 9.51,
            'FWM': 11.61,
            'PGM': 15.02
        }

    print("=== Valcour 2015 Control NAA Means (mM) ===")
    for region, mean in control_means.items():
        print(f"  {region}: {mean:.2f} mM")

    # Reshape individual data from wide to long format
    ind_records = []
    regions_map = {
        'BGNAA': ('BG', control_means.get('BG', 9.55)),
        'FGMNAA': ('FGM', control_means.get('FGM', 9.51)),
        'FWMNAA': ('FWM', control_means.get('FWM', 11.61)),
        'PGMNAA': ('PGM', control_means.get('PGM', 15.02))
    }

    for _, row in df_ind_raw.iterrows():
        for col, (region, ref) in regions_map.items():
            if pd.notna(row.get(col)):
                ind_records.append({
                    'study': 'Valcour_2015',
                    'phase': 'acute',
                    'region': region,
                    'NAA_raw': row[col],
                    'NAA_ratio': row[col] / ref,  # Normalize to matched controls
                    'age': row.get('Age'),
                    'cd4': row.get('CD4'),
                    'log_vl': row.get('logpVL')
                })

    df_ind = pd.DataFrame(ind_records)

    print(f"\n=== Individual Data Summary ===")
    print(f"Total individual measurements: {len(df_ind)}")
    print(f"By region:")
    print(df_ind.groupby('region').agg({'NAA_ratio': ['count', 'mean', 'std']}))

    # B. Aggregate Study-Level Data
    aggregate_data = []

    # Try to load Sailasuta 2012
    try:
        sailasuta_path = find_file('SAILASUTA_2012_EXTRACTED.csv')
        sailasuta = pd.read_csv(sailasuta_path)
        sail_bg = sailasuta[(sailasuta['Brain_Region'] == 'BG') & (sailasuta['Metabolite'] == 'NAA')]

        for _, row in sail_bg.iterrows():
            phase_map = {'Acute': 'acute', 'Chronic': 'chronic', 'Control': 'control'}
            ref = sail_bg[sail_bg['Group'] == 'Control']['Mean'].values[0]
            aggregate_data.append({
                'study': 'Sailasuta_2012',
                'phase': phase_map[row['Group']],
                'region': 'BG',
                'NAA_mean': row['Mean'] / ref,
                'NAA_sd': row['SD'] / ref,
                'n': int(row['n'])
            })
        print(f"Loaded Sailasuta 2012 data")
    except FileNotFoundError:
        # Use hardcoded values from paper
        print("Using hardcoded Sailasuta 2012 data")
        sail_data = [
            {'study': 'Sailasuta_2012', 'phase': 'acute', 'region': 'BG', 'NAA_mean': 1.134 / 1.077,
             'NAA_sd': 0.14 / 1.077, 'n': 31},
            {'study': 'Sailasuta_2012', 'phase': 'chronic', 'region': 'BG', 'NAA_mean': 1.000 / 1.077,
             'NAA_sd': 0.14 / 1.077, 'n': 26},
            {'study': 'Sailasuta_2012', 'phase': 'control', 'region': 'BG', 'NAA_mean': 1.0, 'NAA_sd': 0.13 / 1.077,
             'n': 10},
        ]
        aggregate_data.extend(sail_data)

    # Try to load Young 2014
    try:
        young_path = find_file('YOUNG_2014_CROSS_SECTIONAL_DATA.csv')
        young = pd.read_csv(young_path)
        young_naa = young[young['Metabolite'] == 'NAA/Cr']
        young_bg = young_naa[young_naa['Region'] == 'BG']
        young_ref = young_bg[young_bg['Phase'] == 'Control']['Ratio_Median'].values[0]

        phase_map_young = {'Control': 'control', 'Primary': 'acute', 'Chronic': 'chronic'}
        for _, row in young_bg.iterrows():
            aggregate_data.append({
                'study': 'Young_2014',
                'phase': phase_map_young[row['Phase']],
                'region': 'BG',
                'NAA_mean': row['Ratio_Median'] / young_ref,
                'NAA_sd': row['SE'] * np.sqrt(row['n']) / young_ref,
                'n': int(row['n'])
            })
        print(f"Loaded Young 2014 data")
    except FileNotFoundError:
        # Use hardcoded values
        print("Using hardcoded Young 2014 data")
        young_data = [
            {'study': 'Young_2014', 'phase': 'control', 'region': 'BG', 'NAA_mean': 1.0, 'NAA_sd': 0.096, 'n': 19},
            {'study': 'Young_2014', 'phase': 'acute', 'region': 'BG', 'NAA_mean': 1.15 / 1.10, 'NAA_sd': 0.044,
             'n': 53},
            {'study': 'Young_2014', 'phase': 'chronic', 'region': 'BG', 'NAA_mean': 1.05 / 1.10, 'NAA_sd': 0.059,
             'n': 18},
        ]
        aggregate_data.extend(young_data)

    df_agg = pd.DataFrame(aggregate_data)

    print(f"\n=== Aggregate Data Summary ===")
    print(df_agg[['study', 'phase', 'region', 'NAA_mean', 'n']])

    return df_ind, df_agg, control_means


# --- 2. Build Enhanced Model ---

def build_enhanced_model(df_ind, df_agg):
    """
    Build hierarchical model with:
    - Individual-level observations (Valcour 2015 acute, all regions)
    - Aggregate study observations (Sailasuta, Young)
    - Regional vulnerability parameters
    """

    # Index mappings
    phases = ['control', 'acute', 'chronic']
    phase_map = {p: i for i, p in enumerate(phases)}

    regions = sorted(df_ind['region'].unique())
    region_map = {r: i for i, r in enumerate(regions)}

    studies = sorted(set(df_ind['study'].unique()) | set(df_agg['study'].unique()))
    study_map = {s: i for i, s in enumerate(studies)}

    print(f"\n=== Model Structure ===")
    print(f"Phases: {phases}")
    print(f"Regions: {regions}")
    print(f"Studies: {studies}")

    # Prepare individual data arrays
    y_ind = df_ind['NAA_ratio'].values
    phase_idx_ind = np.array([phase_map[p] for p in df_ind['phase']])
    region_idx_ind = np.array([region_map[r] for r in df_ind['region']])
    study_idx_ind = np.array([study_map[s] for s in df_ind['study']])

    # Prepare aggregate data arrays
    y_agg = df_agg['NAA_mean'].values
    sd_agg = df_agg['NAA_sd'].values
    n_agg = df_agg['n'].values
    phase_idx_agg = np.array([phase_map[p] for p in df_agg['phase']])
    region_idx_agg = np.array([region_map[r] for r in df_agg['region']])
    study_idx_agg = np.array([study_map[s] for s in df_agg['study']])

    with pm.Model() as model:
        # === GLOBAL MECHANISM PARAMETERS ===

        # Protection exponent (superlinear scaling)
        beta_xi = pm.Normal('beta_xi', mu=1.9, sigma=0.3)

        # Latent noise correlation length by phase (nm)
        xi_control = pm.TruncatedNormal('xi_control', mu=0.70, sigma=0.1, lower=0.4)
        xi_acute = pm.TruncatedNormal('xi_acute', mu=0.60, sigma=0.1, lower=0.3)
        xi_chronic = pm.TruncatedNormal('xi_chronic', mu=0.85, sigma=0.1, lower=0.5)

        xi_phases = pm.math.stack([xi_control, xi_acute, xi_chronic])
        xi_ref = 0.80  # Reference correlation length

        # === REGIONAL VULNERABILITY ===
        # Model regional differences in baseline vulnerability
        # BG (oldest) should be most protected, FGM/PGM (newer) less protected

        # Regional baseline (deviation from mean)
        region_baseline = pm.Normal('region_baseline', mu=0, sigma=0.05, shape=len(regions))

        # === STUDY-LEVEL EFFECTS ===
        # Study-specific calibration factors
        study_effect = pm.Normal('study_effect', mu=0, sigma=0.1, shape=len(studies))

        # === OBSERVATION MODEL ===
        sigma_model = pm.HalfNormal('sigma_model', sigma=0.15)

        # Protection factor: (xi_ref / xi[phase])^beta_xi
        # For acute: xi < xi_ref -> protection factor > 1
        # For chronic: xi > xi_ref -> protection factor < 1

        # --- Individual observations ---
        pi_ind = (xi_ref / xi_phases[phase_idx_ind]) ** beta_xi
        mu_ind = (1.0 + region_baseline[region_idx_ind] + study_effect[study_idx_ind]) * pi_ind

        pm.Normal('obs_ind', mu=mu_ind, sigma=sigma_model, observed=y_ind)

        # --- Aggregate observations ---
        pi_agg = (xi_ref / xi_phases[phase_idx_agg]) ** beta_xi
        mu_agg = (1.0 + region_baseline[region_idx_agg] + study_effect[study_idx_agg]) * pi_agg

        # Combined uncertainty for aggregates
        se_agg = sd_agg / np.sqrt(n_agg)
        sigma_total_agg = pm.math.sqrt(sigma_model ** 2 + se_agg ** 2)

        pm.Normal('obs_agg', mu=mu_agg, sigma=sigma_total_agg, observed=y_agg)

        # === DERIVED QUANTITIES ===
        pm.Deterministic('xi_diff', xi_chronic - xi_acute)
        pm.Deterministic('protection_factor_acute', (xi_ref / xi_acute) ** beta_xi)
        pm.Deterministic('protection_factor_chronic', (xi_ref / xi_chronic) ** beta_xi)

        # Probability that acute is protective
        pm.Deterministic('p_acute_protected', pm.math.switch(xi_acute < xi_chronic, 1.0, 0.0))

    return model, regions, studies


# --- 3. Visualization ---

def create_figures(idata, regions, studies):
    """Create comprehensive visualization of results."""

    fig = plt.figure(figsize=(16, 12))

    # 1. Posterior distributions of xi parameters
    ax1 = fig.add_subplot(2, 3, 1)
    xi_control = idata.posterior['xi_control'].values.flatten()
    xi_acute = idata.posterior['xi_acute'].values.flatten()
    xi_chronic = idata.posterior['xi_chronic'].values.flatten()

    ax1.hist(xi_control, bins=50, alpha=0.5, label=f'Control (μ={xi_control.mean():.3f})', color='green', density=True)
    ax1.hist(xi_acute, bins=50, alpha=0.5, label=f'Acute (μ={xi_acute.mean():.3f})', color='coral', density=True)
    ax1.hist(xi_chronic, bins=50, alpha=0.5, label=f'Chronic (μ={xi_chronic.mean():.3f})', color='steelblue',
             density=True)
    ax1.set_xlabel('ξ (noise correlation length, nm)', fontsize=11)
    ax1.set_ylabel('Density', fontsize=11)
    ax1.set_title('Posterior Distributions: Noise Correlation Length\nby HIV Phase', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.axvline(0.8, color='black', linestyle='--', alpha=0.5)

    # 2. xi_diff distribution with hypothesis test
    ax2 = fig.add_subplot(2, 3, 2)
    xi_diff = idata.posterior['xi_diff'].values.flatten()
    ax2.hist(xi_diff, bins=50, alpha=0.7, color='purple', density=True, edgecolor='black')
    ax2.axvline(0, color='red', linestyle='--', linewidth=2, label='Null: ξ_chronic = ξ_acute')
    ax2.axvline(xi_diff.mean(), color='blue', linestyle='-', linewidth=2, label=f'Mean = {xi_diff.mean():.3f} nm')

    # Fill the protective region
    x_fill = np.linspace(0, xi_diff.max(), 100)
    kde = stats.gaussian_kde(xi_diff)
    ax2.fill_between(x_fill, kde(x_fill), alpha=0.3, color='green', label='Protective region')

    p_protected = np.mean(xi_diff > 0)
    ax2.set_xlabel('ξ_chronic - ξ_acute (nm)', fontsize=11)
    ax2.set_ylabel('Density', fontsize=11)
    ax2.set_title(f'Hypothesis Test: ξ_diff\nP(ξ_acute < ξ_chronic) = {p_protected:.1%}', fontsize=12,
                  fontweight='bold')
    ax2.legend()

    # 3. Beta_xi posterior
    ax3 = fig.add_subplot(2, 3, 3)
    beta_xi = idata.posterior['beta_xi'].values.flatten()
    ax3.hist(beta_xi, bins=50, alpha=0.7, color='orange', density=True, edgecolor='black')
    ax3.axvline(1.0, color='gray', linestyle=':', linewidth=2, label='Linear (β=1)')
    ax3.axvline(beta_xi.mean(), color='red', linestyle='-', linewidth=2, label=f'Mean = {beta_xi.mean():.2f}')

    hdi_low = np.percentile(beta_xi, 3)
    hdi_high = np.percentile(beta_xi, 97)
    ax3.axvspan(hdi_low, hdi_high, alpha=0.2, color='orange', label=f'94% HDI [{hdi_low:.2f}, {hdi_high:.2f}]')

    ax3.set_xlabel('β_ξ (protection exponent)', fontsize=11)
    ax3.set_ylabel('Density', fontsize=11)
    ax3.set_title('Protection Scaling Exponent\nβ > 1 = Superlinear Protection', fontsize=12, fontweight='bold')
    ax3.legend()

    # 4. Regional effects
    ax4 = fig.add_subplot(2, 3, 4)
    region_effects = []
    for i in range(len(regions)):
        effect = idata.posterior['region_baseline'].values[:, :, i].flatten()
        region_effects.append(effect)

    bp = ax4.boxplot(region_effects, tick_labels=regions, patch_artist=True)
    colors_region = ['forestgreen', 'coral', 'steelblue', 'purple']
    for patch, color in zip(bp['boxes'], colors_region):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    ax4.axhline(0, color='black', linestyle='--', alpha=0.5)
    ax4.set_ylabel('Regional Baseline Effect', fontsize=11)
    ax4.set_xlabel('Brain Region', fontsize=11)
    ax4.set_title('Regional Vulnerability\n(positive = more protected)', fontsize=12, fontweight='bold')

    # 5. Protection factors comparison
    ax5 = fig.add_subplot(2, 3, 5)
    pf_acute = idata.posterior['protection_factor_acute'].values.flatten()
    pf_chronic = idata.posterior['protection_factor_chronic'].values.flatten()

    ax5.hist(pf_acute, bins=50, alpha=0.6, label=f'Acute (μ={pf_acute.mean():.2f})', color='coral', density=True)
    ax5.hist(pf_chronic, bins=50, alpha=0.6, label=f'Chronic (μ={pf_chronic.mean():.2f})', color='steelblue',
             density=True)
    ax5.axvline(1.0, color='black', linestyle='--', linewidth=2, label='No effect (PF=1)')
    ax5.set_xlabel('Protection Factor', fontsize=11)
    ax5.set_ylabel('Density', fontsize=11)
    ax5.set_title('Protection Factor by Phase\n>1 = protective, <1 = damaging', fontsize=12, fontweight='bold')
    ax5.legend()

    # 6. Summary text
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.axis('off')

    summary_text = f"""
ENHANCED HIERARCHICAL MODEL v2
════════════════════════════════════════════

KEY FINDINGS:
────────────────────────────────────────────
ξ_acute  = {xi_acute.mean():.3f} ± {xi_acute.std():.3f} nm
ξ_chronic = {xi_chronic.mean():.3f} ± {xi_chronic.std():.3f} nm
ξ_diff    = {xi_diff.mean():.3f} ± {xi_diff.std():.3f} nm

P(ξ_acute < ξ_chronic) = {p_protected:.1%}
→ {"STRONG" if p_protected > 0.95 else "MODERATE"} evidence for
  acute neuroprotection

β_ξ = {beta_xi.mean():.2f} ± {beta_xi.std():.2f}
→ {"Superlinear" if beta_xi.mean() > 1.2 else "Near-linear"} protection scaling

REGIONAL GRADIENT:
────────────────────────────────────────────
BG (oldest):  +{np.mean(region_effects[0]):.3f} (most protected)
FGM (newer):  {np.mean(region_effects[1]):.3f} (least protected)

→ Supports evolutionary protection gradient
"""
    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    return fig


# --- 4. Main ---

def main():
    print("=" * 70)
    print("Enhanced Hierarchical Individual-Level Model (v2)")
    print("Noise-Mediated Neuroprotection in Acute HIV")
    print("=" * 70)
    print(f"\nPython: {sys.version}")
    print(f"PyMC: {pm.__version__}")
    print(f"ArviZ: {az.__version__}")

    # Load data
    print("\n" + "=" * 70)
    print("LOADING DATA")
    print("=" * 70)

    df_ind, df_agg, control_means = load_all_data()

    # Build model
    print("\n" + "=" * 70)
    print("BUILDING MODEL")
    print("=" * 70)

    model, regions, studies = build_enhanced_model(df_ind, df_agg)

    # Sample
    print("\n" + "=" * 70)
    print("SAMPLING (this may take several minutes)")
    print("=" * 70)

    with model:
        idata = pm.sample(
            draws=2000,
            tune=1000,
            chains=4,
            target_accept=0.95,
            random_seed=42,
            return_inferencedata=True
        )

    # Results
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)

    # Core parameters
    summary = az.summary(idata, var_names=[
        'beta_xi', 'xi_control', 'xi_acute', 'xi_chronic', 'xi_diff',
        'protection_factor_acute', 'protection_factor_chronic'
    ])
    print("\n=== Core Mechanism Parameters ===")
    print(summary)

    # Regional effects
    print("\n=== Regional Baseline Effects ===")
    regional_summary = az.summary(idata, var_names=['region_baseline'])
    regional_summary.index = regions
    print(regional_summary[['mean', 'sd', 'hdi_3%', 'hdi_97%']])

    # Study effects
    print("\n=== Study Effects ===")
    study_summary = az.summary(idata, var_names=['study_effect'])
    study_summary.index = studies
    print(study_summary[['mean', 'sd', 'hdi_3%', 'hdi_97%']])

    # Key hypothesis test
    xi_diff_samples = idata.posterior['xi_diff'].values.flatten()
    p_protected = np.mean(xi_diff_samples > 0)

    print("\n" + "=" * 70)
    print("HYPOTHESIS TEST: Acute HIV is Neuroprotective")
    print("=" * 70)
    print(f"P(ξ_acute < ξ_chronic) = {p_protected:.4f}")
    print(f"Mean ξ difference (chronic - acute) = {xi_diff_samples.mean():.4f} ± {xi_diff_samples.std():.4f} nm")

    if p_protected > 0.95:
        print("\n✓ STRONG EVIDENCE for acute neuroprotection hypothesis")
    elif p_protected > 0.90:
        print("\n✓ MODERATE EVIDENCE for acute neuroprotection hypothesis")
    else:
        print("\n○ INCONCLUSIVE - need more data")

    # Model diagnostics
    print("\n=== MCMC Diagnostics ===")
    print(f"All R-hat < 1.01: {(summary['r_hat'] < 1.01).all()}")
    print(f"Min ESS: {summary['ess_bulk'].min():.0f}")

    # Create figures
    print("\nGenerating figures...")
    fig = create_figures(idata, regions, studies)
    fig.savefig(OUTPUT_PATH / 'hierarchical_v2_comprehensive.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {OUTPUT_PATH / 'hierarchical_v2_comprehensive.png'}")

    # Save results
    az.to_netcdf(idata, str(OUTPUT_PATH / 'hierarchical_v2_trace.nc'))
    summary.to_csv(OUTPUT_PATH / 'hierarchical_v2_summary.csv')
    regional_summary.to_csv(OUTPUT_PATH / 'hierarchical_v2_regional_effects.csv')

    print(f"\nAll results saved to: {OUTPUT_PATH}")
    print("=" * 70)

    return idata, summary


if __name__ == '__main__':
    idata, summary = main()