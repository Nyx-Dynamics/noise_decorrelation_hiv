import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Output directory
OUTPUT_DIR = Path('results/convergent_evidence')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# EVIDENCE FROM EACH MODEL
# =============================================================================

evidence = {
    # Model v3.6 - Primary hypothesis test (aggregate data)
    'v3.6_Primary': {
        'description': 'Bayesian hierarchical model (aggregate)',
        'n_patients': 143,
        'n_studies': 3,
        'xi_acute_nm': 0.425,
        'xi_acute_sd': 0.065,
        'xi_chronic_nm': 0.790,
        'xi_chronic_sd': 0.065,
        'beta_xi': 2.33,
        'beta_xi_sd': 0.51,
        'P_acute_protected': 0.999,
        'cohens_d': 5.63,
        'prediction_error_pct': 2.0,
        'role': 'PRIMARY',
    },

    # Model v4.0 - Enzyme kinetics validation
    'v4.0_Enzyme': {
        'description': 'Enzyme kinetics mechanistic model',
        'n_patients': 143,
        'n_studies': 3,
        'xi_acute_nm': 0.55,  # From enzyme model
        'xi_acute_sd': 0.08,
        'xi_chronic_nm': 0.82,
        'xi_chronic_sd': 0.10,
        'beta_xi': 1.45,  # |beta| from enzyme model
        'beta_xi_sd': 0.43,
        'P_acute_protected': 0.99,
        'cohens_d': 3.0,
        'prediction_error_pct': 3.0,
        'role': 'VALIDATION (mechanism)',
    },

    # Model v1/v2 - Individual patient validation
    'v2_Individual': {
        'description': 'Hierarchical individual-level model',
        'n_patients': 176,  # Individual measurements
        'n_studies': 3,
        'xi_acute_nm': 0.775,
        'xi_acute_sd': 0.035,
        'xi_chronic_nm': 0.850,
        'xi_chronic_sd': 0.050,
        'beta_xi': 1.75,
        'beta_xi_sd': 0.31,
        'P_acute_protected': 0.924,
        'cohens_d': 1.74,
        'prediction_error_pct': 5.0,
        'role': 'VALIDATION (individual)',
    },
}

# Cross-cohort replication data
cohort_data = {
    'Sailasuta_2012': {
        'location': 'Thailand',
        'acute_n': 31,
        'chronic_n': 26,
        'control_n': 10,
        'acute_naa_ratio': 1.053,  # vs control
        'chronic_naa_ratio': 0.929,
        'p_value': 0.005,
        'finding': 'Acute > Chronic (p=0.005)',
    },
    'Young_2014': {
        'location': 'USA/Uganda',
        'acute_n': 53,
        'chronic_n': 18,
        'control_n': 19,
        'acute_naa_ratio': 1.045,
        'chronic_naa_ratio': 0.955,
        'p_value': 0.003,
        'finding': 'Acute > Chronic (p=0.003)',
    },
    'Valcour_2015': {
        'location': 'Thailand',
        'acute_n': 44,
        'chronic_n': 0,  # No chronic in this study
        'control_n': 28,
        'acute_naa_ratio': 1.090,  # BG region
        'chronic_naa_ratio': None,
        'p_value': None,
        'finding': 'Acute NAA preserved (ratio=1.09)',
    },
}


# =============================================================================
# CREATE SUMMARY TABLE
# =============================================================================

def create_summary_table():
    """Create Table 1: Convergent Evidence Summary"""

    rows = []

    # Model comparison rows
    for model_name, data in evidence.items():
        rows.append({
            'Analysis': model_name.replace('_', ' '),
            'Role': data['role'],
            'N': data['n_patients'],
            'ξ_acute (nm)': f"{data['xi_acute_nm']:.3f} ± {data['xi_acute_sd']:.3f}",
            'ξ_chronic (nm)': f"{data['xi_chronic_nm']:.3f} ± {data['xi_chronic_sd']:.3f}",
            'β_ξ': f"{data['beta_xi']:.2f} ± {data['beta_xi_sd']:.2f}",
            'P(ξ_a < ξ_c)': f"{data['P_acute_protected']:.3f}",
            "Cohen's d": f"{data['cohens_d']:.2f}",
        })

    df = pd.DataFrame(rows)
    return df


def create_cohort_table():
    """Create Table 2: Cross-Cohort Replication"""

    rows = []
    for cohort_name, data in cohort_data.items():
        rows.append({
            'Cohort': cohort_name.replace('_', ' '),
            'Location': data['location'],
            'N (acute)': data['acute_n'],
            'N (chronic)': data['chronic_n'] if data['chronic_n'] else '-',
            'N (control)': data['control_n'],
            'Acute/Control': f"{data['acute_naa_ratio']:.3f}",
            'Chronic/Control': f"{data['chronic_naa_ratio']:.3f}" if data['chronic_naa_ratio'] else '-',
            'Finding': data['finding'],
        })

    df = pd.DataFrame(rows)
    return df


# =============================================================================
# CREATE CONVERGENT EVIDENCE FIGURE
# =============================================================================

def create_convergent_figure():
    """Create publication-quality convergent evidence figure"""

    fig = plt.figure(figsize=(14, 10))

    # Layout: 2x2 grid
    # Top left: ξ estimates across models
    # Top right: β_ξ estimates across models
    # Bottom left: P(acute < chronic) across models
    # Bottom right: Cohort replication summary

    # --- Panel A: ξ estimates ---
    ax1 = fig.add_subplot(2, 2, 1)

    models = list(evidence.keys())
    model_labels = ['v3.6\n(Primary)', 'v4.0\n(Enzyme)', 'v2\n(Individual)']
    x = np.arange(len(models))
    width = 0.35

    xi_acute = [evidence[m]['xi_acute_nm'] for m in models]
    xi_acute_err = [evidence[m]['xi_acute_sd'] for m in models]
    xi_chronic = [evidence[m]['xi_chronic_nm'] for m in models]
    xi_chronic_err = [evidence[m]['xi_chronic_sd'] for m in models]

    bars1 = ax1.bar(x - width / 2, xi_acute, width, yerr=xi_acute_err,
                    label='ξ_acute', color='coral', capsize=5, alpha=0.8)
    bars2 = ax1.bar(x + width / 2, xi_chronic, width, yerr=xi_chronic_err,
                    label='ξ_chronic', color='steelblue', capsize=5, alpha=0.8)

    ax1.set_ylabel('ξ (nm)', fontsize=12)
    ax1.set_title('A. Noise Correlation Length Estimates', fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(model_labels)
    ax1.legend()
    ax1.set_ylim(0, 1.1)

    # Add significance indicator
    for i, m in enumerate(models):
        if evidence[m]['P_acute_protected'] > 0.95:
            ax1.annotate('***', xy=(i, max(xi_chronic[i], xi_acute[i]) + 0.15),
                         ha='center', fontsize=12)

    # --- Panel B: β_ξ estimates ---
    ax2 = fig.add_subplot(2, 2, 2)

    beta_vals = [evidence[m]['beta_xi'] for m in models]
    beta_errs = [evidence[m]['beta_xi_sd'] for m in models]

    colors = ['forestgreen', 'orange', 'purple']
    bars = ax2.bar(x, beta_vals, yerr=beta_errs, capsize=5, color=colors, alpha=0.8)

    ax2.axhline(1.0, color='gray', linestyle='--', linewidth=2, label='Linear (β=1)')
    ax2.set_ylabel('β_ξ (protection exponent)', fontsize=12)
    ax2.set_title('B. Protection Scaling Exponent', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(model_labels)
    ax2.legend()
    ax2.set_ylim(0, 3.5)

    # Annotate superlinear
    for i, (val, err) in enumerate(zip(beta_vals, beta_errs)):
        if val - 1.96 * err > 1.0:
            ax2.annotate('superlinear', xy=(i, val + err + 0.2), ha='center', fontsize=9, style='italic')

    # --- Panel C: Posterior probability ---
    ax3 = fig.add_subplot(2, 2, 3)

    p_vals = [evidence[m]['P_acute_protected'] for m in models]

    bars = ax3.bar(x, p_vals, color=['forestgreen', 'orange', 'purple'], alpha=0.8)
    ax3.axhline(0.95, color='red', linestyle='--', linewidth=2, label='95% threshold')
    ax3.axhline(0.99, color='darkred', linestyle=':', linewidth=2, label='99% threshold')

    ax3.set_ylabel('P(ξ_acute < ξ_chronic)', fontsize=12)
    ax3.set_title('C. Posterior Probability of Acute Protection', fontsize=12, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(model_labels)
    ax3.set_ylim(0.85, 1.02)
    ax3.legend(loc='lower right')

    # Add probability values
    for i, p in enumerate(p_vals):
        ax3.annotate(f'{p:.1%}', xy=(i, p + 0.01), ha='center', fontsize=11, fontweight='bold')

    # --- Panel D: Cross-cohort replication ---
    ax4 = fig.add_subplot(2, 2, 4)

    cohorts = list(cohort_data.keys())
    cohort_labels = [c.replace('_', '\n') for c in cohorts]
    x_cohort = np.arange(len(cohorts))

    acute_ratios = [cohort_data[c]['acute_naa_ratio'] for c in cohorts]
    chronic_ratios = [cohort_data[c]['chronic_naa_ratio'] if cohort_data[c]['chronic_naa_ratio'] else 0 for c in
                      cohorts]

    bars1 = ax4.bar(x_cohort - width / 2, acute_ratios, width, label='Acute/Control',
                    color='coral', alpha=0.8)

    # Only plot chronic where available
    chronic_x = [i for i, c in enumerate(cohorts) if cohort_data[c]['chronic_naa_ratio']]
    chronic_y = [cohort_data[c]['chronic_naa_ratio'] for c in cohorts if cohort_data[c]['chronic_naa_ratio']]
    ax4.bar([x_cohort[i] + width / 2 for i in chronic_x], chronic_y, width,
            label='Chronic/Control', color='steelblue', alpha=0.8)

    ax4.axhline(1.0, color='green', linestyle='--', linewidth=2, label='Healthy baseline')
    ax4.set_ylabel('NAA Ratio (vs Control)', fontsize=12)
    ax4.set_title('D. Cross-Cohort Replication', fontsize=12, fontweight='bold')
    ax4.set_xticks(x_cohort)
    ax4.set_xticklabels(cohort_labels)
    ax4.legend(loc='lower right')
    ax4.set_ylim(0.85, 1.15)

    # Add sample sizes
    for i, c in enumerate(cohorts):
        n_total = cohort_data[c]['acute_n'] + cohort_data[c]['control_n']
        if cohort_data[c]['chronic_n']:
            n_total += cohort_data[c]['chronic_n']
        ax4.annotate(f'n={n_total}', xy=(i, 0.87), ha='center', fontsize=9)

    plt.tight_layout()

    # Add main title
    fig.suptitle('Convergent Evidence for Noise-Mediated Neuroprotection in Acute HIV',
                 fontsize=14, fontweight='bold', y=1.02)

    return fig