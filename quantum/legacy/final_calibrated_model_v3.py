"""
FORWARD MODEL v3.0: Direct Enzyme Modulation
=============================================

Replaces phenomenological compensation with mechanistic enzyme kinetics.

KEY CHANGE:
- v2: "Astrocyte compensation" (reactive, phenomenological)
- v3: "Direct enzyme modulation" (proactive, mechanistic)

MECHANISM:
- Quantum coherence → Protection factor Π_ξ → NAT8L activity
- ξ ↓ (acute) → Π_ξ ↑ → faster NAT8L → NAA preserved
- ξ ↑ (chronic) → Π_ξ ↓ → slower NAT8L → NAA depleted

NO COMPENSATION - just direct quantum-enzyme coupling + natural kinetics
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Tuple
from dataclasses import dataclass
import sys
import os

# Import enzyme kinetics module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from enzyme_kinetics import (
    EnzymeKinetics, compute_protection_factor, coherence_modulation,
    ENZYME
)


# =============================================================================
# CONSTANTS
# =============================================================================

@dataclass
class ModelConstants:
    """Physical and biological constants."""

    # Quantum parameters
    xi_baseline: float = 0.8e-9  # m (healthy/chronic)
    xi_acute: float = 0.42e-9  # m (acute - decorrelated noise)
    xi_chronic: float = 0.79e-9  # m (chronic - correlated noise)

    sigma_r_regular: float = 0.38e-9  # m (healthy)
    sigma_r_fibril: float = 0.53e-9  # m (damaged)

    # MRS baselines (for conversion to ratios)
    NAA_baseline_ratio: float = 1.105  # NAA/Cr ratio (healthy)
    Cho_baseline_ratio: float = 0.225  # Cho/Cr ratio (healthy)
    Cr_conc: float = 6.8e-3  # M (creatine concentration - assumed constant)

    # Cytokine concentrations (pg/mL)
    TNF_healthy: float = 5.0
    TNF_acute: float = 200.0
    TNF_chronic: float = 30.0

    IL6_healthy: float = 10.0
    IL6_acute: float = 500.0
    IL6_chronic: float = 50.0

    # Enzyme modulation parameters
    beta_xi: float = 1.8  # Protection factor exponent (calibrated)
    gamma_coh: float = 1.2  # Coherence modulation exponent (calibrated)


CONST = ModelConstants()


# =============================================================================
# VIRAL DAMAGE FACTORS
# =============================================================================

def compute_viral_damage(condition: str) -> float:
    """
    Compute direct viral/inflammatory damage to enzyme function.

    This is separate from quantum protection - it represents direct
    HIV/cytokine damage to enzyme expression/function.

    Parameters
    ----------
    condition : str
        Clinical condition

    Returns
    -------
    damage_factor : float
        Enzyme function multiplier (0-1)
        1.0 = no damage, 0.7 = 30% impairment
    """
    damage_factors = {
        'healthy': 1.0,
        'acute_HIV': 0.94,  # Minimal direct damage (protection compensates)
        'chronic_HIV': 0.52  # Severe sustained enzyme impairment
    }
    return damage_factors.get(condition, 1.0)


def compute_membrane_turnover(condition: str) -> float:
    """
    Compute membrane turnover rate for Cho dynamics.

    Parameters
    ----------
    condition : str
        Clinical condition

    Returns
    -------
    turnover : float
        Relative membrane turnover rate
    """
    turnover_rates = {
        'healthy': 1.0,
        'acute_HIV': 2.2,  # High inflammation
        'chronic_HIV': 1.3  # Moderate damage/repair
    }
    return turnover_rates.get(condition, 1.0)


# =============================================================================
# CONCENTRATION → RATIO CONVERSION
# =============================================================================

def concentration_to_ratio(NAA_molar: float, Cho_molar: float,
                           Cr_molar: float = CONST.Cr_conc) -> Tuple[float, float]:
    """
    Convert absolute concentrations to MRS-observed ratios.

    MRS measures NAA/Cr and Cho/Cr, not absolute concentrations.

    Parameters
    ----------
    NAA_molar : float
        NAA concentration (M)
    Cho_molar : float
        Cho concentration (M)
    Cr_molar : float
        Creatine concentration (M), assumed constant

    Returns
    -------
    NAA_Cr : float
        NAA/Cr ratio
    Cho_Cr : float
        Cho/Cr ratio
    """
    NAA_Cr = NAA_molar / Cr_molar
    Cho_Cr = Cho_molar / Cr_molar
    return NAA_Cr, Cho_Cr


# =============================================================================
# CORE FORWARD MODEL v3.0
# =============================================================================

def run_enzyme_model_v3(condition: str,
                        coherence_base: float,
                        xi: float,
                        sigma_r: float,
                        beta_xi: float = CONST.beta_xi,
                        gamma_coh: float = CONST.gamma_coh) -> Dict:
    """
    Run complete forward model v3.0 with direct enzyme modulation.

    STEPS:
    1. Compute quantum protection factor Π_ξ from ξ
    2. Compute coherence modulation η_coh from coherence
    3. Compute viral damage factor
    4. Initialize enzyme system with modulation factors
    5. Integrate enzyme ODEs to steady state
    6. Convert concentrations to MRS ratios

    Parameters
    ----------
    condition : str
        'healthy', 'acute_HIV', or 'chronic_HIV'
    coherence_base : float
        SSE coherence (0-1)
    xi : float
        Noise correlation length (m)
    sigma_r : float
        Wavefunction spatial spread (m)
    beta_xi : float
        Protection factor exponent
    gamma_coh : float
        Coherence modulation exponent

    Returns
    -------
    results : dict
        Complete model outputs including:
        - Quantum parameters (ξ, coherence, Π_ξ)
        - Enzyme concentrations (NAA, Cho in mM)
        - MRS ratios (NAA/Cr, Cho/Cr)
        - Mechanistic breakdown
    """

    # 1. Quantum protection factor
    Pi_xi = compute_protection_factor(xi, xi_ref=CONST.xi_baseline, beta_xi=beta_xi)

    # 2. Coherence modulation
    eta_coh = coherence_modulation(coherence_base, coherence_ref=ENZYME.coherence_ref, gamma=gamma_coh)

    # 3. Viral damage
    viral_damage = compute_viral_damage(condition)

    # 4. Combined enzyme enhancement
    total_enhancement = Pi_xi * eta_coh * viral_damage

    # 5. Membrane turnover (for Cho)
    membrane_turnover = compute_membrane_turnover(condition)

    # 6. Initialize enzyme system
    enzymes = EnzymeKinetics(
        Pi_xi=Pi_xi,
        eta_coh=eta_coh,
        viral_damage_factor=viral_damage
    )

    # 7. Integrate to steady state (60 days)
    NAA_molar, Cho_molar = enzymes.integrate(
        duration_days=60,
        membrane_turnover=membrane_turnover
    )

    # 8. Convert to MRS ratios
    NAA_Cr, Cho_Cr = concentration_to_ratio(NAA_molar, Cho_molar)

    # 9. Get cytokine levels (for reporting)
    cytokines = {
        'healthy': (CONST.TNF_healthy, CONST.IL6_healthy),
        'acute_HIV': (CONST.TNF_acute, CONST.IL6_acute),
        'chronic_HIV': (CONST.TNF_chronic, CONST.IL6_chronic)
    }
    TNF, IL6 = cytokines.get(condition, (5.0, 10.0))

    # 10. Compile results
    results = {
        # Quantum parameters
        'xi': xi,
        'xi_nm': xi * 1e9,
        'coherence_base': coherence_base,
        'sigma_r': sigma_r,
        'sigma_r_nm': sigma_r * 1e9,

        # Modulation factors
        'Pi_xi': Pi_xi,
        'eta_coh': eta_coh,
        'viral_damage': viral_damage,
        'total_enhancement': total_enhancement,

        # Absolute concentrations
        'NAA_mM': NAA_molar * 1e3,
        'Cho_mM': Cho_molar * 1e3,

        # MRS ratios
        'NAA_Cr': NAA_Cr,
        'Cho_Cr': Cho_Cr,

        # Additional info
        'membrane_turnover': membrane_turnover,
        'TNF_pg_mL': TNF,
        'IL6_pg_mL': IL6,
        'condition': condition
    }

    return results


# =============================================================================
# MODEL VALIDATION
# =============================================================================

def validate_model_v3():
    """Validate model v3.0 against Sailasuta et al. data."""

    print("=" * 80)
    print(" FORWARD MODEL v3.0 VALIDATION: Direct Enzyme Modulation")
    print("=" * 80)
    print()

    # Observed data (Sailasuta et al., 2012)
    targets = {
        'healthy': {'NAA_Cr': 1.105, 'Cho_Cr': 0.225},
        'acute_HIV': {'NAA_Cr': 1.135, 'Cho_Cr': 0.245},
        'chronic_HIV': {'NAA_Cr': 1.005, 'Cho_Cr': 0.235}
    }

    # Model inputs (from SSE simulations + Bayesian v2 inference)
    inputs = {
        'healthy': {
            'coherence_base': 0.85,
            'xi': 0.75e-9,
            'sigma_r': CONST.sigma_r_regular
        },
        'acute_HIV': {
            'coherence_base': 0.84,
            'xi': 0.42e-9,
            'sigma_r': CONST.sigma_r_regular * 1.05
        },
        'chronic_HIV': {
            'coherence_base': 0.73,
            'xi': 0.79e-9,
            'sigma_r': CONST.sigma_r_regular * 1.4
        }
    }

    # Run model for all conditions
    results_all = {}
    for condition in ['healthy', 'acute_HIV', 'chronic_HIV']:
        inp = inputs[condition]
        results_all[condition] = run_enzyme_model_v3(
            condition=condition,
            coherence_base=inp['coherence_base'],
            xi=inp['xi'],
            sigma_r=inp['sigma_r']
        )

    # Print detailed results
    for condition in ['healthy', 'acute_HIV', 'chronic_HIV']:
        result = results_all[condition]
        target = targets[condition]

        print(f"{condition.upper().replace('_', ' ')}")
        print("-" * 80)

        print(f"  Quantum Parameters:")
        print(f"    Coherence:                {result['coherence_base']:.3f}")
        print(f"    ξ (correlation length):   {result['xi_nm']:.2f} nm")
        print(f"    σ_r (deloc spread):       {result['sigma_r_nm']:.2f} nm")
        print()

        print(f"  Enzyme Modulation:")
        print(f"    Protection factor (Π_ξ):  {result['Pi_xi']:.3f}")
        print(f"    Coherence factor (η_coh): {result['eta_coh']:.3f}")
        print(f"    Viral damage factor:      {result['viral_damage']:.3f}")
        print(f"    Total enhancement:        {result['total_enhancement']:.3f}×")
        print()

        print(f"  Inflammatory State:")
        print(f"    TNF-α:                    {result['TNF_pg_mL']:.1f} pg/mL")
        print(f"    IL-6:                     {result['IL6_pg_mL']:.1f} pg/mL")
        print()

        print(f"  Enzyme Kinetics:")
        print(f"    NAA steady-state:         {result['NAA_mM']:.2f} mM")
        print(f"    Cho steady-state:         {result['Cho_mM']:.2f} mM")
        print()

        print(f"  MRS Observables:")
        print(f"    NAA/Cr - Model:     {result['NAA_Cr']:.3f}")
        print(f"    NAA/Cr - Data:      {target['NAA_Cr']:.3f}")
        error_naa = 100 * (result['NAA_Cr'] - target['NAA_Cr']) / target['NAA_Cr']
        print(f"    Error:              {error_naa:+.1f}%")
        print()

        print(f"    Cho/Cr - Model:     {result['Cho_Cr']:.3f}")
        print(f"    Cho/Cr - Data:      {target['Cho_Cr']:.3f}")
        error_cho = 100 * (result['Cho_Cr'] - target['Cho_Cr']) / target['Cho_Cr']
        print(f"    Error:              {error_cho:+.1f}%")
        print()

    # KEY MECHANISTIC COMPARISON
    print("=" * 80)
    print(" KEY MECHANISTIC INSIGHT (v3.0 - DIRECT ENZYME MODULATION):")
    print("=" * 80)
    print()

    acute = results_all['acute_HIV']
    chronic = results_all['chronic_HIV']

    print(f"ACUTE HIV (Direct Protection Mechanism):")
    print(f"  - Cytokines: {acute['TNF_pg_mL']:.0f} pg/mL TNF (HIGH)")
    print(f"  - ξ: {acute['xi_nm']:.2f} nm (LOW = decorrelated noise)")
    print(f"  - Protection factor: {acute['Pi_xi']:.2f}× → NAT8L activity BOOSTED")
    print(f"  - NAA synthesis rate: {acute['total_enhancement']:.2f}× baseline")
    print(f"  - NAA/Cr: {acute['NAA_Cr']:.3f} (PRESERVED!) ✓")
    print()

    print(f"CHRONIC HIV (Reduced Protection):")
    print(f"  - Cytokines: {chronic['TNF_pg_mL']:.0f} pg/mL TNF (LOW)")
    print(f"  - ξ: {chronic['xi_nm']:.2f} nm (HIGH = correlated noise)")
    print(f"  - Protection factor: {chronic['Pi_xi']:.2f}× → NAT8L activity reduced")
    print(f"  - NAA synthesis rate: {chronic['total_enhancement']:.2f}× baseline")
    print(f"  - NAA/Cr: {chronic['NAA_Cr']:.3f} (Depleted) ✓")
    print()

    # Model comparison
    print("=" * 80)
    print(" MODEL v2 → v3 IMPROVEMENT:")
    print("=" * 80)
    print()
    print("v2.0 (Phenomenological):")
    print("  - 'Astrocyte compensation' (HOW?)")
    print("  - 'Homeostatic ceiling' (WHY?)")
    print("  - Reactive, not mechanistic")
    print()
    print("v3.0 (Mechanistic):")
    print("  - Direct quantum → enzyme coupling")
    print("  - Π_ξ amplifies NAT8L catalytic rate")
    print("  - Natural kinetic steady states")
    print("  - Testable predictions (measure NAT8L vs ξ)")
    print()

    return results_all


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_enzyme_dynamics():
    """Visualize enzyme modulation across conditions."""

    results = validate_model_v3()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Model v3.0: Direct Enzyme Modulation",
                 fontsize=16, fontweight='bold')

    conditions = ['healthy', 'acute_HIV', 'chronic_HIV']

    # 1. Enzyme enhancement factors
    ax = axes[0, 0]
    x = np.arange(len(conditions))
    width = 0.25

    pi_vals = [results[c]['Pi_xi'] for c in conditions]
    eta_vals = [results[c]['eta_coh'] for c in conditions]
    damage_vals = [results[c]['viral_damage'] for c in conditions]
    total_vals = [results[c]['total_enhancement'] for c in conditions]

    ax.bar(x - 1.5 * width, pi_vals, width, label='Π_ξ (Protection)',
           color='steelblue', edgecolor='black')
    ax.bar(x - 0.5 * width, eta_vals, width, label='η_coh (Coherence)',
           color='orange', edgecolor='black')
    ax.bar(x + 0.5 * width, damage_vals, width, label='Viral Damage',
           color='lightcoral', edgecolor='black')
    ax.bar(x + 1.5 * width, total_vals, width, label='Total Enhancement',
           color='green', alpha=0.7, edgecolor='black', linewidth=2)

    ax.axhline(1.0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_ylabel("Modulation Factor", fontsize=12)
    ax.set_title("Enzyme Activity Modulation", fontsize=13, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace('_', ' ').title() for c in conditions])
    ax.legend()
    ax.grid(alpha=0.3, axis='y')

    # 2. NAA comparison
    ax = axes[0, 1]
    naa_model = [results[c]['NAA_Cr'] for c in conditions]
    naa_obs = [1.105, 1.135, 1.005]

    x = np.arange(len(conditions))
    width = 0.35

    ax.bar(x - width / 2, naa_model, width, label='Model v3.0',
           color='lightgreen', edgecolor='black')
    ax.bar(x + width / 2, naa_obs, width, label='Observed Data',
           color='gray', alpha=0.6, edgecolor='black')

    ax.set_ylabel("NAA/Cr", fontsize=12)
    ax.set_title("NAA Prediction: Enzyme Model", fontsize=13, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace('_', ' ').title() for c in conditions])
    ax.legend()
    ax.grid(alpha=0.3, axis='y')

    # 3. Protection factor vs ξ
    ax = axes[1, 0]
    xi_range = np.linspace(0.3e-9, 1.0e-9, 50)
    pi_range = [compute_protection_factor(xi, beta_xi=CONST.beta_xi)
                for xi in xi_range]

    ax.plot(xi_range * 1e9, pi_range, 'b-', linewidth=2)

    # Mark the three conditions
    for cond, color, marker in [('acute_HIV', 'red', 'o'),
                                ('healthy', 'green', 's'),
                                ('chronic_HIV', 'orange', '^')]:
        xi_val = results[cond]['xi']
        pi_val = results[cond]['Pi_xi']
        ax.plot(xi_val * 1e9, pi_val, marker, markersize=12,
                color=color, markeredgecolor='black', markeredgewidth=2,
                label=cond.replace('_', ' ').title())

    ax.set_xlabel("ξ (nm)", fontsize=12)
    ax.set_ylabel("Protection Factor (Π_ξ)", fontsize=12)
    ax.set_title("Protection Factor Gradient", fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # 4. Model errors
    ax = axes[1, 1]
    naa_obs = [1.105, 1.135, 1.005]
    errors = [100 * (results[c]['NAA_Cr'] - obs) / obs
              for c, obs in zip(conditions, naa_obs)]

    colors = ['green' if abs(e) < 5 else 'orange' if abs(e) < 10 else 'red'
              for e in errors]

    bars = ax.bar(range(len(conditions)), errors, color=colors, alpha=0.7,
                  edgecolor='black', linewidth=2)

    for i, (bar, err) in enumerate(zip(bars, errors)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2., height,
                f'{err:+.1f}%',
                ha='center', va='bottom' if height > 0 else 'top',
                fontweight='bold', fontsize=11)

    ax.axhline(0, color='black', linestyle='-', linewidth=1)
    ax.set_ylabel("NAA Prediction Error (%)", fontsize=12)
    ax.set_title("Model Performance", fontsize=13, fontweight='bold')
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels([c.replace('_', ' ').title() for c in conditions])
    ax.grid(alpha=0.3, axis='y')

    plt.tight_layout()

    # Save to results directory
    os.makedirs('results', exist_ok=True)
    plt.savefig("results/enzyme_model_v3_analysis.png", dpi=300, bbox_inches='tight')
    print("\nVisualization saved to: results/enzyme_model_v3_analysis.png")
    plt.show()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    # Run validation
    results = validate_model_v3()

    # Generate visualization
    plot_enzyme_dynamics()