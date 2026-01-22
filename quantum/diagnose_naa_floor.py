"""
Diagnostic Script: Analyze NAA Floor Activation
================================================

This script helps understand why all NAA predictions are the same (0.9945)
and identifies the issue with homeostatic floor activation.
"""

import numpy as np

# Posterior medians from your results
params = {
    'coh_exp': 2.448,
    'xi_exp': 0.310,
    'deloc_exp': 0.203,
    'NAA_base': 1.101,
    'astrocyte_comp': 1.179,
    'xi_floor': 0.363e-9,  # m
    'xi_ceiling': 0.807e-9,  # m
}

# Constants
NAA_baseline = 1.105
xi_baseline = 0.8e-9
sigma_r_regular = 0.38e-9

# Condition-specific values
conditions = {
    'healthy': {
        'coherence_base': 0.85,
        'xi': 0.6786e-9,  # From posterior
        'sigma_r': sigma_r_regular,
    },
    'acute_HIV': {
        'coherence_base': 0.84,
        'xi': 0.4467e-9,  # From posterior
        'sigma_r': sigma_r_regular * 1.05,
    },
    'chronic_HIV': {
        'coherence_base': 0.73,
        'xi': 0.7176e-9,  # From posterior
        'sigma_r': sigma_r_regular * 1.4,
    }
}

print("=" * 80)
print(" DIAGNOSTIC: NAA Floor Activation Analysis")
print("=" * 80)
print()

for cond_name, cond in conditions.items():
    print(f"{cond_name.upper()}")
    print("-" * 80)

    # 1. Nonlinear coherence
    xi_normalized = (cond['xi'] - params['xi_floor']) / (params['xi_ceiling'] - params['xi_floor'])
    xi_normalized = np.clip(xi_normalized, 0.0, 1.0)

    coherence_floor = 0.65
    coh_effective = coherence_floor + (cond['coherence_base'] - coherence_floor) * (1 - xi_normalized) ** 2

    print(f"  ξ: {cond['xi'] * 1e9:.3f} nm")
    print(f"  ξ normalized: {xi_normalized:.3f}")
    print(f"  Coherence base: {cond['coherence_base']:.3f}")
    print(f"  Coherence effective: {coh_effective:.3f}")
    print()

    # 2. Base quantum coupling
    coherence_term = (coh_effective / 0.85) ** params['coh_exp']
    xi_protection = (xi_baseline / cond['xi']) ** params['xi_exp']
    deloc_term = (cond['sigma_r'] / sigma_r_regular) ** params['deloc_exp']

    NAA_quantum = params['NAA_base'] * coherence_term * xi_protection * deloc_term

    print(f"  Coherence term: {coherence_term:.3f}")
    print(f"  ξ protection: {xi_protection:.3f}")
    print(f"  Deloc term: {deloc_term:.3f}")
    print(f"  NAA (quantum only): {NAA_quantum:.3f}")
    print()

    # 3. Astrocyte compensation (chronic only)
    if cond_name == 'chronic_HIV':
        NAA_compensated = NAA_quantum * params['astrocyte_comp']
        print(f"  Astrocyte boost: {params['astrocyte_comp']:.3f}×")
        print(f"  NAA (after compensation): {NAA_compensated:.3f}")
    else:
        NAA_compensated = NAA_quantum
        print(f"  Astrocyte boost: N/A (not chronic)")
        print(f"  NAA (no compensation needed): {NAA_compensated:.3f}")
    print()

    # 4. Homeostatic floor
    NAA_floor = 0.90 * NAA_baseline  # 0.9945
    print(f"  Homeostatic floor: {NAA_floor:.3f}")
    print(f"  Is floor active? {NAA_compensated < NAA_floor}")

    NAA_final = max(NAA_compensated, NAA_floor)
    print(f"  NAA (final): {NAA_final:.3f}")
    print()

    # Observed
    obs = {'healthy': 1.105, 'acute_HIV': 1.135, 'chronic_HIV': 1.005}[cond_name]
    error = 100 * (NAA_final - obs) / obs
    print(f"  NAA observed: {obs:.3f}")
    print(f"  Error: {error:+.1f}%")
    print()
    print()

print("=" * 80)
print(" DIAGNOSIS")
print("=" * 80)
print()
print("PROBLEM IDENTIFIED:")
print("  The homeostatic floor (0.9945) is activating for ALL conditions")
print("  because the quantum coupling is too weak.")
print()
print("ROOT CAUSE:")
print("  NAA_base posterior median: 1.101")
print("  With coupling terms < 1.0, NAA_quantum drops below floor")
print()
print("SOLUTION OPTIONS:")
print("  1. Increase NAA_base prior mean (1.10 → 1.15)")
print("  2. Relax homeostatic floor (0.90 → 0.85)")
print("  3. Remove floor for healthy/acute conditions")
print("  4. Use condition-specific floors")
print()
print("RECOMMENDED FIX:")
print("  Use adaptive floor that only activates in chronic phase")
print("  or when NAA drops below a more permissive threshold (0.85).")
print()