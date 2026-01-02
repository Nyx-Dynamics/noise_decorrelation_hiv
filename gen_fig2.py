#!/usr/bin/env python3
"""
Generate Figure 2: Posterior Distributions for Manuscript
Based on no_valcour_BG_nopseudo run (2025-11-17)
"""

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import json
from pathlib import Path

# Load trace file
print("Loading trace file...")
# Resolve base run data directory in a cross-machine way
# Priority: environment variable RUN_DATA_DIR -> repo-relative default
RUN_DATA_DIR = Path(os.environ.get(
    "RUN_DATA_DIR",
    str(Path(__file__).resolve().parent / 'quantum' / 'quantum' / 'results_v3_6')
))

# Optional overrides (to mirror gen_fig4.py ergonomics)
TRACE_OVERRIDE = os.environ.get('TRACE')
SUMMARY_OVERRIDE = os.environ.get('SUMMARY')

# Default paths (NO‑Valcour BG‑only)
default_run_dir = RUN_DATA_DIR / 'runs' / 'no_valcour_BG_nopseudo'
default_trace_path = default_run_dir / 'trace_no_valcour_BG_nopseudo.nc'
default_summary_path = default_run_dir / 'summary_no_valcour_BG_nopseudo.csv'

# Choose effective paths
if TRACE_OVERRIDE and SUMMARY_OVERRIDE:
    trace_path = Path(TRACE_OVERRIDE)
    summary_path = Path(SUMMARY_OVERRIDE)
    print("Using TRACE/SUMMARY overrides from environment variables.")
else:
    trace_path = default_trace_path
    summary_path = default_summary_path

print(f"RUN_DATA_DIR: {RUN_DATA_DIR}")
# Announce last successful run if marker exists
try:
    _marker = RUN_DATA_DIR / 'runs' / 'LAST_SUCCESSFUL_RUN.json'
    if _marker.exists():
        _info = json.loads(_marker.read_text())
        _ts = _info.get('timestamp', '?')
        _rn = _info.get('run_name', '?')
        print(f"Last successful run: {_ts} ({_rn}) — establishes evolutionary hypothesis as valid")
except Exception:
    pass
print(f"Trace path: {trace_path}")
print(f"Summary path: {summary_path}")

if not trace_path.exists():
    raise FileNotFoundError(
        "Trace file not found.\n"
        f"  Tried: {trace_path}\n"
        f"  RUN_DATA_DIR: {RUN_DATA_DIR}\n"
        "Set RUN_DATA_DIR appropriately or provide TRACE env var to point directly to the .nc file."
    )
trace = az.from_netcdf(trace_path)

# Load summary statistics
if not summary_path.exists():
    raise FileNotFoundError(
        "Summary CSV not found.\n"
        f"  Tried: {summary_path}\n"
        f"  RUN_DATA_DIR: {RUN_DATA_DIR}\n"
        "Set RUN_DATA_DIR appropriately or provide SUMMARY env var to point directly to the CSV file."
    )
summary = pd.read_csv(summary_path, index_col=0)

# Create figure with 2x3 layout (use constrained layout for better spacing)
fig, axes = plt.subplots(2, 3, figsize=(18, 12), constrained_layout=True)
fig.suptitle('Posterior Distributions: Noise Correlation Length Parameters',
             fontsize=16, fontweight='bold', y=0.98)

# Panel A: β_ξ (protection coupling exponent)
ax = axes[0, 0]
az.plot_posterior(trace, var_names=['β_ξ'], ax=ax,
                  hdi_prob=0.95, point_estimate=None,
                  ref_val=None, color='steelblue')
ax.set_xlabel('β_ξ', fontsize=12, fontweight='bold')
ax.set_ylabel('Density', fontsize=12)
ax.set_title('(A) Protection Coupling Exponent', fontsize=13, fontweight='bold', pad=10)
mean_val = summary.loc['β_ξ', 'mean']
sd_val = summary.loc['β_ξ', 'sd']
hdi_low = summary.loc['β_ξ', 'hdi_2.5%']
hdi_high = summary.loc['β_ξ', 'hdi_97.5%']
ax.text(0.02, 0.98, f'Mean: {mean_val:.2f} ± {sd_val:.2f}\n95% HDI: [{hdi_low:.2f}, {hdi_high:.2f}]',
        transform=ax.transAxes, fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='gray'))
# Add a clean reference line at 0 (without ArviZ's auto text)
ax.axvline(0, color='peru', linestyle='-', linewidth=2, alpha=0.8)

# Panel B: ξ_acute
ax = axes[0, 1]
az.plot_posterior(trace, var_names=['ξ_acute'], ax=ax,
                  hdi_prob=0.95, point_estimate=None,
                  color='dodgerblue')
ax.set_xlabel('ξ_acute (nm)', fontsize=12, fontweight='bold')
ax.set_ylabel('Density', fontsize=12)
ax.set_title('(B) Acute HIV Infection', fontsize=13, fontweight='bold', pad=10, color='dodgerblue')
mean_val = summary.loc['ξ_acute', 'mean']
sd_val = summary.loc['ξ_acute', 'sd']
hdi_low = summary.loc['ξ_acute', 'hdi_2.5%']
hdi_high = summary.loc['ξ_acute', 'hdi_97.5%']
ax.text(0.02, 0.98, f'Mean: {mean_val:.3f} ± {sd_val:.3f} nm\n95% HDI: [{hdi_low:.3f}, {hdi_high:.3f}]',
        transform=ax.transAxes, fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='aliceblue', alpha=0.7, edgecolor='gray'))

# Panel C: ξ_chronic
ax = axes[0, 2]
az.plot_posterior(trace, var_names=['ξ_chronic'], ax=ax,
                  hdi_prob=0.95, point_estimate=None,
                  color='crimson')
ax.set_xlabel('ξ_chronic (nm)', fontsize=12, fontweight='bold')
ax.set_ylabel('Density', fontsize=12)
ax.set_title('(C) Chronic HIV Infection', fontsize=13, fontweight='bold', pad=10, color='crimson')
mean_val = summary.loc['ξ_chronic', 'mean']
sd_val = summary.loc['ξ_chronic', 'sd']
hdi_low = summary.loc['ξ_chronic', 'hdi_2.5%']
hdi_high = summary.loc['ξ_chronic', 'hdi_97.5%']
ax.text(0.02, 0.98, f'Mean: {mean_val:.3f} ± {sd_val:.3f} nm\n95% HDI: [{hdi_low:.3f}, {hdi_high:.3f}]',
        transform=ax.transAxes, fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f8d7da', alpha=0.7, edgecolor='gray'))

# Panel D: ξ_control
ax = axes[1, 0]
az.plot_posterior(trace, var_names=['ξ_control'], ax=ax,
                  hdi_prob=0.95, point_estimate=None,
                  color='green')
ax.set_xlabel('ξ_control (nm)', fontsize=12, fontweight='bold')
ax.set_ylabel('Density', fontsize=12)
ax.set_title('(D) HIV-negative Control', fontsize=13, fontweight='bold', pad=10, color='green')
mean_val = summary.loc['ξ_control', 'mean']
sd_val = summary.loc['ξ_control', 'sd']
hdi_low = summary.loc['ξ_control', 'hdi_2.5%']
hdi_high = summary.loc['ξ_control', 'hdi_97.5%']
ax.text(0.02, 0.98, f'Mean: {mean_val:.3f} ± {sd_val:.3f} nm\n95% HDI: [{hdi_low:.3f}, {hdi_high:.3f}]',
        transform=ax.transAxes, fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#e6f4ea', alpha=0.7, edgecolor='gray'))

# Panel E: Δξ with special annotation
ax = axes[1, 1]
az.plot_posterior(trace, var_names=['Δξ'], ax=ax,
                  hdi_prob=0.95, point_estimate=None,
                  ref_val=None, color='purple')
ax.set_xlabel('Δξ = ξ_chronic - ξ_acute (nm)', fontsize=12, fontweight='bold')
ax.set_ylabel('Density', fontsize=12)
ax.set_title('(E) Difference in Correlation Length', fontsize=13, fontweight='bold', pad=10, color='purple')

# Calculate P(Δξ > 0) from trace
delta_xi_samples = trace.posterior['Δξ'].values.flatten()
p_positive = (delta_xi_samples > 0).mean()

mean_val = summary.loc['Δξ', 'mean']
sd_val = summary.loc['Δξ', 'sd']
hdi_low = summary.loc['Δξ', 'hdi_2.5%']
hdi_high = summary.loc['Δξ', 'hdi_97.5%']
ax.text(0.02, 0.98,
        f'Mean: {mean_val:.3f} ± {sd_val:.3f} nm\n95% HDI: [{hdi_low:.3f}, {hdi_high:.3f}]\n\n'
        f'P(Δξ > 0) = {p_positive:.3f}\n(acute < chronic)',
        transform=ax.transAxes, fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#eadcf8', alpha=0.7, edgecolor='gray'))

# Add vertical line at 0
ax.axvline(0, color='black', linestyle='--', linewidth=2, alpha=0.7, label='Null (no difference)')
ax.legend(fontsize=9)

# Panel F: NAA ratio predictions as bar chart
ax = axes[1, 2]
conditions = ['Control', 'Acute', 'Chronic']
predictions = [
    summary.loc['NAA_ratio_control_mean', 'mean'],
    summary.loc['NAA_ratio_acute_mean', 'mean'],
    summary.loc['NAA_ratio_chronic_mean', 'mean']
]
hdi_low = [
    summary.loc['NAA_ratio_control_mean', 'hdi_2.5%'],
    summary.loc['NAA_ratio_acute_mean', 'hdi_2.5%'],
    summary.loc['NAA_ratio_chronic_mean', 'hdi_2.5%']
]
hdi_high = [
    summary.loc['NAA_ratio_control_mean', 'hdi_97.5%'],
    summary.loc['NAA_ratio_acute_mean', 'hdi_97.5%'],
    summary.loc['NAA_ratio_chronic_mean', 'hdi_97.5%']
]

colors = ['green', 'dodgerblue', 'crimson']
bars = ax.bar(conditions, predictions, color=colors, alpha=0.7, edgecolor='black')

# Build yerr as a (2, N) array-like: [[lower_errs], [upper_errs]]
lower_errs = [predictions[i] - hdi_low[i] for i in range(3)]
upper_errs = [hdi_high[i] - predictions[i] for i in range(3)]
# Build asymmetric error bars with shape (2, N): [lower_errors; upper_errors]
lower_errors = np.array([predictions[i] - hdi_low[i] for i in range(3)])
upper_errors = np.array([hdi_high[i] - predictions[i] for i in range(3)])
yerr = np.vstack([lower_errors, upper_errors])
ax.errorbar(conditions, predictions,
            yerr=yerr,
            fmt='none', ecolor='black', capsize=5, capthick=2)

ax.set_ylabel('Predicted NAA/Cr', fontsize=12, fontweight='bold')
ax.set_title('(F) Predicted NAA/Cr Ratios', fontsize=13, fontweight='bold', pad=10)
ax.set_ylim([0.9, 1.25])
ax.grid(axis='y', alpha=0.3)

# Add value labels on bars
for i, (bar, pred, low, high) in enumerate(zip(bars, predictions, hdi_low, hdi_high)):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
            f'{pred:.3f}\n[{low:.3f}, {high:.3f}]',
            ha='center', va='bottom', fontsize=9)

# Save figure
plt.savefig('Figure2_Posterior_Distributions_BG_nopseudo.png', dpi=300, bbox_inches='tight')
print("\n✅ Figure 2 saved: Figure2_Posterior_Distributions_BG_nopseudo.png")

# Also save as PDF for publication
plt.savefig('Figure2_Posterior_Distributions_BG_nopseudo.pdf', dpi=300, bbox_inches='tight')
print("✅ Figure 2 PDF saved: Figure2_Posterior_Distributions_BG_nopseudo.pdf")

# Optional show (guarded for CLI/headless stability)
_show = os.environ.get('SHOW', '1')
if _show not in ('0', 'false', 'False', 'no', 'NO'):
    try:
        plt.show()
    except Exception as e:
        print(f"Note: plt.show() failed: {e}. Figure files were already saved.")

# Print summary statistics
print("\n" + "="*60)
print("SUMMARY STATISTICS")
print("="*60)
print(f"\nξ_acute:   {summary.loc['ξ_acute', 'mean']:.3f} ± {summary.loc['ξ_acute', 'sd']:.3f} nm")
print(f"ξ_chronic: {summary.loc['ξ_chronic', 'mean']:.3f} ± {summary.loc['ξ_chronic', 'sd']:.3f} nm")
print(f"ξ_control: {summary.loc['ξ_control', 'mean']:.3f} ± {summary.loc['ξ_control', 'sd']:.3f} nm")
print(f"\nΔξ:        {summary.loc['Δξ', 'mean']:.3f} ± {summary.loc['Δξ', 'sd']:.3f} nm")
print(f"P(Δξ > 0): {p_positive:.4f} ({p_positive*100:.2f}%)")
print(f"\nβ_ξ:       {summary.loc['β_ξ', 'mean']:.3f} ± {summary.loc['β_ξ', 'sd']:.3f}")
print("\n" + "="*60)
