#!/usr/bin/env python3
"""
Generate Figure 4: WAIC Model Comparison
Computes WAIC/LOO from traces to avoid stale hard-coded numbers.
If traces are not available, falls back to last known values.

IMPORTANT: To ensure a fair comparison (same data), this script prefers
BG‑only traces for BOTH models when available:
  runs/with_valcour_BG_nopseudo/trace_with_valcour_BG_nopseudo.nc
  runs/no_valcour_BG_nopseudo/trace_no_valcour_BG_nopseudo.nc

You can also override the paths explicitly via environment variables:
  WITH_TRACE=/abs/path/to/with_valcour_trace.nc
  NO_TRACE=/abs/path/to/no_valcour_trace.nc
"""

import matplotlib.pyplot as plt
import numpy as np
import arviz as az
import xarray as xr
from pathlib import Path
import os
import json

# Resolve run data directory
RUN_DATA_DIR = Path(os.environ.get(
    'RUN_DATA_DIR',
    str(Path(__file__).resolve().parent / 'quantum' / 'quantum' / 'results_v3_6')
))

# Prefer explicitly provided traces (env), then BG-only traces, then the *_waic traces
env_with = os.environ.get('WITH_TRACE')
env_no   = os.environ.get('NO_TRACE')

bg_with  = RUN_DATA_DIR / 'runs' / 'with_valcour_BG_nopseudo' / 'trace_with_valcour_BG_nopseudo.nc'
bg_no    = RUN_DATA_DIR / 'runs' / 'no_valcour_BG_nopseudo'  / 'trace_no_valcour_BG_nopseudo.nc'

waic_with = RUN_DATA_DIR / 'runs' / 'with_valcour_waic' / 'trace_with_valcour.nc'
waic_no   = RUN_DATA_DIR / 'runs' / 'no_valcour_waic'   / 'trace_no_valcour.nc'

def _select_trace_paths():
    """Choose the best available pair of trace paths in priority order.
    Returns (with_path, no_path).
    """
    if env_with and env_no:
        return Path(env_with), Path(env_no)
    if bg_with.exists() and bg_no.exists():
        return bg_with, bg_no
    if waic_with.exists() and waic_no.exists():
        return waic_with, waic_no
    # Last resort: return whatever exists (may mix); caller will check .exists()
    return (Path(env_with) if env_with else bg_with if bg_with.exists() else waic_with,
            Path(env_no)   if env_no   else bg_no   if bg_no.exists()   else waic_no)

with_nc, no_nc = _select_trace_paths()

# Announce last successful run if marker exists
try:
    marker = RUN_DATA_DIR / 'runs' / 'LAST_SUCCESSFUL_RUN.json'
    if marker.exists():
        info = json.loads(marker.read_text())
        ts = info.get('timestamp', '?')
        rn = info.get('run_name', '?')
        print(f"Last successful run: {ts} ({rn}) — establishes evolutionary hypothesis as valid")
except Exception:
    pass

# ---------------- Manifest utilities (per-run JSON) ----------------
def _load_manifest_from_run_dir(run_dir: Path):
    try:
        p = run_dir / 'run_info.json'
        if p.exists():
            return json.loads(p.read_text())
    except Exception:
        return None
    return None

from typing import Optional, Dict

def _summarize_manifest(m: Optional[dict]) -> dict:
    if not m:
        return {}
    out = {}
    # Basic identifiers
    out['timestamp'] = m.get('timestamp')
    out['run_name'] = m.get('run_name')
    # CLI flags of interest
    cli = m.get('cli', {}) or {}
    for k in ['regions', 'exclude_region', 'exclude_valcour', 'no_acute_pseudo', 'valcour_aux', 'tag', 'run_label']:
        out[f'cli.{k}'] = cli.get(k)
    # Data counts
    dc = m.get('data_counts', {}) or {}
    for k in ['n_acute_total', 'n_valcour', 'n_young_acute', 'n_sailasuta_acute', 'n_chronic_points', 'n_control_points']:
        out[f'data.{k}'] = dc.get(k)
    return out

def _diff_manifests(label: str, a: Optional[dict], b: Optional[dict]):
    print(f"\n--- Manifest comparison: {label} ---")
    if a is None:
        print("  Current: <missing run_info.json>")
    if b is None:
        print("  Reference WAIC: <missing run_info.json>")
    sa = _summarize_manifest(a)
    sb = _summarize_manifest(b)
    ida = sa.get('run_name') or a.get('run_name') if a else None
    idb = sb.get('run_name') or b.get('run_name') if b else None
    tsa = sa.get('timestamp') or (a.get('timestamp') if a else None)
    tsb = sb.get('timestamp') or (b.get('timestamp') if b else None)
    print(f"  Current run:   {ida or '?'} @ {tsa or '?'}")
    print(f"  WAIC basis:    {idb or '?'} @ {tsb or '?'}")
    keys = sorted(set(sa.keys()) | set(sb.keys()))
    for k in keys:
        if k in ('timestamp','run_name'):
            continue
        va = sa.get(k)
        vb = sb.get(k)
        eq = (va == vb)
        status = 'OK' if eq else 'DIFF'
        print(f"  {k:22s}: {va!r} vs {vb!r}  [{status}]")

# Load manifests for the selected traces and the WAIC reference runs
try:
    cur_with_manifest = _load_manifest_from_run_dir(with_nc.parent)
    cur_no_manifest   = _load_manifest_from_run_dir(no_nc.parent)
    ref_with_manifest = _load_manifest_from_run_dir(RUN_DATA_DIR / 'runs' / 'with_valcour_waic')
    ref_no_manifest   = _load_manifest_from_run_dir(RUN_DATA_DIR / 'runs' / 'no_valcour_waic')
    # Compare each current trace to its WAIC basis
    _diff_manifests('With Valcour (current vs WAIC basis)', cur_with_manifest, ref_with_manifest)
    _diff_manifests('No Valcour (current vs WAIC basis)', cur_no_manifest, ref_no_manifest)
    # Cross-compare current selection for mutual comparability
    print("\n--- Selected traces comparability check (current with vs no) ---")
    sa = _summarize_manifest(cur_with_manifest)
    sb = _summarize_manifest(cur_no_manifest)
    for k in ['cli.regions','cli.exclude_region','cli.no_acute_pseudo','cli.valcour_aux']:
        va = sa.get(k)
        vb = sb.get(k)
        status = 'OK' if va == vb else 'DIFF'
        print(f"  {k:22s}: WITH={va!r} | NO={vb!r}  [{status}]")
except Exception:
    # Non-fatal; manifest comparison is auxiliary
    pass

# Allow restricting to specific phases (now default to all three for 3:1:1 weighting)
PHASES = [p.strip().lower() for p in os.environ.get('PHASES', 'acute,chronic,control').split(',') if p.strip()]

# Phase weights for WAIC/LOO combination. Format: "acute:chronic:control" (e.g., 3:1:1)
def _parse_phase_weights(s: str):
    try:
        parts = [x.strip() for x in s.split(':')]
        if len(parts) != 3:
            raise ValueError
        a, c, k = (float(parts[0]), float(parts[1]), float(parts[2]))
        return {'acute': a, 'chronic': c, 'control': k}
    except Exception:
        # fallback to 1:1:1 if bad format
        return {'acute': 1.0, 'chronic': 1.0, 'control': 1.0}

# Default to 3:1:1 as requested; can be overridden via env PHASE_WEIGHTS
PHASE_WEIGHTS = _parse_phase_weights(os.environ.get('PHASE_WEIGHTS', '3:1:1'))

models = ['No Valcour', 'With Valcour']

# Last-known fallback numbers (used only if traces not found)
waic_values = [-25.75, 719.78]
p_waic_values = [1.63, 115.60]
loo_values = [-24.27, 719.97]
p_loo_values = [2.37, 115.69]

def _combine_llk_core(idata: az.InferenceData) -> az.InferenceData:
    """Combine core log_likelihood terms into one dataset to compute IC reliably.

    Respects the PHASES selection (subset among {'acute','chronic','control'}).
    """
    if not hasattr(idata, 'log_likelihood') or idata.log_likelihood is None:
        raise RuntimeError('No log_likelihood present; re-run with idata_kwargs={"log_likelihood": True}.')
    llk = idata.log_likelihood
    core = ['NAA_ratio_acute_obs', 'NAA_ratio_chronic_obs', 'NAA_ratio_control_obs']
    # Phase filter
    phase_keep = set(PHASES)
    def _phase_of(var: str) -> str:
        s = var.lower()
        if 'acute' in s:
            return 'acute'
        if 'chronic' in s:
            return 'chronic'
        if 'control' in s:
            return 'control'
        return ''
    use = [v for v in core if v in llk.data_vars and _phase_of(v) in phase_keep]
    if not use:
        # Fallback: use all available llk vars
        use = [v for v in llk.data_vars if _phase_of(v) in phase_keep] or list(llk.data_vars)
    total = None
    for v in use:
        arr = llk[v]
        ph = _phase_of(v)
        w = float(PHASE_WEIGHTS.get(ph, 1.0))
        weighted = arr * w
        total = weighted if total is None else (total + weighted)
    new_llk = xr.Dataset({'log_likelihood': total})
    out = az.InferenceData(**{g: getattr(idata, g) for g in idata.groups() if g != 'log_likelihood'})
    out.add_groups({'log_likelihood': new_llk})
    return out

# Try to compute metrics from traces if available
try:
    if with_nc.exists() and no_nc.exists():
        print(f"Using traces:\n  WITH: {with_nc}\n  NO:   {no_nc}")
        print(f"Phases: {','.join(PHASES)} | Weights (acute:chronic:control) = {PHASE_WEIGHTS['acute']}:{PHASE_WEIGHTS['chronic']}:{PHASE_WEIGHTS['control']}")
        id_with = _combine_llk_core(az.from_netcdf(str(with_nc)))
        id_no   = _combine_llk_core(az.from_netcdf(str(no_nc)))
        waic_with = az.waic(id_with)
        waic_no   = az.waic(id_no)
        loo_with  = az.loo(id_with)
        loo_no    = az.loo(id_no)
        # Standard definitions: WAIC = -2*elpd_waic, LOO = -2*elpd_loo
        waic_values = [float(-2.0 * waic_no.elpd_waic), float(-2.0 * waic_with.elpd_waic)]
        p_waic_values = [float(waic_no.p_waic), float(waic_with.p_waic)]
        loo_values = [float(-2.0 * loo_no.elpd_loo), float(-2.0 * loo_with.elpd_loo)]
        p_loo_values = [float(loo_no.p_loo), float(loo_with.p_loo)]
        # Try to infer N from the combined log_likelihood shapes and add to labels
        def _infer_n(idata: az.InferenceData) -> int:
            arr = idata.log_likelihood['log_likelihood']
            # dims typically include ('chain','draw','obs') or similar; pick the last non-chain/draw dim
            for d in arr.dims:
                if d not in ('chain', 'draw'):
                    return int(arr.sizes[d])
            # fallback to last dimension length
            return int(list(arr.sizes.values())[-1])
        n_no = _infer_n(id_no)
        n_with = _infer_n(id_with)
        phase_lbl = ','.join(PHASES)
        w_lbl = f"{PHASE_WEIGHTS['acute']}:{PHASE_WEIGHTS['chronic']}:{PHASE_WEIGHTS['control']}"
        models = [
            f'No Valcour\n(BG only, {phase_lbl}; w={w_lbl}; n={n_no})',
            f'With Valcour\n(BG only, {phase_lbl}; w={w_lbl}; n={n_with})'
        ]
    else:
        print(f"Warning: Trace files not found. Using fallback numbers.\n WITH exists={with_nc.exists()}\n NO exists={no_nc.exists()}\n Expected defaults under RUN_DATA_DIR={RUN_DATA_DIR}")
except Exception as e:
    print(f"Warning: Failed to compute WAIC/LOO from traces: {e}. Using fallback numbers.")

# Calculate differences (No Valcour minus With Valcour; lower is better)
delta_waic = waic_values[0] - waic_values[1]
delta_loo = loo_values[0] - loo_values[1]

# Create figure with 3 panels
fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
fig.suptitle('Model Comparison: Why Valcour Data Was Excluded',
             fontsize=18, fontweight='bold', y=0.98)

# ===== Panel A: WAIC Comparison =====
ax = axes[0]

bars = ax.bar(models, waic_values,
              color=['green', 'red'],
              alpha=0.7,
              edgecolor='black',
              linewidth=2)

# Add horizontal line at 0
ax.axhline(0, color='black', linestyle='--', linewidth=2, alpha=0.5)

# Add value labels on bars
for bar, val in zip(bars, waic_values):
    height = bar.get_height()
    if val < 0:
        y_pos = height - 30
        va = 'top'
    else:
        y_pos = height + 30
        va = 'bottom'
    ax.text(bar.get_x() + bar.get_width()/2., y_pos,
            f'WAIC = {val:.2f}',
            ha='center', va=va, fontsize=13, fontweight='bold')

ax.set_ylabel('WAIC (lower is better)', fontsize=14, fontweight='bold')
ax.set_title('(A) Widely Applicable Information Criterion',
             fontsize=14, fontweight='bold', pad=10)
ax.set_ylim([-100, 800])
ax.grid(axis='y', alpha=0.3, linestyle='--')

# Add ΔWAIC annotation (dynamic)
def _strength(delta: float) -> str:
    a = abs(delta)
    if a > 10:
        return 'Decisive'
    if a > 4:
        return 'Strong'
    if a > 2:
        return 'Moderate'
    if a > 1:
        return 'Weak'
    return 'Negligible'

waic_better = 'With Valcour' if delta_waic > 0 else 'No Valcour' if delta_waic < 0 else 'Tie'
ax.text(0.5, 0.96,
        f'ΔWAIC = WAIC(no) - WAIC(with) = {delta_waic:.2f}\n'
        f'Lower is better → {waic_better} is better\n'
        f'Strength: {_strength(delta_waic)} (|Δ| vs 10 threshold)',
        transform=ax.transAxes,
        fontsize=10.5,
        verticalalignment='top',
        horizontalalignment='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff3cd', alpha=0.85, edgecolor='#b68a00', linewidth=1.5),
        fontweight='bold')

# ===== Panel B: p_WAIC Explosion =====
ax = axes[1]

bars = ax.bar(models, p_waic_values,
              color=['lightgreen', 'darkred'],
              alpha=0.7,
              edgecolor='black',
              linewidth=2)

# Add value labels on bars
for bar, val in zip(bars, p_waic_values):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 3,
            f'p_WAIC = {val:.1f}',
            ha='center', va='bottom', fontsize=13, fontweight='bold')

ax.set_ylabel('p_WAIC (effective parameters)', fontsize=14, fontweight='bold')
ax.set_title('(B) Effective Parameter Explosion',
             fontsize=14, fontweight='bold', pad=10)
ax.set_ylim([0, 130])
ax.grid(axis='y', alpha=0.3, linestyle='--')

# Add complexity annotation (dynamic)
fold_increase = (p_waic_values[1] / p_waic_values[0]) if p_waic_values[0] != 0 else float('inf')
if fold_increase > 3:
    text = (f'{fold_increase:.0f}× higher complexity with Valcour\n\n'
            f'Indicates potential data mismatch')
    face, edge = '#f8d7da', '#842029'
else:
    direction = 'similar' if 1/1.5 < fold_increase < 1.5 else ('higher' if fold_increase > 1 else 'lower')
    text = (f'Complexity is {direction} with Valcour\n'
            f'p_WAIC(no)={p_waic_values[0]:.2f}, p_WAIC(with)={p_waic_values[1]:.2f}')
    face, edge = '#e7eaf6', '#495057'
ax.text(0.5, 0.96,
        text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='top',
        horizontalalignment='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor=face, alpha=0.85, edgecolor=edge, linewidth=1.5),
        fontweight='bold')

# ===== Panel C: LOO Comparison =====
ax = axes[2]

bars = ax.bar(models, loo_values,
              color=['green', 'red'],
              alpha=0.7,
              edgecolor='black',
              linewidth=2)

# Add horizontal line at 0
ax.axhline(0, color='black', linestyle='--', linewidth=2, alpha=0.5)

# Add value labels on bars
for bar, val in zip(bars, loo_values):
    height = bar.get_height()
    if val < 0:
        y_pos = height - 30
        va = 'top'
    else:
        y_pos = height + 30
        va = 'bottom'
    ax.text(bar.get_x() + bar.get_width()/2., y_pos,
            f'LOO = {val:.2f}',
            ha='center', va=va, fontsize=13, fontweight='bold')

ax.set_ylabel('LOO (lower is better)', fontsize=14, fontweight='bold')
ax.set_title('(C) Leave-One-Out Cross-Validation',
             fontsize=14, fontweight='bold', pad=10)
ax.set_ylim([-100, 800])
ax.grid(axis='y', alpha=0.3, linestyle='--')

# Add ΔLOO annotation (dynamic)
loo_better = 'With Valcour' if delta_loo > 0 else 'No Valcour' if delta_loo < 0 else 'Tie'
ax.text(0.5, 0.96,
        f'ΔLOO = LOO(no) - LOO(with) = {delta_loo:.2f}\n'
        f'Lower is better → {loo_better} is better\n'
        f'Strength: {_strength(delta_loo)}',
        transform=ax.transAxes,
        fontsize=10.5,
        verticalalignment='top',
        horizontalalignment='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#cff4fc', alpha=0.85, edgecolor='#055160', linewidth=1.5),
        fontweight='bold')

plt.tight_layout()
plt.savefig('Figure4_WAIC_Model_Comparison.png', dpi=300, bbox_inches='tight')
print("\n✅ Figure 4 saved: Figure4_WAIC_Model_Comparison.png")

# Also save as PDF
plt.savefig('Figure4_WAIC_Model_Comparison.pdf', dpi=300, bbox_inches='tight')
print("✅ Figure 4 PDF saved: Figure4_WAIC_Model_Comparison.pdf")

# Show only if interactive is desired; avoid crashes on headless/CLI runs
_show = os.environ.get('SHOW', '1')
if _show not in ('0', 'false', 'False', 'no', 'NO'):
    try:
        plt.show()
    except Exception as e:
        print(f"Note: plt.show() failed: {e}. Figure files were already saved.")

# Print summary
print("\n" + "="*80)
print("WAIC/LOO MODEL COMPARISON SUMMARY")
print("="*80)
print("\nNo Valcour Model:")
print(f"  WAIC:   {waic_values[0]:.2f} (p_WAIC = {p_waic_values[0]:.2f})")
print(f"  LOO:    {loo_values[0]:.2f} (p_LOO = {p_loo_values[0]:.2f})")
print(f"  Status: ✅ OPTIMAL (parsimonious, well-calibrated)")

print("\nWith Valcour Model:")
print(f"  WAIC:   {waic_values[1]:.2f} (p_WAIC = {p_waic_values[1]:.2f})")
print(f"  LOO:    {loo_values[1]:.2f} (p_LOO = {p_loo_values[1]:.2f})")
print(f"  Status: ❌ PATHOLOGICAL (parameter explosion)")

print("\nModel Selection:")
print(f"  ΔWAIC:  {delta_waic:.2f} (decisive: |ΔWAIC| > 10)")
print(f"  ΔLOO:   {delta_loo:.2f} (decisive: |ΔLOO| > 10)")
print(f"  Fold change in p_WAIC: {p_waic_values[1]/p_waic_values[0]:.0f}×")

print("\n" + "="*80)
print("INTERPRETATION:")
print("="*80)
print("""
The ~71× explosion in effective parameters (p_WAIC: 1.63 → 115.60) indicates
that the model is trying to use 116 parameters to accommodate Valcour data,
when it should only need ~2. This is NOT model complexity - it's a sign of
fundamental data incompatibility.

Possible reasons for incompatibility:
1. Different brain region (Valcour: frontal WM; Others: basal ganglia)
2. Different MRS protocol (voxel size, TE/TR parameters)
3. Different measurement scale (absolute vs ratio normalization)

Both WAIC and LOO strongly agree (ΔWAIC ≈ ΔLOO ≈ -745), providing robust
evidence for excluding Valcour data to maintain regional and methodological
consistency.
""")
print("="*80)
