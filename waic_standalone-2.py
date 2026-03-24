#!/usr/bin/env python3
"""
waic_standalone.py
------------------
Standalone WAIC / LOO model comparison for noise_decorrelation_HIV.

Loads existing ArviZ .nc trace files — NO new MCMC sampling.
Writes results to outputs/waic_comparison/waic_results.json
and outputs/waic_comparison/waic_summary.txt

Usage (from repo root):
    # Use .waicenv (Python 3.12, likely has newer ArviZ):
    .waicenv/bin/python waic_standalone.py

    # Or use .venv (Python 3.9):
    .venv/bin/python waic_standalone.py

    # Or system python if arviz is installed:
    python waic_standalone.py
"""

import os
import sys
import json
import warnings
from datetime import datetime
from pathlib import Path

warnings.filterwarnings("ignore")

# ── Dependency check ─────────────────────────────────────────────────────────

def check_deps():
    missing = []
    for pkg in ["arviz", "numpy", "xarray"]:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    if missing:
        print(f"[ERROR] Missing packages: {', '.join(missing)}")
        print("Install with:")
        print(f"  {sys.executable} -m pip install arviz numpy xarray")
        sys.exit(1)

check_deps()

import arviz as az
import numpy as np

print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] ArviZ version: {az.__version__}")
print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] Python: {sys.executable}")

# ── Trace file paths (relative to repo root) ─────────────────────────────────

REPO_ROOT = Path(__file__).parent

TRACES = {
    "v3.6_primary": REPO_ROOT / "results/bayesian_v3_6/trace.nc",
    "v3.6_with_valcour": REPO_ROOT / "evolutionary_framework_repo/results_v3_6/runs/run_20260221_145358/trace_with_valcour.nc",
    "v3.6_no_valcour": REPO_ROOT / "evolutionary_framework_repo/results_v3_6/runs/run_20260221_145603/trace_no_valcour.nc",
    "v4_enzyme": REPO_ROOT / "results/enzyme_v4/trace_v4.nc",
    "v2_individual": REPO_ROOT / "evolutionary_framework_repo/results/bayesian_v2/trace_v2.nc",
}

# ── Output directory ──────────────────────────────────────────────────────────

OUTPUT_DIR = REPO_ROOT / "outputs" / "waic_comparison"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load traces ───────────────────────────────────────────────────────────────

def load_trace(name, path):
    path = Path(path)
    if not path.exists():
        print(f"  [SKIP] {name}: file not found at {path}")
        return None
    size_mb = path.stat().st_size / 1e6
    print(f"  [LOAD] {name}: {path.name} ({size_mb:.1f} MB)")
    try:
        trace = az.from_netcdf(str(path))
        chains = trace.posterior.dims.get("chain", "?")
        draws = trace.posterior.dims.get("draw", "?")
        has_ll = hasattr(trace, "log_likelihood")
        ll_vars = list(trace.log_likelihood.data_vars) if has_ll else []
        ll_status = f"{len(ll_vars)} LL arrays" if ll_vars else "NO log_likelihood"
        print(f"         chains={chains}  draws={draws}  {ll_status}")
        if ll_vars and len(ll_vars) <= 5:
            print(f"         vars: {ll_vars}")
        elif ll_vars:
            print(f"         vars: {ll_vars[:3]} ... (+{len(ll_vars)-3} more)")
        return trace
    except Exception as e:
        print(f"  [FAIL] {name}: {e}")
        return None

print(f"\n[{datetime.now():%Y-%m-%d %H:%M:%S}] Loading traces...")
loaded = {}
for name, path in TRACES.items():
    t = load_trace(name, path)
    if t is not None:
        loaded[name] = t

if len(loaded) < 2:
    print(f"\n[ERROR] Need at least 2 traces for comparison. Found {len(loaded)}.")
    sys.exit(1)

print(f"\n[{datetime.now():%Y-%m-%d %H:%M:%S}] Loaded {len(loaded)} traces: {list(loaded.keys())}")

# ── WAIC computation ──────────────────────────────────────────────────────────

def get_combined_log_likelihood(trace):
    """
    Extract and combine log-likelihood from a trace.

    Handles three cases:
    1. Single log_likelihood array       → use directly
    2. Multiple named arrays             → sum across all observation groups
    3. No log_likelihood group           → return None (cannot compute WAIC)

    Returns an xarray DataArray with dims (chain, draw, obs) or None.
    """
    if not hasattr(trace, "log_likelihood"):
        return None

    ll_ds = trace.log_likelihood
    ll_vars = list(ll_ds.data_vars)

    if len(ll_vars) == 0:
        return None

    if len(ll_vars) == 1:
        # Single array — use directly
        ll = ll_ds[ll_vars[0]]
        # Flatten all obs dims into one
        chain_draw_dims = ["chain", "draw"]
        obs_dims = [d for d in ll.dims if d not in chain_draw_dims]
        if obs_dims:
            ll = ll.stack(obs_point=obs_dims)
        return ll

    # Multiple arrays — stack each to (chain, draw, N_i) then concat along obs
    import xarray as xr
    parts = []
    for var in ll_vars:
        arr = ll_ds[var]
        chain_draw_dims = ["chain", "draw"]
        obs_dims = [d for d in arr.dims if d not in chain_draw_dims]
        if obs_dims:
            arr = arr.stack(obs_point=obs_dims)
        else:
            # Scalar observation — add obs dim
            arr = arr.expand_dims("obs_point")
        parts.append(arr)

    # Rename obs_point coord to avoid conflicts, then concat
    # Use simple integer indexing for concat to avoid dimension mismatch
    parts_clean = []
    for p in parts:
        p = p.drop_vars("obs_point", errors="ignore")
        parts_clean.append(p)
    
    combined = xr.concat(parts_clean, dim="obs_point")
    return combined


def make_idata_with_combined_ll(trace, combined_ll):
    """
    Create a new InferenceData object with a single combined
    log_likelihood variable for WAIC/LOO computation.
    """
    import xarray as xr

    ll_ds = xr.Dataset({"combined": combined_ll})
    new_idata = az.InferenceData(
        posterior=trace.posterior,
        log_likelihood=ll_ds,
    )
    # Carry over observed_data if present
    if hasattr(trace, "observed_data"):
        new_idata.add_groups({"observed_data": trace.observed_data})
    return new_idata


def compute_waic(name, trace):
    """Compute WAIC for a single trace. Returns dict or None."""
    combined_ll = get_combined_log_likelihood(trace)

    if combined_ll is None:
        print(f"    [SKIP] No log_likelihood in {name} — cannot compute WAIC")
        return None

    # Build a clean idata with single combined LL
    idata = make_idata_with_combined_ll(trace, combined_ll)

    try:
        # ArviZ stats (waic, loo) return an ELPDData object.
        # Computing WAIC for a combined LL can be slow for large traces.
        print(f"      [CALC] az.waic(idata, var_name='combined', pointwise=False)...")
        waic_result = az.waic(idata, var_name="combined", pointwise=False)
        # Access via .get() or [key] to avoid potential attribute access issues in some environments
        # ArviZ version 0.17.1 uses 'waic' and 'loo' keys instead of 'elpd_waic' and 'elpd_loo'
        elpd_waic = 0.0
        if "waic" in waic_result:
            elpd_waic = float(waic_result["waic"])
        elif "elpd_waic" in waic_result:
            elpd_waic = float(waic_result["elpd_waic"])
        elif "elpd" in waic_result:
            elpd_waic = float(waic_result["elpd"])
        
        return {
            "model": name,
            "elpd_waic": elpd_waic,
            "se": float(waic_result.get("se", 0)),
            "p_waic": float(waic_result.get("p_waic", 0)),
            "warning": bool(getattr(waic_result, "warning", False)),
            "n_obs": int(combined_ll.sizes.get("obs_point", "?")),
        }
    except Exception as e:
        print(f"    [WARN] WAIC failed for {name}: {e}")
        return None


def compute_loo(name, trace):
    """Compute LOO-CV for a single trace. Returns dict or None."""
    combined_ll = get_combined_log_likelihood(trace)

    if combined_ll is None:
        print(f"    [SKIP] No log_likelihood in {name} — cannot compute LOO")
        return None

    idata = make_idata_with_combined_ll(trace, combined_ll)
    try:
        # Computing LOO can be computationally intensive and may appear to "hang".
        print(f"      [CALC] az.loo(idata, var_name='combined', pointwise=False)...")
        loo_result = az.loo(idata, var_name="combined", pointwise=False)
        
        elpd_loo = 0.0
        if "loo" in loo_result:
            elpd_loo = float(loo_result["loo"])
        elif "elpd_loo" in loo_result:
            elpd_loo = float(loo_result["elpd_loo"])
        elif "elpd" in loo_result:
            elpd_loo = float(loo_result["elpd"])

        return {
            "model": name,
            "elpd_loo": elpd_loo,
            "se": float(loo_result.get("se", 0)),
            "p_loo": float(loo_result.get("p_loo", 0)),
            "warning": bool(getattr(loo_result, "warning", False)),
            "n_obs": int(combined_ll.sizes.get("obs_point", "?")),
        }
    except Exception as e:
        print(f"    [WARN] LOO failed for {name}: {e}")
        return None

print(f"\n[{datetime.now():%Y-%m-%d %H:%M:%S}] Computing WAIC...")
waic_results = {}
for name, trace in loaded.items():
    print(f"  Computing WAIC for {name}...")
    result = compute_waic(name, trace)
    if result:
        waic_results[name] = result
        print(f"    elpd_waic={result['elpd_waic']:.3f}  se={result['se']:.3f}  "
              f"p_waic={result['p_waic']:.3f}")

print(f"\n[{datetime.now():%Y-%m-%d %H:%M:%S}] Computing LOO...")
loo_results = {}
for name, trace in loaded.items():
    print(f"  Computing LOO for {name}...")
    result = compute_loo(name, trace)
    if result:
        loo_results[name] = result
        print(f"    elpd_loo={result['elpd_loo']:.3f}  se={result['se']:.3f}")

# ── Model comparison table ────────────────────────────────────────────────────

def compute_comparison(results_dict, metric_key):
    """Rank models and compute weights from ELPD values."""
    if not results_dict:
        return []

    items = list(results_dict.values())
    # Sort by ELPD (higher = better)
    items.sort(key=lambda x: x[metric_key], reverse=True)

    # Best model ELPD
    best_elpd = items[0][metric_key]

    # Stacking weights via softmax on ELPD differences
    elpds = np.array([x[metric_key] for x in items])
    # Simple Bayesian model weights (proportional to exp(ELPD - max))
    weights_raw = np.exp(elpds - elpds.max())
    weights = weights_raw / weights_raw.sum()

    comparison = []
    for i, (item, w) in enumerate(zip(items, weights)):
        comparison.append({
            "rank": i + 1,
            "model": item["model"],
            metric_key: item[metric_key],
            "se": item["se"],
            "delta_elpd": item[metric_key] - best_elpd,
            "weight": float(w),
            "warning": item.get("warning", False),
        })
    return comparison

waic_comparison = compute_comparison(waic_results, "elpd_waic")
loo_comparison = compute_comparison(loo_results, "elpd_loo")

# ── Save results ──────────────────────────────────────────────────────────────

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

results_payload = {
    "timestamp": timestamp,
    "arviz_version": az.__version__,
    "python_version": sys.version,
    "traces_loaded": list(loaded.keys()),
    "waic": {
        "raw": waic_results,
        "comparison": waic_comparison,
    },
    "loo": {
        "raw": loo_results,
        "comparison": loo_comparison,
    },
    "manuscript_values": {
        "note": "Expected values from manuscript Table S5",
        "full_model_weight": 0.89,
        "linear_model_weight": 0.08,
        "null_model_weight": 0.03,
        "delta_waic_vs_null": 14.15,
    }
}

json_path = OUTPUT_DIR / f"waic_results_{timestamp}.json"
with open(json_path, "w") as f:
    json.dump(results_payload, f, indent=2)
print(f"\n[{datetime.now():%Y-%m-%d %H:%M:%S}] JSON saved: {json_path}")

# ── Human-readable summary ────────────────────────────────────────────────────

summary_lines = [
    "=" * 70,
    "WAIC / LOO MODEL COMPARISON — noise_decorrelation_HIV",
    f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
    f"ArviZ: {az.__version__}",
    "=" * 70,
    "",
    "WAIC COMPARISON (higher ELPD = better fit):",
    f"{'Rank':<5} {'Model':<25} {'ELPD_WAIC':>12} {'SE':>8} {'ΔELPD':>10} {'Weight':>8}",
    "-" * 70,
]
for row in waic_comparison:
    warn = " !" if row["warning"] else ""
    summary_lines.append(
        f"{row['rank']:<5} {row['model']:<25} "
        f"{row['elpd_waic']:>12.3f} {row['se']:>8.3f} "
        f"{row['delta_elpd']:>10.3f} {row['weight']:>8.3f}{warn}"
    )

summary_lines += [
    "",
    "LOO COMPARISON:",
    f"{'Rank':<5} {'Model':<25} {'ELPD_LOO':>12} {'SE':>8} {'ΔELPD':>10} {'Weight':>8}",
    "-" * 70,
]
for row in loo_comparison:
    warn = " !" if row["warning"] else ""
    summary_lines.append(
        f"{row['rank']:<5} {row['model']:<25} "
        f"{row['elpd_loo']:>12.3f} {row['se']:>8.3f} "
        f"{row['delta_elpd']:>10.3f} {row['weight']:>8.3f}{warn}"
    )

summary_lines += [
    "",
    "NOTE: ! = WAIC/LOO warning (Pareto k > 0.7 in LOO, or p_waic > 0.4 in WAIC)",
    "      These do not invalidate results but suggest posterior predictive",
    "      checks should be inspected.",
    "",
    f"Full results: {json_path}",
    "=" * 70,
]

summary_text = "\n".join(summary_lines)
print("\n" + summary_text)

summary_path = OUTPUT_DIR / f"waic_summary_{timestamp}.txt"
with open(summary_path, "w") as f:
    f.write(summary_text)
print(f"\n[{datetime.now():%Y-%m-%d %H:%M:%S}] Summary saved: {summary_path}")
print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] WAIC comparison COMPLETE.")