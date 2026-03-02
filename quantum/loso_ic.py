#!/usr/bin/env python3
"""
loso_ic.py — Leave-One-Study-Out (LOSO) sensitivity analysis.

Iteratively excludes each contributing study and re-fits the Bayesian
v3.6 model to verify that no single cohort disproportionately drives
the overall conclusions.

Studies excluded (one at a time):
  - Valcour 2015      (via --exclude-valcour)
  - Chang 2002        (via --exclude-source)
  - Mohamed 2010      (via --exclude-source)
  - Young 2014        (via --exclude-source)
  - Sailasuta 2016    (via --exclude-source)
  - Sailasuta 2012    (via --exclude-source)

Outputs:
  - loso_summary.csv          Parameter estimates per exclusion
  - SuppFig_loso_analysis.png Three-panel forest plot (ξ_acute, ξ_chronic, β_ξ)

Usage:
    python quantum/quantum/utils/loso_ic.py
    python quantum/quantum/utils/loso_ic.py --draws 1000 --tune 500   # faster
    python quantum/quantum/utils/loso_ic.py --reuse-full              # skip re-running full model

Noise Decorrelation HIV — Demidont AC (2026)
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Studies and the CLI flags needed to exclude each one.
# "label" is the display name for the forest plot (matches supplementary fig).
# "flags" are the extra CLI args passed to bayesian_v3_6_corrected_local.py.
STUDIES = [
    {"label": "Valcour 2015",    "tag": "loso_no_valcour",       "flags": ["--exclude-valcour"]},
    {"label": "Chang 2002",      "tag": "loso_no_chang",         "flags": ["--exclude-source", "Chang"]},
    {"label": "Mohamed 2010",    "tag": "loso_no_mohamed",       "flags": ["--exclude-source", "Mohamed"]},
    {"label": "Young 2014",      "tag": "loso_no_young",         "flags": ["--exclude-source", "Young"]},
    {"label": "Sailasuta 2016",  "tag": "loso_no_sailasuta2016", "flags": ["--exclude-source", "Sailasuta et al. 2016"]},
    {"label": "Sailasuta 2012",  "tag": "loso_no_sailasuta2012", "flags": ["--exclude-source", "Sailasuta_2012"]},
]

# Parameters to extract from the results CSV
PARAMS_OF_INTEREST = ["ξ_acute", "ξ_chronic", "β_ξ"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_project_root() -> Path:
    """Walk up from this file until we find reproduce_all.py or .git."""
    p = Path(__file__).resolve().parent
    for _ in range(10):
        if (p / "reproduce_all.py").exists() or (p / ".git").exists():
            return p
        p = p.parent
    # Fallback: assume CWD
    return Path.cwd()


def find_v36_script(root: Path) -> Path:
    """Locate bayesian_v3_6_corrected_local.py."""
    candidates = [
        root / "quantum" / "quantum" / "bayesian_v3_6_corrected_local.py",
        root / "quantum" / "bayesian_v3_6_corrected_local.py",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(
        "Cannot find bayesian_v3_6_corrected_local.py. "
        f"Searched: {[str(c) for c in candidates]}"
    )


def resolve_output_dir(root: Path) -> Path:
    """
    Determine where results CSVs land.
    Respects NOISE_RESULTS_V36_DIR env var (same as Makefile).
    """
    env_dir = os.environ.get("NOISE_RESULTS_V36_DIR")
    if env_dir:
        return Path(env_dir)
    return root / "results" / "bayesian_v3_6"


def run_model(
    py: str,
    script: Path,
    root: Path,
    tag: str,
    extra_flags: list[str],
    draws: int,
    tune: int,
    chains: int,
) -> bool:
    """Run bayesian_v3_6_corrected_local.py with given flags."""
    cmd = [
        py, str(script),
        "--tag", tag,
        "--draws", str(draws),
        "--tune", str(tune),
        "--chains", str(chains),
    ] + extra_flags

    print(f"  CMD: {' '.join(cmd)}")
    t0 = time.time()

    result = subprocess.run(
        cmd,
        cwd=str(root),
        capture_output=True,
        text=True,
        timeout=3600,
    )

    elapsed = time.time() - t0
    if result.returncode == 0:
        print(f"  ✅ PASS ({elapsed:.1f}s)")
        return True
    else:
        print(f"  ❌ FAIL (exit {result.returncode}, {elapsed:.1f}s)")
        # Print last few lines of stderr for debugging
        if result.stderr:
            for line in result.stderr.strip().split("\n")[-5:]:
                print(f"     {line}")
        return False


def find_results_csv(output_dir: Path, tag: str) -> Path | None:
    """
    Find the results_v3_6_ratio_scale_<tag>.csv file.
    Checks both the output_dir root and the latest run subdirectory.
    """
    # Direct match in output root
    direct = output_dir / f"results_v3_6_ratio_scale_{tag}.csv"
    if direct.exists():
        return direct

    # Search in run subdirectories (sorted newest first)
    runs_dir = output_dir / "runs"
    if runs_dir.exists():
        run_dirs = sorted(runs_dir.iterdir(), reverse=True)
        for rd in run_dirs:
            candidate = rd / f"results_v3_6_ratio_scale_{tag}.csv"
            if candidate.exists():
                return candidate

    # Also check CWD (some configs write there)
    cwd_match = Path(f"results_v3_6_ratio_scale_{tag}.csv")
    if cwd_match.exists():
        return cwd_match

    return None


def find_summary_csv(output_dir: Path, tag: str) -> Path | None:
    """Find the ArviZ summary_<tag>.csv file."""
    direct = output_dir / f"summary_{tag}.csv"
    if direct.exists():
        return direct

    runs_dir = output_dir / "runs"
    if runs_dir.exists():
        for rd in sorted(runs_dir.iterdir(), reverse=True):
            candidate = rd / f"summary_{tag}.csv"
            if candidate.exists():
                return candidate

    cwd_match = Path(f"summary_{tag}.csv")
    if cwd_match.exists():
        return cwd_match

    return None


def parse_results(results_path: Path) -> dict:
    """
    Parse results_v3_6_ratio_scale_*.csv to extract parameter estimates.

    Expected format:
        Parameter,Mean,SD,HDI_2.5%,HDI_97.5%
        ξ_acute,0.607,...
        ξ_chronic,0.808,...
        Δξ,0.201,...
        β_ξ,-2.009,...
        P(ξ_acute < ξ_chronic),0.9105,...
    """
    params = {}
    with open(results_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["Parameter"]
            params[name] = {
                "mean": float(row["Mean"]),
                "sd": float(row["SD"]) if row.get("SD") else 0.0,
                "hdi_low": float(row["HDI_2.5%"]) if row.get("HDI_2.5%") else None,
                "hdi_high": float(row["HDI_97.5%"]) if row.get("HDI_97.5%") else None,
            }
    return params


def parse_summary_for_hdi(summary_path: Path) -> dict:
    """
    Parse the ArviZ summary CSV as a fallback for HDI values.
    First column is unnamed (parameter name).
    """
    params = {}
    with open(summary_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # First column may be unnamed or ''
            name = row.get("") or row.get("Unnamed: 0") or list(row.values())[0]
            try:
                params[name] = {
                    "mean": float(row.get("mean", 0)),
                    "sd": float(row.get("sd", 0)),
                    "hdi_low": float(row.get("hdi_2.5%", 0)),
                    "hdi_high": float(row.get("hdi_97.5%", 0)),
                }
            except (ValueError, TypeError):
                continue
    return params


def load_run_params(output_dir: Path, tag: str) -> dict | None:
    """Load parameter estimates from a completed run."""
    # Try results CSV first (has the key derived quantities)
    rpath = find_results_csv(output_dir, tag)
    if rpath:
        params = parse_results(rpath)
        # If HDI missing, supplement from summary CSV
        for pname in PARAMS_OF_INTEREST:
            if pname in params and params[pname].get("hdi_low") is None:
                spath = find_summary_csv(output_dir, tag)
                if spath:
                    sdata = parse_summary_for_hdi(spath)
                    if pname in sdata:
                        params[pname]["hdi_low"] = sdata[pname]["hdi_low"]
                        params[pname]["hdi_high"] = sdata[pname]["hdi_high"]
        return params

    # Fall back to summary CSV alone
    spath = find_summary_csv(output_dir, tag)
    if spath:
        return parse_summary_for_hdi(spath)

    return None


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_forest_plot(
    full_params: dict,
    loso_results: list[dict],
    output_path: Path,
):
    """
    Generate a 3-panel horizontal forest plot matching the style of
    SuppFig3_loso_analysis.png.

    Each panel shows one parameter (ξ_acute, ξ_chronic, β_ξ).
    Bottom row = Full Model (blue), grey bars = leave-one-out.
    Red dashed vertical line = full model posterior mean.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)

    # Build label list: leave-one-out studies + "Full Model"
    labels = [f"-{r['label']}" for r in loso_results] + ["Full Model"]
    n_rows = len(labels)
    y_pos = np.arange(n_rows)

    for ax, pname in zip(axes, PARAMS_OF_INTEREST):
        # Full model values
        full_mean = full_params[pname]["mean"]
        full_low = full_params[pname].get("hdi_low")
        full_high = full_params[pname].get("hdi_high")

        means = []
        lows = []
        highs = []
        colors = []

        for r in loso_results:
            p = r["params"].get(pname, {})
            m = p.get("mean", np.nan)
            lo = p.get("hdi_low")
            hi = p.get("hdi_high")

            means.append(m)
            # If HDI not available, use mean ± 1.96*sd as approximate 95% CI
            if lo is None or hi is None:
                sd = p.get("sd", 0)
                lo = m - 1.96 * sd
                hi = m + 1.96 * sd
            lows.append(lo)
            highs.append(hi)
            colors.append("0.6")  # grey

        # Append full model
        means.append(full_mean)
        if full_low is not None and full_high is not None:
            lows.append(full_low)
            highs.append(full_high)
        else:
            sd = full_params[pname].get("sd", 0)
            lows.append(full_mean - 1.96 * sd)
            highs.append(full_mean + 1.96 * sd)
        colors.append("#4472C4")  # blue

        means = np.array(means)
        lows = np.array(lows)
        highs = np.array(highs)

        # Error bars: asymmetric [mean - low, high - mean]
        xerr = np.array([means - lows, highs - means])

        ax.barh(
            y_pos, means, xerr=xerr, height=0.6,
            color=colors, edgecolor="0.4", linewidth=0.5,
            capsize=3, error_kw={"linewidth": 1.0},
        )

        # Reference line at full model mean
        ax.axvline(full_mean, color="red", linestyle="--", linewidth=1.2, alpha=0.8)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_xlabel(_param_xlabel(pname))
        ax.tick_params(axis="both", labelsize=9)

        # Clean up
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(
        "Leave-One-Study-Out (LOSO) Sensitivity Analysis",
        fontsize=13, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n✅ Saved LOSO forest plot: {output_path}")


def _param_xlabel(pname: str) -> str:
    labels = {
        "ξ_acute": "ξ_acute",
        "ξ_chronic": "ξ_chronic",
        "β_ξ": "β_ξ",
    }
    return labels.get(pname, pname)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Leave-One-Study-Out sensitivity analysis"
    )
    parser.add_argument(
        "--draws", type=int, default=1500,
        help="MCMC draws per chain (default: 1500)"
    )
    parser.add_argument(
        "--tune", type=int, default=1000,
        help="MCMC tuning steps (default: 1000)"
    )
    parser.add_argument(
        "--chains", type=int, default=4,
        help="Number of MCMC chains (default: 4)"
    )
    parser.add_argument(
        "--reuse-full", action="store_true",
        help="Skip re-running the full model (use existing results)"
    )
    parser.add_argument(
        "--figure-only", action="store_true",
        help="Only regenerate the figure from existing loso_summary.csv"
    )
    args = parser.parse_args()

    root = find_project_root()
    v36_script = find_v36_script(root)
    output_dir = resolve_output_dir(root)
    py = sys.executable

    print("=" * 60)
    print("LEAVE-ONE-STUDY-OUT (LOSO) SENSITIVITY ANALYSIS")
    print("=" * 60)
    print(f"Project root:  {root}")
    print(f"v3.6 script:   {v36_script}")
    print(f"Output dir:    {output_dir}")
    print(f"MCMC config:   draws={args.draws}, tune={args.tune}, chains={args.chains}")
    print()

    # Output paths
    loso_csv = output_dir / "loso_summary.csv"
    loso_fig = output_dir / "SuppFig_loso_analysis.png"

    # ------------------------------------------------------------------
    # Figure-only mode: just re-plot from existing CSV
    # ------------------------------------------------------------------
    if args.figure_only:
        if not loso_csv.exists():
            print(f"❌ Cannot find {loso_csv} — run without --figure-only first.")
            sys.exit(1)
        full_params, loso_results = _load_loso_csv(loso_csv)
        make_forest_plot(full_params, loso_results, loso_fig)
        sys.exit(0)

    # ------------------------------------------------------------------
    # Step 1: Run full model (or reuse existing)
    # ------------------------------------------------------------------
    full_tag = "with_valcour"
    print(f"{'='*60}")
    print(f"[1/{len(STUDIES)+1}] Full model")
    print(f"{'='*60}")

    if args.reuse_full:
        print("  (reusing existing full model results)")
    else:
        ok = run_model(py, v36_script, root, full_tag, [], args.draws, args.tune, args.chains)
        if not ok:
            print("❌ Full model failed — cannot proceed.")
            sys.exit(1)

    full_params = load_run_params(output_dir, full_tag)
    if full_params is None:
        print(f"❌ Cannot find results for tag '{full_tag}' in {output_dir}")
        sys.exit(1)

    print(f"  Full model: ξ_acute={full_params['ξ_acute']['mean']:.3f}, "
          f"ξ_chronic={full_params['ξ_chronic']['mean']:.3f}, "
          f"β_ξ={full_params['β_ξ']['mean']:.3f}")

    # ------------------------------------------------------------------
    # Step 2: Run each leave-one-out
    # ------------------------------------------------------------------
    loso_results = []
    failures = []

    for i, study in enumerate(STUDIES, start=2):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(STUDIES)+1}] Excluding: {study['label']}")
        print(f"{'='*60}")

        ok = run_model(
            py, v36_script, root,
            study["tag"], study["flags"],
            args.draws, args.tune, args.chains,
        )

        if not ok:
            failures.append(study["label"])
            # Still try to load if results exist from a prior run
            params = load_run_params(output_dir, study["tag"])
        else:
            params = load_run_params(output_dir, study["tag"])

        if params is None:
            print(f"  ⚠️  No results found for {study['label']}")
            # Use NaN placeholders so plotting still works
            params = {
                p: {"mean": np.nan, "sd": np.nan, "hdi_low": np.nan, "hdi_high": np.nan}
                for p in PARAMS_OF_INTEREST
            }

        loso_results.append({
            "label": study["label"],
            "tag": study["tag"],
            "params": params,
        })

        # Progress report
        for pname in PARAMS_OF_INTEREST:
            m = params.get(pname, {}).get("mean", np.nan)
            print(f"  {pname} = {m:.3f}" if not np.isnan(m) else f"  {pname} = N/A")

    # ------------------------------------------------------------------
    # Step 3: Save summary CSV
    # ------------------------------------------------------------------
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(loso_csv, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["Excluded_Study", "Tag"]
        for pname in PARAMS_OF_INTEREST:
            header += [f"{pname}_mean", f"{pname}_sd", f"{pname}_hdi_low", f"{pname}_hdi_high"]
        writer.writerow(header)

        # Full model row
        row = ["Full Model", full_tag]
        for pname in PARAMS_OF_INTEREST:
            p = full_params.get(pname, {})
            row += [
                f"{p.get('mean', ''):.4f}" if p.get('mean') is not None else "",
                f"{p.get('sd', ''):.4f}" if p.get('sd') is not None else "",
                f"{p.get('hdi_low', ''):.4f}" if p.get('hdi_low') is not None else "",
                f"{p.get('hdi_high', ''):.4f}" if p.get('hdi_high') is not None else "",
            ]
        writer.writerow(row)

        # LOSO rows
        for r in loso_results:
            row = [r["label"], r["tag"]]
            for pname in PARAMS_OF_INTEREST:
                p = r["params"].get(pname, {})
                m = p.get("mean")
                row += [
                    f"{m:.4f}" if m is not None and not np.isnan(m) else "",
                    f"{p.get('sd', ''):.4f}" if p.get('sd') is not None and not np.isnan(p.get('sd', np.nan)) else "",
                    f"{p.get('hdi_low', ''):.4f}" if p.get('hdi_low') is not None and not np.isnan(p.get('hdi_low', np.nan)) else "",
                    f"{p.get('hdi_high', ''):.4f}" if p.get('hdi_high') is not None and not np.isnan(p.get('hdi_high', np.nan)) else "",
                ]
            writer.writerow(row)

    print(f"\n✅ Saved LOSO summary: {loso_csv}")

    # ------------------------------------------------------------------
    # Step 4: Generate forest plot
    # ------------------------------------------------------------------
    make_forest_plot(full_params, loso_results, loso_fig)

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("LOSO ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"  Studies tested: {len(STUDIES)}")
    print(f"  Failures:       {len(failures)}")

    if failures:
        print(f"  Failed:         {', '.join(failures)}")

    # Quick robustness check: is β_ξ consistently < 0 across all exclusions?
    all_negative = all(
        r["params"].get("β_ξ", {}).get("mean", 0) < 0
        for r in loso_results
        if not np.isnan(r["params"].get("β_ξ", {}).get("mean", np.nan))
    )
    if all_negative:
        print("  ✅ β_ξ < 0 in ALL leave-one-out analyses — robust.")
    else:
        print("  ⚠️  β_ξ sign inconsistent — check individual results.")

    if failures:
        sys.exit(1)


def _load_loso_csv(csv_path: Path):
    """Reload results from a previously saved loso_summary.csv."""
    full_params = {}
    loso_results = []

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            params = {}
            for pname in PARAMS_OF_INTEREST:
                try:
                    params[pname] = {
                        "mean": float(row[f"{pname}_mean"]) if row.get(f"{pname}_mean") else np.nan,
                        "sd": float(row[f"{pname}_sd"]) if row.get(f"{pname}_sd") else np.nan,
                        "hdi_low": float(row[f"{pname}_hdi_low"]) if row.get(f"{pname}_hdi_low") else np.nan,
                        "hdi_high": float(row[f"{pname}_hdi_high"]) if row.get(f"{pname}_hdi_high") else np.nan,
                    }
                except (ValueError, KeyError):
                    params[pname] = {"mean": np.nan, "sd": np.nan, "hdi_low": np.nan, "hdi_high": np.nan}

            if row["Excluded_Study"] == "Full Model":
                full_params = params
            else:
                loso_results.append({
                    "label": row["Excluded_Study"],
                    "tag": row.get("Tag", ""),
                    "params": params,
                })

    return full_params, loso_results


if __name__ == "__main__":
    main()
