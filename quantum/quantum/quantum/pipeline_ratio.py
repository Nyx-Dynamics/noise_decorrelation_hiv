#!/usr/bin/env python3
"""
Ratio triad analysis pipeline

Loads curated ratio datasets with groups (Primary/Acute, Chronic, Control)
and metabolite ratios (NAA/Cr, Cho/Cr). Regions may be present or absent.

Outputs:
- summary.csv: group counts and descriptive stats per metric (and region if available)
- results.csv: tidy normalized long-format data
- manifest.json: run metadata (args, inputs, filters)

This initial version focuses on robust loading/normalization and descriptive
inference (Δ = Acute − Chronic, and an empirical P(Acute < Chronic) via
bootstrap). Hooks for Bayesian modeling can be added later if needed.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Any, List, Optional
import pandas as pd
import numpy as np

from .pipeline_utils import (
    list_curated_files,
    parse_regions,
    write_run_manifest,
)


GROUP_ALIASES = {
    "acute": "Acute",
    "primary": "Acute",
    "chronic": "Chronic",
    "control": "Control",
    "healthy": "Control",
}


def _standardize_ratio_columns(df: pd.DataFrame) -> pd.DataFrame:
    # Try to find group, region, and metric columns
    colmap: Dict[str, str] = {}
    cols = {c.lower(): c for c in df.columns}

    # Group column
    for cand in ("group", "phase", "cohort", "status"):
        if cand in cols:
            colmap[cols[cand]] = "group"
            break

    # Region optional
    for cand in ("region", "roi", "brainregion"):
        if cand in cols:
            colmap[cols[cand]] = "region"
            break

    # Metric columns: prefer explicit NAA/Cr, Cho/Cr
    metric_cols: List[str] = []
    for m in ["naa/cr", "naa_cr", "naacr", "naa_cr_ratio", "naa_to_cr" ]:
        if m in cols:
            metric_cols.append(cols[m])
            colmap[cols[m]] = "NAA_Cr"
            break
    for m in ["cho/cr", "cho_cr", "chocr", "cho_cr_ratio", "cho_to_cr" ]:
        if m in cols:
            metric_cols.append(cols[m])
            colmap[cols[m]] = "Cho_Cr"
            break

    # Fallback: if no obvious metric columns but there are columns named like NAA_Cr etc.
    for m in df.columns:
        ml = m.lower().replace(" ", "").replace("-", "_")
        if ml == "naa_cr" and m not in metric_cols:
            metric_cols.append(m)
            colmap[m] = "NAA_Cr"
        if ml == "cho_cr" and m not in metric_cols:
            metric_cols.append(m)
            colmap[m] = "Cho_Cr"

    out = df.rename(columns=colmap).copy()

    # Normalize group labels
    if "group" in out.columns:
        out["group"] = out["group"].astype(str).str.strip().str.lower().map(GROUP_ALIASES).fillna(out["group"])  # type: ignore
    else:
        out["group"] = np.nan

    # Melt metrics to long format
    long_frames: List[pd.DataFrame] = []
    for metric in ("NAA_Cr", "Cho_Cr"):
        if metric in out.columns:
            tmp = out[[c for c in out.columns if c in {"group", "region", metric}]].copy()
            tmp = tmp.rename(columns={metric: "value"})
            tmp["metric"] = metric
            long_frames.append(tmp)

    if not long_frames:
        # No recognizable metric columns; return empty standardized frame
        return pd.DataFrame(columns=["group", "region", "metric", "value"])  # type: ignore

    long = pd.concat(long_frames, ignore_index=True)
    # Clean types
    long["value"] = pd.to_numeric(long["value"], errors="coerce")
    if "region" in long.columns:
        long["region"] = long["region"].astype(str)
    else:
        long["region"] = "All"

    # Drop rows without group or value
    long = long.dropna(subset=["group", "value"]).copy()
    return long


def _bootstrap_prob_acute_lt_chronic(a: np.ndarray, c: np.ndarray, rng: np.random.Generator, B: int = 10000) -> float:
    if len(a) == 0 or len(c) == 0:
        return np.nan
    a = a[~np.isnan(a)]
    c = c[~np.isnan(c)]
    if len(a) == 0 or len(c) == 0:
        return np.nan
    a_idx = rng.integers(0, len(a), size=(B, len(a)))
    c_idx = rng.integers(0, len(c), size=(B, len(c)))
    a_means = a[a_idx].mean(axis=1)
    c_means = c[c_idx].mean(axis=1)
    return float((a_means < c_means).mean())


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Ratio triad pipeline")
    this_dir = Path(__file__).resolve().parent
    default_ratio = this_dir / "data" / "curated" / "ratio"
    default_results = this_dir / "results" / "ratio"

    parser.add_argument("--ratio-dir", type=str, default=str(default_ratio))
    parser.add_argument("--results-dir", type=str, default=str(default_results))
    parser.add_argument("--metabolites", type=str, choices=["NAA_Cr", "Cho_Cr", "both"], default="both")
    parser.add_argument("--regions", type=str, default="all", help="Region filter (comma) if Region is present")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--save-trace", action="store_true", help="Reserved for future Bayesian modeling")

    args = parser.parse_args(argv)
    ratio_dir = Path(args.ratio_dir)
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    files = list_curated_files(ratio_dir)
    frames: List[pd.DataFrame] = []
    for f in files:
        try:
            if f.suffix.lower() == ".csv":
                df = pd.read_csv(f)
            elif f.suffix.lower() in {".xls", ".xlsx"}:
                df = pd.read_excel(f)
            else:
                continue
            df_long = _standardize_ratio_columns(df)
            if not df_long.empty:
                df_long["study"] = f.stem
                frames.append(df_long)
        except Exception:
            continue

    if frames:
        data = pd.concat(frames, ignore_index=True)
    else:
        data = pd.DataFrame(columns=["study", "group", "region", "metric", "value"])  # type: ignore

    # Filter metrics
    if args.metabolites != "both":
        data = data[data["metric"] == args.metabolites].copy()

    # Filter regions if provided (only applies if Region column exists)
    region_set = parse_regions(args.regions)
    if region_set is not None and "region" in data.columns:
        data = data[data["region"].str.upper().isin(region_set)].copy()

    # Save tidy results
    results_path = results_dir / "results.csv"
    data.to_csv(results_path, index=False)

    # Build summary tables and empirical P(Acute < Chronic)
    rng = np.random.default_rng(args.seed)
    summaries: List[pd.DataFrame] = []

    if not data.empty:
        group_stats = (
            data.groupby(["metric", "region", "group"])  # type: ignore
            ["value"].agg(["count", "mean", "std", "min", "max"]).reset_index()
        )
        group_stats["kind"] = "group_stats"
        summaries.append(group_stats)

        # For each metric/region, compute Δ and probability
        recs: List[Dict[str, Any]] = []
        for (metric, region), sub in data.groupby(["metric", "region"]):  # type: ignore
            acute = sub[sub["group"] == "Acute"]["value"].to_numpy()
            chronic = sub[sub["group"] == "Chronic"]["value"].to_numpy()
            if len(acute) == 0 or len(chronic) == 0:
                p = np.nan
                d = np.nan
            else:
                p = _bootstrap_prob_acute_lt_chronic(acute, chronic, rng)
                d = float(np.nanmean(acute) - np.nanmean(chronic))
            recs.append({
                "metric": metric,
                "region": region,
                "delta_acute_minus_chronic": d,
                "P(acute<chronic)": p,
                "n_acute": int(np.isfinite(acute).sum()),
                "n_chronic": int(np.isfinite(chronic).sum()),
            })
        comp_df = pd.DataFrame(recs)
        comp_df["kind"] = "acute_vs_chronic"
        summaries.append(comp_df)

    summary_df = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    summary_path = results_dir / "summary.csv"
    summary_df.to_csv(summary_path, index=False)

    # Run manifest
    manifest_info: Dict[str, Any] = {
        "pipeline": "ratio",
        "ratio_dir": str(ratio_dir),
        "results_dir": str(results_dir),
        "metabolites": args.metabolites,
        "regions": sorted(list(region_set)) if region_set is not None else "all",
        "inputs": [str(p) for p in files],
        "row_count": int(len(data)),
        "columns": list(data.columns),
    }
    write_run_manifest(results_dir, manifest_info)

    print(f"✅ Ratio pipeline complete. Rows: {len(data)}")
    print(f"   Results: {results_path}")
    print(f"   Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
