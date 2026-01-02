#!/usr/bin/env python3
"""
Absolute (Valcour) analysis pipeline

Loads curated Valcour absolute mM data (NAA, Cho) by Week and Region.
Treats all rows as Acute. Supports time-slice filtering and region selection.

Outputs:
- summary.csv: counts and basic stats by region and week
- results.csv: tidy filtered dataset (one row per observation)
- manifest.json: run metadata (args, inputs, filters)

Note: This initial version focuses on robust loading and summarization, with
hooks for Bayesian modeling to be added in a subsequent step if desired.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Any, List
import pandas as pd

from .pipeline_utils import (
    list_curated_files,
    parse_time_slice,
    parse_regions,
    write_run_manifest,
)


EXPECTED_COLS = {
    "study": ["Study"],
    "subject": ["SubjectID", "Subject", "ID"],
    "week": ["Week", "Weeks", "Valcour_Week"],
    "region": ["Region", "ROI", "BrainRegion"],
    "naa": ["NAA", "Naa"],
    "cho": ["Cho", "CHO"],
    "vl": ["VL", "ViralLoad", "Viral_Load"],
    "logvl": ["logVL", "LogVL", "log10VL", "log_VL"],
}


def _standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    colmap: Dict[str, str] = {}
    low = {c.lower(): c for c in df.columns}
    for key, candidates in EXPECTED_COLS.items():
        for cand in candidates:
            c0 = cand
            cl = cand.lower()
            if cl in low:
                colmap[low[cl]] = key
                break
            if c0 in df.columns:
                colmap[c0] = key
                break
    out = df.rename(columns=colmap).copy()
    # Keep only standardized columns if present
    keep = [c for c in ["study", "subject", "week", "region", "naa", "cho", "vl", "logvl"] if c in out.columns]
    return out[keep]


def _load_absolute_dir(absolute_dir: Path) -> pd.DataFrame:
    files = list_curated_files(absolute_dir)
    frames: List[pd.DataFrame] = []
    for f in files:
        try:
            if f.suffix.lower() == ".csv":
                df = pd.read_csv(f)
            elif f.suffix.lower() in {".xls", ".xlsx"}:
                df = pd.read_excel(f)
            else:
                continue
            df = _standardize_columns(df)
            if "study" not in df.columns:
                df["study"] = f.stem
            frames.append(df)
        except Exception:
            continue
    if not frames:
        return pd.DataFrame(columns=["study", "subject", "week", "region", "naa", "cho", "vl", "logvl"])
    out = pd.concat(frames, ignore_index=True)
    # Coerce types
    if "week" in out.columns:
        out["week"] = pd.to_numeric(out["week"], errors="coerce").astype("Int64")
    for m in ("naa", "cho", "vl", "logvl"):
        if m in out.columns:
            out[m] = pd.to_numeric(out[m], errors="coerce")
    if "region" in out.columns:
        out["region"] = out["region"].astype(str)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Absolute (Valcour) pipeline")
    this_dir = Path(__file__).resolve().parent
    default_abs = this_dir / "data" / "curated" / "absolute"
    default_results = this_dir / "results" / "absolute"

    parser.add_argument("--absolute-dir", type=str, default=str(default_abs))
    parser.add_argument("--results-dir", type=str, default=str(default_results))
    parser.add_argument("--regions", type=str, default="BG", help="Comma list or 'all'")
    parser.add_argument("--time-slice", type=str, default="all")
    parser.add_argument("--include-vl", action="store_true", help="Include VL/logVL in results table if present")
    parser.add_argument("--save-trace", action="store_true", help="Reserved for future Bayesian modeling")

    args = parser.parse_args(argv)

    absolute_dir = Path(args.absolute_dir)
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    region_set = parse_regions(args.regions)
    weeks = parse_time_slice(args.time_slice)

    df = _load_absolute_dir(absolute_dir)

    # Filters
    if region_set is not None and "region" in df.columns:
        df = df[df["region"].str.upper().isin(region_set)].copy()
    if weeks is not None and "week" in df.columns:
        df = df[df["week"].isin(weeks)].copy()

    # Treat all as Acute — annotate for downstream consistency
    df["group"] = "Acute"

    # Prepare outputs
    # results.csv: tidy filtered data with selected columns only
    cols = [c for c in ["study", "subject", "week", "region", "naa", "cho", "vl", "logvl", "group"] if c in df.columns]
    results_path = results_dir / "results.csv"
    df[cols].to_csv(results_path, index=False)

    # summary.csv: counts per region/week and basic stats by region
    summaries: List[pd.DataFrame] = []
    if not df.empty:
        # Counts by region and week
        if "region" in df.columns and "week" in df.columns:
            counts = df.groupby(["region", "week"]).size().reset_index(name="n")
            counts["kind"] = "counts_by_region_week"
            summaries.append(counts)
        # Stats by region for NAA and Cho
        for metric in ["naa", "cho"]:
            if metric in df.columns:
                stats = (
                    df.groupby("region")[metric]
                    .agg(["count", "mean", "std", "min", "max"])
                    .reset_index()
                )
                stats.insert(1, "metric", metric.upper())
                stats["kind"] = "stats_by_region"
                summaries.append(stats)
    summary_df = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    summary_path = results_dir / "summary.csv"
    summary_df.to_csv(summary_path, index=False)

    # Run manifest
    manifest_info: Dict[str, Any] = {
        "pipeline": "absolute",
        "absolute_dir": str(absolute_dir),
        "results_dir": str(results_dir),
        "regions": sorted(list(region_set)) if region_set is not None else "all",
        "time_slice": args.time_slice,
        "inputs": [str(p) for p in list_curated_files(absolute_dir)],
        "row_count": int(len(df)),
        "columns": cols,
    }
    write_run_manifest(results_dir, manifest_info)

    print(f"✅ Absolute pipeline complete. Rows: {len(df)}")
    print(f"   Results: {results_path}")
    print(f"   Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
