"""
Valcour CSF viral load analysis and calibration targets for CNS seeding

Purpose
-------
Read the curated Valcour individual-level dataset, extract CSF viral load
trajectories, compute early detectability metrics, and export calibration
targets for the CNS Monte Carlo seeding model.

Inputs (defaults assume repo layout):
  - quantum/quantum/data/curated/absolute/valcour_2015_full.xlsx

Assumptions
-----------
- Limit of detection (LOD) for CSF: 50 copies/mL.
- Empty/NaN CSF VL or log(CSF VL) entries are treated as TND (target not detected).
- If `cVL` == 1 or `log cVL` == 0 (common encodings for TND), treat as non‑detect.

Outputs
-------
- Tidy CSV: {outdir}/valcour_csf_tidy.csv with columns:
    subject_id, time_weeks, time_days, csf_vl, csf_detect, plasma_vl, arm
- Figures (if --figs):
    {outdir}/csf_spaghetti.(png|pdf)
    {outdir}/csf_detect_by_week.(png|pdf)
    {outdir}/csf_time_to_detect_ecdf.(png|pdf)
- Calibration targets:
    {outdir}/mc_calibration_targets.json
    {outdir}/mc_calibration_targets.csv

CLI
---
python -m quantum.quantum.utils.valcour_csf_vl_analysis \
  --excel quantum/quantum/data/curated/absolute/valcour_2015_full.xlsx \
  --outdir data/analysis_outputs/valcour_csf \
  --figs --pdf

"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Dict, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HOURS_PER_DAY = 24.0
LOD_CSF = 50.0


@dataclass
class Columns:
    subj: str = "id"
    time_weeks: str = "weeks"
    arm: str = "arm"
    csf_vl: str = "cVL"
    csf_log: str = "log cVL"
    plasma_vl: str = "pVL"


def load_valcour_excel(path: str | Path, cols: Columns = Columns()) -> pd.DataFrame:
    """
    Load the curated Excel and return a tidy DataFrame with standardized columns.
    Handles TND encodings.
    """
    path = Path(path)
    # Read first sheet by default
    df = pd.read_excel(path)

    # Heuristics for column names if exact labels differ slightly
    def find_col(preferred: str, alts: Tuple[str, ...]) -> str:
        if preferred in df.columns:
            return preferred
        for a in alts:
            if a in df.columns:
                return a
        # attempt case-insensitive
        lower_map = {c.lower(): c for c in df.columns}
        if preferred.lower() in lower_map:
            return lower_map[preferred.lower()]
        for a in alts:
            if a.lower() in lower_map:
                return lower_map[a.lower()]

        # Fallback for ID if not found: use the first column if no others match
        if preferred == "id" or "subj" in preferred.lower():
            return df.columns[0]

        raise KeyError(f"Could not find column for {preferred} among {df.columns.tolist()}")

    subj_col = find_col(cols.subj, ("subject", "Subject", "ID", "Subj"))
    weeks_col = find_col(cols.time_weeks, ("week", "Week", "weeks_since_enroll", "weeks"))
    arm_col = find_col(cols.arm, ("Arm", "ArmC"))
    csf_vl_col = find_col(cols.csf_vl, ("CSF VL", "CSF_VL", "csf_vl", "cVL"))
    csf_log_col = find_col(cols.csf_log, ("log_cVL", "log_csf", "Log cVL", "logcVL"))
    pvl_col = find_col(cols.plasma_vl, ("plasma", "Plasma VL", "plasma_vl", "pVL", "logpVL"))

    out = pd.DataFrame(
        {
            "subject_id": df[subj_col],
            "time_weeks": df[weeks_col],
            "arm": df[arm_col],
            "csf_vl": df[csf_vl_col],
            "csf_log": df[csf_log_col],
            "plasma_vl": df[pvl_col],
        }
    )

    # Normalize numeric types
    for c in ["time_weeks", "csf_vl", "csf_log", "plasma_vl"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    # TND rules: empty/NaN => non-detect; csf_vl==1 or csf_log==0 also common encodings
    tnd_mask = (
            out["csf_vl"].isna()
            | out["csf_log"].isna()
            | (out["csf_vl"] == 1)
            | (out["csf_log"] == 0)
    )
    # If csf_vl present but < LOD, treat as non‑detect
    below_lod = out["csf_vl"].fillna(0) < LOD_CSF
    csf_detect = ~(tnd_mask | below_lod)
    out["csf_detect"] = csf_detect.astype(int)

    # Time in days (approx 7*weeks where weeks may be fractional)
    out["time_days"] = out["time_weeks"].astype(float) * 7.0

    # Drop rows with no timing information
    out = out.dropna(subset=["time_weeks"]).reset_index(drop=True)
    return out


def proportion_detectable_by_weeks(df: pd.DataFrame, week_bins=(0, 1, 2, 4, 12, 24)) -> pd.DataFrame:
    bins = np.asarray(week_bins, dtype=float)
    # Assign each row to the rightmost bin edge <= time_weeks
    labels = [str(w) for w in bins]
    # Use pandas cut with right=True and include_lowest=True; define bin edges with -inf to each target week
    edges = [-np.inf] + list(bins)
    cats = pd.cut(df["time_weeks"], bins=edges, labels=labels, include_lowest=True, right=True)
    grp = df.assign(bin=cats).groupby("bin", dropna=False)
    out = grp["csf_detect"].agg(["mean", "count"]) \
        .rename(columns={"mean": "prop_detect", "count": "n"}).reset_index()
    # 95% binomial CI (Wald with guard)
    z = 1.96
    p = out["prop_detect"].values
    n = np.maximum(out["n"].values.astype(float), 1.0)
    se = np.sqrt(np.maximum(p * (1 - p), 0) / n)
    out["ci_low"] = np.clip(p - z * se, 0, 1)
    out["ci_high"] = np.clip(p + z * se, 0, 1)
    return out


def time_to_first_detect(df: pd.DataFrame) -> pd.DataFrame:
    # For each subject, find earliest time_days with csf_detect==1
    first = (
        df.sort_values(["subject_id", "time_days"]) \
            .query("csf_detect == 1") \
            .groupby("subject_id", as_index=False)["time_days"].first()
    )
    first = first.rename(columns={"time_days": "time_to_detect_days"})
    return first


def summarize_targets(first_detect: pd.DataFrame, horizons_days=(7, 14, 28)) -> Dict[str, Any]:
    t = first_detect["time_to_detect_days"].values
    out: Dict[str, Any] = {}
    for H in horizons_days:
        key = f"P_CSF_detect_<=_{int(H)}d"
        out[key] = float(np.mean(t <= H)) if len(t) else float("nan")
    if len(t):
        out.update(
            {
                "n_subjects_with_detect": int(np.isfinite(t).sum()),
                "median_time_to_detect_days": float(np.nanmedian(t)),
                "iqr_time_to_detect_days": [
                    float(np.nanpercentile(t, 25)),
                    float(np.nanpercentile(t, 75)),
                ],
            }
        )
    return out


def plot_spaghetti(df: pd.DataFrame, save_base: Path, pdf: bool = False) -> None:
    plt.figure(figsize=(9, 5))
    # Spaghetti plot using log10 scale for VL
    # Use detected values only for visualization; TND shown at y=LOD/2 as markers
    detected = df[df["csf_detect"] == 1]
    tnd = df[df["csf_detect"] == 0]
    for sid, grp in detected.groupby("subject_id"):
        y = np.log10(np.clip(grp["csf_vl"].values, LOD_CSF, None))
        plt.plot(grp["time_days"].values, y, alpha=0.5)
    # Plot TND markers
    if not tnd.empty:
        plt.scatter(tnd["time_days"], np.log10(LOD_CSF / 2.0), c="#999", s=10, label="TND")
    plt.xlabel("Time (days)")
    plt.ylabel("log10 CSF VL (copies/mL)")
    plt.title("Valcour: CSF viral load trajectories (detected values)")
    plt.grid(True, alpha=0.3)
    out_png = save_base.with_suffix(".png")
    plt.tight_layout();
    plt.savefig(out_png, dpi=300);
    plt.close()
    if pdf:
        out_pdf = save_base.with_suffix(".pdf")
        plt.figure(figsize=(9, 5))
        for sid, grp in detected.groupby("subject_id"):
            y = np.log10(np.clip(grp["csf_vl"].values, LOD_CSF, None))
            plt.plot(grp["time_days"].values, y, alpha=0.5)
        if not tnd.empty:
            plt.scatter(tnd["time_days"], np.log10(LOD_CSF / 2.0), c="#999", s=10, label="TND")
        plt.xlabel("Time (days)")
        plt.ylabel("log10 CSF VL (copies/mL)")
        plt.title("Valcour: CSF viral load trajectories (detected values)")
        plt.grid(True, alpha=0.3)
        plt.tight_layout();
        plt.savefig(out_pdf);
        plt.close()


def plot_detect_by_week(prop_df: pd.DataFrame, save_base: Path, pdf: bool = False) -> None:
    plt.figure(figsize=(8, 5))
    x = prop_df["bin"].astype(str)
    y = prop_df["prop_detect"].values
    lo = prop_df["ci_low"].values
    hi = prop_df["ci_high"].values
    plt.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o-", capsize=4)
    plt.ylim(0, 1)
    plt.ylabel("Proportion CSF detectable")
    plt.xlabel("Week bin (≤)")
    plt.title("Valcour: CSF detectability by week bin (LOD 50)")
    plt.grid(True, axis="y", alpha=0.3)
    out_png = save_base.with_suffix(".png")
    plt.tight_layout();
    plt.savefig(out_png, dpi=300);
    plt.close()
    if pdf:
        out_pdf = save_base.with_suffix(".pdf")
        plt.figure(figsize=(8, 5))
        plt.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o-", capsize=4)
        plt.ylim(0, 1)
        plt.ylabel("Proportion CSF detectable")
        plt.xlabel("Week bin (≤)")
        plt.title("Valcour: CSF detectability by week bin (LOD 50)")
        plt.grid(True, axis="y", alpha=0.3)
        plt.tight_layout();
        plt.savefig(out_pdf);
        plt.close()


def plot_time_to_detect_ecdf(first_df: pd.DataFrame, save_base: Path, pdf: bool = False) -> None:
    t = np.sort(first_df["time_to_detect_days"].values)
    y = np.arange(1, len(t) + 1) / len(t) if len(t) else np.array([])
    plt.figure(figsize=(7, 5))
    if len(t):
        plt.step(t, y, where="post")
    plt.xlabel("Time to first CSF detect (days)")
    plt.ylabel("ECDF")
    plt.title("Valcour: ECDF of time to first CSF detect")
    plt.grid(True, alpha=0.3)
    out_png = save_base.with_suffix(".png")
    plt.tight_layout();
    plt.savefig(out_png, dpi=300);
    plt.close()
    if pdf:
        out_pdf = save_base.with_suffix(".pdf")
        plt.figure(figsize=(7, 5))
        if len(t):
            plt.step(t, y, where="post")
        plt.xlabel("Time to first CSF detect (days)")
        plt.ylabel("ECDF")
        plt.title("Valcour: ECDF of time to first CSF detect")
        plt.grid(True, alpha=0.3)
        plt.tight_layout();
        plt.savefig(out_pdf);
        plt.close()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Valcour CSF VL analysis and calibration targets")
    ap.add_argument("--excel", default="quantum/quantum/data/curated/absolute/valcour_2015_full.xlsx",
                    help="Path to Valcour Excel file")
    ap.add_argument("--outdir", default="data/analysis_outputs/valcour_csf",
                    help="Directory to write outputs")
    ap.add_argument("--figs", action="store_true", help="Generate figures")
    ap.add_argument("--pdf", action="store_true", help="Also save PDF copies of figures")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_valcour_excel(args.excel)
    tidy_csv = outdir / "valcour_csf_tidy.csv"
    df.to_csv(tidy_csv, index=False)

    # Descriptives
    prop_df = proportion_detectable_by_weeks(df)
    first_df = time_to_first_detect(df)

    # Export summaries
    prop_csv = outdir / "csf_detect_by_week.csv"
    prop_df.to_csv(prop_csv, index=False)
    first_csv = outdir / "csf_time_to_detect.csv"
    first_df.to_csv(first_csv, index=False)

    targets = summarize_targets(first_df)
    targets_csv = outdir / "mc_calibration_targets.csv"
    targets_json = outdir / "mc_calibration_targets.json"
    # CSV
    pd.DataFrame({"metric": list(targets.keys()), "value": list(targets.values())}).to_csv(targets_csv, index=False)
    # JSON
    with open(targets_json, "w") as f:
        json.dump(targets, f, indent=2)

    # Print key numbers to console
    print("\n=== Valcour CSF detectability calibration targets ===")
    for k in sorted(targets.keys()):
        print(f"{k}: {targets[k]}")
    print(f"Saved tidy CSV: {tidy_csv}")

    if args.figs:
        plot_spaghetti(df, outdir / "csf_spaghetti", pdf=args.pdf)
        plot_detect_by_week(prop_df, outdir / "csf_detect_by_week", pdf=args.pdf)
        plot_time_to_detect_ecdf(first_df, outdir / "csf_time_to_detect_ecdf", pdf=args.pdf)
        print(f"Saved figures under: {outdir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())