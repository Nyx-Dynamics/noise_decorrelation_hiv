#!/usr/bin/env python3
"""
Valcour 2015 absolute-units preparation (patient-level)
=======================================================

Purpose:
- Load Valcour 2015 individual-level CSV/XLSX files (absolute concentrations; mM or mmol/kg)
- Standardize schema to long-form patient-level table
- Keep absolute units for regional modeling pipelines (no conversion to ratios here)
- Emit a concise summary by region and group

Inputs (directory):
- /data/individual/VALCOUR_2015_*.csv (or .xlsx)
  The script is tolerant to multiple source schemas:
  - Long format with columns like: subject_id, group, region, metabolite, value, units
  - Wide format with columns like: NAA, Cho, Cr; will be melted to long format

Outputs:
- quantum/data/curated/valcour_abs.csv          (long-format absolute units)
- data/model_outputs/valcour_abs_summary.csv    (means/SDs by region×group×metabolite)

Notes:
- Units: We keep whatever absolute unit is present but recommend mmol/kg (or mM).
- The script will try to coerce numeric values and drop non-numeric rows.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np


def std_col(col: str) -> str:
    return col.strip().lower().replace(" ", "_").replace("-", "_").replace("/", "_")


def load_table(path: Path) -> pd.DataFrame | None:
    try:
        if path.suffix.lower() == ".csv":
            return pd.read_csv(path)
        elif path.suffix.lower() in {".xls", ".xlsx"}:
            return pd.read_excel(path)
    except Exception as e:
        print(f"[WARN] Failed to load {path}: {e}")
    return None


def detect_schema(df: pd.DataFrame) -> Dict[str, str]:
    cols = {std_col(c): c for c in df.columns}
    mapping: Dict[str, str] = {}

    # Identity/meta
    for key in ["subject_id", "subject", "id"]:
        if key in cols:
            mapping["subject_id"] = cols[key]
            break
    for key in ["group", "phase", "arm", "status"]:
        if key in cols:
            mapping["group"] = cols[key]
            break
    if "region" in cols:
        mapping["region"] = cols["region"]

    # Long-format metabolite columns
    for key in ["metabolite", "analyte", "met"]:
        if key in cols:
            mapping["metabolite"] = cols[key]
            break
    for key in ["value", "conc", "concentration", "amount"]:
        if key in cols:
            mapping["value"] = cols[key]
            break
    for key in ["units", "unit"]:
        if key in cols:
            mapping["units"] = cols[key]
            break

    # Wide-format metabolite columns
    if "naa" in cols:
        mapping["naa"] = cols["naa"]
    if "cho" in cols:
        mapping["cho"] = cols["cho"]
    if "cr" in cols:
        mapping["cr"] = cols["cr"]

    return mapping


def to_long(df: pd.DataFrame, mapping: Dict[str, str]) -> pd.DataFrame:
    df2 = df.copy()
    df2.columns = [std_col(c) for c in df2.columns]

    # Preferred: long format already present
    if "metabolite" in mapping and "value" in mapping:
        long_df = pd.DataFrame({
            "subject_id": df[mapping.get("subject_id")] if "subject_id" in mapping else np.arange(len(df)),
            "group": df[mapping.get("group")] if "group" in mapping else "unknown",
            "region": df[mapping.get("region")] if "region" in mapping else "GLOBAL",
            "metabolite": df[mapping["metabolite"]],
            "value": pd.to_numeric(df[mapping["value"]], errors="coerce"),
            "units": df[mapping["units"]] if "units" in mapping else "mmol/kg",
        })
        long_df = long_df.dropna(subset=["value"])  # keep numeric
        return long_df

    # Fallback: wide format → melt
    wide_cols = [(k, v) for k, v in mapping.items() if k in {"naa", "cho", "cr"}]
    if wide_cols:
        tmp = df[[v for _, v in wide_cols]].copy()
        tmp.columns = [k for k, _ in wide_cols]
        meta = pd.DataFrame({
            "subject_id": df[mapping.get("subject_id")] if "subject_id" in mapping else np.arange(len(df)),
            "group": df[mapping.get("group")] if "group" in mapping else "unknown",
            "region": df[mapping.get("region")] if "region" in mapping else "GLOBAL",
        })
        melted = tmp.melt(var_name="metabolite", value_name="value")
        melted["value"] = pd.to_numeric(melted["value"], errors="coerce")
        long_df = pd.concat([meta.loc[melted.index].reset_index(drop=True), melted.reset_index(drop=True)], axis=1)
        long_df["units"] = "mmol/kg"
        long_df = long_df.dropna(subset=["value"])  # keep numeric
        return long_df

    # If neither pattern found, return empty
    return pd.DataFrame(columns=["subject_id", "group", "region", "metabolite", "value", "units"])


def summarize(long_df: pd.DataFrame) -> pd.DataFrame:
    # Normalize text columns
    for c in ["group", "region", "metabolite", "units"]:
        if c in long_df.columns:
            long_df[c] = long_df[c].astype(str)
    grp_cols = [c for c in ["region", "group", "metabolite", "units"] if c in long_df.columns]
    if not grp_cols:
        return pd.DataFrame()
    g = long_df.groupby(grp_cols)["value"]
    out = g.agg(["count", "mean", "std"]).reset_index()
    return out


def main():
    ap = argparse.ArgumentParser(description="Prepare Valcour 2015 patient-level absolute datasets (long format + summary)")
    ap.add_argument("--individual-dir", required=True, help="Directory containing VALCOUR_2015_*.csv/.xlsx files")
    ap.add_argument("--out-dir", required=True, help="Output directory for curated CSVs (e.g., quantum/data/curated)")
    ap.add_argument("--summary-dir", default=str(Path("/data/model_outputs")), help="Directory for summary CSV")
    args = ap.parse_args()

    in_dir = Path(args.individual_dir).expanduser()
    out_dir = Path(args.out_dir).expanduser()
    summary_dir = Path(args.summary_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    long_parts: List[pd.DataFrame] = []
    files = sorted(list(in_dir.glob("VALCOUR_2015_*")))
    if not files:
        print(f"[WARN] No Valcour files found under {in_dir}")

    for path in files:
        if path.suffix.lower() not in {".csv", ".xls", ".xlsx"}:
            continue
        df = load_table(path)
        if df is None or df.empty:
            print(f"[INFO] Skipping empty/unreadable: {path}")
            continue
        mapping = detect_schema(df)
        long_df = to_long(df, mapping)
        if long_df.empty:
            print(f"[INFO] Could not interpret schema for: {path}")
            continue
        # tag source
        long_df.insert(0, "source", str(path))
        long_parts.append(long_df)
        print(f"[OK] Parsed {path} → rows: {len(long_df)}")

    if long_parts:
        all_long = pd.concat(long_parts, ignore_index=True)
        # Basic cleaning: standardize group labels for Valcour
        def norm_group(x: str) -> str:
            s = str(x).strip().lower()
            if "control" in s or s in {"hiv-", "healthy"}:
                return "control"
            if "acute" in s or "phi" in s or "baseline" in s:
                return "acute"
            if "week" in s or "w24" in s:
                return "acute_week"
            if "chronic" in s:
                return "chronic"
            return x
        all_long["group"] = all_long["group"].map(norm_group)

        # Save curated long-form
        out_long = out_dir / "valcour_abs.csv"
        all_long.to_csv(out_long, index=False)
        print(f"[SAVE] {out_long} | rows: {len(all_long)}")

        # Save summary
        summ = summarize(all_long)
        if not summ.empty:
            out_summ = summary_dir / "valcour_abs_summary.csv"
            summ.to_csv(out_summ, index=False)
            print(f"[SAVE] {out_summ} | rows: {len(summ)}")
    else:
        print("[INFO] No long-form Valcour datasets produced")

    print("\nDone. Review curated CSVs under:")
    print(f"  - {out_dir}")
    print("Summaries under:")
    print(f"  - {summary_dir}")


if __name__ == "__main__":
    main()
