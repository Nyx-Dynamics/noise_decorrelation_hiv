#!/usr/bin/env python3
"""
Patient-level ratios preparation
================================

Purpose:
- Load patient-level extracted datasets from multiple studies (CSV/XLSX) under a root directory
- Detect input scale/reference per file (absolute mM vs metabolite-to-Cr ratios)
- Compute within-study (and optionally within-region) control means
- Produce two curated, unit-consistent patient-level datasets for modeling:
  1) Absolute-to-relative metrics (per metabolite): NAA_rel, Cho_rel, (optional Cr_rel)
  2) Metabolite-to-Cr ratios to relative: NAA_over_Cr_rel, Cho_over_Cr_rel
- Emit a summary report with per-study control means and counts

Outputs (CSV to avoid parquet dependency):
- quantum/data/curated/patient_abs_rel.csv
- quantum/data/curated/patient_naacr_rel.csv
- data/model_outputs/patient_ratios_summary.csv

Notes:
- This script is defensive. It will skip files that do not contain
  the columns needed for a given transformation. It logs what it finds.
- Region-specific control means are used if a region column exists; otherwise,
  study-level control means are used.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple


CONTROL_LABELS = {"control", "healthy", "hiv-", "hiv_neg", "hiv_negative"}
PHASE_MAP = {
    "control": "healthy",
    "healthy": "healthy",
    "hiv-": "healthy",
    "hiv_neg": "healthy",
    "hiv_negative": "healthy",
    "acute": "acute",
    "phi": "acute",
    "primary": "acute",
    "chronic": "chronic",
    "early": "chronic",
    "early_chronic": "chronic",
    "late": "chronic",
    "late_chronic": "chronic",
}


def std_col(col: str) -> str:
    """Normalize a column name for robust matching."""
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


def detect_columns(df: pd.DataFrame) -> Dict[str, str]:
    """Detect key columns and return a mapping from canonical name -> actual name."""
    cols = {std_col(c): c for c in df.columns}
    mapping: Dict[str, str] = {}

    # Identity/meta
    for key in ["subject_id", "subject", "id"]:
        if key in cols:
            mapping["subject_id"] = cols[key]
            break
    for key in ["study", "study_name", "source_study"]:
        if key in cols:
            mapping["study"] = cols[key]
            break
    for key in ["phase", "group", "status"]:
        if key in cols:
            mapping["phase_or_group"] = cols[key]
            break
    if "region" in cols:
        mapping["region"] = cols["region"]

    # Time metadata
    for key in ["days_since_infection", "duration_days", "days"]:
        if key in cols:
            mapping["days_since_infection"] = cols[key]
            break

    # Absolute metabolites
    for key in ["naa_mm", "naa_mmol_kg", "naa"]:
        if key in cols:
            mapping["naa_mM"] = cols[key]
            break
    for key in ["cho_mm", "cho_mmol_kg", "cho"]:
        if key in cols:
            mapping["cho_mM"] = cols[key]
            break
    for key in ["cr_mm", "cr_mmol_kg", "cr"]:
        if key in cols:
            mapping["cr_mM"] = cols[key]
            break

    # Ratios to Cr
    for key in ["naa_cr", "naa_over_cr", "naa_to_cr"]:
        if key in cols:
            mapping["naa_over_cr"] = cols[key]
            break
    for key in ["cho_cr", "cho_over_cr", "cho_to_cr"]:
        if key in cols:
            mapping["cho_over_cr"] = cols[key]
            break

    return mapping


def normalize_phase(val: str) -> str:
    if pd.isna(val):
        return "unknown"
    s = std_col(str(val))
    return PHASE_MAP.get(s, s)


def is_control_label(val: str) -> bool:
    if pd.isna(val):
        return False
    return std_col(str(val)) in CONTROL_LABELS


def compute_control_means(df: pd.DataFrame, group_cols: List[str], value_col: str) -> pd.DataFrame:
    grp = df.groupby(group_cols)[value_col]
    means = grp.mean().rename(f"control_mean_{value_col}")
    return means.reset_index()


def attach_rel(df: pd.DataFrame, value_col: str, control_means: pd.DataFrame, group_cols: List[str], out_col: str) -> pd.DataFrame:
    merged = df.merge(control_means, on=group_cols, how="left")
    with np.errstate(divide='ignore', invalid='ignore'):
        merged[out_col] = merged[value_col] / merged[f"control_mean_{value_col}"]
    merged.drop(columns=[f"control_mean_{value_col}"], inplace=True)
    return merged


def process_file(path: Path, study_id_counter: Dict[str, int]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (abs_rel_rows, naacr_rel_rows, summary_rows)."""
    df = load_table(path)
    if df is None or df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    # Standardize column names
    original_cols = df.columns.tolist()
    df_std = df.copy()
    df_std.columns = [std_col(c) for c in df.columns]

    mapping = detect_columns(df_std)

    # Require subject-level rows
    if "subject_id" not in mapping:
        # Attempt to synthesize subject_id if a column like 'subject' exists
        if "subject" in df_std.columns:
            mapping["subject_id"] = "subject"
        else:
            print(f"[INFO] Skipping (no subject_id): {path}")
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    # Study name
    study_name = df_std[mapping.get("study", None)].iloc[0] if "study" in mapping else path.stem
    study_name = str(study_name)
    if study_name not in study_id_counter:
        study_id_counter[study_name] = len(study_id_counter)
    study_id = study_id_counter[study_name]

    # Phase/group
    if "phase_or_group" in mapping:
        df_std["phase"] = df_std[mapping["phase_or_group"]].map(normalize_phase)
        df_std["is_control"] = df_std[mapping["phase_or_group"]].map(is_control_label)
    else:
        df_std["phase"], df_std["is_control"] = "unknown", False

    # Region handling
    has_region = "region" in mapping
    if has_region:
        df_std["region"] = df_std[mapping["region"]].astype(str)

    group_cols = ["study_name"] + (["region"] if has_region else [])

    df_std["study_name"] = study_name
    df_std["study_id"] = study_id
    df_std["subject_id"] = df_std[mapping["subject_id"]].astype(str)
    if "days_since_infection" in mapping:
        df_std["days_since_infection"] = pd.to_numeric(df_std[mapping["days_since_infection"]], errors="coerce")

    abs_rel_rows = []
    naacr_rel_rows = []
    summary_rows = []

    # Absolute metrics → relative (if present)
    abs_cols = [c for c in ["naa_mM", "cho_mM", "cr_mM"] if c in mapping]
    if abs_cols:
        df_abs = df_std.dropna(subset=[mapping[c] for c in abs_cols], how="all").copy()
        for c in abs_cols:
            col = mapping[c]
            df_abs[col] = pd.to_numeric(df_abs[col], errors="coerce")
        # control means per group
        for c in abs_cols:
            col = mapping[c]
            ctrl = df_abs[df_abs["is_control"] == True]
            if ctrl.empty:
                continue
            means = compute_control_means(ctrl, group_cols, col)
            df_abs = attach_rel(df_abs, col, means, group_cols, out_col=f"{c.replace('_mM','')}_rel")

        keep_cols = [
            "study_name", "study_id", "subject_id", "phase", "is_control",
            "days_since_infection"
        ] + (["region"] if has_region else []) + [
            c for c in ["naa_rel", "cho_rel", "cr_rel"] if c in df_abs.columns
        ]
        abs_rel_rows.append(df_abs[keep_cols])

        # Summary rows for abs
        for c in ["naa_rel", "cho_rel", "cr_rel"]:
            if c in df_abs.columns:
                grp = df_abs.groupby(group_cols + ["phase"]) [c]
                summ = grp.agg(["count", "mean", "std"]).reset_index()
                summ.insert(0, "metric", c)
                summary_rows.append(summ)

    # NAA/Cr (and Cho/Cr) → relative (if present)
    ratio_cols = []
    if "naa_over_cr" in mapping:
        ratio_cols.append((mapping["naa_over_cr"], "naa_over_cr_rel"))
    if "cho_over_cr" in mapping:
        ratio_cols.append((mapping["cho_over_cr"], "cho_over_cr_rel"))

    if ratio_cols:
        df_ratio = df_std.dropna(subset=[c[0] for c in ratio_cols], how="all").copy()
        for c_in, _ in ratio_cols:
            df_ratio[c_in] = pd.to_numeric(df_ratio[c_in], errors="coerce")
        ctrl = df_ratio[df_ratio["is_control"] == True]
        if not ctrl.empty:
            for c_in, c_out in ratio_cols:
                means = compute_control_means(ctrl, group_cols, c_in)
                df_ratio = attach_rel(df_ratio, c_in, means, group_cols, out_col=c_out)

            keep_cols = [
                "study_name", "study_id", "subject_id", "phase", "is_control",
                "days_since_infection"
            ] + (["region"] if has_region else []) + [c_out for _, c_out in ratio_cols]
            naacr_rel_rows.append(df_ratio[keep_cols])

            # Summary rows for ratios
            for _, c_out in ratio_cols:
                grp = df_ratio.groupby(group_cols + ["phase"]) [c_out]
                summ = grp.agg(["count", "mean", "std"]).reset_index()
                summ.insert(0, "metric", c_out)
                summary_rows.append(summ)

    abs_rel_df = pd.concat(abs_rel_rows, ignore_index=True) if abs_rel_rows else pd.DataFrame()
    naacr_rel_df = pd.concat(naacr_rel_rows, ignore_index=True) if naacr_rel_rows else pd.DataFrame()
    summary_df = pd.concat(summary_rows, ignore_index=True) if summary_rows else pd.DataFrame()

    return abs_rel_df, naacr_rel_df, summary_df


def main():
    ap = argparse.ArgumentParser(description="Prepare patient-level relative datasets (absolute→relative and NAA/Cr→relative)")
    ap.add_argument("--extracted-dir", required=True, help="Directory containing extracted patient-level CSV/XLSX files")
    ap.add_argument("--out-dir", required=True, help="Output directory for curated CSVs (e.g., quantum/data/curated)")
    ap.add_argument("--summary-dir", default=str(Path.cwd().parent.parent / "data" / "model_outputs"), help="Directory for summary CSV")
    args = ap.parse_args()

    extracted_dir = Path(args.extracted_dir).expanduser()
    out_dir = Path(args.out_dir).expanduser()
    summary_dir = Path(args.summary_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    all_abs_rel: List[pd.DataFrame] = []
    all_naacr_rel: List[pd.DataFrame] = []
    all_summary: List[pd.DataFrame] = []

    study_id_counter: Dict[str, int] = {}

    files = sorted(list(extracted_dir.glob("**/*")))
    if not files:
        print(f"[WARN] No files found under {extracted_dir}")

    for path in files:
        if path.suffix.lower() not in {".csv", ".xls", ".xlsx"}:
            continue
        abs_rel_df, naacr_rel_df, summary_df = process_file(path, study_id_counter)
        if not abs_rel_df.empty:
            all_abs_rel.append(abs_rel_df)
            print(f"[OK] ABS→REL from {path} | rows: {len(abs_rel_df)}")
        if not naacr_rel_df.empty:
            all_naacr_rel.append(naacr_rel_df)
            print(f"[OK] NAA/Cr→REL from {path} | rows: {len(naacr_rel_df)}")
        if not summary_df.empty:
            # add source filename for traceability
            summary_df.insert(0, "source", str(path))
            all_summary.append(summary_df)

    if all_abs_rel:
        abs_rel = pd.concat(all_abs_rel, ignore_index=True)
        abs_rel.to_csv(out_dir / "patient_abs_rel.csv", index=False)
        print(f"[SAVE] {out_dir / 'patient_abs_rel.csv'} | rows: {len(abs_rel)}")
    else:
        print("[INFO] No ABS→REL datasets produced")

    if all_naacr_rel:
        naacr_rel = pd.concat(all_naacr_rel, ignore_index=True)
        naacr_rel.to_csv(out_dir / "patient_naacr_rel.csv", index=False)
        print(f"[SAVE] {out_dir / 'patient_naacr_rel.csv'} | rows: {len(naacr_rel)}")
    else:
        print("[INFO] No NAA/Cr→REL datasets produced")

    if all_summary:
        summary = pd.concat(all_summary, ignore_index=True)
        # Aggregate to a compact summary per study/region/phase/metric
        group_cols = [c for c in ["study_name", "region", "phase", "metric"] if c in summary.columns]
        if group_cols:
            compact = summary.groupby(group_cols).agg(
                count=("count", "sum"),
                mean=("mean", "mean"),
                std=("std", "mean"),
            ).reset_index()
        else:
            compact = summary
        out_summary = Path("/data/model_outputs") / "patient_ratios_summary.csv"
        compact.to_csv(out_summary, index=False)
        print(f"[SAVE] {out_summary} | rows: {len(compact)}")

    print("\nDone. Review curated CSVs under:")
    print(f"  - {out_dir}")
    print("Summaries under:")
    print(f"  - {summary_dir}")


if __name__ == "__main__":
    main()
