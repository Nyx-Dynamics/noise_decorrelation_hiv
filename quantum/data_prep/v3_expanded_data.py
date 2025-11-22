"""
V3 expanded data loader

- Reads mixed-format files from a source directory and normalizes columns
- Returns (data_dict, df_individual, df_group)
- CLI usage:
    python -m quantum.data_prep.v3_expanded_data --src /path/to/extracted_expanded --out data_dict.json --group group.csv --rows rows.csv

Notes:
- Phases are normalized with preference for the term "control" instead of "healthy".
- Ratios NAA/Cr and Cho/Cr are computed if absolute NAA/Cho/Cr are present.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd

from .common_expanded_loader import load_expanded_dataset


DEFAULT_SRC = "/Users/acdmbpmax/PycharmProjects/noise_decorrelation_HIV/data/extracted_expanded"


def load_v3_expanded_data(src_dir: str = DEFAULT_SRC) -> Tuple[Dict[str, dict], pd.DataFrame, pd.DataFrame]:
    data_dict, df_rows, df_group = load_expanded_dataset(src_dir)
    return data_dict, df_rows, df_group


def _save_optional(obj, path: str | None):
    if not path:
        return
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(obj, (dict, list)):
        p.write_text(json.dumps(obj, indent=2))
    elif isinstance(obj, pd.DataFrame):
        if p.suffix.lower() == ".json":
            p.write_text(obj.to_json(orient="records", indent=2))
        else:
            obj.to_csv(p, index=False)
    else:
        p.write_text(str(obj))


def main():
    ap = argparse.ArgumentParser(description="Load V3 expanded dataset")
    ap.add_argument("--src", default=DEFAULT_SRC, help="Source directory with expanded files")
    ap.add_argument("--out", default=None, help="Optional JSON path to save dict-of-dicts")
    ap.add_argument("--group", default=None, help="Optional CSV/JSON path to save group-level table")
    ap.add_argument("--rows", default=None, help="Optional CSV/JSON path to save individual rows")
    args = ap.parse_args()

    data_dict, df_rows, df_group = load_v3_expanded_data(args.src)

    # Print a short summary and emit JSON if requested
    print(f"Loaded: rows={len(df_rows)}, groups={len(df_group)} from {args.src}")
    # Show phases present
    if "phase" in df_rows.columns:
        print("Phases (rows):", sorted(df_rows["phase"].astype(str).unique()))
    if "study" in df_group.columns:
        print("Studies (groups):", sorted(df_group["study"].astype(str).unique()))

    _save_optional(data_dict, args.out)
    _save_optional(df_group, args.group)
    _save_optional(df_rows, args.rows)

    if not args.out and not args.group and not args.rows:
        # If no outputs requested, print dict-of-dicts to stdout
        print(json.dumps(data_dict, indent=2))


if __name__ == "__main__":
    main()
