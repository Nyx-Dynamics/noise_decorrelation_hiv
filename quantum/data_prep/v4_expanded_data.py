"""
V4 expanded data loader

Differences vs V3:
- Allows aggregation by custom keys (e.g., study × phase × region).
- Exposes CLI flag --by to specify grouping (comma-separated).

CLI usage:
  python -m quantum.data_prep.v4_expanded_data \
      --src /path/to/extracted_expanded \
      --by study,phase,region \
      --out data_dict.json --group group.csv --rows rows.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple, List

import pandas as pd

from .common_expanded_loader import read_tabular_files, normalize_columns, aggregate_group_stats, build_dict_of_dicts


DEFAULT_SRC = "/Users/acdmbpmax/PycharmProjects/noise_decorrelation_HIV/data/extracted_expanded"


def load_v4_expanded_data(src_dir: str = DEFAULT_SRC, aggregate_by: List[str] | None = None) -> Tuple[Dict[str, dict], pd.DataFrame, pd.DataFrame]:
    raw = read_tabular_files(src_dir)
    rows = normalize_columns(raw)
    group = aggregate_group_stats(rows, by=aggregate_by or ["study", "phase"])
    data_dict = build_dict_of_dicts(group)
    return data_dict, rows, group


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
    ap = argparse.ArgumentParser(description="Load V4 expanded dataset with flexible grouping")
    ap.add_argument("--src", default=DEFAULT_SRC, help="Source directory with expanded files")
    ap.add_argument("--by", default="study,phase", help="Comma-separated grouping keys, e.g., 'study,phase,region'")
    ap.add_argument("--out", default=None, help="Optional JSON path to save dict-of-dicts")
    ap.add_argument("--group", default=None, help="Optional CSV/JSON path to save group-level table")
    ap.add_argument("--rows", default=None, help="Optional CSV/JSON path to save individual rows")
    args = ap.parse_args()

    by = [s.strip() for s in args.by.split(",") if s.strip()]
    data_dict, df_rows, df_group = load_v4_expanded_data(args.src, aggregate_by=by)

    print(f"Loaded: rows={len(df_rows)}, groups={len(df_group)} from {args.src}; by={by}")
    if "phase" in df_rows.columns:
        print("Phases (rows):", sorted(df_rows["phase"].astype(str).unique()))
    if "study" in df_group.columns:
        print("Studies (groups):", sorted(df_group["study"].astype(str).unique()))
    if "region" in df_group.columns:
        print("Regions (groups):", sorted(df_group["region"].astype(str).unique()))

    _save_optional(data_dict, args.out)
    _save_optional(df_group, args.group)
    _save_optional(df_rows, args.rows)

    if not args.out and not args.group and not args.rows:
        print(json.dumps(data_dict, indent=2))


if __name__ == "__main__":
    main()
