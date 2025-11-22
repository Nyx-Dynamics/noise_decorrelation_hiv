import json
import re
from pathlib import Path
from typing import Dict, Tuple, List, Optional, Any

import pandas as pd


# -----------------------------------------------------------------------------
# Helpers: IO
# -----------------------------------------------------------------------------

def _read_any(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in [".csv", ".tsv", ".txt"]:
        sep = "," if suffix == ".csv" else "\t"
        try:
            return pd.read_csv(path, sep=sep)
        except Exception:
            # Try flexible engine
            return pd.read_csv(path, sep=sep, engine="python")
    if suffix in [".xlsx", ".xls"]:
        return pd.read_excel(path)
    if suffix in [".json", ".ndjson"]:
        df = pd.read_json(path, lines=(suffix == ".ndjson"))
        # If it's a dict keyed by records, try to normalize
        if df.shape[1] == 1 and isinstance(df.iloc[0, 0], (list, dict)):
            return pd.json_normalize(df.iloc[0, 0])
        return df
    raise ValueError(f"Unsupported file type: {path}")


def read_tabular_files(src_dir: str, include: Optional[List[str]] = None,
                       exclude_hidden: bool = True) -> pd.DataFrame:
    """Read all tabular data files from a directory into a single DataFrame.

    - Supports CSV/TSV/TXT/XLSX/JSON.
    - `include`: optional list of glob patterns to filter files.
    """
    d = Path(src_dir)
    if not d.exists():
        raise FileNotFoundError(f"Source directory does not exist: {src_dir}")

    patterns = include or ["*.csv", "*.tsv", "*.txt", "*.xlsx", "*.xls", "*.json", "*.ndjson"]
    files: List[Path] = []
    for pat in patterns:
        files.extend(d.glob(pat))

    if exclude_hidden:
        files = [p for p in files if not any(part.startswith(".") for part in p.parts)]

    if not files:
        raise FileNotFoundError(f"No data files found in: {src_dir} with patterns {patterns}")

    frames = []
    for f in sorted(files):
        try:
            df = _read_any(f)
            if df is None or df.empty:
                continue
            df["__source_file__"] = str(f)
            frames.append(df)
        except Exception as e:
            # Keep going but note the failure
            frames.append(pd.DataFrame({"__source_file__": [str(f)], "__read_error__": [str(e)]}))
    if not frames:
        raise RuntimeError("No readable files produced a DataFrame.")
    # Concatenate by columns; align mismatched columns
    return pd.concat(frames, axis=0, ignore_index=True, sort=False)


# -----------------------------------------------------------------------------
# Normalization: columns and values
# -----------------------------------------------------------------------------

_CANONICAL_MAP = {
    # study identifiers
    "study": "study",
    "study_id": "study",
    "studyname": "study",
    "publication": "study",
    # publication/study timing
    "pub_year": "publication_year",
    "year": "publication_year",
    "publication_year": "publication_year",
    "publication_date": "publication_date",
    "study_year": "study_year",
    # subject / group
    "subject": "subject_id",
    "subject_id": "subject_id",
    "id": "subject_id",
    "participant": "subject_id",
    "n": "n",
    # phases / groups
    "phase": "phase",
    "group": "phase",
    "condition": "phase",
    "arm": "phase",
    # durations
    "duration": "duration_days",
    "days_since_infection": "duration_days",
    "time_since_infection_days": "duration_days",
    "duration_days": "duration_days",
    # region
    "region": "region",
    "roi": "region",
    # metabolites (absolute)
    "naa": "naa",
    "cho": "cho",
    "cr": "cr",
    # ratios
    "naa_cr": "naa_cr",
    "naa/cr": "naa_cr",
    "naa:cr": "naa_cr",
    "cho_cr": "cho_cr",
    "cho/cr": "cho_cr",
    "cho:cr": "cho_cr",
    # glucose/glutamate ratio if present
    "glu_cr": "glu_cr",
    "glu/cr": "glu_cr",
    "glu:cr": "glu_cr",
    # virologic and immunologic biomarkers
    "vl": "vl_plasma_copies",
    "viral_load": "vl_plasma_copies",
    "vl_plasma": "vl_plasma_copies",
    "plasma_vl": "vl_plasma_copies",
    "copies_vl": "vl_plasma_copies",
    "vl_copies": "vl_plasma_copies",
    "p_vl": "vl_plasma_copies",
    "log_vl": "vl_plasma_log10",
    "log10_vl": "vl_plasma_log10",
    "logvl": "vl_plasma_log10",
    "plasma_log_vl": "vl_plasma_log10",
    "csf_vl": "vl_csf_copies",
    "vl_csf": "vl_csf_copies",
    "csf_viral_load": "vl_csf_copies",
    "csf_log_vl": "vl_csf_log10",
    "log_vl_csf": "vl_csf_log10",
    "cd4": "cd4_abs",
    "total_cd4": "cd4_abs",
    "cd4_abs": "cd4_abs",
    "cd4_percent": "cd4_pct",
    "cd4_%": "cd4_pct",
    "cd4_pct": "cd4_pct",
    "cd4_cd8": "cd4_cd8_ratio",
    "cd4:cd8": "cd4_cd8_ratio",
    "cd4_cd8_ratio": "cd4_cd8_ratio",
    # summary stats
    "naa_mean": "naa_mean",
    "naa_sd": "naa_sd",
    "cho_mean": "cho_mean",
    "cho_sd": "cho_sd",
}

_PHASE_NORMALIZATION = {
    # prefer "control" over "healthy" per user guidance
    "healthy": "control",
    "control": "control",
    "ctrl": "control",
    "hc": "control",
    "acute": "acute",
    "hyperacute": "acute",
    "phi": "acute",  # primary HIV infection
    "early": "acute",
    "chronic": "chronic",
    "art": "art_controlled",
    "art-controlled": "art_controlled",
    "art_controlled": "art_controlled",
    "virally_suppressed": "art_controlled",
    "suppressed": "art_controlled",
}


def _slugify(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", s.strip().lower()).strip("_")


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    # Map columns to canonical names where possible
    mapping = {}
    for c in df.columns:
        c_slug = _slugify(str(c))
        mapping[c] = _CANONICAL_MAP.get(c_slug, c_slug)
    ndf = df.rename(columns=mapping).copy()

    # Phase normalization
    if "phase" in ndf.columns:
        ndf["phase"] = (
            ndf["phase"].astype(str).str.lower().map(lambda x: _PHASE_NORMALIZATION.get(x, x))
        )

    # Compute ratios if raw values exist
    if "naa_cr" not in ndf.columns and {"naa", "cr"}.issubset(ndf.columns):
        with pd.option_context("mode.use_inf_as_na", True):
            ndf["naa_cr"] = pd.to_numeric(ndf["naa"], errors="coerce") / pd.to_numeric(ndf["cr"], errors="coerce")
    if "cho_cr" not in ndf.columns and {"cho", "cr"}.issubset(ndf.columns):
        with pd.option_context("mode.use_inf_as_na", True):
            ndf["cho_cr"] = pd.to_numeric(ndf["cho"], errors="coerce") / pd.to_numeric(ndf["cr"], errors="coerce")

    # Parse/derive publication_year from publication_date if present
    if "publication_year" in ndf.columns:
        ndf["publication_year"] = pd.to_numeric(ndf["publication_year"], errors="coerce")
    if "publication_date" in ndf.columns and "publication_year" not in ndf.columns:
        dt = pd.to_datetime(ndf["publication_date"], errors="coerce")
        ndf["publication_year"] = dt.dt.year

    # Cast numeric summaries if present
    numeric_cast_cols = [
        "naa_mean", "naa_sd", "cho_mean", "cho_sd", "duration_days", "n",
        "vl_plasma_copies", "vl_plasma_log10", "vl_csf_copies", "vl_csf_log10",
        "cd4_abs", "cd4_pct", "cd4_cd8_ratio", "glu_cr",
    ]
    for col in numeric_cast_cols:
        if col in ndf.columns:
            ndf[col] = pd.to_numeric(ndf[col], errors="coerce")

    # Ensure study and phase
    if "study" not in ndf.columns:
        # If missing, derive from source file base name
        if "__source_file__" in ndf.columns:
            ndf["study"] = ndf["__source_file__"].map(lambda p: Path(str(p)).stem)
        else:
            ndf["study"] = "unknown_study"

    if "phase" not in ndf.columns:
        ndf["phase"] = "unknown"

    return ndf


# -----------------------------------------------------------------------------
# Aggregation and dict-of-dicts builder
# -----------------------------------------------------------------------------

def aggregate_group_stats(df: pd.DataFrame, by: Optional[List[str]] = None) -> pd.DataFrame:
    by = by or ["study", "phase"]

    # Choose metrics to aggregate
    metrics = {}
    for metric in ["naa_cr", "cho_cr"]:
        if metric in df.columns:
            metrics[metric] = ["mean", "std", "count"]

    # If precomputed summary is provided, prefer those
    # but still compute count "n" from rows if missing
    agg = None
    if metrics:
        agg = df.groupby(by).agg(metrics)
        # flatten columns
        agg.columns = [f"{m}_{stat}" for m, stat in agg.columns]
        agg = agg.reset_index()

    # Bring forward explicit summary stats if available at group level
    # by grouping and taking first non-null
    summary_cols = [c for c in ["naa_mean", "naa_sd", "cho_mean", "cho_sd", "duration_days", "n"] if c in df.columns]
    if summary_cols:
        firsts = df.groupby(by)[summary_cols].first().reset_index()
        agg = firsts if agg is None else pd.merge(agg, firsts, on=by, how="outer")

    if agg is None:
        # Minimal aggregation: just counts
        counts = df.groupby(by).size().reset_index(name="n")
        agg = counts

    # Ensure "n"
    if "n" not in agg.columns:
        # prefer explicit n from df if present
        if "n" in df.columns:
            n_df = df.groupby(by)["n"].sum().reset_index(name="n")
            agg = pd.merge(agg, n_df, on=by, how="left")
        else:
            # use counts from ratios if present else all rows
            src = "naa_cr_count" if "naa_cr_count" in agg.columns else ("cho_cr_count" if "cho_cr_count" in agg.columns else None)
            if src is not None:
                agg["n"] = agg[src]
            else:
                counts = df.groupby(by).size().reset_index(name="n")
                agg = pd.merge(agg, counts, on=by, how="left")

    # If means/sd missing but ratios available, map
    if "naa_mean" not in agg.columns and "naa_cr_mean" in agg.columns:
        agg["naa_mean"] = agg["naa_cr_mean"]
    if "naa_sd" not in agg.columns and "naa_cr_std" in agg.columns:
        agg["naa_sd"] = agg["naa_cr_std"]
    if "cho_mean" not in agg.columns and "cho_cr_mean" in agg.columns:
        agg["cho_mean"] = agg["cho_cr_mean"]
    if "cho_sd" not in agg.columns and "cho_cr_std" in agg.columns:
        agg["cho_sd"] = agg["cho_cr_std"]

    # Default duration to NaN if absent
    if "duration_days" not in agg.columns:
        agg["duration_days"] = pd.NA

    return agg


def build_dict_of_dicts(group_df: pd.DataFrame) -> Dict[str, Dict]:
    data: Dict[str, Dict] = {}
    # consistent study_id assignment per study label
    study_ids = {s: i for i, s in enumerate(sorted(group_df["study"].astype(str).unique()))}

    for idx, row in group_df.iterrows():
        name = f"{row.get('study', 'study')}_{row.get('phase', 'phase')}"
        data[name] = {
            "phase": str(row.get("phase", "")),
            "n": int(row.get("n", 0)) if pd.notna(row.get("n")) else 0,
            "duration_days": float(row.get("duration_days", float("nan"))) if pd.notna(row.get("duration_days")) else None,
            "NAA_mean": float(row.get("naa_mean", float("nan"))) if pd.notna(row.get("naa_mean")) else None,
            "NAA_sd": float(row.get("naa_sd", float("nan"))) if pd.notna(row.get("naa_sd")) else None,
            "Cho_mean": float(row.get("cho_mean", float("nan"))) if pd.notna(row.get("cho_mean")) else None,
            "Cho_sd": float(row.get("cho_sd", float("nan"))) if pd.notna(row.get("cho_sd")) else None,
            "study_id": int(study_ids.get(str(row.get("study", "")), 0)),
        }
    return data


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

def load_expanded_dataset(src_dir: str, aggregate_by: Optional[List[str]] = None) -> Tuple[Dict[str, Dict], pd.DataFrame, pd.DataFrame]:
    """Load expanded dataset from a directory and return:
    - data_dict: dict-of-dicts compatible with existing loaders
    - df_individual: normalized individual-level DataFrame
    - df_group: aggregated group-level DataFrame

    This function is robust to varying column names and file types.
    """
    raw = read_tabular_files(src_dir)
    norm = normalize_columns(raw)

    by = aggregate_by or ["study", "phase"]
    group = aggregate_group_stats(norm, by=by)

    data_dict = build_dict_of_dicts(group)

    return data_dict, norm, group


def apply_corrections(df: pd.DataFrame, corrections: Optional[Path]) -> pd.DataFrame:
    """Apply optional corrections from a JSON/CSV file.

    Expected schema examples:
    - JSON: {"rename_study": {"Sailasuta2012": "sailasuta_2012"}, "drop_rows": [{"study": "X", "phase": "Y"}]}
    - CSV: columns matching df for patching specific rows.
    """
    if not corrections or not corrections.exists():
        return df
    if corrections.suffix.lower() in [".json"]:
        spec = json.loads(corrections.read_text())
        # rename studies
        for old, new in spec.get("rename_study", {}).items():
            df.loc[df["study"].astype(str) == old, "study"] = new
        # drop rows
        for cond in spec.get("drop_rows", []):
            mask = pd.Series([True] * len(df))
            for k, v in cond.items():
                if k in df.columns:
                    mask &= df[k].astype(str) == str(v)
            df = df.loc[~mask].copy()
        # patch values
        for patch in spec.get("patch", []):
            sel = patch.get("where", {})
            vals = patch.get("set", {})
            if not sel or not vals:
                continue
            mask = pd.Series([True] * len(df))
            for k, v in sel.items():
                if k in df.columns:
                    mask &= df[k].astype(str) == str(v)
            for k, v in vals.items():
                if k in df.columns:
                    df.loc[mask, k] = v
        return df
    # CSV-based patches: merge by keys and prefer patch values where non-null
    patch_df = pd.read_csv(corrections)
    keys = [c for c in ["study", "phase", "region", "subject_id"] if c in df.columns and c in patch_df.columns]
    if not keys:
        return df
    merged = pd.merge(df, patch_df, on=keys, how="left", suffixes=("", "__patch"))
    for c in patch_df.columns:
        if c in keys:
            continue
        patched = f"{c}__patch"
        if patched in merged.columns:
            merged[c] = merged[c].where(merged[patched].isna(), merged[patched])
            merged.drop(columns=[patched], inplace=True)
    return merged
