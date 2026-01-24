"""
Data loaders for ratio-comparison datasets with ART-era classification and
mechanistic Option A (xi → Pi_xi) derivation.

- Supports enzyme_inputs_<ratio>.csv and bayesian_inputs_<ratio>.csv
- Attaches art_era using cutoff year 2006 (<= pre_modern, >=2007 post_modern)
- Computes Pi_xi from xi_estimate_nm with beta_xi=1.89 by default
- Retains enzyme_activity_fold for validation only and computes validation_delta
- Adds derived covariates (log_VL, CD4/CD8 ratio, ART timing deltas) when present

This module is read-only and does not alter input CSVs.
"""
from __future__ import annotations

from pathlib import Path
from typing import Literal, Tuple

import numpy as np
import pandas as pd

# Root for the ratio comparison CSVs
RATIO_DIR = Path('data/extracted_expanded/data_ratios_comparison')

Era = Literal['pre_modern', 'post_modern', 'unknown']

# Fallback mapping for study → publication year (extend as needed)
PAPER_YEAR: dict[str, int] = {
    'Valcour_2015': 2015,
    'Young_2014': 2014,
    'Sailasuta_2012': 2012,
    'Chang_2002': 2002,
}

ART_CUTOFF_YEAR = 2006  # <= pre_modern, >=2007 post_modern


def classify_era(year: float | int | None, cutoff: int = ART_CUTOFF_YEAR) -> Era:
    if pd.isna(year):
        return 'unknown'
    try:
        y = int(year)
    except Exception:
        return 'unknown'
    return 'pre_modern' if y <= cutoff else 'post_modern'


def attach_era(df: pd.DataFrame, cutoff: int = ART_CUTOFF_YEAR) -> pd.DataFrame:
    df = df.copy()
    # Preferred fields: measurement_year / scan_year / publication_year
    year_series = None
    for col in ['measurement_year', 'scan_year', 'publication_year']:
        if col in df.columns:
            year_series = df[col]
            break

    if year_series is None:
        # Map per study if no explicit year column
        if 'study' in df.columns:
            df['publication_year'] = df['study'].map(PAPER_YEAR)
            year_series = df['publication_year']
        else:
            df['publication_year'] = np.nan
            year_series = df['publication_year']

    df['art_era'] = [classify_era(y, cutoff) for y in year_series]
    df['art_era_idx'] = df['art_era'].map({'pre_modern': 0, 'post_modern': 1, 'unknown': -1})
    return df


def add_derived_covariates(df: pd.DataFrame) -> pd.DataFrame:
    """Compute derived covariates when inputs exist. Preserve rows when missing."""
    df = df.copy()
    # pVL → log10 with safe lower bound
    if 'pVL' in df.columns:
        with np.errstate(divide='ignore'):
            df['log_VL'] = np.log10(pd.to_numeric(df['pVL'], errors='coerce').clip(lower=1.0))
    # CD4/CD8 ratio
    if {'CD4', 'CD8'}.issubset(df.columns):
        with np.errstate(divide='ignore', invalid='ignore'):
            c4 = pd.to_numeric(df['CD4'], errors='coerce')
            c8 = pd.to_numeric(df['CD8'], errors='coerce')
            df['CD4_CD8_ratio'] = c4 / c8
    # Time deltas (days) when date columns exist
    date_cols = ['hiv_diagnosis_date', 'art_start_date', 'measurement_date']
    has_dates = [c for c in date_cols if c in df.columns]
    if has_dates:
        for c in has_dates:
            df[c] = pd.to_datetime(df[c], errors='coerce')
        if {'hiv_diagnosis_date', 'art_start_date'}.issubset(df.columns):
            df['time_to_ART_initiation_days'] = (df['art_start_date'] - df['hiv_diagnosis_date']).dt.days
        if {'art_start_date', 'measurement_date'}.issubset(df.columns):
            df['time_since_ART_initiation_days'] = (df['measurement_date'] - df['art_start_date']).dt.days
    return df


def compute_Pi_xi_from_xi_nm(xi_nm: pd.Series, beta_xi: float = 1.89, xi_baseline_nm: float = 0.66) -> pd.Series:
    xi = pd.to_numeric(xi_nm, errors='coerce').astype(float).clip(lower=1e-9)
    return (xi / float(xi_baseline_nm)) ** (-beta_xi)


def load_enzyme_inputs(ratio: str, beta_xi: float = 1.89, xi_baseline_nm: float = 0.66) -> Tuple[pd.DataFrame, Path]:
    """
    Load enzyme inputs for a given ratio (e.g., '3_1_1', '1_2_1').

    Adds columns:
      - art_era, art_era_idx
      - Pi_xi (from xi_estimate_nm via Option A)
      - computed_Pi_xi (alias of Pi_xi)
      - validation_delta vs enzyme_activity_fold (if present)
      - derived covariates (log_VL, CD4_CD8_ratio, timing deltas) when source columns exist
    """
    path = RATIO_DIR / f'enzyme_inputs_{ratio}.csv'
    if not path.exists():
        raise FileNotFoundError(f"Missing enzyme inputs for ratio '{ratio}': {path}")
    df = pd.read_csv(path)
    df = attach_era(df)
    if 'xi_estimate_nm' in df.columns:
        df['Pi_xi'] = compute_Pi_xi_from_xi_nm(df['xi_estimate_nm'], beta_xi=beta_xi, xi_baseline_nm=xi_baseline_nm)
        df['computed_Pi_xi'] = df['Pi_xi']
    else:
        raise ValueError("enzyme_inputs CSV must include 'xi_estimate_nm' column")
    if 'enzyme_activity_fold' in df.columns and 'computed_Pi_xi' in df.columns:
        df['validation_delta'] = (pd.to_numeric(df['computed_Pi_xi'], errors='coerce') -
                                  pd.to_numeric(df['enzyme_activity_fold'], errors='coerce')).abs()
    df = add_derived_covariates(df)
    return df, path


def load_bayesian_inputs(ratio: str) -> Tuple[pd.DataFrame, Path]:
    path = RATIO_DIR / f'bayesian_inputs_{ratio}.csv'
    if not path.exists():
        raise FileNotFoundError(f"Missing bayesian inputs for ratio '{ratio}': {path}")
    df = pd.read_csv(path)
    df = attach_era(df)
    df = add_derived_covariates(df)
    return df, path
