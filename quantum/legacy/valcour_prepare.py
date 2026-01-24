#!/usr/bin/env python3
"""
Valcour Individual Patient Data Preparation (Pipeline 1)

Loads VALCOUR_2015_INDIVIDUAL_PATIENTS.csv (n=62 acute HIV patients)
Outputs: valcour_abs.parquet (absolute mM concentrations, NO conversion)

This is Pipeline 1 - your PRIMARY analysis with highest quality data
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Paths
DATA_DIR = Path(__file__).parent.parent / "data"
INDIVIDUAL_DIR = DATA_DIR / "individual"
PROCESSED_DIR = DATA_DIR / "processed"
PROCESSED_DIR.mkdir(exist_ok=True)


def load_valcour():
    """Load raw Valcour individual patient data"""
    filepath = INDIVIDUAL_DIR / "VALCOUR_2015_INDIVIDUAL_PATIENTS.csv"

    if not filepath.exists():
        raise FileNotFoundError(f"Cannot find: {filepath}")

    df = pd.read_csv(filepath)

    print(f"Loaded {len(df)} patients from Valcour 2015")
    print(f"Columns: {list(df.columns)}")

    return df


def clean_and_validate(df):
    """Clean data and validate quality"""

    print("\n" + "=" * 70)
    print("DATA CLEANING & VALIDATION")
    print("=" * 70)

    # Check for expected columns
    required_cols = ['Age', 'CD4', 'logpVL']
    metabolite_cols = ['FGMNAA', 'FGMCho', 'FWMNAA', 'FWMCho',
                       'BGNAA', 'BGCho', 'PGMNAA', 'PGMCho']

    missing_required = [c for c in required_cols if c not in df.columns]
    if missing_required:
        raise ValueError(f"Missing required columns: {missing_required}")

    # Count missing values per region
    print("\nMissing data by brain region:")
    for region in ['FGM', 'FWM', 'BG', 'PGM']:
        naa_col = f"{region}NAA"
        cho_col = f"{region}Cho"

        if naa_col in df.columns:
            n_missing_naa = df[naa_col].isna().sum()
            n_missing_cho = df[cho_col].isna().sum()
            n_complete = (~df[naa_col].isna() & ~df[cho_col].isna()).sum()

            print(f"  {region}: {n_complete} complete pairs, "
                  f"{n_missing_naa} missing NAA, {n_missing_cho} missing Cho")

    # Patients with complete BG data (primary region)
    df_bg_complete = df[~df['BGNAA'].isna() & ~df['BGCho'].isna()].copy()
    print(f"\nPatients with complete BG data: {len(df_bg_complete)}/{len(df)}")

    # Check for outliers (basic sanity checks)
    print("\nMetabolite value ranges (BG):")
    print(f"  NAA: {df_bg_complete['BGNAA'].min():.2f} - "
          f"{df_bg_complete['BGNAA'].max():.2f} mM")
    print(f"  Cho: {df_bg_complete['BGCho'].min():.2f} - "
          f"{df_bg_complete['BGCho'].max():.2f} mM")

    # Covariate distributions
    print("\nCovariate distributions:")
    print(f"  Age: {df_bg_complete['Age'].mean():.1f} ± "
          f"{df_bg_complete['Age'].std():.1f} years")
    print(f"  CD4: {df_bg_complete['CD4'].mean():.0f} ± "
          f"{df_bg_complete['CD4'].std():.0f} cells/μL")
    print(f"  log VL: {df_bg_complete['logpVL'].mean():.2f} ± "
          f"{df_bg_complete['logpVL'].std():.2f}")

    return df


def prepare_long_format(df):
    """
    Convert wide format to long format for modeling
    Creates one row per patient-region combination
    """

    # Create patient ID if not present
    if 'patient_id' not in df.columns:
        df['patient_id'] = range(len(df))

    # Regions and metabolites
    regions = ['FGM', 'FWM', 'BG', 'PGM']
    region_names = {
        'FGM': 'Frontal_GM',
        'FWM': 'Frontal_WM',
        'BG': 'Basal_Ganglia',
        'PGM': 'Parietal_GM'
    }

    # Build long format
    rows = []
    for _, patient in df.iterrows():
        for region_code, region_name in region_names.items():
            naa_col = f"{region_code}NAA"
            cho_col = f"{region_code}Cho"

            if naa_col in df.columns and pd.notna(patient[naa_col]):
                row = {
                    'patient_id': patient['patient_id'],
                    'age': patient['Age'],
                    'cd4': patient['CD4'],
                    'log_vl': patient['logpVL'],
                    'npz4': patient.get('NPZ4', np.nan),
                    'region': region_name,
                    'NAA_mM': patient[naa_col],
                    'Cho_mM': patient[cho_col] if cho_col in df.columns else np.nan,
                }
                rows.append(row)

    df_long = pd.DataFrame(rows)

    print(f"\nConverted to long format: {len(df_long)} observations "
          f"({len(df_long['patient_id'].unique())} patients × "
          f"{len(df_long['region'].unique())} regions)")

    return df_long


def save_outputs(df, df_long):
    """Save both wide and long format versions"""

    # Wide format (original structure)
    wide_path = PROCESSED_DIR / "valcour_abs_wide.parquet"
    df.to_parquet(wide_path, index=False)
    print(f"\n✓ Saved wide format: {wide_path}")

    # Long format (for hierarchical modeling)
    long_path = PROCESSED_DIR / "valcour_abs_long.parquet"
    df_long.to_parquet(long_path, index=False)
    print(f"✓ Saved long format: {long_path}")

    # Also save as CSV for easy inspection
    csv_path = PROCESSED_DIR / "valcour_abs_long.csv"
    df_long.to_csv(csv_path, index=False)
    print(f"✓ Saved CSV: {csv_path}")

    return wide_path, long_path


def main():
    """Main preparation pipeline"""

    print("=" * 70)
    print("VALCOUR INDIVIDUAL DATA PREPARATION (PIPELINE 1)")
    print("=" * 70)
    print("\nThis script prepares individual patient data (n=62)")
    print("Keeps absolute mM concentrations - NO conversion\n")

    # Load
    df = load_valcour()

    # Clean and validate
    df_clean = clean_and_validate(df)

    # Convert to long format for modeling
    df_long = prepare_long_format(df_clean)

    # Save outputs
    wide_path, long_path = save_outputs(df_clean, df_long)

    # Summary
    print("\n" + "=" * 70)
    print("PREPARATION COMPLETE")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  Wide format:  {wide_path.name}")
    print(f"  Long format:  {long_path.name} ← Use this for hierarchical model")

    print(f"\nData summary:")
    print(f"  Total patients: {len(df_clean)}")
    print(f"  Total observations: {len(df_long)}")
    print(f"  Regions: {sorted(df_long['region'].unique())}")
    print(f"  Metabolites: NAA, Cho (absolute mM)")
    print(f"  Covariates: Age, CD4, log viral load, NPZ4 cognition")

    print(f"\nReady for: models/hierarchical_valcour.py")
    print("=" * 70)


if __name__ == "__main__":
    main()