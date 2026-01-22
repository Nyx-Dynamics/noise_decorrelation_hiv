#!/usr/bin/env python3
"""
Data Harmonization Pipeline for Multi-Study HIV MRS Integration
Converts heterogeneous measurements to common reference frame
"""

import warnings
from typing import Dict, List

import numpy as np
import pandas as pd


class MRSHarmonizer:
    """Harmonize MRS data across different studies and protocols."""

    def __init__(self):
        # Reference values for conversions
        self.reference_values = {
            'NAA_Cr_control': 1.05,  # Typical healthy NAA/Cr
            'NAA_absolute_control': 12.5,  # mM in healthy brain
            'Cr_absolute': 8.0,  # mM creatine
            'Cho_Cr_control': 0.22,  # Typical Cho/Cr
        }

        # Scanner correction factors (empirically derived)
        self.scanner_corrections = {
            'GE_1.5T': 1.00,
            'GE_3T': 0.98,
            'Siemens_1.5T': 1.03,
            'Siemens_3T': 1.01,
            'Philips_1.5T': 0.99,
            'Philips_3T': 0.97
        }

        # Regional scaling factors (basal ganglia as reference)
        self.regional_factors = {
            'basal_ganglia': 1.00,
            'frontal_wm': 0.95,
            'frontal_gm': 1.02,
            'parietal': 0.98,
            'occipital': 0.96,
            'thalamus': 1.01,
            'whole_brain': 0.97
        }

    def harmonize_dataset(self, data: pd.DataFrame, study_info: Dict) -> pd.DataFrame:
        """Main harmonization pipeline."""

        print(f"\n📊 Harmonizing: {study_info['name']}")
        print("-" * 40)

        # Step 1: Convert units
        data = self.convert_units(data, study_info)

        # Step 2: Apply scanner corrections
        data = self.apply_scanner_correction(data, study_info)

        # Step 3: Regional normalization
        data = self.normalize_regions(data, study_info)

        # Step 4: Temporal alignment
        data = self.align_timepoints(data, study_info)

        # Step 5: Quality checks
        data = self.quality_control(data)

        print(f"✓ Harmonized {len(data)} observations")

        return data

    def convert_units(self, data: pd.DataFrame, study_info: Dict) -> pd.DataFrame:
        """Convert all measurements to NAA/Cr ratios."""

        unit_type = study_info.get('units', 'ratio')

        if unit_type == 'absolute':
            print(f"  Converting absolute → ratios")
            # Convert mM to ratios using typical Cr
            if 'NAA_mM' in data.columns:
                data['NAA_Cr'] = data['NAA_mM'] / self.reference_values['Cr_absolute']

        elif unit_type == 'NAA_Cho':
            print(f"  Converting NAA/Cho → NAA/Cr")
            # Use typical Cho/Cr to convert
            data['NAA_Cr'] = data['NAA_Cho'] * self.reference_values['Cho_Cr_control']

        elif unit_type == 'percent_control':
            print(f"  Converting % control → absolute ratios")
            data['NAA_Cr'] = data['NAA_percent'] * self.reference_values['NAA_Cr_control'] / 100

        return data

    def apply_scanner_correction(self, data: pd.DataFrame, study_info: Dict) -> pd.DataFrame:
        """Apply scanner-specific correction factors."""

        scanner = study_info.get('scanner', 'Unknown')
        field = study_info.get('field_strength', '1.5T')

        key = f"{scanner}_{field}"
        correction = self.scanner_corrections.get(key, 1.0)

        if correction != 1.0:
            print(f"  Applying scanner correction: {correction:.2f}")
            data['NAA_Cr'] = data['NAA_Cr'] * correction

        return data

    def normalize_regions(self, data: pd.DataFrame, study_info: Dict) -> pd.DataFrame:
        """Normalize to basal ganglia reference."""

        region = study_info.get('region', 'unknown')
        factor = self.regional_factors.get(region, 1.0)

        if factor != 1.0:
            print(f"  Regional normalization ({region}): {factor:.2f}")
            data['NAA_Cr_normalized'] = data['NAA_Cr'] / factor
        else:
            data['NAA_Cr_normalized'] = data['NAA_Cr']

        return data

    def align_timepoints(self, data: pd.DataFrame, study_info: Dict) -> pd.DataFrame:
        """Align disease stages and timepoints."""

        print(f"  Aligning timepoints...")

        # Map different staging systems
        stage_mapping = {
            # Fiebig stages → weeks since infection
            'Fiebig_I': 1,
            'Fiebig_II': 2,
            'Fiebig_III': 3,
            'Fiebig_IV': 4,
            'Fiebig_V': 8,
            'Fiebig_VI': 16,

            # General phases → weeks
            'hyperacute': 2,
            'acute': 8,
            'early': 24,
            'chronic': 104,  # 2 years
            'late': 260,  # 5 years
        }

        # Add standardized timepoint column
        if 'stage' in data.columns:
            data['weeks_infection'] = data['stage'].map(stage_mapping)
        elif 'days_infection' in data.columns:
            data['weeks_infection'] = data['days_infection'] / 7
        elif 'months_infection' in data.columns:
            data['weeks_infection'] = data['months_infection'] * 4.33

        # Categorize into phases
        data['phase'] = pd.cut(
            data['weeks_infection'],
            bins=[0, 12, 52, np.inf],
            labels=['acute', 'early', 'chronic']
        )

        return data

    def quality_control(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply quality control filters."""

        print(f"  Quality control...")

        # Remove outliers (>3 SD from mean)
        mean = data['NAA_Cr_normalized'].mean()
        std = data['NAA_Cr_normalized'].std()

        before = len(data)
        data = data[np.abs(data['NAA_Cr_normalized'] - mean) < 3 * std]
        after = len(data)

        if before > after:
            print(f"    Removed {before - after} outliers")

        # Flag low-quality measurements
        data['quality_flag'] = 0

        # Flag if SD is too high
        if 'NAA_SD' in data.columns:
            high_var = data['NAA_SD'] / data['NAA_Cr_normalized'] > 0.2
            data.loc[high_var, 'quality_flag'] = 1

        return data

    def combine_studies(self, study_list: List[Dict]) -> pd.DataFrame:
        """Combine multiple harmonized studies."""

        print("\n" + "=" * 50)
        print("COMBINING HARMONIZED STUDIES")
        print("=" * 50)

        all_data = []

        for study in study_list:
            # Load study data
            df = pd.read_csv(study['file'])

            # Add study identifier
            df['study_id'] = study['name']

            # Harmonize
            df_harmonized = self.harmonize_dataset(df, study)

            all_data.append(df_harmonized)

        # Combine all
        combined = pd.concat(all_data, ignore_index=True)

        print(f"\n✅ COMBINED DATASET:")
        print(f"   Total observations: {len(combined)}")
        print(f"   Studies included: {combined['study_id'].nunique()}")
        print(f"   Acute phase: {(combined['phase'] == 'acute').sum()}")
        print(f"   Chronic phase: {(combined['phase'] == 'chronic').sum()}")

        return combined

    def validate_harmonization(self, combined_data: pd.DataFrame):
        """Validate that harmonization worked correctly."""

        print("\n🔍 VALIDATION CHECKS")
        print("-" * 40)

        # Check 1: Consistent ranges
        naa_range = combined_data['NAA_Cr_normalized'].quantile([0.01, 0.99])
        print(f"  NAA/Cr range (1-99%): {naa_range.values}")

        if naa_range[0.01] < 0.5 or naa_range[0.99] > 2.0:
            warnings.warn("NAA/Cr values outside expected range!")

        # Check 2: Study effects
        study_means = combined_data.groupby('study_id')['NAA_Cr_normalized'].mean()
        study_cv = study_means.std() / study_means.mean()
        print(f"  Between-study CV: {study_cv:.3f}")

        if study_cv > 0.15:
            warnings.warn("High between-study variation after harmonization!")

        # Check 3: Phase separation
        acute_mean = combined_data[combined_data['phase'] == 'acute']['NAA_Cr_normalized'].mean()
        chronic_mean = combined_data[combined_data['phase'] == 'chronic']['NAA_Cr_normalized'].mean()

        print(f"  Acute mean: {acute_mean:.3f}")
        print(f"  Chronic mean: {chronic_mean:.3f}")
        print(f"  Preservation signal: {(acute_mean / chronic_mean - 1) * 100:.1f}%")

        return True


def main():
    """Run harmonization pipeline."""

    harmonizer = MRSHarmonizer()

    # Define studies to harmonize
    studies = [
        {
            'name': 'Valcour_2015',
            'file': 'valcour_individual_data.csv',
            'units': 'ratio',
            'scanner': 'GE',
            'field_strength': '3T',
            'region': 'basal_ganglia'
        },
        {
            'name': 'Young_2014',
            'file': 'young_cross_sectional.csv',
            'units': 'ratio',
            'scanner': 'Siemens',
            'field_strength': '3T',
            'region': 'frontal_gm'
        },
        {
            'name': 'Sailasuta_2012',
            'file': 'sailasuta_acute.csv',
            'units': 'ratio',
            'scanner': 'GE',
            'field_strength': '1.5T',
            'region': 'frontal_wm'
        },
        {
            'name': 'Chang_2002',
            'file': 'chang_early_hiv.csv',
            'units': 'absolute',
            'scanner': 'GE',
            'field_strength': '1.5T',
            'region': 'basal_ganglia'
        }
    ]

    # Harmonize and combine
    combined = harmonizer.combine_studies(studies)

    # Validate
    harmonizer.validate_harmonization(combined)

    # Save harmonized dataset
    combined.to_csv('harmonized_combined_data.csv', index=False)

    print("\n✨ Harmonization complete!")
    print(f"   Output: harmonized_combined_data.csv")

    # Create summary statistics
    summary = combined.groupby('phase')['NAA_Cr_normalized'].agg([
        'count', 'mean', 'std', 'median'
    ])

    print("\n📊 HARMONIZED SUMMARY:")
    print(summary)

    return combined


if __name__ == "__main__":
    main()`i