#!/usr/bin/env python3
"""
Data Organization Script for HIV MRS Project
Organizes extracted datasets from Claude project into local Git repo structure

Usage:
    python organize_project_data.py

This will create a structured data directory and move files appropriately.
"""

import os
import shutil
from pathlib import Path

# Define your local Git repo path
# UPDATE THIS PATH TO YOUR ACTUAL REPO LOCATION
LOCAL_REPO = Path(LOCAL_REPO = Path("/Users/acdstudpro/Documents/Github/noise_decorrelation_HIV/data")
)

# Define data directory structure
DATA_DIRS = {
    "raw": LOCAL_REPO / "data" / "raw",  # Original tables/Excel
    "extracted": LOCAL_REPO / "data" / "extracted",  # Study-specific CSVs
    "master": LOCAL_REPO / "data" / "master",  # Integrated databases
    "model_outputs": LOCAL_REPO / "data" / "model_outputs",  # Bayesian results
    "individual": LOCAL_REPO / "data" / "individual",  # Patient-level data (Valcour)
}

# File categorization
FILES_TO_MOVE = {
    "raw": [
        "MRSHIVSuppMatDahmaniRev.xlsx",
        "Table_6.xls",
        "Table3.xls",
        "Table_1_PerinatalHIVInfectionorExposureIsAssociatedWithLowNAcetylaspartateandGlutamateinBasalGangliaatAge9butNot7Years.docx",
        "Table_2_PerinatalHIVInfectionorExposureIsAssociatedWithLowNAcetylaspartateandGlutamateinBasalGangliaatAge9butNot7Years.DOCX",
        "Table_1_Fixel-Based_Analysis_and_Free_Water_Corrected_DTI_Evaluation_of_HIV-Associated_Neurocognitive_Disorders.docx",
        "Table_3_Atypical_Resting-State_Functional_Connectivity_Dynamics_Correlate_With_Early_Cognitive_Dysfunction_in_HIV_Infection.docx",
    ],
    "extracted": [
        "CHANG_2002_EXTRACTED.csv",
        "SAILASUTA_2012_EXTRACTED.csv",
        "SAILASUTA_2012_ACUTE_DATA.csv",
        "SAILASUTA_2016_LONGITUDINAL.csv",
        "YOUNG_2014_CROSS_SECTIONAL_DATA.csv",
        "YOUNG_2014_LONGITUDINAL.csv",
        "CRITICAL_STUDIES_COMPLETE_DATA.csv",
    ],
    "master": [
        "MASTER_HIV_MRS_DATABASE_v2.csv",
        "NAA_DATA_FOR_MODEL.csv",  # This is what v3.6 used
    ],
    "model_outputs": [
        "posterior_predictive.csv",
        "summary_v4.csv",
    ],
    "individual": [
        "VALCOUR_2015_INDIVIDUAL_PATIENTS.csv",
        "VALCOUR_2015_ACUTE_DATA.csv",  # Might be aggregated from individual
        "VALCOUR_2015_DATA_FOR_MASTER.csv",
    ],
}


def create_directory_structure():
    """Create organized data directory structure"""
    print("Creating directory structure...")
    for dir_type, dir_path in DATA_DIRS.items():
        dir_path.mkdir(parents=True, exist_ok=True)
        print(f"  ✓ Created: {dir_path}")


def copy_files_from_project(project_dir="/mnt/project"):
    """Copy files from Claude project to organized structure"""
    project_path = Path(project_dir)

    print("\nCopying files to organized structure...")
    for category, files in FILES_TO_MOVE.items():
        dest_dir = DATA_DIRS[category]
        print(f"\n{category.upper()}:")

        for filename in files:
            src = project_path / filename
            dest = dest_dir / filename

            if src.exists():
                shutil.copy2(src, dest)
                print(f"  ✓ {filename}")
            else:
                print(f"  ✗ {filename} (not found)")


def create_readme():
    """Create README for data directory"""
    readme_content = """# HIV MRS Data Directory

## Structure

### raw/
Original tables and Excel files from published studies
- Contains supplementary materials, raw data tables

### extracted/
Study-specific CSV extractions with group-level statistics
- These are AGGREGATED data (means, SDs, SEs)
- Used for meta-analytic approaches
- Format: Study, Group (Control/Acute/Chronic), n, Metabolite, Mean, SE/SD

### master/
Integrated databases combining multiple studies
- `MASTER_HIV_MRS_DATABASE_v2.csv`: Comprehensive standardized format
- `NAA_DATA_FOR_MODEL.csv`: Input for Bayesian v3.6 model

### individual/
Individual patient-level data (RARE AND PRECIOUS!)
- `VALCOUR_2015_INDIVIDUAL_PATIENTS.csv`: 62 patients, multiple brain regions
- This enables hierarchical modeling with patient random effects
- Cannot be directly compared via WAIC to aggregated data models

### model_outputs/
Results from Bayesian inference
- Posterior predictive samples
- Parameter summaries

## Data Types and Modeling Implications

### Aggregated Data (extracted/)
- **Structure**: Group means ± SE/SD
- **Modeling**: Group-level parameters with measurement error
- **Limitations**: Cannot estimate patient-level heterogeneity
- **Use for**: v3.6-style model (3 group comparison)

### Individual Data (individual/)
- **Structure**: One row per patient
- **Modeling**: Full hierarchical Bayesian with random effects
- **Advantages**: Proper uncertainty quantification, patient heterogeneity
- **Use for**: v4+ models, validation studies

## Critical Note on Model Comparison

⚠️ **WARNING**: Cannot directly compare WAIC between:
- Models fit to aggregated data (3 group means)
- Models fit to individual data (62+ patients)

These have fundamentally different likelihoods and data structures!

## Pipeline Strategy

**Pipeline A: Aggregated Data Analysis**
- Input: extracted/ + master/
- Model: v3.6 approach (3 groups)
- Output: Group-level protection factors

**Pipeline B: Individual Data Analysis**
- Input: individual/
- Model: Hierarchical with patient random effects
- Output: Patient-level predictions, population parameters

**Validation**: Use individual data to validate predictions from aggregated analysis
"""

    readme_path = LOCAL_REPO / "data" / "README.md"
    with open(readme_path, 'w') as f:
        f.write(readme_content)
    print(f"\n✓ Created: {readme_path}")


def create_data_inventory():
    """Create detailed inventory of all datasets"""
    inventory = """# HIV MRS Data Inventory

## Study-by-Study Breakdown

### 1. Sailasuta 2012 (SEARCH 010/RV254)
- **Type**: Aggregated (means ± SD)
- **Groups**: Control (n=10), Acute HIV (n=31, 14d post-infection), Chronic (n=26)
- **Regions**: BG, OGM
- **Metabolites**: NAA, Cho, mI, Glx
- **Format**: Ratio to Creatine
- **Key finding**: NAA preserved in acute phase despite peak viremia

### 2. Valcour 2015 (RV254)
- **Type**: INDIVIDUAL patient data (n=62)
- **Time point**: Acute HIV (Fiebig I-V)
- **Regions**: FGM, FWM, BG, PGM
- **Metabolites**: NAA, Cho (absolute concentrations)
- **Also includes**: Age, CD4, log viral load, NPZ4 cognitive scores
- **Key finding**: +7.7% NAA elevation (p=0.0317) in BG

### 3. Young 2014
- **Type**: Aggregated
- **Groups**: Control (n=19), PHI ~6mo (n=9), Chronic ~10y (n=18)
- **Region**: BG
- **Format**: Ratio to Creatine
- **Key finding**: Intermediate NAA in PHI group

### 4. Chang 2002
- **Type**: Aggregated
- **Groups**: Control (n=25), Early HIV ~2y (n=45)
- **Regions**: BG, FWM, FGM
- **Metabolites**: NAA, Cho, mI
- **Format**: Absolute concentrations (mM)

### 5. Sailasuta 2016
- **Type**: Aggregated longitudinal
- **Groups**: ART-naive baseline, 6mo, 12mo
- **Regions**: BG, FWM
- **Intervention**: Treatment response to cART

### 6. Dahmani 2021
- **Type**: Meta-analysis
- **Acute**: 146 patients, 9 studies
- **Chronic**: 943 patients, 60 studies
- **Result**: No significant NAA change in acute (p=NS), significant decrease in chronic (p<0.001)

## Data Quality Assessment

### High Quality (usable for modeling)
✓ Sailasuta 2012 - Full statistics, multiple metabolites
✓ Valcour 2015 - Individual data (gold standard)
✓ Young 2014 - Clear group definitions
✓ Chang 2002 - Absolute concentrations

### Medium Quality (usable with caveats)
~ Sailasuta 2016 - Missing SD/SE, longitudinal confounds

### Meta-analytic Only
✓ Dahmani 2021 - Validates overall pattern, not for mechanistic modeling

## Current Model Inputs

### v3.6 Model (VALIDATED)
**Input**: `NAA_DATA_FOR_MODEL.csv`
- Control: 8.76 ± 0.79 mM (Chang 2002, BG)
- Acute: 10.0 (Sailasuta 2016, BG, ART-naive)
- Chronic: 7.96 ± 0.91 mM (Chang 2002, BG)

**Result**: P(ξ_acute < ξ_chronic) > 0.999

### Future v4 Model Options
**Option A**: Same aggregated data with enzyme kinetics
**Option B**: Valcour individual data with hierarchical structure
**Option C**: Combined approach (requires careful likelihood specification)

## Missing Data / Opportunities

### What We Need
1. More acute phase individual data (not just Valcour)
2. Longitudinal follow-up: acute → chronic same patients
3. Viral load + metabolite + cognitive scores
4. Regional specificity in acute phase

### Who Might Have It
- Dr. Serena Spudich (Yale) - RV254 extended analyses
- Dr. Sailasuta (?)
- Dr. Valcour - may have additional unpublished data

## Data Access Notes

- Most published data is aggregated (authors report means, not raw data)
- Individual data is rare because:
  * Privacy concerns
  * Not typically required for publication
  * Authors may not share
- Valcour 2015 individual data = GOLD MINE
- Consider reaching out for expanded datasets after first publication
"""

    inventory_path = LOCAL_REPO / "data" / "DATA_INVENTORY.md"
    with open(inventory_path, 'w') as f:
        f.write(inventory)
    print(f"✓ Created: {inventory_path}")


def main():
    """Main execution"""
    print("=" * 60)
    print("HIV MRS Project Data Organization")
    print("=" * 60)

    # Check if local repo exists
    if not LOCAL_REPO.exists():
        print(f"\n❌ ERROR: Local repo not found at {LOCAL_REPO}")
        print("Please update LOCAL_REPO path in the script!")
        return

    # Create structure
    create_directory_structure()

    # Copy files (only if running from within project)
    if Path("/mnt/project").exists():
        copy_files_from_project()
    else:
        print("\n⚠️  Not running in Claude project environment")
        print("   Files not copied - run this script from within your project")

    # Create documentation
    create_readme()
    create_data_inventory()

    print("\n" + "=" * 60)
    print("✓ Data organization complete!")
    print("=" * 60)
    print(f"\nNext steps:")
    print(f"1. cd {LOCAL_REPO}")
    print(f"2. Review data/README.md and data/DATA_INVENTORY.md")
    print(f"3. git add data/")
    print(f"4. git commit -m 'Organize extracted datasets into structured directories'")
    print(f"5. Consider the dual pipeline approach for modeling")


if __name__ == "__main__":
    main()