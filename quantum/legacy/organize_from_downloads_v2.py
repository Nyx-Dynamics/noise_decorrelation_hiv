#!/usr/bin/env python3
"""
Organize HIV MRS data files from Downloads folder
Moves files from ~/Downloads to organized data/ structure
"""

import shutil
from pathlib import Path

# Source and destination
DOWNLOADS = Path.home() / "Downloads"
DATA_DIR = Path.home() / "Documents" / "Github" / "noise_decorrelation_HIV" / "data"

# File categorization - exactly as they appear in Downloads
FILES_TO_ORGANIZE = {
    "raw": [
        "MRS-HIV-SuppMat-Dahmani-Rev.xlsx",  # Note: has hyphens, not underscores
        "Table_6.xls",
        "Table3.xls",
        "Table3-2.xls",  # Appears to be a duplicate
        "valcour_2015.xlsx",  # Additional file
        "Table_1_PerinatalHIVInfectionorExposureIsAssociatedWithLowNAcetylaspartateandGlutamateinBasalGangliaatAge9butNot7Years.docx",
        "Table_2_PerinatalHIVInfectionorExposureIsAssociatedWithLowNAcetylaspartateandGlutamateinBasalGangliaatAge9butNot7Years.DOCX",
        "Table_1_Fixel-Based Analysis and Free Water Corrected DTI Evaluation of HIV-Associated Neurocognitive Disorders.docx",
        "Table_3_Atypical Resting-State Functional Connectivity Dynamics Correlate With Early Cognitive Dysfunction in HIV Infection.docx",
        "Data_Sheet_1_Atypical Resting-State Functional Connectivity Dynamics Correlate With Early Cognitive Dysfunction in HIV Infection.docx",
        "Data_Sheet_2_Atypical Resting-State Functional Connectivity Dynamics Correlate With Early Cognitive Dysfunction in HIV Infection.docx",
        "Data_Sheet_3_Atypical Resting-State Functional Connectivity Dynamics Correlate With Early Cognitive Dysfunction in HIV Infection.docx",
        "Table_1_Brain Volumetric Alterations in Preclinical HIV-Associated Neurocognitive Disorder Using Automatic Brain Quantification and Segmentation Tool.DOCX",
        # URL-encoded filename
        "Data%20Sheet%201_Neuroinflammation%20associated%20with%20proviral%20DNA%20persists%20in%20the%20brain%20of%20virally%20suppressed%20people%20with%20HIV.docx",
        "Table_1_Altered%20regional%20homogeneity%20and%20functional%20connectivity%20of%20brain%20activity%20in%20young%20HIV-infected%20patients%20with%20asymptomatic%20neurocognitive%20imp.DOCX",
    ],
    "extracted": [
        "SAILASUTA_2012_ACUTE_DATA.csv",
        "YOUNG_2014_CROSS_SECTIONAL_DATA.csv",
        "VALCOUR_2015_REGIONAL_SUMMARY.csv",  # Additional file
    ],
    "master": [
        "NAA_DATA_FOR_MODEL.csv",  # KEY FILE for v3.6!
    ],
    "model_outputs": [
        # These may not be in Downloads - that's OK
    ],
    "individual": [
        "VALCOUR_2015_INDIVIDUAL_PATIENTS.csv",  # THE GOLD MINE!
        "VALCOUR_2015_ACUTE_DATA.csv",
        "VALCOUR_2015_DATA_FOR_MASTER.csv",
    ],
    "documentation": [
        # PDF papers and analysis docs
        "valcour_2015.pdf",
        "sailasuta_2012.pdf",
        "sailasuta_2016.pdf",
        "Young_2014.pdf",
        "Change_2002.pdf",
        "Dahmani_2021.pdf",
        "Noise_HIV_NC_v1.pdf",
        "MRS-HIV-SuppMat-Dahmani-Rev.pdf",
        # Analysis documents
        "VALCOUR_2015_COMPLETE_EXTRACTION.md",
        "INTEGRATED_MRS_DATA_EXTRACTION.md",
        "MANUSCRIPT_DATA_EXTRACTION_MASTER.md",
        "ACUTE_HIV_MRS_PAPERS_ANALYSIS.md",
        "INTEGRATED_MASTER_EXTRACTION_V3.md",
        "VALCOUR_S1_DATASET_ANALYSIS.md",
        "YOUNG_2014_FIGURES_DETAILED_ANALYSIS.md",
        "MOHAMED_2010_CHRONIC_HIV_ANALYSIS.md",
        "WAIC_ANALYSIS_SUMMARY.md",
        # Other docs
        "sailasuta_2016_supplement.doc",
    ],
    "figures": [
        "valcour_2015_comprehensive_analysis.png",
        "DATA_EXTRACTION_SUMMARY.png",
        "1763181285273_plot6_comprehensive_summary.png",
        "sailasuta_2012_metabolite_figure.pdf",
    ],
}


def main():
    print("=" * 60)
    print("Organizing HIV MRS Data from Downloads")
    print("=" * 60)

    # Check if Downloads exists
    if not DOWNLOADS.exists():
        print(f"\n❌ Downloads folder not found: {DOWNLOADS}")
        return

    # Create all subdirectories first
    print(f"\nCreating directory structure...")
    for category in FILES_TO_ORGANIZE.keys():
        dest_dir = DATA_DIR / category
        dest_dir.mkdir(parents=True, exist_ok=True)
        print(f"  ✓ {dest_dir.name}/")

    print(f"\nSource: {DOWNLOADS}")
    print(f"Destination: {DATA_DIR}\n")

    moved_count = 0
    missing_count = 0

    for category, files in FILES_TO_ORGANIZE.items():
        if not files:  # Skip empty categories
            continue

        dest_dir = DATA_DIR / category
        print(f"\n{category.upper()}:")

        for filename in files:
            src = DOWNLOADS / filename
            dest = dest_dir / filename

            if src.exists():
                try:
                    # Copy (not move) to be safe - you can delete from Downloads later
                    shutil.copy2(src, dest)
                    print(f"  ✓ {filename}")
                    moved_count += 1
                except Exception as e:
                    print(f"  ✗ {filename} (error copying: {e})")
                    missing_count += 1
            else:
                print(f"  ✗ {filename} (not found)")
                missing_count += 1

    print("\n" + "=" * 60)
    print(f"✓ Organization complete!")
    print(f"  Files copied: {moved_count}")
    print(f"  Files not found: {missing_count}")
    print("=" * 60)

    if missing_count > 0:
        print(f"\n⚠️  {missing_count} files not found in Downloads")
        print("This is OK - you have the essential files!")

    # Note about what they have
    print(f"\n📋 Essential files you have:")
    print(f"  ✓ NAA_DATA_FOR_MODEL.csv (what v3.6 uses)")
    print(f"  ✓ VALCOUR_2015_INDIVIDUAL_PATIENTS.csv (the gold mine!)")
    print(f"  ✓ Study CSVs (Sailasuta, Young, Valcour)")

    print(f"\n📥 Optional files still in Claude project (if needed later):")
    print(f"  • MASTER_HIV_MRS_DATABASE_v2.csv (comprehensive database)")
    print(f"  • CRITICAL_STUDIES_COMPLETE_DATA.csv")
    print(f"  • CHANG_2002_EXTRACTED.csv")
    print(f"  • SAILASUTA_2012_EXTRACTED.csv")
    print(f"  • SAILASUTA_2016_LONGITUDINAL.csv")
    print(f"  • posterior_predictive.csv (model outputs)")
    print(f"  • summary_v4.csv")

    print(f"\nFiles were COPIED (not moved) from Downloads")
    print(f"Original files still in: {DOWNLOADS}")

    print(f"\nNext steps:")
    print(f"1. cd {DATA_DIR.parent}")
    print(f"2. git status  # See what was added")
    print(f"3. git add data/")
    print(f"4. git commit -m 'Add organized extracted datasets and documentation'")
    print(f"5. Review data/README.md for pipeline strategy")


if __name__ == "__main__":
    main()