#!/usr/bin/env python3
"""
Organize HIV MRS Project Files from Claude Project
Based on actual files in /mnt/project directory
"""

import shutil
from pathlib import Path

# Source: Claude project files (user should download these)
# Destination: Local Git repo
REPO_ROOT = Path.home() / "Documents" / "Github" / "noise_decorrelation_HIV"
DATA_DIR = REPO_ROOT / "data"

# Actual files from /mnt/project (verified to exist)
FILES_TO_ORGANIZE = {
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
        "NAA_DATA_FOR_MODEL.csv",
    ],

    "individual": [
        "VALCOUR_2015_INDIVIDUAL_PATIENTS.csv",
        "VALCOUR_2015_ACUTE_DATA.csv",
        "VALCOUR_2015_DATA_FOR_MASTER.csv",
    ],

    "model_outputs": [
        "posterior_predictive.csv",
        "summary_v4.csv",
    ],

    "raw_data": [
        "MRSHIVSuppMatDahmaniRev.xlsx",
        "Table3.xls",
        "Table_6.xls",
        "Table_1_PerinatalHIVInfectionorExposureIsAssociatedWithLowNAcetylaspartateandGlutamateinBasalGangliaatAge9butNot7Years.docx",
        "Table_2_PerinatalHIVInfectionorExposureIsAssociatedWithLowNAcetylaspartateandGlutamateinBasalGangliaatAge9butNot7Years.DOCX",
        "Table_1_Fixel-Based_Analysis_and_Free_Water_Corrected_DTI_Evaluation_of_HIV-Associated_Neurocognitive_Disorders.docx",
        "Table_3_Atypical_Resting-State_Functional_Connectivity_Dynamics_Correlate_With_Early_Cognitive_Dysfunction_in_HIV_Infection.docx",
        "Table_1_Brain_Volumetric_Alterations_in_Preclinical_HIV-Associated_Neurocognitive_Disorder_Using_Automatic_Brain_Quantification_and_Segmentation_Tool.DOCX",
        "Table_1_Altered_20regional_20homogeneity_20and_20functional_20connectivity_20of_20brain_20activity_20in_20young_20HIV-infected_20patients_20with_20asymptomatic_20neurocognitive_20imp.DOCX",
        "Data_20Sheet_201_Neuroinflammation_20associated_20with_20proviral_20DNA_20persists_20in_20the_20brain_20of_20virally_20suppressed_20people_20with_20HIV.docx",
        "Data_Sheet_1_Atypical_Resting-State_Functional_Connectivity_Dynamics_Correlate_With_Early_Cognitive_Dysfunction_in_HIV_Infection.docx",
        "Data_Sheet_2_Atypical_Resting-State_Functional_Connectivity_Dynamics_Correlate_With_Early_Cognitive_Dysfunction_in_HIV_Infection.docx",
        "Data_Sheet_3_Atypical_Resting-State_Functional_Connectivity_Dynamics_Correlate_With_Early_Cognitive_Dysfunction_in_HIV_Infection.docx",
        "MRS-HIV-SuppMat-Dahmani-Rev.docx",
        "CHATGPT_5_0_pro_recs.docx",
    ],

    "documentation": [
        # Markdown analysis files
        "COMPLETE_LITERATURE_SYNTHESIS.md",
        "COMPLETE_EVIDENCE_BASE_SUMMARY.md",
        "NEW_AUTHOR_EXTRACTIONS_SUMMARY.md",
        "MANUSCRIPT_DATA_EXTRACTION_MASTER.md",
        "SAILASUTA_2012_COMPLETE_EXTRACTION.md",
        "VALCOUR_2015_COMPLETE_EXTRACTION.md",
        "FINAL_DATA_EXTRACTION_REPORT.md",
        "DATA_EXTRACTION_INDEX.md",
        "EXTRACTION_MISSION_STATUS_REPORT.md",
        "MASTER_DATA_INVENTORY.md",
        "INTEGRATED_MRS_DATA_EXTRACTION.md",
        "WORKFLOW_GUIDE.md",
        "HIV-Associated_Neurocognitive_Disorders__Evidence_for_the_Acute_Phase_Protective_Paradox.md",
        # Model analysis files
        "BAYESIAN_v3_RESULTS_ANALYSIS.md",
        "bayesian_inference_interpretation.md",
        "chronic_NAA_underprediction_analysis.md",
        "RESULTS_ANALYSIS_AND_FIX.md",
        "advanced_modeling_strategies.md",
        "missing_parameters_literature_mining.md",
        # Integration/framework documents
        "COUPLING_FUNCTIONS_SUMMARY.md",
        "EVOLUTIONARY_FRAMEWORK.md",
        "K_STEP_QUANTUM_INTEGRATION.md",
        "spatial_analysis_comprehensive_report.md",
        "microtubule_analysis_report.md",
        # Project management
        "MASTER_INDEX.md",
        "EXECUTIVE_BRIEFING_FINAL.md",
        "FINAL_INTEGRATED_REPORT.md",
        "PROJECT_REVIEW_AND_INTEGRATION_PLAN.md",
        "README_v2.md",
        "executive_summary.md",
    ],

    "papers": [
        "Bairwa_2016.pdf",
        "Boban_2019.pdf",
        "Bolzenius_2023.pdf",
        "Chaganti2021.pdf",
        "Change_2002.pdf",
        "Chelala_2020.pdf",
        "Cohen_2010.pdf",
        "Dahmani_2021.pdf",
        "Heaps_2015.pdf",
        "Young_2014.pdf",
        "sailasuta_2016.pdf",
        "valcour_2015.pdf",
        "jiv2962.pdf",
        "nihms6573182.pdf",
        "nihms697362.pdf",
        "pone_0070164.pdf",
        "pone_01426002.pdf",
        "Supplementary20materialSupplementary_Figure_1.pdf",
        # Your manuscripts
        "Noise_HIV_NC_v1.pdf",
        "Noise_neuro_v1_1.pdf",
        "manuscript_v2_complete.pdf",
        "Review_of__NoiseMediated_Neuroprotection_in_Acute_HIV__A_Computational_Framework_Proposing_Evoluti.pdf",
    ],

    "figures": [
        "EVIDENCE_PYRAMID_COMPREHENSIVE.png",
        "ULTIMATE_COMPREHENSIVE_ANALYSIS.png",
        "bayesian_inference_results.png",
        "coherence_detailed_analysis.png",
        "hiv_phase_comparison.png",
        "microtubule_analysis_overview.png",
        "spatial_comparison_multirun.png",
        "spatial_quantum_dynamics.png",
        "statistical_comparison_multirun.png",
        "valcour_2015_comprehensive_analysis.png",
        # Model output figures
        "bayesian_v3_6_figures.pdf",
        "v3_6_energy.pdf",
        "v3_6_ppc_az.pdf",
        "v3_6_ppc_pred_vs_obs.pdf",
        "v4_posteriors.pdf",
        "v4_pred_vs_obs.pdf",
        # Screenshots
        "Image_111425_at_4_44_PM.jpg",
        "Image_111425_at_5_32_PM.jpg",
        "Screenshot_20251114_at_4_45_33_PM.png",
        "Screenshot_20251114_at_5_34_06_PM.png",
        "Screenshot_20251114_at_5_34_59_PM.png",
        "Screenshot_20251114_at_5_35_24_PM.png",
    ],
}

# Python code files (keep in repo root or models/ directory)
CODE_FILES = [
    "bayesian_optimization.py",
    "bayesian_optimization_v2.py",
    "bayesian_optimization_v2_final.py",
    "bayesian_enzyme_v4.py",
    "final_calibrated_model.py",
    "final_calibrated_model_v2.py",
    "coupling_functions.py",
    "enzyme_kinetics.py",
    "metabolic_dp_connection.py",
    "open_systems.py",
    "shims.py",
    "model_comparison_WORKING.py",
    "model_comparison_clean.py",
    "model_comparison_fixed.py",
    "diagnose_naa_floor.py",
    "geometry_compare.py",
    "geometry_compare_smoke_test.py",
    "sse_demo_config.py",
    "cli.py",
    "__init__.py",
]


def organize_from_downloads():
    """Organize files assuming user has downloaded them to ~/Downloads"""
    downloads = Path.home() / "Downloads"

    print("=" * 70)
    print("HIV MRS PROJECT ORGANIZATION")
    print("=" * 70)
    print(f"\nSource: {downloads}")
    print(f"Destination: {DATA_DIR}\n")

    # Create directories
    for category in FILES_TO_ORGANIZE.keys():
        (DATA_DIR / category).mkdir(parents=True, exist_ok=True)

    stats = {"copied": 0, "missing": 0}

    # Organize each category
    for category, files in FILES_TO_ORGANIZE.items():
        print(f"\n{category.upper()}:")
        dest_dir = DATA_DIR / category

        for filename in files:
            src = downloads / filename
            dest = dest_dir / filename

            if src.exists():
                try:
                    shutil.copy2(src, dest)
                    print(f"  ✓ {filename}")
                    stats["copied"] += 1
                except Exception as e:
                    print(f"  ✗ {filename} (error: {e})")
                    stats["missing"] += 1
            else:
                print(f"  • {filename} (not in Downloads)")
                stats["missing"] += 1

    # Handle code files
    print(f"\nCODE (to repo root or models/):")
    models_dir = REPO_ROOT / "models"
    models_dir.mkdir(exist_ok=True)

    for filename in CODE_FILES:
        src = downloads / filename
        dest = models_dir / filename

        if src.exists():
            try:
                shutil.copy2(src, dest)
                print(f"  ✓ {filename}")
                stats["copied"] += 1
            except Exception as e:
                print(f"  ✗ {filename}")
                stats["missing"] += 1
        else:
            print(f"  • {filename} (not in Downloads)")
            stats["missing"] += 1

    # Summary
    print(f"\n{'=' * 70}")
    print(f"Total copied: {stats['copied']}")
    print(f"Total missing: {stats['missing']}")
    print(f"{'=' * 70}\n")

    # Next steps
    print("NEXT STEPS:")
    print("1. Create data preparation scripts:")
    print("   - data_prep/patient_ratios_prepare.py")
    print("   - data_prep/valcour_prepare.py")
    print("2. Set up .gitignore (protect individual patient data!)")
    print("3. Begin dual pipeline modeling\n")


if __name__ == "__main__":
    organize_from_downloads()