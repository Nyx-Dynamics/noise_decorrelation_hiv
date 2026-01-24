#!/usr/bin/env python3
"""
Zenodo Data Package Preparation Script for HIV Noise Neuroprotection Study
DOI: 10.5281/zenodo.17512732
Author: A.C. Demidont
Date: November 2024
"""

import argparse
import hashlib
import json
import os
import shutil
from datetime import datetime
from pathlib import Path

# Define Zenodo package structure
ZENODO_STRUCTURE = {
    "title": "Data for: Environmental noise correlation as an evolutionary mechanism for neuroprotection during acute HIV infection",
    "version": "v2.8",
    "doi": "10.5281/zenodo.17512732",
    "authors": [
        {
            "name": "Demidont, A.C.",
            "affiliation": "Nyx Dynamics, LLC",
            "orcid": "0000-0002-9216-8569"
        }
    ],
    "description": "Comprehensive data package supporting manuscript on noise-mediated neuroprotection in acute HIV infection using Bayesian inference on basal ganglia MRS data.",
    "keywords": [
        "HIV neuroprotection",
        "noise correlation",
        "Bayesian inference",
        "basal ganglia",
        "magnetic resonance spectroscopy",
        "NAA metabolism",
        "evolutionary mechanisms"
    ],
    "license": "CC-BY-4.0",
    "related_identifiers": [
        {
            "identifier": "https://github.com/nyxdynamics/hiv-noise-neuroprotection",
            "relation": "isSupplementedBy",
            "scheme": "url"
        }
    ]
}


def organize_zenodo_package(base_dir: Path):
    """Organize project files for Zenodo submission."""

    # Create Zenodo directory structure
    base_dir.mkdir(parents=True, exist_ok=True)

    # Define subdirectories
    directories = {
        "01_manuscript": "Final manuscript PDFs and LaTeX source",
        "02_data_raw": "Original MRS data from published sources",
        "03_data_processed": "Processed data and model inputs",
        "04_model_outputs": "Bayesian inference results and summaries",
        "05_figures": "All figures from manuscript and supplements",
        "06_code": "Analysis code and model implementations",
        "07_cross_validation": "K-fold cross-validation results",
        "08_sensitivity": "Sensitivity analyses and ablation tests",
        "09_supplementary": "Supplementary tables and additional analyses"
    }

    for dir_name, description in directories.items():
        dir_path = base_dir / dir_name
        dir_path.mkdir(parents=True, exist_ok=True)

        # Create README for each directory
        readme_content = f"# {dir_name.replace('_', ' ').title()}\n\n{description}\n"
        (dir_path / "README.md").write_text(readme_content)

    print(f"Created Zenodo package structure at: {base_dir}")
    return base_dir


def create_file_manifest():
    """Create comprehensive file manifest for Zenodo."""

    manifest = {
        "creation_date": datetime.now().isoformat(),
        "package_version": "v2.8",
        "files": {}
    }

    # Categorize files from project directory
    file_categories = {
        "manuscript_files": {
            "pattern": ["noise_v2_8_main.pdf", "noise_v2_8_supplemental.pdf", "*.tex"],
            "description": "Manuscript PDFs and LaTeX source files"
        },
        "model_outputs": {
            "pattern": ["summary_*.csv", "results_v3_6_*.csv", "group_*.csv"],
            "description": "Bayesian model outputs and summaries"
        },
        "figures_main": {
            "pattern": ["forest_hdi_*.png", "posterior_*.png", "ppc_bar_*.png", "trace_plots_*.png"],
            "description": "Main manuscript figures"
        },
        "figures_supplementary": {
            "pattern": ["Supplementary_Figure_*.png", "group_constraints_*.png"],
            "description": "Supplementary figures"
        },
        "cross_validation": {
            "pattern": ["cv_*.png", "cv_*.csv", "*valcour*.csv"],
            "description": "Cross-validation results and Valcour cohort analyses"
        },
        "code_files": {
            "pattern": ["*.py", "Makefile"],
            "description": "Python analysis scripts and workflow automation"
        },
        "tables": {
            "pattern": ["TableS*.csv"],
            "description": "Supplementary tables"
        }
    }

    return manifest, file_categories


def create_zenodo_metadata():
    """Create Zenodo metadata.json file."""

    metadata = {
        "metadata": {
            "title": ZENODO_STRUCTURE["title"],
            "upload_type": "dataset",
            "description": ZENODO_STRUCTURE["description"],
            "creators": ZENODO_STRUCTURE["authors"],
            "keywords": ZENODO_STRUCTURE["keywords"],
            "access_right": "open",
            "license": ZENODO_STRUCTURE["license"],
            "related_identifiers": ZENODO_STRUCTURE["related_identifiers"],
            "version": ZENODO_STRUCTURE["version"],
            "language": "eng",
            "communities": [
                {"identifier": "neuroscience"},
                {"identifier": "hiv-aids"},
                {"identifier": "computational-biology"}
            ],
            "grants": [],
            "subjects": [
                {"term": "Neuroscience", "identifier": "https://id.loc.gov/authorities/subjects/sh85091117"},
                {"term": "HIV infections", "identifier": "https://id.loc.gov/authorities/subjects/sh85060965"},
                {"term": "Bayesian statistical decision theory",
                 "identifier": "https://id.loc.gov/authorities/subjects/sh85012506"}
            ]
        }
    }

    return metadata


def create_data_dictionary():
    """Create comprehensive data dictionary for all CSV files."""

    data_dict = {
        "summary_*.csv": {
            "description": "Parameter summaries from Bayesian inference",
            "columns": {
                "parameter": "Model parameter name",
                "mean": "Posterior mean",
                "sd": "Posterior standard deviation",
                "hdi_3%": "3% highest density interval",
                "hdi_97%": "97% highest density interval",
                "mcse_mean": "Monte Carlo standard error of mean",
                "mcse_sd": "Monte Carlo standard error of SD",
                "ess_bulk": "Bulk effective sample size",
                "ess_tail": "Tail effective sample size",
                "r_hat": "Gelman-Rubin convergence diagnostic"
            }
        },
        "group_inputs_*.csv": {
            "description": "Input data for each analysis group",
            "columns": {
                "group": "HIV infection phase (control/acute/chronic)",
                "NAA_Cr_mean": "Mean NAA/Cr ratio",
                "NAA_Cr_sd": "Standard deviation of NAA/Cr",
                "Cho_Cr_mean": "Mean Cho/Cr ratio",
                "Cho_Cr_sd": "Standard deviation of Cho/Cr",
                "n": "Sample size",
                "source": "Data source study"
            }
        },
        "posterior_predictive_*.csv": {
            "description": "Posterior predictive check results",
            "columns": {
                "group": "HIV infection phase",
                "observed": "Observed NAA/Cr value",
                "predicted_mean": "Predicted mean from model",
                "predicted_sd": "Predicted standard deviation",
                "error_percent": "Percentage prediction error",
                "ppc_pvalue": "Posterior predictive p-value"
            }
        },
        "valcour_*.csv": {
            "description": "Valcour cohort specific analyses",
            "columns": {
                "patient_id": "Anonymized patient identifier",
                "week": "Weeks since infection",
                "NAA_Cr": "NAA/Cr ratio measurement",
                "viral_load": "Plasma viral load (copies/mL)",
                "cd4_count": "CD4+ T-cell count",
                "art_status": "Antiretroviral therapy status"
            }
        }
    }

    return data_dict


def create_readme():
    """Create comprehensive README for Zenodo package."""

    readme = """# Data Package: Environmental Noise Correlation in HIV Neuroprotection

## Overview
This data package contains all data, code, and outputs supporting the manuscript:
"Environmental noise correlation as an evolutionary mechanism for neuroprotection during acute HIV infection"

**DOI:** 10.5281/zenodo.17512732  
**Authors:** A.C. Demidont, DO  
**Affiliation:** Nyx Dynamics, LLC  
**Contact:** acdemidont@nyxdynamics.org  

## Package Contents

### 1. Manuscript Files (`01_manuscript/`)
- `noise_v2_8_main.pdf`: Final main manuscript
- `noise_v2_8_supplemental.pdf`: Supplementary materials
- LaTeX source files for reproduction

### 2. Raw Data (`02_data_raw/`)
Published MRS data from three cohorts:
- Sailasuta et al. 2012 (n=36)
- Young et al. 2014 (n=60)
- Valcour et al. 2013 (n=252)

### 3. Processed Data (`03_data_processed/`)
- Harmonized NAA/Cr and Cho/Cr ratios
- Basal ganglia specific measurements
- Cross-platform normalized values

### 4. Model Outputs (`04_model_outputs/`)
- Bayesian inference results (v3.6 phenomenological model)
- Enzyme kinetics validation (v4.0 model)
- Parameter posterior distributions
- Convergence diagnostics

### 5. Figures (`05_figures/`)
All publication figures including:
- Main manuscript figures (1-4)
- Supplementary figures (S1-S5)
- Diagnostic plots
- Cross-validation visualizations

### 6. Analysis Code (`06_code/`)
- `bayesian_v3_6_corrected_local.py`: Main phenomenological model
- `bayesian_enzyme_v4.py`: Enzyme kinetics validation
- `enzyme_kinetics.py`: NAT8L/ASPA pathway implementation
- `Makefile`: Workflow automation

### 7. Cross-Validation (`07_cross_validation/`)
- 5-fold cross-validation results
- ELPD calculations
- Model complexity metrics
- Held-out predictions

### 8. Sensitivity Analyses (`08_sensitivity/`)
- Prior sensitivity tests
- Ablation studies
- Bootstrap stability analyses
- Outlier robustness checks

### 9. Supplementary Materials (`09_supplementary/`)
- TableS1: Literature review summary
- TableS2: Patient characteristics
- TableS3: Regional brain patterns

## Key Results

### Primary Finding
Noise correlation length (ξ) differs between HIV phases:
- Acute: ξ = 0.61 ± 0.09 nm
- Chronic: ξ = 0.81 ± 0.12 nm
- Control: ξ = 0.51 ± 0.08 nm
- P(ξ_acute < ξ_chronic) = 91.4%

### Protection Factor
β = -2.00 ± 0.60 (superlinear scaling)

### Model Accuracy
- Control: Error +2.6%
- Acute HIV: Error -0.4%
- Chronic HIV: Error -2.2%

### Cross-Validation
ELPD = 0.532 ± 0.069 (p < 0.00003)

## Reproducibility

### Software Requirements
- Python 3.11+
- PyMC 5.12.0
- NumPy 1.24+
- SciPy 1.11+
- Matplotlib 3.7+
- ArviZ 0.16+

### Running the Analysis
```bash
# Install dependencies
pip install -r requirements.txt

# Run main analysis
make bayes-v3-final

# Run cross-validation
make cross-validate

# Generate figures
make figures
```

## Data Sources

All clinical data extracted from published literature:
1. Sailasuta N, et al. (2012) PLoS ONE 7:e49272
2. Young AC, et al. (2014) Neurology 83:1592-1600
3. Valcour VG, et al. (2013) PLoS ONE 8:e70164

## License

This data package is licensed under CC-BY 4.0.
You are free to share and adapt with attribution.

## Citation

If you use this data, please cite:
```
Demidont AC (2024). Data for: Environmental noise correlation as an 
evolutionary mechanism for neuroprotection during acute HIV infection. 
Zenodo. https://doi.org/10.5281/zenodo.17512732
```

## Acknowledgments

We thank the participants in the original clinical studies whose data 
made this analysis possible.

## Version History

- v2.8 (2024-11): Final submission version
- v2.0 (2024-10): Major model revision
- v1.0 (2024-09): Initial development

## Contact

For questions about this data package:
A.C. Demidont, DO
acdemidont@nyxdynamics.org
Nyx Dynamics, LLC
"""

    return readme


def create_requirements_txt():
    """Create requirements.txt for reproducibility."""

    requirements = """# Python requirements for HIV noise neuroprotection analysis
# Generated: November 2024

# Core scientific computing
numpy==1.24.3
scipy==1.11.4
pandas==2.0.3

# Bayesian inference
pymc==5.12.0
arviz==0.16.1
pytensor==2.18.6

# Visualization
matplotlib==3.7.2
seaborn==0.12.2

# Statistical analysis
statsmodels==0.14.0
scikit-learn==1.3.0

# Utilities
tqdm==4.66.1
joblib==1.3.2

# Jupyter support (optional)
jupyter==1.0.0
ipykernel==6.25.2
nbformat==5.9.2

# Documentation
sphinx==7.2.6
sphinx-rtd-theme==1.3.0
"""

    return requirements


def create_citation_file():
    """Create CITATION.cff file for GitHub."""

    citation = """cff-version: 1.2.0
title: Environmental noise correlation as an evolutionary mechanism for neuroprotection during acute HIV infection
message: If you use this software, please cite it as below.
type: software
authors:
  - given-names: A.C.
    family-names: Demidont
    email: acdemidont@nyxdynamics.org
    affiliation: Nyx Dynamics, LLC
    orcid: 'https://orcid.org/0000-0000-0000-0000'
identifiers:
  - type: doi
    value: 10.5281/zenodo.17512732
repository-code: 'https://github.com/nyxdynamics/hiv-noise-neuroprotection'
url: 'https://zenodo.org/record/17512732'
abstract: >-
  Computational framework and data supporting the discovery of
  noise-mediated neuroprotection during acute HIV infection.
  Uses Bayesian inference on magnetic resonance spectroscopy
  data to identify evolutionary adaptive mechanisms protecting
  post-mitotic neurons.
keywords:
  - HIV neuroprotection
  - noise correlation
  - Bayesian inference
  - evolutionary mechanisms
  - basal ganglia
license: CC-BY-4.0
version: 2.8
date-released: '2024-11-20'
"""

    return citation


def generate_checksums(directory):
    """Generate SHA256 checksums for all files."""

    checksums = {}

    for root, dirs, files in os.walk(directory):
        for file in files:
            filepath = Path(root) / file
            if filepath.exists():
                with open(filepath, 'rb') as f:
                    file_hash = hashlib.sha256()
                    while chunk := f.read(8192):
                        file_hash.update(chunk)

                relative_path = filepath.relative_to(directory)
                checksums[str(relative_path)] = file_hash.hexdigest()

    return checksums


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare Zenodo data package")
    default_out = Path("/Users/acdmbpmax/Desktop/zenodo_package")
    parser.add_argument("--out", type=Path, default=default_out, help="Output directory for the Zenodo package")
    parser.add_argument("--checksums", action="store_true", help="Generate CHECKSUMS.sha256.json for all files in the package")
    parser.add_argument("--zip", action="store_true", help="Create a zenodo_package.zip archive in the parent directory of --out")
    return parser.parse_args()


def main():
    """Main execution function."""

    print("=" * 60)
    print("ZENODO DATA PACKAGE PREPARATION")
    print("HIV Noise Neuroprotection Study v2.8")
    print("=" * 60)

    args = parse_args()

    # Create package structure
    base_dir = organize_zenodo_package(args.out)

    # Create metadata files
    metadata = create_zenodo_metadata()
    readme = create_readme()
    requirements = create_requirements_txt()
    citation = create_citation_file()
    data_dict = create_data_dictionary()

    # Save files
    (base_dir / "metadata.json").write_text(json.dumps(metadata, indent=2))
    (base_dir / "README.md").write_text(readme)
    (base_dir / "requirements.txt").write_text(requirements)
    (base_dir / "CITATION.cff").write_text(citation)
    (base_dir / "DATA_DICTIONARY.json").write_text(json.dumps(data_dict, indent=2))

    # Create file manifest (template)
    manifest, categories = create_file_manifest()
    (base_dir / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))

    # Optional checksums
    if args.checksums:
        checksums = generate_checksums(base_dir)
        (base_dir / "CHECKSUMS.sha256.json").write_text(json.dumps(checksums, indent=2))
        print("\n✓ Generated CHECKSUMS.sha256.json")

    # Optional zip
    if args.zip:
        zip_path = base_dir.parent / f"{base_dir.name}.zip"
        if zip_path.exists():
            zip_path.unlink()
        shutil.make_archive(str(zip_path.with_suffix('')), 'zip', root_dir=base_dir)
        print(f"✓ Created archive: {zip_path}")

    print("\n✓ Created metadata files")
    print(f"✓ Package location: {base_dir}")
    print("\nNext steps:")
    print("1. Run organize_zenodo.py to copy v2/v3_6/v4 data into the package")
    print("2. Verify presence of noise_v2_8_main.pdf and noise_v2_8_supplemental.pdf in 01_manuscript/")
    print("3. Generate checksums with --checksums (done if you passed the flag)")
    print("4. Create ZIP archive with --zip (optional)")
    print("5. Upload to Zenodo and update DOI")

    return base_dir


if __name__ == "__main__":
    main()