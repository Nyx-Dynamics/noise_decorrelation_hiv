#!/usr/bin/env python3
"""
Zenodo Data Package Preparation Script for HIV Noise Neuroprotection Study
DOI: 10.5281/zenodo.17512732
Author: A.C. Demidont
Date: November 2024

LOCAL MAC VERSION - Fixed paths for your system
"""

import hashlib
import json
import shutil
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
    "license": "CC-BY-4.0"
}


def organize_zenodo_package():
    """Organize project files for Zenodo submission."""

    # Create Zenodo directory structure in current directory or Desktop
    # Check if we're in the project directory
    current_dir = Path.cwd()

    # Create zenodo_package in the current directory
    base_dir = current_dir / "zenodo_package"

    print(f"Creating Zenodo package at: {base_dir}")
    base_dir.mkdir(exist_ok=True)

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
        dir_path.mkdir(exist_ok=True)

        # Create README for each directory
        readme_content = f"# {dir_name.replace('_', ' ').title()}\n\n{description}\n"
        (dir_path / "README.md").write_text(readme_content)

    print(f"✓ Created Zenodo package structure at: {base_dir}")
    return base_dir


def copy_project_files(base_dir):
    """Copy files from current project to Zenodo structure."""

    current_dir = Path.cwd()
    copied_count = 0

    # File mappings
    file_mappings = {
        "01_manuscript": ["*.pdf", "*.tex"],
        "02_data_raw": ["TableS*.csv"],
        "03_data_processed": ["group_inputs*.csv", "valcour*.csv", "summary*.csv"],
        "04_model_outputs": ["results_v3_6*.csv", "posterior_predictive*.csv", "waic_loo*.csv"],
        "05_figures": ["*.png"],
        "06_code": ["*.py", "Makefile", "requirements.txt"],
        "07_cross_validation": ["cv_*.png", "cv_*.csv"],
        "08_sensitivity": ["*nopseudo*.csv", "*nopseudo*.png"],
        "09_supplementary": ["aux_*.csv", "*waic*.csv"]
    }

    for directory, patterns in file_mappings.items():
        target_dir = base_dir / directory
        for pattern in patterns:
            files = list(current_dir.glob(pattern))
            for file in files:
                if file.is_file() and not file.name.startswith('.'):
                    target_file = target_dir / file.name
                    try:
                        shutil.copy2(file, target_file)
                        print(f"  ✓ Copied {file.name} → {directory}/")
                        copied_count += 1
                    except Exception as e:
                        print(f"  ✗ Could not copy {file.name}: {e}")

    print(f"\n✓ Copied {copied_count} files to Zenodo package")
    return copied_count


def create_metadata_files(base_dir):
    """Create Zenodo metadata files."""

    # Create metadata.json
    metadata = {
        "metadata": {
            "title": ZENODO_STRUCTURE["title"],
            "upload_type": "dataset",
            "description": ZENODO_STRUCTURE["description"],
            "creators": ZENODO_STRUCTURE["authors"],
            "keywords": ZENODO_STRUCTURE["keywords"],
            "access_right": "open",
            "license": ZENODO_STRUCTURE["license"],
            "version": ZENODO_STRUCTURE["version"],
            "language": "eng"
        }
    }

    with open(base_dir / "metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)

    # Create README.md
    readme = f"""# Data Package: Environmental Noise Correlation in HIV Neuroprotection

## Overview
This data package contains all data, code, and outputs supporting the manuscript:
"Environmental noise correlation as an evolutionary mechanism for neuroprotection during acute HIV infection"

**DOI:** {ZENODO_STRUCTURE['doi']}  
**Authors:** {ZENODO_STRUCTURE['authors'][0]['name']}  
**Affiliation:** {ZENODO_STRUCTURE['authors'][0]['affiliation']}  
**ORCID:** {ZENODO_STRUCTURE['authors'][0]['orcid']}  

## Key Results
- Noise correlation length (ξ) differs between HIV phases:
  - Acute: ξ = 0.61 ± 0.09 nm
  - Chronic: ξ = 0.81 ± 0.12 nm
  - P(ξ_acute < ξ_chronic) = 91.4%
- Protection factor: β = -2.00 ± 0.60 (superlinear scaling)
- Model accuracy: <2% error across all conditions

## Package Contents
See individual directory README files for detailed descriptions.

## Reproducibility
All code and data needed to reproduce the analysis are included.
See 06_code/README.md for installation and usage instructions.

## License
{ZENODO_STRUCTURE['license']}

## Citation
If you use this data, please cite:
Demidont AC (2024). Data for: Environmental noise correlation as an 
evolutionary mechanism for neuroprotection during acute HIV infection. 
Zenodo. https://doi.org/{ZENODO_STRUCTURE['doi']}
"""

    with open(base_dir / "README.md", 'w') as f:
        f.write(readme)

    # Create CITATION.cff
    citation = f"""cff-version: 1.2.0
title: {ZENODO_STRUCTURE['title']}
message: If you use this dataset, please cite it as below.
type: dataset
authors:
  - given-names: A.C.
    family-names: Demidont
    email: acdemidont@nyxdynamics.org
    affiliation: {ZENODO_STRUCTURE['authors'][0]['affiliation']}
    orcid: 'https://orcid.org/{ZENODO_STRUCTURE['authors'][0]['orcid']}'
identifiers:
  - type: doi
    value: {ZENODO_STRUCTURE['doi']}
version: {ZENODO_STRUCTURE['version']}
date-released: '2024-11-20'
"""

    with open(base_dir / "CITATION.cff", 'w') as f:
        f.write(citation)

    print("✓ Created metadata files (metadata.json, README.md, CITATION.cff)")


def generate_checksums(directory):
    """Generate SHA256 checksums for all files."""

    checksums = {}
    directory = Path(directory)

    for file_path in directory.rglob('*'):
        if file_path.is_file() and not file_path.name.startswith('.'):
            with open(file_path, 'rb') as f:
                file_hash = hashlib.sha256()
                while chunk := f.read(8192):
                    file_hash.update(chunk)

            relative_path = file_path.relative_to(directory)
            checksums[str(relative_path)] = file_hash.hexdigest()

    # Save checksums
    with open(directory / "checksums.txt", 'w') as f:
        for path, checksum in sorted(checksums.items()):
            f.write(f"{checksum}  {path}\n")

    print(f"✓ Generated checksums for {len(checksums)} files")
    return checksums


def main():
    """Main execution function."""

    print("=" * 60)
    print("ZENODO DATA PACKAGE PREPARATION")
    print("HIV Noise Neuroprotection Study v2.8")
    print("=" * 60)
    print()

    # Check current directory
    current_dir = Path.cwd()
    print(f"Current directory: {current_dir}")
    print()

    # Create package structure
    base_dir = organize_zenodo_package()

    # Create metadata files
    create_metadata_files(base_dir)

    # Copy project files
    print("\nCopying project files...")
    copied_count = copy_project_files(base_dir)

    # Generate checksums
    print("\nGenerating checksums...")
    generate_checksums(base_dir)

    # Calculate package size
    total_size = sum(f.stat().st_size for f in base_dir.rglob('*') if f.is_file())
    size_mb = total_size / (1024 * 1024)

    print("\n" + "=" * 60)
    print("PACKAGE COMPLETE!")
    print("=" * 60)
    print(f"✓ Package location: {base_dir}")
    print(f"✓ Total files: {copied_count}")
    print(f"✓ Package size: {size_mb:.1f} MB")
    print("\nNext steps:")
    print("1. Review the package contents")
    print("2. Add any missing files manually")
    print("3. Create ZIP archive:")
    print(f"   cd {base_dir.parent}")
    print(f"   zip -r zenodo_package.zip zenodo_package/")
    print("4. Upload to Zenodo and verify DOI")
    print("\nRemember to add:")
    print("  - noise_v2_8_main.pdf → 01_manuscript/")
    print("  - noise_v2_8_supplemental.pdf → 01_manuscript/")

    return base_dir


if __name__ == "__main__":
    main()