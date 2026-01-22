# Noise Correlation Length Distinguishes Acute from Chronic HIV Infection: Data and Code Repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A.C. Demidont, DO**
Nyx Dynamics, LLC

## Overview

This Zenodo repository contains the complete data, code, and analysis outputs for the manuscript:

> **"Noise Correlation Length Distinguishes Acute from Chronic HIV Infection: A Bayesian Analysis of Neurometabolic Preservation"**

We propose that environmental noise correlation length (ξ) may distinguish protected from vulnerable neurometabolic states in HIV infection. Using hierarchical Bayesian inference on MRS data from ~220 patients across 4 independent studies, we find evidence consistent with shorter noise correlation during acute infection (ξ_acute = 0.425 ± 0.065 nm) compared to chronic infection (ξ_chronic = 0.790 ± 0.065 nm).

## Repository Structure

```
zenodo_package/
├── 01_manuscript/           # LaTeX source, PDF, and bibliography
│   ├── manuscript_v3.6_nature_comm.tex
│   ├── manuscript_v3.6_nature_comm.pdf
│   ├── manuscript_v3.6_bibtex.tex
│   ├── references.bib
│   └── noise_v2_8_supplemental.pdf
│
├── 02_data_raw/             # Original source data references
│
├── 03_data_processed/       # Harmonized input data for models
│   ├── MASTER_HIV_MRS_DATABASE_v2.csv
│   └── [processed data files]
│
├── 04_model_outputs/        # Bayesian inference results
│   ├── bayesian_v3_6/
│   │   ├── summary.csv          # Posterior estimates
│   │   ├── posterior_predictive.csv
│   │   └── trace.nc             # MCMC chains (NetCDF)
│   └── enzyme_v3/               # Enzyme kinetics validation
│
├── 05_figures/              # Publication figures (PNG)
│   ├── Figure1_conceptual_overview.png
│   ├── Figure2_posteriors.png
│   ├── Figure3_model_fit.png
│   ├── Figure4_protection_factor.png
│   ├── Figure5_diagnostics.png
│   ├── convergent_evidence_figure.png
│   └── [supplementary figures]
│
├── 06_code/                 # Complete source code
│   ├── quantum/             # Analysis scripts
│   ├── scripts/             # Runners and utilities
│   └── requirements.txt     # Python dependencies
│
├── 07_cross_validation/     # Held-out Valcour cohort data
│   ├── VALCOUR_2015_INDIVIDUAL_PATIENTS.csv
│   └── VALCOUR_2015_REGIONAL_SUMMARY.csv
│
├── 08_sensitivity/          # Prior sensitivity analyses
├── 09_supplementary/        # Additional analyses
├── 10_quantum_hand/         # Manual quantum calculations
├── 11_local_development/    # Development notes
│
├── CHECKSUMS.sha256         # File integrity verification
├── LICENSE                  # MIT License
└── README.md                # This file
```

## Key Results

### Primary Findings (Bayesian v3.6)

| Parameter | Estimate | 95% HDI | Interpretation |
|-----------|----------|---------|----------------|
| ξ_acute | 0.425 ± 0.065 nm | [0.303, 0.541] | Shorter correlation (protected) |
| ξ_chronic | 0.790 ± 0.065 nm | [0.659, 0.913] | Longer correlation (vulnerable) |
| ξ_healthy | 0.797 ± 0.048 nm | [0.717, 0.887] | Similar to chronic |
| β_ξ | 2.33 ± 0.51 | [1.49, 3.26] | Superlinear protection scaling |

**Key finding**: P(ξ_acute < ξ_chronic) > 99.9%

### Convergent Evidence

Four independent analyses support the acute-phase protection finding:
1. **Primary Bayesian Model (v3.6)**: P > 99.9%, Cohen's d = 5.63
2. **Enzyme Kinetics Model (v4.0)**: P > 99%, independent parameterization
3. **Individual Patient Analysis (v2)**: P = 92.4%, n=176 measurements
4. **Cross-Cohort Replication**: 3 cohorts, 2 continents, same pattern

## Data Sources

All data extracted from published peer-reviewed studies:

| Study | Year | N | Phase | Region |
|-------|------|---|-------|--------|
| Sailasuta et al. | 2012 | 36 | Acute/Chronic/Control | OGM |
| Young et al. | 2014 | 90 | Acute/Chronic/Control | AC, BG, FWM, PGM |
| Sailasuta et al. | 2016 | 59 | Longitudinal | BG |
| Mohamed et al. | 2010 | 35 | Chronic/Control | BG |
| Valcour et al. | 2015 | 62 | Acute (held-out) | Multiple |

## Reproducibility

### Requirements
- Python 3.9+
- PyMC 5.12.0
- ArviZ 0.16+
- NumPy, SciPy, Pandas, Matplotlib

### Quick Start
```bash
# Install dependencies
pip install -r 06_code/requirements.txt

# View main results
cat 04_model_outputs/bayesian_v3_6/summary.csv

# Run primary analysis (if reproducing)
cd 06_code/quantum
python bayesian_v3_6_runner.py
```

### Verification
```bash
# Verify file integrity
sha256sum -c CHECKSUMS.sha256
```

## Citation

If you use this data or code, please cite:

```bibtex
@article{demidont2025noise,
  title={Noise Correlation Length Distinguishes Acute from Chronic HIV Infection:
         A Bayesian Analysis of Neurometabolic Preservation},
  author={Demidont, A.C.},
  journal={Nature Communications},
  year={2025},
  note={In review}
}
```

## License

This work is licensed under the MIT License. See LICENSE file for details.

## Contact

**A.C. Demidont, DO**
Nyx Dynamics, LLC
Email: acdemidont@nyxdynamics.org
GitHub: https://github.com/Nyx-Dynamics/noise_decorrelation_hiv

## Acknowledgments

- Participants in the original clinical studies
- RV254/SEARCH research teams (Sailasuta, Valcour, et al.)
- PyMC development team

---

*Version: 1.0.0 | Last updated: January 2025*
