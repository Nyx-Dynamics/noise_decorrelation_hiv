# Environmental Noise Structure Predicts Neurometabolic Preservation in HIV Infection

![Python](https://img.shields.io/badge/python-3.9+-blue.svg)
![PyMC](https://img.shields.io/badge/PyMC-5.12-orange.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

**A Bayesian Analysis of Neurometabolic Preservation in HIV Infection**

## Overview

This repository contains code, data, and analyses proposing that environmental noise correlation length (ξ) may distinguish protected from vulnerable neurometabolic states in HIV infection.

### The Clinical Paradox

For 40 years, HIV neuroscience has documented an unexplained paradox:
- **Acute HIV** (peak viral load, cytokine storm): >90% maintain normal cognition with preserved NAA
- **Chronic HIV** (suppressed virus, lower inflammation): 40-50% develop HAND

We propose that **noise correlation structure**, not amplitude, determines neurometabolic outcome.

### Key Findings

| Analysis | Finding | Evidence |
|----------|---------|----------|
| **Primary Bayesian (v3.6)** | ξ_acute < ξ_chronic | P > 99.9%, Cohen's d = 5.63 |
| **Enzyme Kinetics (v4)** | Independent confirmation | P > 99% |
| **Individual Patients (v2)** | Pattern holds at patient level | P = 92.4%, n=176 |
| **Cross-Cohort** | Replicates across 3 studies | 2 continents |

### Parameter Estimates

| Parameter | Estimate | 95% HDI |
|-----------|----------|---------|
| ξ_acute | 0.425 ± 0.065 nm | [0.303, 0.541] |
| ξ_chronic | 0.790 ± 0.065 nm | [0.659, 0.913] |
| ξ_healthy | 0.797 ± 0.048 nm | [0.717, 0.887] |
| β_ξ (protection exponent) | 2.33 ± 0.51 | [1.49, 3.26] |

## Repository Structure

```
noise_decorrelation_HIV/
├── quantum/                 # Analysis scripts & external validation
│   ├── bayesian_v3_6_runner.py
│   ├── bayesian_enzyme_v4.py
│   ├── hierarchical_individual_v2_runner.py
│   └── results/
│       ├── enzyme_v4/       # Enzyme kinetics validation
│       └── regional_v1/     # Regional analysis
│
├── results/                 # Main Bayesian inference outputs
│   └── bayesian_v3_6/       # PRIMARY RESULTS
│       ├── summary.csv
│       ├── posterior_predictive.csv
│       └── trace.nc
│
├── data/
│   ├── extracted/           # Group-level MRS data (13 observations)
│   ├── individual/          # Patient-level data (n=62)
│   └── documentation/       # Data extraction methodology
│
├── zenodo_package/          # Prepared for Zenodo deposition
│
├── requirements.txt
├── LICENSE
└── README.md
```

## Quick Start

### Installation

```bash
git clone https://github.com/Nyx-Dynamics/noise_decorrelation_hiv.git
cd noise_decorrelation_hiv
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### View Results

```bash
# Parameter estimates
cat results/bayesian_v3_6/summary.csv

# Model predictions
cat results/bayesian_v3_6/posterior_predictive.csv
```

### Run Analysis

```bash
cd quantum/
python bayesian_v3_6_runner.py
```

See [QUICKSTART.md](QUICKSTART.md) for detailed reproduction instructions.

## Data Sources

All data extracted from published peer-reviewed studies:

| Study | Year | N | Phase |
|-------|------|---|-------|
| Sailasuta et al. | 2012 | 36 | Acute/Chronic/Control |
| Young et al. | 2014 | 90 | Acute/Chronic/Control |
| Sailasuta et al. | 2016 | 59 | Longitudinal |
| Mohamed et al. | 2010 | 35 | Chronic/Control |
| Valcour et al. | 2015 | 62 | Acute (held-out validation) |

## Technical Details

- **Bayesian Inference**: PyMC 5.12.0 with NUTS sampler
- **Convergence**: R̂ < 1.02, ESS 230-418, 0 divergences
- **Validation**: 5-fold CV on held-out cohort (ELPD = 0.532 ± 0.069)

## Citation

```bibtex
@article{demidont2025noise,
  title={Environmental Noise Structure Predicts Neurometabolic Preservation
         in HIV Infection},
  author={Demidont, A.C.},
  journal={Nature Communications},
  year={2025},
  note={In review}
}
```

## License

MIT License - See [LICENSE](LICENSE) for details.

## Contact

**A.C. Demidont, DO**
Nyx Dynamics, LLC
Email: acdemidont@nyxdynamics.org

---

*"The data require some protective mechanism in acute phase; they do not specifically require quantum coherence. But the framework makes testable predictions."*
