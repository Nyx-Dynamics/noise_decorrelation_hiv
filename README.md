# Noise Correlation Length Distinguishes Neurometabolic Protection from Vulnerability Across HIV Infection Phases

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![PyMC 5.12](https://img.shields.io/badge/PyMC-5.12-orange.svg)](https://www.pymc.io/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18637896.svg)](https://doi.org/10.5281/zenodo.18637896)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2026.703895-b31b1b.svg)](https://doi.org/10.1101/2026.703895)

**A.C. Demidont, DO** · Nyx Dynamics, LLC · [ORCID 0000-0002-9216-8569](https://orcid.org/0000-0002-9216-8569)

---

## The Clinical Paradox

For 40 years, HIV neuroscience has documented an unexplained paradox affecting an estimated 16–20 million people worldwide:

- **Acute HIV** (peak viral load, maximal neuroinflammation): 80–93% maintain normal cognition with preserved NAA
- **Chronic HIV** (suppressed virus, reduced inflammation): 40–50% develop HIV-associated neurocognitive disorders (HAND)

No existing framework explains why neurometabolic injury *increases* as inflammation *decreases*.

## What This Repository Contains

This repository provides the complete code, data, and reproducibility pipeline for a computational framework proposing that environmental noise **correlation length** (ξ), rather than noise amplitude, determines neurometabolic outcomes across HIV infection phases. Four architecturally independent Bayesian models converge on a single finding: shorter noise correlation in acute infection is associated with neurometabolic preservation, while longer correlation in chronic infection is associated with vulnerability.

### Key Results

| Parameter | Estimate | 95% HDI |
|-----------|----------|---------|
| ξ_acute | 0.425 ± 0.065 nm | [0.303, 0.541] |
| ξ_chronic | 0.790 ± 0.065 nm | [0.659, 0.913] |
| β_ξ (coupling exponent) | 2.33 ± 0.51 | [1.49, 3.26] |
| P(ξ_acute < ξ_chronic) | > 0.999 | — |
| Cohen's d | 5.63 | — |

### Convergent Validation

| Model | Approach | Key Finding |
|-------|----------|-------------|
| v3.6 (Primary) | Hierarchical Bayesian, group-level | P(ξ_acute < ξ_chronic) > 99.9% |
| v1.0 (Individual) | Patient-level trajectories (n=62) | P = 95%, 0 divergences |
| v2.0 (Necessity) | Mechanism necessity testing | Null model rejected (WAIC Δ > 10) |
| v4.0 (Enzyme) | Independent enzyme kinetics | P > 99% |

## Repository Structure

```
noise_decorrelation_hiv/
│
├── quantum/                    # Core analysis scripts
│   ├── bayesian_v3_6_runner.py       # Primary hierarchical Bayesian model
│   ├── bayesian_enzyme_v4.py         # Enzyme kinetics validation
│   ├── hierarchical_individual_v1_runner.py  # Individual-level model
│   └── regional_hierarchical_v1.py   # Regional analysis
│
├── data/
│   ├── extracted/              # Group-level MRS data (13 observations, 4 cohorts)
│   ├── individual/             # Patient-level data (n=62, Valcour 2015)
│   └── documentation/          # Data extraction methodology
│
├── results/                    # Inference outputs
│   └── bayesian_v3_6/          # PRIMARY RESULTS
│       ├── summary.csv               # Posterior parameter estimates
│       ├── posterior_predictive.csv   # Model predictions vs observed
│       └── trace.nc                  # Full MCMC trace (NetCDF)
│
├── models/                     # Model specifications
├── figures/                    # Publication figures
├── sections/                   # Manuscript LaTeX sections
│
├── reproduce_all.py            # One-command reproduction pipeline
├── Makefile                    # make reproduce
├── requirements.txt            # Python dependencies
├── CITATION.cff                # Citation metadata
├── CHANGELOG.md                # Version history
└── LICENSE                     # MIT
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

### Reproduce All Results

```bash
python reproduce_all.py
# or
make reproduce
```

### View Primary Results

```bash
cat results/bayesian_v3_6/summary.csv
cat results/bayesian_v3_6/posterior_predictive.csv
```

### Run Individual Models

```bash
cd quantum/
python bayesian_v3_6_runner.py          # Primary model (~5 min)
python bayesian_enzyme_v4.py            # Enzyme kinetics validation
python hierarchical_individual_v1_runner.py  # Individual-level
```

## Data Sources

All data extracted from published peer-reviewed studies. No novel clinical data were collected.

| Study | Year | N | Phase | Region |
|-------|------|---|-------|--------|
| Sailasuta et al. | 2012 | 36 | Acute/Chronic/Control | Basal ganglia, Frontal |
| Young et al. | 2014 | 90 | Acute/Chronic/Control | Basal ganglia, Parietal, Frontal |
| Sailasuta et al. | 2016 | 59 | Longitudinal | Basal ganglia, Frontal |
| Mohamed et al. | 2010 | 35 | Chronic/Control | Frontal white matter |
| Valcour et al. | 2015 | 62 | Acute (held-out validation) | Basal ganglia |

## Technical Validation

### Bayesian Inference (Model v3.6)

- **0 divergent transitions** across all MCMC chains
- **R̂ < 1.02** for all parameters (4 independent chains)
- **ESS 230–418** for primary parameters (ξ, β_ξ)
- **5-fold cross-validation**: ELPD = 0.532 ± 0.069 on held-out data

### Reproducibility

- 12/12 automated tests passing
- Cryptographic checksums for all data files
- Complete `make reproduce` pipeline from raw data to figures
- Zenodo archive: [10.5281/zenodo.18637896](https://doi.org/10.5281/zenodo.18637896)

## Preprint

Demidont, A.C. (2026). Noise Correlation Length Distinguishes Neurometabolic Protection from Vulnerability Across HIV Infection Phases. *bioRxiv*. [https://doi.org/10.1101/2026.703895](https://doi.org/10.1101/2026.703895)

## Citation

```bibtex
@article{Demidont2026_noise_correlation,
  author    = {Demidont, A.C.},
  title     = {Noise Correlation Length Distinguishes Neurometabolic Protection
               from Vulnerability Across {HIV} Infection Phases},
  year      = {2026},
  doi       = {10.1101/2026.703895},
  publisher = {Cold Spring Harbor Laboratory},
  journal   = {bioRxiv}
}
```

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contact

**A.C. Demidont, DO**
Nyx Dynamics, LLC
Email: acdemidont@nyxdynamics.org
ORCID: [0000-0002-9216-8569](https://orcid.org/0000-0002-9216-8569)
