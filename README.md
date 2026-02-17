![Python](https://img.shields.io/badge/python-3.9+-blue.svg)
![PyMC](https://img.shields.io/badge/PyMC-5.x-orange.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

# Noise-Mediated Neuroprotection in Acute HIV

Computational research suggesting a mechanistic solution for the 40-year Acute Protective Paradox in HIV neuroscience.

## Overview

This repository contains code, data, and analyses demonstrating that quantum coherence mechanisms explain why 70–75% of acute HIV patients maintain normal cognition despite massive neuroinflammation, while chronic patients develop progressive cognitive decline with lower inflammation.

### The Paradox

- Acute HIV (highest viral loads, peak inflammation): 70–75% maintain normal cognition
- Chronic HIV (lower inflammation, effective treatment): Progressive cognitive decline

### Our Solution

Inflammatory noise correlation length (ξ) modulates neuronal metabolism via microtubule quantum coherence:
- Shorter correlation lengths (acute HIV) → neuroprotective (noise decorrelation)
- Longer correlation lengths (chronic HIV) → metabolic vulnerability (correlated noise)

### Key Findings

1. Main Analysis (Bayesian v3.6)
   - P(ξ_acute < ξ_chronic) > 0.999
   - Protection exponent β_ξ = 1.89 ± 0.25
   - NAA error < 2%
2. External Validation (Enzyme v4): independent enzyme kinetics approach confirms ξ
3. Individual-Level Validation (Valcour 2015, n=62): significant acute NAA elevation
4. Regional Analysis: evolution-consistent protection differences

---

## Repository Information

- **URL**: [https://github.com/Nyx-Dynamics/noise_decorrelation_hiv](https://github.com/Nyx-Dynamics/noise_decorrelation_hiv)
- **Organization**: Nyx Dynamics LLC

---

## Where is the code?

Use these entry points depending on what you want to do:

- Publication figures (root):
  - gen_fig2.py – Posterior distributions (needs RUN_DATA_DIR pointing to v3.6 runs)
  - gen_fig3.py – Model fit to group constraints (BG)
  - gen_fig4.py – WAIC/LOO comparison from traces (auto-discovers traces)
  - gen_fig5.py – K-fold cross-validation (auto-discovers latest Valcour CV)

- Main Bayesian model (v3.6):
  - quantum/quantum/bayesian_v3_6_corrected_local.py – primary script to reproduce results
  - quantum/quantum/bayesian_v3_6_original.py – original BG-only variant
  - Makefile targets: bayes-v36, bayes-v36-no-valcour, bayes-v36-regions-*, bayes-v36-loro-auto

- Model comparison utilities:
  - quick_waic.py, quick_waic_loo.py – compute WAIC/PSIS-LOO over traces
  - inspect_waic.py – inspect WAIC/LOO interactively

- Calibrated mechanistic models (quantum → MRS):
  - final_calibrated_model.py (root wrapper) → quantum/final_calibrated_model.py
    - python -m final_calibrated_model (validation + ξ dependence plot)
  - final_calibrated_model_v2.py (root and quantum/) – enhanced with compensation

- Enzyme kinetics (external validation, mechanistic):
  - models/bayesian_enzyme_v4.py – PyMC model using enzyme_kinetics forward model
  - enzyme_kinetics.py (root) – lightweight API used by v4

  Run v4 with the larger group-level dataset (BG ratios):

  ```bash
  # From repo root
  python models/bayesian_enzyme_v4.py \
      --use-group-data \
      --data-csv data/extracted/CRITICAL_STUDIES_COMPLETE_DATA.csv \
      --region Basal_Ganglia \
      --studies Young2014_Spudich,Sailasuta2012 \
      --samples 1000 --chains 4 --target-accept 0.99
  ```
  This ingests the larger dataset, builds the v4 enzyme model consistent with v2’s calibration, and saves results under results/enzyme_v4 (trace_v4_group.nc, summary_v4_group.csv). Cho is optional; if absent in the CSV, the model fits NAA-only.

- Ratios comparison pipeline (NAA/Cr, Cho/Cr):
  - quantum/quantum/pipeline_ratio.py – loads curated ratio datasets and produces results/summary

- Data and results locations:
  - data/extracted, data/individual, data/raw – inputs
  - quantum/quantum/results_v3_6 – run artifacts (traces, summaries, figures)
  - results/bayesian_v3_6 – consolidated outputs for quick viewing

See PROJECT_STRUCTURE.md for a complete tree and descriptions.

---

## Quick Start

See QUICKSTART.md for full reproduction steps. Short version:

```bash
# Create and activate environment
python3 -m venv .venv
source .venv/bin/activate

# Install core dependencies
pip install pymc arviz numpy scipy pandas matplotlib seaborn openpyxl
```

View main results:

```bash
cat results/bayesian_v3_6/summary.csv
cat results/bayesian_v3_6/posterior_predictive.csv
```

Run calibrated model(s):

```bash
python -m final_calibrated_model
python final_calibrated_model_v2.py
```

Generate figures (requires RUN_DATA_DIR for v3.6 traces unless default layout is present):

```bash
export RUN_DATA_DIR="$(pwd)/quantum/quantum/results_v3_6"
python gen_fig2.py
python gen_fig4.py
python gen_fig5.py
```

v3.6 helper Make targets (run from repo root):

```bash
make bayes-v36
make bayes-v36-no-valcour
make bayes-v36-regions-bg
make bayes-v36-regions-all
make bayes-v36-loro-auto
make bayes-v36-waic
```

---

## Last Successful Run (marker)

The last successful v3.6 run establishing the evolutionary hypothesis as valid is:

- Timestamp: 20251117_103400
- Run name: with_valcour_BG_nopseudo
- Marker file: quantum/quantum/results_v3_6/runs/LAST_SUCCESSFUL_RUN.json

Most figure scripts will announce this marker on start when RUN_DATA_DIR points to quantum/quantum/results_v3_6 (default). If your artifacts live elsewhere, set:

```bash
export RUN_DATA_DIR="$(pwd)/quantum/quantum/results_v3_6"
```

## Data Sources

Published studies include Sailasuta 2012, Valcour 2015, Young 2014, Chang 2002, and Dahmani 2021. See data/ for files and provenance.

## Statistical Evidence (summary)

- P(ξ_acute < ξ_chronic) > 0.999, β_ξ ≈ 1.89 ± 0.25
- NAA error < 2%; PPC within ~94% HDI

## License

MIT License – see LICENSE.

## Citation

If you use this research or code in your work, please cite it as:

```bibtex
@software{Demidont_Noise-Mediated_Neuroprotection_in_2024,
  author = {Demidont, A.C.},
  title = {{Noise-Mediated Neuroprotection in Acute HIV: A Quantum Coherence Framework}},
  url = {https://github.com/Nyx-Dynamics/noise_decorrelation_hiv},
  version = {1.0.0},
  year = {2024}
}
```

## Contact

A.C. Demidont, DO — Nyx Dynamics LLC

Open an issue on GitHub for questions.

---

“Chaos can protect. Evolution discovered this 500 million years ago.”
