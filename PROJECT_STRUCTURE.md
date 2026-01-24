# Project Structure: Noise-Mediated Neuroprotection in Acute HIV

## Overview
This repository contains code, data, and results for the manuscript demonstrating quantum coherence-based mechanisms explaining the acute protective paradox in HIV-associated neurocognitive disorders.

---

## 🛠 Reproducibility Path (Quick Start)

The most straightforward path to reproduce the results is using the provided `Makefile`:

1. **Environment**: `make venv` followed by `make install`
2. **Full Reproducibility**: `make reproduce` (Runs all models and aggregates results)
3. **Main Analysis (v3.6)**: `make run-main` (or `make smoke-main` for a fast check)
4. **Validation (v4.0)**: `make run-validation` (or `make smoke-validation` for a fast check)
5. **Individual Patient Validation**: `make run-hierarchical` (or `make smoke-hierarchical` for a fast check)
6. **Regional Sensitivity Mapping**: `make run-regional` (or `make smoke-regional` for a fast check)
7. **Figures**: `make figures`

---

## Directory Organization

### 📊 `quantum/` - Core Analysis Framework
The root of `quantum/` contains the primary production scripts used for the manuscript.

- `bayesian_v3_6_runner.py`: **MAIN MODEL**. Bayesian inference on group-level MRS data.
- `bayesian_enzyme_v4.py`: Core mechanistic enzyme kinetics model.
- `enzyme_v4_runner.py`: CLI runner for v4.0 validation.
- `regional_hierarchical_v1.py`: **REGIONAL MODEL**. Evolutionary regional sensitivity analysis.
- `hierarchical_individual_v1_runner.py`: Individual patient-level validation (n=62).
- `legacy_figures.py`: Standardized figure generation for the manuscript.

#### Subdirectories:
- `research/`: Experimental simulations, Stochastic Schrödinger Equation (SSE) models, and DTI simulations.
- `utils/`: Data loaders, environment checks, and shared helper functions.
- `legacy/`: Superseded model versions (v2, v3) and old organization scripts.
- `data_prep/`: Scripts for cleaning and preparing raw data into model-ready CSVs.
- `results/`: MCMC traces (.nc), posterior summaries (.csv), and validation plots.

---

### 📈 `results/` - Model Outputs
Primary statistical inference outputs and historical results.

- `bayesian_v3_6/`: Preserved results from the definitive manuscript run.
- `enzyme_v4/`: Outputs from the mechanistic validation.
- `regional_v1/`: Regional sensitivity maps and evolutionary correlation stats.

---

### 📁 `data/` - Input Data
- `raw/`: Original source files (Excel/Word). **Read-only**.
- `extracted/`: Group-level statistics (n=3) used as primary model inputs.
- `individual/`: Patient-level data (n=62) used for independent validation.
- `master/`: Consolidated MRS databases.

---

### 📦 `archive/` - Legacy Materials
Structured graveyard for old scripts, temporary files, and previous project iterations to keep the workspace clean.

#### Subdirectories:
- `code/`: Superseded experimental scripts, diagnostics, and utility scripts.
- `data/`: Legacy data snapshots and compressed archives.
- `docs/`: Previous manuscript drafts, redlines, and legacy overview figures.
- `env/`: Preserved legacy virtual environments (e.g., `noiseenv`).
- `releases/`: Zenodo packages and deployment artifacts.
- `results/`: Archived model outputs from previous iterations.

---

*Last updated: January 2026*
*Copyright (c) 2026 A.C. Demidont, DO, Nyx Dynamics LLC*
