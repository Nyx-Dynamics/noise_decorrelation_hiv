# Changelog

All notable changes to the **Noise-Mediated Neuroprotection in HIV** project will be documented in this file.

## [Milestone] - 2026-01-01 20:45

### Summary
Stabilized the computational environment and fixed critical model regressions. Successfully validated both the Main Bayesian Model (v3.6) and the External Mechanistic Model (v4.0) after environment recovery.

### Fixed
- **PyMC/PyTensor Regression**: Fixed a `TypeError` in `models/bayesian_enzyme_v4.py` by wrapping the ODE-based forward model in a `pytensor.Op`. This allows the Bayesian v4.0 model to sample using mechanistic enzyme kinetics.
- **Python Environment (noiseenv)**: Recreated and stabilized the **Python 3.9** environment (`noiseenv`). Fixed broken pip installation and resolved dependency conflicts.
- **Dependency Versioning**: Updated `requirements.txt` to `pytensor==2.19.0` to maintain compatibility with `pymc==5.12.0` and Python 3.9.
- **Workflow Automation**: Updated the `Makefile` to include `quantum` in the `PYTHONPATH`, ensuring that submodules are correctly imported during model runs.

### Verified
- **v3.6 Main Analysis**: Confirmed end-to-end execution of the hierarchical Bayesian model (`make run-main`).
- **v4.0 Validation**: Confirmed successful sampling and convergence of the mechanistic enzyme model (`make smoke-validation`).

### Added
- **Requirements Recovery**: Restored `requirements.txt` with exact version pins (`pymc==5.12.0`, `pytensor==2.19.0`, etc.) for reproducibility across machines.
- **Mechanistic Validation (v4.0)**: Introduced independent enzyme kinetics model (`bayesian_enzyme_v4.py`) that confirms $\xi$ estimates within 5% of the main Bayesian model without shared architecture.
- **Hierarchical Bayesian Model (v3.6)**: Primary statistical engine with ART-era hierarchical effects and random study effects, establishing $P(\xi_{acute} < \xi_{chronic}) > 0.999$.
- **Environment Hardening**: Locked project to **Python 3.9** to ensure compatibility with PyMC/PyTensor/ArviZ stack.
- **Simplified Makefile**: Intelligent Python detection and consolidated targets for analysis (`run-main`, `run-validation`), figures, and environment checks (`env-check`).
- **Project Documentation**: Comprehensive organizational summaries, roadmaps, and file inventories to support *Nature Communications* submission.

---

## [v2.0] - 2025-10-18 (Approx)

### Added
- **Enhanced Mechanisms**: Astrocyte compensation, nonlinear $\xi$-coherence coupling with floor (0.42 nm), and homeostatic NAA ceilings.
- **Improved Accuracy**: Reduced chronic NAA prediction error from -16% to +2.1%.
- **Validation Plots**: New 4-panel visualization of compensatory effects.

---

## [v1.0] - 2025-10-10 (Approx)

### Added
- Initial Bayesian optimization for noise decorrelation.
- Basic quantum model mapping correlation length ($\xi$) to NAA observables.
- Initial group-level datasets from Sailasuta, Valcour, and Chang studies.
