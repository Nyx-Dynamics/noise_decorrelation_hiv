# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] — 2026-01-25

### Added
- Primary hierarchical Bayesian model (v3.6) with ART-era effects
- Enzyme kinetics independent validation pathway (v4)
- Individual-level hierarchical model (v1) for Valcour cohort (n = 62)
- Regional hierarchical analysis across 4 brain regions
- 5-model WAIC/LOO comparison framework
- Leave-one-region-out (LORO) sensitivity analysis
- Leave-one-study-out (LOSO) robustness analysis
- K-fold cross-validation on held-out Valcour cohort
- Decoherence baseline simulations (coupled vs uncoupled grid geometries)
- Curated MRS datasets from 4 published studies (~220 patients)
- Comprehensive Makefile with 40+ reproducibility targets
- QUICKSTART.md and PROJECT_STRUCTURE.md documentation
- Full test suite
- MIT license

### Data sources
- Sailasuta et al. 2012 (n = 36)
- Young et al. 2014 (n = 90)
- Sailasuta et al. 2016 (n = 59)
- Mohamed et al. 2010 (n = 35)
- Valcour et al. 2015 (n = 62, held-out validation)

### Key results
- ξ\_acute = 0.425 ± 0.065 nm, ξ\_chronic = 0.790 ± 0.065 nm
- β\_ξ = 2.33 ± 0.51 (superlinear protection scaling)
- P(ξ\_acute < ξ\_chronic) > 0.999
- Cohen's d = 5.63
- WAIC full model weight = 0.89
- Zero divergent transitions, all R̂ < 1.02
