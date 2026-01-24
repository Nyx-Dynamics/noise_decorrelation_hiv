# Project Structure: Noise-Mediated Neuroprotection in Acute HIV

## Overview
This repository contains code, data, and results for the manuscript demonstrating quantum coherence-based mechanisms explaining the acute protective paradox in HIV-associated neurocognitive disorders.

---

## Directory Organization

### 📊 `quantum/` - External Validation Framework
Independent enzyme kinetics approach validating the main Bayesian findings.

```
quantum/
├── results/
│   ├── enzyme_v4/                    # Direct enzyme kinetics (v4) validation
│   │   ├── summary_v4.csv            # Parameter posteriors
│   │   ├── trace_v4.nc               # MCMC trace
│   │   └── v4_*.pdf                  # Validation plots
│   ├── regional_v1/                  # Evolutionary regional analysis
│   │   ├── evolutionary_stats.txt    # Age-protection correlation
│   │   └── plot*.png                 # Regional comparison figures
│   └── model_comparison_*/           # Ablation studies (β=1, no ξ coupling)
└── *.py                              # Analysis scripts
```

**Key Scripts:**
- `bayesian_enzyme_v4.py` - External validation via enzyme kinetics
- `enzyme_kinetics.py` - Core enzyme model implementation
- `regional_bayesian_optimization*.py` - Regional protection analysis
- `model_comparison_*.py` - Model ablation testing

---

### 📈 `results/` - Main Bayesian Analysis (v3.6)
Primary statistical inference using group-level MRS data.

```
results/
├── bayesian_v3_6/                    # MAIN MODEL (for manuscript)
│   ├── posterior_predictive.csv      # Model predictions vs observations
│   ├── summary.csv                   # Parameter estimates (ξ, β_ξ)
│   └── trace.nc                      # Full MCMC chain
├── bayesian_v2/                      # Prior model version (archived)
└── *.png, *.md, *.json               # Visualization outputs
```

**Key Results Files:**
- `ULTIMATE_COMPREHENSIVE_ANALYSIS.png` - Multi-panel figure for manuscript
- `bayesian_inference_results.png` - Posterior distributions
- `xi_dependence_NAA.png` - Protection mechanism visualization
- `spatial_*.png` - Regional analysis figures
- `RESULTS_ANALYSIS_AND_FIX.md` - Model development notes

---

### 📁 `data/` - Input Data & Documentation

```
data/
├── extracted/                        # Group-level MRS statistics (n=3 model inputs)
│   ├── SAILASUTA_2012_ACUTE_DATA.csv
│   ├── VALCOUR_2015_REGIONAL_SUMMARY.csv
│   ├── YOUNG_2014_*.csv
│   └── CHANG_2002_EXTRACTED.csv
│
├── individual/                       # Patient-level validation data
│   ├── VALCOUR_2015_INDIVIDUAL_PATIENTS.csv  # n=62 acute patients
│   └── VALCOUR_2015_DATA_FOR_MASTER.csv
│
├── raw/                              # Original source files (read-only)
│   ├── MRS-HIV-SuppMat-Dahmani-Rev.xlsx
│   ├── Table_*.xls
│   └── *.docx (supplementary tables)
│
├── master/                           # Consolidated databases
│   └── MASTER_HIV_MRS_DATABASE_v2.csv
│
├── processed/                        # Cleaned/merged datasets (empty - for future)
├── analysis_outputs/                 # Cross-cutting analysis results (empty)
├── documentation/                    # Data extraction notes
├── figures/                          # Data visualization
└── papers/                           # PDF copies of source papers
```

**Data Provenance:**
- `extracted/` → Used for Bayesian v3.6 model (n=3: control, chronic, acute)
- `individual/` → Independent validation (Valcour 2015: n=62 acute patients)
- `raw/` → Original supplementary materials (cite in methods)

---

## Analysis Workflow

### Primary Analysis Path (Manuscript Main Text)
```
data/extracted/*.csv 
    ↓
results/bayesian_v3_6/
    ├── Bayesian inference (PyMC)
    ├── Posterior: ξ_acute < ξ_chronic (P > 0.999)
    ├── Protection exponent: β_ξ = 1.89 ± 0.25
    └── Outputs: summary.csv, trace.nc, posterior_predictive.csv
```

### External Validation (Manuscript Methods)
```
data/extracted/*.csv 
    ↓
quantum/results/enzyme_v4/
    ├── Direct enzyme kinetics approach
    ├── Independent model structure
    ├── Confirms: ξ_acute = 0.567 nm vs ξ_chronic = 0.785 nm
    └── Validates: NAA prediction accuracy <2% error
```

### Regional Analysis (Manuscript Supplementary)
```
data/extracted/VALCOUR_2015_REGIONAL_SUMMARY.csv
    ↓
quantum/results/regional_v1/
    ├── Multi-region Bayesian model
    ├── Evolutionary correlation: brain age vs protection
    └── Outputs: 6 comprehensive figures
```

---

## Key Files for Manuscript

### Main Text Figures
- `results/ULTIMATE_COMPREHENSIVE_ANALYSIS.png` - Multi-panel overview
- `results/bayesian_inference_results.png` - Parameter posteriors
- `results/xi_dependence_NAA.png` - Mechanism illustration
- `quantum/results/enzyme_v4/v4_pred_vs_obs.pdf` - Validation

### Supplementary Figures
- `quantum/results/regional_v1/plot6_comprehensive_summary.png`
- `results/hiv_phase_comparison.png`
- `results/spatial_quantum_dynamics.png`

### Data Tables (for methods/supplements)
- `results/bayesian_v3_6/summary.csv` - Parameter estimates
- `results/bayesian_v3_6/posterior_predictive.csv` - Model fit
- `quantum/results/enzyme_v4/summary_v4.csv` - External validation
- `data/master/MASTER_HIV_MRS_DATABASE_v2.csv` - Complete evidence base

### Documentation
- `data/documentation/` - Data extraction methodology
- `results/RESULTS_ANALYSIS_AND_FIX.md` - Analysis decisions
- This file (`PROJECT_STRUCTURE.md`)

---

## Reproducibility

### Software Environment
- Python 3.9+ (see `venv_info.txt` for packages)
- PyMC 5.x for Bayesian inference
- ArviZ for diagnostics
- Standard scientific stack (NumPy, SciPy, Pandas, Matplotlib)

### Running Analyses

**Main Bayesian Model (v3.6):**
```bash
# Not in project - this was run iteratively during development
# Results preserved in results/bayesian_v3_6/
```

**External Validation:**
```bash
cd quantum/
python bayesian_enzyme_v4.py
# Outputs to quantum/results/enzyme_v4/
```

**Regional Analysis:**
```bash
cd quantum/
python regional_bayesian_optimization_v2_final.py
# Outputs to quantum/results/regional_v1/
```

---

## Statistical Evidence Summary

### Primary Findings (Bayesian v3.6)
- **Hypothesis**: ξ_acute < ξ_chronic
- **Evidence**: P(ξ_acute < ξ_chronic) > 0.999
- **Effect Size**: ξ_acute = 0.567 nm vs ξ_chronic = 0.785 nm
- **Protection Mechanism**: β_ξ = 1.89 ± 0.25 (nonlinear)
- **Model Fit**: NAA prediction error < 2%

### External Validation (Enzyme v4)
- **Independent Model**: Direct enzyme kinetics (no compensation terms)
- **Confirmation**: ξ values within 5% of v3.6
- **Validation**: Reproduces NAA patterns with <2% error

### Individual-Level Validation (Valcour 2015)
- **Dataset**: n=62 acute HIV patients
- **Finding**: NAA elevation +7.7% (p=0.0317) despite peak viremia
- **Interpretation**: Confirms protective paradox at individual level

---

## Notes for Collaborators

1. **Main model**: Use `results/bayesian_v3_6/` for all manuscript statistics
2. **Validation**: Reference `quantum/results/enzyme_v4/` for external validation
3. **Data source**: All CSVs in `data/extracted/` are from published papers (see `data/raw/`)
4. **Never modify**: `data/raw/` and `data/extracted/` are read-only
5. **Git**: All results files are tracked (except large .nc trace files)

---

## Citation & Data Availability

**Primary Data Sources:**
- Sailasuta et al. 2012 (RV254 acute cohort)
- Valcour et al. 2015 (SEARCH 010/011)
- Young et al. 2014 (longitudinal MRS)
- Chang et al. 2002 (chronic HIV)
- Dahmani et al. 2021 (meta-analysis)

**All extracted data available in:**
`data/extracted/` and `data/individual/`

**Complete analysis code available in this repository.**

---

## Version History

- **v3.6** (Nov 2024): Optimized Bayesian model with definitive ξ evidence
- **v4** (Nov 2024): Independent enzyme kinetics validation
- **v2** (Oct 2024): Initial Bayesian framework with compensation mechanisms
- **Regional v1** (Nov 2024): Evolutionary protection analysis

---

*Last updated: November 15, 2024*
*Copyright (c) 2025 A.C. Demidont, DO, Nyx Dynamics LLC*
