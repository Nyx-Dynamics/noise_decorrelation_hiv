# Quick Start Guide

## Reproducibility Path

### 0. Quick Reproduce (All Models)
```bash
make reproduce
```
This runs the main Bayesian analysis, enzyme validation, individual validation, and regional mapping, then aggregates all findings into `reproducibility_results/`.

### 1. Clone Repository
```bash
git clone https://github.com/yourusername/noise_decorrelation_HIV.git
cd noise_decorrelation_HIV
```

### 2. Create Virtual Environment
```bash
python3 -m venv .venv
source .venv/bin/activate  # On Mac/Linux
# .venv\Scripts\activate   # On Windows
```

### 3. Install Dependencies
```bash
pip install pymc arviz numpy scipy pandas matplotlib seaborn
pip install -r requirements.txt  # If you create one
```

---

## Key Analyses

### Main Bayesian Model (v3.6)

Results are in `results/bayesian_v3_6/`

**To run:**
```bash
make run-main
```

**To view results:**
```bash
# View parameter estimates
cat results/bayesian_v3_6/summary.csv
```

---

### External Validation (Enzyme v4)

**To run:**
```bash
make run-validation
```

---

### Regional Sensitivity Mapping

**To run:**
```bash
make run-regional
```

---

### Individual Patient Validation (Hierarchical v1)

**To run:**
```bash
make run-hierarchical
```

---

### Generate All Figures

**To run:**
```bash
make figures
```

**Outputs:**
- `quantum/results/regional_v1/regional_summary.csv`
- `quantum/results/regional_v1/evolutionary_stats.txt`
- `quantum/results/regional_v1/plot*.png` (6 figures)

**Expected runtime:** ~10-15 minutes

---

### Model Comparison (Ablation Studies)

**Run model ablation tests:**
```bash
cd quantum/
python model_comparison_clean.py
```

**Outputs:**
- `quantum/results/model_comparison_clean/*.nc` trace files
- Model comparison statistics (WAIC, LOO)

**Tests:**
- Full model (β_ξ estimated, nonlinear)
- Linear model (β_ξ = 1)
- No coupling model (no ξ effect)

---

## Data Access

### View extracted data (used for models)
```bash
ls data/extracted/
cat data/extracted/CRITICAL_STUDIES_COMPLETE_DATA.csv
```

### View individual patient data (for validation)
```bash
ls data/individual/
head data/individual/VALCOUR_2015_INDIVIDUAL_PATIENTS.csv
```

### View master database
```bash
cat data/master/MASTER_HIV_MRS_DATABASE_v2.csv
```

---

## Generate Figures for Manuscript

### Main text figures
```bash
# Figure 1: Comprehensive analysis
cp figures/ULTIMATE_COMPREHENSIVE_ANALYSIS.png manuscript/figures/Figure1.png

# Figure 2: Bayesian inference
cp results/bayesian_inference_results.png manuscript/figures/Figure2.png

# Figure 3: Mechanism illustration
cp figures/xi_dependence_NAA.png manuscript/figures/Figure3.png

# Figure 4: External validation
cp quantum/results/enzyme_v4/v4_pred_vs_obs.pdf manuscript/figures/Figure4.pdf
```

### Supplementary figures
```bash
# Supplementary Figure 1: Regional analysis
cp quantum/results/regional_v1/plot6_comprehensive_summary.png manuscript/figures/FigureS1.png

# Supplementary Figure 2: Phase comparison
cp figures/hiv_phase_comparison.png manuscript/figures/FigureS2.png

# Supplementary Figure 3: Spatial dynamics
cp figures/spatial_quantum_dynamics.png manuscript/figures/FigureS3.png
```

---

## Extract Statistics for Manuscript

### Parameter estimates (Main text)
```bash
# View Bayesian v3.6 results
python -c "
import pandas as pd
df = pd.read_csv('results/bayesian_v3_6/summary.csv')
print(df[['mean', 'sd', 'hdi_3%', 'hdi_97%']])
"
```

### External validation (Methods)
```bash
# View enzyme v4 results
python -c "
import pandas as pd
df = pd.read_csv('quantum/results/enzyme_v4/summary_v4.csv')
print(df[['mean', 'sd', 'hdi_3%', 'hdi_97%']])
"
```

### Model comparison (Supplementary)
```bash
cat quantum/results/model_comparison_clean/waic_comparison.txt
```

---

## Troubleshooting

### "ModuleNotFoundError: No module named 'pymc'"
```bash
source .venv/bin/activate
pip install pymc arviz
```

### "Cannot open trace.nc file"
NetCDF files require `netcdf4` package:
```bash
pip install netcdf4
```

### Long MCMC runtime
Reduce sampling parameters in scripts:
- Change `draws=2000` to `draws=1000`
- Change `tune=1000` to `tune=500`
- Use fewer chains (2 instead of 4)

### Convergence warnings
Check diagnostics:
```python
import arviz as az
trace = az.from_netcdf('results/bayesian_v3_6/trace.nc')
print(az.summary(trace, hdi_prob=0.94))
az.plot_trace(trace)
```

---

## File Organization Reference

```
📦 noise_decorrelation_HIV/
├── 📊 quantum/              # External validation code
│   └── results/            # Enzyme v4, regional, model comparison
├── 📈 results/             # Main Bayesian v3.6 analysis
│   └── bayesian_v3_6/     # PRIMARY RESULTS
├── 📁 data/                # All input data
│   ├── extracted/         # Group-level stats (model input)
│   ├── individual/        # Patient-level data (validation)
│   └── raw/              # Original source files
└── 📝 *.md                # Documentation
```

**See `PROJECT_STRUCTURE.md` for complete details.**

---

## Citation

If you use this code or data, please cite:

```
[Your manuscript citation once published]
```

Data sources:
- Sailasuta et al. 2012 - DOI: XXX
- Valcour et al. 2015 - DOI: XXX  
- Young et al. 2014 - DOI: XXX
- Chang et al. 2002 - DOI: XXX

---

## Contact

For questions about the analysis or code:
- AC, Nyx Dynamics LLC
- Email: [your email]
- GitHub Issues: [repository issues page]

---

*Last updated: January 2026*
