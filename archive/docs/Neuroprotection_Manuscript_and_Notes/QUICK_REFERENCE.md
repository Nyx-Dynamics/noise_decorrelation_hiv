# Enhanced Model v2.0 - Quick Reference Card

## 🚀 Most Common Commands

### Run Enhanced Bayesian Inference
```bash
# Quick test (2 min, 200 samples)
make bayes-v2-smoke

# Standard run (15 min, 3000 samples)
make bayes-v2-run

# Publication quality (25 min, 5000 samples)
make bayes-v2-validate
```

### Run Enhanced Forward Model
```bash
# Validate model vs data
make model-v2-validate

# Generate compensation plots
make model-v2-viz
```

---

## 📁 File Locations

**Code files** (place in project root):
- `bayesian_optimization_v2.py`
- `final_calibrated_model_v2.py`

**Outputs**:
- `results/bayesian_v2/` - Bayesian inference results
- `results/` - Forward model plots

---

## 📊 Check Results

```bash
# View summary
cat results/bayesian_v2/results_summary.txt

# View posterior statistics
cat results/bayesian_v2/summary_v2.csv

# View predictions
cat results/bayesian_v2/posterior_predictive_v2.csv

# View plots (macOS)
open results/bayesian_v2/compensatory_mechanisms.png
open results/enhanced_model_compensation.png
```

---

## ✅ Success Indicators

**After running `make bayes-v2-run`, look for:**
- ✓ P(ξ_acute < ξ_chronic) = 1.000
- ✓ Chronic NAA error: ~+2% (was -16%)
- ✓ All R-hat < 1.05 (convergence)
- ✓ All ESS > 400 (sufficient samples)
- ✓ Astrocyte compensation: ~1.18

---

## 🔧 Customize Parameters

```bash
# More samples
make bayes-v2-run BAYES_V2_DRAWS=5000 BAYES_V2_TUNE=2000

# Higher acceptance rate (if divergences)
make bayes-v2-run BAYES_V2_TARGET_ACCEPT=0.95

# Different random seed
make bayes-v2-run BAYES_V2_SEED=999
```

---

## 🆘 Troubleshooting

### Problem: ImportError for pymc
```bash
# Install dependencies
pip install pymc>=5.0.0 arviz pandas pytensor
```

### Problem: Low ESS
```bash
# Increase samples
make bayes-v2-run BAYES_V2_DRAWS=5000 BAYES_V2_TUNE=2000
```

### Problem: Divergences
```bash
# Increase target acceptance
make bayes-v2-run BAYES_V2_TARGET_ACCEPT=0.95
```

### Problem: Chronic NAA still wrong
```bash
# Check astrocyte compensation parameter
python -c "
import arviz as az
idata = az.from_netcdf('results/bayesian_v2/trace_v2.nc')
az.plot_posterior(idata, var_names=['astrocyte_comp'])
"
```

---

## 📖 Full Documentation

- **INDEX.md** - Master guide to all files
- **QUICKSTART.md** - Detailed quick start
- **IMPLEMENTATION_SUMMARY.md** - What was implemented
- **README_v2.md** - Complete technical docs
- **venv_info.txt** - All available commands

---

## 💡 One-Line Commands

```bash
# Check if Bayesian stack is installed
make check-bayes-env

# Run enhanced Bayesian inference
make bayes-v2-run

# Validate enhanced forward model
make model-v2-validate

# View help
make help
```

---

## 📈 Expected Output

```
P(ξ_acute < ξ_chronic) = 1.0000 ✓

Posterior Medians:
  astrocyte_comp:  1.182 (18% boost)
  xi_floor_nm:     0.42 nm
  xi_ceil_nm:      0.79 nm

Predictions:
  Chronic NAA error: +2.1% ✓
  (Previously: -16.0%)
```

---

## 🎯 Next Steps After Success

1. ✅ Verify convergence (R-hat, ESS)
2. ✅ Check chronic NAA error < 5%
3. ✅ Review compensation mechanisms
4. 📝 Start manuscript writing
5. 🔬 Plan validation experiments

---

**For complete instructions, see INDEX.md or venv_info.txt**
