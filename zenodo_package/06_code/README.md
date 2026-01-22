# Noise-Mediated Neuroprotection in Acute HIV


**Computational Study suggests mechanism for the 40-year Acute Protective Paradox in HIV Neuroscience**

## Overview

This repository contains code, data, and analyses for demonstrating that quantum coherence mechanisms explain why 70-75% of acute HIV patients maintain normal cognition despite massive neuroinflammation, while chronic patients develop progressive decline.

### Key Findings

- **Main Analysis** (Bayesian v3.6): Definitive evidence that noise correlation length (ξ) differs between acute and chronic HIV (P > 0.999)
- **External Validation** (Enzyme v4): Independent validation via direct enzyme kinetics
- **Regional Analysis**: Evolutionary pattern showing older brain regions have optimal protection
- **Clinical Impact**: First mechanistic explanation for well-established clinical paradox

### Repository Structure

📊 **quantum/** - External validation (enzyme kinetics, regional analysis)  
📈 **results/** - Main Bayesian analysis (v3.6)  
📁 **data/** - All input data (group-level + individual patient data)  

**See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for complete details**

## Quick Start

See [QUICKSTART.md](QUICKSTART.md) for instructions to reproduce all analyses.

## Citation

[Manuscript in preparation for Nature Communications]

## Technical Details

Includes quantum coherence simulator for microtubular structures with stochastic Schrödinger equation (SSE) modeling of Tegmark-style dephasing as spatially-varying, cytokine-modulated noise.
=======
![Python](https://img.shields.io/badge/python-3.9+-blue.svg)
![PyMC](https://img.shields.io/badge/PyMC-5.x-orange.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

**Computational Research Suggesting a Mechanistic Solution for the 40-year Acute Protective Paradox in HIV Neuroscience**

## Overview

This repository contains code, data, and analyses demonstrating that quantum coherence mechanisms explain why 70-75% of acute HIV patients maintain normal cognition despite massive neuroinflammation, while chronic patients develop progressive cognitive decline with lower inflammation levels.

### The Paradox

For 40 years, HIV neuroscience has documented a clinical paradox:
- **Acute HIV** (highest viral loads, peak inflammation): 70-75% maintain normal cognition
- **Chronic HIV** (lower inflammation, effective treatment): Progressive cognitive decline

Despite extensive epidemiological data across 135+ studies and >5,000 patients, no mechanistic explanation existed until now.

### Our Solution

We propose that inflammatory **noise correlation length** (ξ) modulates neuronal metabolism through microtubule quantum coherence:
- **Shorter correlation lengths** (acute HIV) → Paradoxically neuroprotective
- **Longer correlation lengths** (chronic HIV) → Metabolic vulnerability

### Key Findings

1. **Main Analysis (Bayesian v3.6)**: 
   - Definitive evidence: P(ξ_acute < ξ_chronic) > 0.999
   - Protection exponent: β_ξ = 1.89 ± 0.25 (nonlinear mechanism)
   - Model prediction accuracy: <2% error for NAA

2. **External Validation (Enzyme v4)**:
   - Independent enzyme kinetics approach
   - Confirms ξ values within 5% of main model
   - No shared model structure with v3.6

3. **Individual-Level Validation**:
   - n=62 acute HIV patients (Valcour 2015)
   - Statistically significant NAA elevation (+7.7%, p=0.0317)
   - Validates paradox at individual patient level

4. **Regional/Evolutionary Analysis**:
   - Older brain structures (500M years) show optimal protection
   - Newer regions (frontal cortex) maximally vulnerable
   - Suggests evolutionary-scale optimization

## Repository Structure

```
📦 noise_decorrelation_HIV/
├── 📊 quantum/              # External validation & regional analyses
│   └── results/
│       ├── enzyme_v4/      # Independent enzyme kinetics validation
│       ├── regional_v1/    # Evolutionary protection analysis
│       └── model_comparison_*/ # Ablation studies
│
├── 📈 results/              # Main Bayesian inference
│   └── bayesian_v3_6/      # PRIMARY RESULTS FOR MANUSCRIPT
│
└── 📁 data/
    ├── extracted/          # Group-level MRS data (n=3 model inputs)
    ├── individual/         # Patient-level validation data
    └── raw/               # Original source materials
```

**See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for complete documentation.**

## Quick Start

### Setup
```bash
git clone https://github.com/TheonlyqueenAC/noise_decorrelation_HIV.git
cd noise_decorrelation_HIV
python3 -m venv .venv
source .venv/bin/activate
pip install pymc arviz numpy scipy pandas matplotlib seaborn
```

### View Main Results
```bash
# Parameter estimates (Bayesian v3.6)
cat results/bayesian_v3_6/summary.csv

# Model predictions vs observations
cat results/bayesian_v3_6/posterior_predictive.csv

# Figures
open results/ULTIMATE_COMPREHENSIVE_ANALYSIS.png
```

### Run External Validation
```bash
cd quantum/
python bayesian_enzyme_v4.py
# Outputs to quantum/results/enzyme_v4/
```

**See [QUICKSTART.md](QUICKSTART.md) for complete reproduction instructions.**

## Data Sources

All data extracted from published studies:
- **Sailasuta et al. 2012** - RV254 acute HIV cohort MRS
- **Valcour et al. 2015** - SEARCH 010/011 longitudinal study
- **Young et al. 2014** - Perinatal HIV exposure MRS
- **Chang et al. 2002** - Chronic HIV neuroimaging
- **Dahmani et al. 2021** - Meta-analysis (135+ studies)

Complete data available in `data/` directory. See `data/README.md` for details.

## Statistical Evidence

### Primary Hypothesis Test
- **H₀**: ξ_acute ≥ ξ_chronic (no protective mechanism)
- **H₁**: ξ_acute < ξ_chronic (protective paradox)
- **Evidence**: P(H₁ | data) > 0.999
- **Bayes Factor**: BF₁₀ > 1000 (decisive evidence)

### Effect Sizes
- ξ_acute = 0.567 ± 0.068 nm
- ξ_chronic = 0.785 ± 0.073 nm
- Difference: 27.8% shorter in acute (95% HDI: 18.5% - 36.2%)
- Protection exponent: β_ξ = 1.89 ± 0.25

### Model Validation
- NAA prediction error: <2% across all conditions
- Posterior predictive checks: All observations within 94% HDI
- External validation: Independent model confirms findings

## Technical Implementation

### Quantum Coherence Simulator
- Stochastic Schrödinger equation (SSE) framework
- Tegmark-style dephasing with spatially-varying noise
- Cytokine-modulated correlation length
- Microtubule structure modeling

### Bayesian Inference
- PyMC probabilistic programming
- NUTS sampler (4 chains, 2000 draws)
- Comprehensive convergence diagnostics (R̂ < 1.01)
- ArviZ for posterior analysis

### Enzyme Kinetics Integration
- Direct coupling to aspartoacylase activity
- No phenomenological compensation terms
- Independent validation of quantum effects

## Manuscript Status

**Target Journal**: Nature Communications (Impact Factor: 16.6)

**Positioning**: First mechanistic explanation for well-established 40-year clinical paradox, not speculative quantum biology.

**Testable Predictions**: 5 experimental designs with specific costs and timelines.

## Citation

Manuscript in preparation for Nature Communications.

If you use this code or data, please cite:
```
[Citation will be added upon publication]
```

## License

MIT License - See LICENSE file for details.

## Contact

**A.C. Demidont, DO**  
Nyx Dynamics LLC  
Infectious Diseases Physician & Computational Researcher

For questions or collaboration inquiries:
- GitHub Issues: [Open an issue](https://github.com/TheonlyqueenAC/noise_decorrelation_HIV/issues)
- Email: [Your email if you want to include it]

## Acknowledgments

Data from the RV254/SEARCH studies (Sailasuta, Valcour et al.)  
Statistical methodology guidance from MIT MicroMasters program  
Computational infrastructure: [Your institution/resources]

---

*"Chaos can protect. Evolution discovered this 500 million years ago."*
>>>>>>> d68e1f5 (docs: Add comprehensive project documentation)
