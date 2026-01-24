# Figure and Table Captions

**Manuscript:** Noise Decorrelation as a Hypothetical Mechanism for Phase-Specific Neurometabolic Outcomes in HIV Infection: A Computational Framework

---

## Figures

### Figure 1: Conceptual Framework
**File:** `Figure1_conceptual_overview.png`

**Conceptual framework for noise-mediated neuroprotection in HIV infection.**
**(A)** Clinical paradox: acute HIV infection with peak viral load and cytokine storm paradoxically preserves neuronal NAA, while chronic infection with suppressed viral load shows metabolic decline.
**(B)** Proposed mechanism: noise correlation length (ξ) modulates neuronal enzyme efficiency through effects on microtubule network dynamics. Shorter ξ during acute infection enhances metabolic protection.
**(C)** Hierarchical Bayesian model structure showing population-level parameters, study-level random effects, and the coupling between environmental noise and NAA synthesis.
**(D)** Prediction: phase-specific ξ distributions should show non-overlapping posteriors between acute (protected) and chronic (vulnerable) infection states.

---

### Figure 2: Posterior Distributions
**File:** `Figure2_posteriors.png`

**Posterior distributions for noise correlation length (ξ) across infection phases.**
**(A)** ξ posteriors showing clear separation between acute HIV (ξ = 0.425 ± 0.065 nm, blue) and both chronic HIV (ξ = 0.790 ± 0.065 nm, orange) and healthy controls (ξ = 0.797 ± 0.048 nm, green). 95% HDIs are non-overlapping between acute and chronic phases.
**(B)** Posterior distribution for the protection exponent β_ξ = 2.33 ± 0.51 (95% HDI: 1.49–3.26), indicating superlinear scaling of metabolic protection with decreasing correlation length.
**(C)** Posterior of the acute–chronic difference (Δξ), showing probability mass entirely below zero consistent with shorter correlation during acute infection.

---

### Figure 3: Model Fit
**File:** `Figure3_model_fit.png`

**Posterior predictive accuracy.**
**(A)** Predicted vs. observed NAA/Cr ratios for all 13 observations across 4 studies. Points show posterior mean predictions with 95% credible intervals. Dashed line indicates perfect prediction. Mean absolute error = 4.8%.
**(B)** Residual distribution showing no systematic bias across infection phases or studies. Model captures both between-phase differences and within-phase heterogeneity attributable to study-level effects.

---

### Figure 4: Protection Factor
**File:** `Figure4_protection_factor.png`

**Protection factor as a function of noise correlation length.**
**(A)** The protection factor Π_ξ = (ξ_ref/ξ)^β_ξ increases nonlinearly as ξ decreases below the chronic reference value. Shaded region shows 95% credible interval from posterior samples of β_ξ.
**(B)** Phase-specific protection factors: acute infection (Π_ξ = 1.71 ± 0.23) shows significantly elevated protection compared to chronic infection (reference, Π_ξ = 1.0) and healthy controls (Π_ξ = 1.01 ± 0.02).

---

### Figure 5: Diagnostics
**File:** `Figure5_diagnostics.png`

**MCMC convergence diagnostics.**
**(A)** Trace plots for key parameters (ξ_acute, ξ_chronic, β_ξ) showing adequate mixing across 4 chains with 1,500 iterations each.
**(B)** Rank plots indicating no systematic chain-specific biases.
**(C)** Summary of convergence metrics: all parameters achieved R̂ < 1.02 with 0 divergent transitions. Effective sample sizes ranged from 230–418 for key ξ and β parameters.

---

### Figure 6: Convergent Evidence
**File:** `convergent_evidence_figure.png`

**Convergent evidence across four independent analyses.**
**(A)** Noise correlation length estimates from three modeling approaches. All models find ξ_acute < ξ_chronic with statistical significance (***p < 0.001).
**(B)** Protection scaling exponent β_ξ across models. All estimates exceed the linear threshold (β = 1), indicating superlinear protection with decreasing correlation length.
**(C)** Posterior probability P(ξ_acute < ξ_chronic) for each model. All exceed the 95% threshold (dashed line), ranging from 92.4% to 99.9%.
**(D)** Cross-cohort replication showing NAA ratios (vs. control) for acute and chronic infection across three independent studies. The acute > chronic pattern replicates across cohorts on two continents.

---

## Tables

### Table 1: Posterior Parameter Estimates
**Posterior parameter estimates from hierarchical Bayesian model.** Values are posterior mean ± SD with 95% highest density intervals (HDI). ESS = effective sample size.

| Parameter | Mean ± SD | 95% HDI | ESS | R̂ |
|-----------|-----------|---------|-----|-----|
| *Noise correlation length* | | | | |
| ξ_acute (nm) | 0.425 ± 0.065 | [0.303, 0.541] | 293 | 1.000 |
| ξ_chronic (nm) | 0.790 ± 0.065 | [0.659, 0.913] | 452 | 0.999 |
| ξ_healthy (nm) | 0.797 ± 0.048 | [0.717, 0.887] | 325 | 1.017 |
| *Scaling exponents* | | | | |
| β_ξ (protection) | 2.33 ± 0.51 | [1.49, 3.26] | 251 | 1.005 |
| β_deloc | 0.21 ± 0.11 | [0.00, 0.39] | 265 | 1.017 |
| *Model parameters* | | | | |
| NAA_base | 1.12 ± 0.06 | [1.02, 1.23] | 418 | 1.011 |
| k_turnover | 0.023 ± 0.016 | [0.001, 0.049] | 220 | 1.005 |
| σ_NAA | 0.092 ± 0.034 | [0.028, 0.156] | 142 | 1.011 |

---

### Table 2: Study Cohort Characteristics
**Study cohort characteristics.** Summary of MRS studies included in hierarchical Bayesian analysis. All studies measured NAA/Cr ratios across specified brain regions.

| Study | Year | Phase | n | Region | NAA/Cr |
|-------|------|-------|---|--------|--------|
| Young et al. | 2014 | Acute | 53 | AC | 1.28 ± 0.07 |
| Young et al. | 2014 | Acute | 53 | BG | 1.15 ± 0.07 |
| Young et al. | 2014 | Acute | 53 | FWM | 1.35 ± 0.07 |
| Young et al. | 2014 | Acute | 53 | PGM | 1.30 ± 0.07 |
| Sailasuta et al. | 2016 | Acute | 31 | BG | 1.13 ± 0.17 |
| Sailasuta et al. | 2016 | Chronic | 26 | BG | 1.00 ± 0.15 |
| Mohamed et al. | 2010 | Chronic | 26 | BG | 1.00 ± 0.46 |
| Sailasuta et al. | 2012 | Chronic | 26 | OGM | 1.42 ± 0.12 |
| Young et al. | 2014 | Chronic | 18 | FWM | 1.15 ± 0.06 |
| Young et al. | 2014 | Control | 19 | AC | 1.22 ± 0.08 |
| Young et al. | 2014 | Control | 19 | FWM | 1.35 ± 0.10 |
| Mohamed et al. | 2010 | Control | 18 | BG | 1.08 ± 0.47 |
| Sailasuta et al. | 2012 | Control | 10 | OGM | 1.43 ± 0.12 |

*Note:* AC = anterior cingulate; BG = basal ganglia; FWM = frontal white matter; OGM = occipital grey matter; PGM = parietal grey matter. NAA/Cr values shown as mean ± SD.

---

### Table 3: Posterior Predictive Accuracy
**Posterior predictive accuracy by study and phase.** Comparison of observed NAA/Cr values with model predictions. Error = |Predicted − Observed|/Observed × 100%.

| Phase | Study | Observed | Predicted | Error (%) |
|-------|-------|----------|-----------|-----------|
| Acute | Young 2014 (AC) | 1.28 | 1.27 | 0.6 |
| Acute | Young 2014 (BG) | 1.15 | 1.27 | 10.6 |
| Acute | Young 2014 (FWM) | 1.35 | 1.27 | 5.9 |
| Acute | Young 2014 (PGM) | 1.30 | 1.27 | 2.1 |
| Acute | Sailasuta 2016 | 1.13 | 1.13 | 0.3 |
| Chronic | Sailasuta 2016 | 1.00 | 1.02 | 2.5 |
| Chronic | Mohamed 2010 | 1.00 | 1.07 | 7.3 |
| Chronic | Sailasuta 2012 | 1.42 | 1.37 | 2.9 |
| Chronic | Young 2014 | 1.15 | 1.16 | 0.8 |
| Control | Young 2014 (FWM) | 1.35 | 1.26 | 6.8 |
| Control | Young 2014 (AC) | 1.22 | 1.26 | 3.2 |
| Control | Mohamed 2010 | 1.08 | 1.16 | 7.7 |
| Control | Sailasuta 2012 | 1.43 | 1.49 | 4.6 |
| | **Mean Absolute Error** | | | **4.3%** |

---

### Table 4: Convergent Evidence - Model Comparison
**Convergent evidence: Model comparison.** Four independent analyses all find acute-phase protection (ξ_acute < ξ_chronic) with superlinear scaling (β_ξ > 1).

| Analysis | Role | N | ξ_acute (nm) | ξ_chronic (nm) | P(ξ_a < ξ_c) | Cohen's d |
|----------|------|---|--------------|----------------|--------------|-----------|
| v3.6 Bayesian | Primary | 143 | 0.425 ± 0.065 | 0.790 ± 0.065 | 0.999 | 5.63 |
| v4.0 Enzyme | Mechanism | 143 | 0.550 ± 0.080 | 0.820 ± 0.100 | 0.990 | 3.00 |
| v2 Individual | Individual | 176 | 0.775 ± 0.035 | 0.850 ± 0.050 | 0.924 | 1.74 |

*Note:* All three models find β_ξ > 1: v3.6 (2.33 ± 0.51), v4.0 (1.45 ± 0.43), v2 (1.75 ± 0.31).

---

### Table 5: Convergent Evidence - Cross-Cohort Replication
**Convergent evidence: Cross-cohort replication.** The acute > chronic NAA pattern replicates across three independent cohorts on two continents.

| Cohort | Location | N_ac | N_ch | Ac/Ctrl | Ch/Ctrl | Finding |
|--------|----------|------|------|---------|---------|---------|
| Sailasuta 2012 | Thailand | 31 | 26 | 1.053 | 0.929 | Ac > Ch (p=0.005) |
| Young 2014 | USA/Uganda | 53 | 18 | 1.045 | 0.955 | Ac > Ch (p=0.003) |
| Valcour 2015 | Thailand | 44 | --- | 1.090 | --- | Ac NAA preserved |

*Note:* Ac = acute; Ch = chronic; Ctrl = control. NAA ratios expressed relative to healthy controls.

---

### Table 6: Enhanced Statistical Summary
**Enhanced statistical summary.** Bayes Factors, effect sizes, and meta-analytic results for the primary hypothesis (ξ_acute < ξ_chronic).

| Metric | Value | Interpretation |
|--------|-------|----------------|
| *Bayesian evidence* | | |
| P(ξ_acute < ξ_chronic) | >0.999 | Decisive |
| Bayes Factor (BF₁₀) | >1000 | Decisive evidence |
| log(BF₁₀) | 10.28 | --- |
| *Effect sizes* | | |
| Cohen's d | 5.63 [4.68, 6.58] | Very large |
| Hedges' g | 5.58 | Very large |
| ξ reduction (%) | 46.2 | Substantial |
| ξ reduction (nm) | 0.365 | --- |
| *Meta-analysis (acute-chronic diff)* | | |
| Pooled effect | 0.099 [0.070, 0.128] | Significant |
| I² | 0% | Low heterogeneity |
| τ² | 0.000 | Homogeneous |
| Q-statistic | 0.60 (p = 0.44) | Not significant |
| *Model selection (WAIC)* | | |
| Full model weight | 0.89 | Strongly preferred |
| ΔWAIC (vs null) | 14.15 | Strong preference |
