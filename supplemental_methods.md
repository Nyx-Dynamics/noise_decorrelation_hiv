# Supplemental Methods: Expanded Dataset (3:1:1) — Noise Decorrelation HIV Models

This document records the exact environment, commands, run IDs, diagnostics, and figure-generation steps used to reproduce the expanded-dataset analyses (3:1:1 primary) for the phenomenological (v3.6, v3) and forward mechanistic ODE (enzyme v3) models. It is written to ensure end-to-end reproducibility from a clean checkout on Apple Silicon (arm64) Macs.

## Environment (Apple Silicon, Conda-forge)

- Hardware/OS: Apple Silicon (arm64) macOS
- Python: 3.12 (arm64)
- Distribution: Miniforge/Conda-forge

Create the environment:

```bash
brew install --cask miniforge
conda init zsh
exec zsh

conda create -n hiv-models -c conda-forge \
  python=3.12 numpy scipy pandas matplotlib pymc arviz pytensor \
  openblas llvm-openmp -y
conda activate hiv-models
```

Set threading for stable performance:

```bash
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
```

## Data and Era Definition

- Primary expanded dataset (phase-balance 3:1:1) files live under:
  - `data/extracted_expanded/data_ratios_comparison/bayesian_inputs_3_1_1.csv`
  - `data/extracted_expanded/data_ratios_comparison/enzyme_inputs_3_1_1.csv`

- Era classification (pre vs post modern ART):
  - `pre_modern`: year ≤ 2006
  - `post_modern`: year ≥ 2007
  - If year is missing, rows are labeled `unknown`.

- The loader attaches `art_era` via a fallback chain (measurement year → publication year → study→year mapping). In this project, a one-off preprocessing step filled missing `publication_year` for known studies (e.g., `Chang_2002` → 2002, `Young_2014` → 2014). This enabled pre/post runs.

Note: rows with `art_era = unknown` are included only in `--era both` analyses; they are excluded from `pre_modern`/`post_modern` splits to avoid imputing era when year information is unavailable.

## Models and Roles

- v3.6 (phenomenological): primary Bayesian model used for manuscript inference and figures. Better posterior predictive fit vs v3.
- v3 (phenomenological): optional parity model; included for completeness in sensitivity or historical comparisons.
- enzyme v3 ODE (mechanistic forward): produces steady-state NAA/Cr and Cho/Cr by phase using the ξ → Π_ξ mapping (Option A), with ART era treated as a covariate for filtering/stratification of inputs, not as a mechanistic modifier.

All Bayesian variants employ Option A mapping exclusively:
- Mechanistic driver: `xi_estimate_nm`
- Derived modulation: `Π_ξ = (xi / xi_baseline)^(−β_ξ)`
- `enzyme_activity_fold` columns are retained for validation only (not used in the forward modulation).

## Primary Converged Run (Both Eras): v3.6, 3:1:1

- Command (example):

```bash
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
python -m quantum.bayesian_v3_6_runner --ratio 3_1_1 --era both \
  --draws 1000 --tune 1000 --chains 4 --seed 2025
```

- Run ID used for figures: 
  - `results/bayesian_v3_6/3_1_1/both/20251123T223826Z_a8c794cc/`

- Diagnostics (from ArviZ):
  - max R-hat ≈ 1.0064 (≤ 1.01)
  - min ESS bulk ≈ 1079
  - min ESS tail ≈ 1216
  - Divergences: none at final settings (or negligible after modest tuning)

- Key outputs used:
  - `trace_v3_6.nc` (InferenceData NetCDF)
  - `posterior_predictive.csv` (observed vs posterior predicted means)
  - `run_manifest.json` (per-run provenance)

- Manuscript figures generated (see Figure Generation section):
  - `forest_hdi_3_1_1_both_v3_6.png`
  - `posterior_beta_xi_kde_3_1_1_both_v3_6.png`
  - `posterior_xi_control_kde_3_1_1_both_v3_6.png`
  - `posterior_xi_acute_kde_3_1_1_both_v3_6.png`
  - `posterior_xi_chronic_kde_3_1_1_both_v3_6.png`
  - `ppd_overlay_3_1_1_both_v3_6.png`
  - `era_effect_violin_3_1_1_both_v3_6.png`

## Pre/Post Era Splits (Sensitivity Note)

- After filling publication years by study, the dataset contained:
  - `pre_modern`: 1 row (Chang 2002)
  - `post_modern`: 27 rows
  - `unknown`: 2 rows

- With `pre_modern` only (N=1), era-specific inference is weakly identified. Short-tuning MCMC exhibited many divergences (as expected). Longer tuning and higher target acceptance (e.g., `--tune 6000`, `--target-accept 0.97`) mitigate but do not fully eliminate divergences due to the single observation.

- We therefore report the primary inference and all manuscript figures from the converged **both-eras** v3.6 run, with ART era treated as a covariate in the combined model. Pre/post panels may be shown as supplementary illustrations with appropriate caveats.

## Forward Mechanistic ODE (enzyme v3)

- Command (3:1:1, both):

```bash
python -m quantum.run_enzyme_v3 --ratio 3_1_1 --era both
```

- Run ID used for ODE figure:
  - `results/enzyme_v3/3_1_1/both/20251123T200507Z_1bdc3701/`
- Figure produced (manuscript):
  - `ode_phase_ratios_3_1_1_both.png`

## Figure Generation (Manuscript Set)

A dedicated generator produces the legacy-style figure set into `figures/figures/` based on the latest (or specified) run.

- Both eras (v3.6, 3:1:1):

```bash
python -m quantum.legacy_figures \
  --model bayesian_v3_6 --ratio 3_1_1 --era both \
  --run-id 20251123T223826Z_a8c794cc \
  --suffix 3_1_1_both_v3_6 --outdir figures/figures
```

- Pre / Post (optional supplemental panels):

```bash
python -m quantum.legacy_figures \
  --model bayesian_v3_6 --ratio 3_1_1 --era pre_modern \
  --suffix 3_1_1_pre_v3_6 --outdir figures/figures

python -m quantum.legacy_figures \
  --model bayesian_v3_6 --ratio 3_1_1 --era post_modern \
  --suffix 3_1_1_post_v3_6 --outdir figures/figures
```

- Forward ODE figure:

```bash
python -m quantum.legacy_figures \
  --model enzyme_v3_ode --ratio 3_1_1 --era both \
  --suffix 3_1_1_both_ode --outdir figures/figures
```

(Optionally add `--mirror-root "/Users/.../runs/run_YYYYMMDD_HHMMSS"` to copy figures into the iCloud manuscript folder.)

## Quick Diagnostics Reproduction for the Converged Run

```bash
python - <<'PY'
import arviz as az
idata = az.from_netcdf('results/bayesian_v3_6/3_1_1/both/20251123T223826Z_a8c794cc/trace_v3_6.nc')
print('max R-hat:', az.rhat(idata).to_array().max().values)
print('min ESS bulk:', az.ess(idata, method='bulk').to_array().min().values)
print('min ESS tail:', az.ess(idata, method='tail').to_array().min().values)
PY
```

## Makefile Targets (Convenience)

To regenerate the manuscript figures with your active interpreter:

```bash
make -f Makefile.txt PY="$(which python)" figs-v3_6-3-1-1-both
make -f Makefile.txt PY="$(which python)" figs-ode-3-1-1-both
# Optional supplemental panels (require pre/post runs to exist)
make -f Makefile.txt PY="$(which python)" figs-v3_6-3-1-1-pre
make -f Makefile.txt PY="$(which python)" figs-v3_6-3-1-1-post
```

For mirroring figures to your iCloud results tree in the same call:

```bash
make -f Makefile.txt PY="$(which python)" \
  LEGACY_OUTPUT_ROOT="/Users/acdmbpmax/Library/Mobile Documents/com~apple~CloudDocs/PycharmProjects - Noise Decorrelation/studio_results_v3_6/results_v3_6/runs/run_20251116_134612" \
  figs-v3_6-3-1-1-both
```

## Reproducibility and Provenance

- Every run writes a per-run `run_manifest.json` capturing CLI args, environment versions, git info, inputs, and outputs.
- All outputs are written into collision-safe, timestamped directories (no overwrites of prior results).
- Figures are regenerated from the exact run folder used for the manuscript via a dedicated CLI.

## Rationale for Using the Combined (Both-Eras) Run in Main Text

Splitting by era left the pre_modern group with a single study-level observation. Small-N strata are inherently unstable for fully Bayesian models with multiple hierarchical components, and NUTS divergences at standard tuning are expected. The combined both-eras model, with ART era as a covariate, converged cleanly (R-hat ≤ 1.01 and strong ESS) and supports the study’s claims without relying on underpowered era splits. The pre/post runs are documented here for transparency and are available as supplemental panels with appropriate caveats.

---

This `supplemental_methods.md` is intended to be included in the manuscript/Zenodo package alongside the run manifests, figures, and (optionally) a `figures_manifest.json` that maps each figure to its originating run_id and CLI.
