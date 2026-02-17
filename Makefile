.PHONY: help run-local run-corr validate smoke tegmark-compare phase-sweep demo-config demo-run viz-coherence viz-summary viz-phase-sweep viz-kernel viz-latest-coherence viz-latest-summary viz-latest-kernel viz-geometry compare-geometry mc-run mc-viz interpret-run interpret-phase interpret-mc mc-smoke bayes-v2-run bayes-v2-smoke bayes-v2-validate model-v2-validate model-v2-viz bayes-v36 bayes-v36-original bayes-v36-no-valcour bayes-v36-validate bayes-v36-aux-time bayes-v36-aux-time-vl bayes-v36-aux-plasma bayes-v36-aux-csf bayes-v36-aux-both bayes-v36-overlay bayes-v36-regions-bg bayes-v36-regions-all bayes-v36-waic bayes-v36-loro-auto bayes-v36-ic-compare bayes-v36-loso bayes-v36-valcour-cv

PY=python

help:
	@echo "Microtubule_Simulation — terminal commands"
	@echo ""
	@echo "=== BAYESIAN v3.6 (BG ratios; VL-free primary) ==="
	@echo "  make bayes-v36                 # Run v3.6 with Valcour included (timestamped run folder)"
	@echo "  make bayes-v36-original       # Run preserved ORIGINAL v3.6 (BG-only; no multi-region; VL-free primary)"
	@echo "  make bayes-v36-no-valcour      # Run v3.6 excluding Valcour acute individuals"
	@echo "  make bayes-v36-overlay TRACE=path  # Re-run v3.6 and overlay densities vs a given NetCDF trace"
	@echo "  make bayes-v36-validate        # Validation (p-values + Fig.2) and exit"
	@echo "  make bayes-v36-regions-bg      # Run v3.6 with curated group constraints restricted to BG only (--regions BG)"
	@echo "  make bayes-v36-regions-all     # Run v3.6 with curated group constraints from all regions (--regions all)"
	@echo "  make bayes-v36-waic            # Compute WAIC/LOO for traces (TRACE_LIST=path1.nc,path2.nc or auto)"
	@echo "  make bayes-v36-loro-auto       # Leave-one-region-out sweep (BG,FWM,FGM,AC) on curated constraints"
	@echo "  make bayes-v36-ic-compare TRACE_WITH=... TRACE_NO=...  # Compare IC with SE (core, harmonized)"
	@echo "  make bayes-v36-loso            # Leave-one-source-out IC sweep (curated sources)"
	@echo "  make bayes-v36-valcour-cv K=5  # K-fold CV on Valcour acute (predictive alignment)"
	@echo "  make bayes-v36-aux-time        # Run Valcour auxiliary (time-only), decoupled"
	@echo "  make bayes-v36-aux-time-vl     # Run Valcour auxiliary (time+VL), decoupled"
	@echo "  make bayes-v36-aux-plasma      # Run Valcour auxiliary (plasma VL arm), decoupled"
	@echo "  make bayes-v36-aux-csf         # Run Valcour auxiliary (CSF VL arm), decoupled"
	@echo "  make bayes-v36-aux-both        # Run Valcour auxiliary (plasma+CSF arms), decoupled"
	@echo ""
	@echo "=== ENHANCED MODEL v2.0 (NEW) ==="
	@echo "  make bayes-v2-run              # Run enhanced Bayesian inference (3000 samples)"
	@echo "  make bayes-v2-smoke            # Quick test of enhanced Bayesian inference"
	@echo "  make bayes-v2-validate         # Full validation run (5000 samples for publication)"
	@echo "  make model-v2-validate         # Validate enhanced forward model vs data"
	@echo "  make model-v2-viz              # Run forward model and generate compensation plots"
	@echo ""
	@echo "=== ORIGINAL TARGETS ==="
	@echo "  make venv                # Create .venv (Python virtual environment)"
	@echo "  make install             # Install requirements into .venv"
	@echo "  make run-local           # Run SSE_local with default small grid"
	@echo "  make run-corr            # Run SSE_correlated with xi=0.8"
	@echo "  make validate            # Run Extra/sse_validation.py"
	@echo "  make smoke               # Run Extra/sse_smoke_test.py"
	@echo "  make tegmark-compare SUMMARY=path DELTA_X=1e-9  # Compare to Tegmark estimate"
	@echo "  make phase-sweep [K=8 phases='none art_controlled chronic acute'] [MODE=SSE_local] [XI=0.8]"
	@echo "  make demo-config         # Print recommended demo config"
	@echo "  make demo-run            # Run small-grid demo"
	@echo ""
	@echo "=== VISUALIZATION ==="
	@echo "  make viz-coherence CSV=path      # Plot SSE coherence time series from CSV"
	@echo "  make viz-summary SUMMARY=path     # Plot final |psi|^2 and Gamma overlay from summary/NPZ"
	@echo "  make viz-phase-sweep              # Plot variance bands from sse_phase_sweep_summary.json"
	@echo "  make viz-kernel KERNEL=path       # Preview SSE correlated kernel from NPZ"
	@echo "  make viz-geometry CSV=path        # Plot reg vs fib coherence from CSV (+Δ panel)"
	@echo "  make viz-latest-coherence         # Auto-plot coherence for the latest *_sse_coherence.csv"
	@echo "  make viz-latest-summary           # Auto-plot summary overlay for the latest *_summary.json"
	@echo ""
	@echo "=== MONTE CARLO & INTERPRETATION ==="
	@echo "  make compare-geometry SUMMARY=path  # Run geometry comparison (fits + report) for a summary JSON"
	@echo "  make mc-run                      # Run a small Monte Carlo ensemble (override vars)"
	@echo "  make mc-viz                      # Visualize Monte Carlo aggregate summary JSON"
	@echo "  make mc-smoke                    # Run a tiny MC smoke test"
	@echo "  make interpret-run SUMMARY=path   # Create Markdown/JSON interpretation for a single run"
	@echo "  make interpret-phase              # Interpret phase-sweep aggregate JSON"
	@echo "  make interpret-mc                 # Interpret Monte Carlo aggregate JSON"
	@echo ""
	@echo "=== ORIGINAL BAYESIAN ==="
	@echo "  make check-bayes-env    # Check PyMC/ArviZ installation"
	@echo "  make bayes-run          # Run original Bayesian inference"
	@echo "  make bayes-smoke        # Quick test of original Bayesian inference"
	@echo ""
	@echo "=== UTILITIES ==="
	@echo "  make commit-msg         # Generate git commit message"
	@echo "  make pip-freeze         # Export exact package versions"
	@echo "  make clean-venv         # Remove virtual environment"

# =============================================================================
# BAYESIAN v3.6 (BG RATIOS; VL-FREE PRIMARY) — NEW
# =============================================================================

V36_SCRIPT=quantum/quantum/bayesian_v3_6_corrected_local.py
V36_ORIGINAL=quantum/quantum/bayesian_v3_6_original.py
WAIC_HELPER=quantum/quantum/utils/waic_loo_helper.py
IC_COMPARE=compare_ic_with_se.py
LOSO_HELPER=quantum/quantum/utils/loso_ic.py
VALCOUR_CV=quantum/quantum/utils/valcour_kfold_cv.py

# Default run with Valcour included (Acute pool includes Valcour 0/4/12/24)
bayes-v36:
	$(PY) $(V36_SCRIPT) \
		--tag with_valcour \
		--plot-densities --plots-extended

# Preserved ORIGINAL (O2): BG-only, no multi-region augmentation, VL-free primary
bayes-v36-original:
	$(PY) $(V36_ORIGINAL) \
		--tag original \
		--plot-densities --plots-extended

# Ablation: exclude Valcour individuals from Acute pool
bayes-v36-no-valcour:
	$(PY) $(V36_SCRIPT) \
		--exclude-valcour \
		--tag no_valcour \
		--plot-densities --plots-extended

# Regions scope: restrict curated group constraints to BG only
bayes-v36-regions-bg:
	$(PY) $(V36_SCRIPT) \
		--regions BG \
		--tag regions_bg \
		--plot-densities --plots-extended

# Regions scope: include curated group constraints from all regions (default)
bayes-v36-regions-all:
	$(PY) $(V36_SCRIPT) \
		--regions all \
		--tag regions_all \
		--plot-densities --plots-extended

# Compute WAIC/LOO from one or more trace NetCDF files.
# Usage:
#   make bayes-v36-waic TRACE_LIST=path1.nc,path2.nc[,path3.nc]
# If TRACE_LIST is empty, auto-discovers up to the 5 most recent traces under results_v3_6/runs/.
bayes-v36-waic:
	@if [ -z "$(TRACE_LIST)" ]; then \
		$(PY) $(WAIC_HELPER) --auto 5; \
	else \
		$(PY) $(WAIC_HELPER) --traces $(TRACE_LIST); \
	fi

# IC compare with SE (harmonized core)
# Usage:
#   make bayes-v36-ic-compare TRACE_WITH=path/to/with.nc TRACE_NO=path/to/no.nc
bayes-v36-ic-compare:
	@if [ -z "$(TRACE_WITH)" ] || [ -z "$(TRACE_NO)" ]; then \
		echo "ERROR: Provide TRACE_WITH=... and TRACE_NO=..."; exit 1; \
	fi
	$(PY) $(IC_COMPARE) --with-trace $(TRACE_WITH) --no-trace $(TRACE_NO)

# Leave-one-source-out (curated sources) — orchestrates runs and IC comparison
bayes-v36-loso:
	$(PY) $(LOSO_HELPER)

# Valcour K-fold CV (acute hold-out predictive alignment)
# Usage: make bayes-v36-valcour-cv K=5 SEED=123
K?=5
SEED?=42
bayes-v36-valcour-cv:
	$(PY) $(VALCOUR_CV) --kfold $(K) --seed $(SEED)

# Leave-one-region-out sweep for curated group constraints (adjunctive, multi-region):
# Regions in priority order: BG, FWM, FGM, AC (CG maps to AC in code)
bayes-v36-loro-auto:
	@for R in BG FWM FGM AC; do \
		echo "[LORO] Excluding region $$R"; \
		$(PY) $(V36_SCRIPT) \
			--regions all \
			--exclude-region $$R \
			--tag loro_$$R \
			--plot-densities --plots-extended; \
	done

# Validation pass (non-Bayesian): p-values and Fig.2-style plots, then exit
bayes-v36-validate:
	$(PY) $(V36_SCRIPT) \
		--validate-valcour --fig2 --week-pvals \
		--tag valcour_validation

# Auxiliary model (fully decoupled), time-only
bayes-v36-aux-time:
	$(PY) $(V36_SCRIPT) \
		--tag valcour_time \
		--valcour-aux time

# Auxiliary model (fully decoupled), time+VL (falls back to time if VL unavailable)
bayes-v36-aux-time-vl:
	$(PY) $(V36_SCRIPT) \
		--tag valcour_time_vl \
		--valcour-aux time+vl

# Auxiliary model (fully decoupled), plasma VL arm
bayes-v36-aux-plasma:
	$(PY) $(V36_SCRIPT) \
		--tag valcour_plasma \
		--valcour-aux plasma

# Auxiliary model (fully decoupled), CSF VL arm
bayes-v36-aux-csf:
	$(PY) $(V36_SCRIPT) \
		--tag valcour_csf \
		--valcour-aux csf

# Auxiliary model (fully decoupled), plasma + CSF arms
bayes-v36-aux-both:
	$(PY) $(V36_SCRIPT) \
		--tag valcour_both \
		--valcour-aux both

# Re-run v3.6 (with Valcour) and overlay densities vs a provided NetCDF trace
# Usage: make bayes-v36-overlay TRACE=quantum/quantum/results_v3_6/runs/<run_no_valcour>/trace_no_valcour.nc
bayes-v36-overlay:
	@if [ -z "$(TRACE)" ]; then \
		echo "Usage: make bayes-v36-overlay TRACE=path/to/trace_no_valcour.nc"; \
		echo "Hint: run 'make bayes-v36-no-valcour' first to generate the ablation trace."; \
		exit 2; \
	fi
	$(PY) $(V36_SCRIPT) \
		--tag with_valcour \
		--plot-densities --plots-extended \
		--compare-trace $(TRACE)

# =============================================================================
# ENHANCED MODEL v2.0 TARGETS (NEW)
# =============================================================================

# Enhanced Bayesian inference with compensatory mechanisms
BAYES_V2_DRAWS?=3000
BAYES_V2_TUNE?=1500
BAYES_V2_CHAINS?=4
BAYES_V2_TARGET_ACCEPT?=0.92
BAYES_V2_SEED?=42

bayes-v2-run:
	@echo "==================================================================="
	@echo " Running Enhanced Bayesian Inference v2.0"
	@echo " - Astrocyte compensation parameter"
	@echo " - Nonlinear ξ-coherence coupling with floor"
	@echo " - Homeostatic NAA ceiling"
	@echo "==================================================================="
	@echo ""
	$(PY) quantum/bayesian_optimization_v2.py \
		--draws $(BAYES_V2_DRAWS) \
		--tune $(BAYES_V2_TUNE) \
		--chains $(BAYES_V2_CHAINS) \
		--target-accept $(BAYES_V2_TARGET_ACCEPT) \
		--seed $(BAYES_V2_SEED)
	@echo ""
	@echo "✓ Inference complete. Results saved to results/bayesian_v2/"
	@echo ""
	@echo "Check convergence with:"
	@echo "  - All R-hat < 1.05? (in summary_v2.csv)"
	@echo "  - All ESS > 400? (in summary_v2.csv)"
	@echo "  - Chronic NAA error < 5%? (in results_summary.txt)"

bayes-v2-smoke:
	@echo "Running quick smoke test of enhanced Bayesian inference..."
	$(PY) quantum/bayesian_optimization_v2.py \
		--draws 200 \
		--tune 200 \
		--chains 2 \
		--target-accept 0.90 \
		--seed 999
	@echo ""
	@echo "✓ Smoke test complete. Check results/bayesian_v2/ for outputs."

bayes-v2-validate:
	@echo "==================================================================="
	@echo " Running FULL VALIDATION (Publication Quality)"
	@echo " This will take 20-30 minutes..."
	@echo "==================================================================="
	@echo ""
	$(PY) quantum/bayesian_optimization_v2.py \
		--draws 5000 \
		--tune 2000 \
		--chains 4 \
		--target-accept 0.95 \
		--seed 42 \
		--plot
	@echo ""
	@echo "✓ Full validation complete!"
	@echo ""
	@echo "Publication checklist:"
	@echo "  □ P(ξ_acute < ξ_chronic) > 0.95"
	@echo "  □ Chronic NAA error < ±5%"
	@echo "  □ All R-hat < 1.05"
	@echo "  □ All ESS > 400"
	@echo "  □ Astrocyte compensation: 1.15-1.25"
	@echo ""
	@echo "See results/bayesian_v2/results_summary.txt for details"

# Enhanced forward model validation
model-v2-validate:
	@echo "==================================================================="
	@echo " Validating Enhanced Forward Model v2.0"
	@echo "==================================================================="
	@echo ""
	$(PY) -c "from quantum.final_calibrated_model_v2 import validate_model_v2; validate_model_v2()"
	@echo ""
	@echo "✓ Validation complete"

model-v2-viz:
	@echo "==================================================================="
	@echo " Generating Compensation Mechanism Visualizations"
	@echo "==================================================================="
	@echo ""
	$(PY) -c "from quantum.final_calibrated_model_v2 import validate_model_v2, plot_compensation_effects; validate_model_v2(); plot_compensation_effects()"
	@echo ""
	@echo "✓ Plots saved to results/enhanced_model_compensation.png"

# =============================================================================
# ORIGINAL TARGETS (UNCHANGED)
# =============================================================================

# Basic runs via Python CLI
run-local:
	$(PY) -m quantum.cli --mode SSE_local --hiv_phase acute --N_r 36 --N_z 36 --dt 0.01 --time_steps 120 --frames_to_save 12 --rng_seed 1234

run-corr:
	$(PY) -m quantum.cli --mode SSE_correlated --xi 0.8 --hiv_phase acute --N_r 36 --N_z 36 --dt 0.01 --time_steps 120 --frames_to_save 12 --rng_seed 1234

# Validation helpers
validate:
	PYTHONPATH=. $(PY) Extra/sse_validation.py

smoke:
	PYTHONPATH=. $(PY) Extra/sse_smoke_test.py

# Analytical comparison (requires SUMMARY path)
SUMMARY?=
DELTA_X?=1e-9
tegmark-compare:
	@if [ -z "$(SUMMARY)" ]; then echo "ERROR: Provide SUMMARY=path/to/*_summary.json"; exit 1; fi
	$(PY) Extra/tegmark_compare.py $(SUMMARY) --delta_x $(DELTA_X)

# Phase sweep (override variables as needed)
K?=8
phases?=none art_controlled chronic acute
MODE?=SSE_local
XI?=0.8
N_R?=36
N_Z?=36
DT?=0.01
time_steps?=120
frames_to_save?=12
phase-sweep:
	PYTHONPATH=. $(PY) Extra/sse_phase_sweep.py --K $(K) --phases $(phases) --mode $(MODE) --xi $(XI) --N_r $(N_R) --N_z $(N_Z) --dt $(DT) --time_steps $(time_steps) --frames_to_save $(frames_to_save)

# Demo config and run
demo-config:
	PYTHONPATH=. $(PY) Extra/sse_demo_config.py

demo-run:
	PYTHONPATH=. $(PY) Extra/sse_demo_config.py --run

# Visualization targets
viz-coherence:
	@if [ -z "$(CSV)" ]; then echo "ERROR: Provide CSV=path/to/*_sse_coherence.csv"; exit 1; fi
	PYTHONPATH=. $(PY) Extra/sse_visualize.py coherence --csv $(CSV)

viz-summary:
	@if [ -z "$(SUMMARY)" ]; then echo "ERROR: Provide SUMMARY=path/to/*_summary.json"; exit 1; fi
	PYTHONPATH=. $(PY) Extra/sse_visualize.py summary --summary $(SUMMARY)

viz-phase-sweep:
	PYTHONPATH=. $(PY) Extra/sse_visualize.py phase --summary-json datafiles/sse_phase_sweep_summary.json

viz-kernel:
	@if [ -z "$(KERNEL)" ]; then echo "ERROR: Provide KERNEL=path/to/*_sse_kernel.npz"; exit 1; fi
	PYTHONPATH=. $(PY) Extra/sse_visualize.py kernel --kernel $(KERNEL)

viz-geometry:
	@if [ -z "$(CSV)" ]; then echo "ERROR: Provide CSV=path/to/*_sse_coherence.csv"; exit 1; fi
	PYTHONPATH=. $(PY) Extra/sse_visualize.py geometry --csv $(CSV)

compare-geometry:
	@if [ -z "$(SUMMARY)" ]; then echo "ERROR: Provide SUMMARY=path/to/*_summary.json"; exit 1; fi
	# Attempt to infer CSV from SUMMARY; pass if present
	@CSV_INFER=$${SUMMARY%_summary.json}_sse_coherence.csv; \
	if [ -f "$$CSV_INFER" ]; then \
	  echo "Using coherence CSV: $$CSV_INFER"; \
	  PYTHONPATH=. $(PY) Extra/geometry_compare.py run --summary $(SUMMARY) --csv $$CSV_INFER; \
	else \
	  PYTHONPATH=. $(PY) Extra/geometry_compare.py run --summary $(SUMMARY); \
	fi

viz-latest-coherence:
	@if ls datafiles/*_sse_coherence.csv >/dev/null 2>&1; then \
		CSV=$$(ls -t datafiles/*_sse_coherence.csv | head -n1); \
		echo "Using latest CSV: $$CSV"; \
		PYTHONPATH=. $(PY) Extra/sse_visualize.py coherence --csv $$CSV; \
	else echo "No *_sse_coherence.csv found under datafiles/"; exit 1; fi

viz-latest-summary:
	@if ls datafiles/*_summary.json >/dev/null 2>&1; then \
		SUMMARY=$$(ls -t datafiles/*_summary.json | head -n1); \
		echo "Using latest SUMMARY: $$SUMMARY"; \
		PYTHONPATH=. $(PY) Extra/sse_visualize.py summary --summary $$SUMMARY; \
	else echo "No *_summary.json found under datafiles/"; exit 1; fi

viz-latest-kernel:
	@if ls datafiles/*_sse_kernel.npz >/dev/null 2>&1; then \
		KERNEL=$$(ls -t datafiles/*_sse_kernel.npz | head -n1); \
		echo "Using latest KERNEL: $$KERNEL"; \
		PYTHONPATH=. $(PY) Extra/sse_visualize.py kernel --kernel $$KERNEL; \
	else echo "No *_sse_kernel.npz found under datafiles/. Note: kernels are only saved for SSE_correlated runs."; exit 1; fi

# Monte Carlo analytics
MC_N?=16
MC_MODE?=SSE_local
MC_PHASES?=none art_controlled chronic acute
MC_G0_MIN?=0.03
MC_G0_MAX?=0.07
MC_A_MIN?=0.08
MC_A_MAX?=0.12
MC_XI_MIN?=0.4
MC_XI_MAX?=1.2
MC_N_R?=36
MC_N_Z?=36
MC_DT?=0.01
MC_STEPS?=120
MC_FRAMES?=12
MC_BASE_SEED?=3000
MC_OUT?=datafiles/sse_mc_summary.json

mc-run:
	PYTHONPATH=. $(PY) Extra/sse_mc_analytics.py run --N $(MC_N) --mode $(MC_MODE) --phases $(MC_PHASES) \
		--Gamma0_min $(MC_G0_MIN) --Gamma0_max $(MC_G0_MAX) --alpha_min $(MC_A_MIN) --alpha_max $(MC_A_MAX) \
		--xi_min $(MC_XI_MIN) --xi_max $(MC_XI_MAX) --N_r $(MC_N_R) --N_z $(MC_N_Z) --dt $(MC_DT) \
		--time_steps $(MC_STEPS) --frames_to_save $(MC_FRAMES) --base_seed $(MC_BASE_SEED) --out $(MC_OUT)

mc-viz:
	PYTHONPATH=. $(PY) Extra/sse_mc_analytics.py viz --summary $(MC_OUT)

mc-smoke:
	PYTHONPATH=. $(PY) Extra/sse_mc_smoke_test.py

# Interpretation targets
interpret-run:
	@if [ -z "$(SUMMARY)" ]; then echo "ERROR: Provide SUMMARY=path/to/*_summary.json"; exit 1; fi
	PYTHONPATH=. $(PY) Extra/sse_interpret.py run --summary $(SUMMARY)

interpret-phase:
	PYTHONPATH=. $(PY) Extra/sse_interpret.py phase --summary-json datafiles/sse_phase_sweep_summary.json

interpret-mc:
	PYTHONPATH=. $(PY) Extra/sse_interpret.py mc --summary datafiles/sse_mc_summary.json

# Virtual environment management
VENV_DIR:=.venv
PY_VENV:=$(VENV_DIR)/bin/python
PIP_VENV:=$(VENV_DIR)/bin/pip

venv:
	@echo "Creating virtual environment in $(VENV_DIR)..."
	python3 -m venv $(VENV_DIR)
	@echo "Activate with: source $(VENV_DIR)/bin/activate (Linux/macOS)"
	@echo "On Windows PowerShell: .\\$(VENV_DIR)\\Scripts\\Activate.ps1"

install: venv
	@echo "Installing requirements into $(VENV_DIR)..."
	$(PIP_VENV) install --upgrade pip
	$(PIP_VENV) install -r requirements.txt
	@echo "Done. Use 'source $(VENV_DIR)/bin/activate' before running make targets if your default PY isn't the venv."

pip-freeze:
	$(PIP_VENV) freeze > requirements.lock.txt
	@echo "Locked dependencies written to requirements.lock.txt"

clean-venv:
	rm -rf $(VENV_DIR)
	@echo "Removed $(VENV_DIR)"


# Bayesian stack helpers (original)
check-bayes-env:
	PYTHONPATH=. $(PY) Extra/check_bayes_env.py

bayes-run:
	$(PY) -m quantum.bayesian_optimization --draws 2000 --tune 1000 --chains 4 --target-accept 0.9

bayes-smoke:
	$(PY) -m quantum.bayesian_optimization --draws 200 --tune 200 --chains 2 --target-accept 0.9 --seed 123


# Commit message helper
commit-msg:
	PYTHONPATH=. $(PY) Extra/commit_message.py
