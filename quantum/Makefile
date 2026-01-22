.PHONY: help run-local run-corr validate smoke tegmark-compare phase-sweep demo-config demo-run viz-coherence viz-summary viz-phase-sweep viz-kernel viz-latest-coherence viz-latest-summary viz-latest-kernel viz-geometry compare-geometry mc-run mc-viz interpret-run interpret-phase interpret-mc mc-smoke bayes-v2-run bayes-v2-smoke bayes-v2-validate model-v2-validate model-v2-viz

PY=python

help:
	@echo "Microtubule_Simulation — terminal commands"
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
	$(PY) -m quantum.bayesian_optimization_v2 \
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
	$(PY) -m quantum.bayesian_optimization_v2 \
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
	$(PY) -m quantum.bayesian_optimization_v2 \
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
	$(PY) -c "from final_calibrated_model_v2 import validate_model_v2, plot_compensation_effects; validate_model_v2(); plot_compensation_effects()"
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
