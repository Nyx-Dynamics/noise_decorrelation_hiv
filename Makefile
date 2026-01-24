# =============================================================================
# Makefile for Noise-Mediated Neuroprotection in HIV
# =============================================================================

# Default Python interpreter
PY ?= python3
ifneq ($(wildcard .venv/bin/python),)
  PY = .venv/bin/python
else ifneq ($(wildcard noiseenv/bin/python3.9),)
  PY = noiseenv/bin/python3.9
endif

# Environment
VENV_DIR = .venv
VENV_PY = $(VENV_DIR)/bin/python

.PHONY: help venv install clean-venv run-main run-validation run-hierarchical run-regional run-all smoke-main smoke-validation smoke-hierarchical smoke-regional figures report env-check

help:
	@echo "--- Noise-Mediated Neuroprotection Project ---"
	@echo ""
	@echo "Environment:"
	@echo "  make venv             # Create virtual environment (.venv)"
	@echo "  make install          # Install dependencies into .venv"
	@echo "  make clean-venv       # Remove virtual environment"
	@echo "  make env-check        # Check Python version and dependencies"
	@echo ""
	@echo "Main Analysis (v3.6):"
	@echo "  make run-main         # Run main Bayesian v3.6 analysis"
	@echo "  make smoke-main       # Quick test of main Bayesian analysis"
	@echo ""
	@echo "External Validation (v4):"
	@echo "  make run-validation   # Run Enzyme v4 mechanistic validation"
	@echo "  make smoke-validation # Quick test of Enzyme v4 validation"
	@echo ""
	@echo "Individual Patient Validation:"
	@echo "  make run-hierarchical # Run hierarchical individual-level model"
	@echo "  make smoke-hierarchical # Quick test of individual-level model"
	@echo ""
	@echo "Regional Sensitivity Mapping:"
	@echo "  make run-regional     # Run regional hierarchical model"
	@echo "  make smoke-regional   # Quick test of regional model"
	@echo ""
	@echo "Outputs & Reporting:"
	@echo "  make figures          # Generate all manuscript figures"
	@echo "  make report           # Summarize the project status"
	@echo "  make run-all          # Run main, validation, hierarchical, regional, and figures"
	@echo ""
	@echo "Using Python: $(PY)"
	@echo "To override: make PY=python3.9 <target>"

env-check:
	@echo "Checking Python interpreter..."
	@$(PY) --version || echo "Error: $(PY) not found"
	@echo "Checking dependencies..."
	@$(PY) -c "import pymc, pytensor, pandas, numpy, arviz; print('All dependencies found.')" || echo "Error: Missing dependencies. Run 'make install' or check your environment."

# --- Environment ---

venv:
	@echo "Creating virtual environment..."
	python3 -m venv $(VENV_DIR)
	@echo "Done. Activate with: source $(VENV_DIR)/bin/activate"

install:
	@echo "Installing dependencies..."
	$(VENV_PY) -m pip install --upgrade pip
	$(VENV_PY) -m pip install -r requirements.txt
	@echo "Done."

clean-venv:
	rm -rf $(VENV_DIR)
	@echo "Virtual environment removed."

# --- Main Analysis v3.6 ---

run-main:
	@echo "Running Main Bayesian Analysis v3.6..."
	$(PY) -m quantum.bayesian_v3_6_runner --ratio 3_1_1 --era both
	@echo "Results saved to results/bayesian_v3_6/"

smoke-main:
	@echo "Running Smoke Test for Main Analysis..."
	$(PY) -m quantum.bayesian_v3_6_runner --ratio 3_1_1 --era both --draws 100 --tune 100 --chains 2
	@echo "Smoke test complete."

# --- External Validation v4 ---

run-validation:
	@echo "Running Enzyme v4 Mechanistic Validation..."
	$(PY) -m quantum.enzyme_v4_runner --ratio 3_1_1 --era both
	@echo "Results saved to results/enzyme_v4/"

smoke-validation:
	@echo "Running Smoke Test for Enzyme v4..."
	$(PY) -m quantum.enzyme_v4_runner --ratio 3_1_1 --era both --draws 100 --tune 100 --chains 2
	@echo "Smoke test complete."

# --- Individual Patient Validation ---

run-hierarchical:
	@echo "Running Hierarchical Individual-Level Analysis..."
	$(PY) -m quantum.hierarchical_individual_v1_runner
	@echo "Results saved to results/hierarchical_individual_v1/"

smoke-hierarchical:
	@echo "Running Smoke Test for Hierarchical Individual Model..."
	$(PY) -m quantum.hierarchical_individual_v1_runner --draws 100 --tune 100 --chains 2
	@echo "Smoke test complete."

# --- Regional Sensitivity Mapping ---

run-regional:
	@echo "Running Regional Hierarchical Analysis..."
	$(PY) -m quantum.regional_hierarchical_v1
	@echo "Results saved to results/regional_hierarchical_v1/"

smoke-regional:
	@echo "Running Smoke Test for Regional Model..."
	$(PY) -m quantum.regional_hierarchical_v1 --draws 100 --tune 100 --chains 2
	@echo "Smoke test complete."

# --- Figures ---

figures:
	@echo "Generating Manuscript Figures..."
	$(PY) -m quantum.legacy_figures --model bayesian_v3_6 --ratio 3_1_1 --era both --suffix 3_1_1_both_v3_6 --outdir figures/figures
	@echo "Figures generated in figures/figures/"

# --- Reporting ---

report:
	@echo "Project: Noise-Mediated Neuroprotection in HIV"
	@echo "Status: Analysis Complete (v3.6 Main + v4.0 Validation)"
	@echo "Key Files:"
	@echo "  - Main Model: quantum/bayesian_v3_6_runner.py"
	@echo "  - Validation Model: quantum/bayesian_enzyme_v4.py"
	@echo "  - Data Inventory: data/DATA_INVENTORY.md"
	@echo "  - Structure: PROJECT_STRUCTURE.md"

# --- Aggregates ---

run-all: run-main run-validation run-hierarchical run-regional figures
	@echo "Full project workflow completed."

reproduce:
	@echo "Running full reproducibility orchestrator..."
	$(PY) -m quantum.reproduce_all
	@echo "Reproducibility results available in reproducibility_results/"
