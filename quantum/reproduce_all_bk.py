#!/usr/bin/env python3
"""
reproduce_all.py — Single-command reproduction of all manuscript results.

Noise Decorrelation as a Hypothetical Mechanism for Phase-Specific
Neurometabolic Outcomes in HIV Infection: A Computational Framework

Demidont AC (2026). bioRxiv doi: 10.64898/2026.02.10.703895

Usage:
    python reproduce_all.py              # Full reproduction (~2-3 hours)
    python reproduce_all.py --quick      # Smoke test (~5 minutes)
    python reproduce_all.py --core       # Core results only (~45 minutes)

Requirements:
    pip install -r requirements.txt
    (PyMC 5.12+, ArviZ 0.16+, Python 3.9+)
"""

import subprocess
import sys
import time
import argparse
from pathlib import Path
from datetime import datetime


ROOT = Path(__file__).resolve().parent
LOG_FILE = ROOT / "reproducibility_results" / "reproduce_all_log.txt"


def log(msg, log_fh=None):
    """Print and optionally log a message."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    if log_fh:
        log_fh.write(line + "\n")
        log_fh.flush()


def run_step(name, cmd, log_fh=None, cwd=None):
    """Run a step, report pass/fail, return success bool."""
    log(f"{'='*60}", log_fh)
    log(f"STEP: {name}", log_fh)
    log(f"CMD:  {cmd}", log_fh)
    log(f"{'='*60}", log_fh)

    t0 = time.time()
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=cwd or str(ROOT),
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour max per step
        )
        elapsed = time.time() - t0

        if result.returncode == 0:
            log(f"  PASS ({elapsed:.1f}s)", log_fh)
            if log_fh and result.stdout:
                lines = result.stdout.strip().split("\n")
                for line in lines[-20:]:
                    log_fh.write(f"    {line}\n")
            return True
        else:
            log(f"  FAIL (exit {result.returncode}, {elapsed:.1f}s)", log_fh)
            if result.stderr:
                for line in result.stderr.strip().split("\n")[-10:]:
                    log(f"    STDERR: {line}", log_fh)
            return False

    except subprocess.TimeoutExpired:
        log(f"  TIMEOUT after 3600s", log_fh)
        return False
    except Exception as e:
        log(f"  ERROR: {e}", log_fh)
        return False


def main():
    parser = argparse.ArgumentParser(description="Reproduce all manuscript results")
    parser.add_argument("--quick", action="store_true",
                        help="Smoke test only (~5 min)")
    parser.add_argument("--core", action="store_true",
                        help="Core results only (~45 min)")
    args = parser.parse_args()

    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    py = sys.executable
    results = []

    with open(LOG_FILE, "w") as log_fh:
        log("=" * 60, log_fh)
        log("NOISE DECORRELATION HIV — FULL REPRODUCTION SUITE", log_fh)
        log(f"Mode: {'quick' if args.quick else 'core' if args.core else 'full'}", log_fh)
        log(f"Python: {py}", log_fh)
        log(f"Working dir: {ROOT}", log_fh)
        log("=" * 60, log_fh)

        # ============================================================
        # QUICK MODE — smoke tests only
        # ============================================================
        if args.quick:
            steps = [
                ("Bayesian v2 smoke test",
                 f"{py} quantum/bayesian_optimization_v2.py "
                 f"--draws 200 --tune 200 --chains 2 "
                 f"--target-accept 0.90 --seed 999"),
                ("SSE smoke test",
                 f"PYTHONPATH=. {py} quantum/sse_smoke_test.py"),
            ]
            for name, cmd in steps:
                ok = run_step(name, cmd, log_fh)
                results.append((name, ok))

        else:
            # ========================================================
            # CORE STEPS (always run in core and full modes)
            # ========================================================

            # Step 1: Primary Bayesian v3.6 (Figures 1-2, Table 1)
            results.append(("Bayesian v3.6 (primary)", run_step(
                "Primary Bayesian v3.6 (with Valcour)",
                f"{py} quantum/quantum/bayesian_v3_6_corrected_local.py "
                f"--tag with_valcour --plot-densities --plots-extended",
                log_fh,
            )))

            # Step 2: Enzyme kinetics v4 (independent validation)
            results.append(("Enzyme kinetics v4", run_step(
                "Enzyme kinetics v4 (independent validation)",
                f"{py} quantum/bayesian_enzyme_v4.py",
                log_fh,
            )))

            # Step 3: Enhanced Bayesian v2 (publication quality)
            results.append(("Bayesian v2 validate", run_step(
                "Enhanced Bayesian v2 (publication quality)",
                f"{py} quantum/bayesian_optimization_v2.py "
                f"--draws 5000 --tune 2000 --chains 4 "
                f"--target-accept 0.95 --seed 42 --plot",
                log_fh,
            )))

            # Step 4: WAIC/LOO model comparison (Figure 6)
            results.append(("WAIC comparison", run_step(
                "WAIC/LOO 5-model comparison",
                f"{py} quantum/quantum/utils/waic_loo_helper.py --auto 5",
                log_fh,
            )))

            # Step 5: Validation pass (p-values, Fig.2)
            results.append(("Validation", run_step(
                "Validation pass (p-values + Fig.2)",
                f"{py} quantum/quantum/bayesian_v3_6_corrected_local.py "
                f"--validate-valcour --fig2 --week-pvals "
                f"--tag valcour_validation",
                log_fh,
            )))

            if not args.core:
                # ====================================================
                # FULL MODE — sensitivity and robustness analyses
                # ====================================================

                # Step 6: Ablation (exclude Valcour)
                results.append(("No-Valcour ablation", run_step(
                    "Ablation: exclude Valcour cohort",
                    f"{py} quantum/quantum/bayesian_v3_6_corrected_local.py "
                    f"--exclude-valcour --tag no_valcour "
                    f"--plot-densities --plots-extended",
                    log_fh,
                )))

                # Step 7: Basal ganglia only
                results.append(("BG-only", run_step(
                    "Region restriction: basal ganglia only",
                    f"{py} quantum/quantum/bayesian_v3_6_corrected_local.py "
                    f"--regions BG --tag regions_bg "
                    f"--plot-densities --plots-extended",
                    log_fh,
                )))

                # Step 8: Leave-one-region-out
                results.append(("LORO sweep", run_step(
                    "Leave-one-region-out sweep (BG, FWM, FGM, AC)",
                    "make bayes-v36-loro-auto",
                    log_fh,
                )))

                # Step 9: Leave-one-study-out
                results.append(("LOSO", run_step(
                    "Leave-one-study-out",
                    f"{py} quantum/quantum/utils/loso_ic.py",
                    log_fh,
                )))

                # Step 10: K-fold cross-validation
                results.append(("Valcour 5-fold CV", run_step(
                    "5-fold CV on Valcour cohort",
                    f"{py} quantum/quantum/utils/valcour_kfold_cv.py "
                    f"--kfold 5 --seed 42",
                    log_fh,
                )))

                # Step 11: Decoherence baselines
                results.append(("Decoherence (uncoupled)", run_step(
                    "Decoherence baseline (uncoupled grid)",
                    f"{py} -m quantum.cli --mode SSE_local "
                    f"--hiv_phase acute --N_r 36 --N_z 36 "
                    f"--dt 0.01 --time_steps 120 --frames_to_save 12 "
                    f"--rng_seed 1234",
                    log_fh,
                )))

                results.append(("Decoherence (coupled)", run_step(
                    "Decoherence baseline (coupled, ξ=0.8)",
                    f"{py} -m quantum.cli --mode SSE_correlated "
                    f"--xi 0.8 --hiv_phase acute --N_r 36 --N_z 36 "
                    f"--dt 0.01 --time_steps 120 --frames_to_save 12 "
                    f"--rng_seed 1234",
                    log_fh,
                )))

        # ============================================================
        # SUMMARY
        # ============================================================
        log("", log_fh)
        log("=" * 60, log_fh)
        log("REPRODUCTION SUMMARY", log_fh)
        log("=" * 60, log_fh)

        passed = sum(1 for _, ok in results if ok)
        total = len(results)

        for name, ok in results:
            status = "PASS" if ok else "FAIL"
            log(f"  [{status}] {name}", log_fh)

        log("", log_fh)
        log(f"Result: {passed}/{total} steps passed", log_fh)
        log(f"Log:    {LOG_FILE}", log_fh)

        if passed == total:
            log("All manuscript results reproduced successfully.", log_fh)
        else:
            log("Some steps failed. Check log for details.", log_fh)
            sys.exit(1)


if __name__ == "__main__":
    main()
