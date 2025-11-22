### Project Summary — Microtubule_Simulation (as of 2025-10-17 18:47)

Purpose
- Simulate quantum wavefunction dynamics in microtubules with cylindrical geometry and explore environmental (Tegmark‑style) decoherence.
- Implement Option B: stochastic Schrödinger equation (SSE) dephasing, with cytokine/HIV‑phase coupling and optional spatial correlations.

Executive highlights
- Geometry: cylindrical (r, z) grids with regular and Fibonacci variants; proper cylindrical normalization in outputs/metrics.
- Open quantum system (SSE, Option B):
  - SSE_local (uncorrelated, per‑site pure‑dephasing noise).
  - SSE_correlated (Gaussian‑kernel spatial correlations, correlation length ξ).
- Data/visualization: Γ_map and kernel NPZs, SSE coherence CSV, overlay figures, and a visualization helper.
- Reproducibility and stability: rng_seed support; dt·Γ guard with optional clipping; convergence checklist and small‑grid demo.

Key milestones delivered
1) Open‑systems scaffolding
- quantum/open_systems.py: build_dephasing_map, sse_dephasing_step, sse_dephasing_step_correlated, build_gaussian_kernel, coherence_metrics_sse.
- quantum/shims.py: stable imports for cylindrical‑aware metrics from legacy dirs.

2) Simulator wiring and config
- quantum/microtubule_quantum_coherence_full_simulation.py:
  - Config keys: dephasing_model ∈ {"none","SSE_local","SSE_correlated"}, temperature_K, ionic_strength_M, dielectric_rel, gamma_scale_alpha, corr_length_xi, rng_seed, dt_gamma_guard_max.
  - Per‑step SSE application after Hamiltonian update; cylindrical renormalization; Γ_map rebuild as cytokines evolve.

3) CLI and Makefile entry points
- CLI: python -m quantum.cli (mode, hiv_phase, grid/time, corr_length_xi, rng_seed, Γ params, guard).
- Makefile targets: run-local, run-corr, validate, smoke, tegmark-compare, phase-sweep, demo-config, demo-run, viz-* helpers.

4) Validation and analysis helpers (Extra/)
- sse_validation.py — quick local vs correlated runs with diagnostics.
- sse_smoke_test.py — sanity checks on norms and metric finiteness.
- tegmark_compare.py — OOM analytical Γ_Tegmark vs simulated decay fit.
- sse_phase_sweep.py — multi‑trajectory phase sweep; mean and 10–90% variance bands; half‑life stats per phase.
- sse_demo_config.py — small‑grid demo + guard status report.
- sse_visualize.py — plot coherence time series, final overlays, phase bands, and kernels.

5) Documentation
- README.md — project overview, terminal quickstart, visualization usage.
- project docs/Option_B_SSE_Quickstart.md — enabling and using SSE modes; outputs.
- project docs/SSE_Convergence_Checklist.md — stability/convergence guidance.
- Project Roadmap — scope, plan, and multiple dated status updates (latest: 2025‑10‑17 18:24).
- Extra/venv_info.txt — end‑to‑end venv + execution and visualization guide.

Current capabilities (at a glance)
- Run SSE_local or SSE_correlated simulations; save Γ_map and kernel NPZs, SSE coherence CSV, summary JSON (config, flags, metrics), and overlay PNGs.
- Visualize results via helper targets: viz-latest-summary, viz-latest-coherence, viz-phase-sweep, viz-latest-kernel.
- Compare simulated decoherence with a Tegmark‑style analytical estimate (order‑of‑magnitude guidance).

How to run (short examples)
- Local SSE: python -m quantum.cli --mode SSE_local --hiv_phase acute --N_r 36 --N_z 36 --dt 0.01 --time_steps 120 --frames_to_save 12 --rng_seed 1234
- Correlated SSE (ξ=0.8): python -m quantum.cli --mode SSE_correlated --xi 0.8 --hiv_phase acute --N_r 36 --N_z 36 --dt 0.01 --time_steps 120 --frames_to_save 12 --rng_seed 1234
- Makefile: make run-local | make run-corr | make validate | make smoke | make tegmark-compare SUMMARY=... DELTA_X=1e-9

Primary outputs
- datafiles/:
  - {run}_summary.json — configs, timestamps, legacy and SSE metrics, flags (Gamma_map_sse_final_present, SSE_kernel_present, dt_gamma_guard_triggered).
  - {run}_sse_coherence.csv — time, sse_coherence_reg, sse_coherence_fib.
  - {run}_arrays.npz — final ψ arrays and grids; {run}_gamma_map_sse.npz — Γ_map; {run}_sse_kernel.npz — kernel (correlated runs only).
- figures/: {run}_overlay.png; phase sweep: sse_phase_sweep.png.
- Mirrored to ~/Desktop/microtubule_simulation/ when possible.

Stability and convergence
- dt·Γ guard: configurable via dt_gamma_guard_max (default 0.2), clips Γ_map and records a flag if exceeded.
- Convergence guide: project docs/SSE_Convergence_Checklist.md; small‑grid demo in Extra/sse_demo_config.py.

Roadmap status (abridged)
- Items 1–10: Completed (scaffolding, wiring, visualization/outputs, validation helpers, smoke test, documentation).
- 11: Ongoing (monitor stability guard across diverse configs; minor doc polish as needed).
- 12: Pending (final submission/closure after confirmation).

Repository map (SSE‑related)
- quantum/: open_systems.py, shims.py, microtubule_quantum_coherence_full_simulation.py, cli.py
- Extra/: sse_validation.py, sse_smoke_test.py, sse_phase_sweep.py, tegmark_compare.py, sse_demo_config.py, sse_visualize.py, venv_info.txt
- project docs/: Option_B_SSE_Quickstart.md, SSE_Convergence_Checklist.md, recommendations_tegmark_extension.md, sse_legacy_dirs_review.md, PROJECT_SUMMARY_2025-10-17.md (this file)
- Root: README.md, Project Roadmap, Makefile, requirements.txt

Next steps (suggested)
- Close Roadmap Item 12 after final confirmation; keep monitoring dt·Γ guard under larger grids and correlated runs.
- Optional: refine Γ_map physical calibration vs analytical estimates; add unit annotations for SI alignment if needed.

Reproducibility snapshot
- Python deps: requirements.txt (numpy, scipy, matplotlib). venv workflow documented in Extra/venv_info.txt.
- Determinism: set rng_seed in configs/CLI.

Citations
- M. Tegmark, Phys. Rev. E 61, 4194 (2000); Breuer & Petruccione (2002); Haken–Strobl–Reineker (1972).