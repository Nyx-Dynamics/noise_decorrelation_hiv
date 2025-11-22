### Option B (SSE) Quickstart

This Quickstart shows how to enable and use the SSE (stochastic Schrödinger equation) dephasing options in the Microtubule_Simulation, how to interpret outputs, and how to run the included validation helpers.

Prerequisites
- Python environment with numpy, matplotlib, and scipy installed.
- This repository checked out locally.

Enabling SSE in the simulator
- Local (uncorrelated) SSE:
  - config={"dephasing_model": "SSE_local"}
- Spatially correlated SSE:
  - config={"dephasing_model": "SSE_correlated", "corr_length_xi": 0.8}

Relevant configuration parameters
- hiv_phase: one of {"none","acute","art_controlled","chronic"} controlling cytokine scaling of Γ.
- Gamma_0, alpha_c: baseline and cytokine coupling for Γ_map construction.
- temperature_K, ionic_strength_M, dielectric_rel: environmental context (used for documentation/exports; analytical comparisons).
- gamma_scale_alpha: global multiplier on cytokine coupling in Γ_map.
- corr_length_xi: spatial correlation length (grid units) for SSE_correlated.
- rng_seed: optional numpy RNG seed for reproducibility of stochastic steps.
- dt_gamma_guard_max: stability threshold; if max(Γ)*dt exceeds this value, Γ_map is clipped and a warning is logged.

Running a short SSE simulation
- Example (Python):
  from quantum.microtubule_quantum_coherence_full_simulation import MicrotubuleQuantumSimulator
  sim = MicrotubuleQuantumSimulator(config={
      'dephasing_model': 'SSE_local',
      'hiv_phase': 'acute',
      'rng_seed': 42,
      'N_r': 40, 'N_z': 40, 'dt': 0.01, 'time_steps': 80, 'frames_to_save': 8,
  })
  sim.run_simulation()
  sim.save_data()

Outputs and where to find them
- Summary JSON in data_dir and Desktop mirror, includes:
  - timestamps, variance/coherence time-series (legacy), and SSE metrics when enabled:
    - sse_coherence_reg/fib, sse_entropy_reg/fib, sse_variance_reg/fib
  - flags: Gamma_map_sse_final_present, SSE_kernel_present, dt_gamma_guard_triggered
  - config (with rng_seed)
- Arrays NPZ: final |ψ|^2 for both grids, cytokine_field, and Γ_map_sse (when enabled). Kernel NPZ saved for SSE_correlated.
- CSV: {run}_sse_coherence.csv with columns [time, sse_coherence_reg, sse_coherence_fib] when SSE is enabled.
- Figure: {run}_overlay.png — overlay of final |ψ|^2 (reg) with Γ_map_sse.

Validation helpers
- Extra/sse_validation.py: runs SSE_local then SSE_correlated and prints basic diagnostics.
  - python Extra/sse_validation.py
- Extra/tegmark_compare.py: compares simulated coherence decay with an order-of-magnitude Tegmark-style estimate.
  - python Extra/tegmark_compare.py path/to/{run}_summary.json --delta_x 1e-9

Smoke test
- Extra/sse_smoke_test.py: runs a very short SSE_local simulation and checks that norms remain ~1 and metrics are finite.
  - python Extra/sse_smoke_test.py

Notes and limitations
- Γ_map units are currently relative to the simulator’s internal time/space scales; the comparison script provides order-of-magnitude guidance rather than precise physical calibration.
- The dt·Γ guard helps maintain numerical stability; if you see frequent guard warnings, reduce dt, Gamma_0, alpha_c, or gamma_scale_alpha.
- For reproducible stochastic behavior, set rng_seed.
