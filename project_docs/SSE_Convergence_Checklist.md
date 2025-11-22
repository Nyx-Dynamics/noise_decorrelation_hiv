### SSE Convergence Checklist and Stability Guide (Roadmap item 8)

Purpose
- Provide a lightweight, practical procedure to verify numerical stability and convergence for SSE runs.
- Define recommended small‑grid demo settings and acceptance criteria.

Key stability guard
- dt·Γ guard (implemented in simulator via dt_gamma_guard_max): aim for max(Γ)*dt < 0.2 by default.
- If triggered frequently:
  - Reduce dt, or
  - Reduce Gamma_0 / alpha_c / gamma_scale_alpha, or
  - Increase spatial/temporal resolution gradually and retest.

Recommended small‑grid demo config (fast, for smoke/convergence checks)
- N_r = 36, N_z = 36
- dt = 0.01, time_steps = 120, frames_to_save = 12
- hiv_phase = 'acute' (stress‑test) and 'art_controlled' (moderate)
- dephasing_model = 'SSE_local' (use 'SSE_correlated' with corr_length_xi ≈ 0.8 for correlation tests)
- rng_seed = fixed (e.g., 1234) for reproducibility

Convergence procedure
1) Time step refinement (fixed grid)
- Start with dt = 0.02, then halve to 0.01 and 0.005 keeping other params fixed.
- For each dt, record sse_coherence_reg(t) and compute γ_fit (Extra/tegmark_compare.py).
- Acceptance: relative change in γ_fit between successive dt refinements < 10% and dt·Γ guard not triggered or rarely.

2) Spatial refinement (fixed dt once time step converged)
- Increase N_r, N_z by ~50% (e.g., 36→54) keeping dt fixed.
- Compare γ_fit and late‑time coherence; monitor memory/time costs.
- Acceptance: relative change in γ_fit < 10%; no systematic drift in late‑time coherence.

3) Correlation sensitivity (optional)
- If using SSE_correlated, vary corr_length_xi across {0.4, 0.8, 1.2} (grid units).
- Expect broadened variance bands with larger ξ; verify stability and dt·Γ guard status.

Acceptance criteria (minimum)
- Norms remain within 0.9–1.1 throughout (final norms checked in smoke tests).
- γ_fit stable under dt and modest grid refinement (< 10% change).
- No persistent dt·Γ guard clipping at the final chosen settings.

Artifacts to save
- Summary JSONs, CSV of sse_coherence, overlay figures per run.
- Phase sweep outputs (Extra/sse_phase_sweep.py) showing mean and 10–90% bands per phase.

Tips
- Use rng_seed to make runs reproducible when comparing parameter changes.
- For speed, prefer SSE_local for convergence testing; add SSE_correlated later as a sensitivity check.
- If coherence shows non‑monotonic early transients, fit γ on the later monotonic segment (built into helpers).
