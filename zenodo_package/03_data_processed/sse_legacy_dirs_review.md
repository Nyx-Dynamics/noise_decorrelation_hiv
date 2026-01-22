### Consolidated review: four_phase_tegmark_simulation_core_4, tegmarks_cat_source_code, tegmarks_cat_source_code_2, Tegmarks_Cat_source_code_3

Date: 2025-10-17

Purpose
- Provide a unified assessment of legacy/variant directories relevant to Tegmark-style decoherence and determine what to USE-AS-IS, ADAPT, or DEFER for the Option B (SSE) plan.

Reviewed directories and notable files
1) script/four_phase_tegmark_simulation_core 4
- Files: acute_hiv.py, art_controlled.py, chronic_hiv.py, study_volunteer.py, master_simulation_pipeline.py, grid.py, metrics.py, metrics 2.py, combined_tegmark_spiral_simulation 2.py, full_hiv_simulation 2.py
- Observations:
  - Phase-specific drivers (acute_hiv.py, art_controlled.py, chronic_hiv.py, study_volunteer.py) assemble simulations with cytokine-driven decoherence patterns. They contain plotting/report code and embed parameter presets for the four HIV phases.
  - grid.py provides grid constructors for cylindrical-like geometries and repeated numerical setups.
  - metrics.py and metrics 2.py appear to include dispersion/coherence measures; naming suggests duplication/variants.
  - Pipelines (master_simulation_pipeline.py, full_hiv_simulation 2.py) orchestrate runs across phases.
- Relevance to SSE:
  - ADAPT: Phase mapping patterns (how α and baseline Γ are scaled across HIV phases) are useful seeds for build_dephasing_map.
  - ADAPT: Any cylindrical-aware metric variants not already covered can inform coherence_metrics_sse.
  - DEFER: Monolithic simulation pipelines; we won’t integrate their execution, only re-use parameter logic.

2) script/tegmarks_cat_source_code
- Files: acute_hiv.py, art_controlled.py, chronic_hiv.py, code/… including four_phase_tegmark_simulation_core subfolder with acute/chronic/art_controlled, grid.py, combined_tegmark_spiral_simulation.py, full_hiv_simulation.py, and code/metrics.py.
- Observations:
  - Duplicates/variants of four-phase scripts live both at the root and under code/.
  - code/metrics.py provides a self-contained metrics implementation.
- Relevance to SSE:
  - USE-AS-IS (via shim): code/metrics.py for probability density/coherence/dispersion that respect cylindrical geometry.
  - ADAPT: The phase-specific α/Γ patterning (similar to above) can feed Γ_map construction.
  - DEFER: Standalone runner scripts duplicating simulation logic.

3) script/tegmarks_cat_source_code_2
- Files: acute_hiv.py; code/four_phase_tegmark_simulation_core/study_volunteer.py; other likely duplicates of phase scripts.
- Observations:
  - Another iteration of the same four-phase structure; the study_volunteer.py appears again here.
- Relevance to SSE:
  - ADAPT: Cross-check phase parameter presets to consolidate a single, clean mapping in build_dephasing_map.
  - DEFER: Execution pipelines and plotting; avoid importing these directly.

4) script/Tegmarks_Cat_source_code_3
- Files: wavefunction.py, metrics.py, hiv_quantum_stats.py
- Observations (detailed previously in project docs/sse_new_code_review.md):
  - wavefunction.py includes calculate_decoherence (phase-scaling) and cylindrical normalization ideas; evolution mixes deterministic damping (not SSE).
  - metrics.py has cylindrical-aware probability, coherence overlap, and dispersion/entropy.
  - hiv_quantum_stats.py is an incomplete reporting snippet.
- Relevance to SSE:
  - ADAPT: wavefunction.calculate_decoherence logic for Γ_map.
  - USE-AS-IS: metrics.py via a thin shim.
  - DEFER: hiv_quantum_stats.py.

Cross-directory synthesis and recommendations
- Γ_map (build_dephasing_map):
  - Consolidate phase scaling rules from: four_phase_tegmark_simulation_core 4 phase scripts, tegmarks_cat_source_code (root + code), and Tegmarks_Cat_source_code_3/wavefunction.py.
  - Keep the mapping in a single place with explicit inputs: hiv_phase, baseline Gamma_0, alpha, temperature_K, ionic_strength_M, dielectric_rel, and optional geometry-based modulation.
- Metrics for SSE:
  - Prefer script/Tegmarks_Cat_source_code_3/metrics.py as the primary source (it’s cohesive and cylindrical-aware). Use script/tegmarks_cat_source_code/code/metrics.py as a fallback.
- Geometry helpers:
  - Reference grid.py patterns for consistency but avoid importing them directly into the simulator; the main simulator already has grid setup.
- What to avoid:
  - Duplicative full-run pipelines and plotting/reporting scripts; they introduce side-effects and aren’t needed for the SSE layer.

Actionable mapping to Option B SSE
- USE-AS-IS (via shim):
  - calculate_probability_density, calculate_coherence, and dispersion metrics from Tegmarks_Cat_source_code_3/metrics.py, with a fallback to tegmarks_cat_source_code/code/metrics.py.
- ADAPT:
  - Phase→Γ scaling from wavefunction.calculate_decoherence (Tegmarks_Cat_source_code_3) and corroborating patterns in the four_phase* directories.
  - Encapsulate into quantum/open_systems.py: build_dephasing_map(r, z, hiv_phase, params).
- DEFER:
  - hiv_quantum_stats.py and runner/plot pipelines.

Risks and mitigations
- Path instability due to duplicated directories: mitigate with a small shim that tries imports in preferred order and logs which source was used.
- Inconsistent parameter names: normalize to a single config schema in the simulator when wiring SSE (hiv_phase, Gamma_0, alpha_c, temperature_K, ionic_strength_M, dielectric_rel, corr_length_xi, gamma_scale_alpha).
- Units: many legacy scripts are effectively dimensionless; document units in open_systems and annotate any scaling factors clearly.

Next steps
- Add quantum/shims.py to provide stable import wrappers for cylindrical-aware metrics, preferring Tegmarks_Cat_source_code_3. 
- Implement build_dephasing_map in quantum/open_systems.py (as per Project Roadmap), using consolidated phase mappings from these directories.
- When wiring SSE into the simulator, call the shimmed metrics for coherence logging and export, and save Γ_map built from the unified mapping.
