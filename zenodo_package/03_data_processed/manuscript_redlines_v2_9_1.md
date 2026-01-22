# v2.9.1 Redlines (Abstract, Results, Discussion)

This file applies the agreed redlines with canonical mapping, tempered language, and final numbers from the primary run `20251125T160846Z_ce3b0657`.

Final primary run (locked)
- Run ID: `20251125T160846Z_ce3b0657` (log‑ξ parameterization, `β_ξ ≥ 0`, 0 divergences)
- Key numbers: `β_ξ = 1.929 ± 0.247` (94% HDI `[1.464, 2.386]`); `ξ_acute = 0.786` (94% HDI `[0.697, 0.875]`); `ξ_chronic = 0.826` (94% HDI `[0.735, 0.923]`); `P(ξ_acute < ξ_chronic) = 1.00`.
- Canonical mapping: `Π_ξ(ξ) = (ξ_ref/ξ)^{β_ξ}` with `β_ξ > 0`.
- Phrase throughout: “shorter `ξ` ⇒ larger `Π_ξ` ⇒ more protection; `|β_ξ| ≈ 2`.”

---

## Abstract (results sentence)

Replace the current strong/mechanistic phrasing with:

“Our hierarchical Bayesian model infers a shorter effective correlation length in acute relative to chronic infection with posterior probability `P(ξ_acute < ξ_chronic) = 1.00`. A noise‑dependent protection factor `Π_ξ(ξ) = (ξ_ref/ξ)^{β_ξ}` with `β_ξ = 1.93` [1.46, 2.39] is supported over simpler variants by information criteria and cross‑validation, and the model achieves sub‑MRS‑precision prediction error across cohorts.”

Optional tail sentence to temper mechanism:

“These results offer a quantitative framework consistent with preserved NAA during acute HIV and generate specific testable hypotheses about noise‑mediated neuroprotection, while not establishing a specific biophysical mechanism.”

---

## Results (primary numbers, mapping, and comparisons)

Insert or update the Results summary to read:

- Primary inference (v3.6; 3:1:1; both eras; final run `20251125T160846Z_ce3b0657`) yields `β_ξ = 1.93` (mean 1.929, SD 0.247; 94% HDI [1.464, 2.386]).
- Phase‑specific effective parameter estimates: `ξ_acute = 0.786` (94% HDI [0.697, 0.875]) and `ξ_chronic = 0.826` (94% HDI [0.735, 0.923]); `P(ξ_acute < ξ_chronic) = 1.00`.
- We adopt the canonical mapping `Π_ξ(ξ) = (ξ_ref/ξ)^{β_ξ}` with `β_ξ > 0`; under this parameterization, shorter `ξ` implies larger `Π_ξ` (greater protection).
- Model comparison (WAIC/LOO; Table X) favours the full noise‑dependent model over constant‑`ξ` and `β_ξ = 0` variants with improved posterior predictive accuracy.
- Diagnostics (Supplementary Table Sx): all reported parameters exhibit R̂ ≈ 1.00 and large effective sample sizes; no divergences were observed in the final run.

Replace any older values (e.g., ~0.91 for `P(ξ_acute < ξ_chronic)`) with the numbers above; move historical values to a short sensitivity note if retained.

### Identifiability and ablation (short statements to accompany tables)

- “Across ablation variants, WAIC/LOO strongly favours the full model with noise‑dependent protection (Table X), consistent with improved posterior predictive accuracy.”
- “Across ratio configurations, the inferred exponent remains near `|β_ξ| ≈ 2`, and the acute–chronic separation is maintained with high posterior probability (Table Y).”
- “Analyses restricted to three global phase means cannot distinguish these mechanisms; the evidence for phase‑specific noise‑dependent protection emerges only when combining individual and group‑level data across studies.”

---

## Discussion (explicit caveat; hypothesis framing)

Use the following language to temper mechanistic claims while preserving the proposed hypothesis:

- “We emphasise that `ξ` is an inferred latent parameter constrained by MRS data and model structure; it should be interpreted as an effective phenomenological measure of noise structure, not yet a directly measured physical correlation length in vivo.”
- “If enzyme catalysis involves significant tunnelling contributions, a `|β_ξ| ≈ 2` scaling is compatible with simple tunnelling models; however, our data do not directly measure tunnelling or enzyme structure, and alternative mechanisms remain plausible.”
- “Accordingly, we frame the quantum/microtubule component as a mechanistic hypothesis consistent with the inferred scaling, and we outline specific experiments to test it.”

---

## Tables (Main) — compact templates

These tables are intended for the Main text. Place extended versions (per‑config diagnostics, SEs, and per‑condition PPC metrics) in the Supplement.

### Table X. Model ablation (Main)

Columns: Variant | ΔWAIC vs Full (negative favours Full) | ELPD_diff ± SE (optional) | PPC error (%; optional)

Rows: Full (phase‑ξ; `β_ξ` free) | Constant‑`ξ` (`β_ξ` free) | No protection (`β_ξ = 0`).

### Table Y. Ratio configuration comparison (Main)

Columns: Configuration | `β_ξ` mean [94% HDI] | `P(ξ_acute < ξ_chronic)` | PPC error (%; optional)

Rows: 3:1:1 (both eras) [primary], and additional configs (pre‑modern; post‑modern; optional other schemes) as space allows.

---

## Supplement — extended tables and diagnostics

Include full WAIC/LOO outputs (ELPD, pWAIC, SEs), per‑condition PPC metrics, full diagnostics (R̂, ESS, BFMI, divergences), and any legacy content clearly marked as superseded.
