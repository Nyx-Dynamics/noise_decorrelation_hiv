import argparse
import json
import os
from datetime import datetime
from typing import Dict, Any, Optional

import numpy as np

"""
Data interpretation helper for SSE outputs.

Generates concise, human-readable interpretations (Markdown) and a
machine-readable JSON summary from existing artifacts:
- Per-run summary JSON (+ optional SSE coherence CSV)
- Phase-sweep aggregate JSON (Extra/sse_phase_sweep.py)
- Monte Carlo aggregate JSON (Extra/sse_mc_analytics.py)

Outputs are saved to results/ and mirrored to
~/Desktop/microtubule_simulation/results when possible.
"""

k_B = 1.380649e-23
hbar = 1.054571817e-34
NA = 6.02214076e23
E_CHARGE = 1.602176634e-19


def ensure_dir(path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)


def mirror_to_desktop(path: str):
    try:
        desk_dir = os.path.expanduser('~/Desktop/microtubule_simulation/results')
        os.makedirs(desk_dir, exist_ok=True)
        import shutil
        shutil.copyfile(path, os.path.join(desk_dir, os.path.basename(path)))
    except Exception:
        pass


# --- Shared helpers ---

def fit_exponential_decay(times, coherences) -> Optional[float]:
    """Robust exponential fit C(t) ~ C0 * exp(-gamma * t). Returns gamma or None."""
    t = np.asarray(times, dtype=float)
    c = np.asarray(coherences, dtype=float)
    n = min(t.size, c.size)
    if n < 5:
        return None
    t = t[:n]
    c = c[:n]
    mask = np.isfinite(c) & (c > 0)
    if not np.any(mask):
        return None
    t = t[mask]
    c = c[mask]
    if t.size < 5:
        return None
    c_mon = np.minimum.accumulate(c)
    try:
        y = np.log(c_mon)
        A = np.vstack([t, np.ones_like(t)]).T
        sol, *_ = np.linalg.lstsq(A, y, rcond=None)
        slope = sol[0]
        gamma = float(max(0.0, -slope))
        return gamma
    except Exception:
        return None


def debye_length_m(temperature_K: float, ionic_strength_M: float, dielectric_rel: float) -> float:
    eps0 = 8.8541878128e-12
    I = max(1e-6, float(ionic_strength_M))
    I_m3 = 1000.0 * I
    return float(np.sqrt((dielectric_rel * eps0 * k_B * temperature_K) / (2.0 * (NA * (E_CHARGE ** 2)) * I_m3)))


def tegmark_gamma_estimate(delta_x: float,
                           temperature_K: float,
                           ionic_strength_M: float,
                           dielectric_rel: float,
                           scale_A: float = 1e-3) -> float:
    lamD = debye_length_m(temperature_K, ionic_strength_M, dielectric_rel)
    ratio = float(delta_x) / max(lamD, 1e-12)
    return float(scale_A * (k_B * temperature_K / hbar) * (ratio ** 2))


# --- Per-run interpretation ---

def interpret_run(summary_json: str,
                  delta_x: float = 1e-9,
                  scale_A: float = 1e-3,
                  out_prefix: Optional[str] = None) -> Dict[str, Any]:
    with open(summary_json, 'r') as f:
        summary = json.load(f)

    cfg = summary.get('config', {})
    times = np.array(summary.get('timestamps', []), dtype=float)
    sse_coh = summary.get('sse_coherence_reg') or []
    legacy_coh = summary.get('coherence_reg') or []
    coherences = sse_coh if sse_coh else legacy_coh

    gamma_fit = fit_exponential_decay(times, coherences)
    t_half = float(np.log(2.0) / gamma_fit) if (gamma_fit is not None and gamma_fit > 0) else None

    # Analytical Tegmark OOM estimate
    T = float(cfg.get('temperature_K', 310.0))
    I = float(cfg.get('ionic_strength_M', 0.15))
    epsr = float(cfg.get('dielectric_rel', 80.0))
    gamma_T = tegmark_gamma_estimate(delta_x, T, I, epsr, scale_A)

    guard = bool(summary.get('dt_gamma_guard_triggered', False))
    kernel_present = bool(summary.get('SSE_kernel_present', False))
    gamma_present = bool(summary.get('Gamma_map_sse_final_present', False))

    # Simple traffic-light interpretation for stability
    if guard:
        stability_flag = 'warn'  # dt*Gamma clipping occurred
    else:
        stability_flag = 'ok'

    # Relative comparison text
    compare_txt = None
    if gamma_fit is not None and gamma_fit > 0:
        ratio = gamma_fit / max(gamma_T, 1e-30)
        if ratio < 0.1:
            compare_txt = 'Simulated decay slower than Tegmark OOM (ratio < 0.1); coherence more robust.'
        elif ratio > 10:
            compare_txt = 'Simulated decay faster than Tegmark OOM (ratio > 10); strong dephasing regime.'
        else:
            compare_txt = 'Simulated decay within an order of magnitude of Tegmark estimate.'

    findings = {
        'run_summary_path': summary_json,
        'hiv_phase': cfg.get('hiv_phase', 'unknown'),
        'dephasing_model': cfg.get('dephasing_model', 'unknown'),
        'corr_length_xi': cfg.get('corr_length_xi', 0.0),
        'rng_seed': cfg.get('rng_seed'),
        'gamma_fit': gamma_fit,
        't_half': t_half,
        'gamma_Tegmark': gamma_T,
        'delta_x_used': delta_x,
        'stability_flag': stability_flag,
        'dt_gamma_guard_triggered': guard,
        'Gamma_map_present': gamma_present,
        'kernel_present': kernel_present,
        'comparison_text': compare_txt,
    }

    # Compose Markdown
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    lines = []
    lines.append(f"# Data Interpretation — Single Run\n")
    lines.append(f"Generated: {ts}\n")
    lines.append("")
    lines.append(f"- HIV phase: {findings['hiv_phase']}")
    lines.append(f"- Dephasing model: {findings['dephasing_model']} (xi={findings['corr_length_xi']})")
    if gamma_fit is not None:
        lines.append(f"- Simulated decay rate γ_fit ≈ {gamma_fit:.3e} [1/sim time]")
        if t_half is not None:
            lines.append(f"- Coherence half-life t1/2 ≈ {t_half:.3g} sim time units")
    else:
        lines.append("- γ_fit: insufficient data to fit (non-monotonic or too short)")
    lines.append(f"- Tegmark OOM Γ_T ≈ {gamma_T:.3e} s^-1 (Δx={delta_x:g} m, T={T} K, I={I} M, ε_r={epsr})")
    if compare_txt:
        lines.append(f"- Comparison: {compare_txt}")
    lines.append(f"- Stability: {'CLIPPED (guard triggered)' if guard else 'OK'}")
    lines.append(f"- Artifacts: Γ_map={gamma_present}, kernel={kernel_present}")

    # Write outputs
    base = os.path.basename(summary_json).replace('_summary.json', '')
    out_base = out_prefix if out_prefix else os.path.join('results', base)
    md_path = f"{out_base}_interpretation.md"
    json_path = f"{out_base}_interpretation.json"
    ensure_dir(md_path)
    ensure_dir(json_path)
    with open(md_path, 'w') as f:
        f.write("\n".join(lines) + "\n")
    with open(json_path, 'w') as f:
        json.dump(findings, f, indent=2)
    mirror_to_desktop(md_path)
    mirror_to_desktop(json_path)
    print('Wrote:', md_path)
    print('Wrote:', json_path)
    return findings


# --- Phase-sweep interpretation ---

def interpret_phase_sweep(summary_json: str, out_prefix: Optional[str] = None) -> Dict[str, Any]:
    with open(summary_json, 'r') as f:
        data = json.load(f)
    phases = data.get('phases') or list((data.get('stats_per_phase') or {}).keys())
    stats = data.get('stats_per_phase', {})
    if not phases or not stats:
        raise ValueError('Invalid phase-sweep summary JSON: missing phases/stats_per_phase')

    # Build findings per phase
    findings: Dict[str, Any] = {
        'source': summary_json,
        'per_phase': {},
    }
    lines = ["# Data Interpretation — Phase Sweep", ""]
    for ph in phases:
        entry = stats.get(ph)
        if not entry:
            continue
        gamma_mean = entry.get('gamma_mean')
        gamma_std = entry.get('gamma_std')
        t_half_mean = entry.get('t_half_mean')
        t_half_std = entry.get('t_half_std')
        n_trajs = entry.get('n_trajs')
        findings['per_phase'][ph] = {
            'gamma_mean': gamma_mean,
            'gamma_std': gamma_std,
            't_half_mean': t_half_mean,
            't_half_std': t_half_std,
            'n_trajs': n_trajs,
        }
        lines.append(f"## {ph}")
        if gamma_mean is not None:
            lines.append(f"- γ_fit mean±std: {gamma_mean:.3e} ± {gamma_std:.3e} (n={n_trajs})")
        else:
            lines.append("- Insufficient data for γ_fit")
        if t_half_mean is not None:
            lines.append(f"- t1/2 mean±std: {t_half_mean:.3g} ± {t_half_std:.3g} (sim units)")
        lines.append("")

    # Simple qualitative ranking by half-life if available
    ranked = [
        (ph, findings['per_phase'][ph].get('t_half_mean')) for ph in phases
        if findings['per_phase'][ph].get('t_half_mean') is not None
    ]
    ranked.sort(key=lambda x: x[1], reverse=True)
    if ranked:
        lines.append("### Relative robustness (by t1/2 mean, highest first)")
        for ph, th in ranked:
            lines.append(f"- {ph}: t1/2 ≈ {th:.3g}")

    base = os.path.splitext(os.path.basename(summary_json))[0].replace('_summary', '')
    out_base = out_prefix if out_prefix else os.path.join('results', f'{base}_phase')
    md_path = f"{out_base}_interpretation.md"
    json_path = f"{out_base}_interpretation.json"
    ensure_dir(md_path)
    ensure_dir(json_path)
    with open(md_path, 'w') as f:
        f.write("\n".join(lines) + "\n")
    with open(json_path, 'w') as f:
        json.dump(findings, f, indent=2)
    mirror_to_desktop(md_path)
    mirror_to_desktop(json_path)
    print('Wrote:', md_path)
    print('Wrote:', json_path)
    return findings


# --- Monte Carlo interpretation ---

def interpret_mc(summary_json: str, out_prefix: Optional[str] = None) -> Dict[str, Any]:
    with open(summary_json, 'r') as f:
        data = json.load(f)

    # Expect schema from Extra/sse_mc_analytics.py; be robust if keys differ.
    per_phase = data.get('per_phase') or data.get('stats_per_phase') or {}
    findings: Dict[str, Any] = {'source': summary_json, 'per_phase': {}}

    lines = ["# Data Interpretation — Monte Carlo Ensemble", ""]
    for ph, entry in per_phase.items():
        gamma_ci = entry.get('gamma_95ci') or entry.get('gamma_ci95')
        t_half_ci = entry.get('t_half_95ci') or entry.get('t_half_ci95')
        gamma_mean = entry.get('gamma_mean')
        t_half_mean = entry.get('t_half_mean')
        findings['per_phase'][ph] = {
            'gamma_mean': gamma_mean,
            'gamma_95ci': gamma_ci,
            't_half_mean': t_half_mean,
            't_half_95ci': t_half_ci,
        }
        lines.append(f"## {ph}")
        if gamma_mean is not None:
            lines.append(f"- γ_fit mean: {gamma_mean:.3e}")
        if gamma_ci:
            lo, hi = gamma_ci
            lines.append(f"- γ_fit 95% CI: [{lo:.3e}, {hi:.3e}]")
        if t_half_mean is not None:
            lines.append(f"- t1/2 mean: {t_half_mean:.3g} (sim units)")
        if t_half_ci:
            lo, hi = t_half_ci
            lines.append(f"- t1/2 95% CI: [{lo:.3g}, {hi:.3g}]")
        lines.append("")

    # Optional: quick sensitivity note if chronic >> none in gamma
    try:
        g_none = per_phase.get('none', {}).get('gamma_mean')
        g_acute = per_phase.get('acute', {}).get('gamma_mean')
        if g_none is not None and g_acute is not None:
            ratio = float(g_acute) / max(float(g_none), 1e-30)
            if ratio > 2.0:
                lines.append(f"Observation: acute shows >2x faster decay than none (γ ratio ≈ {ratio:.1f}).")
    except Exception:
        pass

    base = os.path.splitext(os.path.basename(summary_json))[0].replace('_summary', '')
    out_base = out_prefix if out_prefix else os.path.join('results', f'{base}_mc')
    md_path = f"{out_base}_interpretation.md"
    json_path = f"{out_base}_interpretation.json"
    ensure_dir(md_path)
    ensure_dir(json_path)
    with open(md_path, 'w') as f:
        f.write("\n".join(lines) + "\n")
    with open(json_path, 'w') as f:
        json.dump(findings, f, indent=2)
    mirror_to_desktop(md_path)
    mirror_to_desktop(json_path)
    print('Wrote:', md_path)
    print('Wrote:', json_path)
    return findings


# --- CLI ---

def main(argv=None):
    p = argparse.ArgumentParser(description='Data interpretation helper for SSE outputs')
    sub = p.add_subparsers(dest='cmd', required=True)

    p_run = sub.add_parser('run', help='Interpret a single run summary JSON')
    p_run.add_argument('--summary', required=True, help='Path to *_summary.json')
    p_run.add_argument('--delta_x', type=float, default=1e-9, help='Δx for Tegmark OOM estimate (meters)')
    p_run.add_argument('--scale_A', type=float, default=1e-3, help='Scale factor for Tegmark OOM estimate')
    p_run.add_argument('--out_prefix', type=str, default=None, help='Optional output prefix (path without extension)')

    p_phase = sub.add_parser('phase', help='Interpret phase-sweep summary JSON')
    p_phase.add_argument('--summary-json', default='datafiles/sse_phase_sweep_summary.json')
    p_phase.add_argument('--out_prefix', type=str, default=None)

    p_mc = sub.add_parser('mc', help='Interpret Monte Carlo aggregate summary JSON')
    p_mc.add_argument('--summary', default='datafiles/sse_mc_summary.json')
    p_mc.add_argument('--out_prefix', type=str, default=None)

    args = p.parse_args(argv)
    if args.cmd == 'run':
        interpret_run(args.summary, args.delta_x, args.scale_A, args.out_prefix)
    elif args.cmd == 'phase':
        interpret_phase_sweep(args.summary_json, args.out_prefix)
    elif args.cmd == 'mc':
        interpret_mc(args.summary, args.out_prefix)


if __name__ == '__main__':
    main()
