import argparse
import json
import os
from typing import List, Dict, Any

import matplotlib.pyplot as plt
import numpy as np

"""
Phase sweep and multi-trajectory validation for SSE (Roadmap item 6)

- Runs multiple SSE trajectories per HIV phase and aggregates coherence decay statistics.
- Computes exponential-decay fits per trajectory and reports half-life distributions.
- Optionally annotates plots with a simple Tegmark-style Γ_Tegmark estimate.

Usage examples:
  python Extra/sse_phase_sweep.py --K 8 --phases none art_controlled chronic acute \
      --mode SSE_local --N_r 36 --N_z 36 --dt 0.01 --time_steps 120 --xi 0.8

Outputs (saved under data_dir and figures_dir via simulator plus a phase-sweep figure):
  - {run}_phase_sweep_summary.json (aggregate stats per phase)
  - {run}_phase_sweep.csv (optional per-trajectory metrics)
  - {run}_phase_sweep.png (mean coherence with 10–90% bands per phase)

Notes:
- Uses simulator’s internal saving and also returns the last base_filename per run; this script writes an
  additional aggregate JSON and figure to project figures_dir and Desktop mirror.
- To keep changes minimal, this helper instantiates the simulator repeatedly with different rng_seed.
"""

k_B = 1.380649e-23
hbar = 1.054571817e-34
NA = 6.02214076e23
e = 1.602176634e-19

def debye_length_m(temperature_K: float, ionic_strength_M: float, dielectric_rel: float) -> float:
    eps0 = 8.8541878128e-12
    I = max(1e-6, ionic_strength_M)
    I_m3 = 1000.0 * I
    return float(np.sqrt((dielectric_rel * eps0 * k_B * temperature_K) / (2.0 * (NA * (e ** 2)) * I_m3)))


def tegmark_gamma_estimate(delta_x: float,
                           temperature_K: float,
                           ionic_strength_M: float,
                           dielectric_rel: float,
                           scale_A: float = 1e-3) -> float:
    lamD = debye_length_m(temperature_K, ionic_strength_M, dielectric_rel)
    ratio = delta_x / max(lamD, 1e-12)
    return float(scale_A * (k_B * temperature_K / hbar) * (ratio ** 2))


def fit_exponential_decay(times, coherences):
    """Robust exponential decay fit C(t) ~ C0*exp(-gamma*t).

    Aligns times and coherences to the same length, filters non-finite/non-positive
    values, and fits on the monotonic (cumulative-min) segment.
    Returns gamma (>=0) or None if insufficient data.
    """
    # Convert to arrays
    t = np.asarray(times)
    c = np.asarray(coherences)

    # Align lengths defensively
    n = min(t.size, c.size)
    if n < 5:
        return None
    t = t[:n]
    c = c[:n]

    # Filter to valid values
    mask = np.isfinite(c) & (c > 0)
    if not np.any(mask):
        return None
    t = t[mask]
    c = c[mask]
    if t.size < 5:
        return None

    # Enforce monotonic non-increasing segment to avoid early transients
    c_mon = np.minimum.accumulate(c)

    try:
        y = np.log(c_mon)
        A = np.vstack([t, np.ones_like(t)]).T
        sol, *_ = np.linalg.lstsq(A, y, rcond=None)
        slope = sol[0]
        gamma_fit = -slope
        return float(max(0.0, gamma_fit))
    except Exception:
        return None


def run_trajectory(phase: str, mode: str, xi: float, base_seed: int, traj_idx: int,
                   N_r: int, N_z: int, dt: float, time_steps: int, frames_to_save: int,
                   Gamma_0: float, alpha_c: float,
                   temperature_K: float, ionic_strength_M: float, dielectric_rel: float) -> Dict[str, Any]:
    from Legacy.microtubule_quantum_coherence_full_simulation import MicrotubuleQuantumSimulator

    seed = int(base_seed + traj_idx)
    cfg = {
        'N_r': N_r,
        'N_z': N_z,
        'dt': dt,
        'time_steps': time_steps,
        'frames_to_save': frames_to_save,
        'hiv_phase': phase,
        'dephasing_model': mode,
        'corr_length_xi': xi,
        'Gamma_0': Gamma_0,
        'alpha_c': alpha_c,
        'temperature_K': temperature_K,
        'ionic_strength_M': ionic_strength_M,
        'dielectric_rel': dielectric_rel,
        'rng_seed': seed,
    }
    sim = MicrotubuleQuantumSimulator(config=cfg)
    sim.run_simulation()
    base = sim.save_data()
    times = np.array(sim.simulation_timestamps, dtype=float)
    # Prefer SSE coherence if present
    coh = np.array(sim.sse_coherence_reg if sim.sse_coherence_reg else sim.coherence_measure_reg, dtype=float)
    gamma_fit = fit_exponential_decay(times, coh)
    t_half = float(np.log(2.0) / gamma_fit) if (gamma_fit is not None and gamma_fit > 0) else None
    return {
        'phase': phase,
        'seed': seed,
        'timestamps': times.tolist(),
        'coherence': coh.tolist(),
        'gamma_fit': gamma_fit,
        't_half': t_half,
        'run_id': base,
    }


def aggregate_phase(results: List[Dict[str, Any]]):
    # Align by the minimal length across both timestamps and coherence arrays per trajectory
    if not results:
        return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float), {
            'gamma_mean': None, 'gamma_std': None, 't_half_mean': None, 't_half_std': None, 'n_trajs': 0,
        }
    min_len = min(min(len(r['timestamps']), len(r['coherence'])) for r in results)
    # Defensive: if extremely short, return minimal aligned arrays to avoid plotting errors
    if min_len < 2:
        times = np.array(results[0]['timestamps'][:min_len], dtype=float)
        mean = np.array(results[0]['coherence'][:min_len], dtype=float)
        p10 = mean.copy()
        p90 = mean.copy()
        gammas = np.array([r['gamma_fit'] for r in results if r.get('gamma_fit') is not None], dtype=float)
        thalf = np.array([r['t_half'] for r in results if r.get('t_half') is not None], dtype=float)
        stats = {
            'gamma_mean': float(np.nanmean(gammas)) if gammas.size else None,
            'gamma_std': float(np.nanstd(gammas)) if gammas.size else None,
            't_half_mean': float(np.nanmean(thalf)) if thalf.size else None,
            't_half_std': float(np.nanstd(thalf)) if thalf.size else None,
            'n_trajs': int(len(results)),
        }
        return times, mean, p10, p90, stats
    # Build aligned matrices
    times = np.array(results[0]['timestamps'][:min_len], dtype=float)
    mat = np.vstack([np.array(r['coherence'][:min_len], dtype=float) for r in results])
    mean = np.nanmean(mat, axis=0)
    p10 = np.nanpercentile(mat, 10, axis=0)
    p90 = np.nanpercentile(mat, 90, axis=0)
    gammas = np.array([r['gamma_fit'] for r in results if r.get('gamma_fit') is not None], dtype=float)
    thalf = np.array([r['t_half'] for r in results if r.get('t_half') is not None], dtype=float)
    stats = {
        'gamma_mean': float(np.nanmean(gammas)) if gammas.size else None,
        'gamma_std': float(np.nanstd(gammas)) if gammas.size else None,
        't_half_mean': float(np.nanmean(thalf)) if thalf.size else None,
        't_half_std': float(np.nanstd(thalf)) if thalf.size else None,
        'n_trajs': int(len(results)),
    }
    return times, mean, p10, p90, stats


def make_plot(phases, aggregates, outpath_png: str, annotate: Dict[str, float] = None):
    plt.figure(figsize=(8, 5))
    colors = {
        'none': '#4CAF50',
        'art_controlled': '#2196F3',
        'chronic': '#FF9800',
        'acute': '#F44336',
    }
    for ph in phases:
        times, mean, p10, p90, _ = aggregates[ph]
        c = colors.get(ph, None)
        label = f"{ph}"
        plt.plot(times, mean, label=label, color=c)
        plt.fill_between(times, p10, p90, color=c, alpha=0.2)
    plt.xlabel('time (sim units)')
    plt.ylabel('coherence overlap (reg)')
    plt.title('SSE coherence: mean and 10–90% bands across HIV phases')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath_png, dpi=150)
    try:
        desktop_png = os.path.expanduser(os.path.join('~/Desktop/microtubule_simulation/figures', os.path.basename(outpath_png)))
        plt.savefig(desktop_png, dpi=150)
    except Exception:
        pass
    plt.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--K', type=int, default=8, help='Number of trajectories per phase')
    p.add_argument('--phases', nargs='+', default=['none', 'art_controlled', 'chronic', 'acute'])
    p.add_argument('--mode', type=str, default='SSE_local', choices=['SSE_local', 'SSE_correlated'])
    p.add_argument('--xi', type=float, default=0.0, help='Correlation length for SSE_correlated')
    p.add_argument('--N_r', type=int, default=36)
    p.add_argument('--N_z', type=int, default=36)
    p.add_argument('--dt', type=float, default=0.01)
    p.add_argument('--time_steps', type=int, default=120)
    p.add_argument('--frames_to_save', type=int, default=12)
    p.add_argument('--Gamma_0', type=float, default=0.05)
    p.add_argument('--alpha_c', type=float, default=0.1)
    p.add_argument('--temperature_K', type=float, default=310.0)
    p.add_argument('--ionic_strength_M', type=float, default=0.15)
    p.add_argument('--dielectric_rel', type=float, default=80.0)
    p.add_argument('--delta_x', type=float, default=1e-9, help='For Tegmark Γ estimate annotation')
    p.add_argument('--scale_A', type=float, default=1e-3)
    p.add_argument('--base_seed', type=int, default=1000)
    args = p.parse_args()

    aggregates = {}
    per_phase_results: Dict[str, List[Dict[str, Any]]] = {}

    # Run trajectories
    for phase in args.phases:
        results = []
        for k in range(args.K):
            r = run_trajectory(
                phase=phase,
                mode=args.mode,
                xi=args.xi,
                base_seed=args.base_seed,
                traj_idx=k,
                N_r=args.N_r,
                N_z=args.N_z,
                dt=args.dt,
                time_steps=args.time_steps,
                frames_to_save=args.frames_to_save,
                Gamma_0=args.Gamma_0,
                alpha_c=args.alpha_c,
                temperature_K=args.temperature_K,
                ionic_strength_M=args.ionic_strength_M,
                dielectric_rel=args.dielectric_rel,
            )
            results.append(r)
        per_phase_results[phase] = results
        aggregates[phase] = aggregate_phase(results)

    # Prepare outputs
    # Aggregate JSON
    summary = {
        'config': vars(args),
        'phases': args.phases,
        'stats_per_phase': {},
        'tegmark_gamma': tegmark_gamma_estimate(args.delta_x, args.temperature_K, args.ionic_strength_M, args.dielectric_rel, args.scale_A),
    }
    for ph in args.phases:
        times, mean, p10, p90, stats = aggregates[ph]
        summary['stats_per_phase'][ph] = {
            'times': times.tolist(),
            'mean': mean.tolist(),
            'p10': p10.tolist(),
            'p90': p90.tolist(),
            **stats,
        }

    # Write JSON next to project data_dir
    data_dir = 'datafiles'
    os.makedirs(data_dir, exist_ok=True)
    sweep_json = os.path.join(data_dir, 'sse_phase_sweep_summary.json')
    with open(sweep_json, 'w') as f:
        json.dump(summary, f, indent=2)
    # Mirror to Desktop
    try:
        desktop_json = os.path.expanduser('~/Desktop/microtubule_simulation/datafiles/sse_phase_sweep_summary.json')
        os.makedirs(os.path.dirname(desktop_json), exist_ok=True)
        with open(desktop_json, 'w') as f:
            json.dump(summary, f, indent=2)
    except Exception:
        pass

    # Plot with variance bands
    figures_dir = 'figures'
    os.makedirs(figures_dir, exist_ok=True)
    out_png = os.path.join(figures_dir, 'sse_phase_sweep.png')
    make_plot(args.phases, aggregates, out_png)

    print('Phase sweep complete.')
    for ph in args.phases:
        _, _, _, _, stats = aggregates[ph]
        print(f"  {ph}: n={stats['n_trajs']}, t1/2 mean±std = {stats['t_half_mean']:.3g} ± {stats['t_half_std']:.3g} (sim units)" if stats['t_half_mean'] is not None else f"  {ph}: insufficient data for half-life")
    print('Saved:', sweep_json, 'and', out_png)


if __name__ == '__main__':
    main()
