import argparse
import json
import os
from typing import List, Dict, Any, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

"""
Monte Carlo analytics helper for SSE-enabled simulator.

Capabilities (subcommands):
- run: Execute an ensemble of simulations with sampled parameters; save per-trial results and an aggregate JSON.
- summarize: Aggregate existing per-trial JSON/CSV into a single summary JSON.
- viz: Plot 95% confidence bands and distribution of fitted decay rates/half-lives from the aggregate JSON.

Notes:
- Reuses fit_exponential_decay from Extra/sse_phase_sweep.py to avoid code duplication differences.
- Runs simulator per trial with seeded RNG based on a base_seed and trial index for reproducibility.
- Outputs are written under datafiles/ and figures/ and mirrored to ~/Desktop/microtubule_simulation/ when possible.

Usage examples:
  # Run a small MC ensemble (N=16) with local SSE
  PYTHONPATH=. python Extra/sse_mc_analytics.py run --N 16 --mode SSE_local \
      --N_r 36 --N_z 36 --dt 0.01 --time_steps 120 --frames_to_save 12 \
      --Gamma0_min 0.03 --Gamma0_max 0.07 --alpha_min 0.08 --alpha_max 0.12 \
      --phases none art_controlled chronic acute --base_seed 2000

  # Visualize aggregated results
  PYTHONPATH=. python Extra/sse_mc_analytics.py viz --summary datafiles/sse_mc_summary.json
"""


# ---- Utilities reused (import lightweight functions) ----

def fit_exponential_decay(times, coherences) -> Optional[float]:
    # Local import to avoid circulars
    import numpy as _np
    t = _np.asarray(times)
    c = _np.asarray(coherences)
    n = min(t.size, c.size)
    if n < 5:
        return None
    t = t[:n]
    c = c[:n]
    mask = _np.isfinite(c) & (c > 0)
    if not _np.any(mask):
        return None
    t = t[mask]
    c = c[mask]
    if t.size < 5:
        return None
    c_mon = _np.minimum.accumulate(c)
    try:
        y = _np.log(c_mon)
        A = _np.vstack([t, _np.ones_like(t)]).T
        sol, *_ = _np.linalg.lstsq(A, y, rcond=None)
        slope = sol[0]
        gamma_fit = -slope
        return float(max(0.0, gamma_fit))
    except Exception:
        return None


# ---- Sampling ----

def sample_params(N: int,
                  phases: List[str],
                  Gamma0_range: Tuple[float, float],
                  alpha_range: Tuple[float, float],
                  xi_range: Tuple[float, float],
                  mode: str,
                  rng: np.random.Generator) -> List[Dict[str, Any]]:
    samples: List[Dict[str, Any]] = []
    for i in range(N):
        phase = phases[i % len(phases)]
        g0 = rng.uniform(Gamma0_range[0], Gamma0_range[1])
        a = rng.uniform(alpha_range[0], alpha_range[1])
        if 'correlated' in mode.lower():
            xi = rng.uniform(xi_range[0], xi_range[1])
        else:
            xi = 0.0
        samples.append({'hiv_phase': phase, 'Gamma_0': float(g0), 'alpha_c': float(a), 'xi': float(xi)})
    return samples


# ---- Core execution ----

def run_trial(idx: int,
              cfg_base: Dict[str, Any],
              sample: Dict[str, Any],
              base_seed: int) -> Dict[str, Any]:
    from Legacy.microtubule_quantum_coherence_full_simulation import MicrotubuleQuantumSimulator

    cfg = dict(cfg_base)
    cfg.update({
        'hiv_phase': sample['hiv_phase'],
        'Gamma_0': sample['Gamma_0'],
        'alpha_c': sample['alpha_c'],
        'corr_length_xi': sample.get('xi', 0.0),
        'rng_seed': int(base_seed + idx),
    })
    sim = MicrotubuleQuantumSimulator(config=cfg)
    sim.run_simulation()
    run_id = sim.save_data()

    times = np.array(sim.simulation_timestamps, dtype=float)
    coh = np.array(sim.sse_coherence_reg if sim.sse_coherence_reg else sim.coherence_measure_reg, dtype=float)
    gamma_fit = fit_exponential_decay(times, coh)
    t_half = float(np.log(2.0) / gamma_fit) if (gamma_fit is not None and gamma_fit > 0) else None

    return {
        'idx': idx,
        'config': cfg,
        'sample': sample,
        'timestamps': times.tolist(),
        'coherence_reg': coh.tolist(),
        'gamma_fit': gamma_fit,
        't_half': t_half,
        'run_id': run_id,
    }


def aggregate_trials(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    gammas = np.array([r['gamma_fit'] for r in results if r.get('gamma_fit') is not None], dtype=float)
    thalf = np.array([r['t_half'] for r in results if r.get('t_half') is not None], dtype=float)
    summary = {
        'n_trials': int(len(results)),
        'gamma_mean': float(np.nanmean(gammas)) if gammas.size else None,
        'gamma_std': float(np.nanstd(gammas)) if gammas.size else None,
        'gamma_ci95': (
            float(np.nanpercentile(gammas, 2.5)) if gammas.size else None,
            float(np.nanpercentile(gammas, 97.5)) if gammas.size else None,
        ),
        't_half_mean': float(np.nanmean(thalf)) if thalf.size else None,
        't_half_std': float(np.nanstd(thalf)) if thalf.size else None,
        't_half_ci95': (
            float(np.nanpercentile(thalf, 2.5)) if thalf.size else None,
            float(np.nanpercentile(thalf, 97.5)) if thalf.size else None,
        ),
    }
    return summary


# ---- Visualization ----

def viz_summary(mc_json: str, out_png: Optional[str] = None):
    with open(mc_json, 'r') as f:
        data = json.load(f)
    gammas = np.array([r['gamma_fit'] for r in data.get('results', []) if r.get('gamma_fit') is not None], dtype=float)
    thalf = np.array([r['t_half'] for r in data.get('results', []) if r.get('t_half') is not None], dtype=float)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    if gammas.size:
        ax[0].hist(gammas, bins=20, color='#1565C0', alpha=0.8)
        ax[0].set_title('γ_fit distribution')
        ax[0].set_xlabel('γ_fit (sim units)')
        ax[0].set_ylabel('count')
    if thalf.size:
        ax[1].hist(thalf, bins=20, color='#EF6C00', alpha=0.8)
        ax[1].set_title('Half-life distribution')
        ax[1].set_xlabel('t1/2 (sim units)')
        ax[1].set_ylabel('count')
    plt.tight_layout()

    if out_png is None:
        out_png = os.path.join('figures', 'sse_mc_distributions.png')
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=150)
    try:
        desk_dir = os.path.expanduser('~/Desktop/microtubule_simulation/figures')
        os.makedirs(desk_dir, exist_ok=True)
        import shutil
        shutil.copyfile(out_png, os.path.join(desk_dir, os.path.basename(out_png)))
    except Exception:
        pass
    plt.close()
    print('Saved MC viz to', out_png)


# ---- CLI ----

def main(argv=None):
    p = argparse.ArgumentParser(description='Monte Carlo analytics for SSE simulator')
    sub = p.add_subparsers(dest='cmd', required=True)

    # run
    pr = sub.add_parser('run', help='Run a Monte Carlo ensemble')
    pr.add_argument('--N', type=int, default=16, help='Number of trials')
    pr.add_argument('--mode', type=str, default='SSE_local', choices=['SSE_local', 'SSE_correlated'], help='Dephasing model')
    pr.add_argument('--phases', nargs='+', default=['none', 'art_controlled', 'chronic', 'acute'])
    pr.add_argument('--Gamma0_min', type=float, default=0.03)
    pr.add_argument('--Gamma0_max', type=float, default=0.07)
    pr.add_argument('--alpha_min', type=float, default=0.08)
    pr.add_argument('--alpha_max', type=float, default=0.12)
    pr.add_argument('--xi_min', type=float, default=0.4)
    pr.add_argument('--xi_max', type=float, default=1.2)
    pr.add_argument('--N_r', type=int, default=36)
    pr.add_argument('--N_z', type=int, default=36)
    pr.add_argument('--dt', type=float, default=0.01)
    pr.add_argument('--time_steps', type=int, default=120)
    pr.add_argument('--frames_to_save', type=int, default=12)
    pr.add_argument('--base_seed', type=int, default=3000)
    pr.add_argument('--out', type=str, default='datafiles/sse_mc_summary.json')

    # summarize
    ps = sub.add_parser('summarize', help='Summarize existing per-trial JSON files into an aggregate JSON')
    ps.add_argument('--glob', type=str, default='datafiles/*_summary.json', help='Glob for per-trial summary JSONs')
    ps.add_argument('--out', type=str, default='datafiles/sse_mc_summary.json')

    # viz
    pv = sub.add_parser('viz', help='Visualize aggregate MC results')
    pv.add_argument('--summary', type=str, default='datafiles/sse_mc_summary.json')
    pv.add_argument('--out', type=str, default=None)

    args = p.parse_args(argv)

    if args.cmd == 'run':
        rng = np.random.default_rng(args.base_seed)
        samples = sample_params(
            N=args.N,
            phases=args.phases,
            Gamma0_range=(args.Gamma0_min, args.Gamma0_max),
            alpha_range=(args.alpha_min, args.alpha_max),
            xi_range=(args.xi_min, args.xi_max),
            mode=args.mode,
            rng=rng,
        )
        cfg_base = {
            'dephasing_model': args.mode,
            'N_r': args.N_r,
            'N_z': args.N_z,
            'dt': args.dt,
            'time_steps': args.time_steps,
            'frames_to_save': args.frames_to_save,
        }
        results: List[Dict[str, Any]] = []
        for i, samp in enumerate(samples):
            print(f"[MC] Trial {i+1}/{args.N}: phase={samp['hiv_phase']} G0={samp['Gamma_0']:.4f} alpha={samp['alpha_c']:.4f} xi={samp['xi']:.3f}")
            r = run_trial(i, cfg_base, samp, base_seed=args.base_seed)
            results.append(r)
        agg = aggregate_trials(results)
        payload = {
            'config': {
                'mode': args.mode,
                'phases': args.phases,
                'N': args.N,
                'Gamma0_range': [args.Gamma0_min, args.Gamma0_max],
                'alpha_range': [args.alpha_min, args.alpha_max],
                'xi_range': [args.xi_min, args.xi_max],
                'grid_time': {'N_r': args.N_r, 'N_z': args.N_z, 'dt': args.dt, 'time_steps': args.time_steps, 'frames_to_save': args.frames_to_save},
                'base_seed': args.base_seed,
            },
            'summary': agg,
            'results': results,
        }
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, 'w') as f:
            json.dump(payload, f, indent=2)
        # Mirror to Desktop
        try:
            out2 = os.path.expanduser('~/Desktop/microtubule_simulation/datafiles/sse_mc_summary.json')
            os.makedirs(os.path.dirname(out2), exist_ok=True)
            with open(out2, 'w') as f:
                json.dump(payload, f, indent=2)
        except Exception:
            pass
        print('MC run complete. Aggregate written to', args.out)
        return

    if args.cmd == 'summarize':
        import glob
        files = glob.glob(args.glob)
        results: List[Dict[str, Any]] = []
        for fp in sorted(files):
            try:
                with open(fp, 'r') as f:
                    js = json.load(f)
                # Pull sse_coherence_reg and timestamps if present
                times = js.get('timestamps', [])
                coh = js.get('sse_coherence_reg', []) or js.get('coherence_reg', [])
                gamma = fit_exponential_decay(times, coh)
                th = float(np.log(2.0) / gamma) if (gamma is not None and gamma > 0) else None
                results.append({'src': fp, 'gamma_fit': gamma, 't_half': th})
            except Exception:
                continue
        agg = aggregate_trials(results)
        payload = {'from_glob': args.glob, 'summary': agg, 'results': results}
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, 'w') as f:
            json.dump(payload, f, indent=2)
        print('Summarized', len(results), 'files into', args.out)
        return

    if args.cmd == 'viz':
        viz_summary(args.summary, args.out)
        return


if __name__ == '__main__':
    main()
