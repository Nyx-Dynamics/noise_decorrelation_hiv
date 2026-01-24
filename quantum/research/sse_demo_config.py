import json
from datetime import datetime

import numpy as np

"""
Small‑grid demo config and quick stability check (Roadmap item 8)

Usage:
  python Extra/sse_demo_config.py --run

- Prints a recommended config dict for quick tests.
- With --run, executes a brief simulation and reports norms and dt·Γ guard status.
"""

import argparse

def recommended_config():
    return {
        'N_r': 36,
        'N_z': 36,
        'dt': 0.01,
        'time_steps': 60,
        'frames_to_save': 6,
        'hiv_phase': 'art_controlled',
        'dephasing_model': 'SSE_local',
        'rng_seed': 1234,
        'dt_gamma_guard_max': 0.2,
    }


def run_demo(cfg):
    from Legacy.microtubule_quantum_coherence_full_simulation import MicrotubuleQuantumSimulator
    sim = MicrotubuleQuantumSimulator(config=cfg)
    sim.run_simulation()
    run_id = sim.save_data()
    norm_reg = float(np.sum(np.abs(sim.Psi_reg) ** 2 * sim.r[:, None]) * sim.dr * sim.dz)
    norm_fib = float(np.sum(np.abs(sim.Psi_fib) ** 2 * sim.r[:, None]) * sim.dr * sim.dz)
    guard = bool(getattr(sim, '_dt_gamma_guard_triggered', False))
    print('Demo run complete:', run_id)
    print(f'  Final norms (reg, fib): {norm_reg:.6f}, {norm_fib:.6f}')
    print(f'  dt·Γ guard triggered: {guard}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', action='store_true', help='Run a short demo simulation')
    args = p.parse_args()
    cfg = recommended_config()
    print('Recommended small‑grid SSE config:')
    print(json.dumps(cfg, indent=2))
    if args.run:
        print('Running demo at', datetime.now().isoformat())
        run_demo(cfg)


if __name__ == '__main__':
    main()
