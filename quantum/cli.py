import argparse
import json
import sys
from typing import Optional


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Run Microtubule Quantum Simulator from terminal (SSE-enabled)."
    )
    # Modes and phase
    p.add_argument('--dephasing_model', '--mode', dest='dephasing_model',
                   default='SSE_local', choices=['none', 'SSE_local', 'SSE_correlated'],
                   help='Dephasing model to use')
    p.add_argument('--hiv_phase', default='none', choices=['none', 'art_controlled', 'chronic', 'acute'],
                   help='HIV/inflammation phase controlling Γ scaling')

    # Grid and time
    p.add_argument('--N_r', type=int, default=36, help='Radial grid points')
    p.add_argument('--N_z', type=int, default=36, help='Axial grid points')
    p.add_argument('--dt', type=float, default=0.01, help='Time step')
    p.add_argument('--time_steps', type=int, default=120, help='Number of time steps')
    p.add_argument('--frames_to_save', type=int, default=12, help='Frames to save')

    # SSE correlation and RNG
    p.add_argument('--corr_length_xi', '--xi', dest='corr_length_xi', type=float, default=0.0,
                   help='Correlation length for SSE_correlated (grid units)')
    p.add_argument('--rng_seed', type=int, default=None, help='RNG seed for reproducibility')

    # Environment / Γ map params
    p.add_argument('--Gamma_0', type=float, default=0.05, help='Baseline decoherence rate')
    p.add_argument('--alpha_c', type=float, default=0.1, help='Cytokine coupling for Γ')
    p.add_argument('--gamma_scale_alpha', type=float, default=1.0, help='Global multiplier on cytokine coupling')
    p.add_argument('--temperature_K', type=float, default=310.0, help='Temperature (K)')
    p.add_argument('--ionic_strength_M', type=float, default=0.15, help='Ionic strength (M)')
    p.add_argument('--dielectric_rel', type=float, default=80.0, help='Relative dielectric constant')

    # Stability guard
    p.add_argument('--dt_gamma_guard_max', type=float, default=0.2,
                   help='Guard threshold for max(Γ)*dt; Γ map is clipped if exceeded')

    # Output dirs (optional overrides)
    p.add_argument('--output_dir', type=str, default=None, help='Output directory for results')
    p.add_argument('--figures_dir', type=str, default=None, help='Figures output directory')
    p.add_argument('--data_dir', type=str, default=None, help='Data output directory')

    return p.parse_args(argv)


def run_from_args(args) -> Optional[str]:
    # Local import to avoid heavy import on help
    from Legacy.microtubule_quantum_coherence_full_simulation import MicrotubuleQuantumSimulator

    cfg = {
        'dephasing_model': args.dephasing_model,
        'hiv_phase': args.hiv_phase,
        'N_r': args.N_r,
        'N_z': args.N_z,
        'dt': args.dt,
        'time_steps': args.time_steps,
        'frames_to_save': args.frames_to_save,
        'corr_length_xi': args.corr_length_xi,
        'rng_seed': args.rng_seed,
        'Gamma_0': args.Gamma_0,
        'alpha_c': args.alpha_c,
        'gamma_scale_alpha': args.gamma_scale_alpha,
        'temperature_K': args.temperature_K,
        'ionic_strength_M': args.ionic_strength_M,
        'dielectric_rel': args.dielectric_rel,
        'dt_gamma_guard_max': args.dt_gamma_guard_max,
    }
    if args.output_dir:
        cfg['output_dir'] = args.output_dir
    if args.figures_dir:
        cfg['figures_dir'] = args.figures_dir
    if args.data_dir:
        cfg['data_dir'] = args.data_dir

    sim = MicrotubuleQuantumSimulator(config=cfg)
    sim.run_simulation()
    result = sim.save_data()

    # save_data() may or may not return an identifier; handle both
    try:
        run_id = str(result) if result is not None else None
    except Exception:
        run_id = None

    # Print a concise terminal summary
    print("Simulation finished.")
    if run_id:
        print(f"Run ID: {run_id}")
    print("Key config:")
    print(json.dumps({
        'dephasing_model': args.dephasing_model,
        'hiv_phase': args.hiv_phase,
        'N_r': args.N_r,
        'N_z': args.N_z,
        'dt': args.dt,
        'time_steps': args.time_steps,
        'corr_length_xi': args.corr_length_xi,
        'rng_seed': args.rng_seed,
    }, indent=2))
    return run_id


def main(argv=None):
    args = parse_args(argv)
    run_from_args(args)


if __name__ == '__main__':
    sys.exit(main())
