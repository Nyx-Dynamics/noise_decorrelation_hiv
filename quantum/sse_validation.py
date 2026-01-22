from datetime import datetime

import numpy as np


# Minimal validation helper for SSE_local vs SSE_correlated
# Runs two short simulations with small grids and prints basic metrics.


def run_short(mode: str, xi: float = 0.0, seed: int = 42):
    from Legacy.microtubule_quantum_coherence_full_simulation import MicrotubuleQuantumSimulator

    config = {
        'N_r': 40,
        'N_z': 40,
        'dt': 0.01,
        'time_steps': 80,
        'frames_to_save': 8,
        'hiv_phase': 'acute',
        'dephasing_model': mode,
        'corr_length_xi': xi,
        'gamma_scale_alpha': 1.0,
        'temperature_K': 310.0,
        'ionic_strength_M': 0.15,
        'dielectric_rel': 80.0,
    }

    sim = MicrotubuleQuantumSimulator(config=config)
    # Optional: set numpy RNG seed for reproducibility within this process
    np.random.seed(seed)

    results = sim.run_simulation()
    base = sim.save_data()

    # Quick scalar diagnostics
    final_prob_reg = float(np.sum(np.abs(sim.Psi_reg) ** 2 * sim.r[:, None]) * sim.dr * sim.dz)
    final_prob_fib = float(np.sum(np.abs(sim.Psi_fib) ** 2 * sim.r[:, None]) * sim.dr * sim.dz)
    coh_reg = results['coherence_reg'][-1] if results['coherence_reg'] else None
    coh_fib = results['coherence_fib'][-1] if results['coherence_fib'] else None

    print(f"Mode={mode} xi={xi} -> final_norms (reg,fib)=({final_prob_reg:.6f},{final_prob_fib:.6f})")
    print(f"  final coherence (reg,fib)=({coh_reg:.4f},{coh_fib:.4f})  run_id={base}")

    return {
        'mode': mode,
        'xi': xi,
        'final_norm_reg': final_prob_reg,
        'final_norm_fib': final_prob_fib,
        'final_coh_reg': coh_reg,
        'final_coh_fib': coh_fib,
        'run_id': base,
    }


def main():
    print("SSE validation helper starting at", datetime.now().isoformat())
    out_local = run_short('SSE_local', xi=0.0)
    out_corr = run_short('SSE_correlated', xi=0.8)

    # Simple comparative message
    print("\nComparison summary:")
    print(f"  SSE_local coherence(reg)={out_local['final_coh_reg']:.4f}")
    print(f"  SSE_corr  coherence(reg)={out_corr['final_coh_reg']:.4f}  (xi={out_corr['xi']})")


if __name__ == '__main__':
    main()
