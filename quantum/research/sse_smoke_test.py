import numpy as np

"""
Lightweight smoke test for SSE_local mode.
Run: python Extra/sse_smoke_test.py
"""

def run_smoke():
    from Legacy.microtubule_quantum_coherence_full_simulation import MicrotubuleQuantumSimulator
    cfg = {
        'N_r': 24,
        'N_z': 24,
        'dt': 0.01,
        'time_steps': 30,
        'frames_to_save': 5,
        'hiv_phase': 'acute',
        'dephasing_model': 'SSE_local',
        'rng_seed': 1234,
    }
    sim = MicrotubuleQuantumSimulator(config=cfg)
    sim.run_simulation()
    sim.save_data()

    # Norm checks
    norm_reg = float(np.sum(np.abs(sim.Psi_reg) ** 2 * sim.r[:, None]) * sim.dr * sim.dz)
    norm_fib = float(np.sum(np.abs(sim.Psi_fib) ** 2 * sim.r[:, None]) * sim.dr * sim.dz)
    assert 0.9 < norm_reg < 1.1, f"Norm (reg) out of bounds: {norm_reg}"
    assert 0.9 < norm_fib < 1.1, f"Norm (fib) out of bounds: {norm_fib}"

    # If SSE metrics are present, ensure they are finite
    if getattr(sim, 'sse_coherence_reg', None):
        assert np.isfinite(sim.sse_coherence_reg[-1]), "Invalid coherence value"

    print("SSE smoke test passed. Final norms:", norm_reg, norm_fib)


if __name__ == '__main__':
    run_smoke()
