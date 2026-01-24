import json
import os

"""
Lightweight smoke test for Monte Carlo analytics.
Runs a very small ensemble (N=4) and checks that an aggregate JSON is written
with finite means when possible, and that distributions can be plotted.

Usage:
  PYTHONPATH=. python Extra/sse_mc_smoke_test.py
"""


def run_smoke():
    from quantum.sse_mc_analytics import main as mc_main

    # Run a tiny ensemble
    out = 'datafiles/sse_mc_summary_smoke.json'
    mc_main(['run', '--N', '4', '--mode', 'SSE_local', '--N_r', '24', '--N_z', '24', '--dt', '0.01',
             '--time_steps', '60', '--frames_to_save', '6', '--out', out])

    assert os.path.exists(out), f"Aggregate JSON not found: {out}"
    with open(out, 'r') as f:
        data = json.load(f)
    summ = data.get('summary', {})
    # Means may be None for very noisy tiny sets; just ensure keys present
    assert 'gamma_mean' in summ and 't_half_mean' in summ, 'Missing keys in MC summary'

    # Try visualization (should save a PNG)
    from quantum.sse_mc_analytics import viz_summary
    fig_out = 'figures/sse_mc_distributions_smoke.png'
    viz_summary(out, fig_out)
    assert os.path.exists(fig_out), 'MC viz figure was not saved'

    print('MC smoke test passed: summary at', out, 'figure at', fig_out)


if __name__ == '__main__':
    run_smoke()
