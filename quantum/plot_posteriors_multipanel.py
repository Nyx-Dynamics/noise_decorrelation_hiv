#!/usr/bin/env python3
"""
Make a compact multi‑panel posterior figure from a v3.6 run.

Panels (default):
- KDEs (or hist) for beta_xi, xi_acute, xi_chronic, xi_control
- Optional pair plot for (xi_acute, xi_chronic) with divergences overlay if available

Usage examples (from project root):

# Use final manuscript run (3:1:1, both eras)
python -m quantum.plot_posteriors_multipanel \
  --run-dir results/bayesian_v3_6/3_1_1/both/20251125T160846Z_ce3b0657 \
  --out figures/figures/posteriors_multipanel_3_1_1_both_v3_6.png

# Or pass a trace path directly
python -m quantum.plot_posteriors_multipanel \
  --trace results/bayesian_v3_6/3_1_1/both/20251125T160846Z_ce3b0657/trace_v3_6.nc
"""
from __future__ import annotations

import argparse
from pathlib import Path

import arviz as az
import numpy as np
import matplotlib.pyplot as plt


def find_vars(idata: az.InferenceData):
    post = idata.posterior
    # beta_xi
    beta = post.get('beta_xi')
    # Flexible names for xi across versions
    xi_a = post.get('xi_nm_acute', post.get('xi_acute_nm'))
    xi_c = post.get('xi_nm_chronic', post.get('xi_chronic_nm'))
    xi_ctrl = post.get('xi_nm_control', post.get('xi_healthy_nm'))
    return beta, xi_a, xi_c, xi_ctrl


def to_1d(x):
    if x is None:
        return None
    v = x.values if hasattr(x, 'values') else x
    return np.asarray(v).reshape(-1)


def hdi94(x):
    try:
        h = az.hdi(x, hdi_prob=0.94)
        return float(h[0]), float(h[1])
    except Exception:
        q = np.quantile(x, [0.03, 0.97])
        return float(q[0]), float(q[1])


def make_plot(idata: az.InferenceData, out_path: Path, title_suffix: str = ""):
    beta_raw, xi_a_raw, xi_c_raw, xi_ctrl_raw = find_vars(idata)
    beta = to_1d(beta_raw)
    xi_a = to_1d(xi_a_raw)
    xi_c = to_1d(xi_c_raw)
    xi_ctrl = to_1d(xi_ctrl_raw)

    # Figure layout: 2x2 KDEs + optional pair plot below
    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.1], hspace=0.35)

    # Panel A: beta_xi
    ax1 = fig.add_subplot(gs[0, 0])
    if beta is not None:
        az.plot_kde(beta, ax=ax1)
        m, sd = float(np.mean(beta)), float(np.std(beta))
        lo, hi = hdi94(beta)
        ax1.set_title(fr"$\\beta_\\xi$  (mean={m:.2f}, SD={sd:.2f}; 94% HDI [{lo:.2f}, {hi:.2f}])")
    else:
        ax1.text(0.5, 0.5, 'beta_xi missing', ha='center', va='center')
        ax1.set_title(r"$\beta_\xi$")

    # Panel B: xi_control
    ax2 = fig.add_subplot(gs[0, 1])
    if xi_ctrl is not None:
        az.plot_kde(xi_ctrl, ax=ax2)
        m, sd = float(np.mean(xi_ctrl)), float(np.std(xi_ctrl))
        lo, hi = hdi94(xi_ctrl)
        ax2.set_title(fr"$\\xi_{{control}}$ (mean={m:.3f} nm; 94% HDI [{lo:.3f}, {hi:.3f}])")
    else:
        ax2.text(0.5, 0.5, 'xi_control missing', ha='center', va='center')
        ax2.set_title(r"$\xi_{control}$")

    # Panel C: xi_acute
    ax3 = fig.add_subplot(gs[1, 0])
    if xi_a is not None:
        az.plot_kde(xi_a, ax=ax3)
        m, sd = float(np.mean(xi_a)), float(np.std(xi_a))
        lo, hi = hdi94(xi_a)
        ax3.set_title(fr"$\\xi_{{acute}}$ (mean={m:.3f} nm; 94% HDI [{lo:.3f}, {hi:.3f}])")
    else:
        ax3.text(0.5, 0.5, 'xi_acute missing', ha='center', va='center')
        ax3.set_title(r"$\xi_{acute}$")

    # Panel D: xi_chronic
    ax4 = fig.add_subplot(gs[1, 1])
    if xi_c is not None:
        az.plot_kde(xi_c, ax=ax4)
        m, sd = float(np.mean(xi_c)), float(np.std(xi_c))
        lo, hi = hdi94(xi_c)
        ax4.set_title(fr"$\\xi_{{chronic}}$ (mean={m:.3f} nm; 94% HDI [{lo:.3f}, {hi:.3f}])")
    else:
        ax4.text(0.5, 0.5, 'xi_chronic missing', ha='center', va='center')
        ax4.set_title(r"$\xi_{chronic}$")

    # Panel E: pair (xi_acute, xi_chronic) with divergences overlay if present
    ax5 = fig.add_subplot(gs[2, :])
    try:
        az.plot_pair(
            idata,
            var_names=[
                'xi_nm_acute' if 'xi_nm_acute' in idata.posterior else 'xi_acute_nm',
                'xi_nm_chronic' if 'xi_nm_chronic' in idata.posterior else 'xi_chronic_nm',
            ],
            kind='kde',
            divergences=True,
            ax=ax5,
        )
        ax5.set_title(r"$\xi_{acute}$ vs $\xi_{chronic}$ (pair plot; divergences overlay)")
    except Exception:
        # Fallback: scatter from arrays if available
        if xi_a is not None and xi_c is not None:
            ax5.scatter(xi_a, xi_c, s=5, alpha=0.2)
            ax5.set_xlabel(r"$\xi_{acute}$ (nm)")
            ax5.set_ylabel(r"$\xi_{chronic}$ (nm)")
            ax5.set_title(r"$\xi_{acute}$ vs $\xi_{chronic}$")
        else:
            ax5.text(0.5, 0.5, 'pair plot unavailable', ha='center', va='center')

    if title_suffix:
        fig.suptitle(f"Posterior summaries {title_suffix}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Wrote {out_path}")


def main():
    p = argparse.ArgumentParser(description='Multi‑panel posterior figure for v3.6 run')
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--trace', help='Path to trace_v3_6.nc')
    g.add_argument('--run-dir', help='Path to run directory containing trace_v3_6.nc')
    p.add_argument('--out', default='figures/figures/posteriors_multipanel_3_1_1_both_v3_6.png', help='Output PNG path')
    args = p.parse_args()

    if args.trace:
        trace_path = Path(args.trace)
        title_suffix = f"\n{trace_path.parent.name}"
    else:
        rd = Path(args.run_dir)
        trace_path = rd / 'trace_v3_6.nc'
        title_suffix = f"\n{rd.name}"

    idata = az.from_netcdf(trace_path)
    make_plot(idata, Path(args.out), title_suffix=title_suffix)


if __name__ == '__main__':
    main()
