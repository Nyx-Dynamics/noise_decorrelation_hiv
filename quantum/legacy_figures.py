"""
Legacy-style figure generator for manuscript outputs.

Usage examples:
  python -m quantum.legacy_figures --model bayesian_v3_6 --ratio 3_1_1 --era both \
    --suffix 3_1_1_both_v3_6 --outdir figures/figures

  python -m quantum.legacy_figures --model enzyme_v3_ode --ratio 3_1_1 --era both \
    --suffix 3_1_1_both_ode --outdir figures/figures

- If --run-id is omitted, the latest run directory is selected automatically.
- Optional --mirror-root will copy the generated figures to an external folder (e.g., iCloud runs/*).
- This module intentionally focuses on figure generation only and does not change core model logic.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import arviz as az  # type: ignore
except Exception:  # pragma: no cover
    az = None  # type: ignore


def parse_args():
    p = argparse.ArgumentParser(description="Generate legacy-style figures for manuscript")
    p.add_argument("--model", required=True, choices=["bayesian_v3_6", "bayesian_v3", "enzyme_v3_ode"])
    p.add_argument("--ratio", required=True, choices=["3_1_1", "1_1_1", "1_2_1", "2_1"])
    p.add_argument("--era", required=True, choices=["pre_modern", "post_modern", "both"])
    p.add_argument("--run-id", default=None, help="Specific run_id to use; if omitted, picks latest run")
    p.add_argument("--suffix", default=None, help="Suffix used in output filenames (e.g., 3_1_1_both_v3_6)")
    p.add_argument("--outdir", default="figures/figures", help="Directory to write figures into")
    p.add_argument("--mirror-root", default=None, help="Optional folder to mirror generated figures to")
    return p.parse_args()


def _find_latest_run(base: Path) -> Optional[Path]:
    if not base.exists():
        return None
    runs: List[Path] = [p for p in base.iterdir() if p.is_dir()]
    if not runs:
        return None
    runs.sort()
    return runs[-1]


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _mirror_dir(src: Path, dst: Path) -> None:
    import shutil
    _ensure_dir(dst)
    # Copy only files from src into dst (flat); do not delete existing
    for f in src.glob("*.png"):
        shutil.copy2(str(f), str(dst / f.name))


def _fig_ppd_overlay(pp_df: pd.DataFrame, out: Path, title: str) -> None:
    phases = ["Control", "Acute", "Chronic"]
    obs = pp_df.groupby("Phase")["NAA_mean_observed"].mean()
    pp = pp_df.groupby("Phase")["NAA_mean_predicted"].mean() if "NAA_mean_predicted" in pp_df.columns else None
    x = np.arange(len(phases)); w = 0.35
    plt.figure(figsize=(6, 4))
    plt.bar(x - w/2, [obs.get(p, np.nan) for p in phases], width=w, label="Observed")
    if pp is not None:
        plt.bar(x + w/2, [pp.get(p, np.nan) for p in phases], width=w, label="Posterior Predicted")
    plt.xticks(x, phases)
    plt.ylabel("NAA/Cr")
    plt.title(title)
    plt.legend()
    plt.savefig(out, dpi=200, bbox_inches="tight"); plt.close()


def _fig_forest_hdi(post, out: Path, title: str) -> None:
    if az is None:
        return
    xi_baseline_nm = 0.66
    baseline_ratio = 2.2
    beta = post["beta_xi"]
    xi_ctrl = post["xi_nm_control"]; xi_acute = post["xi_nm_acute"]; xi_chron = post["xi_nm_chronic"]
    Pi_ctrl = (xi_ctrl / xi_baseline_nm) ** (-beta)
    Pi_acute = (xi_acute / xi_baseline_nm) ** (-beta)
    Pi_chron = (xi_chron / xi_baseline_nm) ** (-beta)
    means = [
        (baseline_ratio * Pi_ctrl).stack(draws=("chain", "draw")).values,
        (baseline_ratio * Pi_acute).stack(draws=("chain", "draw")).values,
        (baseline_ratio * Pi_chron).stack(draws=("chain", "draw")).values,
    ]
    az.plot_forest(means, model_names=["Control", "Acute", "Chronic"], hdi_prob=0.95)
    plt.title(title)
    plt.savefig(out, dpi=200, bbox_inches="tight"); plt.close()


def _fig_param_kdes(post, outdir: Path, suffix: str) -> None:
    if az is None:
        return
    items = [
        (post['beta_xi'], 'β_ξ', f'posterior_beta_xi_kde_{suffix}.png'),
        (post['xi_nm_control'], 'ξ_control (nm)', f'posterior_xi_control_kde_{suffix}.png'),
        (post['xi_nm_acute'],   'ξ_acute (nm)',   f'posterior_xi_acute_kde_{suffix}.png'),
        (post['xi_nm_chronic'], 'ξ_chronic (nm)', f'posterior_xi_chronic_kde_{suffix}.png'),
    ]
    for var, label, name in items:
        arr = var.stack(draws=("chain", "draw")).values
        az.plot_kde(arr)
        plt.title(label); plt.xlabel(label)
        plt.savefig(outdir / name, dpi=200, bbox_inches="tight"); plt.close()


def _fig_era_effect(post, out: Path, title: str) -> None:
    if az is None or 'era_effect' not in post:
        return
    era = post['era_effect']
    pre = era.sel(era_effect_dim_0=0).stack(draws=("chain", "draw")).values
    postv = era.sel(era_effect_dim_0=1).stack(draws=("chain", "draw")).values
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.violinplot([pre, postv], positions=[0, 1], showmeans=True)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['Pre-2007', 'Post-2007'])
    ax.set_title(title)
    plt.savefig(out, dpi=200, bbox_inches='tight'); plt.close()


def generate_bayesian_figs(model: str, base: Path, suffix: str, mirror_root: Optional[Path]) -> None:
    outdir = Path("figures/figures")
    _ensure_dir(outdir)

    # Load trace/PPD
    trace_path = base / ("trace_v3_6.nc" if model == "bayesian_v3_6" else "trace.nc")
    idata = az.from_netcdf(trace_path) if az is not None else None
    post = idata.posterior if idata is not None else None

    pp_path = base / 'posterior_predictive.csv'
    pp_df = pd.read_csv(pp_path) if pp_path.exists() else pd.DataFrame()

    # PPD overlay
    _fig_ppd_overlay(pp_df, outdir / f"ppd_overlay_{suffix}.png", f"Posterior predictive — {suffix}")

    # Forest HDI and parameter KDEs
    if post is not None:
        _fig_forest_hdi(post, outdir / f"forest_hdi_{suffix}.png", f"Posterior NAA/Cr by phase — {suffix}")
        _fig_param_kdes(post, outdir, suffix)
        _fig_era_effect(post, outdir / f"era_effect_violin_{suffix}.png", f"ART era effect — {suffix}")

    # Optional mirror
    if mirror_root is not None:
        _mirror_dir(outdir, mirror_root / 'figures')


def generate_ode_figs(base: Path, suffix: str, mirror_root: Optional[Path]) -> None:
    outdir = Path("figures/figures"); _ensure_dir(outdir)
    pred_path = base / 'predictions.csv'
    if not pred_path.exists():
        return
    pred = pd.read_csv(pred_path)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].bar(pred['Phase'], pred['NAA_Cr']); ax[0].set_title(f'ODE NAA/Cr by Phase ({suffix})'); ax[0].set_ylabel('NAA/Cr')
    ax[1].bar(pred['Phase'], pred['Cho_Cr']); ax[1].set_title(f'ODE Cho/Cr by Phase ({suffix})')
    plt.tight_layout(); plt.savefig(outdir / f"ode_phase_ratios_{suffix}.png", dpi=200, bbox_inches='tight'); plt.close()
    if mirror_root is not None:
        _mirror_dir(outdir, mirror_root / 'figures')


def main():
    args = parse_args()

    if args.model in ("bayesian_v3_6", "bayesian_v3"):
        model_dir = Path("results") / args.model / args.ratio / args.era
        run_dir = model_dir / args.run_id if args.run_id else _find_latest_run(model_dir)
        if run_dir is None:
            raise SystemExit(f"No runs found under {model_dir}")
        suffix = args.suffix or f"{args.ratio}_{args.era}_{args.model}"
        mirror_root = Path(args.mirror_root) if args.mirror_root else None
        generate_bayesian_figs(args.model, run_dir, suffix=suffix, mirror_root=mirror_root)
        print(f"Figures generated for {args.model} at {run_dir}")
    else:
        # enzyme_v3_ode
        model_dir = Path("results/enzyme_v3") / args.ratio / args.era
        run_dir = model_dir / args.run_id if args.run_id else _find_latest_run(model_dir)
        if run_dir is None:
            raise SystemExit(f"No enzyme_v3 runs found under {model_dir}")
        suffix = args.suffix or f"{args.ratio}_{args.era}_ode"
        mirror_root = Path(args.mirror_root) if args.mirror_root else None
        generate_ode_figs(run_dir, suffix=suffix, mirror_root=mirror_root)
        print(f"ODE figures generated for {run_dir}")


if __name__ == "__main__":
    main()
