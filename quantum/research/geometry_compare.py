import argparse
import json
import os
from typing import Optional, Dict, Any, Tuple

import numpy as np

"""
Geometry comparison helper (Fibonacci vs. uniform/regular).

Capabilities:
- run: consume a single *_summary.json (and optionally its *_sse_coherence.csv),
       compute γ_fit and half-life for reg and fib, and produce Markdown + JSON
       reports highlighting differences.
- batch: aggregate across multiple summaries via a glob pattern; compute mean±std
         and 95% CI for γ_fit and half-life by geometry.
- viz: plot reg vs fib coherence from a CSV (or from summary if CSV not present),
       including an optional Δ(t) = reg − fib panel.

Outputs are saved under results/ and figures/ and mirrored to
~/Desktop/microtubule_simulation/ when possible.
"""


def ensure_dir(path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)


def mirror(path: str, subdir: str):
    try:
        d = os.path.expanduser(os.path.join('~/Desktop/microtubule_simulation', subdir))
        os.makedirs(d, exist_ok=True)
        import shutil
        shutil.copyfile(path, os.path.join(d, os.path.basename(path)))
    except Exception:
        pass


# --- Fitting utility (robust) ---

def fit_exponential_decay(times, coherences) -> Optional[float]:
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


# --- Core comparison helpers ---

def load_summary(summary_path: str) -> Dict[str, Any]:
    with open(summary_path, 'r') as f:
        return json.load(f)


def load_coherence_from_csv(csv_path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    import csv
    t, reg, fib = [], [], []
    with open(csv_path, 'r') as f:
        r = csv.DictReader(f)
        for row in r:
            t.append(float(row['time']))
            reg.append(float(row['sse_coherence_reg']))
            fib.append(float(row['sse_coherence_fib']))
    return np.array(t), np.array(reg), np.array(fib)


def compute_geom_fits(times: np.ndarray, reg: np.ndarray, fib: np.ndarray) -> Dict[str, Any]:
    g_reg = fit_exponential_decay(times, reg)
    g_fib = fit_exponential_decay(times, fib)
    th_reg = float(np.log(2.0) / g_reg) if (g_reg and g_reg > 0) else None
    th_fib = float(np.log(2.0) / g_fib) if (g_fib and g_fib > 0) else None
    diff_g = (g_fib - g_reg) if (g_reg is not None and g_fib is not None) else None
    ratio_g = (g_fib / g_reg) if (g_reg not in (None, 0) and g_fib is not None) else None
    diff_th = (th_fib - th_reg) if (th_reg is not None and th_fib is not None) else None
    ratio_th = (th_fib / th_reg) if (th_reg not in (None, 0) and th_fib is not None) else None
    return {
        'gamma_reg': g_reg,
        'gamma_fib': g_fib,
        't_half_reg': th_reg,
        't_half_fib': th_fib,
        'gamma_diff_fib_minus_reg': diff_g,
        'gamma_ratio_fib_over_reg': ratio_g,
        't_half_diff_fib_minus_reg': diff_th,
        't_half_ratio_fib_over_reg': ratio_th,
    }


def bootstrapped_delta(times: np.ndarray, reg: np.ndarray, fib: np.ndarray, B: int = 500, seed: int = 1234) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)
    n = min(reg.size, fib.size, times.size)
    if n < 5:
        return {'delta_mean': None, 'delta_ci95': None}
    reg = reg[:n]
    fib = fib[:n]
    delta = fib - reg
    idx = np.arange(n)
    boots = []
    for _ in range(B):
        sample_idx = rng.choice(idx, size=n, replace=True)
        boots.append(float(np.nanmean(delta[sample_idx])))
    boots = np.array(boots)
    return {
        'delta_mean': float(np.nanmean(delta)),
        'delta_ci95': [float(np.nanpercentile(boots, 2.5)), float(np.nanpercentile(boots, 97.5))],
    }


def compare_single(summary_path: str, out_prefix: Optional[str] = None, csv_hint: Optional[str] = None) -> Dict[str, Any]:
    summary = load_summary(summary_path)
    base = os.path.basename(summary_path).replace('_summary.json', '')
    data_dir = os.path.dirname(summary_path)
    # Prefer CSV for SSE coherence; fallback to arrays from summary
    csv_path = csv_hint or os.path.join(data_dir, f'{base}_sse_coherence.csv')
    if os.path.exists(csv_path):
        times, reg, fib = load_coherence_from_csv(csv_path)
    else:
        # Fallback to time series in summary (SSE or legacy)
        times = np.array(summary.get('timestamps', []), dtype=float)
        reg = np.array(summary.get('sse_coherence_reg') or summary.get('coherence_reg') or [], dtype=float)
        fib = np.array(summary.get('sse_coherence_fib') or summary.get('coherence_fib') or [], dtype=float)
    fits = compute_geom_fits(times, reg, fib)
    delta_stats = bootstrapped_delta(times, reg, fib)

    findings = {
        'source_summary': summary_path,
        'config': summary.get('config', {}),
        'fits': fits,
        'delta_stats': delta_stats,
    }

    # Write Markdown + JSON
    out_base = out_prefix if out_prefix else os.path.join('results', f'{base}_geom')
    md = f"{out_base}_comparison.md"
    js = f"{out_base}_comparison.json"
    ensure_dir(md)
    ensure_dir(js)

    lines = []
    lines.append('# Geometry comparison — Fibonacci vs. regular')
    lines.append('')
    lines.append(f"Source summary: {summary_path}")
    lines.append('')
    gr = fits.get('gamma_reg'); gf = fits.get('gamma_fib')
    tr = fits.get('t_half_reg'); tf = fits.get('t_half_fib')
    lines.append(f"- γ_fit (reg) = {gr:.3e}" if gr is not None else "- γ_fit (reg) = n/a")
    lines.append(f"- γ_fit (fib) = {gf:.3e}" if gf is not None else "- γ_fit (fib) = n/a")
    if gr is not None and gf is not None and gr != 0:
        lines.append(f"- γ ratio (fib/reg) = {gf/gr:.3f}")
    lines.append(f"- t1/2 (reg) = {tr:.3g}" if tr is not None else "- t1/2 (reg) = n/a")
    lines.append(f"- t1/2 (fib) = {tf:.3g}" if tf is not None else "- t1/2 (fib) = n/a")
    if tr is not None and tf is not None and tr != 0:
        lines.append(f"- t1/2 ratio (fib/reg) = {tf/tr:.3f}")
    dm = delta_stats.get('delta_mean'); dci = delta_stats.get('delta_ci95')
    if dm is not None and dci:
        lines.append(f"- Mean coherence gap Δ(t)=fib−reg ≈ {dm:.3e} (95% CI [{dci[0]:.3e}, {dci[1]:.3e}])")
    lines.append('')
    lines.append('Interpretation:')
    if gf is not None and gr is not None and gf < gr:
        lines.append('- Faster decay on Fibonacci than regular (γ_fib < γ_reg).')
    elif gf is not None and gr is not None and gf > gr:
        lines.append('- Slower decay on Fibonacci than regular (γ_fib > γ_reg), consistent with prolonged decoherence hypothesis if sustained across runs.')
    else:
        lines.append('- Insufficient data to compare decay rates robustly.')

    with open(md, 'w') as f:
        f.write("\n".join(lines) + "\n")
    with open(js, 'w') as f:
        json.dump(findings, f, indent=2)
    mirror(md, 'results')
    mirror(js, 'results')
    print('Wrote:', md)
    print('Wrote:', js)
    return findings


def aggregate_batch(glob_pattern: str, out_json: str = 'results/geometry_compare_aggregate.json') -> Dict[str, Any]:
    import glob
    files = glob.glob(glob_pattern)
    rows = []
    for fp in files:
        try:
            s = load_summary(fp)
            times = np.array(s.get('timestamps', []), dtype=float)
            reg = np.array(s.get('sse_coherence_reg') or s.get('coherence_reg') or [], dtype=float)
            fib = np.array(s.get('sse_coherence_fib') or s.get('coherence_fib') or [], dtype=float)
            fits = compute_geom_fits(times, reg, fib)
            rows.append({'src': fp, **fits})
        except Exception:
            continue
    def agg(key: str):
        vals = np.array([r[key] for r in rows if r.get(key) is not None], dtype=float)
        if vals.size == 0:
            return None, None, [None, None]
        mean = float(np.nanmean(vals)); std = float(np.nanstd(vals))
        ci = [float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))]
        return mean, std, ci
    summary = {
        'n_files': len(files),
        'gamma_reg': agg('gamma_reg'),
        'gamma_fib': agg('gamma_fib'),
        't_half_reg': agg('t_half_reg'),
        't_half_fib': agg('t_half_fib'),
        'gamma_ratio_fib_over_reg': agg('gamma_ratio_fib_over_reg'),
        't_half_ratio_fib_over_reg': agg('t_half_ratio_fib_over_reg'),
        'rows': rows,
    }
    ensure_dir(out_json)
    with open(out_json, 'w') as f:
        json.dump(summary, f, indent=2)
    mirror(out_json, 'results')
    print('Aggregated geometry comparison written to', out_json)
    return summary


# --- Visualization ---

def plot_geometry(csv_path: str, out_png: Optional[str] = None, with_delta: bool = True):
    import matplotlib.pyplot as plt
    times, reg, fib = load_coherence_from_csv(csv_path)
    fig, ax = plt.subplots(2 if with_delta else 1, 1, figsize=(8, 5), sharex=True)
    if with_delta:
        ax0 = ax[0]; ax1 = ax[1]
    else:
        ax0 = ax
    ax0.plot(times, reg, label='regular', color='#1565C0')
    ax0.plot(times, fib, label='fibonacci', color='#EF6C00')
    ax0.set_ylabel('coherence overlap')
    ax0.set_title('Geometry coherence comparison')
    ax0.legend()
    if with_delta:
        delta = fib - reg
        ax1.plot(times, delta, color='#6A1B9A')
        ax1.axhline(0.0, color='gray', lw=0.8)
        ax1.set_xlabel('time (sim units)')
        ax1.set_ylabel('Δ(t) fib−reg')
    else:
        ax0.set_xlabel('time (sim units)')
    plt.tight_layout()
    if out_png is None:
        base = os.path.basename(csv_path).replace('_sse_coherence.csv', '')
        out_png = os.path.join('figures', f'{base}_geometry_compare.png')
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=150)
    mirror(out_png, 'figures')
    plt.close()
    print('Saved geometry comparison plot to', out_png)


# --- CLI ---

def main(argv=None):
    p = argparse.ArgumentParser(description='Geometry comparison helper (Fibonacci vs. regular)')
    sub = p.add_subparsers(dest='cmd', required=True)

    pr = sub.add_parser('run', help='Compare geometries for a single run summary JSON')
    pr.add_argument('--summary', required=True, help='Path to *_summary.json')
    pr.add_argument('--csv', default=None, help='Optional path to *_sse_coherence.csv (if not provided, inferred)')
    pr.add_argument('--out_prefix', default=None, help='Optional output prefix for reports')

    pb = sub.add_parser('batch', help='Aggregate geometry comparison across multiple summaries (glob)')
    pb.add_argument('--glob', default='datafiles/*_summary.json')
    pb.add_argument('--out', default='results/geometry_compare_aggregate.json')

    pv = sub.add_parser('viz', help='Plot coherence reg vs fib from CSV (with optional Δ(t) panel)')
    pv.add_argument('--csv', required=True, help='Path to *_sse_coherence.csv')
    pv.add_argument('--out', default=None, help='Output PNG path')
    pv.add_argument('--no-delta', action='store_true', help='Disable Δ(t) panel')

    args = p.parse_args(argv)
    if args.cmd == 'run':
        compare_single(args.summary, args.out_prefix, csv_hint=args.csv)
        return
    if args.cmd == 'batch':
        aggregate_batch(args.glob, args.out)
        return
    if args.cmd == 'viz':
        plot_geometry(args.csv, args.out, with_delta=(not args.no_delta))
        return


if __name__ == '__main__':
    main()
