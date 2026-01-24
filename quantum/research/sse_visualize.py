import argparse
import json
import os
from typing import Optional, Dict, Any

import matplotlib.pyplot as plt
import numpy as np

"""
Lightweight visualization helper for SSE outputs.

Capabilities (select via subcommands):
- coherence: Plot coherence time series from CSV (reg vs fib).
- summary: Plot final |psi|^2 heatmaps and overlay with Gamma map using saved NPZ files;
           falls back to existing overlay PNG if NPZs are missing.
- phase: Plot mean and 10–90% variance bands from sse_phase_sweep_summary.json.
- kernel: Preview SSE correlated kernel from *_sse_kernel.npz.

Usage examples:
  # 1) Coherence time series from CSV
  python Extra/sse_visualize.py coherence --csv datafiles/microtubule_simulation_..._sse_coherence.csv

  # 2) Summary overlay from *_summary.json
  python Extra/sse_visualize.py summary --summary datafiles/microtubule_simulation_..._summary.json

  # 3) Phase-sweep bands
  python Extra/sse_visualize.py phase --summary-json datafiles/sse_phase_sweep_summary.json

  # 4) Kernel preview
  python Extra/sse_visualize.py kernel --kernel datafiles/microtubule_simulation_..._sse_kernel.npz

Outputs are written under figures/ by default and mirrored to ~/Desktop/microtubule_simulation/figures when possible.
"""


def ensure_dirs(path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)


def mirror_to_desktop(fig_path: str):
    try:
        desktop_dir = os.path.expanduser('~/Desktop/microtubule_simulation/figures')
        os.makedirs(desktop_dir, exist_ok=True)
        dst = os.path.join(desktop_dir, os.path.basename(fig_path))
        import shutil
        shutil.copyfile(fig_path, dst)
    except Exception:
        pass


# ---- Coherence time series ----

def plot_coherence(csv_path: str, out: Optional[str] = None):
    import csv
    times = []
    reg = []
    fib = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row['time']))
            reg.append(float(row['sse_coherence_reg']))
            fib.append(float(row['sse_coherence_fib']))
    times = np.array(times)
    reg = np.array(reg)
    fib = np.array(fib)

    plt.figure(figsize=(7, 4))
    plt.plot(times, reg, label='reg', color='#1565C0')
    plt.plot(times, fib, label='fib', color='#EF6C00')
    plt.xlabel('time (sim units)')
    plt.ylabel('SSE coherence overlap')
    plt.title('SSE coherence over time')
    plt.legend()
    plt.tight_layout()

    if out is None:
        figures_dir = 'figures'
        os.makedirs(figures_dir, exist_ok=True)
        base = os.path.basename(csv_path).replace('_sse_coherence.csv', '')
        out = os.path.join(figures_dir, f'{base}_sse_coherence.png')
    ensure_dirs(out)
    plt.savefig(out, dpi=150)
    mirror_to_desktop(out)
    plt.close()
    print('Saved coherence plot to', out)


# ---- Summary / NPZ overlays ----

def infer_base_from_summary(summary_path: str) -> str:
    name = os.path.basename(summary_path)
    if name.endswith('_summary.json'):
        return name[:-len('_summary.json')]
    return os.path.splitext(name)[0]


def load_npz_or_none(path: str) -> Optional[Dict[str, Any]]:
    try:
        return dict(np.load(path))
    except Exception:
        return None


def first_present(mapping: Dict[str, Any], keys) -> Optional[Any]:
    """Return the first present value for a list of keys in a mapping (no array truth checks)."""
    for k in keys:
        if k in mapping:
            return mapping[k]
    return None


def plot_summary(summary_path: str, out_prefix: Optional[str] = None):
    with open(summary_path, 'r') as f:
        summary = json.load(f)
    base = infer_base_from_summary(summary_path)
    data_dir = os.path.dirname(summary_path)
    # Expected NPZ files
    arrays_npz = os.path.join(data_dir, f'{base}_arrays.npz')
    gamma_npz = os.path.join(data_dir, f'{base}_gamma_map_sse.npz')

    arrays = load_npz_or_none(arrays_npz)
    gamma = load_npz_or_none(gamma_npz)

    figures_dir = 'figures'
    os.makedirs(figures_dir, exist_ok=True)

    # 1) If arrays available, plot final |psi|^2 (reg, fib)
    if arrays is not None:
        # Try to fetch arrays with flexible keys (avoid array truth-value evaluation)
        psi_reg = first_present(arrays, ['Psi_reg_final', 'psi_reg_final', 'psi_reg'])
        psi_fib = first_present(arrays, ['Psi_fib_final', 'psi_fib_final', 'psi_fib'])
        R = first_present(arrays, ['R'])
        Z = first_present(arrays, ['Z'])
        r = first_present(arrays, ['r'])
        z = first_present(arrays, ['z'])
        extent = None
        if r is not None and z is not None:
            extent = [float(np.min(r)), float(np.max(r)), float(np.min(z)), float(np.max(z))]
        elif R is not None and Z is not None:
            extent = [float(np.min(R)), float(np.max(R)), float(np.min(Z)), float(np.max(Z))]

        if psi_reg is not None:
            plt.figure(figsize=(6, 4))
            plt.imshow(np.abs(psi_reg) ** 2, origin='lower', aspect='auto', cmap='viridis', extent=extent)
            plt.colorbar(label='|psi|^2 (reg)')
            plt.xlabel('r')
            plt.ylabel('z')
            plt.title('|psi|^2 (regular)')
            out = os.path.join(figures_dir, f'{base}_psi_reg.png') if not out_prefix else f'{out_prefix}_psi_reg.png'
            ensure_dirs(out)
            plt.tight_layout()
            plt.savefig(out, dpi=150)
            mirror_to_desktop(out)
            plt.close()
            print('Saved', out)
        if psi_fib is not None:
            plt.figure(figsize=(6, 4))
            plt.imshow(np.abs(psi_fib) ** 2, origin='lower', aspect='auto', cmap='viridis', extent=extent)
            plt.colorbar(label='|psi|^2 (fib)')
            plt.xlabel('r')
            plt.ylabel('z')
            plt.title('|psi|^2 (fibonacci)')
            out = os.path.join(figures_dir, f'{base}_psi_fib.png') if not out_prefix else f'{out_prefix}_psi_fib.png'
            ensure_dirs(out)
            plt.tight_layout()
            plt.savefig(out, dpi=150)
            mirror_to_desktop(out)
            plt.close()
            print('Saved', out)

    # 2) Overlay: |psi|^2 (reg) with Gamma map
    if arrays is not None and gamma is not None:
        psi_reg = first_present(arrays, ['Psi_reg_final', 'psi_reg_final', 'psi_reg'])
        r = first_present(arrays, ['r'])
        z = first_present(arrays, ['z'])
        Gamma = first_present(gamma, ['Gamma_map_sse', 'Gamma'])
        if psi_reg is not None and Gamma is not None:
            extent = None
            if r is not None and z is not None:
                extent = [float(np.min(r)), float(np.max(r)), float(np.min(z)), float(np.max(z))]
            plt.figure(figsize=(6, 4))
            im0 = plt.imshow(np.abs(psi_reg) ** 2, origin='lower', aspect='auto', cmap='viridis', extent=extent)
            im1 = plt.imshow(Gamma, origin='lower', aspect='auto', cmap='inferno', alpha=0.35, extent=extent)
            plt.xlabel('r')
            plt.ylabel('z')
            plt.title('Final |psi|^2 (reg) with Gamma overlay')
            cbar0 = plt.colorbar(im0, fraction=0.046, pad=0.04)
            cbar0.set_label('|psi|^2')
            cbar1 = plt.colorbar(im1, fraction=0.046, pad=0.12)
            cbar1.set_label('Gamma')
            out = os.path.join(figures_dir, f'{base}_overlay.png') if not out_prefix else f'{out_prefix}_overlay.png'
            ensure_dirs(out)
            plt.tight_layout()
            plt.savefig(out, dpi=150)
            mirror_to_desktop(out)
            plt.close()
            print('Saved', out)
    elif gamma is None:
        # Try to reference existing overlay saved by simulator
        possible_overlay = os.path.join('figures', f'{base}_overlay.png')
        if os.path.exists(possible_overlay):
            print('Overlay already exists at', possible_overlay)
        else:
            print('Gamma NPZ not found; skipping overlay generation.')


# ---- Phase sweep visualization ----

def plot_phase_sweep(summary_json: str, out: Optional[str] = None):
    with open(summary_json, 'r') as f:
        data = json.load(f)
    phases = data.get('phases') or list((data.get('stats_per_phase') or {}).keys())
    stats = data.get('stats_per_phase', {})
    if not phases or not stats:
        raise ValueError('Invalid phase-sweep summary JSON: missing phases/stats_per_phase')

    plt.figure(figsize=(8, 5))
    colors = {
        'none': '#4CAF50',
        'art_controlled': '#2196F3',
        'chronic': '#FF9800',
        'acute': '#F44336',
    }
    for ph in phases:
        entry = stats.get(ph)
        if not entry:
            continue
        times = np.array(entry['times'])
        mean = np.array(entry['mean'])
        p10 = np.array(entry['p10'])
        p90 = np.array(entry['p90'])
        c = colors.get(ph, None)
        plt.plot(times, mean, label=ph, color=c)
        plt.fill_between(times, p10, p90, color=c, alpha=0.2)
    plt.xlabel('time (sim units)')
    plt.ylabel('coherence overlap (reg)')
    plt.title('SSE coherence: mean and 10–90% bands across HIV phases')
    plt.legend()
    plt.tight_layout()

    if out is None:
        out = os.path.join('figures', 'sse_phase_sweep_viz.png')
    ensure_dirs(out)
    plt.savefig(out, dpi=150)
    mirror_to_desktop(out)
    plt.close()
    print('Saved phase-sweep plot to', out)


# ---- Kernel preview ----

def _find_latest_kernel(default_dir: str = 'datafiles') -> Optional[str]:
    try:
        candidates = [os.path.join(default_dir, f) for f in os.listdir(default_dir) if f.endswith('_sse_kernel.npz')]
        if not candidates:
            return None
        candidates.sort(key=lambda p: os.path.getmtime(p), reverse=True)
        return candidates[0]
    except Exception:
        return None


def plot_kernel(kernel_npz: str, out: Optional[str] = None):
    # If provided path does not exist, try to auto-detect latest kernel under datafiles/
    if not os.path.exists(kernel_npz):
        auto = _find_latest_kernel('datafiles')
        if auto and os.path.exists(auto):
            print(f"Provided kernel not found: {kernel_npz}. Using latest detected kernel: {auto}")
            kernel_npz = auto
        else:
            raise FileNotFoundError(
                f"Kernel NPZ not found at {kernel_npz}. No *_sse_kernel.npz found under datafiles/. "
                f"Note: kernels are only saved for SSE_correlated runs. Try: make viz-latest-kernel"
            )

    data = np.load(kernel_npz)
    # Try common keys without evaluating array truthiness
    if isinstance(data, np.lib.npyio.NpzFile):
        mapping = {k: data[k] for k in data.files}
    else:
        # Fallback, though np.load with NPZ returns NpzFile
        mapping = data
    K = first_present(mapping, ['sse_kernel', 'kernel'])
    if K is None:
        # If single array without key
        if isinstance(data, np.lib.npyio.NpzFile) and len(list(data.files)) == 1:
            K = data[list(data.files)[0]]
        else:
            raise ValueError('Kernel array not found in NPZ (expected key sse_kernel or kernel)')

    K = np.array(K, dtype=float)
    plt.figure(figsize=(5, 4))
    plt.imshow(K, origin='lower', cmap='magma')
    plt.colorbar(label='kernel weight')
    plt.title('SSE correlated kernel')
    plt.tight_layout()

    if out is None:
        base = os.path.basename(kernel_npz).replace('_sse_kernel.npz', '')
        out = os.path.join('figures', f'{base}_sse_kernel.png')
    ensure_dirs(out)
    plt.savefig(out, dpi=150)
    mirror_to_desktop(out)
    plt.close()
    print('Saved kernel image to', out)


# ---- CLI ----

def main(argv=None):
    p = argparse.ArgumentParser(description='SSE visualization helper')
    sub = p.add_subparsers(dest='cmd', required=True)

    p_coh = sub.add_parser('coherence', help='Plot coherence time series from CSV')
    p_coh.add_argument('--csv', required=True, help='Path to *_sse_coherence.csv')
    p_coh.add_argument('--out', default=None, help='Output PNG path')

    p_sum = sub.add_parser('summary', help='Plot final |psi|^2 and Gamma overlay from *_summary.json and NPZ files')
    p_sum.add_argument('--summary', required=True, help='Path to *_summary.json')
    p_sum.add_argument('--out_prefix', default=None, help='Optional prefix for output file names')

    p_phase = sub.add_parser('phase', help='Plot phase-sweep mean and 10–90% bands')
    p_phase.add_argument('--summary-json', default='datafiles/sse_phase_sweep_summary.json', help='Path to phase-sweep summary JSON')
    p_phase.add_argument('--out', default=None, help='Output PNG path')

    p_kern = sub.add_parser('kernel', help='Plot SSE kernel from NPZ')
    p_kern.add_argument('--kernel', required=True, help='Path to *_sse_kernel.npz')
    p_kern.add_argument('--out', default=None, help='Output PNG path')

    # Geometry comparison (reg vs fib) from CSV
    p_geom = sub.add_parser('geometry', help='Plot reg vs fib coherence from *_sse_coherence.csv (with optional Δ panel)')
    p_geom.add_argument('--csv', required=True, help='Path to *_sse_coherence.csv')
    p_geom.add_argument('--out', default=None, help='Output PNG path')
    p_geom.add_argument('--no-delta', action='store_true', help='Disable Δ(t)=fib−reg panel')

    args = p.parse_args(argv)

    if args.cmd == 'coherence':
        plot_coherence(args.csv, args.out)
    elif args.cmd == 'summary':
        plot_summary(args.summary, args.out_prefix)
    elif args.cmd == 'phase':
        plot_phase_sweep(args.summary_json, args.out)
    elif args.cmd == 'kernel':
        plot_kernel(args.kernel, args.out)
    elif args.cmd == 'geometry':
        # Lazy import to avoid heavy deps if not used
        from quantum.geometry_compare import plot_geometry as _plot_geometry
        _plot_geometry(args.csv, args.out, with_delta=(not args.no_delta))


if __name__ == '__main__':
    main()
