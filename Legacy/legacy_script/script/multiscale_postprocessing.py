import os
import json
import csv
from dataclasses import dataclass, asdict
from typing import List, Optional, Dict, Any, Tuple
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime


"""
Multiscale post-processing pipeline (Levels 2–5)

Purpose:
- Consume quantum simulation outputs (SSE coherence, variance, config) produced by
  quantum/microtubule_quantum_coherence_full_simulation.py
- Map quantum coherence metrics to transport efficiency (Level 2), mitochondrial
  ATP production (Level 3), NAA synthesis (Level 4), and a simple MRS proxy for
  NAA/Cr ratio (Level 5).

Inputs:
- <base>_summary.json (required): contains config, timestamps, coherence/variance series
- <base>_sse_coherence.csv (optional): time, sse_coherence_reg, sse_coherence_fib

Outputs:
- <base>_multiscale_summary.json: compiled downstream metrics and parameters
- <base>_multiscale_timeseries.csv: time-series of η, J, ATP, NAA, NAA/Cr
- <base>_multiscale_overview.png: quicklook plot of key series

Usage:
- python Extra/script/multiscale_postprocessing.py --summary <path_to_summary.json> [--xi 0.8] [--xi_acute 0.4]

Notes:
- Parameters are placeholders calibrated to be numerically reasonable and
  should be refined against experimental data.
- The ξ (correlation length) can modulate the effective coherence contribution.
"""


@dataclass
class TransportParams:
    eta_0: float = 0.5           # base transport efficiency (quantum-independent)
    alpha: float = 0.6           # weight for delocalization factor
    beta: float = 0.8            # weight for coherence factor
    sigma_r_baseline: float = 0.38  # nm, typical regular lattice delocalization proxy
    coherence_baseline: float = 1.0  # dimensionless baseline


@dataclass
class SubstrateParams:
    D_eff: float = 0.5           # um^2/s, effective cytosolic diffusion
    tau_halflife: float = 600.0  # s, degradation half-life scale (placeholder)
    compartment_factor: float = 1.0  # unitless, delivery compartment scaling
    distance: float = 10.0       # um, soma→dendrite effective path


@dataclass
class MetabolismParams:
    # Linear saturation mapping S -> ATP for simplicity
    k_ATP: float = 1.0           # scaling constant to map S to ATP rate (a.u.)
    ATP_max: float = 5.0         # cap for ATP rate (a.u.)


@dataclass
class NAAParams:
    Km_asp: float = 0.5          # mM
    Km_acCoA: float = 0.1        # mM
    Vmax: float = 10.0           # μmol/min/g protein (relative scale)
    K_ATP: float = 1.0           # a.u., ATP activation constant
    asp_conc: float = 2.0        # mM, placeholder pool
    acCoA_conc: float = 0.2      # mM, placeholder pool


@dataclass
class MRSParams:
    Cr_baseline: float = 1.0     # arbitrary concentration units
    T2_NAA_ms: float = 260.0     # ms, rough at 3T
    T2_Cr_ms: float = 180.0      # ms, rough at 3T
    TE_ms: float = 35.0          # ms (PRESS)


@dataclass
class PipelineParams:
    transport: TransportParams = TransportParams()
    substrate: SubstrateParams = SubstrateParams()
    metabolism: MetabolismParams = MetabolismParams()
    naa: NAAParams = NAAParams()
    mrs: MRSParams = MRSParams()
    xi_healthy_nm: float = 0.8
    xi_acute_nm: float = 0.4
    xi_sensitivity: float = 0.5  # how strongly ξ modulates coherence contribution (0..1)


def load_summary(summary_path: str) -> Dict[str, Any]:
    with open(summary_path, 'r') as f:
        return json.load(f)


def load_coherence_csv(csv_path: str) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    if not os.path.exists(csv_path):
        return None
    times, c_reg, c_fib = [], [], []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row['time']))
            c_reg.append(float(row['sse_coherence_reg']))
            c_fib.append(float(row['sse_coherence_fib']))
    return np.array(times), np.array(c_reg), np.array(c_fib)


def xi_modulated_coherence(coherence: np.ndarray, xi_nm: float, params: PipelineParams) -> np.ndarray:
    # Normalize ξ between acute and healthy, scale coherence contribution accordingly
    xi_min = min(params.xi_healthy_nm, params.xi_acute_nm)
    xi_max = max(params.xi_healthy_nm, params.xi_acute_nm)
    if xi_max <= xi_min:
        return coherence
    xi_norm = (xi_nm - xi_min) / (xi_max - xi_min)
    # Map: higher ξ → more correlated dephasing → stronger beneficial coherence factor
    scale = 1.0 - params.xi_sensitivity + params.xi_sensitivity * xi_norm
    return coherence * scale


def transport_efficiency(coherence: np.ndarray, sigma_r: np.ndarray, tp: TransportParams) -> np.ndarray:
    # Delocalization factor and coherence factor
    deloc = np.divide(sigma_r, tp.sigma_r_baseline, out=np.ones_like(sigma_r), where=(tp.sigma_r_baseline > 0))
    coh = np.divide(coherence, tp.coherence_baseline, out=np.ones_like(coherence), where=(tp.coherence_baseline > 0))
    eta = tp.eta_0 * (1.0 + tp.alpha * deloc) * (1.0 + tp.beta * coh)
    return np.clip(eta, 0.0, 1.0)


def motor_flux(eta: np.ndarray, flux_0: float = 1.0) -> np.ndarray:
    # Simple proportional mapping to cargo flux
    return flux_0 * eta


def substrate_delivery(flux: np.ndarray, sp: SubstrateParams) -> np.ndarray:
    tau_diff = (sp.distance ** 2) / (2.0 * sp.D_eff + 1e-12)
    # Convert halflife to exponential decay constant
    lam = np.log(2.0) / max(sp.tau_halflife, 1e-6)
    degradation = np.exp(-lam * tau_diff)
    S = flux * degradation * sp.compartment_factor
    return S


def atp_rate(S: np.ndarray, mp: MetabolismParams) -> np.ndarray:
    ATP = mp.k_ATP * S
    return np.minimum(ATP, mp.ATP_max)


def naa_synthesis(ATP: np.ndarray, np_params: NAAParams) -> np.ndarray:
    ATP_factor = (ATP ** 2) / (ATP ** 2 + np_params.K_ATP ** 2)
    v = np_params.Vmax * ATP_factor
    # Incorporate substrate saturation as constant pools (can be time-dependent in future)
    v *= (np_params.asp_conc / (np_params.Km_asp + np_params.asp_conc))
    v *= (np_params.acCoA_conc / (np_params.Km_acCoA + np_params.acCoA_conc))
    return v


def mrs_proxy(NAA: np.ndarray, mrs: MRSParams) -> np.ndarray:
    # Simple T2-weighted signal model at TE
    TE = mrs.TE_ms
    E2_NAA = np.exp(-TE / max(mrs.T2_NAA_ms, 1e-6))
    E2_Cr = np.exp(-TE / max(mrs.T2_Cr_ms, 1e-6))
    # Assume Cr ~ baseline constant pool for proxy
    Cr = np.full_like(NAA, fill_value=mrs.Cr_baseline)
    S_NAA = NAA * E2_NAA
    S_Cr = Cr * E2_Cr
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.divide(S_NAA, S_Cr, out=np.zeros_like(S_NAA), where=(S_Cr > 0))
    return ratio


def compute_sigma_r_from_variance(variance_series: List[float]) -> np.ndarray:
    # σ_r is sqrt of variance along r; input variance may already incorporate volume element
    arr = np.asarray(variance_series, dtype=float)
    arr = np.maximum(arr, 0.0)
    return np.sqrt(arr)


def run_pipeline(summary_path: str, xi_nm: Optional[float] = None, params: PipelineParams = PipelineParams()) -> Dict[str, Any]:
    summary = load_summary(summary_path)
    config = summary.get('config', {})
    timestamps = np.array(summary.get('timestamps', []), dtype=float)

    # Prefer SSE coherence if present; else fallback to generic coherence
    c_reg = np.array(summary.get('sse_coherence_reg', summary.get('coherence_reg', [])), dtype=float)
    c_fib = np.array(summary.get('sse_coherence_fib', summary.get('coherence_fib', [])), dtype=float)

    # Delocalization proxies from SSE variance if present; else compute from variance_reg
    if 'sse_variance_reg' in summary and summary['sse_variance_reg']:
        sigma_r_reg = compute_sigma_r_from_variance(summary['sse_variance_reg'])
        sigma_r_fib = compute_sigma_r_from_variance(summary['sse_variance_fib'])
    else:
        sigma_r_reg = compute_sigma_r_from_variance(summary.get('variance_reg', []))
        sigma_r_fib = compute_sigma_r_from_variance(summary.get('variance_fib', []))

    # Ensure equal length
    n = min(len(timestamps), len(c_reg), len(c_fib), len(sigma_r_reg), len(sigma_r_fib))
    timestamps = timestamps[:n]
    c_reg, c_fib = c_reg[:n], c_fib[:n]
    sigma_r_reg, sigma_r_fib = sigma_r_reg[:n], sigma_r_fib[:n]

    # Select ξ scenario
    xi_value = xi_nm if xi_nm is not None else float(config.get('corr_length_xi', params.xi_healthy_nm))

    # Apply ξ modulation to coherence
    c_reg_eff = xi_modulated_coherence(c_reg, xi_value, params)
    c_fib_eff = xi_modulated_coherence(c_fib, xi_value, params)

    # Level 2: transport efficiency and cargo flux
    eta_reg = transport_efficiency(c_reg_eff, sigma_r_reg, params.transport)
    eta_fib = transport_efficiency(c_fib_eff, sigma_r_fib, params.transport)
    J_reg = motor_flux(eta_reg)
    J_fib = motor_flux(eta_fib)

    # Level 2.5: substrate at mitochondria
    S_reg = substrate_delivery(J_reg, params.substrate)
    S_fib = substrate_delivery(J_fib, params.substrate)

    # Level 3: ATP production rate
    ATP_reg = atp_rate(S_reg, params.metabolism)
    ATP_fib = atp_rate(S_fib, params.metabolism)

    # Level 4: NAA synthesis (rate proxy proportional to concentration; integrate crudely)
    vNAA_reg = naa_synthesis(ATP_reg, params.naa)
    vNAA_fib = naa_synthesis(ATP_fib, params.naa)
    # Cumulative proxy for [NAA]
    # Use simple Euler integration with dt inferred from timestamps if available
    if len(timestamps) > 1:
        dt = np.diff(timestamps, prepend=timestamps[0])
    else:
        dt = np.ones_like(vNAA_reg)
    NAA_reg = np.cumsum(vNAA_reg * dt)
    NAA_fib = np.cumsum(vNAA_fib * dt)

    # Level 5: MRS proxy ratios
    NAA_over_Cr_reg = mrs_proxy(NAA_reg, params.mrs)
    NAA_over_Cr_fib = mrs_proxy(NAA_fib, params.mrs)

    result = {
        'config_used': config,
        'params': asdict(params),
        'xi_applied_nm': xi_value,
        'timestamps': timestamps.tolist(),
        'eta_reg': eta_reg.tolist(),
        'eta_fib': eta_fib.tolist(),
        'J_reg': J_reg.tolist(),
        'J_fib': J_fib.tolist(),
        'S_reg': S_reg.tolist(),
        'S_fib': S_fib.tolist(),
        'ATP_reg': ATP_reg.tolist(),
        'ATP_fib': ATP_fib.tolist(),
        'NAA_reg': NAA_reg.tolist(),
        'NAA_fib': NAA_fib.tolist(),
        'NAA_over_Cr_reg': NAA_over_Cr_reg.tolist(),
        'NAA_over_Cr_fib': NAA_over_Cr_fib.tolist(),
        'final_ratios': {
            'reg': float(NAA_over_Cr_reg[-1]) if len(NAA_over_Cr_reg) else None,
            'fib': float(NAA_over_Cr_fib[-1]) if len(NAA_over_Cr_fib) else None,
        }
    }
    return result


def save_outputs(summary_path: str, base_result: Dict[str, Any]) -> Tuple[str, str, str]:
    # Derive base filename and directories similar to simulator
    base_dir = os.path.dirname(summary_path)
    base_name = os.path.basename(summary_path).replace('_summary.json', '')

    # Save JSON
    out_json_proj = os.path.join(base_dir, f"{base_name}_multiscale_summary.json")
    with open(out_json_proj, 'w') as f:
        json.dump(base_result, f, indent=2)

    out_json_desktop = os.path.expanduser(f"~/Desktop/microtubule_simulation/datafiles/{base_name}_multiscale_summary.json")
    os.makedirs(os.path.dirname(out_json_desktop), exist_ok=True)
    with open(out_json_desktop, 'w') as f:
        json.dump(base_result, f, indent=2)

    # Save CSV time series
    out_csv_proj = os.path.join(base_dir, f"{base_name}_multiscale_timeseries.csv")
    headers = ['time','eta_reg','eta_fib','J_reg','J_fib','S_reg','S_fib','ATP_reg','ATP_fib','NAA_reg','NAA_fib','NAA_over_Cr_reg','NAA_over_Cr_fib']
    rows = zip(base_result['timestamps'], base_result['eta_reg'], base_result['eta_fib'], base_result['J_reg'], base_result['J_fib'], base_result['S_reg'], base_result['S_fib'], base_result['ATP_reg'], base_result['ATP_fib'], base_result['NAA_reg'], base_result['NAA_fib'], base_result['NAA_over_Cr_reg'], base_result['NAA_over_Cr_fib'])
    with open(out_csv_proj, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)

    out_csv_desktop = os.path.expanduser(f"~/Desktop/microtubule_simulation/datafiles/{base_name}_multiscale_timeseries.csv")
    with open(out_csv_desktop, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)

    # Quick overview plot
    out_fig_proj = os.path.join(base_dir, f"{base_name}_multiscale_overview.png")
    out_fig_desktop = os.path.expanduser(f"~/Desktop/microtubule_simulation/figures/{base_name}_multiscale_overview.png")
    os.makedirs(os.path.dirname(out_fig_desktop), exist_ok=True)
    try:
        t = np.array(base_result['timestamps'])
        fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
        axes[0].plot(t, base_result['eta_reg'], label='η reg')
        axes[0].plot(t, base_result['eta_fib'], label='η fib', ls='--')
        axes[0].set_ylabel('Transport η')
        axes[0].legend()

        axes[1].plot(t, base_result['ATP_reg'], label='ATP reg')
        axes[1].plot(t, base_result['ATP_fib'], label='ATP fib', ls='--')
        axes[1].set_ylabel('ATP rate (a.u.)')
        axes[1].legend()

        axes[2].plot(t, base_result['NAA_over_Cr_reg'], label='NAA/Cr reg')
        axes[2].plot(t, base_result['NAA_over_Cr_fib'], label='NAA/Cr fib', ls='--')
        axes[2].set_ylabel('NAA/Cr (proxy)')
        axes[2].set_xlabel('time (s)')
        axes[2].legend()

        fig.suptitle('Multiscale Post-Processing Overview')
        fig.tight_layout()
        fig.savefig(out_fig_proj, dpi=150)
        fig.savefig(out_fig_desktop, dpi=150)
        plt.close(fig)
    except Exception:
        # Non-fatal if plotting fails
        pass

    return out_json_proj, out_csv_proj, out_fig_proj


def _cli():
    import argparse
    parser = argparse.ArgumentParser(description='Multiscale post-processing pipeline (Levels 2–5)')
    parser.add_argument('--summary', required=True, help='Path to <base>_summary.json produced by the simulator')
    parser.add_argument('--xi', type=float, default=None, help='Correlation length ξ (nm) to apply in coherence modulation')
    parser.add_argument('--xi_acute', type=float, default=0.4, help='Acute ξ (nm) reference')
    parser.add_argument('--xi_healthy', type=float, default=0.8, help='Healthy ξ (nm) reference')
    parser.add_argument('--xi_sensitivity', type=float, default=0.5, help='Weight of ξ effect (0..1)')
    args = parser.parse_args()

    params = PipelineParams()
    params.xi_acute_nm = args.xi_acute
    params.xi_healthy_nm = args.xi_healthy
    params.xi_sensitivity = args.xi_sensitivity

    result = run_pipeline(args.summary, xi_nm=args.xi, params=params)
    out_json, out_csv, out_fig = save_outputs(args.summary, result)
    print(f"Saved: {out_json}\nSaved: {out_csv}\nSaved: {out_fig}")


if __name__ == '__main__':
    _cli()
