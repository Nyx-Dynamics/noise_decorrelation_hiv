"""
Option B (SSE) minimal open-systems scaffolding for Microtubule_Simulation.

Provides:
- build_dephasing_map(r, z, cytokine_field, hiv_phase, params) -> Gamma_map
- sse_dephasing_step(psi, Gamma_map, dt, rng=None) -> psi_new (local/uncorrelated)
- sse_dephasing_step_correlated(psi, Gamma_map, dt, kernel=None, xi=None, dr=None, dz=None, rng=None) -> psi_new
- build_gaussian_kernel(dr, dz, N_z, N_r, xi) -> 2D normalized kernel
- coherence_metrics_sse(psi, psi0, R_grid, Z_grid, dr, dz) -> dict

Notes:
- Physiological defaults are used where parameters are not provided:
  temperature_K=310.0, ionic_strength_M=0.15, dielectric_rel=80.0.
- Γ_map currently adapts legacy phase scaling patterns and cytokine coupling; units are
  effectively relative rates consistent with existing simulator usage. Future work can
  introduce SI-consistent scaling tied to Tegmark-style estimates.

References (see project docs for full citations): Tegmark (2000); Haken–Strobl–Reineker; Breuer–Petruccione.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Dict, Any

import numpy as np

from .shims import (
    calculate_probability_density,
    calculate_coherence,
    calculate_dispersion_metrics,
)


@dataclass
class EnvParams:
    # Physiological/environmental parameters
    temperature_K: float = 310.0
    ionic_strength_M: float = 0.15
    dielectric_rel: float = 80.0
    # Dephasing parameters
    Gamma_0: float = 0.05
    alpha_c: float = 0.1
    gamma_scale_alpha: float = 1.0
    corr_length_xi: float = 0.0  # not used in local SSE


_PHASE_MULTIPLIERS = {
    # These multipliers adapt patterns seen in legacy calculate_decoherence logic
    # Acute stronger scaling; chronic persistent; art moderate; none minimal.
    "acute": dict(alpha=0.28 * 4.0, base=1.0),
    "chronic": dict(alpha=0.9, base=1.1),
    "art_controlled": dict(alpha=0.8, base=1.0),
    "none": dict(alpha=1.0, base=1.0),
}


def build_dephasing_map(r: np.ndarray,
                        z: np.ndarray,
                        cytokine_field: np.ndarray,
                        hiv_phase: str,
                        params: Optional[Dict[str, Any]] = None) -> np.ndarray:
    """
    Construct a spatial dephasing-rate map Γ(r, z) informed by cytokine field and HIV phase.

    Args:
        r: 1D radial grid (size N_r)
        z: 1D axial grid (size N_z)
        cytokine_field: 2D array (N_z x N_r) with values in [0,1]
        hiv_phase: one of {"none","acute","art_controlled","chronic"}
        params: optional dict that may include {Gamma_0, alpha_c, gamma_scale_alpha,
                temperature_K, ionic_strength_M, dielectric_rel}

    Returns:
        Gamma_map: 2D array (N_z x N_r)
    """
    p = EnvParams(**(params or {}))

    # Normalize inputs
    hiv = hiv_phase if hiv_phase in _PHASE_MULTIPLIERS else "none"
    mult = _PHASE_MULTIPLIERS[hiv]

    # Effective alpha scaled by phase and optional global scale
    alpha_eff = p.alpha_c * mult["alpha"] * p.gamma_scale_alpha
    base_eff = p.Gamma_0 * mult["base"]

    # Cytokine-coupled dephasing (legacy-compatible exponential form)
    # Γ = Γ0' * exp(alpha_eff * C)
    Gamma_map = base_eff * np.exp(alpha_eff * cytokine_field)

    # Clamp to reasonable bounds to avoid numerical blowups
    Gamma_map = np.clip(Gamma_map, 0.0, 10.0)
    return Gamma_map


def sse_dephasing_step(psi: np.ndarray,
                       Gamma_map: np.ndarray,
                       dt: float,
                       rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """
    Apply a local SSE (phase-noise) dephasing step to the wavefunction.

    Ito/Euler–Maruyama discretization for pure dephasing unraveling:
    psi <- psi * exp(i * sqrt(2*Γ*dt) * ξ - 0.5*Γ*dt), with ξ ~ N(0,1) i.i.d.

    This preserves norm in expectation; we renormalize numerically each step.

    Args:
        psi: complex wavefunction (N_z x N_r)
        Gamma_map: dephasing rates (N_z x N_r)
        dt: time step
        rng: optional numpy Generator for reproducibility

    Returns:
        psi_new: updated wavefunction
    """
    if rng is None:
        rng = np.random.default_rng()

    # Ensure shapes compatible
    if psi.shape != Gamma_map.shape:
        raise ValueError("psi and Gamma_map must have the same shape")

    # Draw standard normal noise field
    xi = rng.standard_normal(size=psi.shape)

    # Compute phase-noise factor
    phase_std = np.sqrt(np.maximum(0.0, 2.0 * Gamma_map * dt))
    phase = phase_std * xi

    # Deterministic damping term for pure dephasing unraveling
    damp = np.exp(-0.5 * Gamma_map * dt)

    # Apply multiplicative update
    psi_new = psi * damp * np.exp(1j * phase)

    # Renormalize to unit probability (using cylindrical-aware approach handled by caller)
    return psi_new


def coherence_metrics_sse(psi: np.ndarray,
                          psi0: np.ndarray,
                          R_grid: np.ndarray,
                          Z_grid: np.ndarray,
                          dr: float,
                          dz: float) -> Dict[str, Any]:
    """
    Compute basic coherence/dispersion diagnostics using cylindrical-aware metrics.

    Returns: dict with keys {coherence_overlap, variance, entropy, kurtosis}
    """
    # Probability densities with Jacobian baked by metrics when using R_grid in overlap
    prob_current = calculate_probability_density(psi, R_grid)
    prob_initial = calculate_probability_density(psi0, R_grid)

    coherence_overlap = calculate_coherence(prob_current, prob_initial, dr, dz, R_grid)
    dispersion = calculate_dispersion_metrics(psi, R_grid, Z_grid, dr, dz)

    return {
        "coherence_overlap": float(coherence_overlap),
        **dispersion,
    }


# --- Correlated SSE support ---

def build_gaussian_kernel(dr: float, dz: float, N_z: int, N_r: int, xi: float) -> np.ndarray:
    """
    Build a separable 2D Gaussian kernel approximating spatial correlation with length ξ.
    The kernel is normalized to sum to 1. Uses same physical length for r and z.
    If xi <= 0, returns a 1x1 kernel (delta), i.e., no correlation.
    """
    if xi is None or xi <= 0:
        return np.array([[1.0]], dtype=float)

    # Convert correlation length (space units) to standard deviations in index units
    sigma_r = max(1.0, xi / max(dr, 1e-12))
    sigma_z = max(1.0, xi / max(dz, 1e-12))

    # Limit kernel extents to 3 sigma each side
    half_r = int(min(N_r // 2, np.ceil(3.0 * sigma_r)))
    half_z = int(min(N_z // 2, np.ceil(3.0 * sigma_z)))

    r_idx = np.arange(-half_r, half_r + 1)
    z_idx = np.arange(-half_z, half_z + 1)

    kr = np.exp(-0.5 * (r_idx / sigma_r) ** 2)
    kz = np.exp(-0.5 * (z_idx / sigma_z) ** 2)

    kernel = np.outer(kz, kr)
    s = kernel.sum()
    if s > 0:
        kernel = kernel / s
    return kernel.astype(float)


def _fft_convolve2d(field: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """
    Fast 2D convolution via FFT with wrap-around (circular) boundary conditions.
    Shapes: field (N_z x N_r), kernel (kz x kr). Kernel is zero-padded to field size.
    """
    Nz, Nr = field.shape
    kz, kr = kernel.shape
    # Zero-pad kernel to field size with wrap-around center
    pad = np.zeros_like(field, dtype=float)
    # Place kernel at top-left; FFT assumes (0,0) at origin; roll to center after FFT mult
    pad[:kz, :kr] = kernel
    F_field = np.fft.rfft2(field)
    F_kern = np.fft.rfft2(np.fft.ifftshift(pad))
    conv = np.fft.irfft2(F_field * F_kern, s=field.shape)
    return conv.real


def sse_dephasing_step_correlated(psi: np.ndarray,
                                   Gamma_map: np.ndarray,
                                   dt: float,
                                   kernel: Optional[np.ndarray] = None,
                                   xi: Optional[float] = None,
                                   dr: Optional[float] = None,
                                   dz: Optional[float] = None,
                                   rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """
    Apply an SSE dephasing step with spatially correlated phase noise.

    The correlated standard normal field ξ_corr is obtained by filtering white noise with a
    Gaussian kernel (normalized to unit variance afterwards). Then use the same phase update:
        psi <- psi * exp(i * sqrt(2*Γ*dt) * ξ_corr - 0.5*Γ*dt)

    If kernel is None and xi>0 with dr,dz provided, a kernel is built on-the-fly.
    If xi<=0, falls back to local (uncorrelated) SSE.
    """
    if rng is None:
        rng = np.random.default_rng()

    if psi.shape != Gamma_map.shape:
        raise ValueError("psi and Gamma_map must have the same shape")

    Nz, Nr = psi.shape

    if kernel is None:
        if xi is None or xi <= 0 or dr is None or dz is None:
            return sse_dephasing_step(psi, Gamma_map, dt, rng=rng)
        kernel = build_gaussian_kernel(dr, dz, Nz, Nr, xi)

    # White noise then filter to impose spatial correlation
    white = rng.standard_normal(size=psi.shape)
    xi_corr = _fft_convolve2d(white, kernel)

    # Normalize to unit variance to ensure correct scaling
    std = np.std(xi_corr)
    if std > 0:
        xi_corr = xi_corr / std

    phase_std = np.sqrt(np.maximum(0.0, 2.0 * Gamma_map * dt))
    phase = phase_std * xi_corr

    damp = np.exp(-0.5 * Gamma_map * dt)
    psi_new = psi * damp * np.exp(1j * phase)
    return psi_new
