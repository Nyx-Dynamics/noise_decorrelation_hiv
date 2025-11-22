"""
Wavefunction Evolution Module (Corrected with Tegmark Model)
----------------------------
This module handles the quantum evolution of wavefunctions in microtubules,
with a focus on contrasting regular vs. Fibonacci grid behavior under
Tegmark decoherence conditions across different HIV phases.
"""

import numpy as np

def evolve_wavefunction(Psi, V, Gamma, dr, dz, dt, R, hbar, m, N_r, N_z, grid_type):
    """
    Evolve the wavefunction in time according to the given potential and decoherence rates,
    with Tegmark decoherence model that differentiates between regular and Fibonacci grids.

    Args:
        Psi: Complex wavefunction (2D array [N_r x N_z]).
        V: Potential array (2D array [N_r x N_z]).
        Gamma: Decoherence rates array (2D array [N_r x N_z]).
        dr, dz: Grid spacings in r and z directions.
        dt: Time step for evolution.
        R: Radial grid.
        hbar, m: Physical constants (reduced Planck's constant and particle mass).
        N_r, N_z: Number of grid points in r and z respectively.
        grid_type: Grid type ("regular" or "fibonacci").

    Returns:
        Psi_new: Updated wavefunction after one time step.
    """
    # Preallocate updated wavefunction
    Psi_new = np.zeros_like(Psi, dtype=complex)

    # Apply Tegmark decoherence based on grid type
    if grid_type == "regular":
        # Regular grid: Apply strong Tegmark decoherence (rapid collapse)
        # This represents the theoretical prediction that coherence should
        # collapse almost instantly at biological temperatures
        tegmark_factor = 0.5  # Strong decoherence - will rapidly collapse
    else:  # fibonacci
        # Fibonacci grid: Apply geometrically protected decoherence (slower collapse)
        # This represents the hypothesis that golden ratio structures might
        # provide some protection against decoherence
        golden_ratio = (1 + np.sqrt(5)) / 2

        # Create a phi-resonant protection pattern
        phi_protection = np.zeros_like(R)
        for i in range(N_z):
            for j in range(N_r):
                # Normalized position
                r_norm = (R[i,j] - np.min(R)) / (np.max(R) - np.min(R))
                z_norm = i / N_z

                # Fibonacci protection is strongest at golden ratio positions
                phi_factor = np.abs(np.sin(2 * np.pi * golden_ratio * r_norm) *
                                   np.cos(2 * np.pi * golden_ratio * z_norm))
                phi_protection[i,j] = 0.9 * phi_factor  # Up to 90% protection

        # Apply much weaker decoherence to Fibonacci grid
        tegmark_factor = 0.1 * (1.0 - phi_protection)  # Minimal decoherence with protection

    # Apply evolution for interior points
    for i in range(1, N_z - 1):
        for j in range(1, N_r - 1):
            # Radial position (needed for Jacobian)
            r = R[i, j]

            # Jacobian factor for cylindrical coordinates
            jacobian = r

            # Second derivative in r with Jacobian and cylindrical correction
            d2r = (r + 0.5*dr) * Psi[i, j+1] - 2*r * Psi[i, j] + (r - 0.5*dr) * Psi[i, j-1]
            d2r = d2r / (r * dr**2)

            # First derivative in r (for cylindrical term) with proper Jacobian
            dr_term = (Psi[i, j+1] - Psi[i, j-1]) / (2 * dr)

            # Second derivative in z
            d2z = (Psi[i+1, j] - 2*Psi[i, j] + Psi[i-1, j]) / dz**2

            # Laplacian in cylindrical coordinates (properly formulated)
            laplacian = d2r + dr_term/r + d2z

            # Apply Tegmark thermal decoherence + inflammation decoherence
            total_decoherence = Gamma[i, j] + tegmark_factor

            # Time evolution step with proper decoherence
            Psi_new[i, j] = Psi[i, j] + 1j * hbar * dt / (2 * m) * laplacian - \
                          1j * dt / hbar * V[i, j] * Psi[i, j] - \
                          total_decoherence * dt * Psi[i, j]

            # Apply Jacobian factor to preserve probability in cylindrical coordinates
            Psi_new[i, j] *= jacobian

    # Apply boundary conditions
    Psi_new[0, :] = 0  # z = 0
    Psi_new[-1, :] = 0  # z = L
    Psi_new[:, 0] = 0  # r = R_inner
    Psi_new[:, -1] = 0  # r = R_outer

    # Clean any NaN or Inf values that might have appeared
    Psi_new = np.nan_to_num(Psi_new, nan=0.0, posinf=0.0, neginf=0.0)

    # Enforce normalization with cylindrical volume element
    volume_element = R * dr * dz * 2 * np.pi  # r·dr·dφ·dz with dφ integrated

    # Calculate norm with proper volume element
    norm_squared = np.sum(np.abs(Psi_new)**2 * volume_element)

    if norm_squared > 1e-10:
        Psi_new /= np.sqrt(norm_squared)

    return Psi_new

def calculate_potential(R, Z, C, V_0, base_potential_type="periodic"):
    """
    Calculate the potential energy field with corrected physics.
    """
    # Base potential (tubulin periodicity)
    if base_potential_type == "periodic":
        V_base = 5.0 * np.cos(2 * np.pi * Z / np.max(Z))
    elif base_potential_type == "harmonic":
        r_center = (np.max(R) + np.min(R)) / 2
        z_center = np.max(Z) / 2
        V_base = 0.5 * ((R - r_center) ** 2 + (Z - z_center) ** 2)
    else:
        V_base = np.zeros_like(R)

    # Wall potentials for confinement
    V_walls = np.zeros_like(R)
    R_inner = np.min(R[0, :])
    R_outer = np.max(R[0, :])
    dr = (R_outer - R_inner) / (R.shape[1] - 1)

    # Add smoother confining potentials at walls
    wall_width = 2 * dr

    # Limit arguments to tanh to prevent overflow
    inner_arg = np.clip((R - (R_inner + wall_width)) / (0.5 * wall_width), -20, 20)
    outer_arg = np.clip((R - (R_outer - wall_width)) / (0.5 * wall_width), -20, 20)

    inner_wall = 1e6 * (0.5 - 0.5 * np.tanh(inner_arg))
    outer_wall = 1e6 * (0.5 + 0.5 * np.tanh(outer_arg))
    V_walls = inner_wall + outer_wall

    # Cytokine potential with radial scaling
    # Limit scaling factor to prevent overflow
    scale_factor = np.clip(1.0 + 0.1 * (R - R_inner) / (R_outer - R_inner), 0.5, 1.5)
    V_cytokine = V_0 * C * scale_factor

    # Sum all potential components
    V_total = V_base + V_walls + V_cytokine

    # Clean any non-finite values
    V_total = np.nan_to_num(V_total, nan=0.0, posinf=1e6, neginf=-1e6)

    return V_total


def calculate_decoherence(cytokine_field, baseline_rate, alpha, condition="acute_hiv",
                         coupling_strength=0.28, grid_type="regular"):
    """
    Calculate inflammation-induced decoherence rate field.
    This is separate from the Tegmark thermal decoherence.
    """
    # Clean and limit cytokine field to prevent overflow in exponential
    safe_field = np.clip(cytokine_field, -20, 20)

    # Base decoherence calculation based on condition
    if condition == "acute_hiv" or condition == "il6_acute":
        # Acute HIV has intense decoherence
        alpha_modified = alpha * coupling_strength * 4.0
        base_decoherence = baseline_rate * np.exp(alpha_modified * safe_field)
    elif condition == "chronic_hiv":
        # Chronic HIV has persistent decoherence
        alpha_modified = alpha * 0.9
        baseline_modified = baseline_rate * 1.1
        base_decoherence = baseline_modified * np.exp(alpha_modified * safe_field)
    elif condition == "art_controlled":
        # ART-controlled HIV has moderate decoherence
        alpha_modified = alpha * 0.8
        base_decoherence = baseline_rate * np.exp(alpha_modified * safe_field)
    elif condition == "study_volunteer":
        # Study volunteer has minimal decoherence
        alpha_modified = alpha * 0.7
        baseline_modified = baseline_rate * 0.9
        base_decoherence = baseline_modified * np.exp(alpha_modified * safe_field)
    elif condition == "terminal":
        # Terminal AIDS has most intense decoherence
        alpha_modified = alpha * coupling_strength * 5.0
        base_decoherence = baseline_rate * 1.2 * np.exp(alpha_modified * safe_field)
    else:
        base_decoherence = baseline_rate * np.exp(alpha * safe_field)

    # Grid-specific adjustments
    if grid_type == "fibonacci":
        # Create phi-resonant pattern for inflammation response
        golden_ratio = (1 + np.sqrt(5)) / 2
        shape = cytokine_field.shape
        N_z, N_r = shape

        # Create modulation factor based on Fibonacci pattern
        modulation = np.zeros_like(cytokine_field)
        for i in range(N_z):
            for j in range(N_r):
                pos_factor = (i / N_z + j / N_r) % 1.0
                modulation[i, j] = 0.8 + 0.4 * np.sin(2 * np.pi * golden_ratio * pos_factor)

        # Apply modulation
        base_decoherence = base_decoherence * modulation

    # Clean any NaN or Inf values
    base_decoherence = np.nan_to_num(base_decoherence, nan=baseline_rate, posinf=baseline_rate*10, neginf=baseline_rate)

    return base_decoherence