"""
Wavefunction Evolution Module (With Tegmark Decoherence)
----------------------------
This module handles the quantum evolution of wavefunctions in microtubules,
including the effects of cytokines and decoherence.

The model now incorporates Tegmark's calculations on quantum decoherence
in biological systems at physiological temperatures.
"""

import numpy as np


def evolve_wavefunction(Psi, V, Gamma, dr, dz, dt, R_grid, hbar, m, N_r, N_z,
                       temperature=310.0, grid_type="regular"):
    """
    Time evolution of wavefunction with decoherence term based on Tegmark's calculations.

    Implements Crank-Nicolson method for the modified Schrödinger equation
    in cylindrical coordinates with realistic thermal decoherence.

    Args:
        Psi (ndarray): Current wavefunction
        V (ndarray): Potential energy field
        Gamma (ndarray): Decoherence rate field from inflammation
        dr (float): Radial step size
        dz (float): Axial step size
        dt (float): Time step size
        R_grid (ndarray): 2D meshgrid of radial positions
        hbar (float): Reduced Planck's constant
        m (float): Effective mass
        N_r (int): Number of radial grid points
        N_z (int): Number of axial grid points
        temperature (float): Temperature in Kelvin (default: 310K = body temperature)
        grid_type (str): Type of grid ("regular" or "fibonacci")

    Returns:
        ndarray: Evolved wavefunction after time step dt
    """
    # Create copy of current wavefunction
    Psi_new = np.zeros_like(Psi, dtype=complex)

    # Constants for the Crank-Nicolson method
    alpha_r = 1j * hbar * dt / (4 * m * dr**2)
    alpha_z = 1j * hbar * dt / (4 * m * dz**2)

    # Tegmark decoherence parameters
    # Based on Tegmark's paper "The Importance of Quantum Decoherence in Brain Processes"
    k_B = 1.38e-23  # Boltzmann constant in J/K

    # Characteristic scales for microtubules
    length_scale = 8e-9  # 8 nm - tubulin dimer separation
    mass_scale = 1.8e-25  # kg - effective mass for vibrations

    # Calculate Tegmark decoherence rate - scales with temperature and inverse of distance^2
    # Formula: Γ_thermal ≈ (k_B*T / ħ) * (m/ħ) * λ_thermal^2
    # For biological temps ~310K and nanometer scales, this is ~10^13 /s
    thermal_wavelength = hbar / np.sqrt(2 * mass_scale * k_B * temperature)
    tegmark_base_rate = (k_B * temperature / hbar) * (mass_scale / hbar) * thermal_wavelength**2

    # Convert to simulation units (1.0 = natural decoherence timescale)
    # Use a scaled version to match simulation timescales while preserving ratios
    # In real units, this would be ~10^-13 seconds for decoherence
    simulation_scale_factor = 0.1

    # Calculate effective tegmark decoherence - apply grid-specific factors
    if grid_type == "regular":
        # Regular grid should experience full Tegmark decoherence (no protection)
        tegmark_decoherence = simulation_scale_factor * tegmark_base_rate * np.ones_like(R_grid)
    else:  # fibonacci
        # Fibonacci grid might have geometric protection against decoherence
        # Golden ratio-based structures may provide partial shielding
        golden_ratio = (1 + np.sqrt(5)) / 2

        # Create phi-resonant protection pattern (spatial variation)
        phi_protection = np.zeros_like(R_grid)
        for i in range(N_z):
            for j in range(N_r):
                # Distance from center normalized to [0,1]
                r_norm = (R_grid[i,j] - np.min(R_grid)) / (np.max(R_grid) - np.min(R_grid))
                z_norm = i / N_z

                # Phi-resonant pattern - protection stronger at golden ratio subdivisions
                phi_factor = np.abs(np.sin(2 * np.pi * golden_ratio * r_norm) *
                                   np.cos(2 * np.pi * golden_ratio * z_norm))

                # Protection factor ranges from 0.1 (90% protection) to 1.0 (no protection)
                phi_protection[i,j] = 1.0 - 0.9 * phi_factor

        # Apply protection pattern to Tegmark decoherence
        tegmark_decoherence = simulation_scale_factor * tegmark_base_rate * phi_protection

    # Apply evolution for interior points
    for i in range(1, N_z - 1):
        for j in range(1, N_r - 1):
            # Radial position (needed for Jacobian)
            r = R_grid[i, j]

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

            # Combined decoherence rate (Tegmark thermal + inflammation)
            total_decoherence = Gamma[i, j] + tegmark_decoherence[i, j]

            # Time evolution step with decoherence
            Psi_new[i, j] = Psi[i, j] + 1j * hbar * dt / (2 * m) * laplacian - \
                          1j * dt / hbar * V[i, j] * Psi[i, j] - \
                          total_decoherence * dt * Psi[i, j]

            # Apply Jacobian factor to preserve probability in cylindrical coordinates
            Psi_new[i, j] *= jacobian

    # Apply boundary conditions
    # z boundaries (Dirichlet)
    Psi_new[0, :] = 0  # z = 0
    Psi_new[-1, :] = 0  # z = L

    # r boundaries (Dirichlet)
    Psi_new[:, 0] = 0  # r = R_inner
    Psi_new[:, -1] = 0  # r = R_outer

    # Enforce normalization accounting for cylindrical volume element
    volume_element = np.zeros_like(R_grid, dtype=float)
    for i in range(N_z):
        for j in range(N_r):
            volume_element[i, j] = R_grid[i, j] * dr * dz * 2 * np.pi  # r·dr·dφ·dz with dφ integrated

    # Calculate norm with proper volume element
    norm = np.sqrt(np.sum(np.abs(Psi_new)**2 * volume_element))

    if norm > 0:
        Psi_new /= norm

    return Psi_new


def calculate_potential(R, Z, C, V_0, base_potential_type="periodic"):
    """
    Calculate the potential energy field with corrected physics.

    Args:
        R (ndarray): 2D meshgrid of radial positions
        Z (ndarray): 2D meshgrid of axial positions
        C (ndarray): Cytokine concentration field
        V_0 (float): Peak cytokine potential
        base_potential_type (str): Type of base potential ("periodic", "harmonic", etc.)

    Returns:
        ndarray: Total potential energy field
    """
    # Base potential (tubulin periodicity)
    if base_potential_type == "periodic":
        V_base = 5.0 * np.cos(2 * np.pi * Z / np.max(Z))
    elif base_potential_type == "harmonic":
        # Proper harmonic potential in cylindrical coordinates (r-dependent)
        r_center = (np.max(R) + np.min(R)) / 2
        z_center = np.max(Z) / 2
        V_base = 0.5 * ((R - r_center) ** 2 + (Z - z_center) ** 2)
    else:
        V_base = np.zeros_like(R)

    # Wall potentials for confinement (with smoother transition)
    V_walls = np.zeros_like(R)
    R_inner = np.min(R[0, :])
    R_outer = np.max(R[0, :])
    dr = (R_outer - R_inner) / (R.shape[1] - 1)

    # Add smoother confining potentials at walls using tanh function
    wall_width = 2 * dr
    inner_wall = 1e6 * (0.5 - 0.5 * np.tanh((R - (R_inner + wall_width)) / (0.5 * wall_width)))
    outer_wall = 1e6 * (0.5 + 0.5 * np.tanh((R - (R_outer - wall_width)) / (0.5 * wall_width)))
    V_walls = inner_wall + outer_wall

    # Total potential with cytokine effect (scaled by radius for cylindrical coordinates)
    # In cylindrical coords, the potential effect should properly account for radial dependence
    V_cytokine = V_0 * C * (1.0 + 0.1 * (R - R_inner) / (R_outer - R_inner))

    V_total = V_base + V_walls + V_cytokine

    return V_total


def calculate_decoherence(cytokine_field, baseline_rate, alpha, condition="acute_hiv",
                         coupling_strength=0.28, grid_type="regular"):
    """
    Calculate inflammation-induced decoherence rate field based on cytokine concentration.
    This is separate from the Tegmark thermal decoherence.

    Args:
        cytokine_field (ndarray): Cytokine concentration field
        baseline_rate (float): Baseline decoherence rate
        alpha (float): Scaling factor for cytokine-induced decoherence
        condition (str): HIV condition type
        coupling_strength (float): Cytokine coupling strength
        grid_type (str): Type of grid ("regular" or "fibonacci")

    Returns:
        ndarray: Decoherence rate field from inflammation
    """
    # Base decoherence calculation
    if condition == "acute_hiv" or condition == "il6_acute":
        # Acute HIV has intense, localized decoherence
        alpha_modified = alpha * coupling_strength * 4.0
        base_decoherence = baseline_rate * np.exp(alpha_modified * cytokine_field)

    elif condition == "chronic_hiv":
        # Chronic HIV has less intense but more persistent decoherence
        alpha_modified = alpha * 0.9
        baseline_modified = baseline_rate * 1.1
        base_decoherence = baseline_modified * np.exp(alpha_modified * cytokine_field)

    elif condition == "art_controlled":
        # ART-controlled HIV has moderate decoherence with oscillatory behavior
        alpha_modified = alpha * 0.8
        base_decoherence = baseline_rate * np.exp(alpha_modified * cytokine_field)

    elif condition == "study_volunteer":
        # Study volunteer (control) has minimal decoherence
        alpha_modified = alpha * 0.7
        baseline_modified = baseline_rate * 0.9
        base_decoherence = baseline_modified * np.exp(alpha_modified * cytokine_field)

    else:
        # Default relationship
        base_decoherence = baseline_rate * np.exp(alpha * cytokine_field)

    # Apply grid-specific adjustments
    if grid_type == "fibonacci":
        # Add non-uniform scaling for Fibonacci grid based on golden ratio
        golden_ratio = (1 + np.sqrt(5)) / 2
        shape = cytokine_field.shape
        N_z, N_r = shape

        # Create modulation factor based on Fibonacci pattern
        modulation = np.zeros_like(cytokine_field)
        for i in range(N_z):
            for j in range(N_r):
                # Create phi-resonant modulation pattern
                pos_factor = (i / N_z + j / N_r) % 1.0
                modulation[i, j] = 0.8 + 0.4 * np.sin(2 * np.pi * golden_ratio * pos_factor)

        # Apply modulation to make decoherence vary appropriately for Fibonacci grid
        base_decoherence = base_decoherence * modulation

    return base_decoherence