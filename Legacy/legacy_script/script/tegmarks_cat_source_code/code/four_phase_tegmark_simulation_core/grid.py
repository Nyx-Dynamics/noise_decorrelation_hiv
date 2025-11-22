"""
Grid Generation Module (Corrected)
---------------------
This module handles the generation of regular and Fibonacci-scaled grids
for microtubule quantum simulations.

The Fibonacci scaling is crucial for modeling event horizon-like boundaries
that form in microtubules in response to perturbations.
"""

import numpy as np


def generate_regular_grid(R_inner, R_outer, L, N_r, N_z):
    """
    Generate a regular grid in cylindrical coordinates.

    Args:
        R_inner (float): Inner radius of microtubule
        R_outer (float): Outer radius of microtubule
        L (float): Axial length of microtubule
        N_r (int): Number of radial grid points
        N_z (int): Number of axial grid points

    Returns:
        tuple: (r, z, R, Z, dr, dz) where:
            r (ndarray): 1D array of radial positions
            z (ndarray): 1D array of axial positions
            R (ndarray): 2D meshgrid of radial positions
            Z (ndarray): 2D meshgrid of axial positions
            dr (float): Radial step size
            dz (float): Axial step size
    """
    # Create 1D grids
    r = np.linspace(R_inner, R_outer, N_r)
    z = np.linspace(0, L, N_z)

    # Calculate step sizes
    dr = (R_outer - R_inner) / (N_r - 1)
    dz = L / (N_z - 1)

    # Create 2D meshgrid
    R, Z = np.meshgrid(r, z)

    return r, z, R, Z, dr, dz


def generate_fibonacci_sequence(n):
    """
    Generate Fibonacci sequence up to n terms.

    Args:
        n (int): Number of terms to generate

    Returns:
        ndarray: Fibonacci sequence
    """
    # Initialize with first two Fibonacci numbers
    fib = [0, 1]

    # Generate remaining terms
    for i in range(2, n):
        fib.append(fib[i - 1] + fib[i - 2])

    return np.array(fib)


def generate_fibonacci_grid(R_inner, R_outer, L, N_r, N_z):
    """
    Generate a Fibonacci-scaled grid in cylindrical coordinates.

    The Fibonacci scaling provides non-uniform spacing that follows
    the Fibonacci sequence pattern, which has been hypothesized to
    better preserve quantum coherence in biological systems.

    Args:
        R_inner (float): Inner radius of microtubule
        R_outer (float): Outer radius of microtubule
        L (float): Axial length of microtubule
        N_r (int): Number of radial grid points
        N_z (int): Number of axial grid points

    Returns:
        tuple: (r_fib, z_fib, R_fib, Z_fib, dr, dz) where:
            r_fib (ndarray): 1D array of Fibonacci-scaled radial positions
            z_fib (ndarray): 1D array of Fibonacci-scaled axial positions
            R_fib (ndarray): 2D meshgrid of Fibonacci-scaled radial positions
            Z_fib (ndarray): 2D meshgrid of Fibonacci-scaled axial positions
            dr (float): Effective average radial step size (for normalization)
            dz (float): Effective average axial step size (for normalization)
    """
    # Get regular grid for comparison and step sizes
    r, z, _, _, dr, dz = generate_regular_grid(R_inner, R_outer, L, N_r, N_z)

    # Generate Fibonacci sequence with padding to ensure we have enough elements
    fib_seq = generate_fibonacci_sequence(max(N_r, N_z) + 10)

    # Remove the leading zeros to avoid singularities in the grid
    fib_seq = fib_seq[fib_seq > 0]

    # Extract the appropriate number of elements
    if len(fib_seq) < N_r:
        # Pad with additional Fibonacci numbers if needed
        while len(fib_seq) < N_r:
            fib_seq = np.append(fib_seq, fib_seq[-1] + fib_seq[-2])

    fib_r = fib_seq[-(N_r):]  # Use last N_r elements

    # Same for z-direction
    if len(fib_seq) < N_z:
        # Pad with additional Fibonacci numbers if needed
        while len(fib_seq) < N_z:
            fib_seq = np.append(fib_seq, fib_seq[-1] + fib_seq[-2])

    fib_z = fib_seq[-(N_z):]  # Use last N_z elements

    # Normalize sequences to physical dimensions with improved scaling
    # Ensure min value is not 0 to avoid singularities in cylindrical coordinates
    min_fib_r = np.min(fib_r)
    max_fib_r = np.max(fib_r)
    fib_r_scaled = R_inner + (fib_r - min_fib_r) / (max_fib_r - min_fib_r) * (R_outer - R_inner)

    min_fib_z = np.min(fib_z)
    max_fib_z = np.max(fib_z)
    fib_z_scaled = (fib_z - min_fib_z) / (max_fib_z - min_fib_z) * L

    # Create 2D Fibonacci grid with correct ordering
    R_fib = np.zeros((N_z, N_r))
    Z_fib = np.zeros((N_z, N_r))

    for i in range(N_z):
        for j in range(N_r):
            R_fib[i, j] = fib_r_scaled[j]
            Z_fib[i, j] = fib_z_scaled[i]

    return fib_r_scaled, fib_z_scaled, R_fib, Z_fib, dr, dz


def initialize_wavefunction(grid_type, R, Z, R_fib, Z_fib, R_inner, R_outer, L, dr, dz):
    """
    Initialize the wavefunction on regular or Fibonacci grid with proper normalization
    in cylindrical coordinates.

    Args:
        grid_type (str): Either 'regular' or 'fibonacci'
        R (ndarray): 2D meshgrid of regular radial positions
        Z (ndarray): 2D meshgrid of regular axial positions
        R_fib (ndarray): 2D meshgrid of Fibonacci-scaled radial positions
        Z_fib (ndarray): 2D meshgrid of Fibonacci-scaled axial positions
        R_inner (float): Inner radius of microtubule
        R_outer (float): Outer radius of microtubule
        L (float): Axial length of microtubule
        dr (float): Radial step size
        dz (float): Axial step size

    Returns:
        ndarray: Initial wavefunction
    """
    # Common parameters
    sigma_z = L / 10  # Width of initial Gaussian
    z0 = L / 2  # Center position

    if grid_type == 'regular':
        # Initialize on regular grid
        Psi = np.exp(-0.5 * ((Z - z0) / sigma_z) ** 2) * np.sin(np.pi * (R - R_inner) / (R_outer - R_inner))

        # Create volume element for cylindrical coordinates
        volume_element = R * dr * dz * 2 * np.pi  # r·dr·dφ·dz with dφ integrated

    elif grid_type == 'fibonacci':
        # Initialize on Fibonacci grid
        Psi = np.exp(-0.5 * ((Z_fib - z0) / sigma_z) ** 2) * np.sin(np.pi * (R_fib - R_inner) / (R_outer - R_inner))

        # Create volume element for cylindrical coordinates
        volume_element = R_fib * dr * dz * 2 * np.pi  # Using average dr for simplicity
    else:
        raise ValueError(f"Unknown grid type: {grid_type}")

    # Normalize wavefunction with proper volume element
    norm = np.sqrt(np.sum(np.abs(Psi) ** 2 * volume_element))

    if norm > 0:
        Psi /= norm

    return Psi