"""
Physics Metrics Module (Corrected)
---------------------
This module calculates various quantum metrics for analyzing
wavefunction evolution in microtubules, including coherence,
variance, entropy, and event horizon detection.
"""

import numpy as np
from scipy import integrate


def calculate_probability_density(wavefunction, R_grid=None):
    """
    Calculate probability density from wavefunction,
    accounting for cylindrical coordinates.

    Args:
        wavefunction (ndarray): Complex wavefunction
        R_grid (ndarray, optional): 2D meshgrid of radial positions
                                   for Jacobian correction

    Returns:
        ndarray: Probability density
    """
    # Basic probability density
    prob = np.abs(wavefunction) ** 2

    # Apply Jacobian correction if R_grid is provided
    if R_grid is not None:
        # In cylindrical coordinates, probability density includes r factor
        # from the volume element r·dr·dφ·dz
        prob = prob * R_grid

    return prob


def calculate_coherence(psi_current, psi_initial, dr, dz, R_grid):
    """
    Calculate quantum coherence as correlation with initial state,
    with proper normalization in cylindrical coordinates.

    Args:
        psi_current (ndarray): Current probability density
        psi_initial (ndarray): Initial probability density
        dr (float): Radial step size
        dz (float): Axial step size
        R_grid (ndarray): 2D meshgrid of radial positions for Jacobian

    Returns:
        float: Coherence measure (1.0 = perfect coherence)
    """
    # Create volume element
    volume_element = R_grid * dr * dz * 2 * np.pi  # r·dr·dφ·dz with dφ integrated

    # Calculate overlap between current and initial states with proper volume element
    coherence = np.sum(np.sqrt(psi_initial * psi_current) * volume_element)

    # Normalize coherence to range [0, 1]
    initial_norm = np.sqrt(np.sum(psi_initial * volume_element))
    current_norm = np.sqrt(np.sum(psi_current * volume_element))

    if initial_norm > 0 and current_norm > 0:
        coherence /= (initial_norm * current_norm)

    return coherence


def calculate_dispersion_metrics(wavefunction, R_grid, Z_grid, dr, dz):
    """
    Calculate comprehensive dispersion metrics for wavefunction,
    properly accounting for cylindrical coordinates.

    Args:
        wavefunction (ndarray): Complex wavefunction
        R_grid (ndarray): 2D meshgrid of radial positions
        Z_grid (ndarray): 2D meshgrid of axial positions
        dr (float): Radial step size
        dz (float): Axial step size

    Returns:
        dict: Dictionary of dispersion metrics
    """
    # Create volume element
    volume_element = R_grid * dr * dz * 2 * np.pi  # r·dr·dφ·dz with dφ integrated

    # Calculate probability density with Jacobian correction
    prob = np.abs(wavefunction) ** 2

    # Normalize probability density with volume element
    total_prob = np.sum(prob * volume_element)
    if total_prob > 0:
        prob = prob / total_prob

    # 1. Calculate expectation values (means)
    r_mean = np.sum(R_grid * prob * volume_element)
    z_mean = np.sum(Z_grid * prob * volume_element)

    # 2. Calculate variance (standard dispersion measure)
    var_r = np.sum((R_grid - r_mean) ** 2 * prob * volume_element)
    var_z = np.sum((Z_grid - z_mean) ** 2 * prob * volume_element)
    total_variance = var_r + var_z

    # 3. Calculate entropy (information-theoretic measure of dispersion)
    # Add small epsilon to avoid log(0)
    epsilon = 1e-10
    # Correct entropy calculation for continuous probability distribution
    entropy = -np.sum(prob * np.log(prob + epsilon) * volume_element)

    # 4. Calculate kurtosis (measure of the "tailedness" of the distribution)
    # Higher values indicate more extreme deviations
    if var_r > 0:
        r_kurtosis = np.sum((R_grid - r_mean) ** 4 * prob * volume_element) / (var_r ** 2)
    else:
        r_kurtosis = 0

    if var_z > 0:
        z_kurtosis = np.sum((Z_grid - z_mean) ** 4 * prob * volume_element) / (var_z ** 2)
    else:
        z_kurtosis = 0

    total_kurtosis = r_kurtosis + z_kurtosis

    return {
        'variance': total_variance,
        'entropy': entropy,
        'kurtosis': total_kurtosis,
        'var_r': var_r,
        'var_z': var_z,
        'r_mean': r_mean,
        'z_mean': z_mean
    }


def calculate_event_horizon(psi, r, coherence_level, grid_type="regular", R_grid=None,
                            min_threshold=0.05, max_threshold=0.2):
    """
    Calculate event horizon boundary based on probability density threshold,
    accounting for proper cylindrical coordinates.

    Args:
        psi (ndarray): Probability density
        r (ndarray): 1D array of radial positions
        coherence_level (float): Current coherence level (used for adaptive threshold)
        grid_type (str): Grid type ("regular" or "fibonacci")
        R_grid (ndarray, optional): 2D meshgrid of radial positions for Jacobian
        min_threshold (float): Minimum threshold ratio
        max_threshold (float): Maximum threshold ratio

    Returns:
        ndarray: Radial positions of event horizon for each axial position
    """
    N_z, N_r = psi.shape

    # Apply Jacobian correction if R_grid is provided
    if R_grid is not None:
        # Correct probability density for cylindrical coordinates
        psi_corrected = psi * R_grid
    else:
        psi_corrected = psi

    current_max_prob = np.max(psi_corrected)

    # Set different thresholds for different grid types to ensure meaningful comparison
    if grid_type == "regular":
        threshold_ratio = max(min_threshold, min(max_threshold, coherence_level * 0.22))
    else:  # fibonacci
        # Use same threshold calculation for both grid types - corrected to ensure proper comparison
        threshold_ratio = max(min_threshold, min(max_threshold, coherence_level * 0.22))

    threshold = threshold_ratio * current_max_prob

    # Calculate horizon for each z position
    horizon_r_indices = np.zeros(N_z, dtype=int)

    for i in range(N_z):
        # Find where probability drops below threshold
        # Moving from outer boundary inward
        for j in range(N_r - 1, 0, -1):
            if psi_corrected[i, j] > threshold:
                horizon_r_indices[i] = j
                break

    # Convert indices to radial positions
    return r[horizon_r_indices]


def calculate_integrated_metrics(time_array, coherence_array):
    """
    Calculate integrated metrics over time.

    Args:
        time_array (ndarray): Array of time points
        coherence_array (ndarray): Array of coherence values

    Returns:
        dict: Dictionary of integrated metrics
    """
    # Check for valid inputs
    if len(time_array) != len(coherence_array) or len(time_array) == 0:
        return {
            'integrated_coherence': 0.0,
            'half_life': 0.0,
            'initial_coherence': 0.0 if len(coherence_array) == 0 else coherence_array[0],
            'final_coherence': 0.0 if len(coherence_array) == 0 else coherence_array[-1]
        }

    # Area under the coherence curve
    integrated_coherence = integrate.trapezoid(coherence_array, time_array)

    # Calculate coherence half-life
    initial_coherence = coherence_array[0]
    final_coherence = coherence_array[-1]
    half_decay_value = initial_coherence - (initial_coherence - final_coherence) / 2

    # Find the time at which coherence drops to half-decay value
    half_life = time_array[-1]  # Default to final time if never reaches half
    for i, coh in enumerate(coherence_array):
        if coh <= half_decay_value:
            if i > 0:
                # Linear interpolation for more precise half-life
                t1, t2 = time_array[i - 1], time_array[i]
                c1, c2 = coherence_array[i - 1], coherence_array[i]
                half_life = t1 + (half_decay_value - c1) * (t2 - t1) / (c2 - c1)
            else:
                half_life = time_array[i]
            break

    return {
        'integrated_coherence': integrated_coherence,
        'half_life': half_life,
        'initial_coherence': initial_coherence,
        'final_coherence': final_coherence
    }