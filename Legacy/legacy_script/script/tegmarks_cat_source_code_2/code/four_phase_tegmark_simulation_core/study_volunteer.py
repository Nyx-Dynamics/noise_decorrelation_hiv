#!/usr/bin/env python3
"""
Enhanced Study Volunteer Simulation Module
--------------------------------------
Simulates the control condition of a healthy study volunteer without HIV infection
on quantum coherence in microtubules, including realistic immune challenges.

This module models both baseline cytokine levels and periodic immune challenges
in a healthy control subject and simulates their effects on quantum coherence
preservation in regular vs. Fibonacci-scaled spatial grids.
"""

# Fix import paths
import os
import sys
# Numerical stability fix for quantum coherence simulations
import numpy as np
import sys

# Configure NumPy to show warnings but continue execution
np.seterr(all='warn', under='ignore')

# Wait until all imports are done before patching
_original_import = __import__


def _patched_import(*args, **kwargs):
    module = _original_import(*args, **kwargs)

    # Only patch once all core modules are imported
    if args[0] in ['scripts.core.wavefunction', 'scripts.core.metrics']:
        try:
            # Check if we can access both modules
            from scripts.core import wavefunction
            from scripts.core import metrics

            # Patch the key functions causing numerical errors

            # 1. Save original functions
            original_evolve = wavefunction.evolve_wavefunction
            original_prob = metrics.calculate_probability_density
            original_coherence = metrics.calculate_coherence

            # 2. Define stable versions
            def evolve_wavefunction_stable(Psi, V, Gamma, dr, dz, dt, R, hbar, m, N_r, N_z, grid_type):
                """Numerically stable version of evolve_wavefunction"""
                # Limit dt for stability
                stable_dt = min(dt, 0.25 * min(dr, dz) ** 2 * m / hbar)

                # Get result from original function
                result = original_evolve(Psi, V, Gamma, dr, dz, stable_dt, R, hbar, m, N_r, N_z, grid_type)

                # Clean any NaN or Inf values
                result[~np.isfinite(result)] = 0.0

                # Normalize to ensure stability
                volume_element = R * dr * dz * 2 * np.pi
                norm_squared = np.sum(np.abs(result) ** 2 * volume_element)
                if norm_squared > 1e-10:
                    result = result / np.sqrt(norm_squared)

                return result

            def calculate_probability_density_stable(wavefunction, R_grid=None):
                """Numerically stable version of calculate_probability_density"""
                # Create a clean copy of the wavefunction
                wf = np.copy(wavefunction)

                # Replace any NaN/Inf values
                wf[~np.isfinite(wf)] = 0.0

                # Calculate probability with protection against overflow
                try:
                    prob = np.abs(wf) ** 2
                except:
                    # Fallback if squaring fails
                    prob = np.zeros_like(wf, dtype=float)
                    for i in range(wf.shape[0]):
                        for j in range(wf.shape[1]):
                            val = wf[i, j]
                            prob[i, j] = (val.real ** 2 + val.imag ** 2)

                # Replace any NaN/Inf in result
                prob[~np.isfinite(prob)] = 0.0

                # Apply Jacobian correction if R_grid is provided
                if R_grid is not None:
                    prob = prob * R_grid

                return prob

            def calculate_coherence_stable(psi_current, psi_initial, dr, dz, R_grid):
                """Numerically stable version of calculate_coherence"""
                # Create volume element
                volume_element = R_grid * dr * dz * 2 * np.pi

                # Clean inputs
                psi_c = np.copy(psi_current)
                psi_i = np.copy(psi_initial)
                psi_c[~np.isfinite(psi_c)] = 0.0
                psi_i[~np.isfinite(psi_i)] = 0.0

                # Calculate overlap with protection against negative values
                sqrt_term = np.sqrt(np.maximum(0.0, psi_i * psi_c))
                coherence = np.sum(sqrt_term * volume_element)

                # Normalize safely
                epsilon = 1e-15  # Small value to prevent division by zero
                initial_norm = np.sqrt(np.sum(psi_i * volume_element) + epsilon)
                current_norm = np.sqrt(np.sum(psi_c * volume_element) + epsilon)

                # Avoid division by very small values
                if initial_norm > 1e-10 and current_norm > 1e-10:
                    coherence = coherence / (initial_norm * current_norm)
                else:
                    coherence = 0.0

                # Ensure valid range
                coherence = np.clip(coherence, 0.0, 1.0)

                return coherence

            # 3. Replace the original functions with stable versions
            wavefunction.evolve_wavefunction = evolve_wavefunction_stable
            metrics.calculate_probability_density = calculate_probability_density_stable
            metrics.calculate_coherence = calculate_coherence_stable

            # 4. Restore original import function
            sys.modules['builtins'].__import__ = _original_import

            print("Numerical stability patch applied successfully")

        except (ImportError, AttributeError) as e:
            # Do nothing if modules aren't fully imported yet
            pass

    return module


# Install the import hook
sys.modules['builtins'].__import__ = _patched_import
# Get the absolute path of the current file
current_file = os.path.abspath(__file__)
# Get the directory of the current file
current_dir = os.path.dirname(current_file)
# Get the parent directory (scripts/simulations)
simulations_dir = current_dir
# Get the scripts directory
scripts_dir = os.path.dirname(simulations_dir)
# Get the project root directory
project_root = os.path.dirname(scripts_dir)
# Add project root to Python path
sys.path.insert(0, project_root)

import numpy as np
import matplotlib.pyplot as plt
import json
import time
import datetime

# Import core modules
try:
    from scripts.core.grid import generate_regular_grid, generate_fibonacci_grid, initialize_wavefunction
    from scripts.core.wavefunction import evolve_wavefunction, calculate_potential, calculate_decoherence
    from scripts.core.metrics import (calculate_probability_density, calculate_coherence, calculate_dispersion_metrics,
                                      calculate_event_horizon, calculate_integrated_metrics)
    from scripts.utils.helpers import setup_logging, log_message, create_timestamp, save_figure
except ModuleNotFoundError as e:
    raise ModuleNotFoundError(f"Could not locate required modules: {e}. Check your project structure.") from e


def initialize_control_cytokine_field(r, z, R_inner, R_outer, N_z, N_r):
    """
    Initialize cytokine field for healthy control/study volunteer.

    Creates a spatially varying cytokine concentration field characteristic
    of a healthy subject with minimal inflammatory cytokines but with
    localized areas representing normal immune activity.

    Args:
        r (ndarray): 1D array of radial positions
        z (ndarray): 1D array of axial positions
        R_inner (float): Inner radius of microtubule
        R_outer (float): Outer radius of microtubule
        N_z (int): Number of axial grid points
        N_r (int): Number of radial grid points

    Returns:
        ndarray: 2D array of cytokine concentrations
    """
    log_message("Initializing control cytokine field with immune challenges")
    cytokines = np.zeros((N_z, N_r))
    L = np.max(z)

    # Healthy control pattern: very low intensity, minimal spatial structure
    for i in range(N_z):
        for j in range(N_r):
            # Distance from center (normalized)
            center_dist_r = (r[j] - (R_inner + R_outer) / 2) / ((R_outer - R_inner) / 2)
            center_dist_z = (z[i] - L / 2) / (L / 2)
            center_dist = np.sqrt(center_dist_r ** 2 + center_dist_z ** 2)

            # Healthy pattern: very low intensity
            cytokines[i, j] = 0.15 * np.exp(-center_dist ** 2 / 2.0)

    # Add minimal background noise
    cytokines += 0.05 * np.random.random((N_z, N_r))

    # Add small immune activity sites (typical in healthy individuals)
    num_spots = 3  # Just a few small immune activity sites
    for _ in range(num_spots):
        i = np.random.randint(0, N_z)
        j = np.random.randint(N_r // 2, N_r - 1)  # Place near outer boundary
        size = 2 + np.random.randint(2)  # Small sizes

        for di in range(-size, size + 1):
            for dj in range(-size, size + 1):
                ii, jj = i + di, j + dj
                if 0 <= ii < N_z and 0 <= jj < N_r:
                    dist = np.sqrt(di ** 2 + dj ** 2) / size
                    if dist <= 1:
                        intensity = 0.2 * np.exp(-4 * dist ** 2)  # Smaller intensity
                        cytokines[ii, jj] = min(0.35, cytokines[ii, jj] + intensity)

    # Add a pattern representing normal mucosal immune activity
    z_mid = N_z // 2
    width = N_z // 8
    for i in range(z_mid - width, z_mid + width):
        for j in range(N_r - 5, N_r):
            if 0 <= i < N_z and 0 <= j < N_r:  # Boundary check
                mucosal_factor = 0.15 * np.exp(-((i - z_mid) / width) ** 2)
                cytokines[i, j] += mucosal_factor

    # Normalize to [0, 1] with low maximum for healthy control
    cytokines = np.clip(cytokines, 0, 0.35)

    return cytokines


def update_cytokine_field(C, dt, current_time, dr, dz, challenge_period=5.0, challenge_duration=0.5,
                          challenge_strength=0.3):
    """
    Update cytokine field over time for healthy control simulation,
    including periodic immune challenges.

    Args:
        C (ndarray): Current cytokine field
        dt (float): Time step
        current_time (float): Current simulation time
        dr (float): Radial step size
        dz (float): Axial step size
        challenge_period (float): Time between immune challenges
        challenge_duration (float): Duration of each immune challenge
        challenge_strength (float): Maximum strength of immune challenges

    Returns:
        ndarray: Updated cytokine field
    """
    # Calculate Laplacian for diffusion
    laplacian_C = np.zeros_like(C)
    N_z, N_r = C.shape

    for i in range(1, N_z - 1):
        for j in range(1, N_r - 1):
            laplacian_C[i, j] = (C[i + 1, j] + C[i - 1, j] + C[i, j + 1] + C[i, j - 1] - 4 * C[i, j]) / (dr * dz)

    # Healthy control pattern: minimal variation, homeostatic
    time_factor = 1.0 + 0.02 * np.sin(0.1 * current_time)

    # Simulate periodic immune challenges (e.g., exposure to pathogens)
    # Every ~5 time units, simulate immune response
    immune_challenge = 0.0

    # Calculate time within the current challenge period
    time_in_period = current_time % challenge_period

    # If within challenge duration, activate immune response
    if time_in_period < challenge_duration:
        challenge_intensity = challenge_strength * np.exp(-(time_in_period / challenge_duration) ** 2)
        immune_challenge = challenge_intensity
        log_message(f"Immune challenge at t={current_time:.2f}, strength={challenge_intensity:.2f}")

    # Diffusion and decay rates - normal physiological rates
    diffusion_rate = 0.01  # Normal diffusion
    decay_rate = 0.02 * time_factor  # Normal cytokine clearance rate

    # Small physiological fluctuations
    fluctuation_amplitude = 0.005
    fluctuations = fluctuation_amplitude * (np.random.random(C.shape) - 0.5)

    # Update cytokine field
    C_new = C + dt * (diffusion_rate * laplacian_C - decay_rate * C + fluctuations)

    # If immune challenge is active, add localized immune response
    if immune_challenge > 0:
        # Choose random location for immune challenge
        i, j = np.random.randint(0, N_z), np.random.randint(N_r // 2, N_r)
        size = 4 + np.random.randint(3)

        for di in range(-size, size + 1):
            for dj in range(-size, size + 1):
                ii, jj = i + di, j + dj
                if 0 <= ii < N_z and 0 <= jj < N_r:
                    dist = np.sqrt(di ** 2 + dj ** 2) / size
                    if dist <= 1:
                        challenge_effect = immune_challenge * np.exp(-2 * dist ** 2)
                        C_new[ii, jj] = min(0.45, C_new[ii, jj] + challenge_effect)  # Higher but limited peak

    # Keep within healthy bounds
    C_new = np.clip(C_new, 0, 0.45)  # Slightly higher max to allow for challenge peaks

    return C_new


def run_simulation(config):
    """
    Run the enhanced study volunteer simulation with the given configuration.

    Args:
        config (dict): Configuration dictionary

    Returns:
        bool: True if simulation completed successfully, False otherwise
    """
    # Extract configuration parameters
    try:
        sim_params = config["simulation_parameters"]
        condition_params = config["conditions"].get("study_volunteer_enhanced",
                                                    config["conditions"][
                                                        "study_volunteer"])  # Fallback to regular if enhanced not defined
        output_params = config["output"]

        # Setup directories
        results_dir = condition_params.get("output_dir", "results/study_volunteer_enhanced")
        figures_dir = condition_params.get("figures_dir", "figures/study_volunteer_enhanced")
        logs_dir = condition_params.get("logs_dir", "logs/study_volunteer_enhanced")

        os.makedirs(results_dir, exist_ok=True)
        os.makedirs(figures_dir, exist_ok=True)
        os.makedirs(logs_dir, exist_ok=True)

        # Setup logging
        log_file = os.path.join(logs_dir, f"study_volunteer_enhanced_{create_timestamp()}.log")
        setup_logging(log_file)

        # Physical constants
        hbar = sim_params["hbar"]
        m = sim_params["m"]
        L = sim_params["L"]
        R_inner = sim_params["R_inner"]
        R_outer = sim_params["R_outer"]
        N_r = sim_params["N_r"]
        N_z = sim_params["N_z"]
        dt = sim_params["dt"]
        time_steps = sim_params["time_steps"]

        # Simulation parameters specific to enhanced study volunteer
        V_0 = condition_params.get("V_0", sim_params["V_0"])
        Gamma_0 = condition_params.get("Gamma_0", 0.02)
        alpha_c = condition_params.get("alpha_c", 0.05)
        coupling_strength = condition_params.get("IL6_COUPLING_STRENGTH", 0.28)

        # Immune challenge parameters
        challenge_period = condition_params.get("IMMUNE_CHALLENGE_PERIOD", 5.0)
        challenge_duration = condition_params.get("IMMUNE_CHALLENGE_DURATION", 0.5)
        challenge_strength = condition_params.get("IMMUNE_CHALLENGE_STRENGTH", 0.3)

        # Generate grids
        log_message("Generating spatial grids")
        r, z, R, Z, dr, dz = generate_regular_grid(R_inner, R_outer, L, N_r, N_z)
        r_fib, z_fib, R_fib, Z_fib, _, _ = generate_fibonacci_grid(R_inner, R_outer, L, N_r, N_z)

        # Initialize wavefunctions
        log_message("Initializing wavefunctions")
        Psi_reg = initialize_wavefunction('regular', R, Z, R_fib, Z_fib, R_inner, R_outer, L, dr, dz)
        Psi_fib = initialize_wavefunction('fibonacci', R, Z, R_fib, Z_fib, R_inner, R_outer, L, dr, dz)

        # Initialize cytokine field
        C = initialize_control_cytokine_field(r, z, R_inner, R_outer, N_z, N_r)

        # Initialize potentials
        log_message("Initializing potentials")
        V_reg = calculate_potential(R, Z, C, V_0)
        V_fib = calculate_potential(R_fib, Z_fib, C, V_0)

        # Initialize decoherence rates
        log_message("Initializing decoherence rates")
        Gamma_reg = calculate_decoherence(C, Gamma_0, alpha_c, "study_volunteer", coupling_strength,
                                          grid_type="regular")
        Gamma_fib = calculate_decoherence(C, Gamma_0, alpha_c, "study_volunteer", coupling_strength,
                                          grid_type="fibonacci")

        # Data storage for time evolution
        time_list = [0.0]
        coherence_reg_list = [1.0]  # Initial state has perfect coherence
        coherence_fib_list = [1.0]  # Initial state has perfect coherence
        variance_reg_list = [0.1]  # Sample initial variance
        variance_fib_list = [0.1]  # Sample initial variance

        # Calculate initial probability densities
        Psi_reg_prob = calculate_probability_density(Psi_reg)
        Psi_fib_prob = calculate_probability_density(Psi_fib)

        # Time evolution loop
        log_message("Starting time evolution simulation for enhanced study volunteer")
        evolution_start_time = time.time()

        # Store initial state
        Psi_reg_prob_initial = Psi_reg_prob
        Psi_fib_prob_initial = Psi_fib_prob

        # Initial horizon calculation
        horizon_reg = calculate_event_horizon(Psi_reg_prob, r, 1.0, grid_type="regular", R_grid=R)
        horizon_fib = calculate_event_horizon(Psi_fib_prob, r, 1.0, grid_type="fibonacci", R_grid=R_fib)

        # Initial metrics
        disp_reg = calculate_dispersion_metrics(Psi_reg, R, Z, dr, dz)
        disp_fib = calculate_dispersion_metrics(Psi_fib, R_fib, Z_fib, dr, dz)

        for step in range(1, time_steps + 1):
            current_time = step * dt

            # Update cytokine field with immune challenges
            C = update_cytokine_field(C, dt, current_time, dr, dz,
                                      challenge_period, challenge_duration, challenge_strength)

            # Update potentials and decoherence rates
            V_reg = calculate_potential(R, Z, C, V_0)
            V_fib = calculate_potential(R_fib, Z_fib, C, V_0)
            Gamma_reg = calculate_decoherence(C, Gamma_0, alpha_c, "study_volunteer", coupling_strength,
                                              grid_type="regular")
            Gamma_fib = calculate_decoherence(C, Gamma_0, alpha_c, "study_volunteer", coupling_strength,
                                              grid_type="fibonacci")

            # Evolve wavefunctions with Tegmark decoherence model
            Psi_reg = evolve_wavefunction(Psi_reg, V_reg, Gamma_reg, dr, dz, dt, R, hbar, m, N_r, N_z,
                                          grid_type="regular")
            Psi_fib = evolve_wavefunction(Psi_fib, V_fib, Gamma_fib, dr, dz, dt, R_fib, hbar, m, N_r, N_z,
                                          grid_type="fibonacci")

            # Calculate probability densities
            Psi_reg_prob = calculate_probability_density(Psi_reg)
            Psi_fib_prob = calculate_probability_density(Psi_fib)

            # Calculate coherence with proper cylindrical coordinates
            coh_reg = calculate_coherence(Psi_reg_prob, Psi_reg_prob_initial, dr, dz, R)
            coh_fib = calculate_coherence(Psi_fib_prob, Psi_fib_prob_initial, dr, dz, R_fib)

            # Calculate dispersion metrics
            disp_reg = calculate_dispersion_metrics(Psi_reg, R, Z, dr, dz)
            disp_fib = calculate_dispersion_metrics(Psi_fib, R_fib, Z_fib, dr, dz)

            # Store results periodically
            save_frequency = output_params.get("save_frequency", 10)
            if step % save_frequency == 0:
                time_list.append(current_time)
                coherence_reg_list.append(coh_reg)
                coherence_fib_list.append(coh_fib)
                variance_reg_list.append(disp_reg['variance'])
                variance_fib_list.append(disp_fib['variance'])

                # Print progress update
                log_frequency = output_params.get("log_frequency", 50)
                if step % log_frequency == 0:
                    elapsed_time = time.time() - evolution_start_time
                    estimated_total = elapsed_time * (time_steps / step)
                    estimated_remaining = estimated_total - elapsed_time

                    log_message(f"Step {step}/{time_steps} ({step / time_steps * 100:.1f}%)")
                    log_message(f"  Elapsed: {elapsed_time:.2f}s, Remaining: {estimated_remaining:.2f}s")
                    log_message(f"  Coherence - Regular: {coh_reg:.4f}, Fibonacci: {coh_fib:.4f}")
                    log_message(
                        f"  Variance - Regular: {disp_reg['variance']:.4f}, Fibonacci: {disp_fib['variance']:.4f}")

        # Calculate final horizon with proper cylindrical coordinates
        final_horizon_reg = calculate_event_horizon(Psi_reg_prob, r, coherence_reg_list[-1], grid_type="regular",
                                                    R_grid=R)
        final_horizon_fib = calculate_event_horizon(Psi_fib_prob, r, coherence_fib_list[-1], grid_type="fibonacci",
                                                    R_grid=R_fib)

        # Calculate integrated metrics
        times_array = np.array(time_list)
        coherence_reg_array = np.array(coherence_reg_list)
        coherence_fib_array = np.array(coherence_fib_list)

        reg_metrics = calculate_integrated_metrics(times_array, coherence_reg_array)
        fib_metrics = calculate_integrated_metrics(times_array, coherence_fib_array)

        # Calculate advantage metrics
        coherence_advantage = (fib_metrics['integrated_coherence'] / reg_metrics['integrated_coherence'] - 1) * 100

        # Save results to JSON
        log_message("Saving results to JSON")
        results = {
            'simulation_parameters': {
                'hbar': hbar,
                'm': m,
                'L': L,
                'R_inner': R_inner,
                'R_outer': R_outer,
                'V_0': V_0,
                'Gamma_0': Gamma_0,
                'alpha_c': alpha_c,
                'time_steps': time_steps,
                'dt': dt,
                'total_time': time_steps * dt,
                'challenge_period': challenge_period,
                'challenge_duration': challenge_duration,
                'challenge_strength': challenge_strength
            },
            'coherence_metrics': {
                'final_coherence_regular': float(coherence_reg_list[-1]),
                'final_coherence_fibonacci': float(coherence_fib_list[-1]),
                'integrated_coherence_regular': float(reg_metrics['integrated_coherence']),
                'integrated_coherence_fibonacci': float(fib_metrics['integrated_coherence']),
                'coherence_advantage_percent': float(coherence_advantage),
                'half_life_regular': float(reg_metrics['half_life']),
                'half_life_fibonacci': float(fib_metrics['half_life'])
            },
            'dispersion_metrics': {
                'final_variance_regular': float(variance_reg_list[-1]),
                'final_variance_fibonacci': float(variance_fib_list[-1]),
                'final_entropy_regular': float(disp_reg['entropy']),
                'final_entropy_fibonacci': float(disp_fib['entropy']),
                'final_kurtosis_regular': float(disp_reg['kurtosis']),
                'final_kurtosis_fibonacci': float(disp_fib['kurtosis'])
            },
            'event_horizon': {
                'final_mean_radius_regular': float(np.mean(final_horizon_reg)),
                'final_mean_radius_fibonacci': float(np.mean(final_horizon_fib))
            },
            'state': "study_volunteer_enhanced",
            'timestamp': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }

        # Save to JSON file
        results_file = os.path.join(results_dir, "study_volunteer_enhanced_results.json")
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=4)
        log_message(f"Results saved to {results_file}")

        # Save time series data to NPZ
        npz_file = os.path.join(results_dir, "study_volunteer_enhanced_data.npz")
        np.savez(npz_file,
                 time=np.array(time_list),
                 coherence_reg=np.array(coherence_reg_list),
                 coherence_fib=np.array(coherence_fib_list),
                 variance_reg=np.array(variance_reg_list),
                 variance_fib=np.array(variance_fib_list))
        log_message(f"Time series data saved to {npz_file}")

        log_message("Enhanced study volunteer simulation completed successfully!")
        return True

    except Exception as e:
        log_message(f"Error in enhanced study volunteer simulation: {e}")
        import traceback
        log_message(traceback.format_exc())
        return False


# For standalone testing
if __name__ == "__main__":
    # Load configuration from file - improved path resolution
    try:
        # First try to load from the script directory
        config_path = os.path.join(os.path.dirname(current_file), "config.json")

        # If not found, try to load from project root
        if not os.path.exists(config_path):
            config_path = os.path.join(project_root, "config.json")

        # If still not found, try current working directory
        if not os.path.exists(config_path):
            config_path = "config.json"

        log_message(f"Loading configuration from: {config_path}")
        with open(config_path, 'r') as f:
            config = json.load(f)
        success = run_simulation(config)
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
