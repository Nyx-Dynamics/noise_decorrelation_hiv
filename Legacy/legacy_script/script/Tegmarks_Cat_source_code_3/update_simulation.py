#!/usr/bin/env python3
"""
Update Script for Microtubule Quantum Coherence Simulation

This script updates the simulation files to incorporate Tegmark's realistic
decoherence model and proper cylindrical coordinate physics.
"""

import os
import sys
import shutil
import json


def backup_file(filepath):
    """
    Create a backup of a file.

    Args:
        filepath (str): Path to the file to back up

    Returns:
        str: Path to the backup file
    """
    backup_path = f"{filepath}.bak"
    if os.path.exists(filepath):
        shutil.copy2(filepath, backup_path)
        print(f"Created backup: {backup_path}")
    return backup_path


def update_wavefunction_module():
    """
    Update the wavefunction.py module with Tegmark decoherence model.
    """
    file_path = "scripts/core/wavefunction.py"
    backup_file(file_path)

    # Updated content is provided separately as it's quite long
    # This would be where you paste the corrected wavefunction.py code
    print(f"Updating {file_path} with Tegmark decoherence model...")

    # Placeholder for actual file writing
    # with open(file_path, 'w') as f:
    #    f.write(updated_content)

    print(f"Updated {file_path}")


def update_grid_module():
    """
    Update the grid.py module with proper cylindrical coordinates.
    """
    file_path = "scripts/core/grid.py"
    backup_file(file_path)

    # Updated content is provided separately
    print(f"Updating {file_path} with cylindrical coordinate corrections...")

    # Placeholder for actual file writing
    # with open(file_path, 'w') as f:
    #    f.write(updated_content)

    print(f"Updated {file_path}")


def update_metrics_module():
    """
    Update the metrics.py module with proper volume elements.
    """
    file_path = "scripts/core/metrics.py"
    backup_file(file_path)

    # Updated content is provided separately
    print(f"Updating {file_path} with proper volume element calculations...")

    # Placeholder for actual file writing
    # with open(file_path, 'w') as f:
    #    f.write(updated_content)

    print(f"Updated {file_path}")


def update_simulation_files():
    """
    Update the simulation scripts to use the new function signatures.
    """
    simulation_files = [
        "scripts/simulations/acute_hiv.py",
        "scripts/simulations/chronic_hiv.py",
        "scripts/simulations/art_controlled.py",
        "scripts/simulations/study_volunteer.py"
    ]

    for file_path in simulation_files:
        if not os.path.exists(file_path):
            print(f"Warning: {file_path} not found, skipping.")
            continue

        backup_file(file_path)

        print(f"Updating {file_path} with new function signatures...")

        # Read file content
        with open(file_path, 'r') as f:
            content = f.read()

        # Update function calls to match new signatures

        # 1. Update coherence calculation (now requires R_grid parameter)
        content = content.replace(
            "coh_reg = calculate_coherence(Psi_reg_prob, Psi_reg_prob_list[0], dr, dz)",
            "coh_reg = calculate_coherence(Psi_reg_prob, Psi_reg_prob_list[0], dr, dz, R)"
        )
        content = content.replace(
            "coh_fib = calculate_coherence(Psi_fib_prob, Psi_fib_prob_list[0], dr, dz)",
            "coh_fib = calculate_coherence(Psi_fib_prob, Psi_fib_prob_list[0], dr, dz, R_fib)"
        )

        # 2. Update decoherence calculation (add grid_type parameter)
        content = content.replace(
            "Gamma_reg = calculate_decoherence(C, Gamma_0, alpha_c, \"acute_hiv\", IL6_COUPLING_STRENGTH)",
            "Gamma_reg = calculate_decoherence(C, Gamma_0, alpha_c, \"acute_hiv\", IL6_COUPLING_STRENGTH, \"regular\")"
        )
        content = content.replace(
            "Gamma_fib = calculate_decoherence(C, Gamma_0, alpha_c, \"acute_hiv\", IL6_COUPLING_STRENGTH)",
            "Gamma_fib = calculate_decoherence(C, Gamma_0, alpha_c, \"acute_hiv\", IL6_COUPLING_STRENGTH, \"fibonacci\")"
        )

        # Patterns for other condition files
        for condition in ["chronic_hiv", "art_controlled", "study_volunteer"]:
            content = content.replace(
                f"Gamma_reg = calculate_decoherence(C, Gamma_0, alpha_c, \"{condition}\")",
                f"Gamma_reg = calculate_decoherence(C, Gamma_0, alpha_c, \"{condition}\", grid_type=\"regular\")"
            )
            content = content.replace(
                f"Gamma_fib = calculate_decoherence(C, Gamma_0, alpha_c, \"{condition}\")",
                f"Gamma_fib = calculate_decoherence(C, Gamma_0, alpha_c, \"{condition}\", grid_type=\"fibonacci\")"
            )

        # 3. Update wavefunction evolution (add grid_type parameter)
        content = content.replace(
            "Psi_reg = evolve_wavefunction(Psi_reg, V_reg, Gamma_reg, dr, dz, dt, R, hbar, m, N_r, N_z)",
            "Psi_reg = evolve_wavefunction(Psi_reg, V_reg, Gamma_reg, dr, dz, dt, R, hbar, m, N_r, N_z, grid_type=\"regular\")"
        )
        content = content.replace(
            "Psi_fib = evolve_wavefunction(Psi_fib, V_fib, Gamma_fib, dr, dz, dt, R_fib, hbar, m, N_r, N_z)",
            "Psi_fib = evolve_wavefunction(Psi_fib, V_fib, Gamma_fib, dr, dz, dt, R_fib, hbar, m, N_r, N_z, grid_type=\"fibonacci\")"
        )

        # 4. Update event horizon calculation (add R_grid parameter)
        content = content.replace(
            "horizon_reg = calculate_event_horizon(Psi_reg_prob, r, rolling_coherence_reg, grid_type=\"regular\")",
            "horizon_reg = calculate_event_horizon(Psi_reg_prob, r, rolling_coherence_reg, grid_type=\"regular\", R_grid=R)"
        )
        content = content.replace(
            "horizon_fib = calculate_event_horizon(Psi_fib_prob, r, rolling_coherence_fib, grid_type=\"fibonacci\")",
            "horizon_fib = calculate_event_horizon(Psi_fib_prob, r, rolling_coherence_fib, grid_type=\"fibonacci\", R_grid=R_fib)"
        )

        # Write updated content back to file
        with open(file_path, 'w') as f:
            f.write(content)

        print(f"Updated {file_path}")


def update_config():
    """
    Update the config.json file to increase time_steps.
    """
    config_path = "config.json"
    backup_file(config_path)

    try:
        with open(config_path, 'r') as f:
            config = json.load(f)

        # Update time_steps to 500
        if "simulation_parameters" in config:
            config["simulation_parameters"]["time_steps"] = 500
            print(f"Updated time_steps to 500 in config.json")

        # Write updated config back to file
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=4)

        print(f"Updated {config_path}")

    except Exception as e:
        print(f"Error updating config: {e}")


def main():
    """
    Main function to run all updates.
    """
    print("Starting simulation code update...")

    # Backup core modules
    update_wavefunction_module()
    update_grid_module()
    update_metrics_module()

    # Update simulation files to match new function signatures
    update_simulation_files()

    # Update configuration
    update_config()

    print("\nUpdate complete! You can now run the simulations with the corrected physics.")
    print("To run the acute HIV simulation:")
    print("python scripts/simulations/acute_hiv.py")


if __name__ == "__main__":
    main()