#!/usr/bin/env python3

import argparse
import random
import statistics
import os
import pandas as pd
from statistics import mean
import json

# -----------------------
# ARGUMENTS
# -----------------------
def load_real_data(subject_group):
    base_path = os.path.expanduser("~/Desktop/Microtubule_Simulation/comparison_results/phase_coherence_data")
    file_path = os.path.join(base_path, subject_group, f"{subject_group}_results.json")

    try:
        with open(file_path, "r") as f:
            data = json.load(f)
            return data
    except FileNotFoundError:
        print(f"[WARN] Real data for '{subject_group}' not found.")
        return None

parser = argparse.ArgumentParser(description="Monte Carlo Microtubule Simulation")
parser.add_argument("input_file", help="Simulation config (template-based .txt)")
parser.add_argument("--type", required=True, help="Subject group (e.g., chronic_hiv)")
parser.add_argument("--grid", choices=["fibonacci", "uniform", "both"], default="both", help="Grid type")
args = parser.parse_args()

# Placeholder phase info (normally replaced by analyzer)
phase_info = {
    "subject_group": args.type,
    "phase_name": "unknown_phase",
}

# Simulated placeholder config — normally parsed from file
# Define all default parameters (baseline values)
parameters = {
    "num_steps": 100,
    "current_phase": 0,
    "nucleation_growth_rate": 0.8,
    "nucleation_catastrophe_prob": 0.05,
    "elongation_growth_rate": 1.2,
    "elongation_catastrophe_prob": 0.03,
    "steady_state_growth_rate": 1.0,
    "steady_state_catastrophe_prob": 0.07,
    "catastrophe_growth_rate": 0.5,
    "catastrophe_catastrophe_prob": 0.15,
    "fib_sequence": [1, 1, 2, 3, 5, 8, 13, 21, 34, 55],
    "growth_adjustment": 1.0,
    "stability_adjustment": 1.0,
    "catastrophe_prob": 0.1,  # this will be overwritten below if real data is found
}

# Inject values from .json if available
real_data = load_real_data(args.type)

if real_data:
    parameters["catastrophe_prob"] = 0.25 if real_data["final_variance_fibonacci"] > 1.2 else 0.12
    parameters["growth_adjustment"] = real_data["final_coherence_fibonacci"] / max(real_data["final_coherence_regular"], 1e-6)
    parameters["stability_adjustment"] = 1.0 / (real_data["final_variance_fibonacci"] + 1e-6)

    print(f"[INFO] Injected real-world data for {args.type}")
    print(f"        catastrophe_prob = {parameters['catastrophe_prob']}")
    print(f"        growth_adjustment = {parameters['growth_adjustment']:.4f}")
    print(f"        stability_adjustment = {parameters['stability_adjustment']:.4f}")
else:
    print(f"[INFO] Using default parameters for {args.type}")

# -----------------------
# FUNCTION DEFINITIONS
# -----------------------
def apply_catastrophe(mt_length):
    return max(0, mt_length - (mt_length * random.uniform(0.3, 0.5)))

def get_fib_step(idx, growth_rate):
    return parameters["fib_sequence"][idx % len(parameters["fib_sequence"])] * growth_rate * 0.1

def get_uniform_step(growth_rate):
    return growth_rate * 0.1

def calculate_stability(history):
    return (1.0 / (statistics.pvariance(history) + 1e-6)) * parameters["stability_adjustment"]

def count_drops(history):
    return sum(1 for i in range(len(history) - 1) if history[i] - history[i + 1] > 0.5)

def find_recovery_end(history, start):
    return next((j for j in range(start + 1, len(history)) if history[j] >= history[start - 1]), len(history))

def find_recovery_periods(history):
    return [
        (i, find_recovery_end(history, i)) for i in range(len(history) - 1)
        if i > 0 and history[i] < history[i - 1] and find_recovery_end(history, i) > i
    ]

def calculate_recovery_time(history):
    periods = find_recovery_periods(history)
    return mean([j - i for i, j in periods]) if periods else 0

# -----------------------
# PHASE CONFIG
# -----------------------
phase = parameters["current_phase"]
growth_rate = (
    parameters["nucleation_growth_rate"] if phase == 0 else
    parameters["elongation_growth_rate"] if phase == 1 else
    parameters["steady_state_growth_rate"] if phase == 2 else
    parameters["catastrophe_growth_rate"]
) * parameters["growth_adjustment"]

catastrophe_prob = (
    parameters["nucleation_catastrophe_prob"] if phase == 0 else
    parameters["elongation_catastrophe_prob"] if phase == 1 else
    parameters["steady_state_catastrophe_prob"] if phase == 2 else
    parameters["catastrophe_catastrophe_prob"]
)

# -----------------------
# SIMULATION EXECUTION
# -----------------------
def run_model(grid_type):
    mt_length = 1
    history = [mt_length]

    for i in range(parameters["num_steps"]):
        if random.random() < catastrophe_prob:
            mt_length = apply_catastrophe(mt_length)
        else:
            step = get_fib_step(i, growth_rate) if grid_type == "fibonacci" else get_uniform_step(growth_rate)
            mt_length += step
        history.append(mt_length)

    return {
        "final_length": history[-1],
        "max_length": max(history),
        "growth_rate": history[-1] / parameters["num_steps"],
        "stability": calculate_stability(history),
        "catastrophe_frequency": count_drops(history) / len(history),
        "recovery_time": calculate_recovery_time(history),
    }

# -----------------------
# RUN & EXPORT
# -----------------------
subject = phase_info["subject_group"]
grid_modes = ["fibonacci", "uniform"] if args.grid == "both" else [args.grid]

for grid_type in grid_modes:
    print(f"[RUN] {subject} | {grid_type}")

    results = run_model(grid_type)

    # Export path
    export_dir = os.path.expanduser(
        f"~/Desktop/Microtubule_Simulation/comparison_results/phase_coherence_data/{subject}/{grid_type}"
    )
    os.makedirs(export_dir, exist_ok=True)

    filename_base = f"{subject}_{grid_type}"
    csv_path = os.path.join(export_dir, f"{filename_base}_mt_results.csv")
    xlsx_path = os.path.join(export_dir, f"{filename_base}_mt_results.xlsx")

    # Add metadata
    results["subject_group"] = subject
    results["grid_type"] = grid_type
    results["phase"] = "unspecified"  # to be updated by analyzer

    # Export
    df = pd.DataFrame([results])
    df.to_csv(csv_path, index=False)
    df.to_excel(xlsx_path, index=False)
    print(f"[EXPORT] Saved to: {csv_path}")