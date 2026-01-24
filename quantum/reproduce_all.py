#!/usr/bin/env python3
"""
Reproducibility Orchestrator for Noise-Mediated Neuroprotection in HIV.
Runs all core models, aggregates statistical results, and centralizes visualizations.
"""

import os
import subprocess
import pandas as pd
import numpy as np
import arviz as az
from pathlib import Path
from datetime import datetime
import json
import matplotlib.pyplot as plt
import shutil

# --- Configuration ---
MODELS = {
    "Main_v3.6": {
        "module": "quantum.bayesian_v3_6_runner",
        "args": ["--ratio", "3_1_1", "--era", "both", "--draws", "1000", "--tune", "1000"],
        "out_base": "results/bayesian_v3_6/3_1_1/both"
    },
    "Enzyme_v4": {
        "module": "quantum.enzyme_v4_runner",
        "args": ["--ratio", "3_1_1", "--era", "both", "--draws", "1000", "--tune", "1000"],
        "out_base": "results/enzyme_v4/3_1_1/both"
    },
    "Individual_v1": {
        "module": "quantum.hierarchical_individual_v1_runner",
        "args": ["--draws", "1000", "--tune", "1000"],
        "out_base": "results/hierarchical_individual_v1"
    },
    "Regional_v1": {
        "module": "quantum.regional_hierarchical_v1",
        "args": ["--draws", "1000", "--tune", "1000"],
        "out_base": "results/regional_hierarchical_v1"
    }
}

REPRO_DIR = Path("reproducibility_results")

def run_model(name, config):
    print(f"\n>>> Running {name}...")
    cmd = [os.sys.executable, "-m", config["module"]] + config["args"]
    subprocess.run(cmd, check=True)
    
    # Find the latest run directory for this model
    base = Path(config["out_base"])
    if not base.exists():
        return None
        
    runs = sorted([d for d in base.iterdir() if d.is_dir()])
    if not runs:
        # Some runners might output directly to out_base
        return base
    return runs[-1]

def aggregate_stats(run_dirs):
    print("\n>>> Aggregating statistics...")
    summary_rows = []
    
    # 1. Main v3.6
    main_dir = run_dirs.get("Main_v3.6")
    if main_dir:
        trace_path = main_dir / "trace_v3_6.nc"
        if trace_path.exists():
            idata = az.from_netcdf(trace_path)
            post = idata.posterior
            summary_rows.append({
                "Model": "Main_v3.6",
                "beta_xi": float(post["beta_xi"].mean()),
                "xi_acute": float(post["xi_nm_acute"].mean()),
                "xi_chronic": float(post["xi_nm_chronic"].mean()),
                "p_protected": float((post["xi_nm_acute"] < post["xi_nm_chronic"]).mean())
            })

    # 2. Enzyme v4
    enz_dir = run_dirs.get("Enzyme_v4")
    if enz_dir:
        summary_path = enz_dir / "summary_v4.csv"
        if summary_path.exists():
            df = pd.read_csv(summary_path, index_col=0)
            summary_rows.append({
                "Model": "Enzyme_v4",
                "beta_xi": df.loc["beta_xi", "mean"],
                "xi_acute": df.loc["xi_acute", "mean"],
                "xi_chronic": df.loc["xi_chronic", "mean"],
                "p_protected": "N/A" # Mechanistic validation
            })

    # 3. Individual v1
    ind_dir = run_dirs.get("Individual_v1")
    if ind_dir:
        summary_path = ind_dir / "summary.csv"
        if summary_path.exists():
            df = pd.read_csv(summary_path, index_col=0)
            summary_rows.append({
                "Model": "Individual_v1",
                "beta_xi": df.loc["beta_xi", "mean"],
                "xi_acute": df.loc["xi_acute", "mean"],
                "xi_chronic": df.loc["xi_chronic", "mean"],
                "p_protected": "N/A"
            })
            
    # 4. Regional v1
    reg_dir = run_dirs.get("Regional_v1")
    if reg_dir:
        summary_path = reg_dir / "summary.csv"
        if summary_path.exists():
            df = pd.read_csv(summary_path, index_col=0)
            summary_rows.append({
                "Model": "Regional_v1",
                "beta_xi": df.loc["beta_xi_region[0]", "mean"], # BG region
                "xi_acute": df.loc["xi_acute_region[0]", "mean"],
                "xi_chronic": "N/A",
                "p_protected": "N/A"
            })

    master_df = pd.DataFrame(summary_rows)
    master_df.to_csv(REPRO_DIR / "master_summary.csv", index=False)
    print(f"Master summary saved to {REPRO_DIR / 'master_summary.csv'}")
    return master_df

def collect_visuals(run_dirs):
    print("\n>>> Collecting visualizations...")
    viz_dir = REPRO_DIR / "visualizations"
    viz_dir.mkdir(exist_ok=True)
    
    # 1. Main v3.6 (using legacy_figures)
    main_dir = run_dirs.get("Main_v3.6")
    if main_dir:
        subprocess.run([
            os.sys.executable, "-m", "quantum.legacy_figures",
            "--model", "bayesian_v3_6", "--ratio", "3_1_1", "--era", "both",
            "--outdir", str(viz_dir)
        ])

    # Generate Publication Figures (2-5)
    print("\n>>> Generating publication-ready figures (2-5)...")
    subprocess.run([os.sys.executable, "quantum/research/hierarchical_fig_2_5.py"])
    # Copy generated pub figures to visualizations folder
    for f in Path("figures").glob("Figure*.png"):
        shutil.copy(f, viz_dir / f.name)

    # 2. Enzyme v4
    enz_dir = run_dirs.get("Enzyme_v4")
    if enz_dir:
        for f in enz_dir.glob("*.pdf"):
            shutil.copy(f, viz_dir / f.name)
        for f in enz_dir.glob("*.png"):
            shutil.copy(f, viz_dir / f.name)

    # 3. Individual/Regional
    for name in ["Individual_v1", "Regional_v1"]:
        d = run_dirs.get(name)
        if d:
            for f in d.glob("fig*.png"): # Copy our new descriptive filenames
                shutil.copy(f, viz_dir / f"{name}_{f.name}")

def main():
    REPRO_DIR.mkdir(exist_ok=True)
    print("="*60)
    print("NOISE-MEDIATED NEUROPROTECTION REPRODUCIBILITY SUITE")
    print("="*60)
    
    run_dirs = {}
    for name, config in MODELS.items():
        try:
            run_dir = run_model(name, config)
            if run_dir:
                run_dirs[name] = run_dir
                print(f"Completed {name}. Results in {run_dir}")
        except Exception as e:
            print(f"Error running {name}: {e}")
            
    aggregate_stats(run_dirs)
    collect_visuals(run_dirs)
    
    print("\n" + "="*60)
    print("REPRODUCIBILITY RUN COMPLETE")
    print(f"Unified results: {REPRO_DIR}")
    print("="*60)

if __name__ == "__main__":
    main()
