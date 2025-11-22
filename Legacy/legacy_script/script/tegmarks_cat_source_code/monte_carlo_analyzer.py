import os
import subprocess
import time
import json
import pandas as pd

timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
# ---- Settings ----
SUBJECT_GROUPS = ["art_controlled", "chronic_hiv", "acute_hiv", "study_volunteer"]
PHASES = {
    "nucleation": 0,
    "elongation": 1,
    "steady_state": 2,
    "catastrophe": 3
}
GRID_MODE = "both"  # "fibonacci", "uniform", or "both"

# Set your interpreter location (relative to script/)
INTERPRETER_PATH = os.path.join(os.path.dirname(__file__), "monte_carlo_interpreter.py")
PROJECT_ROOT = os.path.expanduser("~/Desktop/Microtubule_Simulation")
TEMP_DIR = os.path.join(PROJECT_ROOT, "script", f"temp_sim_files_{time.strftime('%Y%m%d_%H%M%S')}")

os.makedirs(TEMP_DIR, exist_ok=True)

# ---- Run All Simulations ----
for group in SUBJECT_GROUPS:
    for phase_name, phase_id in PHASES.items():
        print("=" * 80)
        print(f"▶ Subject: {group} | Phase: {phase_name}")

        # Create input file content
        sim_txt = f"""# Auto-generated simulation input
PHASE_ID = {phase_id}
PHASE_NAME = {phase_name}
SUBJECT_GROUP = {group}
TIMESTAMP = {time.strftime("%Y-%m-%d %H:%M:%S")}
"""

        # Write .txt file
        sim_file = os.path.join(TEMP_DIR, f"{group}_{phase_name}_input.txt")
        with open(sim_file, "w") as f:
            f.write(sim_txt)

        # Build command
        command = [
            "python3",
            INTERPRETER_PATH,
            sim_file,
            "--type", group,
            "--grid", GRID_MODE
        ]

        print(f"⏳ Running: {' '.join(command)}")

        # Run the interpreter
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"✅ SUCCESS: {group} | {phase_name}")
        else:
            print(f"❌ ERROR: {group} | {phase_name}")
            print(result.stderr)

# Export metrics
metric_rows = []
grid_types = ["fibonacci", "uniform"]

for group in SUBJECT_GROUPS:
    for grid in grid_types:
        json_file_path = os.path.join("script/Core Scripts/comparison_results", f"{group}_results.json")
        with open(json_file_path, 'r') as json_file:
            results = json.load(json_file)

        key_map = {
            "fibonacci": "fibonacci",
            "uniform": "regular"
        }

        m = results.get("results", results)
        if not m:
            print(f"⚠️ Warning: No metrics found in {json_file_path}")
            continue
            
        try:
            final_coherence = m[f"final_coherence_{key_map[grid]}"]
            half_life = m[f"half_life_{key_map[grid]}"]
            auc = m[f"auc_{key_map[grid]}"]
        except KeyError as e:
            print(f"❌ Missing key {e} in {json_file_path} for grid: {grid}")
            continue

        initial_coherence = 1.0
        coherence_loss = initial_coherence - final_coherence
        coherence_loss_percent = (coherence_loss / initial_coherence) * 100

        row = {
            'timestamp': timestamp,
            'grid': grid,
            'hiv_phase': 'combined',  # or infer from file structure if needed
            'initial_coherence': initial_coherence,
            'final_coherence': final_coherence,
            'coherence_loss': coherence_loss,
            'coherence_loss_percent': coherence_loss_percent,
            'half_life': half_life,
            'auc': auc
        }
        metric_rows.append(row)

df_metrics = pd.DataFrame(metric_rows)
df_metrics.to_csv('il6_metrics.csv', index=False)
print("✅ Saved IL-6 metrics to il6_metrics.csv")

# Export summary-level advantage data
adv_rows = []

for group in SUBJECT_GROUPS:
    json_file_path = os.path.join("script/Core Scripts/comparison_results", f"{group}_results.json")
    with open(json_file_path, 'r') as json_file:
        m = json.load(json_file)

    try:
        coherence_diff = m["final_coherence_fibonacci"] - m["final_coherence_regular"]
        half_life_diff = m["half_life_fibonacci"] - m["half_life_regular"]
        auc_diff = m["auc_fibonacci"] - m["auc_regular"]

        coherence_ratio = m["final_coherence_fibonacci"] / m["final_coherence_regular"]
        half_life_ratio = m["half_life_fibonacci"] / m["half_life_regular"]
        auc_ratio = m["auc_fibonacci"] / m["auc_regular"]

        row = {
            'timestamp': timestamp,
            'hiv_phase': group,
            'coherence_diff': coherence_diff,
            'auc_diff': auc_diff,
            'half_life_diff': half_life_diff,
            'coherence_ratio': coherence_ratio,
            'half_life_ratio': half_life_ratio,
            'auc_ratio': auc_ratio
        }
        adv_rows.append(row)
    except KeyError as e:
        print(f"❌ Missing key {e} in {json_file_path}")
        continue

df_adv = pd.DataFrame(adv_rows)
df_adv.to_csv('il6_advantage.csv', index=False)
print("✅ Saved Fibonacci advantage data to il6_advantage.csv")
