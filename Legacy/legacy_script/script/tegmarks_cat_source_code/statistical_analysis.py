# Re-load the data files after execution state reset
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
import ace_tools as tools

# Define file paths
file_paths = {
    "chronic": "/mnt/data/chronic_hiv_data.npz",
    "acute": "/mnt/data/acute_hiv_data.npz",
    "art_controlled": "/mnt/data/art_controlled_data.npz",
    "non_infected": "/mnt/data/non_infected_data.npz",
}

# Load the datasets
datasets = {key: np.load(path) for key, path in file_paths.items()}

# Extract coherence and variance measures
time = datasets["chronic"]["time"]
coherence = {key: datasets[key]["coherence_fib"] for key in datasets}
variance = {key: datasets[key]["variance_fib"] for key in datasets}

# Perform statistical analysis
# Compare coherence between conditions using paired t-tests
stat_results = {}
conditions = list(coherence.keys())

for i in range(len(conditions)):
    for j in range(i + 1, len(conditions)):
        t_stat, p_val = stats.ttest_rel(coherence[conditions[i]], coherence[conditions[j]])
        stat_results[f"{conditions[i]} vs {conditions[j]}"] = (t_stat, p_val)

# Summarize results in a DataFrame
stat_df = pd.DataFrame.from_dict(stat_results, orient="index", columns=["T-Statistic", "P-Value"])
tools.display_dataframe_to_user(name="Statistical Analysis Results", dataframe=stat_df)

# Visualization of Coherence and Variance Trends
plt.figure(figsize=(12, 6))

# Coherence Plot
plt.subplot(1, 2, 1)
for key in coherence:
    plt.plot(time, coherence[key], label=key.capitalize(), linestyle="-", marker="o")
plt.xlabel("Time")
plt.ylabel("Coherence Measure")
plt.title("Coherence Over Time")
plt.legend()
plt.grid()

# Variance Plot
plt.subplot(1, 2, 2)
for key in variance:
    plt.plot(time, variance[key], label=key.capitalize(), linestyle="-", marker="o")
plt.xlabel("Time")
plt.ylabel("Spatial Variance")
plt.title("Variance Over Time")
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()