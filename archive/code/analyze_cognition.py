
import pandas as pd
import numpy as np
import statsmodels.api as sm
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_cognition():
    base_path = Path(__file__).resolve().parent
    ind_path = base_path / 'data/individual/VALCOUR_2015_INDIVIDUAL_PATIENTS.csv'
    
    if not ind_path.exists():
        print("Data not found.")
        return

    df = pd.read_csv(ind_path)
    # Focus on BGNAA as it's the primary region for the acute surge
    # Clean data: need both BGNAA and NPZ4
    df_clean = df.dropna(subset=['BGNAA', 'NPZ4'])
    
    print(f"Analyzing {len(df_clean)} patients with both BG NAA and NPZ4 scores.")
    
    # Calculate NAA ratio (using healthy baseline 8.76 mM as per guidelines)
    df_clean['NAA_ratio'] = df_clean['BGNAA'] / 8.76
    
    # Correlation
    corr = df_clean['NAA_ratio'].corr(df_clean['NPZ4'])
    print(f"Correlation (NAA Ratio vs NPZ4): {corr:.4f}")
    
    # Regression
    X = sm.add_constant(df_clean['NAA_ratio'])
    model = sm.OLS(df_clean['NPZ4'], X).fit()
    print(model.summary())
    
    # Visualization
    plt.figure(figsize=(10, 6))
    sns.regplot(x='NAA_ratio', y='NPZ4', data=df_clean)
    plt.title('Cognitive Performance (NPZ4) vs Neuronal Protection (NAA Ratio)')
    plt.xlabel('NAA Ratio (Acute Patient / Healthy Control)')
    plt.ylabel('Cognitive Score (NPZ4)')
    plt.axvline(1.0, color='red', linestyle='--', label='Healthy Baseline')
    plt.legend()
    plt.savefig('results/cognition_correlation.png')
    print("Visual saved to results/cognition_correlation.png")

    # Save stats
    with open('results/cognition_stats.txt', 'w') as f:
        f.write(f"Correlation (NAA Ratio vs NPZ4): {corr:.4f}\n")
        f.write(model.summary().as_text())

if __name__ == "__main__":
    analyze_cognition()
