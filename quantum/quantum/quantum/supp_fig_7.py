import matplotlib.pyplot as plt
import numpy as np

# Data for Cross-System Validation
systems = ['Tegmark Limit', 'FMO Complex', 'Cryptochrome', 'HIV (Chronic)', 'HIV (Acute/BG)']
xi_values = [5.0, 2.0, 0.7, 0.79, 0.425]  # nm
xi_errors = [0.5, 0.3, 0.2, 0.065, 0.065]

# Data for Evolutionary Signal (from regional_summary.csv)
regions = ['Basal Ganglia', 'Parietal GM', 'Frontal WM', 'Frontal GM']
age_my = [500, 300, 200, 50]  # Million Years
xi_acute = [0.50, 0.59, 0.55, 0.63]

plt.figure(figsize=(12, 6))

# Panel A: Universal Correlation Scales
plt.subplot(1, 2, 1)
colors = ['gray', 'forestgreen', 'royalblue', 'darkorange', 'crimson']
plt.barh(systems, xi_values, xerr=xi_errors, color=colors, alpha=0.7)
plt.axvline(x=1.0, color='black', linestyle='--', label='Sub-nanometer Threshold')
plt.xlabel('Correlation Length ξ (nm)')
plt.title('A. Universal Correlation Scales')
plt.grid(axis='x', linestyle=':', alpha=0.6)

# Panel B: Evolutionary Optimization Signal
plt.subplot(1, 2, 2)
plt.scatter(age_my, xi_acute, s=150, c='crimson', edgecolors='black', label='Regional Data')

# Exponential Fit from evolutionary_stats.txt: xi = 0.53 + 0.13 * exp(-0.0057 * age)
x_fit = np.linspace(0, 550, 100)
y_fit = 0.53 + 0.13 * np.exp(-0.0057 * x_fit)
plt.plot(x_fit, y_fit, 'k--', label='Evolutionary Asymptote')

for i, txt in enumerate(regions):
    plt.annotate(txt, (age_my[i], xi_acute[i]), xytext=(5, 5), textcoords='offset points')

plt.xlabel('Phylogenetic Age (Million Years)')
plt.ylabel('Inferred ξ_acute (nm)')
plt.title('B. Evolutionary Optimization Signal')
plt.legend()

plt.tight_layout()
plt.savefig('SuppFig7.png', dpi=300)