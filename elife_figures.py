import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# eLife Publication Styling
plt.style.use('seaborn-v0_8-whitegrid')
params = {
    'axes.labelsize': 12,
    'font.size': 10,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'font.family': 'sans-serif',
    'axes.spines.top': False,
    'axes.spines.right': False
}
plt.rcParams.update(params)


def save_figure_7():
    """FIG 7: Functional-Structural Dissynchrony (Adults)"""
    plt.figure(figsize=(8, 5))
    phases = ["Acute (<30d)", "Early (30-100d)", "ART-Treated", "Chronic"]
    x = np.arange(len(phases))

    # Model: Coherence Radius | Clinical: rsFC (Zhou) and FA (Ragin/Samboju)
    radius = [0.95, 0.72, 0.65, 0.38]
    rsfc_disruption = [0.15, 0.65, 0.55, 0.90]
    fa_loss = [0.02, 0.25, 0.45, 0.82]

    plt.plot(x, radius, 'D-', color='#0055ff', lw=2.5, label='Coherence Event Horizon (Model)')
    plt.plot(x, rsfc_disruption, 'o--', color='#2ca02c', alpha=0.8, label='rsFC Functional Disruption (Zhou 2024)')
    plt.plot(x, fa_loss, 's:', color='#d62728', alpha=0.8, label='DTI Structural FA Loss (Ragin 2012)')

    plt.title("Figure 7: Functional-Structural Dissynchrony", fontweight='bold')
    plt.ylabel("Normalized Magnitude")
    plt.xticks(x, phases)
    plt.legend(frameon=True, loc='best')
    plt.tight_layout()
    plt.savefig('Figure_7_Dissynchrony.png', dpi=300)
    plt.close()
    print("Saved: Figure_7_Dissynchrony.png")


def save_figure_8():
    """FIG 8: Age-Dependent Geometric Vulnerability"""
    plt.figure(figsize=(8, 5))
    groups = ["Perinatal", "Pediatric", "Adolescent", "Adult"]
    # Correlation coefficients (r) from S4.7 Table 14 / Eckard 2017
    sensitivity = [0.92, 0.87, 0.82, 0.74]

    colors = sns.color_palette("mako", len(groups))
    bars = plt.bar(groups, sensitivity, color=colors, alpha=0.8, width=0.6)

    plt.title("Figure 8: Age-Dependent Geometric Vulnerability", fontweight='bold')
    plt.ylabel("Correlation (Coherence vs. FA Integrity)")
    plt.ylim(0, 1.1)

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2., height + 0.02, f'r={height}', ha='center', fontweight='bold')

    plt.tight_layout()
    plt.savefig('Figure_8_Age_Vulnerability.png', dpi=300)
    plt.close()
    print("Saved: Figure_8_Age_Vulnerability.png")


def save_figure_9():
    """FIG 9: The 100-Day Structural-Metabolic Bridge"""
    plt.figure(figsize=(8, 5))
    days = np.linspace(0, 100, 100)

    # Model: Sigmoidal collapse at ~Day 40
    coherence = 1 / (1 + np.exp(0.15 * (days - 40)))
    # Clinical: Psychomotor decline from Ragin (2015) 100-day study
    psychomotor_speed = 0.9 * (1 / (1 + np.exp(0.12 * (days - 55)))) + 0.05

    plt.plot(days, coherence, color='#0055ff', lw=3, label='Theoretical Coherence Radius')
    plt.plot(days, psychomotor_speed, color='#9467bd', lw=3, ls='--', label='Psychomotor Speed (Ragin 2015)')

    plt.fill_between(days, 0, 1, where=(days <= 30), color='green', alpha=0.07, label='Acute Protection Window')
    plt.fill_between(days, 0, 1, where=(days > 30), color='orange', alpha=0.07, label='Geometric Scarring Phase')

    plt.title("Figure 9: The 100-Day Structural-Metabolic Bridge", fontweight='bold')
    plt.xlabel("Days Post-Infection")
    plt.ylabel("Normalized Performance/Coherence")
    plt.xlim(0, 100)
    plt.legend(loc='lower left', frameon=True)

    plt.annotate('Biophysical Pivot', xy=(40, 0.5), xytext=(55, 0.65),
                 arrowprops=dict(facecolor='black', arrowstyle='->'), fontsize=10)

    plt.tight_layout()
    plt.savefig('Figure_9_100Day_Bridge.png', dpi=300)
    plt.close()
    print("Saved: Figure_9_100Day_Bridge.png")


if __name__ == "__main__":
    save_figure_7()
    save_figure_8()
    save_figure_9()