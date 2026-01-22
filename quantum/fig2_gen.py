#!/usr/bin/env python3
"""
Figure 1: Conceptual Overview
==============================
Graphical abstract showing the relationship between:
- HIV infection phases (acute -> chronic)
- Noise correlation length (xi)
- Quantum coherence protection mechanism
- NAA metabolite preservation

For Nature Communications submission.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, FancyArrowPatch
from matplotlib.patches import ArrowStyle
import os

# Set up publication quality
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

# Color palette (colorblind-friendly)
COLORS = {
    'acute': '#E64B35',      # Red
    'chronic': '#4DBBD5',    # Teal
    'healthy': '#00A087',    # Green
    'hiv': '#7E6148',        # Brown
    'neuron': '#3C5488',     # Navy
    'microtubule': '#8491B4',# Slate
    'quantum': '#F39B7F',    # Salmon
    'arrow': '#333333',      # Dark gray
}


def create_figure1():
    """Create the main conceptual overview figure."""

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # =========================================================================
    # PANEL A: Infection Timeline (Top Left)
    # =========================================================================
    ax = axes[0, 0]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')
    ax.text(-0.05, 1.05, 'A', fontsize=14, fontweight='bold', transform=ax.transAxes)

    # Timeline arrow
    ax.annotate('', xy=(95, 20), xytext=(5, 20),
                arrowprops=dict(arrowstyle='->', color=COLORS['arrow'], lw=2))
    ax.text(50, 12, 'Time since infection', ha='center', fontsize=10)

    # Phase boxes
    phases = [
        ('Acute\n(0-6 mo)', 20, COLORS['acute']),
        ('Chronic\n(>6 mo)', 55, COLORS['chronic']),
        ('Stable\n(years)', 85, COLORS['healthy']),
    ]

    for label, x, color in phases:
        box = FancyBboxPatch((x-12, 30), 24, 20,
                             boxstyle="round,pad=0.02,rounding_size=2",
                             facecolor=color, edgecolor='black', alpha=0.8, lw=1.5)
        ax.add_patch(box)
        ax.text(x, 40, label, ha='center', va='center', fontsize=9, fontweight='bold', color='white')

    # Xi progression curve
    x_curve = np.linspace(8, 92, 100)
    xi_curve = 0.4 + 0.4 * (1 - np.exp(-(x_curve - 8) / 35))
    y_curve = 60 + 30 * (xi_curve - 0.4) / 0.4
    ax.plot(x_curve, y_curve, color=COLORS['quantum'], lw=3)
    ax.fill_between(x_curve, 60, y_curve, alpha=0.2, color=COLORS['quantum'])

    # Labels on curve
    ax.text(18, 65, r'$\xi_{acute}$' + '\n0.42 nm', fontsize=9, color=COLORS['acute'],
            fontweight='bold', ha='center')
    ax.text(75, 88, r'$\xi_{chronic}$' + '\n0.79 nm', fontsize=9, color=COLORS['chronic'],
            fontweight='bold', ha='center')

    ax.set_title('HIV Infection Phases and Noise Correlation Recovery', fontsize=11, fontweight='bold', pad=10)

    # =========================================================================
    # PANEL B: Microtubule Schematic (Top Right)
    # =========================================================================
    ax = axes[0, 1]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')
    ax.text(-0.05, 1.05, 'B', fontsize=14, fontweight='bold', transform=ax.transAxes)

    # Draw microtubule structure
    mt_y = 70
    n_dimers = 10

    for i in range(n_dimers):
        x = 10 + i * 8
        # Alpha tubulin (top)
        circle_a = Circle((x, mt_y + 4), 3, facecolor=COLORS['microtubule'],
                          edgecolor='black', lw=0.5)
        ax.add_patch(circle_a)
        # Beta tubulin (bottom)
        circle_b = Circle((x, mt_y - 4), 3, facecolor=COLORS['neuron'],
                          edgecolor='black', lw=0.5)
        ax.add_patch(circle_b)

    ax.text(50, mt_y + 15, 'Neuronal Microtubule', ha='center', fontsize=10, fontweight='bold')
    ax.text(92, mt_y + 4, r'$\alpha$', fontsize=8, color=COLORS['microtubule'])
    ax.text(92, mt_y - 4, r'$\beta$', fontsize=8, color=COLORS['neuron'])

    # Chronic: long correlation wave
    x_wave = np.linspace(10, 90, 150)
    y_chronic = 45 + 8 * np.sin(x_wave * 0.08) * np.exp(-((x_wave-50)/60)**2)
    ax.plot(x_wave, y_chronic, color=COLORS['chronic'], lw=2.5, label='Chronic')
    ax.fill_between(x_wave, 45, y_chronic, alpha=0.2, color=COLORS['chronic'])

    # Correlation length arrow for chronic
    ax.annotate('', xy=(75, 46), xytext=(25, 46),
                arrowprops=dict(arrowstyle='<->', color=COLORS['chronic'], lw=2))
    ax.text(50, 38, r'$\xi_{chronic} \approx 0.79$ nm', ha='center', fontsize=9,
            color=COLORS['chronic'], fontweight='bold')

    # Acute: short correlation (noisy)
    y_acute = 20 + 4 * np.sin(x_wave * 0.25) * np.exp(-((x_wave-50)/25)**2)
    ax.plot(x_wave, y_acute, color=COLORS['acute'], lw=2.5, label='Acute')

    # Add noise to acute
    for offset in [-1.5, 1.5]:
        y_noise = y_acute + offset + np.random.normal(0, 0.5, len(x_wave))
        ax.plot(x_wave, y_noise, color=COLORS['acute'], lw=0.5, alpha=0.4)

    # Correlation length arrow for acute
    ax.annotate('', xy=(60, 21), xytext=(40, 21),
                arrowprops=dict(arrowstyle='<->', color=COLORS['acute'], lw=2))
    ax.text(50, 10, r'$\xi_{acute} \approx 0.42$ nm', ha='center', fontsize=9,
            color=COLORS['acute'], fontweight='bold')

    ax.set_title('Quantum Coherence in Microtubules', fontsize=11, fontweight='bold', pad=10)

    # =========================================================================
    # PANEL C: Protection Mechanism (Bottom Left)
    # =========================================================================
    ax = axes[1, 0]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')
    ax.text(-0.05, 1.05, 'C', fontsize=14, fontweight='bold', transform=ax.transAxes)

    # Flow diagram boxes
    boxes = [
        (20, 85, 'HIV\nInfection', COLORS['hiv']),
        (50, 85, 'Acute\nInflammation', COLORS['acute']),
        (80, 85, 'Mitochondrial\nStress', '#E67E22'),
        (20, 55, 'Reduced\nATP', '#9B59B6'),
        (50, 55, r'Low $\xi$', COLORS['acute']),
        (80, 55, 'Reduced\nCoherence', COLORS['quantum']),
        (50, 25, 'NAA\nPreserved', COLORS['healthy']),
    ]

    for x, y, label, color in boxes:
        box = FancyBboxPatch((x-12, y-8), 24, 16,
                             boxstyle="round,pad=0.02,rounding_size=1",
                             facecolor=color, edgecolor='black', alpha=0.85, lw=1)
        ax.add_patch(box)
        ax.text(x, y, label, ha='center', va='center', fontsize=8,
                fontweight='bold', color='white')

    # Arrows connecting boxes
    arrow_style = dict(arrowstyle='->', color=COLORS['arrow'], lw=1.5)
    ax.annotate('', xy=(38, 85), xytext=(32, 85), arrowprops=arrow_style)
    ax.annotate('', xy=(68, 85), xytext=(62, 85), arrowprops=arrow_style)
    ax.annotate('', xy=(80, 77), xytext=(80, 63), arrowprops=arrow_style)
    ax.annotate('', xy=(20, 77), xytext=(20, 63), arrowprops=arrow_style)
    ax.annotate('', xy=(38, 55), xytext=(32, 55), arrowprops=arrow_style)
    ax.annotate('', xy=(68, 55), xytext=(62, 55), arrowprops=arrow_style)
    ax.annotate('', xy=(50, 47), xytext=(50, 33), arrowprops=arrow_style)

    # Protection equation box
    eq_box = FancyBboxPatch((60, 15), 35, 20,
                            boxstyle="round,pad=0.02",
                            facecolor='wheat', edgecolor='black', alpha=0.9, lw=1)
    ax.add_patch(eq_box)
    ax.text(77.5, 25, r'$\Pi_\xi \propto \xi^{\beta_\xi}$' + '\n' + r'$\beta_\xi = 2.33 \pm 0.51$',
            ha='center', va='center', fontsize=10)

    ax.set_title('Neuroprotection Mechanism', fontsize=11, fontweight='bold', pad=10)

    # =========================================================================
    # PANEL D: Key Findings Bar Chart (Bottom Right)
    # =========================================================================
    ax = axes[1, 1]
    ax.text(-0.05, 1.05, 'D', fontsize=14, fontweight='bold', transform=ax.transAxes)

    # Bar chart data
    conditions = ['Acute', 'Chronic', 'Healthy']
    xi_means = [0.425, 0.790, 0.797]
    xi_errors = [0.065, 0.065, 0.048]
    colors = [COLORS['acute'], COLORS['chronic'], COLORS['healthy']]

    x_pos = np.arange(len(conditions))
    bars = ax.bar(x_pos, xi_means, yerr=xi_errors, capsize=5, color=colors,
                  edgecolor='black', lw=1.5, alpha=0.85)

    # Value labels on bars
    for i, (bar, mean) in enumerate(zip(bars, xi_means)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + xi_errors[i] + 0.02,
                f'{mean:.3f}', ha='center', fontsize=10, fontweight='bold')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(conditions, fontsize=11)
    ax.set_ylabel(r'Noise Correlation Length $\xi$ (nm)', fontsize=11)
    ax.set_ylim(0, 1.0)

    # Significance bracket
    ax.plot([0, 0, 1, 1], [0.88, 0.90, 0.90, 0.88], color='black', lw=1.5)
    ax.text(0.5, 0.91, '***', ha='center', fontsize=14)
    ax.text(0.5, 0.95, 'P > 0.999', ha='center', fontsize=9)

    # Statistics text box
    stats_text = ("Key Statistics:\n"
                  r"$\bullet$ P($\xi_{acute} < \xi_{chronic}$) > 0.999" + "\n"
                  r"$\bullet$ Cohen's d = 5.63" + "\n"
                  r"$\bullet$ 86% increase acute$\rightarrow$chronic")

    ax.text(0.98, 0.02, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3))

    ax.set_title('Noise Correlation by HIV Phase', fontsize=11, fontweight='bold', pad=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # =========================================================================
    # Save figure
    # =========================================================================
    plt.tight_layout()

    output_dir = '/Users/acdmbpmax/Desktop/noise canonical/figures'
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, 'Figure1_conceptual_overview.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Figure saved to: {output_path}")

    pdf_path = output_path.replace('.png', '.pdf')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"PDF saved to: {pdf_path}")

    plt.close()


if __name__ == '__main__':
    create_figure1()
    print("\nFigure 1 generation complete!")
