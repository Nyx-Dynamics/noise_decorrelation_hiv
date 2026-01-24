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
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle
from matplotlib.patches import ConnectionPatch
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap

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

    fig = plt.figure(figsize=(12, 9))

    # Create grid layout
    # Panel A: Infection timeline and xi progression
    # Panel B: Microtubule quantum coherence schematic
    # Panel C: NAA protection mechanism
    # Panel D: Key findings summary

    # =========================================================================
    # PANEL A: Infection Timeline (Top Left)
    # =========================================================================
    ax_a = fig.add_axes([0.05, 0.55, 0.42, 0.40])
    ax_a.set_xlim(0, 100)
    ax_a.set_ylim(0, 100)
    ax_a.axis('off')
    ax_a.text(0, 97, 'A', fontsize=14, fontweight='bold', transform=ax_a.transAxes)

    # Timeline arrow
    ax_a.annotate('', xy=(95, 25), xytext=(5, 25),
                  arrowprops=dict(arrowstyle='->', color=COLORS['arrow'], lw=2))
    ax_a.text(50, 18, 'Time since infection', ha='center', fontsize=10)

    # Phase boxes
    phases = [
        ('Acute\n(0-6 mo)', 20, COLORS['acute'], 'High viral load\nInflammation'),
        ('Chronic\n(>6 mo)', 55, COLORS['chronic'], 'Controlled viremia\nART initiated'),
        ('Stable\n(years)', 85, COLORS['healthy'], 'Viral suppression\nImmune reconstitution'),
    ]

    for label, x, color, desc in phases:
        box = FancyBboxPatch((x-12, 35), 24, 25,
                             boxstyle="round,pad=0.03,rounding_size=2",
                             facecolor=color, edgecolor='black', alpha=0.8, lw=1.5)
        ax_a.add_patch(box)
        ax_a.text(x, 47, label, ha='center', va='center', fontsize=10, fontweight='bold', color='white')

        # Description below
        ax_a.text(x, 68, desc, ha='center', va='top', fontsize=8, color='gray')

    # Xi progression curve
    x_curve = np.linspace(5, 95, 100)
    # Acute phase: low xi, then recovery
    xi_curve = 0.4 + 0.4 * (1 - np.exp(-(x_curve - 5) / 30))
    ax_a.plot(x_curve, 85 + 10 * (xi_curve - 0.4) / 0.4, color=COLORS['quantum'],
              lw=3, label=r'$\xi$ (noise correlation)')

    ax_a.text(20, 88, r'$\xi_{acute}$', fontsize=11, color=COLORS['acute'], fontweight='bold')
    ax_a.text(75, 94, r'$\xi_{chronic}$', fontsize=11, color=COLORS['chronic'], fontweight='bold')

    ax_a.set_title('HIV Infection Phases and Noise Correlation Recovery', fontsize=11, fontweight='bold', pad=5)

    # =========================================================================
    # PANEL B: Microtubule Schematic (Top Right)
    # =========================================================================
    ax_b = fig.add_axes([0.53, 0.55, 0.42, 0.40])
    ax_b.set_xlim(0, 100)
    ax_b.set_ylim(0, 100)
    ax_b.axis('off')
    ax_b.text(0, 97, 'B', fontsize=14, fontweight='bold', transform=ax_b.transAxes)

    # Draw microtubule structure
    mt_y = 50
    mt_length = 80
    mt_x_start = 10

    # Tubulin dimers as circles
    n_dimers = 12
    dimer_spacing = mt_length / n_dimers

    for i in range(n_dimers):
        x = mt_x_start + i * dimer_spacing
        # Alpha tubulin
        circle_a = Circle((x, mt_y + 5), 3, facecolor=COLORS['microtubule'],
                          edgecolor='black', lw=0.5)
        ax_b.add_patch(circle_a)
        # Beta tubulin
        circle_b = Circle((x, mt_y - 5), 3, facecolor=COLORS['neuron'],
                          edgecolor='black', lw=0.5)
        ax_b.add_patch(circle_b)

    ax_b.text(50, mt_y + 20, 'Neuronal Microtubule', ha='center', fontsize=10, fontweight='bold')
    ax_b.text(95, mt_y + 5, r'$\alpha$', fontsize=8, color=COLORS['microtubule'])
    ax_b.text(95, mt_y - 5, r'$\beta$', fontsize=8, color=COLORS['neuron'])

    # Quantum coherence illustration
    # Wavy lines showing coherent oscillations
    x_wave = np.linspace(15, 85, 200)

    # Healthy/chronic: long correlation
    y_wave_healthy = mt_y - 25 + 5 * np.sin(x_wave * 0.1) * np.exp(-((x_wave-50)/50)**2)
    ax_b.plot(x_wave, y_wave_healthy, color=COLORS['chronic'], lw=2, alpha=0.8)
    ax_b.text(50, mt_y - 35, r'Chronic: $\xi \approx 0.79$ nm', ha='center',
              fontsize=9, color=COLORS['chronic'])

    # Acute: short correlation (more noise)
    y_wave_acute = mt_y - 50 + 3 * np.sin(x_wave * 0.3) * np.exp(-((x_wave-50)/20)**2)
    for offset in np.linspace(-2, 2, 5):
        ax_b.plot(x_wave, y_wave_acute + offset * 0.5, color=COLORS['acute'],
                  lw=1, alpha=0.4)
    ax_b.text(50, mt_y - 60, r'Acute: $\xi \approx 0.42$ nm', ha='center',
              fontsize=9, color=COLORS['acute'])

    # Coherence length arrows
    ax_b.annotate('', xy=(70, mt_y - 22), xytext=(30, mt_y - 22),
                  arrowprops=dict(arrowstyle='<->', color=COLORS['chronic'], lw=2))
    ax_b.annotate('', xy=(60, mt_y - 48), xytext=(40, mt_y - 48),
                  arrowprops=dict(arrowstyle='<->', color=COLORS['acute'], lw=2))

    ax_b.set_title('Quantum Coherence in Microtubules', fontsize=11, fontweight='bold', pad=5)

    # =========================================================================
    # PANEL C: Protection Mechanism (Bottom Left)
    # =========================================================================
    ax_c = fig.add_axes([0.05, 0.08, 0.42, 0.40])
    ax_c.set_xlim(0, 100)
    ax_c.set_ylim(0, 100)
    ax_c.axis('off')
    ax_c.text(0, 97, 'C', fontsize=14, fontweight='bold', transform=ax_c.transAxes)

    # Mechanism flow diagram
    # HIV -> Inflammation -> Mitochondrial dysfunction -> Reduced ATP -> Lower xi
    # BUT: Lower xi -> Reduced coherence -> Paradoxical NAA preservation

    steps = [
        (15, 80, 'HIV\nInfection', COLORS['hiv']),
        (40, 80, 'Acute\nInflammation', COLORS['acute']),
        (65, 80, 'Mitochondrial\nStress', '#E67E22'),
        (15, 45, 'Reduced\nATP', '#9B59B6'),
        (40, 45, r'Low $\xi$' + '\n(Decorrelated)', COLORS['acute']),
        (65, 45, 'Reduced\nCoherence', COLORS['quantum']),
        (40, 15, 'NAA\nPreserved', COLORS['healthy']),
    ]

    for x, y, label, color in steps:
        box = FancyBboxPatch((x-10, y-8), 20, 16,
                             boxstyle="round,pad=0.02,rounding_size=1",
                             facecolor=color, edgecolor='black', alpha=0.8, lw=1)
        ax_c.add_patch(box)
        ax_c.text(x, y, label, ha='center', va='center', fontsize=8,
                  fontweight='bold', color='white')

    # Arrows
    arrows = [
        ((25, 80), (30, 80)),
        ((50, 80), (55, 80)),
        ((65, 72), (65, 53)),
        ((15, 72), (15, 53)),
        ((25, 45), (30, 45)),
        ((50, 45), (55, 45)),
        ((40, 37), (40, 23)),
    ]

    for start, end in arrows:
        ax_c.annotate('', xy=end, xytext=start,
                      arrowprops=dict(arrowstyle='->', color=COLORS['arrow'], lw=1.5))

    # Protection factor equation
    ax_c.text(85, 45, r'$\Pi_\xi \propto \xi^{\beta_\xi}$' + '\n' + r'$\beta_\xi \approx 2.3$',
              ha='center', va='center', fontsize=11,
              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    ax_c.set_title('Neuroprotection Mechanism', fontsize=11, fontweight='bold', pad=5)

    # =========================================================================
    # PANEL D: Key Findings (Bottom Right)
    # =========================================================================
    ax_d = fig.add_axes([0.53, 0.08, 0.42, 0.40])
    ax_d.set_xlim(0, 100)
    ax_d.set_ylim(0, 100)
    ax_d.axis('off')
    ax_d.text(0, 97, 'D', fontsize=14, fontweight='bold', transform=ax_d.transAxes)

    # Bar chart of xi values
    conditions = ['Acute', 'Chronic', 'Healthy']
    xi_means = [0.425, 0.790, 0.797]
    xi_errors = [0.065, 0.065, 0.048]
    colors = [COLORS['acute'], COLORS['chronic'], COLORS['healthy']]

    bar_width = 15
    bar_positions = [20, 50, 80]

    for i, (pos, mean, err, color) in enumerate(zip(bar_positions, xi_means, xi_errors, colors)):
        # Scale for display
        bar_height = mean * 60
        rect = Rectangle((pos - bar_width/2, 25), bar_width, bar_height,
                         facecolor=color, edgecolor='black', lw=1)
        ax_d.add_patch(rect)

        # Error bar
        ax_d.plot([pos, pos], [25 + bar_height - err*60, 25 + bar_height + err*60],
                  color='black', lw=2)
        ax_d.plot([pos-3, pos+3], [25 + bar_height - err*60, 25 + bar_height - err*60],
                  color='black', lw=2)
        ax_d.plot([pos-3, pos+3], [25 + bar_height + err*60, 25 + bar_height + err*60],
                  color='black', lw=2)

        # Value label
        ax_d.text(pos, 25 + bar_height + 8, f'{mean:.3f}', ha='center', fontsize=9, fontweight='bold')

        # Condition label
        ax_d.text(pos, 18, conditions[i], ha='center', fontsize=10, fontweight='bold')

    # Y-axis label
    ax_d.text(5, 55, r'$\xi$ (nm)', ha='center', va='center', fontsize=10, rotation=90)

    # Significance brackets
    ax_d.plot([20, 20, 50, 50], [78, 82, 82, 78], color='black', lw=1)
    ax_d.text(35, 84, '***', ha='center', fontsize=12)
    ax_d.text(35, 88, 'P > 0.999', ha='center', fontsize=8)

    # Key statistics box
    stats_text = (
        "Key Statistics:\n"
        r"$\bullet$ P($\xi_{acute} < \xi_{chronic}$) > 0.999" + "\n"
        r"$\bullet$ Cohen's d = 5.6 (very large)" + "\n"
        r"$\bullet$ 86% recovery in chronic phase" + "\n"
        r"$\bullet$ $\beta_\xi$ = 2.33 $\pm$ 0.51"
    )
    ax_d.text(50, 5, stats_text, ha='center', va='bottom', fontsize=8,
              bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3))

    ax_d.set_title('Noise Correlation by HIV Phase', fontsize=11, fontweight='bold', pad=5)

    # =========================================================================
    # Save figure
    # =========================================================================
    output_path = '/Users/acdmbpmax/Desktop/noise canonical/figures/Figure1_conceptual_overview.png'
    import os
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Figure saved to: {output_path}")

    # Also save as PDF for publication
    pdf_path = output_path.replace('.png', '.pdf')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"PDF saved to: {pdf_path}")

    plt.close()


if __name__ == '__main__':
    create_figure1()
    print("\nFigure 1 generation complete!")
