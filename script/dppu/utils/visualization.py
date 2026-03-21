"""
Visualization Utility
=====================

Standardized plot styles and helper functions for DPPUv2 paper figures.
"""

import matplotlib.pyplot as plt
import seaborn as sns

def set_style():
    """Set global matplotlib style for high-quality figures."""
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    plt.rcParams.update({
        'font.family': 'serif',
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'figure.figsize': (8, 6),
        'lines.linewidth': 2.0,
        'figure.dpi': 300
    })

def save_plot(filename, output_dir="output"):
    """Save plot with standard settings."""
    import os
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, filename)
    plt.tight_layout()
    plt.savefig(path, bbox_inches='tight')
    print(f"Saved figure: {path}")
