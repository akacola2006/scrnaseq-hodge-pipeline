#!/usr/bin/env python3
"""
Phase 12: Create NVU Integrated Model Figure

Creates a comprehensive visualization of the 4-node NVU model
showing evidence levels for different relationships.

Author: Claude Code
Date: 2025-11-24
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Output directory
OUTPUT_DIR = 'results/phase12_nvu_integrated'

print("=" * 80)
print("Phase 12: Creating NVU Integrated Model Figure")
print("=" * 80)
print()

# Create figure with two subplots (Control vs ALS)
fig = plt.figure(figsize=(18, 10))

# Define node positions (x, y)
node_positions = {
    'Vascular': (0.15, 0.7),
    'Glia': (0.35, 0.3),
    'Upper': (0.65, 0.7),
    'VAT1L': (0.85, 0.3)
}

# Define colors for each node
node_colors = {
    'Vascular': '#87CEEB',  # Sky blue
    'Glia': '#98D8C8',      # Mint green
    'Upper': '#FFB347',     # Orange
    'VAT1L': '#DDA0DD'      # Plum
}

# Evidence levels for edges
# Format: (source, target, evidence_level, description)
# Evidence levels: 'high', 'moderate', 'exploratory', 'broken_ALS'

edges_control = [
    ('Glia', 'Upper', 'high', 'ρ=-0.78**\n(preserved)'),
    ('Upper', 'VAT1L', 'moderate', 'ρ=+0.32\n(ns)'),
    ('Vascular', 'Glia', 'exploratory', 'PT先行\n(cell-level)'),
    ('Vascular', 'Upper', 'exploratory', 'PT先行\n(cell-level)'),
]

edges_ALS = [
    ('Glia', 'Upper', 'broken_ALS', 'ρ=-0.07\n(BROKEN)'),
    ('Upper', 'VAT1L', 'exploratory', 'ρ=+0.28\n(weak)'),
    ('Vascular', 'Glia', 'exploratory', 'PT先行\n(不明)'),
    ('Vascular', 'Upper', 'exploratory', 'PT先行\n(不明)'),
]

# Edge style mapping
edge_styles = {
    'high': {'color': 'darkgreen', 'width': 4, 'style': '-', 'alpha': 1.0},
    'moderate': {'color': 'orange', 'width': 3, 'style': '-', 'alpha': 0.8},
    'exploratory': {'color': 'gray', 'width': 2, 'style': '--', 'alpha': 0.6},
    'broken_ALS': {'color': 'red', 'width': 2, 'style': ':', 'alpha': 0.5}
}

def draw_nvu_network(ax, title, edges, node_positions, node_colors):
    """Draw NVU network on given axes."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=18, fontweight='bold', pad=20)

    # Draw edges first (so they appear behind nodes)
    for source, target, evidence_level, label in edges:
        src_pos = node_positions[source]
        tgt_pos = node_positions[target]

        style = edge_styles[evidence_level]

        # Create arrow
        arrow = FancyArrowPatch(
            src_pos, tgt_pos,
            arrowstyle='->' if evidence_level != 'broken_ALS' else '-',
            connectionstyle='arc3,rad=0.2',
            color=style['color'],
            linewidth=style['width'],
            linestyle=style['style'],
            alpha=style['alpha'],
            mutation_scale=30,
            zorder=1
        )
        ax.add_patch(arrow)

        # Add edge label
        mid_x = (src_pos[0] + tgt_pos[0]) / 2
        mid_y = (src_pos[1] + tgt_pos[1]) / 2
        # Offset label slightly
        offset_x = 0.02 * (tgt_pos[1] - src_pos[1])  # Perpendicular offset
        offset_y = -0.02 * (tgt_pos[0] - src_pos[0])

        ax.text(mid_x + offset_x, mid_y + offset_y, label,
                fontsize=9, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                         edgecolor=style['color'], alpha=0.8),
                zorder=2)

    # Draw nodes on top
    for node_name, (x, y) in node_positions.items():
        # Node box
        box = FancyBboxPatch(
            (x - 0.08, y - 0.05), 0.16, 0.10,
            boxstyle='round,pad=0.01',
            facecolor=node_colors[node_name],
            edgecolor='black',
            linewidth=3,
            zorder=3
        )
        ax.add_patch(box)

        # Node label
        ax.text(x, y, node_name, fontsize=14, fontweight='bold',
                ha='center', va='center', zorder=4)

# Create Control subplot
ax1 = plt.subplot(1, 2, 1)
draw_nvu_network(ax1, 'Control: Preserved NVU Coupling',
                 edges_control, node_positions, node_colors)

# Create ALS subplot
ax2 = plt.subplot(1, 2, 2)
draw_nvu_network(ax2, 'ALS: NVU Coupling Breakdown',
                 edges_ALS, node_positions, node_colors)

# Add legend
legend_elements = [
    mpatches.Patch(facecolor='darkgreen', edgecolor='black', label='High confidence (strong stat. evidence)'),
    mpatches.Patch(facecolor='orange', edgecolor='black', label='Moderate (suggestive but weak stats)'),
    mpatches.Patch(facecolor='gray', edgecolor='black', label='Exploratory (cell-level or weak patient-level)'),
    mpatches.Patch(facecolor='red', edgecolor='black', label='Broken in ALS (lost coupling)')
]

fig.legend(handles=legend_elements, loc='lower center', ncol=2,
          fontsize=11, frameon=True, fancybox=True, shadow=True,
          bbox_to_anchor=(0.5, -0.02))

# Add overall title
fig.suptitle('Phase 12: Integrated NVU Model - Control vs ALS',
             fontsize=20, fontweight='bold', y=0.98)

# Add caption
caption_text = """
Key Findings:
• Control: Strong Glia ↔ Upper coupling (ρ_partial = -0.78, p = 0.002) suggests preserved homeostatic regulation
• ALS: Complete breakdown of Glia-Upper coupling (ρ_partial = -0.07, ns), indicating loss of NVU coordination
• Vascular involvement: Early PT_dpt at cell level, but patient-level causality unclear (requires validation)
• VAT1L: Fragile component with weak coupling to Upper in both groups

Evidence Levels: High (patient-level, FDR-corrected p < 0.05) | Moderate (p < 0.1 or LiNGAM) | Exploratory (cell-level PT or weak patient stats)
"""

fig.text(0.5, 0.06, caption_text, fontsize=10, ha='center', va='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

plt.tight_layout(rect=[0, 0.15, 1, 0.96])

# Save figure
output_file = f'{OUTPUT_DIR}/NVU_phase12_integrated_model.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"  Saved: {output_file}")

plt.close()

print()
print("=" * 80)
print("Phase 12 NVU Integrated Figure Complete!")
print("=" * 80)
