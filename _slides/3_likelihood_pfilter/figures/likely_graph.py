import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib import font_manager
from scipy.stats import beta
from pathlib import Path

# Output path: next to this script, in A_pfilter/figures/
out_dir = Path(__file__).resolve().parent
pdf_path = out_dir / "likely_graph.pdf"

# Clean sans-serif labels (not hand-drawn). First available wins.
_preferred = ["Helvetica", "Arial", "Helvetica Neue", "DejaVu Sans"]
_available = {f.name for f in font_manager.fontManager.ttflist}
_family = next((n for n in _preferred if n in _available), "DejaVu Sans")
hand_font = font_manager.FontProperties(family=_family)
print(f"Using label font: {_family}")

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "hatch.linewidth": 1.15,
})

fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.25))
fig.subplots_adjust(left=0.035, right=0.985, bottom=0.09, top=0.96, wspace=0.12)

red = "#e01818"
blue = "#1616d8"

def sketched(line, scale=1.2, length=85, randomness=1.6):
    return line   # clean lines: no hand-drawn sketch effect

def draw_axes(ax):
    ax.set_xlim(-0.08, 1.06)
    ax.set_ylim(-0.09, 1.08)
    ax.axis("off")

    x_arrow = FancyArrowPatch(
        (0, 0), (1.03, 0), arrowstyle="->", mutation_scale=24,
        linewidth=4.1, color="black", shrinkA=0, shrinkB=0,
        capstyle="round", joinstyle="round"
    )
    y_arrow = FancyArrowPatch(
        (0, 0), (0, 1.03), arrowstyle="->", mutation_scale=24,
        linewidth=4.1, color="black", shrinkA=0, shrinkB=0,
        capstyle="round", joinstyle="round"
    )
    ax.add_patch(x_arrow)
    ax.add_patch(y_arrow)

    ax.text(-0.055, 1.02, "case", fontproperties=hand_font, fontsize=28,
            ha="center", va="bottom", rotation=8)
    ax.text(0.99, -0.045, "time", fontproperties=hand_font, fontsize=24,
            ha="right", va="top", rotation=-5)
    ax.text(-0.035, 0.015, "0", fontproperties=hand_font, fontsize=20,
            ha="right", va="bottom")
    ax.text(0.012, -0.05, "0", fontproperties=hand_font, fontsize=20,
            ha="left", va="top")

def outer_curve(x, a, b, height):
    y = beta.pdf(x, a, b)
    y = y / np.nanmax(y) * height
    y[(x <= 0) | (x >= 1)] = 0
    return y

def add_panel(ax, color, model_label, outer_ab, outer_height, data_xy,
              inner_xy, label_xy, pointer_xy, data_label_xy, data_pointer_xy):
    draw_axes(ax)

    x = np.linspace(0.001, 0.999, 500)
    y = outer_curve(x, outer_ab[0], outer_ab[1], outer_height)

    # Diagonal shaded model region
    fill = ax.fill_between(
        x, 0, y, facecolor="none", edgecolor=color,
        hatch="/", alpha=0.1, linewidth=0.0, zorder=1
    )
    # Outer model envelope
    outer, = ax.plot(x, y, color=color, linewidth=3.0, alpha=0.3,
                     solid_capstyle="round", zorder=2)
    sketched(outer, scale=1.3, length=75, randomness=1.8)

    # Inner colored model trajectory
    inner_x, inner_y = zip(*inner_xy)
    inner, = ax.plot(inner_x, inner_y, color=color, linewidth=3.0, alpha=0.5,
                     solid_capstyle="round", solid_joinstyle="round", zorder=3)
    sketched(inner, scale=1.1, length=65, randomness=1.5)

    # Observed data trajectory
    data_x, data_y = zip(*data_xy)
    data, = ax.plot(data_x, data_y, color="black", linewidth=5.0,
                    solid_capstyle="round", solid_joinstyle="round", zorder=4)
    sketched(data, scale=0.9, length=90, randomness=1.15)

    # Model label and pointer
    ax.text(label_xy[0], label_xy[1], model_label, color=color,
            fontproperties=hand_font, fontsize=25, ha="center", va="center",
            rotation=7 if color == red else 3)
    pointer, = ax.plot(
        [label_xy[0] - 0.03, pointer_xy[0]],
        [label_xy[1] - 0.08, pointer_xy[1]],
        color=color, linewidth=4.0, solid_capstyle="round", zorder=5
    )
    sketched(pointer, scale=0.8, length=70, randomness=1.1)

    # Data label and leader
    ax.text(data_label_xy[0], data_label_xy[1], "data", color="black",
            fontproperties=hand_font, fontsize=25, ha="left", va="center",
            rotation=-3)
    leader, = ax.plot(
        [data_label_xy[0] - 0.02, data_pointer_xy[0]],
        [data_label_xy[1] - 0.02, data_pointer_xy[1]],
        color="black", linewidth=3.5, solid_capstyle="round", zorder=5
    )
    sketched(leader, scale=0.75, length=75, randomness=1.0)

# Left panel
left_data = [
    (0.00, 0.00), (0.10, 0.075), (0.14, 0.12), (0.18, 0.05),
    (0.27, 0.13), (0.37, 0.09), (0.49, 0.49), (0.54, 0.44),
    (0.68, 0.66), (0.73, 0.37), (0.80, 0.39), (0.82, 0.28),
    (0.83, 0.09), (0.94, 0.00)
]
left_inner = [
    (0.00, 0.00), (0.10, 0.10), (0.20, 0.24), (0.31, 0.42),
    (0.42, 0.55), (0.50, 0.50), (0.58, 0.32), (0.64, 0.20),
    (0.72, 0.16), (0.78, 0.08), (0.90, 0.03), (0.97, 0.00)
]
add_panel(
    axes[0], red, "model 1", (2.35, 3.0), 0.98,
    left_data, left_inner,
    label_xy=(0.71, 0.94), pointer_xy=(0.55, 0.78),
    data_label_xy=(0.82, 0.48), data_pointer_xy=(0.73, 0.45)
)

# Right panel
right_data = [
    (0.00, 0.00), (0.11, 0.11), (0.15, 0.05), (0.24, 0.13),
    (0.34, 0.09), (0.48, 0.55), (0.54, 0.49), (0.64, 0.62),
    (0.79, 0.75), (0.84, 0.37), (0.90, 0.42), (0.91, 0.09),
    (0.99, 0.00)
]
right_inner = [
    (0.00, 0.00), (0.12, 0.10), (0.25, 0.20), (0.36, 0.32),
    (0.49, 0.45), (0.62, 0.58), (0.74, 0.68), (0.82, 0.56),
    (0.87, 0.35), (0.91, 0.18), (0.96, 0.06), (1.00, 0.00)
]
add_panel(
    axes[1], blue, "model 2", (2.0, 1.75), 0.92,
    right_data, right_inner,
    label_xy=(0.31, 0.91), pointer_xy=(0.45, 0.75),
    data_label_xy=(0.88, 0.48), data_pointer_xy=(0.84, 0.44)
)

fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.04)

plt.close(fig)

print(f"Created:\n{pdf_path}")
