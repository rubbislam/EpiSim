import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib import font_manager
from pathlib import Path

# Output path: next to this script, in A_pfilter/figures/
out_dir = Path(__file__).resolve().parent
pdf_path = out_dir / "pfilter_graph.pdf"

# Clean sans-serif labels (not hand-drawn). First available wins.
_preferred = ["Helvetica", "Arial", "Helvetica Neue", "DejaVu Sans"]
_available = {f.name for f in font_manager.fontManager.ttflist}
_family = next((n for n in _preferred if n in _available), "DejaVu Sans")
hand = font_manager.FontProperties(family=_family)
print(f"Using label font: {_family}")

fig, ax = plt.subplots(figsize=(15.2, 5.1))
ax.set_xlim(0, 18.2)
ax.set_ylim(0, 6.0)
ax.axis("off")
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

black = "black"
gray = "#cfcfcf"
blue = "#087cf0"
green = "#00845f"
purple = "#7e66ee"

def sketch(artist, scale=1.0, length=80, randomness=1.2):
    return artist   # clean lines: no hand-drawn sketch effect

def arrow(start, end, color=black, lw=3.5, mutation=22, z=5):
    patch = FancyArrowPatch(
        start, end, arrowstyle="->", mutation_scale=mutation,
        linewidth=lw, color=color, shrinkA=0, shrinkB=0,
        capstyle="round", joinstyle="round", zorder=z
    )
    ax.add_patch(patch)
    sketch(patch, 0.9, 90, 1.0)
    return patch

def line(xs, ys, color=black, lw=3.0, alpha=1.0, z=3):
    obj, = ax.plot(
        xs, ys, color=color, linewidth=lw, alpha=alpha,
        solid_capstyle="round", solid_joinstyle="round", zorder=z
    )
    sketch(obj, 0.75, 85, 1.05)
    return obj

def txt(x, y, s, size=22, color=black, ha="center", va="center", rotation=0):
    return ax.text(
        x, y, s, fontproperties=hand, fontsize=size, color=color,
        ha=ha, va=va, rotation=rotation
    )

# Section titles
txt(2.65, 5.70, r"simulate given $\theta$", size=28)
txt(9.05, 5.70, "weight + resample + compute", size=25)
txt(14.85, 5.70, "continue", size=28)

# -----------------------------
# Left panel: simulate particles
# -----------------------------
x0, x1 = 0.65, 5.35
y0, y1 = 1.45, 5.05

arrow((x0, y0), (x1, y0), lw=3.8)
arrow((x0, y0), (x0, y1), lw=3.8)

txt(0.20, 4.85, "X", size=31)
txt(5.25, 1.10, "t", size=27)
txt(1.35, 1.03, r"$t_{n-1}$", size=23)
txt(3.09, 1.03, r"$t_n$", size=23)
txt(1.35, 0.42, "current\nstate", size=21)

# Tick marks. t_n sits at the same RELATIVE x within this panel as t_n does
# in the right ("continue") panel, so the two panels share one time grid.
line([1.35, 1.35], [1.45, 1.72], lw=3.5)
line([3.09, 3.09], [1.45, 1.72], lw=3.5)

# Particle states and simulated trajectories
start_x = 1.35
end_x = 3.09     # trajectories end at t_n (compressed to match the right panel)
start_ys = [4.05, 3.62, 3.28, 2.94, 2.30]

paths_left = [
    ([1.35, 1.60, 1.85, 2.15, 2.48, 2.75, 3.05, 3.30, 3.58, 3.95],
     [4.05, 3.98, 4.13, 4.02, 4.22, 4.05, 4.42, 4.10, 4.72, 4.42]),
    ([1.35, 1.65, 1.95, 2.30, 2.62, 2.88, 3.18, 3.48, 3.72, 3.95],
     [3.62, 3.55, 3.66, 3.52, 3.57, 3.42, 3.88, 3.98, 3.86, 3.78]),
    ([1.35, 1.60, 1.90, 2.20, 2.55, 2.85, 3.12, 3.45, 3.70, 3.95],
     [3.28, 3.12, 3.25, 3.05, 3.15, 2.88, 3.38, 3.50, 3.58, 3.42]),
    ([1.35, 1.62, 1.88, 2.18, 2.46, 2.75, 3.05, 3.33, 3.65, 3.95],
     [2.94, 2.77, 2.40, 2.28, 2.48, 2.30, 2.62, 2.48, 2.74, 2.84]),
    ([1.35, 1.67, 1.96, 2.26, 2.56, 2.90, 3.18, 3.48, 3.72, 3.95],
     [2.30, 2.34, 2.14, 2.04, 2.22, 2.38, 2.55, 2.40, 2.60, 2.48]),
]

# Compress the left-panel trajectories horizontally so they end at t_n=end_x.
L_scale = (end_x - start_x) / (3.95 - 1.35)
for xs, ys in paths_left:
    line([start_x + (x - 1.35) * L_scale for x in xs], ys, color=gray, lw=3.8, z=2)

ax.scatter([start_x] * 5, start_ys, s=145, color=black, zorder=6)
ax.scatter([end_x] * 5, [p[1][-1] for p in paths_left], s=145, color=gray, zorder=5)

txt(2.00, 4.85, r"$f_{X_n\mid X_{n-1}}$", size=25, color="#7a7a7a")

# Purple transition arrow
arrow((5.60, 3.10), (6.25, 3.10), color=purple, lw=4.0, mutation=24)

# -----------------------------------------
# Middle: weight, resample, compute estimate
# -----------------------------------------
row_y = [4.42, 3.78, 3.38, 2.84, 2.18]
labels = [r"$X_{n,1}$", r"$X_{n,2}$", r"$X_{n,3}$", r"$X_{n,4}$", r"$X_{n,5}$"]
weights = [r"$w_{n,1}$", r"$w_{n,2}$", r"$w_{n,3}$", r"$w_{n,4}$", r"$w_{n,5}$"]

for y, lab, wlab in zip(row_y, labels, weights):
    ax.text(7.30, y, lab, fontsize=22, color=gray, ha="right", va="center")
    ax.scatter([7.48], [y], s=120, color=gray, zorder=5)
    line([7.72, 8.98], [y, y], color="#7dc5ff", lw=4.0, alpha=0.75, z=2)
    arrow((8.55, y), (9.12, y), color="#7dc5ff", lw=3.0, mutation=18, z=4)
    ax.text(9.22, y, wlab, fontsize=21, color=blue, ha="left", va="center")

# Blue likelihood curve
curve_y = np.linspace(1.82, 4.78, 160)
curve_x = 8 + 0.55 * np.exp(-((curve_y - 3.40) / 0.85) ** 2)
line(curve_x, curve_y, color=blue, lw=4.0, z=4)
ax.text(8.45, 1.35, r"$f_{Y_n\mid X_n}$", fontsize=27, color=blue,
        ha="center", va="center")

# Green resampling arrow
arrow((10.05, 3.36), (10.72, 3.36), color=green, lw=4.0, mutation=24)

# Resampled particles: duplicated states are represented by repeated dots
resample_x = [10.96, 11.28]
# Resampled swarm: the high-weight middle particles survive and are
# duplicated, drawn at the SAME heights as their sources X_{n,2..4}
# (row_y[1:4] = 3.78, 3.38, 2.84) so the correspondence is legible.
resample_points = [
    (10.96, 3.78), (11.28, 3.78),   # two copies of X_{n,2}
    (10.96, 3.38), (11.28, 3.38),   # two copies of X_{n,3}
    (11.12, 2.84),                  # one copy of X_{n,4}
]
for x, y in resample_points:
    ax.scatter([x], [y], s=120, color=gray, zorder=5)

# Likelihood estimate
arrow((11.12, 2.35), (11.12, 1.65), color=green, lw=4.0, mutation=24)
ax.text(11.12, 1.06, r"$\widehat{\mathcal{L}}_n(\theta)$",
        fontsize=29, color=green, ha="center", va="center")

# Purple arrow to continuation
arrow((12.05, 3.10), (12.70, 3.10), color=purple, lw=4.0, mutation=24)

# ------------------------
# Right panel: continuation
# ------------------------
rx0, rx1 = 13.18, 18.00
ry0, ry1 = 1.45, 5.05
arrow((rx0, ry0), (rx1, ry0), lw=3.8)
arrow((rx0, ry0), (rx0, ry1), lw=3.8)

txt(12.75, 4.85, "X", size=31)
txt(17.90, 1.10, "t", size=27)
txt(13.85, 1.03, r"$t_{n-1}$", size=23)
txt(15.68, 1.03, r"$t_n$", size=23)
txt(17.62, 1.03, r"$t_{n+1}$", size=23)
txt(15.68, 0.42, "current\nstate", size=21)

line([13.85, 13.85], [1.45, 1.72], lw=3.5)
line([15.68, 15.68], [1.45, 1.72], lw=3.5)
line([17.62, 17.62], [1.45, 1.72], lw=3.5)

# Reuse the five initial trajectories, shifted and slightly compressed
shift = 12.50
scale_x = (15.68 - 13.85) / (3.95 - 1.35)
for xs, ys in paths_left:
    xs2 = [13.85 + (x - 1.35) * scale_x for x in xs]
    line(xs2, ys, color=gray, lw=3.8, z=2)

ax.scatter([13.85] * 5, start_ys, s=145, color=black, zorder=6)

# Resampled states at t_n: two upper duplicates, two middle duplicates, one lower
tn_states = [3.82, 3.47, 3.47, 2.76, 3.82]
# Draw slightly offset duplicate dots so multiplicity is visible
tn_xs = [15.64, 15.72, 15.64, 15.68, 15.72]
ax.scatter(tn_xs, tn_states, s=145, color=black, zorder=7)

# Continuation paths from resampled particles
cont_paths = [
    ([15.68, 16.05, 16.42, 16.78, 17.18, 17.62],
     [3.82, 3.92, 4.05, 4.34, 4.58, 4.78]),
    ([15.68, 16.02, 16.35, 16.72, 17.16, 17.62],
     [3.82, 3.66, 3.84, 4.02, 4.22, 4.38]),
    ([15.68, 16.05, 16.40, 16.78, 17.17, 17.62],
     [3.47, 3.56, 3.63, 3.58, 3.64, 3.66]),
    ([15.68, 16.06, 16.38, 16.76, 17.14, 17.62],
     [3.47, 3.18, 3.00, 2.84, 2.72, 2.64]),
    ([15.68, 16.03, 16.42, 16.80, 17.16, 17.62],
     [2.76, 2.62, 2.56, 2.45, 2.31, 2.18]),
]

for xs, ys in cont_paths:
    line(xs, ys, color=gray, lw=3.8, z=2)
    ax.scatter([xs[-1]], [ys[-1]], s=145, color=gray, zorder=5)

fig.subplots_adjust(left=0.01, right=0.995, bottom=0.02, top=0.995)
fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.04)
plt.close(fig)

print(f"Created:\n{pdf_path}")
