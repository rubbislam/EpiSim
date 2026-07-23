import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib import font_manager
from pathlib import Path

# Output path: next to this script, in A_pfilter/figures/
out_dir = Path(__file__).resolve().parent
pdf_path = out_dir / "directsim_graph.pdf"

# Clean sans-serif labels (not hand-drawn). First available wins.
# Arial is preferred over Helvetica because macOS Helvetica is a .ttc with no
# separately addressable bold face, so weight="bold" would render as regular.
_preferred = ["Arial", "Helvetica", "Helvetica Neue", "DejaVu Sans"]
_available = {f.name for f in font_manager.fontManager.ttflist}
_family = next((n for n in _preferred if n in _available), "DejaVu Sans")
hand = font_manager.FontProperties(family=_family)
print(f"Using label font: {_family}")

fig, ax = plt.subplots(figsize=(15.8, 7.2))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")
ax.set_xlim(0, 20)
ax.set_ylim(0, 10)
ax.axis("off")

black = "#000000"
gray = "#d0d0d0"
blue = "#087cf0"
red = "#d00000"
green = "#96d9c0"

def sketched(artist, scale=0.9, length=90, randomness=1.15):
    return artist   # clean lines: no hand-drawn sketch effect

def draw_line(xs, ys, color=black, lw=3.6, alpha=1.0, z=3):
    obj, = ax.plot(
        xs, ys, color=color, linewidth=lw, alpha=alpha,
        solid_capstyle="round", solid_joinstyle="round", zorder=z
    )
    sketched(obj)
    return obj

def draw_arrow(start, end, color=black, lw=3.8, mutation=23, z=5):
    obj = FancyArrowPatch(
        start, end, arrowstyle="->", mutation_scale=mutation,
        linewidth=lw, color=color, shrinkA=0, shrinkB=0,
        capstyle="round", joinstyle="round", zorder=z
    )
    ax.add_patch(obj)
    sketched(obj, scale=0.8, length=100, randomness=1.05)
    return obj

def htext(x, y, s, size=24, color=black, ha="center", va="center",
          rotation=0, z=8, weight="normal"):
    fp = hand.copy()
    fp.set_weight(weight)
    return ax.text(
        x, y, s, fontproperties=fp, fontsize=size, color=color,
        ha=ha, va=va, rotation=rotation, zorder=z
    )

# Titles
htext(2.2, 9.55, "Expected", size=30, weight="bold")
htext(14.15, 9.48, "common cases", size=30, weight="bold")

# Blue log-likelihood expression (simplified: sum over observations j)
ax.text(
    1.5, 8.15,
    r"$\hat{\ell}=\sum_{j}\log f_{Y\mid X}(y_j^{*}\mid x_j^{*})$",
    fontsize=26, color=blue, ha="left", va="center"
)

def draw_axes(x0, x1, y0, y1):
    draw_arrow((x0, y0), (x1, y0), color=black, lw=3.8, mutation=22)
    draw_arrow((x0, y0), (x0, y1), color=black, lw=3.8, mutation=22)
    htext(x0, y1 + 0.28, "case", size=27, ha="center", va="bottom")
    htext(x1 - 0.15, y0 - 0.58, "time", size=26)
    htext(x0 - 0.33, y0 + 0.28, "0", size=22)
    htext(x0 + 0.24, y0 - 0.48, "0", size=22)

# Left panel axes
lx0, lx1 = 1.05, 9.35
ly0, ly1 = 1.15, 8.55
draw_axes(lx0, lx1, ly0, ly1)

# Right panel axes
rx0, rx1 = 11.55, 19.65
ry0, ry1 = 1.15, 8.55
draw_axes(rx0, rx1, ry0, ry1)

# Common x-coordinates
xL = np.array([1.05, 1.85, 2.25, 3.00, 3.72, 4.42, 5.15, 5.55, 6.35, 7.75, 8.38, 8.82, 8.92])
xR = xL - lx0 + rx0

# Black observed data path
data_y = np.array([1.18, 1.70, 1.34, 1.92, 1.58, 4.80, 4.48, 5.15, 6.10, 7.15, 3.55, 3.84, 1.55])

draw_line(xL, data_y, color=black, lw=4.0, z=5)
ax.scatter(xL, data_y, s=135, color=black, zorder=7)

draw_line(xR, data_y, color=black, lw=4.0, z=5)
ax.scatter(xR, data_y, s=135, color=black, zorder=7)

# Gray simulated path in the expected panel
sim_y = np.array([1.22, 1.18, 1.18, 3.72, 3.05, 3.08, 5.70, 7.55, 8.75, 7.40, 5.05, 1.18, 1.65])
draw_line(xL, sim_y, color=green, lw=4.2, z=2)
ax.scatter(xL, sim_y, s=120, color=green, zorder=4)

# Gray baseline simulated observation dots in the right panel
baseline_y = np.full_like(xR, 1.18, dtype=float)
draw_line(xR, baseline_y, color=red, lw=2.5, z=1)
ax.scatter(xR, baseline_y, s=105, color=red, zorder=3)

# Labels in left panel
htext(4.45, 6.55, "sim\n(X)", size=26, color=green, ha="center")
# draw_line([4.98, 5.55], [6.30, 6.15], color=green, lw=3.0, z=4)

htext(6.92, 3.88, "data\n(y*)", size=23, color=black, ha="center", rotation=-4)
draw_line([7.42, 7.76], [4.38, 4.62], color=black, lw=2.8, z=6)

# Labels in right panel
htext(17.10, 3.87, "data\n(y*)", size=23, color=black, ha="center", rotation=-4)
draw_line([17.62, 17.92], [4.40, 4.62], color=black, lw=2.8, z=6)

# Red "common case" trajectory
red_x = np.array([11.58, 12.38, 12.74, 13.22, 13.95, 14.95])
red_y = np.array([1.95, 2.80, 3.55, 4.35, 5.95, 8.65])
draw_line(red_x, red_y, color=red, lw=4.3, z=5)
ax.scatter(red_x, red_y, s=155, color=red, zorder=7)

htext(17.00, 7.90, "inefficient!", size=28, color=red)

fig.subplots_adjust(left=0.015, right=0.995, bottom=0.02, top=0.985)

fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.04)
plt.close(fig)

print(f"Created:\n{pdf_path}")
