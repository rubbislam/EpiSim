import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyArrow, Ellipse, Polygon
from matplotlib import font_manager
from pathlib import Path

# Output path: next to this script, in figures/
out_dir = Path(__file__).resolve().parent
pdf_path = out_dir / "mif_graph.pdf"

# Clean sans-serif labels (not hand-drawn). First available wins.
# Arial is preferred over Helvetica because macOS Helvetica is a .ttc with no
# separately addressable bold face, so weight="bold" would render as regular.
_preferred = ["Arial", "Helvetica", "Helvetica Neue", "DejaVu Sans"]
_available = {f.name for f in font_manager.fontManager.ttflist}
_family = next((n for n in _preferred if n in _available), "DejaVu Sans")
hand = font_manager.FontProperties(family=_family)
print(f"Using label font: {_family}")

fig, ax = plt.subplots(figsize=(16, 9))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")
ax.set_xlim(0, 16)
ax.set_ylim(0, 9)
ax.axis("off")

black = "#000000"
gray = "#cfcfcf"
light_gray = "#e7e7e7"
blue = "#087df0"
red = "#d00000"

def sketch(artist, scale=0.75, length=90, randomness=1.05):
    return artist   # clean lines: no hand-drawn sketch effect

def line(xs, ys, color=black, lw=3.4, alpha=1.0, z=4):
    obj, = ax.plot(xs, ys, color=color, linewidth=lw, alpha=alpha,
                   solid_capstyle="round", solid_joinstyle="round", zorder=z)
    sketch(obj)
    return obj

def arrow(start, end, color=black, lw=3.4, mutation=21, z=5, style="->", alpha=1.0):
    obj = FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=mutation,
                          linewidth=lw, color=color, alpha=alpha,
                          shrinkA=0, shrinkB=0, capstyle="round",
                          joinstyle="round", zorder=z)
    ax.add_patch(obj)
    sketch(obj, 0.7, 95, 1.0)
    return obj

def flow_arrow(start, end, width=0.34, hw=0.72, hl=0.46):
    # A pale-grey BLOCK arrow (shaft + head) connecting the loop panels.
    x0, y0 = start
    x1, y1 = end
    a = FancyArrow(x0, y0, x1 - x0, y1 - y0, width=width,
                   head_width=hw, head_length=hl, length_includes_head=True,
                   facecolor=light_gray, edgecolor="none", zorder=0)
    ax.add_patch(a)
    return a

def htext(x, y, s, size=22, color=black, ha="center", va="center", rotation=0, z=8):
    return ax.text(x, y, s, fontproperties=hand, fontsize=size, color=color,
                   ha=ha, va=va, rotation=rotation, zorder=z)

def mtext(x, y, s, size=22, color=black, ha="center", va="center", z=8):
    return ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, zorder=z)

def axes_panel(x0, y0, width=4.1, height=2.75, ticks=(0.85, 2.65), tick_labels=(r"$t_0$", r"$t_1$")):
    arrow((x0, y0), (x0 + width, y0), lw=3.2, mutation=18)
    arrow((x0, y0), (x0, y0 + height), lw=3.2, mutation=18)
    mtext(x0 - 0.22, y0 + height - 0.02, r"$\theta$", size=25)
    mtext(x0 + width - 0.03, y0 - 0.25, r"$t$", size=23)
    for dx, lab in zip(ticks, tick_labels):
        line([x0 + dx, x0 + dx], [y0, y0 + 0.15], lw=3.0)
        mtext(x0 + dx, y0 - 0.28, lab, size=22)
    return x0, y0, width, height

def colored_tuple(x, y, left_symbol, left_sub, theta_sub, red_left=False, size=25):
    # Places: ( left , theta )
    mtext(x, y, "(", size=size, color=black, ha="left")
    mtext(x+0.18, y, rf"${left_symbol}_{{{left_sub}}}^{{1:5}}$", size=size,
          color=red if red_left else black, ha="left")
    mtext(x+0.88, y, ",", size=size, color=black, ha="left")
    mtext(x+1.06, y, rf"$\theta_{{{theta_sub}}}^{{1:5}}$", size=size,
          color=red if not red_left else black, ha="left")
    mtext(x+1.94, y, ")", size=size, color=black, ha="left")

# Pale grey block arrows connecting the panels around the loop. All the same
# length (0.72), each centred in the gap between the panels it joins.
# (The top-left return into "initialize" is drawn by left_top below.)
flow_arrow((5.89, 7.0), (6.61, 7.0))     # initialize -> perturb  (centred in gap)
flow_arrow((10.89, 7.0), (11.61, 7.0))   # perturb -> simulate    (centred in gap)
flow_arrow((13.75, 5.11), (13.75, 4.39))   # simulate -> filter+resample (down, centred)
flow_arrow((11.76, 2.40), (11.04, 2.40))   # filter+resample -> perturb (left)
flow_arrow((6.78, 2.40), (6.06, 2.40))     # perturb -> simulate (left)

# Left loop return: bent arrows, same arm width (0.34) and head (0.72 x 0.46)
# as the straight flow arrows.
left_top = Polygon([(0.22,5.15),(0.56,5.15),(0.56,6.68),(1.19,6.68),
                    (1.19,6.49),(1.65,6.85),(1.19,7.21),(1.19,7.02),
                    (0.22,7.02)], closed=True, facecolor=light_gray,
                   edgecolor="none", zorder=0)
left_bottom = Polygon([(1.55,1.45),(1.55,1.79),(0.64,1.79),(0.64,2.64),
                       (0.83,2.64),(0.47,3.10),(0.11,2.64),(0.30,2.64),
                       (0.30,1.45)], closed=True, facecolor=light_gray,
                      edgecolor="none", zorder=0)
ax.add_patch(left_top)
ax.add_patch(left_bottom)
for yy in [4.40, 4.10, 3.80]:
    ax.scatter([0.42], [yy], s=85, color=gray, zorder=2)

# Panel positions
p1 = (1.75, 5.75)
p2 = (6.75, 5.75)     # evenly spaced between p1 and p3 -> equal gaps
p3 = (11.75, 5.75)
p5 = (7.05, 1.05)
p6 = (1.80, 1.05)

# Top-left
axes_panel(*p1, width=4.0, height=2.65, ticks=(0.45,), tick_labels=(r"$t_0$",))
htext(2.95, 8.52, "initialize", size=28)
colored_tuple(4.02, 8.50, "x", "0", "0", red_left=True, size=24)

init_x = p1[0] + 0.45
init_ys = np.linspace(p1[1] + 0.65, p1[1] + 2.05, 5)
ax.scatter([init_x]*5, init_ys, s=115, color=black, zorder=7)
ell = Ellipse((init_x, np.mean(init_ys)), width=0.55, height=1.85,
              edgecolor=blue, facecolor="none", linewidth=3.2, zorder=6)
ax.add_patch(ell)
sketch(ell)
htext(3.20, 7.75, "particles", size=25, color=blue)

# Top-middle
axes_panel(*p2, width=4.0, height=2.65, ticks=(0.45, 2.05),
           tick_labels=(r"$t_0$", r"$t_1$"))
htext(7.90, 8.52, "perturb", size=28)
colored_tuple(8.72, 8.50, "x", "0", "1", red_left=False, size=24)

pert0_x = p2[0] + 0.45
pert0_ys = [p2[1]+2.05, p2[1]+1.72, p2[1]+1.45, p2[1]+1.20, p2[1]+0.57]
ax.scatter([pert0_x]*5, pert0_ys, s=115, color=black, zorder=7)

# Top-right
axes_panel(*p3, width=4.0, height=2.65, ticks=(0.45, 2.05),
           tick_labels=(r"$t_0$", r"$t_1$"))
htext(12.82, 8.52, "simulate", size=28)
colored_tuple(13.62, 8.50, "x", "1", "1", red_left=True, size=24)

sim0_x = p3[0] + 0.45
sim1_x = p3[0] + 2.05
sim0_ys = [p3[1]+2.00, p3[1]+1.67, p3[1]+1.41, p3[1]+1.15, p3[1]+0.62]
sim1_ys = [p3[1]+2.00, p3[1]+1.67, p3[1]+1.40, p3[1]+1.15, p3[1]+0.62]
ax.scatter([sim0_x]*5, sim0_ys, s=105, color=black, zorder=7)
for y0, y1 in zip(sim0_ys, sim1_ys):
    arrow((sim0_x+0.12, y0), (sim1_x-0.08, y1),
          color=gray, lw=3.3, mutation=18, z=4)
ax.scatter([sim1_x]*5, sim1_ys, s=105, color=gray, zorder=6)
mtext(12.95, 8, r"$f_{x_1\mid x_0}$", size=25, color=gray)

# Bottom-right: filter + resample
htext(13.85, 4.02, "filter + resample", size=30)
# The whole resampling diagram is shifted DOWN by 0.4 relative to the
# "filter + resample" title (which stays put).
row_y = [2.95, 2.62, 2.25, 1.78, 1.28]
for i, y in enumerate(row_y, start=1):
    mtext(12.80, y, rf"$(x_1^{{{i}}},\theta_1^{{{i}}})$", size=17,
          color=gray, ha="right")
    ax.scatter([13.20], [y], s=105, color=gray, zorder=5)

cy = np.linspace(1.05, 3.25, 160)
# Peak the weight curve at the TOP rows so X_1^1 and X_1^2 are heaviest.
cx = 13.5 + 0.35*np.exp(-((cy-2.78)/0.5)**2)
line(cx, cy, color=blue, lw=3.5, z=6)
mtext(13.7, 3.5, r"$f_{y_1\mid x_1}$", size=24, color=blue)

# Resampled: X_1^1 is duplicated (two copies) and X_1^3 is dropped;
# X_1^2, X_1^4 and X_1^5 each survive once.
out_positions = [
    (14.15, 2.95), (14.55, 2.95),   # X_1^1: two copies
    (14.15, 2.62),                   # X_1^2: one copy
    (14.15, 1.78),                   # X_1^4: one copy
    (14.15, 1.28),                   # X_1^5: one copy
]
for y in row_y:
    arrow((13.28, y), (14.00, y), color="#7ec5ff", lw=2.6, mutation=15, z=4)
for x, y in out_positions:
    ax.scatter([x], [y], s=105, color=blue, zorder=7)
mtext(13.70, 0.66, r"$(x_{1}^{{1:5}^\prime},\,\theta_{1}^{{1:5}^\prime})$",
      size=25, color=blue)

# Bottom-middle
axes_panel(*p5, width=4.0, height=2.65, ticks=(0.45, 2.05, 3.45),
           tick_labels=(r"$t_0$", r"$t_1$", r"$t_2$"))
htext(8.12, 4.02, "perturb", size=28)
colored_tuple(8.82, 4.00, "x", "1", "2", red_left=False, size=24)

x_t0 = p5[0]+0.45
x_t1 = p5[0]+2.05
ys_t0 = [p5[1]+2.02, p5[1]+1.66, p5[1]+1.38, p5[1]+1.12, p5[1]+0.50]
ys_t1 = [p5[1]+2.18, p5[1]+1.92, p5[1]+1.66, p5[1]+0.88, p5[1]+0.47]
ax.scatter([x_t0]*5, ys_t0, s=110, color=black, zorder=7)
ax.scatter([x_t1]*5, ys_t1, s=110, color=black, zorder=7)

# Bottom-left
axes_panel(*p6, width=4.0, height=2.65, ticks=(0.45, 2.05, 3.45),
           tick_labels=(r"$t_0$", r"$t_1$", r"$t_2$"))
htext(2.78, 4.02, "simulate", size=28)
colored_tuple(3.62, 4.00, "x", "2", "2", red_left=True, size=24)

b0 = p6[0] + 0.45
b1 = p6[0] + 2.05
b2 = p6[0] + 3.45
b0_ys = [p6[1]+1.98, p6[1]+1.63, p6[1]+1.37, p6[1]+1.12, p6[1]+0.55]
b1_ys = [p6[1]+2.10, p6[1]+1.82, p6[1]+1.57, p6[1]+0.95, p6[1]+0.53]
# Grey dots at t_2 sit at the SAME heights as their black t_1 dots, so the
# f_{x_2|x_1} arrows are perfectly horizontal.
b2_ys = b1_ys
ax.scatter([b0]*5, b0_ys, s=105, color=black, zorder=7)
ax.scatter([b1]*5, b1_ys, s=105, color=black, zorder=7)
for y0, y1 in zip(b1_ys, b2_ys):
    arrow((b1+0.12, y0), (b2-0.08, y1),
          color=gray, lw=3.3, mutation=18, z=4)
ax.scatter([b2]*5, b2_ys, s=105, color=gray, zorder=6)
mtext(4.5, 3.4, r"$f_{x_2\mid x_1}$", size=24, color=gray)

fig.subplots_adjust(left=0.01, right=0.995, bottom=0.015, top=0.99)
fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.03)
plt.close(fig)

print(f"Created:\n{pdf_path}")
