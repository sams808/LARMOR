"""make_logo.py — generate LARMOR brand assets into assets/:
  larmor_logo.png    (512×512, window/taskbar icon source)
  larmor_splash.png  (760×440, startup splash)
  larmor.ico         (multi-size Windows icon, for the exe + title bar)

Design: a nuclear spin precessing about B0 (the Larmor precession the app is
named for) traces a helix whose projection radiates an NMR lineshape — the
second-order quadrupolar powder pattern LARMOR was built to fit. Physics in,
resolved spectrum out. Regenerate any time: python make_logo.py
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")

INK = "#0f1a24"          # deep slate-navy ground
TEAL = "#12b0a0"         # LARMOR accent (matches the app's #0e7c86 family)
TEAL_HI = "#4fd1c5"
AMBER = "#f2b134"
CORAL = "#e2705a"
MIST = "#d7e3e6"


def _quad_lineshape(x):
    """A stylized second-order quadrupolar central-transition powder pattern:
    the two horns and the tail LARMOR fits every day."""
    def horn(x0, w):
        return 1.0 / (1.0 + ((x - x0) / w) ** 2)
    y = 0.95 * horn(-0.15, 0.05) + 0.75 * horn(0.28, 0.09)
    y += 0.25 * np.exp(-((x + 0.5) / 0.35) ** 2)     # broad shoulder
    y[x > 0.42] *= np.exp(-(x[x > 0.42] - 0.42) * 6)  # sharp high edge
    return y / y.max()


def _draw_mark(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect("equal")
    ax.axis("off")

    # --- B0 axis (vertical) ---
    ax.add_patch(FancyArrowPatch((2.4, 1.4), (2.4, 8.8), arrowstyle="-|>",
                                 mutation_scale=18, lw=2.2, color=MIST,
                                 zorder=2))
    ax.text(2.4, 9.2, "B₀", color=MIST, fontsize=13, ha="center",
            va="bottom", fontweight="bold")

    # --- Larmor precession helix about B0 ---
    t = np.linspace(0, 4 * np.pi, 400)
    z = np.linspace(2.0, 8.2, t.size)
    rad = 1.15
    xh = 2.4 + rad * np.sin(t)
    # fake perspective: squash x a touch and fade with height
    for i in range(t.size - 1):
        a = 0.35 + 0.6 * (i / t.size)
        ax.plot(xh[i:i + 2], z[i:i + 2], color=TEAL_HI, lw=2.4, alpha=a,
                solid_capstyle="round", zorder=3)
    # the spin vector (magnetic moment) at the top of the helix
    ax.add_patch(FancyArrowPatch((2.4, 5.1), (xh[-1], z[-1]),
                                 arrowstyle="-|>", mutation_scale=16,
                                 lw=2.6, color=AMBER, zorder=4))

    # --- radiated spectrum on the right ---
    xs = np.linspace(-1, 1, 400)
    ys = _quad_lineshape(xs)
    px = 4.8 + (xs + 1) * 2.35        # 4.8 .. 9.5
    py = 2.6 + ys * 4.6               # baseline 2.6
    ax.plot(px, py, color=TEAL, lw=3.0, solid_capstyle="round", zorder=4)
    ax.fill_between(px, 2.6, py, color=TEAL, alpha=0.16, zorder=3)
    # a few colored "fit component" ticks under the pattern
    for xc, col in ((5.4, CORAL), (6.9, AMBER), (8.3, TEAL_HI)):
        ax.plot([xc, xc], [2.6, 2.6 + 0.9], color=col, lw=2.6,
                solid_capstyle="round", alpha=0.9, zorder=4)

    # baseline
    ax.plot([4.7, 9.6], [2.6, 2.6], color=MIST, lw=1.2, alpha=0.5, zorder=2)


def make_logo(path, px=512):
    fig = plt.figure(figsize=(px / 100, px / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    bg = FancyBboxPatch((0.25, 0.25), 9.5, 9.5,
                        boxstyle="round,pad=0.02,rounding_size=1.7",
                        facecolor=INK, edgecolor="none")
    ax.add_patch(bg)
    _draw_mark(ax)
    fig.savefig(path, transparent=True)
    plt.close(fig)


def make_splash(path, w=860, h=440):
    fig = plt.figure(figsize=(w / 100, h / 100), dpi=100)
    fig.patch.set_facecolor(INK)
    ax = fig.add_axes([0.01, 0.10, 0.44, 0.86])
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    _draw_mark(ax)
    fig.text(0.50, 0.63, "LARMOR", color="white", fontsize=52,
             fontweight="bold", family="DejaVu Sans", va="center")
    fig.text(0.505, 0.45, "Lineshape Analysis & Refinement\n"
                          "for Magnetic-resonance Of solids",
             color="#8fb8b4", fontsize=13, va="center")
    fig.text(0.04, 0.045, "an open successor to dmfit — fits with uncertainties",
             color="#5a7370", fontsize=10.5)
    fig.savefig(path, facecolor=INK)
    plt.close(fig)


def make_ico(png_path, ico_path):
    from PIL import Image

    img = Image.open(png_path)
    img.save(ico_path, sizes=[(16, 16), (24, 24), (32, 32), (48, 48),
                              (64, 64), (128, 128), (256, 256)])


if __name__ == "__main__":
    os.makedirs(ASSETS, exist_ok=True)
    logo = os.path.join(ASSETS, "larmor_logo.png")
    make_logo(logo)
    make_splash(os.path.join(ASSETS, "larmor_splash.png"))
    try:
        make_ico(logo, os.path.join(ASSETS, "larmor.ico"))
    except Exception as exc:      # Pillow optional
        print("icon skipped:", exc)
    print("assets written to", ASSETS)
