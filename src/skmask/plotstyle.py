"""One visual style for every figure: recessive chrome, thin marks, validated colours.

Colours are the first two slots of a categorical palette checked for colour-vision deficiency
(adjacent CVD Delta E >= 8, OKLab x100) with the dataviz validator; text never takes a series colour.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834"]   # slot 1 blue, slot 2 orange


def apply():
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": ["Segoe UI", "system-ui", "DejaVu Sans", "sans-serif"],
        "font.size": 9,
        "text.color": INK,
        "axes.labelcolor": INK_SECONDARY,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "grid.linestyle": "-",
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK_SECONDARY,
        "ytick.labelcolor": INK_SECONDARY,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "lines.linewidth": 2.0,
        "lines.solid_capstyle": "round",
    })


def marker(color, shape="o", size=7):
    """Keyword arguments for a data marker with a surface-coloured ring."""
    return dict(marker=shape, markersize=size, markerfacecolor=color, markeredgecolor=SURFACE,
                markeredgewidth=1.5, linestyle="none")
