"""Figure for R5: which mask fires when a defect that is not its own is the only one present.

Reads results/cross_defect_false_positives/cross_defect.json. One panel per sensor; a row per mask
and a column per defect switched on. The area of each dot is the fraction of trials in which that
mask fired, the number beside it is that fraction, and the colour says whether the fraction is
consistent with alpha or above it. The diagonal, where the mask meets its own defect, is drawn as a
hollow square: firing there is the mask working, not a false positive.

Run:  python analysis/figure_cross_defect.py
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from skmask import plotstyle as ps
from skmask.provenance import write_sidecar

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "cross_defect_false_positives" / "cross_defect.json"
OUTPUT = ROOT / "results" / "cross_defect_false_positives" / "cross_defect.png"
SENSORS = {"deep_underground": "deep underground", "shallow_underground": "shallow underground",
           "surface_lab": "surface laboratory"}
MASKS = {"hot_columns": "hot columns\n& pixels", "cti": "charge-transfer\ntrails", "halo": "halo",
         "serial": "serial-register\nhits", "low_energy_clusters": "low-energy\nclusters"}
GROUPS = {"hot_columns_pixels": "hot columns\n& pixels", "cti": "charge-transfer\ntrails",
          "serial": "serial-register\nhits", "low_energy_clusters": "low-energy\nclusters",
          "halo": "halo"}
MAX_AREA = 210.0          # area of a dot at rate 1


def main():
    ps.apply()
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    alpha = data["alpha"]
    sensors = [s for s in SENSORS if s in data["per_sensor"]]
    masks = list(MASKS)

    fig, axes = plt.subplots(1, len(sensors), figsize=(10.5, 3.6), sharey=True, constrained_layout=True)
    for ax, sensor in zip(axes, sensors):
        groups = [g for g in GROUPS if g in data["per_sensor"][sensor]]
        for x, group in enumerate(groups):
            rows = data["per_sensor"][sensor][group]
            for y, mask in enumerate(masks):
                row = rows[mask]
                if row["own_defect"]:
                    ax.scatter([x], [y], s=70, facecolors="none", edgecolors=ps.AXIS, linewidths=1.2)
                    continue
                if not row["trials"]:
                    continue
                fires = row["verdict"] == "fires on this other defect"
                colour = ps.SERIES[1] if fires else ps.SERIES[0]
                rate = row["rate"]
                ax.scatter([x], [y], s=max(MAX_AREA * rate, 12), color=colour,
                           edgecolors=ps.SURFACE, linewidths=1.2, zorder=3)
                if rate > 0:
                    ax.text(x + 0.22, y, f"{rate:.2f}", va="center", ha="left", fontsize=7,
                            color=ps.INK_SECONDARY)
        ax.set_title(SENSORS[sensor], loc="left", fontsize=10, color=ps.INK)
        ax.set_xticks(range(len(groups)), [GROUPS[g] for g in groups], fontsize=7)
        ax.set_xlim(-0.6, len(groups) - 0.2)
        ax.grid(visible=False)
        ax.tick_params(length=0)

    axes[0].set_yticks(range(len(masks)), [MASKS[m] for m in masks], fontsize=7)
    axes[0].set_ylim(len(masks) - 0.5, -0.5)
    axes[0].set_ylabel("mask", color=ps.INK_SECONDARY)
    handles = [Line2D([], [], **ps.marker(ps.SERIES[0], "o", 7)),
               Line2D([], [], **ps.marker(ps.SERIES[1], "o", 7)),
               Line2D([], [], marker="o", markersize=8, markerfacecolor="none",
                      markeredgecolor=ps.AXIS, linestyle="none")]
    labels = [f"consistent with α = {alpha}", "fires on this other defect", "its own defect"]
    axes[-1].legend(handles, labels, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8,
                    labelcolor=ps.INK_SECONDARY, handletextpad=0.4)
    fig.supxlabel(f"the only defect switched on; dot area and number are the fraction of trials in which "
                  f"the mask fired ({data['n_runs_completed']} runs per configuration)",
                  fontsize=9, color=ps.INK_SECONDARY)
    fig.savefig(OUTPUT, dpi=200, bbox_inches="tight")
    write_sidecar(OUTPUT, __file__, inputs=[INPUT], notes="figure of R5")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
