"""Figure for R3: how often each adaptive mask fires with its target absent, against alpha.

Reads results/null_false_positive_rates/null_rates.json (the presets with every target defect
switched off). One panel per sensor; one row per mask, with its Clopper-Pearson 95 % interval and
the trial unit (a stack of images for the stack-calibrated masks, one image for the per-image ones).

Run:  python analysis/figure_null_rates.py
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt

from skmask import plotstyle as ps
from skmask.provenance import write_sidecar

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "null_false_positive_rates" / "null_rates.json"
OUTPUT = ROOT / "results" / "null_false_positive_rates" / "null_rates.png"
SENSORS = {"deep_underground": "deep underground", "shallow_underground": "shallow underground",
           "surface_lab": "surface laboratory"}
LABELS = {"hot_columns": "hot columns & pixels", "cti": "charge-transfer trails", "halo": "halo",
          "serial": "serial-register hits", "low_energy_clusters": "low-energy clusters"}


def main():
    ps.apply()
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    alpha = data["alpha"]
    sensors = [s for s in SENSORS if s in data["per_sensor"]]
    masks = [m for m in LABELS if m in data["per_sensor"][sensors[0]]]

    fig, axes = plt.subplots(1, len(sensors), figsize=(9.6, 2.9), sharey=True, sharex=True,
                             constrained_layout=True)
    x_max = max(row["ci95"][1] for s in sensors for row in
                (data["per_sensor"][s][m] for m in masks) if row["trials"])
    for ax, sensor in zip(axes, sensors):
        rows = data["per_sensor"][sensor]
        for i, mask in enumerate(masks):
            row = rows[mask]
            if not row["trials"]:
                continue
            lo, hi = row["ci95"]
            ax.plot([lo, hi], [i, i], color=ps.SERIES[0], linewidth=2, solid_capstyle="round")
            ax.plot([row["rate"]], [i], **ps.marker(ps.SERIES[0], "o", 6))
            ax.text(hi + 0.01 * x_max, i, f"{row['fired']}/{row['trials']}", va="center",
                    color=ps.INK_SECONDARY, fontsize=7)
        ax.axvline(alpha, color=ps.INK_SECONDARY, linewidth=1)
        ax.set_title(SENSORS[sensor], loc="left", fontsize=10, color=ps.INK)
        ax.grid(axis="y", visible=False)
        ax.set_xlim(-0.03 * x_max * 1.35, x_max * 1.35)
    units = {m: data["per_sensor"][sensors[0]][m]["unit"] for m in masks}
    axes[0].set_yticks(range(len(masks)), [f"{LABELS[m]}\n(per {units[m]})" for m in masks])
    axes[0].set_ylim(len(masks) - 0.4, -0.6)
    fig.supxlabel(f"fraction of trials without the defect in which the mask fired, with 95 % intervals; "
                  f"the vertical line is α = {alpha}", fontsize=9, color=ps.INK_SECONDARY)
    fig.savefig(OUTPUT, dpi=200)
    write_sidecar(OUTPUT, __file__, inputs=[INPUT], notes="figure of R3")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
