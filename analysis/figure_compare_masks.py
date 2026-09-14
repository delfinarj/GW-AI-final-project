"""Figure for R4: each mask's figure of merit relative to the oracle tuned on the same sensor.

Reads results/compare_masks/compare_masks.json. Small multiples, one panel per sensor, one row per
mask (and the six combined). The oracle is the reference line at 1; blue circles are the adaptive
masks, orange markers the fixed masks transplanted from each of the other two sensors (marker shape
names the source, in the legend). One x-axis for all panels.

Run:  python analysis/figure_compare_masks.py
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from skmask import plotstyle as ps
from skmask.provenance import write_sidecar

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "compare_masks" / "compare_masks.json"
OUTPUT = ROOT / "results" / "compare_masks" / "compare_masks.png"
SENSORS = {"deep_underground": "deep underground", "shallow_underground": "shallow underground",
           "surface_lab": "surface laboratory"}
MASKS = {"hot": "hot columns & pixels", "cti": "charge-transfer trails", "halo": "halo",
         "serial": "serial-register hits", "lec": "low-energy clusters", "muon": "muon tracks",
         "all": "all six combined"}
SHAPES = {"deep_underground": "^", "shallow_underground": "s", "surface_lab": "D"}


def relative(rows, preset, mask, label):
    oracle_key = "oracle" if mask != "all" else f"fixed_tuned_on_{preset}"
    oracle = rows[mask][oracle_key]["fom"]
    value = rows[mask][label]["fom"]
    return value / oracle if oracle > 0 else float("nan")


def main():
    ps.apply()
    data = json.loads(INPUT.read_text(encoding="utf-8"))["per_sensor"]
    presets = [p for p in SENSORS if p in data]
    fig, axes = plt.subplots(1, len(presets), figsize=(9.6, 3.8), sharey=True, sharex=True,
                             constrained_layout=True)
    masks = list(MASKS)
    for ax, preset in zip(axes, presets):
        rows = data[preset]
        ax.axvline(1.0, color=ps.INK_SECONDARY, linewidth=1)
        for i, mask in enumerate(masks):
            for source in presets:
                if source == preset:
                    continue
                label = f"transplant_from_{source}" if mask != "all" else f"fixed_tuned_on_{source}"
                ax.plot([relative(rows, preset, mask, label)], [i], **ps.marker(ps.SERIES[1], SHAPES[source], 6))
            ax.plot([relative(rows, preset, mask, "adaptive")], [i], **ps.marker(ps.SERIES[0], "o", 7))
        ax.set_title(SENSORS[preset], loc="left", fontsize=10, color=ps.INK)
        ax.grid(axis="y", visible=False)
        ax.set_xlim(0, 1.3)
    axes[0].set_yticks(range(len(masks)), [MASKS[m] for m in masks])
    axes[0].invert_yaxis()
    fig.supxlabel("figure of merit relative to the oracle (1 = fixed mask tuned with truth on this sensor)",
                  fontsize=9, color=ps.INK_SECONDARY)
    handles = [Line2D([], [], **ps.marker(ps.SERIES[0], "o", 7), label="adaptive (no truth, no retuning)")]
    handles += [Line2D([], [], **ps.marker(ps.SERIES[1], SHAPES[s], 6), label=f"fixed, tuned on {SENSORS[s]}")
                for s in presets]
    fig.legend(handles=handles, loc="outside upper center", ncols=4, fontsize=8)
    fig.savefig(OUTPUT, dpi=200)
    write_sidecar(OUTPUT, __file__, inputs=[INPUT], notes="figure of R4")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
