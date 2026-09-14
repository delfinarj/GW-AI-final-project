"""Figure for R4: each mask's figure of merit relative to the oracle tuned on the same sensor.

Reads every results/compare_masks/**/compare_masks.json (one per independent seed). Small multiples,
one panel per sensor, one row per mask (and the six combined). The oracle is the reference line at 1.
Within a row each series is offset vertically so no marker hides another: blue circles are the
adaptive masks, orange markers the fixed masks transplanted from each of the other two sensors
(marker shape names the source). Markers are the median over seeds, whiskers the min-max range.

Run:  python analysis/figure_compare_masks.py
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from skmask import plotstyle as ps
from skmask.provenance import write_sidecar

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "compare_masks"
OUTPUT = RESULTS / "compare_masks.png"
SENSORS = {"deep_underground": "deep underground", "shallow_underground": "shallow underground",
           "surface_lab": "surface laboratory"}
MASKS = {"hot": "hot columns & pixels", "cti": "charge-transfer trails", "halo": "halo",
         "serial": "serial-register hits", "lec": "low-energy clusters", "muon": "muon tracks",
         "all": "all six combined"}
SHAPES = {"deep_underground": "^", "shallow_underground": "s", "surface_lab": "D"}
SIZES = {"^": 7, "s": 6, "D": 5.5}


def relative(rows, preset, mask, label):
    oracle_key = "oracle" if mask != "all" else f"fixed_tuned_on_{preset}"
    oracle = rows[mask][oracle_key]["fom"]
    return rows[mask][label]["fom"] / oracle if oracle > 0 else np.nan


def series_labels(preset, mask, presets):
    """(key in the JSON, colour, marker shape, legend name) for every series of one row."""
    out = [("adaptive", ps.SERIES[0], "o", "adaptive")]
    for source in presets:
        if source != preset:
            key = f"transplant_from_{source}" if mask != "all" else f"fixed_tuned_on_{source}"
            out.append((key, ps.SERIES[1], SHAPES[source], source))
    return out


def main():
    ps.apply()
    inputs = sorted(RESULTS.glob("**/compare_masks.json"))
    runs = [json.loads(p.read_text(encoding="utf-8"))["per_sensor"] for p in inputs]
    presets = [p for p in SENSORS if p in runs[0]]
    masks = list(MASKS)
    fig, axes = plt.subplots(1, len(presets), figsize=(9.6, 4.2), sharey=True, sharex=True,
                             constrained_layout=True)
    offsets = [-0.22, 0.0, 0.22]
    summary = {}
    for ax, preset in zip(axes, presets):
        ax.axvline(1.0, color=ps.INK_SECONDARY, linewidth=1)
        for i, mask in enumerate(masks):
            for (key, color, shape, name), dy in zip(series_labels(preset, mask, presets), offsets):
                values = np.array([relative(run[preset], preset, mask, key) for run in runs])
                y = i + dy
                if len(values) > 1:
                    ax.plot([np.nanmin(values), np.nanmax(values)], [y, y], color=color, linewidth=2)
                ax.plot([np.nanmedian(values)], [y], **ps.marker(color, shape, SIZES.get(shape, 7)))
                summary[f"{preset}/{mask}/{name}"] = {"median": float(np.nanmedian(values)),
                                                      "min": float(np.nanmin(values)),
                                                      "max": float(np.nanmax(values))}
        ax.set_title(SENSORS[preset], loc="left", fontsize=10, color=ps.INK)
        ax.grid(axis="y", visible=False)
        ax.set_xlim(-0.03, 1.3)
    axes[0].set_yticks(range(len(masks)), [MASKS[m] for m in masks])
    axes[0].set_ylim(len(masks) - 0.5, -0.5)
    fig.supxlabel(f"figure of merit relative to the oracle (1 = fixed mask tuned with truth on this sensor); "
                  f"median and range over {len(runs)} seed{'s' if len(runs) > 1 else ''}",
                  fontsize=9, color=ps.INK_SECONDARY)
    handles = [Line2D([], [], **ps.marker(ps.SERIES[0], "o", 7), label="adaptive (no truth, no retuning)")]
    handles += [Line2D([], [], **ps.marker(ps.SERIES[1], SHAPES[s], SIZES[SHAPES[s]]),
                       label=f"fixed, tuned on {SENSORS[s]}") for s in presets]
    fig.legend(handles=handles, loc="outside upper center", ncols=4, fontsize=8)
    fig.savefig(OUTPUT, dpi=200)
    summary_path = RESULTS / "relative_fom_summary.json"
    summary_path.write_text(json.dumps({"n_seeds": len(runs), "inputs": [str(p.relative_to(ROOT)) for p in inputs],
                                        "relative_fom": summary}, indent=2), encoding="utf-8")
    write_sidecar(OUTPUT, __file__, inputs=inputs, notes=f"figure of R4 over {len(runs)} seeds")
    write_sidecar(summary_path, __file__, inputs=inputs, notes="relative figure of merit, median and range over seeds")
    print(f"wrote {OUTPUT} and {summary_path} from {len(runs)} runs")


if __name__ == "__main__":
    main()
