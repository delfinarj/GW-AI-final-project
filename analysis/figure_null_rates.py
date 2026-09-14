"""Figure for R3: how often each adaptive mask fires with its target absent, against alpha.

Reads results/null_false_positive_rates/null_rates.json. One series (the measured rate with its
Clopper-Pearson 95 % interval), one reference line at alpha.

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
LABELS = {"hot_columns": "hot columns", "cti": "charge-transfer trails", "halo": "halo",
          "serial": "serial-register hits", "low_energy_clusters": "low-energy clusters"}
ALPHA = 0.01


def main():
    ps.apply()
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    names = [n for n in LABELS if n in data]
    fig, ax = plt.subplots(figsize=(6.4, 2.9), constrained_layout=True)
    for i, name in enumerate(names):
        row = data[name]
        lo, hi = row["ci95"]
        ax.plot([lo, hi], [i, i], color=ps.SERIES[0], linewidth=2, solid_capstyle="round")
        ax.plot([row["rate"]], [i], **ps.marker(ps.SERIES[0]))
        ax.text(hi + 0.002, i, f"{row['fired']}/{row['runs']}", va="center", color=ps.INK_SECONDARY, fontsize=8)
    ax.axvline(ALPHA, color=ps.INK_SECONDARY, linewidth=1)
    ax.text(ALPHA, len(names) - 0.4, f" alpha = {ALPHA}", color=ps.INK_SECONDARY, fontsize=8, va="bottom")
    ax.set_yticks(range(len(names)), [LABELS[n] for n in names])
    ax.invert_yaxis()
    ax.set_ylim(len(names) - 0.3, -0.7)
    ax.set_xlim(0, max(max(data[n]["ci95"][1] for n in names) * 1.25, ALPHA * 3))
    ax.set_xlabel("fraction of defect-free runs in which the mask fired (95 % interval)")
    ax.grid(axis="y", visible=False)
    ax.set_title("Adaptive masks on sensors without their defect", loc="left", fontsize=10, color=ps.INK)
    fig.savefig(OUTPUT, dpi=200)
    write_sidecar(OUTPUT, __file__, inputs=[INPUT], notes="figure of R3")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
