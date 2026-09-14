"""Figure for R1: single-electron density of the public SENSEI release against exposure.

Reads results/release_rate/release_rate.json (written by reproduce_release_rate.py) and draws the
fitted densities with their uncertainties, the weighted straight line, and the rate compared with
the published value.

Run:  python analysis/figure_release_rate.py
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from skmask import plotstyle as ps
from skmask.provenance import write_sidecar

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "release_rate" / "release_rate.json"
OUTPUT = ROOT / "results" / "release_rate" / "release_rate.png"
ROWS_PER_SUPERPIXEL = 32


def main():
    r = json.loads(INPUT.read_text(encoding="utf-8"))
    rows = r["per_exposure"]
    exposure = np.array([row["exposure_days"] for row in rows])
    density = np.array([row["mu"] for row in rows])
    error = np.array([row["mu_err"] for row in rows])
    slope = r["rate_e_per_pix_day"] * ROWS_PER_SUPERPIXEL
    intercept = r["intercept_e_per_superpix"]

    ps.apply()
    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    grid = np.linspace(0, exposure.max() * 1.05, 100)
    ax.plot(grid, (intercept + slope * grid) * 1e4, color=ps.SERIES[1], label="weighted straight line")
    for x, y, err in zip(exposure, density, error):
        ax.plot([x, x], [(y - err) * 1e4, (y + err) * 1e4], color=ps.SERIES[0], linewidth=2)
    ax.plot(exposure, density * 1e4, **ps.marker(ps.SERIES[0]), label="fitted 1e density per exposure")
    ax.set_xlabel("exposure (days)")
    ax.set_ylabel(r"1e density ($10^{-4}$ e$^-$ / superpixel)")
    ax.set_title("Public SENSEI SNOLAB release, reproduced", loc="left", fontsize=10, color=ps.INK)
    text = (f"this work: ({r['rate_e_per_pix_day'] * 1e5:.2f} $\\pm$ {r['rate_err'] * 1e5:.2f})"
            f" $\\times 10^{{-5}}$ e$^-$/pix/day\n"
            f"published: ({r['published_rate'] * 1e5:.2f} $\\pm$ {r['published_err'] * 1e5:.2f})"
            f" $\\times 10^{{-5}}$ (arXiv:2410.18716)")
    ax.text(0.03, 0.97, text, transform=ax.transAxes, va="top", fontsize=9)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax.set_xlim(left=0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(OUTPUT, dpi=200)
    write_sidecar(OUTPUT, __file__, inputs=[INPUT], notes="figure of R1")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
