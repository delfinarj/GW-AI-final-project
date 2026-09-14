"""Which bit of the public release's `mask` branch is which mask?

The release README names the eight masks selected by 0x067d (neighbour / 1-pixel cluster cut,
bleeding, high-energy halo, crosstalk, noisy row, edge, bad pixel, bad column) but does not give
the bit of each. The hypothesis tested here is that they appear in bit order:

    0x1 neighbour, 0x4 bleeding, 0x8 halo, 0x10 crosstalk, 0x20 noisy row, 0x40 edge,
    0x200 bad pixel, 0x400 bad column

Each mask leaves a geometric signature that does not depend on its parameters, so the hypothesis
is checked bit by bit against the signature its name predicts:

- bad column: masked superpixels lie in columns masked on (nearly) all active rows;
- noisy row: masked superpixels lie in rows masked on (nearly) all columns;
- edge: masked superpixels are concentrated at the borders of the active area;
- neighbour: every masked superpixel touches a superpixel above 0.5 e-;
- halo: masked superpixels are close to superpixels above 100 e- in every direction;
- bleeding: masked superpixels have a >100 e- superpixel upstream (smaller x in the row or
  smaller y in the column), i.e. the distribution is one-sided;
- crosstalk and bad pixel: small, isolated; not uniquely identifiable by geometry, reported only.

Run:  python analysis/check_release_mask_bits.py
"""
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

from skmask import sensei_public as sp
from skmask.provenance import write_sidecar

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "release_mask_bits"
HYPOTHESIS = {0x1: "neighbour", 0x4: "bleeding", 0x8: "halo", 0x10: "crosstalk",
              0x20: "noisy row", 0x40: "edge", 0x200: "bad pixel", 0x400: "bad column"}
EXPOSURE_S = 21600    # the file in which the most bits occur


def images(data):
    """Split the flat tree into per-image (20, 3200) arrays, keyed by RUNID."""
    out = {}
    for run in np.unique(data["RUNID"]):
        sel = data["RUNID"] == run
        charge = np.zeros((sp.N_ROWS, sp.N_COLUMNS))
        mask = np.zeros((sp.N_ROWS, sp.N_COLUMNS), dtype=np.int64)
        charge[data["y"][sel], data["x"][sel]] = data["ePix"][sel]
        mask[data["y"][sel], data["x"][sel]] = data["mask"][sel]
        out[int(run)] = (charge, mask)
    return out


def signatures(charge, mask, bit):
    active = np.zeros_like(mask, dtype=bool)
    active[1:sp.ACTIVE_ROWS + 1, :sp.ACTIVE_COLUMNS] = True
    masked = ((mask & bit) != 0) & active
    n = int(masked.sum())
    if n == 0:
        return None
    rows_active = sp.ACTIVE_ROWS
    full_columns = masked.sum(axis=0) >= 0.9 * rows_active
    full_rows = masked.sum(axis=1) >= 0.9 * sp.ACTIVE_COLUMNS
    ys, xs = np.nonzero(masked)

    border = np.minimum.reduce([xs, sp.ACTIVE_COLUMNS - 1 - xs, (ys - 1) * sp.ROWS_PER_SUPERPIXEL,
                                (sp.ACTIVE_ROWS - ys) * sp.ROWS_PER_SUPERPIXEL])

    touching = ndimage.maximum_filter(charge * active, size=3) > 0.5
    bright = (charge > 100) & active
    if bright.any():
        # physical distance in pixels: one superpixel is 1 column wide and 32 rows tall
        distance = ndimage.distance_transform_edt(~bright, sampling=(sp.ROWS_PER_SUPERPIXEL, 1))[ys, xs]
        upstream_row = np.maximum.accumulate(bright, axis=1)
        upstream_col = np.maximum.accumulate(bright, axis=0)
        upstream = (upstream_row | upstream_col)[ys, xs]
        downstream_row = np.maximum.accumulate(bright[:, ::-1], axis=1)[:, ::-1]
        downstream_col = np.maximum.accumulate(bright[::-1, :], axis=0)[::-1, :]
        downstream_only = ((downstream_row | downstream_col) & ~(upstream_row | upstream_col))[ys, xs]
    else:
        distance = np.full(n, np.inf)
        upstream = np.zeros(n, dtype=bool)
        downstream_only = np.zeros(n, dtype=bool)

    return {
        "n_masked": n,
        "frac_in_full_columns": float(full_columns[xs].mean()),
        "frac_in_full_rows": float(full_rows[ys].mean()),
        "median_distance_to_border_pix": float(np.median(border)),
        "frac_touching_charge_above_0.5e": float(touching[ys, xs].mean()),
        "median_distance_to_100e_pix": float(np.median(distance)),
        "frac_with_100e_upstream": float(upstream.mean()),
        "frac_with_100e_only_downstream": float(downstream_only.mean()),
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = sp.load(EXPOSURE_S, extra=("RUNID",))
    per_image = images(data)
    summary = {}
    for bit, name in HYPOTHESIS.items():
        rows = [s for s in (signatures(c, m, bit) for c, m in per_image.values()) if s]
        if not rows:
            summary[hex(bit)] = {"hypothesis": name, "images_with_bit": 0}
            continue
        keys = [k for k in rows[0] if k != "n_masked"]
        summary[hex(bit)] = {
            "hypothesis": name,
            "images_with_bit": len(rows),
            "n_masked_total": int(sum(r["n_masked"] for r in rows)),
            **{k: float(np.average([r[k] for r in rows], weights=[r["n_masked"] for r in rows]))
               for k in keys},
        }
        s = summary[hex(bit)]
        print(f"{hex(bit):>6} {name:<11} images {s['images_with_bit']:>2}  masked {s['n_masked_total']:>8}  "
              f"fullcol {s['frac_in_full_columns']:.2f}  fullrow {s['frac_in_full_rows']:.2f}  "
              f"border {s['median_distance_to_border_pix']:>7.1f}  touch>0.5e {s['frac_touching_charge_above_0.5e']:.2f}  "
              f"d100e {s['median_distance_to_100e_pix']:>7.1f}  up {s['frac_with_100e_upstream']:.2f}  "
              f"downonly {s['frac_with_100e_only_downstream']:.2f}")

    output = OUT_DIR / "mask_bit_signatures.json"
    output.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_sidecar(output, __file__, inputs=[sp.DATA_DIR / sp.EXPOSURES_S[EXPOSURE_S]],
                  parameters={"exposure_s": EXPOSURE_S, "hypothesis": {hex(k): v for k, v in HYPOTHESIS.items()},
                              "full_line_fraction": 0.9, "bright_threshold_e": 100, "touch_threshold_e": 0.5},
                  notes="Bit-to-name mapping of the public release, tested by geometric signature")


if __name__ == "__main__":
    main()
