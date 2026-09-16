"""R3. How often does each adaptive mask fire when its target is absent?

Every adaptive mask sets its threshold with a Bonferroni-corrected tail probability alpha, so where
its defect is absent it should fire on a fraction ~alpha of independent trials. This script measures
that fraction on the three presets of the project with every defect switched off (no hot columns or
pixels, no charge-transfer trails, no serial-register hits, no low-energy clusters, no halo), while
keeping what cannot be switched off in a real sensor: dark current, spurious charge, the injected
signal, muon tracks and high-energy deposits. One simulation therefore serves all five masks.

Two rates are reported, because the masks work at different granularities:

- per run (a stack of images, the unit a stack-calibrated mask decides on): hot columns, CTI, halo;
- per image: serial-register rows, low-energy clusters.

The check is two-sided: the Clopper-Pearson 95 % interval must contain alpha. An interval entirely
below alpha means the mask is more conservative than it claims, which is also worth knowing.

Run:  python analysis/null_false_positive_rates.py [n_runs]
"""
import json
import sys
from pathlib import Path

import numpy as np

from skmask import masks as M
from skmask.estimate import electrons_from_image
from skmask.presets import PRESETS
from skmask.provenance import write_sidecar
from skmask.simulate import simulate
from skmask.stats import clopper_pearson

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "null_false_positive_rates"
ALPHA = 0.01
IMAGES_PER_RUN = 2
SEED = 20260915

# Every defect the masks look for, switched off; everything unavoidable kept.
DEFECTS_OFF = dict(n_hot_columns=0, hot_column_e_per_pix=0.0, n_hot_pixels=0, hot_pixel_e=0.0,
                   cti_h_prob=0.0, cti_v_prob=0.0, serial_hits_per_image=0.0,
                   lowE_clusters_per_image=0.0, halo_yield_per_e=0.0)
CLEAN_PRESETS = {name: sensor.with_(**DEFECTS_OFF) for name, sensor in PRESETS.items()}
PER_RUN = ("hot_columns", "cti", "halo")
PER_IMAGE = ("serial", "low_energy_clusters")


def main(n_runs):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    fired = {name: {k: 0 for k in PER_RUN + PER_IMAGE} for name in CLEAN_PRESETS}
    trials = {name: {k: 0 for k in PER_RUN + PER_IMAGE} for name in CLEAN_PRESETS}
    uncalibrated = {name: 0 for name in CLEAN_PRESETS}

    for run in range(n_runs):
        for name, sensor in CLEAN_PRESETS.items():
            stack = [electrons_from_image(simulate(sensor, rng).measured)[0] for _ in range(IMAGES_PER_RUN)]

            length_h, length_v, _ = M.adaptive_cti_lengths(stack, alpha=ALPHA)
            cti = [M.cti_mask(e, length_h, length_v) for e in stack]
            radius, info = M.adaptive_halo_radius(stack, alpha=ALPHA, exclude=cti)
            halo = [c | M.halo_mask(e, radius) for c, e in zip(cti, stack)]
            columns = M.adaptive_hot_columns(stack, alpha=ALPHA, exclude=halo)
            pixels = M.adaptive_hot_pixels(stack, alpha=ALPHA, exclude=halo)
            hot = [M.column_mask(e.shape, columns) | pixels for e in stack]

            fired[name]["cti"] += (length_h > 0) or (length_v > 0)
            fired[name]["hot_columns"] += (len(columns) > 0) or bool(pixels.any())
            if info["calibrated"]:
                fired[name]["halo"] += radius > 0
                trials[name]["halo"] += 1
            else:
                uncalibrated[name] += 1
            trials[name]["cti"] += 1
            trials[name]["hot_columns"] += 1

            for e, h, c in zip(stack, hot, cti):
                serial = M.adaptive_serial_rows(e, alpha=ALPHA, exclude=h | c)
                lec = M.adaptive_low_energy_cluster_mask(
                    e, alpha=ALPHA, exclude=h | c | M.row_mask(e.shape, serial))
                fired[name]["serial"] += len(serial) > 0
                fired[name]["low_energy_clusters"] += bool(lec.any())
                trials[name]["serial"] += 1
                trials[name]["low_energy_clusters"] += 1
        if (run + 1) % 5 == 0:
            print(f"{run + 1} runs done", flush=True)

    results = {"alpha": ALPHA, "n_runs": n_runs, "images_per_run": IMAGES_PER_RUN, "per_sensor": {}}
    for name in CLEAN_PRESETS:
        rows = {}
        for mask in PER_RUN + PER_IMAGE:
            k, n = fired[name][mask], trials[name][mask]
            lo, hi = clopper_pearson(k, n) if n else (float("nan"), float("nan"))
            rows[mask] = {"unit": "run" if mask in PER_RUN else "image", "fired": int(k), "trials": int(n),
                          "rate": (k / n) if n else None, "ci95": [lo, hi],
                          "contains_alpha": bool(lo <= ALPHA <= hi) if n else None,
                          "below_alpha": bool(hi < ALPHA) if n else None}
            print(f"{name:<20} {mask:<20} {k:>3}/{n:<4} rate {rows[mask]['rate'] if n else float('nan'):.3f}  "
                  f"CI [{lo:.3f}, {hi:.3f}]  {'contains' if rows[mask]['contains_alpha'] else ('below' if rows[mask]['below_alpha'] else 'ABOVE')} alpha")
        rows["halo_uncalibrated_runs"] = uncalibrated[name]
        results["per_sensor"][name] = rows

    output = OUT_DIR / "null_rates.json"
    output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    write_sidecar(output, __file__, seed=SEED,
                  parameters={"alpha": ALPHA, "n_runs": n_runs, "images_per_run": IMAGES_PER_RUN,
                              "defects_off": DEFECTS_OFF,
                              "presets": {k: v.__dict__ for k, v in CLEAN_PRESETS.items()}},
                  notes="false-positive rates of the adaptive masks on the project's presets with every target defect switched off")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 50)
