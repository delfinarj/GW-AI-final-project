"""R7. How often does the muon mask fire when there are no tracks at all?

The other false-positive tests cannot ask this. R3 switches every defect off but keeps muon tracks and
high-energy deposits, because they cannot be switched off in a real sensor and because the halo is
generated from their charge; R5 switches the defects on one at a time and keeps tracks on for the same
reason; and on the public release the muon mask is not run, since an image binning 32 physical rows
per superpixel has no tracks of the shape the mask predicts. So the sixth mask has never been asked
the question the other five were asked.

Here every defect is off and, in addition, the muon flux and the high-energy background are zero.
What is left is dark current, spurious charge, the injected signal and read noise. Any pixel the muon
mask flags is a false positive, and unlike the other masks its threshold is not a tail probability but
a physical statement: a cluster long enough, straight enough and carrying about the charge a
minimum-ionising particle leaves crossing the sensor. The rate is therefore not expected to equal
alpha; it is whatever coincidences of dark current can imitate a track, which is what is measured.

Run:  python analysis/muon_mask_null_rate.py [n_runs]
"""
import json
import sys
from pathlib import Path

import numpy as np

from skmask import masks as M
from skmask.estimate import electrons_from_image
from skmask.presets import PRESETS
from skmask.provenance import write_sidecar
from skmask.simulate import diffusion_sigma_um, simulate
from skmask.stats import clopper_pearson

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "muon_mask_null_rate"
IMAGES_PER_RUN = 2
SEED = 20260917
# Everything the masks look for, plus the two sources of tracks, switched off.
ALL_OFF = dict(n_hot_columns=0, hot_column_e_per_pix=0.0, n_hot_pixels=0, hot_pixel_e=0.0,
               cti_h_prob=0.0, cti_v_prob=0.0, serial_hits_per_image=0.0,
               lowE_clusters_per_image=0.0, halo_yield_per_e=0.0,
               muon_flux_per_cm2_day=0.0, highE_dru=0.0)


def main(n_runs):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {"n_runs": n_runs, "images_per_run": IMAGES_PER_RUN, "seed": SEED, "per_sensor": {}}
    for index, (name, sensor) in enumerate(PRESETS.items()):
        quiet = sensor.with_(**ALL_OFF)
        rng = np.random.default_rng([SEED, index])   # not hash(name): that varies between processes
        fired, masked, images = 0, [], 0
        for _ in range(n_runs):
            for _ in range(IMAGES_PER_RUN):
                electrons, _ = electrons_from_image(simulate(quiet, rng).measured)
                mask = M.adaptive_muon_mask(electrons, quiet.pixel_um, quiet.thickness_um,
                                            float(diffusion_sigma_um(quiet, quiet.thickness_um)))
                fired += bool(mask.any())
                masked.append(float(mask.mean()))
                images += 1
        lo, hi = clopper_pearson(fired, images)
        results["per_sensor"][name] = {"images": images, "fired": int(fired),
                                       "rate": fired / images, "ci95": [lo, hi],
                                       "largest_fraction_masked": float(max(masked)),
                                       "mean_fraction_masked": float(np.mean(masked))}
        print(f"{name:<20} {fired:>3}/{images:<4} rate {fired / images:.3f}  CI [{lo:.3f}, {hi:.3f}]  "
              f"largest fraction masked {max(masked):.2e}", flush=True)

    output = OUT_DIR / "muon_null.json"
    output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    write_sidecar(output, __file__, seed=SEED,
                  parameters={"n_runs": n_runs, "images_per_run": IMAGES_PER_RUN, "all_off": ALL_OFF},
                  notes="how often the adaptive muon mask fires on sensors with no tracks at all")
    print(f"wrote {output}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
