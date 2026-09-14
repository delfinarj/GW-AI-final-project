"""R3. How often does each adaptive mask fire when its target is absent?

Every adaptive mask sets its threshold with a Bonferroni-corrected tail probability alpha, so on
data with no defect it should fire in at most a fraction ~alpha of independent runs. This script
measures that fraction over many independent realisations of defect-free sensors and reports it
with a Clopper-Pearson 95 % interval.

"Fires" means, per run (a stack of images of one sensor):
- hot columns: at least one column flagged;
- CTI: a non-zero length in either direction;
- halo: a non-zero radius (runs where the halo could not be calibrated are counted separately);
- serial register: at least one row flagged in any image of the stack;
- low-energy clusters: any pixel masked in any image of the stack.

Run:  python analysis/null_false_positive_rates.py [n_runs]
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import beta

from skmask import masks as M
from skmask.events import to_electrons
from skmask.provenance import write_sidecar
from skmask.simulate import Sensor, simulate

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "null_false_positive_rates"
ALPHA = 0.01
IMAGES_PER_RUN = 3
SEED = 20260914

BASE = dict(noise_e=0.12, exposure_days=1.0, dark_e_per_pix_day=5e-4, signal_e_per_pix=5e-4)
# Defect-free sensors; the bright deposits give halo and CTI tests something to trigger on.
PLAIN = Sensor(name="plain", nx=800, ny=400, **BASE)
BRIGHT = Sensor(name="bright", nx=800, ny=400, highE_dru=7e4, highE_keV=(2.0, 7.0), **BASE)


def clopper_pearson(k, n, level=0.95):
    lo = beta.ppf((1 - level) / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - (1 - level) / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def electrons_stack(sensor, rng):
    images = [simulate(sensor, rng) for _ in range(IMAGES_PER_RUN)]
    return [to_electrons(im.measured, sensor.noise_e, im.total.mean()) for im in images]


def main(n_runs):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    fired = {"hot_columns": 0, "cti": 0, "halo": 0, "serial": 0, "low_energy_clusters": 0}
    halo_uncalibrated = 0
    for run in range(n_runs):
        plain = electrons_stack(PLAIN, rng)
        bright = electrons_stack(BRIGHT, rng)
        fired["hot_columns"] += len(M.adaptive_hot_columns(plain, alpha=ALPHA)) > 0
        fired["serial"] += any(len(M.adaptive_serial_rows(e, alpha=ALPHA)) for e in plain)
        fired["low_energy_clusters"] += any(M.adaptive_low_energy_cluster_mask(e, alpha=ALPHA).any() for e in plain)
        h, v, _ = M.adaptive_cti_lengths(bright, alpha=ALPHA)
        fired["cti"] += (h > 0) or (v > 0)
        radius, info = M.adaptive_halo_radius(bright, alpha=ALPHA, max_radius=60)
        if not info["calibrated"]:
            halo_uncalibrated += 1
        else:
            fired["halo"] += radius > 0
        if (run + 1) % 10 == 0:
            print(f"{run + 1} runs: {fired}")

    results = {}
    for name, k in fired.items():
        n = n_runs - (halo_uncalibrated if name == "halo" else 0)
        lo, hi = clopper_pearson(k, n)
        results[name] = {"fired": int(k), "runs": int(n), "rate": k / n, "ci95": [lo, hi],
                         "consistent_with_alpha": lo <= ALPHA}
        print(f"{name:<20} fired {k:>3}/{n}  rate {k / n:.3f}  95% CI [{lo:.3f}, {hi:.3f}]  "
              f"{'consistent with' if lo <= ALPHA else 'ABOVE'} alpha = {ALPHA}")
    results["halo_uncalibrated_runs"] = halo_uncalibrated

    output = OUT_DIR / "null_rates.json"
    output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    write_sidecar(output, __file__, parameters={"alpha": ALPHA, "n_runs": n_runs,
                                                "images_per_run": IMAGES_PER_RUN,
                                                "plain": PLAIN.__dict__, "bright": BRIGHT.__dict__},
                  seed=SEED, results={k: v["rate"] for k, v in results.items() if isinstance(v, dict)})


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 200)
