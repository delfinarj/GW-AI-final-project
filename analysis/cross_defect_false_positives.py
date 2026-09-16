"""R5. Does an adaptive mask fire on a defect that is not the one it looks for?

R3 measures how often each adaptive mask fires when *every* defect is switched off, so a mask that
fires on another mask's defect passes it unseen: R3 cannot tell "this mask is well behaved" from
"this mask is firing on something else that is also absent". That gap was visible in R4, where the
adaptive low-energy-cluster mask kept only ~0.85 of the signal on the surface sensor even in a seed
with no low-energy clusters at all.

This script closes it by switching the defects on ONE AT A TIME. For each sensor and each of the five
switchable defects it simulates a sensor with that defect at its preset value and the other four at
zero, runs the same calibration chain as R3, and records which masks fire. A mask that fires when its
own defect is off and another one is on is firing on the wrong thing, and the rate is measured.

What cannot be switched off is kept in every configuration, exactly as in R3: dark current, spurious
charge, the injected signal, muon tracks and high-energy deposits. Muons in particular stay on
because the halo is generated from deposited charge: with no tracks there is no halo to switch on.

The reference for "should not fire" is the alpha the mask is built with: alpha for four of them, and
twice alpha for the hot-column row, which fires when either of two alpha-level procedures does (a
column or a pixel). A mask is called to fire on another defect when a Clopper-Pearson lower bound
lies above that rate. The bound is taken at a confidence corrected for the number of cells tested,
because the same question is asked of every mask in every configuration; without that correction,
with 60 cells at 5 %, about one and a half cells would be declared to fire by chance alone. Both
intervals are written out: the plain 95 % one, and the corrected one the verdict uses.

Run:  python analysis/cross_defect_false_positives.py [n_runs]

Each run draws from its own generator, seeded by (SEED, run index), so a run does not depend on the
runs before it. The result file is rewritten after every completed run and carries the counts it was
built from, so an interrupted run resumes from it and lands on the same numbers as an uninterrupted
one; `n_runs_completed` says how many runs each rate is based on. This matters on the machine these
were produced on, which killed one attempt with a segmentation fault (see PROVENANCE.md).
"""
import json
import sys
from pathlib import Path

import numpy as np

from skmask import masks as M
from skmask.estimate import electrons_from_image
from skmask.presets import PRESETS
from skmask.provenance import STATE_AT_START, write_sidecar
from skmask.simulate import simulate
from skmask.stats import clopper_pearson

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "cross_defect_false_positives"
ALPHA = 0.01
IMAGES_PER_RUN = 2
SEED = 20260916

# The same switches R3 uses to turn every defect off; a configuration restores one group of them.
DEFECTS_OFF = dict(n_hot_columns=0, hot_column_e_per_pix=0.0, n_hot_pixels=0, hot_pixel_e=0.0,
                   cti_h_prob=0.0, cti_v_prob=0.0, serial_hits_per_image=0.0,
                   lowE_clusters_per_image=0.0, halo_yield_per_e=0.0)
GROUPS = {"hot_columns_pixels": ("n_hot_columns", "hot_column_e_per_pix", "n_hot_pixels", "hot_pixel_e"),
          "cti": ("cti_h_prob", "cti_v_prob"),
          "serial": ("serial_hits_per_image",),
          "low_energy_clusters": ("lowE_clusters_per_image",),
          "halo": ("halo_yield_per_e",)}
# Which defect each mask is meant to catch; anything else it fires on is the wrong thing.
MASK_TARGET = {"hot_columns": "hot_columns_pixels", "cti": "cti", "halo": "halo",
               "serial": "serial", "low_energy_clusters": "low_energy_clusters"}
PER_RUN = ("hot_columns", "cti", "halo")
PER_IMAGE = ("serial", "low_energy_clusters")
FAMILY = 0.05          # family-wise error the corrected interval controls, over the cells tested
# The rate each row should not exceed when its own defect is absent. The hot-column row fires when
# either the column procedure or the pixel procedure does, and each carries its own alpha.
NOMINAL = {"hot_columns": 2 * ALPHA, "cti": ALPHA, "halo": ALPHA, "serial": ALPHA,
           "low_energy_clusters": ALPHA}


def configurations():
    """(sensor name, defect group, sensor with only that group switched on)."""
    for name, sensor in PRESETS.items():
        off = sensor.with_(**DEFECTS_OFF)
        for group, keys in GROUPS.items():
            on = {k: getattr(sensor, k) for k in keys}
            if all(not v for v in on.values()):
                continue        # this preset does not have this defect: nothing to switch on
            yield name, group, off.with_(**on)


def one_run(sensor, rng):
    """The R3 calibration chain on one stack; per mask, whether it fired and how much it masked."""
    stack = [electrons_from_image(simulate(sensor, rng).measured)[0] for _ in range(IMAGES_PER_RUN)]

    length_h, length_v, _ = M.adaptive_cti_lengths(stack, alpha=ALPHA)
    cti = [M.cti_mask(e, length_h, length_v) for e in stack]
    radius, info = M.adaptive_halo_radius(stack, alpha=ALPHA, exclude=cti)
    halo = [c | M.halo_mask(e, radius) for c, e in zip(cti, stack)]
    columns = M.adaptive_hot_columns(stack, alpha=ALPHA, exclude=halo)
    pixels = M.adaptive_hot_pixels(stack, alpha=ALPHA, exclude=halo)
    hot = [M.column_mask(e.shape, columns) | pixels for e in stack]

    out = {"cti": ((length_h > 0) or (length_v > 0), float(np.mean([c.mean() for c in cti]))),
           "hot_columns": ((len(columns) > 0) or bool(pixels.any()),
                           float(np.mean([h.mean() for h in hot]))),
           "halo": ((radius > 0, float(np.mean([(h & ~c).mean() for h, c in zip(halo, cti)])))
                    if info["calibrated"] else None)}
    serial_out, lec_out = [], []
    for e, h, c in zip(stack, hot, cti):
        rows = M.adaptive_serial_rows(e, alpha=ALPHA, exclude=h | c)
        row_mask = M.row_mask(e.shape, rows)
        lec = M.adaptive_low_energy_cluster_mask(e, alpha=ALPHA, exclude=h | c | row_mask)
        serial_out.append((len(rows) > 0, float(row_mask.mean())))
        lec_out.append((bool(lec.any()), float(lec.mean())))
    out["serial"] = serial_out
    out["low_energy_clusters"] = lec_out
    return out


def raw_key(name, group):
    return f"{name}|{group}"


def smallest_firing_count(n, level, nominal):
    """The fewest firings at which this cell could be called to fire: its detection floor."""
    for k in range(1, n + 1):
        if clopper_pearson(k, n, level)[0] > nominal:
            return k
    return None


def summarise(cells, n_runs_done, n_runs_target):
    tested = sum(1 for (name, group), cell in cells.items() for mask in PER_RUN + PER_IMAGE
                 if MASK_TARGET[mask] != group and cell["trials"][mask])
    level = 1.0 - FAMILY / max(tested, 1)
    per_sensor = {}
    for (name, group), cell in cells.items():
        rows = {}
        for mask in PER_RUN + PER_IMAGE:
            k, n = cell["fired"][mask], cell["trials"][mask]
            lo, hi = clopper_pearson(k, n) if n else (float("nan"), float("nan"))
            lo_family = clopper_pearson(k, n, level)[0] if n else float("nan")
            own = MASK_TARGET[mask] == group
            fractions = cell["fraction"][mask]
            floor = smallest_firing_count(n, level, NOMINAL[mask]) if n and not own else None
            rows[mask] = {"unit": "run" if mask in PER_RUN else "image",
                          "own_defect": own, "fired": int(k), "trials": int(n),
                          "rate": (k / n) if n else None, "ci95": [lo, hi],
                          "nominal_rate": NOMINAL[mask], "lower_bound_family_corrected": lo_family,
                          "smallest_count_that_would_fire": floor,
                          "at_detection_floor": bool(floor is not None and k == floor),
                          "median_masked_fraction": float(np.median(fractions)) if fractions else None,
                          "verdict": ("its own defect" if own else
                                      ("fires on this other defect" if n and lo_family > NOMINAL[mask] else
                                       "consistent with alpha" if n else None))}
        rows["halo_uncalibrated_runs"] = cell["uncalibrated"]
        per_sensor.setdefault(name, {})[group] = rows
    raw = {raw_key(n, g): {"fired": c["fired"], "trials": c["trials"], "fraction": c["fraction"],
                           "uncalibrated": c["uncalibrated"]} for (n, g), c in cells.items()}
    return {"alpha": ALPHA, "family_wise_error": FAMILY, "cells_tested": tested,
            "confidence_of_the_corrected_bound": level, "nominal_rate_per_mask": NOMINAL,
            "n_runs_completed": n_runs_done, "n_runs_target": n_runs_target,
            "images_per_run": IMAGES_PER_RUN, "seed": SEED,
            "images_of_a_run_are_not_independent_trials": True,
            "per_sensor": per_sensor, "counts_this_was_built_from": raw}


def sidecar_parameters(n_runs):
    return {"alpha": ALPHA, "family_wise_error": FAMILY, "n_runs": n_runs,
            "images_per_run": IMAGES_PER_RUN, "defects_off": DEFECTS_OFF,
            "groups": {k: list(v) for k, v in GROUPS.items()}, "mask_target": MASK_TARGET,
            "nominal_rate_per_mask": NOMINAL}


def main(n_runs):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    configs = list(configurations())
    print(f"{len(configs)} configurations: " + ", ".join(f"{n}/{g}" for n, g, _ in configs), flush=True)
    cells = {(n, g): {"fired": {m: 0 for m in PER_RUN + PER_IMAGE},
                      "trials": {m: 0 for m in PER_RUN + PER_IMAGE},
                      "fraction": {m: [] for m in PER_RUN + PER_IMAGE},
                      "uncalibrated": 0} for n, g, _ in configs}
    output = OUT_DIR / "cross_defect.json"
    first_run = 0
    segments = []

    if output.exists():                       # continue an interrupted run, run by run
        previous = json.loads(output.read_text(encoding="utf-8"))
        counts = previous.get("counts_this_was_built_from", {})
        same = (previous.get("seed") == SEED and previous.get("alpha") == ALPHA
                and previous.get("images_per_run") == IMAGES_PER_RUN
                and set(counts) == {raw_key(n, g) for n, g, _ in configs})
        segments = previous.get("code_commits", [])
        if same and previous.get("n_runs_completed", 0) >= n_runs:
            print(f"{previous['n_runs_completed']} runs already done; rewriting the summary from the "
                  f"counts they left", flush=True)
            for (name, group), cell in cells.items():
                cell.update(counts[raw_key(name, group)])
            results = summarise(cells, previous["n_runs_completed"], previous["n_runs_completed"])
            results["code_commits"] = segments
            output.write_text(json.dumps(results, indent=2), encoding="utf-8")
            write_sidecar(output, __file__, seed=SEED, parameters=sidecar_parameters(n_runs),
                          notes="firing rates of the adaptive masks with one defect switched on at a time")
            return
        if same and previous["n_runs_completed"] < n_runs:
            for (name, group), cell in cells.items():
                cell.update(counts[raw_key(name, group)])
            first_run = previous["n_runs_completed"]
            print(f"resuming after {first_run} completed runs", flush=True)
        else:
            print("the file present was made with other settings; starting over", flush=True)

    # every stretch of running records the code it ran under, so a resumed file cannot claim one commit
    commits = segments + [{"from_run": first_run, "commit": STATE_AT_START["commit"],
                           "dirty": STATE_AT_START["dirty"]}]

    for run in range(first_run, n_runs):
        rng = np.random.default_rng([SEED, run])   # this run does not depend on the runs before it
        for name, group, sensor in configs:
            cell, res = cells[(name, group)], one_run(sensor, rng)
            for mask in PER_RUN:
                if res[mask] is None:               # the halo could not be calibrated in this run
                    cell["uncalibrated"] += 1
                    continue
                fired, fraction = res[mask]
                cell["fired"][mask] += bool(fired)
                cell["trials"][mask] += 1
                cell["fraction"][mask].append(fraction)
            for mask in PER_IMAGE:
                for fired, fraction in res[mask]:
                    cell["fired"][mask] += bool(fired)
                    cell["trials"][mask] += 1
                    cell["fraction"][mask].append(fraction)
        partial = summarise(cells, run + 1, n_runs)
        partial["code_commits"] = commits
        output.write_text(json.dumps(partial, indent=2), encoding="utf-8")
        write_sidecar(output, __file__, seed=SEED, parameters=sidecar_parameters(n_runs),
                      notes="firing rates of the adaptive masks with one defect switched on at a time "
                            f"({run + 1} of {n_runs} runs)")
        print(f"{run + 1}/{n_runs} runs done", flush=True)

    results = summarise(cells, n_runs, n_runs)
    results["code_commits"] = commits
    output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    for name, groups in results["per_sensor"].items():
        for group, rows in groups.items():
            for mask in PER_RUN + PER_IMAGE:
                r = rows[mask]
                if not r["own_defect"] and r["trials"]:
                    print(f"{name:<20} only {group:<20} {mask:<20} {r['fired']:>3}/{r['trials']:<4} "
                          f"rate {r['rate']:.3f}  {r['verdict']}", flush=True)
    write_sidecar(output, __file__, seed=SEED, parameters=sidecar_parameters(n_runs),
                  notes="firing rates of the adaptive masks with one defect switched on at a time")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
