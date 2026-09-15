"""R4. Do the masks work on sensors they were not tuned on?

For each of the three presets (src/skmask/presets.py) two independent stacks are simulated: a
calibration stack and a test stack. Every mask is then evaluated on the test stack in three forms:

- ORACLE: the fixed mask with parameters chosen on the calibration stack of the same sensor using
  the simulator's truth (the best a hand-tuned mask can do on that sensor);
- TRANSPLANT: the fixed mask with the oracle parameters of a *different* sensor, unchanged;
- ADAPTIVE: calibrated on the calibration stack of the same sensor without any truth.

Figure of merit for choosing oracle parameters (CHOICE): S / sqrt(S + B), with S the surviving
injected signal events and B the surviving single-electron events of the mask's target origin
(pixels of muon origin for the muon mask). It rewards removing the target while keeping exposure,
and it is the significance of a uniform signal over that background.

A fourth form, NO MASK, is evaluated as well: the figure of merit of doing nothing. It is the level
any mask must beat to be worth applying, and it separates "this transplanted mask does harm" from
"this transplanted mask simply does not help".

Reported per sensor, per mask, per form: target removed, clean pixels kept, signal efficiency, FoM;
and the same for the six masks combined.

Run:  python analysis/compare_masks_across_sensors.py [n_images_per_stack]
"""
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

from skmask import masks as M
from skmask.estimate import electrons_from_image
from skmask.events import clusters, pixel_origin, single_electron_events
from skmask.presets import PRESETS
from skmask.provenance import write_sidecar
from skmask.simulate import diffusion_sigma_um, simulate

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "compare_masks"
SEED = 20260915
ALPHA = 0.01

TARGETS = {"hot": ("hot_column", "hot_pixel"), "cti": ("cti",), "halo": ("halo",),
           "serial": ("serial",), "lec": ("lowE_cluster",), "muon": ("muon",)}

GRIDS = {
    "hot": [dict(factor=f, min_images=m) for f, m in itertools.product([1.2, 1.5, 2, 3, 5, 8, 12], [1, 2, 3, 4, 99])],
    "cti": [dict(length_h=h, length_v=v)
            for h, v in itertools.product([0, 10, 25, 50, 100, 200, 300, 400], [0, 5, 10, 25, 50, 100, 200])],
    "halo": [dict(radius=r) for r in [0, 5, 10, 20, 30, 45, 60, 90, 120, 160, 200, 260]],
    "serial": [dict(window=w, min_charged=m)
               for w, m in itertools.product([5, 10, 20, 50, 100], [2, 3, 4, 5, 6, 8, 12]) if m <= w],
    "lec": [dict(radius=r, min_neighbours=m)
            for r, m in itertools.product([2, 5, 10, 20, 40, 80, 120], [1, 2, 3, 4, 6, 10, 15])],
    "muon": [dict(min_charge=q, min_length=l, dilation=d)
             for q, l, d in itertools.product([50, 100, 500, 2000, 10000, 50000], [2, 3, 10, 30, 100], [0, 2, 4, 8])],
}


class Frame:
    """What an evaluation needs from one simulated image, without the per-origin charge maps."""

    def __init__(self, image):
        sensor = image.sensor
        self.electrons, self.charge_fit = electrons_from_image(image.measured)
        events = single_electron_events(self.electrons)
        origin = pixel_origin(image.charge)
        self.events_of = {name: np.flatnonzero(events & (origin == name)) for name in image.charge}
        self.muon_pixels = np.flatnonzero(image.charge["muon"] > 0)
        background = sum(image.charge[n] for n in ("hot_column", "hot_pixel", "muon", "highE", "halo",
                                                   "serial", "lowE_cluster", "cti"))
        self.clean = np.flatnonzero(background == 0)
        self.labels_table = clusters(self.electrons)
        self.shape = self.electrons.shape


def simulate_frames(sensor, n, rng):
    return [Frame(simulate(sensor, rng)) for _ in range(n)]


def evaluate(masks, frames, mask_name):
    """Aggregate metrics of one mask per frame over a stack."""
    target_total = target_removed = signal_total = signal_kept = clean_total = clean_kept = 0
    for mask, frame in zip(masks, frames):
        flat = mask.ravel()
        if mask_name == "muon":
            idx = frame.muon_pixels
        elif mask_name == "all":
            idx = np.concatenate([frame.events_of[o] for t in TARGETS.values() for o in t if o != "muon"])
        else:
            idx = np.concatenate([frame.events_of[o] for o in TARGETS[mask_name]])
        target_total += len(idx)
        target_removed += flat[idx].sum()
        s = frame.events_of["signal"]
        signal_total += len(s)
        signal_kept += (~flat[s]).sum()
        clean_total += len(frame.clean)
        clean_kept += (~flat[frame.clean]).sum()
    surviving = target_total - target_removed
    return {"target": int(target_total),
            "target_removed": float(target_removed / target_total) if target_total else None,
            "clean_kept": float(clean_kept / clean_total),
            "signal_efficiency": float(signal_kept / signal_total) if signal_total else None,
            "fom": float(signal_kept / np.sqrt(signal_kept + surviving)) if signal_kept + surviving else 0.0}


def fixed_masks(name, params, frames):
    stack = [f.electrons for f in frames]
    if name == "hot":
        columns = M.fixed_hot_columns(stack, params["factor"])
        pixels = M.fixed_hot_pixels(stack, params["min_images"])
        return [M.column_mask(f.shape, columns) | pixels for f in frames]
    if name == "cti":
        return [M.cti_mask(f.electrons, params["length_h"], params["length_v"]) for f in frames]
    if name == "halo":
        return [M.halo_mask(f.electrons, params["radius"]) for f in frames]
    if name == "serial":
        return [M.row_mask(f.shape, M.fixed_serial_rows(f.electrons, params["window"], params["min_charged"]))
                for f in frames]
    if name == "lec":
        return [M.fixed_clustered_events_mask(f.electrons, params["radius"], params["min_neighbours"]) for f in frames]
    if name == "muon":
        return [M.fixed_muon_mask(f.electrons, params["min_charge"], params["min_length"], params["dilation"],
                                  labels_table=f.labels_table) for f in frames]
    raise ValueError(name)


def adaptive_calibration(sensor, frames):
    """Everything an adaptive mask learns from a calibration stack, without truth."""
    # order: CTI -> halo -> hot columns/pixels (see the module docstring of skmask.masks)
    stack = [f.electrons for f in frames]
    length_h, length_v, _ = M.adaptive_cti_lengths(stack, alpha=ALPHA)
    cti = [M.cti_mask(f.electrons, length_h, length_v) for f in frames]
    radius, info = M.adaptive_halo_radius(stack, alpha=ALPHA, exclude=cti)
    halo = [c | M.halo_mask(f.electrons, radius) for c, f in zip(cti, frames)]
    columns = M.adaptive_hot_columns(stack, alpha=ALPHA, exclude=halo)
    pixels = M.adaptive_hot_pixels(stack, alpha=ALPHA, exclude=halo)
    return {"hot_columns": columns.tolist(), "hot_pixels": pixels, "length_h": length_h, "length_v": length_v,
            "halo_radius": radius, "halo_calibrated": bool(info["calibrated"]),
            "halo_at_limit": bool(info["at_limit"]),
            "sigma_back_um": float(diffusion_sigma_um(sensor, sensor.thickness_um))}


def adaptive_masks(sensor, calibration, frames):
    out = {name: [] for name in TARGETS}
    for f in frames:
        hot = M.column_mask(f.shape, calibration["hot_columns"]) | calibration["hot_pixels"]
        cti = M.cti_mask(f.electrons, calibration["length_h"], calibration["length_v"])
        halo = M.halo_mask(f.electrons, calibration["halo_radius"])
        serial = M.row_mask(f.shape, M.adaptive_serial_rows(f.electrons, alpha=ALPHA, exclude=hot | cti | halo))
        lec = M.adaptive_low_energy_cluster_mask(f.electrons, alpha=ALPHA, exclude=hot | cti | halo | serial)
        muon = M.adaptive_muon_mask(f.electrons, sensor.pixel_um, sensor.thickness_um,
                                    calibration["sigma_back_um"], labels_table=f.labels_table)
        for name, m in zip(TARGETS, (hot, cti, halo, serial, lec, muon)):
            out[name].append(m)
    return out


def main(n_images, seed=SEED, out_dir=OUT_DIR):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    started = time.time()
    frames = {}
    for name, sensor in PRESETS.items():
        frames[name] = {"cal": simulate_frames(sensor, n_images, rng), "test": simulate_frames(sensor, n_images, rng)}
        print(f"simulated {name} ({time.time() - started:.0f} s)")

    oracle = {}
    for preset in PRESETS:
        oracle[preset] = {}
        for mask_name, grid in GRIDS.items():
            scores = [evaluate(fixed_masks(mask_name, p, frames[preset]["cal"]), frames[preset]["cal"], mask_name)["fom"]
                      for p in grid]
            oracle[preset][mask_name] = grid[int(np.argmax(scores))]
        print(f"oracle {preset}: {oracle[preset]} ({time.time() - started:.0f} s)")

    results = {"oracle_parameters": oracle, "adaptive_calibration": {}, "per_sensor": {}}
    for preset, sensor in PRESETS.items():
        test = frames[preset]["test"]
        calibration = adaptive_calibration(sensor, frames[preset]["cal"])
        adaptive = adaptive_masks(sensor, calibration, test)
        results["adaptive_calibration"][preset] = {k: v for k, v in calibration.items() if k != "hot_pixels"} | \
            {"n_hot_pixels": int(calibration["hot_pixels"].sum())}
        rows = {}
        combined = {"adaptive": [np.zeros(f.shape, dtype=bool) for f in test]}
        for tuned_on in PRESETS:
            combined[f"fixed_tuned_on_{tuned_on}"] = [np.zeros(f.shape, dtype=bool) for f in test]
        empty = [np.zeros(f.shape, dtype=bool) for f in test]
        for mask_name in TARGETS:
            rows[mask_name] = {"adaptive": evaluate(adaptive[mask_name], test, mask_name),
                               "no_mask": evaluate(empty, test, mask_name)}
            combined["adaptive"] = [c | m for c, m in zip(combined["adaptive"], adaptive[mask_name])]
            for tuned_on in PRESETS:
                masks = fixed_masks(mask_name, oracle[tuned_on][mask_name], test)
                label = "oracle" if tuned_on == preset else f"transplant_from_{tuned_on}"
                rows[mask_name][label] = evaluate(masks, test, mask_name)
                key = f"fixed_tuned_on_{tuned_on}"
                combined[key] = [c | m for c, m in zip(combined[key], masks)]
        rows["all"] = {label: evaluate(m, test, "all") for label, m in combined.items()}
        rows["all"]["no_mask"] = evaluate(empty, test, "all")
        results["per_sensor"][preset] = rows
        print(f"\n=== {preset} ===")
        for mask_name, forms in rows.items():
            for label, r in forms.items():
                removed = "   -  " if r["target_removed"] is None else f"{r['target_removed']:.3f}"
                print(f"{mask_name:<7} {label:<34} target {r['target']:>6} removed {removed}  clean kept "
                      f"{r['clean_kept']:.3f}  signal eff {r['signal_efficiency']:.3f}  FoM {r['fom']:.1f}")

    output = out_dir / "compare_masks.json"
    output.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
    write_sidecar(output, __file__, parameters={"n_images_per_stack": n_images, "alpha": ALPHA, "grids": GRIDS,
                                                "presets": {k: v.__dict__ for k, v in PRESETS.items()}},
                  seed=seed, notes="oracle = fixed mask tuned with truth on the same sensor's calibration stack")
    print(f"done in {time.time() - started:.0f} s")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("n_images", nargs="?", type=int, default=4, help="images per calibration and test stack")
    parser.add_argument("--seed", type=int, default=SEED, help="random seed of the whole run")
    parser.add_argument("--out", default=str(OUT_DIR), help="output directory")
    args = parser.parse_args()
    main(args.n_images, args.seed, args.out)
