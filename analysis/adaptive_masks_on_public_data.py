"""R6. What the self-calibrating masks do on a real sensor: the public SENSEI SNOLAB release.

Everything else in this project is measured on the simulator. This runs the adaptive masks, unchanged,
on the public release images and asks two questions a simulation cannot answer:

1. Does the calibration procedure produce sensible constants on a sensor nobody tuned it for? The
   release sensor is a fourth sensor, real, with a geometry none of the presets has: 3200 x 20
   superpixels, each binning 32 physical rows, of which 3072 x 16 are active.
2. Where the release publishes its own mask, do the two agree? The release's masks were made with
   other definitions and other parameters, so disagreement is not by itself an error; what the
   comparison shows is whether the procedure lands on the same pixels a human-tuned mask did.

The bits of the release mask are the hypothesis tested in R2 (`analysis/check_release_mask_bits.py`),
and the pairing of a bit with one of our masks is written here, not derived.

The muon mask is not run: it uses the sensor's thickness, pixel size and back-surface diffusion to
predict the charge and length of a track, and in an image that bins 32 rows into one superpixel a
track is neither straight nor of the predicted length.

Run:  python analysis/adaptive_masks_on_public_data.py
"""
import json
from pathlib import Path

import numpy as np

from skmask import masks as M
from skmask import sensei_public as sp
from skmask.estimate import electrons_from_image
from skmask.provenance import write_sidecar

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "adaptive_on_public"
ALPHA = 0.01
# A superpixel bins 32 physical rows, so a step along y is 32 times a step along x. Distances for the
# halo are measured in units of the column pitch; without this the halo is a disc in superpixels,
# which on a 16-row image reaches the top and bottom of the frame and swallows whole columns.
SAMPLING = (float(sp.ROWS_PER_SUPERPIXEL), 1.0)
# Our mask -> the release bit it should correspond to, under the R2 hypothesis. Written, not derived.
COUNTERPART = {"hot_columns_pixels": 0x400 | 0x200, "cti": 0x4, "halo": 0x8, "serial_rows": 0x20}
BIT_NAMES = {0x400 | 0x200: "bad column | bad pixel", 0x4: "bleeding", 0x8: "halo", 0x20: "noisy row"}


def active_images(exposure_s):
    """Per image of one exposure, the active area as (charge, release mask) arrays."""
    data = sp.load(exposure_s, extra=("RUNID",))
    out = []
    for run in np.unique(data["RUNID"]):
        sel = data["RUNID"] == run
        charge = np.zeros((sp.N_ROWS, sp.N_COLUMNS))
        release = np.zeros((sp.N_ROWS, sp.N_COLUMNS), dtype=np.int64)
        charge[data["y"][sel], data["x"][sel]] = data["ePix"][sel]
        release[data["y"][sel], data["x"][sel]] = data["mask"][sel]
        rows = slice(1, sp.ACTIVE_ROWS + 1)
        columns = slice(0, sp.ACTIVE_COLUMNS)
        out.append((charge[rows, columns], release[rows, columns]))
    return out


def overlap(ours, theirs):
    """How much of each mask the other one also masks."""
    both = int(np.count_nonzero(ours & theirs))
    n_ours, n_theirs = int(np.count_nonzero(ours)), int(np.count_nonzero(theirs))
    return {"ours": n_ours, "theirs": n_theirs, "both": both,
            "fraction_of_ours_also_theirs": (both / n_ours) if n_ours else None,
            "fraction_of_theirs_also_ours": (both / n_theirs) if n_theirs else None}


def run_exposure(exposure_s):
    pairs = active_images(exposure_s)
    stack, release = [], [r for _, r in pairs]
    noise, density = [], []
    for charge, _ in pairs:
        electrons, info = electrons_from_image(charge)
        stack.append(electrons)
        noise.append(float(info["noise_e"]))
        density.append(float(info["density"]))

    length_h, length_v, cti_profiles = M.adaptive_cti_lengths(stack, alpha=ALPHA)
    cti = [M.cti_mask(e, length_h, length_v) for e in stack]
    radius, halo_info = M.adaptive_halo_radius(stack, alpha=ALPHA, exclude=cti, sampling=SAMPLING)
    halo = [M.halo_mask(e, radius, sampling=SAMPLING) for e in stack]
    before_hot = [c | h for c, h in zip(cti, halo)]
    columns = M.adaptive_hot_columns(stack, alpha=ALPHA, exclude=before_hot)
    pixels = M.adaptive_hot_pixels(stack, alpha=ALPHA, exclude=before_hot)
    hot = [M.column_mask(e.shape, columns) | pixels for e in stack]

    # Is a flagged column anomalous on its own terms? Its rate of charged low-charge pixels against
    # the rate of the columns that were not flagged, the same counts the calibration works with.
    charged = np.zeros(stack[0].shape[1])
    valid_pixels = np.zeros_like(charged)
    for electrons, excluded in zip(stack, before_hot):
        charged += (M.low_charge_occupancy(electrons) & ~excluded).sum(axis=0)
        valid_pixels += (~excluded).sum(axis=0)
    flagged = np.zeros(charged.size, dtype=bool)
    flagged[list(columns)] = True
    common_rate = float(charged[~flagged].sum() / max(valid_pixels[~flagged].sum(), 1.0))
    flagged_rates = (charged[flagged] / np.maximum(valid_pixels[flagged], 1.0)).tolist()
    hot_column_evidence = {"common_rate_charged_pixels": common_rate,
                           "flagged_column_rates": [float(r) for r in flagged_rates],
                           "ratio_to_common": [float(r / common_rate) for r in flagged_rates] if common_rate else None}

    serial, lec = [], []
    for e, h, c in zip(stack, hot, cti):
        rows = M.adaptive_serial_rows(e, alpha=ALPHA, exclude=h | c)
        serial.append(M.row_mask(e.shape, rows))
        lec.append(M.adaptive_low_energy_cluster_mask(e, alpha=ALPHA, exclude=h | c | serial[-1]))

    ours = {"hot_columns_pixels": hot, "cti": cti, "halo": halo, "serial_rows": serial,
            "low_energy_clusters": lec}
    fractions = {name: float(np.mean([m.mean() for m in masks])) for name, masks in ours.items()}
    # the masks overlap, so what is actually removed is their union, not the sum of the fractions
    union = [np.logical_or.reduce([masks[i] for masks in ours.values()]) for i in range(len(stack))]
    union_fraction = float(np.mean([u.mean() for u in union]))

    # The loudest columns before any mask is applied, and whether the procedure ended up flagging
    # them: a column can be hidden from the hot-column calibration by an earlier mask.
    raw_charged = np.zeros(stack[0].shape[1])
    for electrons in stack:
        raw_charged += M.low_charge_occupancy(electrons).sum(axis=0)
    raw_rate = raw_charged / (len(stack) * stack[0].shape[0])
    flagged_set = {int(c) for c in columns}
    loudest_columns = [{"column": int(c), "charged_pixel_rate": float(raw_rate[c]),
                        "flagged_by_us": int(c) in flagged_set,
                        "in_release_bad_columns": bool(((release[0][:, c] & COUNTERPART["hot_columns_pixels"]) != 0).any())}
                       for c in np.argsort(raw_rate)[::-1][:8]]
    comparison = {}
    for name, bit in COUNTERPART.items():
        per_image = [overlap(m, (r & bit) != 0) for m, r in zip(ours[name], release)]
        kept = [o["fraction_of_ours_also_theirs"] for o in per_image
                if o["fraction_of_ours_also_theirs"] is not None]
        theirs_seen = [o["fraction_of_theirs_also_ours"] for o in per_image
                       if o["fraction_of_theirs_also_ours"] is not None]
        comparison[name] = {
            "release_bit": BIT_NAMES[bit],
            "median_fraction_of_ours_also_theirs": (float(np.median(kept)) if kept else None),
            "median_fraction_of_theirs_also_ours": (float(np.median(theirs_seen)) if theirs_seen else None),
            "median_pixels_ours": float(np.median([o["ours"] for o in per_image])),
            "median_pixels_theirs": float(np.median([o["theirs"] for o in per_image]))}

    # what the trail calibration had to work with: the release blinds hits, so triggers are few
    block, max_distance = 5, 400
    threshold = ALPHA / (2 * (max_distance // block))
    triggers = [int(np.count_nonzero(e >= M.TRIGGER_E)) for e in stack]
    cti_detail = {}
    for axis in ("h", "v"):
        profile = cti_profiles[axis]
        best = int(np.argmin(profile["p"]))
        cti_detail[axis] = {"smallest_p": float(profile["p"][best]),
                            "at_distance_pix": int((best + 1) * block),
                            "n_downstream": int(profile["n_down"][best]),
                            "n_upstream": int(profile["n_up"][best]),
                            "bonferroni_threshold": threshold}

    return {"exposure_s": exposure_s, "images": len(stack), "shape": list(stack[0].shape),
            "trigger_pixels": {"total": int(sum(triggers)), "median_per_image": float(np.median(triggers))},
            "hot_column_evidence": hot_column_evidence,
            "cti_detail": cti_detail,
            "measured_noise_e": {"median": float(np.median(noise)), "min": min(noise), "max": max(noise)},
            "measured_density_1e": {"median": float(np.median(density)), "min": min(density), "max": max(density)},
            "constants_chosen": {"cti_length_h": int(length_h), "cti_length_v": int(length_v),
                                 "halo_radius": int(radius), "halo_calibrated": bool(halo_info["calibrated"]),
                                 "hot_columns": [int(c) for c in columns],
                                 "hot_pixels": int(np.count_nonzero(pixels))},
            "masked_fraction": fractions,
            "masked_fraction_union": union_fraction,
            "loudest_columns": loudest_columns,
            "release_mask_fraction": float(np.mean([((r & sp.MASK_1E_ANALYSIS) != 0).mean() for r in release])),
            "against_release": comparison}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {"alpha": ALPHA, "per_exposure": {}}
    for exposure_s in sorted(sp.EXPOSURES_S):
        print(f"exposure {exposure_s} s", flush=True)
        r = run_exposure(exposure_s)
        results["per_exposure"][str(exposure_s)] = r
        print(f"  {r['images']} images, noise {r['measured_noise_e']['median']:.3f} e, "
              f"1e density {r['measured_density_1e']['median']:.2e}", flush=True)
        print(f"  chose {r['constants_chosen']}", flush=True)
        ratios = r["hot_column_evidence"]["ratio_to_common"] or []
        if ratios:
            print(f"  flagged columns carry {min(ratios):.0f}-{max(ratios):.0f} times the charged-pixel rate "
                  f"of the columns left alone", flush=True)
        print(f"  {r['trigger_pixels']['total']} trigger pixels in all, median "
              f"{r['trigger_pixels']['median_per_image']:.0f} per image; smallest trail p-value "
              f"h {r['cti_detail']['h']['smallest_p']:.2g}, v {r['cti_detail']['v']['smallest_p']:.2g} "
              f"(needs < {r['cti_detail']['h']['bonferroni_threshold']:.2g})", flush=True)
        missed = [c for c in r["loudest_columns"] if not c["flagged_by_us"] and c["charged_pixel_rate"] > 0.01]
        if missed:
            print(f"  loud columns NOT flagged: {[(c['column'], round(c['charged_pixel_rate'], 3)) for c in missed]}",
                  flush=True)
        for name, c in r["against_release"].items():
            print(f"  {name:<20} masks {r['masked_fraction'][name]:.4f} of the image; of ours "
                  f"{c['median_fraction_of_ours_also_theirs']} is also '{c['release_bit']}'", flush=True)

    output = OUT_DIR / "adaptive_on_public.json"
    output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    write_sidecar(output, __file__,
                  inputs=[sp.DATA_DIR / name for name in sp.EXPOSURES_S.values()],
                  parameters={"alpha": ALPHA, "sampling_rows_columns": list(SAMPLING),
                              "counterpart_bits": {k: hex(v) for k, v in COUNTERPART.items()}},
                  notes="the adaptive masks run unchanged on the public SENSEI release, against its own mask bits")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
