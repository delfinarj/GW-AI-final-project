"""Everything the short summaries say, computed once and shared by every language they say it in.

The five-page summary exists in English and in Spanish. If each generator read the result files for
itself the two could drift apart with nobody noticing, which is exactly the failure this project is
built to avoid, so both import `collect()` from here: the numbers and the table *rows* are shared,
and only the prose and the table *headers* belong to a language.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results"


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sci(x, digits=1):
    """Scientific notation that stays on one line inside a narrow table column."""
    mantissa, exponent = f"{x:.{digits}e}".split("e")
    return f"<span class='nb'>{mantissa}&times;10<sup>{int(exponent)}</sup></span>"


def collect(sensor_names, mask_names, group_names):
    """Every value the summaries quote. The three dictionaries carry the language's own labels."""
    from skmask.presets import PRESETS

    head = load(RES / "report_numbers.json")
    null = load(RES / "null_false_positive_rates" / "null_rates.json")
    muon = load(RES / "muon_mask_null_rate" / "muon_null.json")
    cross = load(RES / "cross_defect_false_positives" / "cross_defect.json")
    public = load(RES / "adaptive_on_public" / "adaptive_on_public.json")
    rate = load(RES / "release_rate" / "release_rate.json")

    sensor_rows = [[sensor_names[n], f"{s.nx}&times;{s.ny}", f"{s.thickness_um:.0f}", f"{s.noise_e:.2f}",
                    f"{s.exposure_days * 24:.0f} h", sci(s.dark_e_per_pix_day), f"{s.muon_flux_per_cm2_day:.2g}"]
                   for n, s in PRESETS.items()]

    null_cells = [v for rows in null["per_sensor"].values() for v in rows.values() if isinstance(v, dict)]
    null_tested = [v for v in null_cells if v["trials"]]

    fires = []
    for sensor, groups in cross["per_sensor"].items():
        for group, rows in groups.items():
            for mask, v in rows.items():
                if isinstance(v, dict) and not v["own_defect"] and v["trials"] and v["verdict"].startswith("fires"):
                    fires.append((v["rate"], v["median_masked_fraction"] or 0.0, sensor, group, mask, v))
    fires.sort(reverse=True)

    exposures = [public["per_exposure"][k] for k in sorted(public["per_exposure"], key=int)]
    last = exposures[-1]
    loud = last["loudest_columns"]

    return {
        "head": head, "null": null, "muon": muon, "cross": cross, "public": public, "rate": rate,
        "tr": head["transplant"], "ad": head["adaptive"], "held": head["held_out"],
        "sensor_rows": sensor_rows,
        "null_tested": null_tested,
        "null_fired": sum(v["fired"] for v in null_tested),
        "muon_fired": sum(v["fired"] for v in muon["per_sensor"].values()),
        "muon_images": sum(v["images"] for v in muon["per_sensor"].values()),
        "fires": fires,
        "run_survivors": [v for *_, v in fires
                          if v.get("counted_per_run_instead", {}).get("verdict", "").startswith("fires")],
        "fire_rows": [[sensor_names[s], group_names[g], mask_names[m], f"{v['fired']}/{v['trials']}",
                       f"{r:.2f}", f"{f:.4f}"] for r, f, s, g, m, v in fires],
        "worst_fire": fires[0] if fires else None,
        "exposures": exposures, "last": last,
        "public_rows": [[f"{r['exposure_s'] / 3600:.0f} h", r["images"], r["trigger_pixels"]["total"],
                         f"{r['measured_noise_e']['median']:.3f}", sci(r["measured_density_1e"]["median"]),
                         f"{len(r['constants_chosen']['hot_columns'])}", f"{r['constants_chosen']['halo_radius']}",
                         f"{r['masked_fraction_union']:.4f}", f"{r['release_mask_fraction']:.3f}"]
                        for r in exposures],
        "n_loud_flagged": next((i for i, c in enumerate(loud) if not c["flagged_by_us"]), len(loud)),
        "precision": [r["against_release"]["hot_columns_pixels"]["median_fraction_of_ours_also_theirs"]
                      for r in exposures],
        "recall": [r["against_release"]["hot_columns_pixels"]["median_fraction_of_theirs_also_ours"]
                   for r in exposures],
        "ratios": [x for r in exposures for x in (r["hot_column_evidence"]["ratio_to_common"] or [])],
        "r1_noise": [x["noise"] for x in rate["per_exposure"]],
        "noise_medians": [r["measured_noise_e"]["median"] for r in exposures],
        "inputs": [RES / "report_numbers.json", RES / "null_false_positive_rates" / "null_rates.json",
                   RES / "muon_mask_null_rate" / "muon_null.json",
                   RES / "cross_defect_false_positives" / "cross_defect.json",
                   RES / "adaptive_on_public" / "adaptive_on_public.json",
                   RES / "release_rate" / "release_rate.json"],
    }
