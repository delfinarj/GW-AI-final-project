"""Independent arithmetic checks for the report's highest-impact claims.

This script does not run the three long simulations.  It recomputes summaries
from their committed JSON, checks sidecar hashes using the repository's stated
LF-normalized convention, and checks public-data mask bits directly.

Run with:
    uv run python verification/recompute_claims_20260925.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import uproot
from scipy.stats import beta


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
TEXT_SUFFIXES = {".json", ".md", ".py", ".html", ".css", ".txt", ".yml", ".yaml", ".toml", ".cfg", ".C"}
MIN_TARGET = 20
MARGIN = 0.01


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def recorded_hash(path: Path) -> str:
    data = path.read_bytes()
    if path.suffix in TEXT_SUFFIXES:
        data = data.replace(bytes((13, 10)), bytes((10,)))
    return hashlib.sha256(data).hexdigest()


def check_sidecars():
    sidecars = list(RESULTS.glob("**/*.provenance.json")) + list((ROOT / "report").glob("*.provenance.json"))
    output_bad, input_bad = [], []
    for sidecar in sidecars:
        record = read_json(sidecar)
        output = ROOT / record["output"]
        actual = recorded_hash(output)
        if actual != record["output_sha256"]:
            output_bad.append(
                {"sidecar": str(sidecar.relative_to(ROOT)), "output": record["output"],
                 "expected": record["output_sha256"], "actual": actual}
            )
        for item in record.get("inputs", []):
            actual = recorded_hash(ROOT / item["path"])
            if actual != item["sha256"]:
                input_bad.append(
                    {"sidecar": str(sidecar.relative_to(ROOT)), "input": item["path"],
                     "expected": item["sha256"], "actual": actual}
                )
    return {"sidecars": len(sidecars), "output_mismatches": output_bad, "input_mismatches": input_bad}


def frame_rows(seed, sensor, mask):
    forms = seed["per_sensor"][sensor][mask]
    oracle_key = "oracle"
    out = {}
    for name, row in forms.items():
        label = name.replace("transplant_from_", "")
        out[label] = row
    return {"target": forms[oracle_key]["target"], "forms": out}


def headline(indices, seeds):
    cases = {}
    for sensor in seeds[0]["per_sensor"]:
        for mask in seeds[0]["per_sensor"][sensor]:
            if mask == "all":
                continue
            rows = [frame_rows(seeds[i], sensor, mask) for i in indices]
            valid = [row for row in rows if row["target"] >= MIN_TARGET]
            if len(valid) >= (3 if len(indices) == 5 else 1):
                cases[(sensor, mask)] = valid

    transplant_harm = []
    transplant_same = []
    adaptive_harm = []
    adaptive_benefit = []
    adaptive_rel = {}
    adaptive_worst_seed = {}
    transplant_beats = []
    for key, rows in cases.items():
        oracle = np.array([row["forms"]["oracle"]["fom"] for row in rows])
        no_mask = np.array([row["forms"]["no_mask"]["fom"] for row in rows])
        adaptive = np.array([row["forms"]["adaptive"]["fom"] for row in rows])
        adaptive_rel[key] = float(np.median(adaptive / oracle))
        adaptive_worst_seed[key] = float(np.min(adaptive / oracle))
        if np.all(adaptive < no_mask * (1 - MARGIN)):
            adaptive_harm.append(key)
        if np.all(adaptive > no_mask * (1 + MARGIN)):
            adaptive_benefit.append(key)
        sources = [name for name in rows[0]["forms"] if name not in ("oracle", "adaptive", "no_mask")]
        source_medians = {}
        for source in sources:
            values = np.array([row["forms"][source]["fom"] for row in rows])
            source_medians[source] = float(np.median(values))
            if np.all(values < no_mask * (1 - MARGIN)):
                transplant_harm.append((*key, source))
            elif np.all(np.abs(values - no_mask) <= no_mask * MARGIN):
                transplant_same.append((*key, source))
        best = max(source_medians, key=source_medians.get)
        best_values = np.array([row["forms"][best]["fom"] for row in rows])
        if np.all(best_values > adaptive):
            transplant_beats.append((*key, best))
    worst = min(adaptive_rel, key=adaptive_rel.get)
    return {
        "adaptive_cases": len(cases),
        "transplant_cases": sum(
            len([name for name in rows[0]["forms"] if name not in ("oracle", "adaptive", "no_mask")])
            for rows in cases.values()
        ),
        "transplant_harmful": len(transplant_harm),
        "transplant_no_difference": len(transplant_same),
        "adaptive_harmful": len(adaptive_harm),
        "adaptive_better": len(adaptive_benefit),
        "adaptive_median_of_oracle": float(np.median(list(adaptive_rel.values()))),
        "adaptive_worst_case": list(worst),
        "adaptive_worst_median": adaptive_rel[worst],
        "adaptive_worst_seed": adaptive_worst_seed[worst],
        "transplant_beats_adaptive": len(transplant_beats),
    }


def cp_interval(k, n, level=0.95):
    lo = beta.ppf((1 - level) / 2, k, n - k + 1) if k else 0.0
    hi = beta.ppf(1 - (1 - level) / 2, k + 1, n - k) if k < n else 1.0
    return [float(lo), float(hi)]


def check_intervals():
    null = read_json(RESULTS / "null_false_positive_rates" / "null_rates.json")
    max_null_delta = 0.0
    for rows in null["per_sensor"].values():
        for name, row in rows.items():
            if not isinstance(row, dict) or "ci95" not in row:
                continue
            recomputed = cp_interval(row["fired"], row["trials"])
            max_null_delta = max(max_null_delta, *(abs(a - b) for a, b in zip(recomputed, row["ci95"])))

    cross = read_json(RESULTS / "cross_defect_false_positives" / "cross_defect.json")
    level = cross["confidence_of_the_corrected_bound"]
    max_cross_delta = 0.0
    reported_cross_fires = 0
    paired_recheck = []
    for sensor, groups in cross["per_sensor"].items():
        for defect, rows in groups.items():
            if not isinstance(rows, dict):
                continue
            raw = cross["counts_this_was_built_from"][f"{sensor}|{defect}"]
            for mask in ("hot_columns", "cti", "halo", "serial", "low_energy_clusters"):
                row = rows[mask]
                if row["trials"]:
                    lo = cp_interval(row["fired"], row["trials"], level)[0]
                    max_cross_delta = max(max_cross_delta, abs(lo - row["lower_bound_family_corrected"]))
                if not row["own_defect"] and row["verdict"].startswith("fires"):
                    reported_cross_fires += 1
                    if row["unit"] == "image":
                        fired = [value > 0 for value in raw["fraction"][mask]]
                        run_any = sum(any(fired[i : i + 2]) for i in range(0, len(fired), 2))
                        run_lo = cp_interval(run_any, 20, level)[0]
                        paired_recheck.append(
                            {"sensor": sensor, "defect": defect, "mask": mask,
                             "image_fired": row["fired"], "image_trials": row["trials"],
                             "run_any_fired": run_any, "run_trials": 20,
                             "family_corrected_lower": run_lo,
                             "exceeds_union_bound_null_0.02": run_lo > 0.02}
                        )
    return {
        "null_max_abs_interval_difference": max_null_delta,
        "cross_max_abs_family_lower_difference": max_cross_delta,
        "reported_cross_firing_cells": reported_cross_fires,
        "paired_run_level_recheck": paired_recheck,
    }


def check_public_columns():
    public = read_json(RESULTS / "adaptive_on_public" / "adaptive_on_public.json")
    exposures = [public["per_exposure"][key] for key in sorted(public["per_exposure"], key=int)]
    sets = [set(row["constants_chosen"]["hot_columns"]) for row in exposures]
    selected = sorted(sets[-1])
    path = RESULTS.parent / "data" / "public" / "sensei_snolab_binned_1" / "RELEASE_hits_blinded_EXP72000_13.root"
    with uproot.open(path) as root_file:
        arrays = root_file["calPixTree"].arrays(["x", "y", "mask"], library="np")
    active = (arrays["x"] < 3072) & (arrays["y"] > 0) & (arrays["y"] <= 16)
    x, mask = arrays["x"][active], arrays["mask"][active]
    fractions = {column: float(np.mean((mask[x == column] & 0x400) != 0)) for column in selected + [1659]}
    release_columns = [column for column in range(3072) if np.mean((mask[x == column] & 0x400) != 0) >= 0.9]
    return {
        "counts": [len(group) for group in sets],
        "nested": all(sets[i] <= sets[i + 1] for i in range(len(sets) - 1)),
        "selected_20h": selected,
        "bad_column_bit_0x400_fraction": fractions,
        "release_full_bad_columns": len(release_columns),
        "ours_over_release_full_columns": len(selected) / len(release_columns),
    }


def main():
    seed_files = sorted((RESULTS / "compare_masks").glob("seed_*/compare_masks.json"))
    seeds = [read_json(path) for path in seed_files]
    result = {
        "sidecars": check_sidecars(),
        "headline_all_five": headline(range(5), seeds),
        "headline_held_out": headline((3, 4), seeds),
        "intervals": check_intervals(),
        "public_columns": check_public_columns(),
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
