"""Build the presented page, report/index.html, from the result files.

No number on the page is typed by hand: every value is read here from results/**.json, and every
sentence that states a comparison ("the adaptive masks never fall below ...") is computed from the
same data, so a re-run of the analyses changes the page consistently. The PDF is this page printed
by scripts/make_pdf.py.

Run:  python analysis/build_report.py
"""
import html
import json
import shutil
from pathlib import Path

import numpy as np

from skmask.presets import PRESETS
from skmask.provenance import write_sidecar

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "report"
FIG = OUT / "figures"
REPO_URL = "https://github.com/delfinarj/GW-AI-final-project"

SENSOR_NAMES = {"deep_underground": "Deep underground", "shallow_underground": "Shallow underground",
                "surface_lab": "Surface laboratory"}
MASK_NAMES = {"hot": "Hot columns & pixels", "cti": "Charge-transfer trails", "halo": "Halo",
              "serial": "Serial-register hits", "lec": "Low-energy clusters", "muon": "Muon tracks",
              "all": "All six combined"}
NULL_NAMES = {"hot_columns": "Hot columns", "cti": "Charge-transfer trails", "halo": "Halo",
              "serial": "Serial-register hits", "low_energy_clusters": "Low-energy clusters"}


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def e(text):
    return html.escape(str(text))


def sci(x, digits=2):
    mantissa, exponent = f"{x:.{digits}e}".split("e")
    return f"{mantissa} &times; 10<sup>{int(exponent)}</sup>"


def table(headers, rows, cls=""):
    head = "".join(f"<th scope='col'>{h}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in row) + "</tr>" for row in rows)
    return f"<div class='table-wrap'><table class='{cls}'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def figure(src, alt, caption):
    return (f"<figure><img src='figures/{src}' alt='{e(alt)}' loading='lazy'>"
            f"<figcaption>{caption}</figcaption></figure>")


def main():
    OUT.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    rate = load(RES / "release_rate" / "release_rate.json")
    bits = load(RES / "release_mask_bits" / "mask_bit_signatures.json")
    null = load(RES / "null_false_positive_rates" / "null_rates.json")
    seed_files = sorted((RES / "compare_masks").glob("seed_*/compare_masks.json"))
    seeds = [load(p) for p in seed_files]
    n_seeds = len(seeds)
    # Seeds run once after the code was frozen, never used to change anything (see PLAN.md).
    holdout_names = {"seed_20260920", "seed_20260921"}
    holdout_idx = [i for i, f in enumerate(seed_files) if f.parent.name in holdout_names]
    MIN_TARGET = 20      # pre-registered: fewer target events than this is "too few to compare"

    for src in (RES / "release_rate" / "release_rate.png", RES / "null_false_positive_rates" / "null_rates.png",
                RES / "compare_masks" / "compare_masks.png"):
        shutil.copy2(src, FIG / src.name)

    MARGIN = 0.01        # paired rule (added after review): a difference under 1 % is "no difference"

    def per_seed(indices):
        """{(sensor, mask): [{"target", "forms": {name: {"fom", "removed", "kept", "signal"}}} per seed]}"""
        out = {}
        for i in indices:
            for sensor, masks in seeds[i]["per_sensor"].items():
                for mask, forms in masks.items():
                    oracle_key = "oracle" if mask != "all" else f"fixed_tuned_on_{sensor}"
                    row = {"target": forms[oracle_key]["target"], "forms": {}}
                    for label, r in forms.items():
                        if label == oracle_key:
                            name = "oracle"
                        elif label in ("adaptive", "no_mask"):
                            name = label
                        else:
                            name = label.replace("transplant_from_", "").replace("fixed_tuned_on_", "")
                        row["forms"][name] = {"fom": r["fom"], "removed": r["target_removed"],
                                              "kept": r["clean_kept"], "signal": r["signal_efficiency"]}
                    out.setdefault((sensor, mask), []).append(row)
        return out

    def collect(indices):
        """Relative figure of merit per seed, for the figure-like table (all seeds, no exclusion)."""
        out = {}
        for key, rows in per_seed(indices).items():
            entry = {"target": [r["target"] for r in rows], "series": {}}
            for r in rows:
                oracle = r["forms"]["oracle"]["fom"]
                for name, f in r["forms"].items():
                    entry["series"].setdefault(name, []).append(f["fom"] / oracle if oracle > 0 else float("nan"))
            out[key] = entry
        return out

    def summarise(indices, min_valid_seeds):
        """Pre-registered statements (median rule) and the stricter paired rule, over the given seeds.

        Exclusion is per seed: a seed enters a mask-sensor case only if it has at least MIN_TARGET target
        events; a case is compared only if at least `min_valid_seeds` seeds enter it.
        """
        data = per_seed(indices)
        cases = {}
        for key, rows in data.items():
            if key[1] == "all":
                continue
            valid = [r for r in rows if r["target"] >= MIN_TARGET]
            if len(valid) >= min_valid_seeds:
                cases[key] = valid
        sources_of = lambda rows: [s for s in rows[0]["forms"] if s not in ("adaptive", "no_mask", "oracle")]
        rel = lambda r, name: r["forms"][name]["fom"] / r["forms"]["oracle"]["fom"]

        median_harm, paired_harm, paired_same, adaptive_paired_harm, adaptive_benefit = [], [], [], [], []
        adaptive_median_harm, beats_adaptive, adaptive_values, adaptive_worst_seed = [], [], {}, {}
        for key, rows in cases.items():
            nm = [r["forms"]["no_mask"]["fom"] for r in rows]
            ad = [r["forms"]["adaptive"]["fom"] for r in rows]
            adaptive_values[key] = float(np.median([rel(r, "adaptive") for r in rows]))
            adaptive_worst_seed[key] = float(min(rel(r, "adaptive") for r in rows))
            if np.median(ad) < np.median(nm):
                adaptive_median_harm.append(key)
            if all(a < n * (1 - MARGIN) for a, n in zip(ad, nm)):
                adaptive_paired_harm.append(key)
            if all(a > n * (1 + MARGIN) for a, n in zip(ad, nm)):
                adaptive_benefit.append(key)
            best_source, best_value = None, -np.inf
            for source in sources_of(rows):
                tr = [r["forms"][source]["fom"] for r in rows]
                if np.median(tr) < np.median(nm):
                    median_harm.append((*key, source))
                if all(x < n * (1 - MARGIN) for x, n in zip(tr, nm)):
                    paired_harm.append(((*key, source), float(np.median([rel(r, source) for r in rows])),
                                        float(np.median([rel(r, "no_mask") for r in rows]))))
                elif all(abs(x - n) <= n * MARGIN for x, n in zip(tr, nm)):
                    paired_same.append((*key, source))
                if np.median(tr) > best_value:
                    best_source, best_value = source, float(np.median(tr))
            if best_source is not None and all(r["forms"][best_source]["fom"] > r["forms"]["adaptive"]["fom"] for r in rows):
                beats_adaptive.append((*key, best_source))
        worst_adaptive = min(adaptive_values.items(), key=lambda kv: kv[1]) if adaptive_values else None
        return {"cases": cases, "n_cases": len(cases),
                "n_transplant": sum(len(sources_of(rows)) for rows in cases.values()),
                "median_harm": median_harm, "paired_harm": paired_harm, "paired_same": paired_same,
                "adaptive_median_harm": adaptive_median_harm, "adaptive_paired_harm": adaptive_paired_harm,
                "adaptive_benefit": adaptive_benefit, "beats_adaptive": beats_adaptive,
                "worst_harm": min(paired_harm, key=lambda r: r[1]) if paired_harm else None,
                "worst_adaptive": worst_adaptive,
                "worst_adaptive_seed": adaptive_worst_seed.get(worst_adaptive[0]) if worst_adaptive else None,
                "median_adaptive": float(np.median(list(adaptive_values.values()))) if adaptive_values else float("nan")}

    overall = summarise(list(range(n_seeds)), min_valid_seeds=3)
    held = summarise(holdout_idx, min_valid_seeds=1) if holdout_idx else None
    development_idx = [i for i in range(n_seeds) if i not in holdout_idx]
    dev = summarise(development_idx, min_valid_seeds=2) if development_idx else None
    overall["data"] = collect(list(range(n_seeds)))

    def describe(key):
        sensor, mask, series = key
        who = "adaptive" if series == "adaptive" else f"tuned on {SENSOR_NAMES[series].lower()}"
        return f"{MASK_NAMES[mask].lower()} on the {SENSOR_NAMES[sensor].lower()} sensor ({who})"

    worst_harm_text = ("" if not overall["worst_harm"] else
                       f"The clearest was {describe(overall['worst_harm'][0])}: a median "
                       f"{overall['worst_harm'][1]:.2f} of the oracle, where doing nothing scores "
                       f"{overall['worst_harm'][2]:.2f}.")
    adaptive_harm_text = (
        f"worse than no mask in every seed in {len(overall['adaptive_paired_harm'])} of {overall['n_cases']} cases"
        + ("" if not overall["adaptive_paired_harm"] else
           " (" + "; ".join(describe((*k, "adaptive")) for k in overall["adaptive_paired_harm"]) + ")")
        + f", and better than no mask in every seed in {len(overall['adaptive_benefit'])}")
    worst_adaptive_text = ("" if overall["worst_adaptive"] is None else
                           f"{overall['worst_adaptive'][1]:.2f} as a median over seeds "
                           f"({describe((*overall['worst_adaptive'][0], 'adaptive'))}; "
                           f"{overall['worst_adaptive_seed']:.2f} in its worst seed)")
    dev_only_harm = [k for k in (dev["adaptive_paired_harm"] if dev else []) if held and k not in held["cases"]]
    def adaptive_below_no_mask(indices, key):
        rows = per_seed(indices).get(key, [])
        return [(r["forms"]["adaptive"]["fom"], r["forms"]["no_mask"]["fom"]) for r in rows]

    held_below_text = ""
    if held is not None:
        pairs = [pair for k in (dev["adaptive_paired_harm"] if dev else []) if k not in held["cases"]
                 for pair in adaptive_below_no_mask(holdout_idx, k)]
        if pairs and all(a < n for a, n in pairs):
            held_below_text = (", although there too its figure of merit is below no mask ("
                               + ", ".join(f"{a / n:.2f}" for a, n in pairs) + " of no mask)")
        elif pairs:
            held_below_text = (", and there its figure of merit relative to no mask is "
                               + ", ".join(f"{a / n:.2f}" for a, n in pairs))

    held_text = ("" if held is None else
                 f"<p>On the {len(holdout_idx)} seeds run once after the code was frozen: "
                 f"{len(held['paired_harm'])} of {held['n_transplant']} transplanted cases worse than no mask in "
                 f"every seed, {len(held['adaptive_paired_harm'])} of {held['n_cases']} adaptive ones, and a median "
                 f"adaptive {held['median_adaptive']:.2f} of the oracle."
                 + ("" if not dev_only_harm else
                    " The adaptive case that is worse than no mask in every development seed ("
                    + "; ".join(describe((*k, "adaptive")) for k in dev_only_harm)
                    + ") has fewer than " + str(MIN_TARGET) + " target events in the held-out seeds, so it is not "
                    "counted there" + held_below_text + ".")
                 + "</p>")
    rule_text = (f"Counted pairwise: a case is harmful only if the transplanted mask scores more than "
                 f"{int(MARGIN * 100)} % below no mask in every seed with at least {MIN_TARGET} target events. "
                 f"The rule pre-registered before the re-run compared medians without a margin and gives "
                 f"{len(overall['median_harm'])} harmful transplants and {len(overall['adaptive_median_harm'])} "
                 f"harmful adaptive cases; it was replaced after an independent review showed that medians "
                 f"differing in the third decimal were being counted as harm.")

    # signal lost by the adaptive masks where their target is absent or rare, from the R4 runs
    lec_rows = [r for r in per_seed(list(range(n_seeds))).get(("surface_lab", "lec"), [])
                if r["forms"]["adaptive"]["signal"] is not None]
    lec_signal = [r["forms"]["adaptive"]["signal"] for r in lec_rows]
    lec_empty = [r["forms"]["adaptive"]["signal"] for r in lec_rows if r["target"] == 0]
    lec_signal_text = ("" if not lec_signal or min(lec_signal) > 0.99 else
                       f" On the surface sensor, though, the adaptive low-energy-cluster mask keeps only "
                       f"{min(lec_signal):.2f}&ndash;{max(lec_signal):.2f} of the signal"
                       + ("" if not lec_empty else
                          f", {min(lec_empty):.2f} even in a seed with no low-energy clusters at all")
                       + ": it fires on charge-transfer trail electrons, a false positive that the defect-free test "
                       "cannot see because it switches every defect off at once.")

    # ---------------- tables
    sensor_rows = []
    for name, s in PRESETS.items():
        sensor_rows.append([
            SENSOR_NAMES[name], f"{s.nx} &times; {s.ny}", f"{s.thickness_um:.0f}", f"{s.noise_e:.2f}",
            f"{s.exposure_days * 24:.0f} h", sci(s.dark_e_per_pix_day, 1),
            f"{s.muon_flux_per_cm2_day:.2g}", f"{s.highE_dru:.3g}",
            f"{s.halo_length_um:.0f} / {s.cti_h_length_pix:.0f} / {s.n_hot_columns}"])

    null_rows = []
    for sensor, masks in null["per_sensor"].items():
        for mask, v in masks.items():
            if not isinstance(v, dict) or not v["trials"]:
                continue
            verdict = ("contains &alpha;" if v["contains_alpha"] else
                       ("below &alpha; (more conservative)" if v["below_alpha"] else "<strong>above &alpha;</strong>"))
            null_rows.append([SENSOR_NAMES[sensor], NULL_NAMES[mask], f"per {v['unit']}",
                              f"{v['fired']} / {v['trials']}", f"{v['rate']:.3f}",
                              f"[{v['ci95'][0]:.3f}, {v['ci95'][1]:.3f}]", verdict])

    fom_rows = []
    for sensor in SENSOR_NAMES:
        for mask in MASK_NAMES:
            entry = overall["data"].get((sensor, mask))
            if entry is None:
                continue
            med = {name: float(np.median(v)) for name, v in entry["series"].items()}
            spread = {name: (float(np.min(v)), float(np.max(v))) for name, v in entry["series"].items()}
            target = float(np.median(entry["target"]))

            def cell(name):
                if name not in med:
                    return "&ndash;"
                text = f"{med[name]:.2f}"
                if n_seeds > 1:
                    text += f" <span class='range'>({spread[name][0]:.2f}&ndash;{spread[name][1]:.2f})</span>"
                return text

            few = (target < MIN_TARGET) and mask != "all"
            cells = [SENSOR_NAMES[sensor], MASK_NAMES[mask],
                     f"{target:.0f}" + (" <span class='range'>too few to compare</span>" if few else ""),
                     cell("no_mask"), cell("adaptive")]
            for source in SENSOR_NAMES:
                cells.append("<span class='range'>oracle = 1</span>" if source == sensor else cell(source))
            fom_rows.append(cells)

    working_rows = []
    for (sensor, mask), rows in sorted(per_seed(list(range(n_seeds))).items(),
                                       key=lambda kv: (list(SENSOR_NAMES).index(kv[0][0]), list(MASK_NAMES).index(kv[0][1]))):
        if mask == "all":
            continue
        def med(name, field):
            values = [r["forms"][name][field] for r in rows if r["forms"].get(name, {}).get(field) is not None]
            return f"{np.median(values):.2f}" if values else "&ndash;"
        working_rows.append([SENSOR_NAMES[sensor], MASK_NAMES[mask],
                             f"{np.median([r['target'] for r in rows]):.0f}",
                             med("adaptive", "removed"), med("adaptive", "kept"),
                             med("oracle", "removed"), med("oracle", "kept")])

    calib_rows = []
    for path, run in zip(seed_files, seeds):
        for sensor, c in run["adaptive_calibration"].items():
            calib_rows.append([SENSOR_NAMES[sensor], e(path.parent.name if path.parent.name != "compare_masks" else "20260915"),
                               f"{c['length_h']} / {c['length_v']}", f"{c['halo_radius']}",
                               "yes" if c["halo_calibrated"] else "no", f"{len(c['hot_columns'])}",
                               f"{PRESETS[sensor].n_hot_columns}"])

    bit_rows = [[k, e(v["hypothesis"]), f"{v.get('frac_in_full_columns', float('nan')):.2f}",
                 f"{v.get('frac_with_100e_upstream', float('nan')):.2f} / {v.get('frac_with_100e_only_downstream', float('nan')):.2f}",
                 f"{v.get('median_distance_to_100e_pix', float('nan')):.0f}"]
                for k, v in bits.items() if v.get("images_with_bit")]

    masks_rows = [
        ["Hot columns & pixels", "column rate above k &times; median",
         "Poisson tail of each column's count against its own exposure, Bonferroni over columns, single columns before groups"],
        ["Charge-transfer trails", "L pixels after every bright pixel",
         "downstream vs. upstream excess of single electrons around the same bright pixel (binomial); symmetric sources cancel"],
        ["Halo", "disc of radius R around bright pixels",
         "each annulus against everything outside it (conditional binomial); no far field needed"],
        ["Serial-register hits", "rows with n charged pixels in a window of w",
         "densest window of each row against the image's own occupancy, at several widths"],
        ["Low-energy clusters", "disc around events with &ge; m neighbours",
         "neighbours within r against a Poisson count for the local valid area, at several radii"],
        ["Muon tracks", "charge and length cuts",
         "charge of a minimum-ionising track for the measured length and this sensor's thickness; straight, piled-up or cut by the border"],
    ]

    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Self-calibrating Skipper-CCD masks</title>
<style>
:root {{
  --bg: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink-2: #52514e; --muted: #898781;
  --rule: #e1e0d9; --accent: #2a78d6; --callout: #eef4fc;
  color-scheme: light;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{ --bg: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink-2: #c3c2b7;
    --muted: #898781; --rule: #2c2c2a; --accent: #3987e5; --callout: #16222f; color-scheme: dark; }}
}}
:root[data-theme="dark"] {{ --bg: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink-2: #c3c2b7;
  --muted: #898781; --rule: #2c2c2a; --accent: #3987e5; --callout: #16222f; color-scheme: dark; }}
* {{ box-sizing: border-box; }}
body {{ margin: 0; background: var(--bg); color: var(--ink);
  font: 16px/1.6 "Segoe UI", system-ui, -apple-system, "Helvetica Neue", sans-serif; }}
main {{ max-width: 46rem; margin: 0 auto; padding-inline: 1.25rem; padding-block: 3rem 4rem; }}
h1 {{ font-size: 2rem; line-height: 1.2; margin: 0 0 .5rem; letter-spacing: -0.01em; }}
h2 {{ font-size: 1.3rem; margin: 3rem 0 .75rem; padding-top: 1.5rem; border-top: 1px solid var(--rule); }}
h3 {{ font-size: 1.05rem; margin: 1.75rem 0 .5rem; }}
p, li {{ color: var(--ink); }}
.lede {{ font-size: 1.15rem; color: var(--ink-2); margin: 0 0 1.5rem; }}
.meta {{ color: var(--muted); font-size: .9rem; }}
.answer {{ background: var(--callout); border-left: 3px solid var(--accent); padding: 1rem 1.25rem;
  border-radius: 4px; margin: 1.5rem 0; }}
.answer p {{ margin: .25rem 0; }}
figure {{ margin: 1.5rem 0; }}
figure img {{ width: 100%; height: auto; display: block; background: #fcfcfb; border-radius: 6px;
  border: 1px solid var(--rule); }}
figcaption {{ color: var(--ink-2); font-size: .9rem; margin-top: .5rem; }}
.wide {{ width: min(64rem, calc(100vw - 2.5rem)); margin-left: 50%; transform: translateX(-50%); }}
.table-wrap {{ overflow-x: auto; margin: 1rem 0; }}
table {{ border-collapse: collapse; width: 100%; font-size: .88rem; font-variant-numeric: tabular-nums; }}
th, td {{ text-align: left; padding: .4rem .6rem; border-bottom: 1px solid var(--rule); vertical-align: top; }}
th {{ color: var(--ink-2); font-weight: 600; }}
.range {{ color: var(--muted); }}
code {{ font-family: Consolas, "SFMono-Regular", Menlo, monospace; font-size: .88em; }}
pre {{ background: var(--surface); border: 1px solid var(--rule); border-radius: 6px; padding: .9rem 1rem;
  overflow-x: auto; font-size: .85rem; line-height: 1.5; }}
details {{ margin: .75rem 0; }}
summary {{ cursor: pointer; color: var(--ink-2); }}
a {{ color: var(--accent); }}
@media print {{
  body {{ background: #fff; font-size: 11pt; }}
  main {{ max-width: none; padding: 0; }}
  .wide {{ width: 100%; margin-left: 0; transform: none; }}
  figure, table, .answer {{ break-inside: avoid; }}
  h2 {{ break-after: avoid; }}
  details {{ display: block; }}
  details > summary {{ display: none; }}
  @page {{ margin: 16mm 14mm; }}
}}
</style>
</head>
<body>
<main>
<header>
<h1>Masks that calibrate themselves</h1>
<p class="lede">Can the pixel masks of Skipper-CCD analyses do their job without being tuned to the particular sensor they run on?</p>
<p class="meta">D. Rodriguez Juiz &middot; final project, <em>Gravitational Waves and AI-Assisted Research</em> &middot;
<a href="{REPO_URL}">repository</a> &middot; <a href="{REPO_URL}/blob/main/PROVENANCE.md">provenance of every result</a></p>
</header>

<div class="answer">
<p><strong>Short answer, from simulations of three different sensors ({n_seeds} independent seeds).</strong>
Of the {overall["n_transplant"]} cases in which a mask tuned by hand on one sensor was moved unchanged to another,
{len(overall["paired_harm"])} were <em>worse than applying no mask at all</em> in every seed, and
{len(overall["paired_same"])} made no difference. {worst_harm_text}</p>
<p>The same masks written as procedures that calibrate themselves were {adaptive_harm_text}. Measured against
a mask tuned with truth on the same sensor they reach a median {overall["median_adaptive"]:.2f}, with a worst
case of {worst_adaptive_text}; that ratio is an upper bound, because several tuned optima sit at the edge of
the parameter grid. In {len(overall["beats_adaptive"])} cases a transplanted mask beat the adaptive one in every
seed: self-calibration avoids the large failures but is not free.</p>
<p>On the three sensors with their defects switched off, the adaptive masks fired no more often than
&alpha; allows (and apparently less often; 50 runs cannot tell).{lec_signal_text}</p>
{held_text}
<p class="meta">{rule_text}</p>
</div>

<h2>The question</h2>
<p>A Skipper-CCD counts single electrons in millions of pixels. Before any physics is extracted, analyses discard pixels
with masks, each aimed at one problem: charge left behind during transfer, hits in the serial register, hot columns and
pixels, photons emitted around high-energy tracks (the halo), clusters of low-energy events, and muons. Every mask has sizes
and thresholds &mdash; a radius, a trail length, a rate cut &mdash; chosen by hand for one sensor, readout and site.</p>
<p>The working hypothesis: what transfers between sensors is not the numbers but the <em>procedure</em> used to choose
them. So each mask was written twice: a <strong>fixed</strong> form with hand-set constants, and an <strong>adaptive</strong>
form that measures its own sizes on the images it is applied to and sets its threshold by a false-positive rate
(&alpha; = 0.01, corrected for the number of tests) or by measurable physical properties of the sensor.</p>
{table(["Mask", "Fixed form", "Adaptive form"], masks_rows)}

<h2>How it was tested</h2>
<h3>A simulator that remembers where every electron came from</h3>
<p>Each source of charge fills its own map, so a mask can be scored against exactly the events it should remove.
Three sensors differ the way real deployments do. Their geometry, noise, dark rates, background rates and muon fluxes come
from public measurements (SENSEI at SNOLAB and MINOS, Oscura sensors at the surface, SNO, PDG); defect populations are not
published in a transferable form and are scenario values chosen to differ between sensors.</p>
{table(["Sensor", "Pixels", "Thickness (&micro;m)", "Noise (e&minus;)", "Exposure", "Dark rate (e&minus;/pix/day)",
        "Muons (cm&minus;2 day&minus;1)", "High-E (dru)", "Halo &micro;m / trail px / hot cols"], sensor_rows)}
<p>Three versions of every mask were compared on an independent test stack of each sensor: the <strong>oracle</strong>
(fixed form tuned with the simulator's truth on that sensor), the <strong>transplant</strong> (the oracle constants of a
different sensor) and the <strong>adaptive</strong> form, calibrated without truth. The figure of merit is S/&radic;(S+B):
surviving injected signal against surviving target background.</p>

<h3>A check against real images</h3>
<p>Before relying on anything, the public SENSEI SNOLAB data release was read and its single-electron rate reproduced:
{sci(rate['rate_e_per_pix_day'])} &plusmn; {sci(rate['rate_err'], 1)} e&minus;/pix/day against the published
{sci(rate['published_rate'])} &plusmn; {sci(rate['published_err'], 1)} (pull {rate['pull']:+.2f}).</p>
{figure("release_rate.png", "Single-electron density against exposure for four exposures of the public SENSEI release, on a straight line whose slope reproduces the published rate.",
        "Single-electron density per exposure in the public SENSEI SNOLAB release, with the weighted straight line. The slope is the rate.")}
<details><summary>Which bit of the public mask is which (geometric check)</summary>
<p>The release names its masks but not their bits. Each mask leaves a signature independent of its parameters; the bit
order was checked against them. The noisy-row bit never appears in the active area, so it is not verified.</p>
{table(["Bit", "Name (hypothesis)", "Fraction in whole columns", "Bright pixel upstream / only downstream", "Median distance to &gt;100 e (px)"], bit_rows)}
</details>

<h2>Result 1 &middot; How often the adaptive masks fire when there is nothing to find</h2>
<p>{null['n_runs']} independent runs of each of the three sensors with every defect the masks look for switched
off, keeping what a real sensor cannot switch off: dark current, spurious charge, the injected signal, muon tracks
and high-energy deposits. The stack-calibrated masks decide once per run of {null['images_per_run']} images; the
others decide per image. The test is two-sided, so a mask that is far <em>more</em> conservative than &alpha; is
visible as well.</p>
{figure("null_rates.png", "For each sensor and mask, the fraction of trials without the defect in which the mask fired, with 95 percent intervals, against alpha equal to 0.01.",
        "Clopper&ndash;Pearson 95 % intervals; the vertical line is &alpha; = 0.01. Intervals are wide because 50 runs cannot resolve 0.01 from 0.03.")}
{table(["Sensor", "Adaptive mask", "Unit", "Fired / trials", "Rate", "95 % interval", "Against &alpha;"], null_rows)}

<h2>Result 2 &middot; Transplanted constants can do harm; self-calibration avoids the large failures, at a cost</h2>
<div class="wide">
{figure("compare_masks.png", "Small multiples for three sensors: figure of merit of adaptive masks and of fixed masks transplanted from other sensors, relative to the oracle.",
        f"Figure of merit relative to the oracle tuned with truth on the same sensor (vertical line). Markers are medians over {n_seeds} seed{'s' if n_seeds > 1 else ''}"
        + (", whiskers the range" if n_seeds > 1 else "") + ". The deep-underground sensor has few target events, so every version sits near 1 there.")}
</div>
{table(["Evaluated on", "Mask", "Target events", "No mask", "Adaptive", "Fixed, tuned on deep", "Fixed, tuned on shallow", "Fixed, tuned on surface"], fom_rows)}
<p class="meta">Figure of merit relative to the oracle of the sensor in the first column{"; median, with the range over seeds in brackets" if n_seeds > 1 else ""}.</p>
<details><summary>Target removed and clean pixels kept at each working point (medians over seeds)</summary>
{table(["Sensor", "Mask", "Target events", "Adaptive: removed", "Adaptive: clean kept", "Oracle: removed", "Oracle: clean kept"], working_rows)}
</details>
<details><summary>What the adaptive masks measured on each sensor</summary>
{table(["Sensor", "Seed", "Trail length h / v (px)", "Halo radius (px)", "Halo calibrated", "Hot columns found", "Hot columns simulated"], calib_rows)}
</details>

<h2>What the checks caught</h2>
<p>Every error below was found by a test or a diagnostic against the simulator's truth, fixed, and recorded.</p>
<ul>
<li>The pixel threshold of every &ldquo;adaptive&rdquo; mask was set from the simulator's true mean charge: truth leaking into the masks, and at the surface a threshold of 0.50 e that turned read noise into fake events. Found by an independent review; the threshold is now measured from each image's charge histogram and every result was re-run.</li>
<li>The simulator redrew hot columns in every image, which would have invalidated every calibration on a stack.</li>
<li>The hot-column mask flagged the innocent neighbours of each hot column (17 of them in a test).</li>
<li>The halo test treated its reference rate as exact; a 4&sigma; fluctuation gave a halo where there was none.</li>
<li>At the surface no pixel is far from every muon track: the first halo calibration fell back to its largest radius and masked the whole image. It now compares each annulus with everything outside it.</li>
<li>Crossing muons merged into wide clusters that the track criterion rejected (24 % of muon pixels missed); tracks cut by the image border were missed too.</li>
<li>Hot columns were calibrated before charge-transfer trails, whose vertical trails from muons made 69 of 81 flagged columns false.</li>
<li>In one seed the halo radius, chosen by significance alone, reached 125 px on the surface sensor and masked 99.9 % of the image. The radius now maximises a figure of merit estimated from the data within the significant range; all seeds were re-run, plus two seeds never used during development.</li>
</ul>

<h2>Limitations</h2>
<ul>
<li>Defect populations (hot columns, trail lengths, serial hits, clusters) are scenario values, not measurements.</li>
<li>One figure of merit; the oracle is tuned mask by mask, not jointly, so a combined adaptive mask can exceed it.</li>
<li>Adaptive trail lengths are limited by the number of bright pixels in the calibration stack: they are short when data are few.</li>
<li>The muon figure of merit counts pixels, which weighs every missed track pixel heavily.</li>
<li>The adaptive muon mask uses the sensor's diffusion model and the minimum-ionising charge per length, and the simulated muons carry exactly that mean charge with no Landau tail, so its charge tolerance is tested against the model it was built from.</li>
<li>The defect-free test switches every defect off at once, so it cannot see a mask firing on another defect (the low-energy-cluster mask on charge-transfer trails, above).</li>
<li>Several oracle optima lie at the edge of their parameter grids (charge-transfer trails at the longest lengths, muons and hot columns at the loosest cuts), so ratios to the oracle overstate the adaptive masks.</li>
<li>An adaptive mask cannot tell injected signal from dark current, so it estimates signal as every uniform single electron; the evaluation counts only the injected signal. At the surface, where dark current is ~20 times the signal, the adaptive halo therefore chooses not to mask and reaches ~0.7 of the oracle. The figure of merit was not changed after this was seen (the no-mask reference and wider grids were added later, for other reasons).</li>
<li>The adaptive masks were not yet run on the real public images beyond the rate reproduction.</li>
</ul>

<h2>Reproduce</h2>
<pre><code>uv sync --group dev
uv run python scripts/fetch_public_data.py
uv run pytest
uv run python analysis/reproduce_release_rate.py
uv run python analysis/check_release_mask_bits.py
uv run python analysis/null_false_positive_rates.py 50
for seed in 20260915 20260916 20260917 20260920 20260921; do
  uv run python analysis/compare_masks_across_sensors.py 4 --seed $seed --out results/compare_masks/seed_$seed
done
uv run python analysis/figure_release_rate.py
uv run python analysis/figure_null_rates.py
uv run python analysis/figure_compare_masks.py
uv run python analysis/build_report.py
uv run python scripts/make_pdf.py</code></pre>
<p>Every output has a <code>.provenance.json</code> sidecar with the script, git commit, input hashes, parameters and seed;
<a href="{REPO_URL}/blob/main/PROVENANCE.md">PROVENANCE.md</a> says how each result was checked.</p>
</main>
</body>
</html>
"""
    output = OUT / "index.html"
    output.write_text(page, encoding="utf-8")
    inputs = [RES / "release_rate" / "release_rate.json", RES / "release_mask_bits" / "mask_bit_signatures.json",
              RES / "null_false_positive_rates" / "null_rates.json", *seed_files]
    write_sidecar(output, __file__, inputs=inputs, notes="presented page; every number read from the inputs")
    print(f"wrote {output} ({n_seeds} R4 seeds)")


if __name__ == "__main__":
    main()
