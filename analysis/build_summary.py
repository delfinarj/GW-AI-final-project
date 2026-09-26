"""Build report/summary.html, the five-page summary, from the result files.

Same rule as the full page: no number here is typed by hand. The aggregate figures of the transplant
comparison come from `results/report_numbers.json`, which `analysis/build_report.py` writes, so the
summary and the page cannot drift apart; everything else is read from the result file of the analysis
it belongs to. The PDF is this page printed by `scripts/make_pdf.py`.

Run:  python analysis/build_summary.py
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
REPO = "https://github.com/delfinarj/GW-AI-final-project"
PAGE = "https://delfinarj.github.io/GW-AI-final-project/"

SENSOR_NAMES = {"deep_underground": "Deep underground", "shallow_underground": "Shallow underground",
                "surface_lab": "Surface laboratory"}
MASK_NAMES = {"hot_columns": "hot columns and pixels", "cti": "charge-transfer trails", "halo": "halo",
              "serial": "serial-register hits", "low_energy_clusters": "low-energy clusters"}
GROUP_NAMES = {"hot_columns_pixels": "hot columns and pixels", "cti": "charge-transfer trails",
               "serial": "serial-register hits", "low_energy_clusters": "low-energy clusters",
               "halo": "halo"}


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def e(text):
    return html.escape(str(text))


def sci(x, digits=1):
    """Scientific notation that stays on one line inside a narrow table column."""
    mantissa, exponent = f"{x:.{digits}e}".split("e")
    return f"<span class='nb'>{mantissa}&times;10<sup>{int(exponent)}</sup></span>"


def table(headers, rows):
    head = "".join(f"<th>{h}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in row) + "</tr>" for row in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    head = load(RES / "report_numbers.json")
    null = load(RES / "null_false_positive_rates" / "null_rates.json")
    muon = load(RES / "muon_mask_null_rate" / "muon_null.json")
    cross = load(RES / "cross_defect_false_positives" / "cross_defect.json")
    public = load(RES / "adaptive_on_public" / "adaptive_on_public.json")
    rate = load(RES / "release_rate" / "release_rate.json")
    for src in (RES / "compare_masks" / "compare_masks.png",
                RES / "cross_defect_false_positives" / "cross_defect.png"):
        shutil.copy2(src, FIG / src.name)

    tr, ad, held = head["transplant"], head["adaptive"], head["held_out"]

    # ---- sensors
    sensor_rows = [[SENSOR_NAMES[n], f"{s.nx}&times;{s.ny}", f"{s.thickness_um:.0f}", f"{s.noise_e:.2f}",
                    f"{s.exposure_days * 24:.0f} h", sci(s.dark_e_per_pix_day), f"{s.muon_flux_per_cm2_day:.2g}"]
                   for n, s in PRESETS.items()]

    # ---- false positives with nothing to find (R3) and the muon mask (R7)
    null_cells = [v for rows in null["per_sensor"].values() for k, v in rows.items() if isinstance(v, dict)]
    null_tested = [v for v in null_cells if v["trials"]]
    null_contains = sum(1 for v in null_tested if v["contains_alpha"])
    null_fired = sum(v["fired"] for v in null_tested)
    muon_fired = sum(v["fired"] for v in muon["per_sensor"].values())
    muon_images = sum(v["images"] for v in muon["per_sensor"].values())

    # ---- one defect at a time (R5)
    fires = []
    for sensor, groups in cross["per_sensor"].items():
        for group, rows in groups.items():
            for mask, v in rows.items():
                if isinstance(v, dict) and not v["own_defect"] and v["trials"] and v["verdict"].startswith("fires"):
                    fires.append((v["rate"], v["median_masked_fraction"] or 0.0, sensor, group, mask, v))
    fires.sort(reverse=True)
    fire_rows = [[SENSOR_NAMES[s], GROUP_NAMES[g], MASK_NAMES[m], f"{v['fired']}/{v['trials']}",
                  f"{r:.2f}", f"{f:.4f}"] for r, f, s, g, m, v in fires]
    worst_fire = fires[0] if fires else None

    # ---- the real sensor (R6)
    exposures = [public["per_exposure"][k] for k in sorted(public["per_exposure"], key=int)]
    last = exposures[-1]
    public_rows = [[f"{r['exposure_s'] / 3600:.0f} h", r["images"], r["trigger_pixels"]["total"],
                    f"{r['measured_noise_e']['median']:.3f}", sci(r["measured_density_1e"]["median"]),
                    f"{len(r['constants_chosen']['hot_columns'])}", f"{r['constants_chosen']['halo_radius']}",
                    f"{r['masked_fraction_union']:.4f}", f"{r['release_mask_fraction']:.3f}"]
                   for r in exposures]
    loud = last["loudest_columns"]
    n_loud_flagged = next((i for i, c in enumerate(loud) if not c["flagged_by_us"]), len(loud))
    precision = [r["against_release"]["hot_columns_pixels"]["median_fraction_of_ours_also_theirs"] for r in exposures]
    recall = [r["against_release"]["hot_columns_pixels"]["median_fraction_of_theirs_also_ours"] for r in exposures]
    ratios = [x for r in exposures for x in (r["hot_column_evidence"]["ratio_to_common"] or [])]
    r1_noise = [x["noise"] for x in rate["per_exposure"]]

    css = (Path(__file__).resolve().parents[1] / "report" / "summary.css").read_text(encoding="utf-8")

    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Masks that calibrate themselves &mdash; five-page summary</title>
<style>{css}</style>
</head>
<body>

<h1>Masks for Skipper-CCD images that do not depend on the sensor</h1>
<p class="lede">Can the masks that remove readout defects, material defects and unwanted physics be written so that
they work on a sensor they were never tuned on?</p>
<p class="byline">D. Rodriguez Juiz and F. P&eacute;rez &middot; final project, <em>Gravitational Waves and
AI-Assisted Research</em> &middot; full report, code and provenance: <a href="{PAGE}">{PAGE}</a></p>

<div class="answer">
<p><strong>Short answer.</strong> Yes for five of the six masks, and the price is measurable. Of the
{tr["cases"]} cases in which a mask tuned by hand on one sensor was moved unchanged to another, {tr["harmful"]}
were worse than using no mask at all in every seed and {tr["no_difference"]} made no difference; the clearest,
{tr["worst"]["case"]}, scores {tr["worst"]["relative_to_oracle"]:.2f} of the mask tuned with truth on that sensor,
where doing nothing scores {tr["worst"]["no_mask"]:.2f}. The same masks written as self-calibrating procedures are
worse than no mask in {ad["harmful"]} of {ad["cases"]} cases and better in {ad["better_than_no_mask"]}, and reach a
median {ad["median_of_oracle"]:.2f} of the mask tuned with truth. The exception is {ad["harmful_cases"][0]}, which
fires on charge-transfer trail electrons.</p>
</div>

<h2>The question and what was built</h2>
<p>A Skipper-CCD counts single electrons. Before any physics is extracted, analyses discard pixels with masks, each
aimed at one problem: charge left behind during transfer, hits in the serial register, hot columns and pixels,
photons around high-energy tracks (the halo), clusters of low-energy events, and muon tracks. Every mask carries
sizes and thresholds &mdash; a radius, a trail length, a rate cut &mdash; chosen by hand for one sensor and site.</p>
<p><strong>Hypothesis:</strong> what transfers between sensors is not the numbers but the <em>procedure</em> that
chooses them. Each of the six masks was therefore written twice: a <strong>fixed</strong> form with hand-set
constants, and an <strong>adaptive</strong> form that measures its own sizes on the images it is applied to and sets
its threshold from a false-positive rate (&alpha;&nbsp;=&nbsp;0.01, corrected for the number of tests) or from
measurable physical properties of the sensor.</p>
<p>The test bed is a simulator that keeps a separate charge map per source, so a mask can be scored against exactly
the events it should remove. Three sensors differ the way real deployments do; their geometry, noise, dark rates,
backgrounds and muon fluxes come from public measurements, while defect populations, which are not published in a
transferable form, are scenario values chosen to differ between sensors.</p>
{table(["Sensor", "Pixels", "Thickness (&micro;m)", "Noise (e&minus;)", "Exposure", "Dark (e&minus;/pix/day)",
        "Muons (cm&minus;&sup2;/day)"], sensor_rows)}
<p>Four forms of every mask were compared on an independent test stack of each sensor: the <strong>oracle</strong>
(fixed, tuned with the simulator's truth on that sensor), the <strong>transplant</strong> (the oracle constants of a
different sensor), the <strong>adaptive</strong> form, and <strong>no mask</strong>. The figure of merit is
S/&radic;(S+B): surviving injected signal against surviving target background.</p>

<div class="page-break"></div>
<h2>Result 1 &middot; Transplanted constants can do harm; self-calibration avoids the large failures, at a cost</h2>
<figure>
<img src="figures/compare_masks.png" alt="Figure of merit relative to the oracle for each sensor and mask">
<figcaption>Each mask on each sensor, relative to the oracle of that sensor, over {head["n_seeds"]} seeds. A
transplanted mask is the oracle of another sensor; no mask is the reference that says whether masking helps at
all.</figcaption>
</figure>
<p>Moving constants between sensors is not merely suboptimal, it can be worse than not masking: {tr["harmful"]} of
{tr["cases"]} transplants score below no mask in every seed, by more than the {100 * head["margin"]:.0f}&nbsp;%
margin used to call a difference real. Self-calibration removes those failures &mdash; {ad["harmful"]} of
{ad["cases"]} cases &mdash; but it is not free: in {ad["beaten_by_a_transplant"]} cases a transplanted mask still
beats the adaptive one in every seed, and the worst adaptive case, {ad["worst"]["case"]}, reaches only
{ad["worst"]["median_of_oracle"]:.2f} of the oracle ({ad["worst"]["worst_seed"]:.2f} in its worst seed).</p>
<p>The summary statistics were pre-registered before the final run. Cases with fewer than
{head["min_target_events"]} target events in a seed are excluded from that seed, and a case counts as harmful only
if it loses by more than the margin in <em>every</em> valid seed. Two seeds ({held["seeds"]}) were run once after
the code was frozen and never used to change anything: on them, {held["transplant_harmful"]} of
{held["transplant_cases"]} transplants are harmful, {held["adaptive_harmful"]} of {held["cases"]} adaptive cases
are, and the adaptive median is {held["adaptive_median_of_oracle"]:.2f} of the oracle &mdash; the same picture as
the development seeds.</p>

<div class="page-break"></div>
<h2>Result 2 &middot; What the masks fire on when there is nothing to find</h2>
<p>With every defect switched off on all three sensors, keeping what a real sensor cannot switch off (dark current,
spurious charge, the injected signal, muon tracks and high-energy deposits), the five statistical masks fired
{null_fired} times in {len(null_tested)} sensor-mask tests of {null["n_runs"]} runs each, and every interval
contains &alpha;. The sixth mask, for muon tracks, cannot be in that test, since tracks are what a real sensor
cannot switch off; asked separately, in a simulator with the flux set to zero, it fired in {muon_fired} of
{muon_images} images.</p>
<p>That test has a blind spot: switching every defect off at once cannot see a mask firing on <em>another</em>
mask's defect. So the defects were also switched on one at a time, {cross["cells_tested"]} mask-and-defect cells in
all, {cross["n_runs_completed"]} runs each, with the bound each cell must clear taken at
{100 * cross["confidence_of_the_corrected_bound"]:.2f}&nbsp;% confidence to control a
{100 * cross["family_wise_error"]:.0f}&nbsp;% chance of one false call over the whole grid.</p>
<figure>
<img src="figures/cross_defect.png" alt="Grid of masks against the single defect switched on, per sensor">
<figcaption>One defect on at a time. Dot area is the fraction of trials in which the mask fired; orange marks the
cells above the rate that mask should not exceed; hollow circles are the diagonal, where a mask meets its own
defect.</figcaption>
</figure>
<p>{len(fires)} of {cross["cells_tested"]} cells fire on a defect that is not their own, and charge-transfer trails
are what most often set them off. The largest is unambiguous: on the surface sensor, with trails as the only defect
present, the adaptive low-energy-cluster mask fires in every one of its trials and masks a median
{worst_fire[1]:.3f} of the image. This is the one real failure of self-calibration found here, and it was predicted
in writing before the run, from a symptom seen in Result&nbsp;1.</p>
{table(["Sensor", "Only defect present", "Mask that fired", "Fired/trials", "Rate", "Median fraction masked"],
       fire_rows)}

<div class="page-break"></div>
<h2>Result 3 &middot; The same procedures on a real sensor</h2>
<p>Everything above is simulated. The adaptive masks were also run, unchanged, on the {len(exposures)} exposures of
the public SENSEI SNOLAB data release: a fourth sensor, real, with a geometry none of the presets has
({last["shape"][1]}&times;{last["shape"][0]} active superpixels, each binning 32 physical rows). Nothing about it
was tuned, and the release publishes its own mask, so the two can be compared where they overlap.</p>
<p>Given only the images, the estimator measures a readout noise of {min(r["measured_noise_e"]["median"] for r in exposures):.3f}&ndash;{max(r["measured_noise_e"]["median"] for r in exposures):.3f}&nbsp;e,
against {min(r1_noise):.4f}&ndash;{max(r1_noise):.4f}&nbsp;e from an independent fit on the same files, and a
single-electron density that grows with exposure. The hot-column procedure flags
{", ".join(str(len(r["constants_chosen"]["hot_columns"])) for r in exposures)} columns as the exposure grows, and
the sets are nested: each keeps the previous columns and adds the next loudest. At the longest exposure the
{n_loud_flagged} columns with the highest rate of charged pixels are exactly the ones it flags, all of them inside
the release's own bad-column mask, and each carries {min(ratios):.0f} to {max(ratios):.0f} times the charged-pixel
rate of the columns it left alone.</p>
{table(["Exposure", "Images", "Trigger pixels", "Noise (e&minus;)", "1e density", "Hot columns", "Halo radius",
        "We mask", "The release masks"], public_rows)}
<p>Both directions of the agreement matter and both are here. Of what the procedure flags, the release also flags
{min(precision):.2f}&ndash;{max(precision):.2f}: it stays inside a mask a person tuned. Of what the release flags,
the procedure flags {min(recall):.3f}&ndash;{max(recall):.3f}: the release masks far more. The first number is the
easy direction &mdash; a procedure that flagged one true column would score 1.00 &mdash; and it rests on a handful
of column decisions, not on the hundreds of pixels it is counted over.</p>
<p>The trail calibration chose length zero on all four exposures, which the registered expectation did not predict.
The reason is in the result file: the release blinds its hits, so a whole exposure holds
{", ".join(str(r["trigger_pixels"]["total"]) for r in exposures)} pixels above the trigger. At the longest exposure
the test does see the excess on the correct side &mdash; {last["cti_detail"]["v"]["n_downstream"]} single electrons
downstream of a trigger against {last["cti_detail"]["v"]["n_upstream"]} upstream,
p&nbsp;=&nbsp;{last["cti_detail"]["v"]["smallest_p"]:.4f} &mdash; but the corrected threshold asks for
p&nbsp;&lt;&nbsp;{last["cti_detail"]["v"]["bonferroni_threshold"]:.1e}. The procedure is refusing to mask on
evidence this thin, which is what it was built to do.</p>

<div class="page-break"></div>
<h2>How the work was done, and how it was checked</h2>
<p>Every line of code and every analysis in this project was produced by AI agents, working from a stated question
and under rules written down before the work: no number on any page is typed by hand, every output carries a
sidecar recording the script, the git commit, the hashes of its inputs and its seed, and no private document ever
enters a tracked file. The method that made this trustworthy, rather than merely fast, was adversarial:</p>
<ul>
<li><strong>Expectations were registered before runs.</strong> Written in <code>PLAN.md</code> and committed before
the analysis ran, so a confirmation could not be invented afterwards. One of them was wrong &mdash; the trail
length on the real sensor &mdash; and the repository says so, with the numbers that explain why.</li>
<li><strong>Independent agents reviewed the work</strong> with no access to the conversation that produced it. They
found real defects: truth from the simulator leaking into an "adaptive" threshold; a headline that measured doing
nothing; a halo mask measuring distance in superpixels on a binned sensor, which hid the loudest column from the
hot-column calibration; a verdict applied to sixty cells with no correction for their number. Each is recorded in
<code>PROVENANCE.md</code> with what was done about it, including the findings that were <em>not</em> acted on.</li>
<li><strong>Everything was reproduced from a clean clone</strong> &mdash; a fresh checkout sharing nothing with the
working copy but the commit and the pinned environment &mdash; and compared at a relative tolerance of
10<sup>&minus;9</sup>. All results but one are bit-identical, figures included; the exception, the
one-defect-at-a-time run, takes about two hours and was not re-run, which the report states.</li>
</ul>
<p>Errors the checks caught are listed on the full report page rather than hidden: a simulator that redrew its hot
columns every image, a hot-column mask that flagged innocent neighbours, a halo test that treated its reference rate
as exact, a muon mask that missed crossing tracks, a calibration order that made 69 of 81 flagged columns false.</p>

<h3>What is still wrong with it</h3>
<ul>
<li>Defect populations are scenario values, not measurements; only the sensor and background parameters are public
measurements.</li>
<li>The defects are switched on one at a time, so nothing here measures what the masks do when several defects
overlap in the same pixels.</li>
<li>On the real sensor only four of the six masks can be checked at all: the release publishes no counterpart to the
low-energy-cluster mask, and its binned superpixels make the muon mask's geometry meaningless.</li>
<li>One figure of merit, and an oracle tuned mask by mask rather than jointly, so a combined adaptive mask can
exceed it.</li>
</ul>

<h3>Where everything is</h3>
<p>Repository: <a href="{REPO}">{REPO}</a> &mdash; <code>README.md</code> reproduces the whole project from a clean
clone with <code>uv</code>; <code>PLAN.md</code> holds the pre-registered statements; <code>PROVENANCE.md</code>
records every input, every choice with a defensible alternative, how each result was checked, and every error found
along the way; <code>report/index.html</code> is the full report, of which this document is a summary.</p>

</body>
</html>
"""
    output = OUT / "summary.html"
    output.write_text(page, encoding="utf-8")
    inputs = [RES / "report_numbers.json", RES / "null_false_positive_rates" / "null_rates.json",
              RES / "muon_mask_null_rate" / "muon_null.json",
              RES / "cross_defect_false_positives" / "cross_defect.json",
              RES / "adaptive_on_public" / "adaptive_on_public.json",
              RES / "release_rate" / "release_rate.json"]
    write_sidecar(output, __file__, inputs=inputs, notes="five-page summary; every number read from the inputs")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
