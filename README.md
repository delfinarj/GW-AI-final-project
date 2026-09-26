# Sensor-independent pixel masks for Skipper-CCD images

Final project for the course *Gravitational Waves and AI-Assisted Research*
([course repository](https://github.com/matiaszaldarriaga/GW-AI-course),
[project brief](https://matiaszaldarriaga.github.io/GW-AI-course/final-project.html)).

## The question

Skipper-CCD analyses discard pixels with masks, each aimed at one problem: deficient readout
(charge-transfer inefficiency, serial-register hits), material defects (hot columns and pixels),
and physics other than the one under study (halo around high-energy deposits, low-energy
clusters, muon tracks).

**Can each of these masks be made efficient at its own problem without depending on the
particular Skipper-CCD it runs on?**

The working answer tested here: a mask stops depending on the sensor when it is a *procedure that
calibrates itself from the images*, with thresholds set by a controlled false-positive rate or by
measurable physical properties, instead of constants tuned by hand. See [`PLAN.md`](PLAN.md).

## What was found

- Moving hand-tuned constants to another sensor is not merely suboptimal: in 5 of 28 cases the
  transplanted mask scored **below using no mask at all**, in every seed.
- The self-calibrating form avoids those failures (1 of 14 cases below no mask) and reaches a median
  **0.98** of a mask tuned with the simulator's truth on that same sensor &mdash; but it is not free:
  in 8 cases a transplanted mask still beat it in every seed.
- The one real failure: on the surface sensor the adaptive low-energy-cluster mask fires on
  charge-transfer trail electrons. With trails as the only defect present it fires in 40 of 40 images
  and masks 12 % of the image.
- On the real public SENSEI release, which nothing here was tuned for, the hot-column procedure flags
  the seven loudest columns of the longest exposure and nothing else, all of them inside the mask the
  collaboration published.

The full argument, with the figures and what is wrong with it, is the report page; the five-page
version is [`report/summary.pdf`](report/summary.pdf).

## How it is tested

- A simulator ([`src/skmask/simulate.py`](src/skmask/simulate.py)) that records the origin of every
  electron, run for three sensors built from public parameters
  ([`src/skmask/presets.py`](src/skmask/presets.py)): deep underground, shallow underground, and a
  surface laboratory.
- Six masks, each in a fixed and an adaptive form ([`src/skmask/masks.py`](src/skmask/masks.py)).
- The public SENSEI SNOLAB data release as a check against real images.

## Reproduce everything

The commands below assume a POSIX shell (Git Bash, WSL, macOS or Linux); only the `for` loop over
seeds is shell-specific. They need [uv](https://docs.astral.sh/uv/) &mdash;
`curl -LsSf https://astral.sh/uv/install.sh | sh`, or `irm https://astral.sh/uv/install.ps1 | iex` on
Windows &mdash; which installs Python 3.11 itself if it is missing. The whole set takes about seven
hours, almost all of it in R4 and R5; each line says what it costs.

```
# keep the committed numbers to compare against: every analysis overwrites results/ in place
git archive HEAD results | tar -x -C ../reference

# environment: uv installs the exact versions in uv.lock, in seconds and without a solver
uv sync --group dev
uv run python scripts/fetch_public_data.py       # downloads and SHA-256-verifies the public data
uv run pytest                                    # simulator, estimation, events and masks
uv run python analysis/reproduce_release_rate.py        # R1
uv run python analysis/check_release_mask_bits.py       # R2
uv run python analysis/adaptive_masks_on_public_data.py     # R6 (the masks on the real release)
uv run python analysis/null_false_positive_rates.py 50  # R3 (about 40 min)
uv run python analysis/muon_mask_null_rate.py 100        # R7 (about 15 min)
uv run python analysis/cross_defect_false_positives.py 20   # R5 (about 2 h; add --resume to continue one that was killed)
# R4: three seeds used while developing and two run once after the code was frozen
# (about 70 minutes per seed, so about six hours for the five)
for seed in 20260915 20260916 20260917 20260920 20260921; do
  uv run python analysis/compare_masks_across_sensors.py 4 --seed $seed       --out results/compare_masks/seed_$seed
done
uv run python analysis/figure_release_rate.py           # figures
uv run python analysis/figure_null_rates.py
uv run python analysis/figure_cross_defect.py
uv run python analysis/figure_compare_masks.py
uv run python analysis/build_report.py                  # report/index.html, numbers read from results/
uv run python analysis/build_summary.py                 # report/summary.html, the five-page version
uv run python scripts/make_pdf.py                       # report/report.pdf (needs Chrome, Chromium or Edge)
uv run python scripts/make_pdf.py report/summary.html report/summary.pdf
uv run python scripts/compare_results.py ../reference/results results   # did the numbers come back?
```

`compare_results.py` exits non-zero if any result moved by more than a relative 1e-9, and reports the
largest difference per file. Re-running R1, R2, R6 and R7 in a clean clone reproduces the committed
files bit for bit; R3, R4 and R5 have been checked the same way and are recorded in
[`PROVENANCE.md`](PROVENANCE.md).

### Which result is which

| Label | Script | Where it appears on the page |
|---|---|---|
| R1 | `analysis/reproduce_release_rate.py` | "A check against real images", in *How it was tested* |
| R2 | `analysis/check_release_mask_bits.py` | the mask-bit table, in *How it was tested* |
| R3 | `analysis/null_false_positive_rates.py` | Result 1, the table |
| R7 | `analysis/muon_mask_null_rate.py` | Result 1, the paragraph after the table |
| R5 | `analysis/cross_defect_false_positives.py` | Result 2 |
| R4 | `analysis/compare_masks_across_sensors.py` | Result 3 |
| R6 | `analysis/adaptive_masks_on_public_data.py` | Result 4 |

Without uv, `environment.yml` builds the same environment with conda (slower, and conda resolves
versions itself, so small numerical differences are possible).

The presented page is [`report/index.html`](report/index.html) and the PDF
[`report/report.pdf`](report/report.pdf).

Every script writes its output under `results/` with a `.provenance.json` sidecar (script, git
commit, inputs with hashes, parameters, seed). [`PROVENANCE.md`](PROVENANCE.md) lists every result
with where it came from and how it was checked, including the errors the checks caught.

## Layout

| Path | What |
|---|---|
| `src/skmask/` | simulator, presets, thresholding and scoring (`events.py`), masks, public-data reader, provenance helper |
| `tests/` | checks against answers derived independently of the code |
| `analysis/` | one script per result, plus the two that build the pages |
| `results/` | outputs and their provenance sidecars (see [`results/README.md`](results/README.md)) |
| `scripts/` | data download, PDF printing, and `compare_results.py` |
| `report/` | the generated page, the two PDFs and their figures |
| `verification/` | logs of the clean-clone reproductions, the evidence behind the claims in `PROVENANCE.md` |
| `index.html` | the front door served by GitHub Pages |
| `PLAN.md` | the plan, and the statements pre-registered before each run |
| `PROVENANCE.md` | every result: where it came from, how it was checked, what went wrong |
| `METHOD.md` | how the work was done with agents, and what they got wrong |
| `AGENTS.md` | instructions for agents working in this repository |

## Data

Only public data and simulations are used. `data/` is not versioned; the public release
([sensei-skipper/DataReleases](https://github.com/sensei-skipper/DataReleases), CC-BY-4.0) is
downloaded by `scripts/fetch_public_data.py`. Private material used for background reading never
enters the repository.
