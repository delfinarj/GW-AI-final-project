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

## How it is tested

- A simulator ([`src/skmask/simulate.py`](src/skmask/simulate.py)) that records the origin of every
  electron, run for three sensors built from public parameters
  ([`src/skmask/presets.py`](src/skmask/presets.py)): deep underground, shallow underground, and a
  surface laboratory.
- Six masks, each in a fixed and an adaptive form ([`src/skmask/masks.py`](src/skmask/masks.py)).
- The public SENSEI SNOLAB data release as a check against real images.

## Reproduce everything

```
# environment: uv installs the exact versions in uv.lock, in seconds and without a solver
uv sync --group dev
uv run python scripts/fetch_public_data.py       # downloads and SHA-256-verifies the public data
uv run pytest                                    # simulator, estimation, events and masks
uv run python analysis/reproduce_release_rate.py        # R1
uv run python analysis/check_release_mask_bits.py       # R2
uv run python analysis/null_false_positive_rates.py 50  # R3 (about 40 min)
# R4: three seeds used while developing and two run once after the code was frozen
for seed in 20260915 20260916 20260917 20260920 20260921; do
  uv run python analysis/compare_masks_across_sensors.py 4 --seed $seed       --out results/compare_masks/seed_$seed
done
uv run python analysis/figure_release_rate.py           # figures
uv run python analysis/figure_null_rates.py
uv run python analysis/figure_compare_masks.py
uv run python analysis/build_report.py                  # report/index.html, numbers read from results/
uv run python scripts/make_pdf.py                       # report/report.pdf (needs Chrome, Chromium or Edge)
uv run python scripts/compare_results.py <a reference results/> results  # optional: compare with committed numbers
```

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
| `analysis/` | one script per result |
| `results/` | outputs and their provenance sidecars |
| `scripts/` | data download |
| `AGENTS.md` | instructions for agents working in this repository |

## Data

Only public data and simulations are used. `data/` is not versioned; the public release
([sensei-skipper/DataReleases](https://github.com/sensei-skipper/DataReleases), CC-BY-4.0) is
downloaded by `scripts/fetch_public_data.py`. Private material used for background reading never
enters the repository.
