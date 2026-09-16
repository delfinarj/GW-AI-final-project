# Provenance

One entry per result that appears in the deliverables, and per input the results depend on.
Each entry says where the result came from and how it was checked.

## Inputs

### Public SENSEI SNOLAB data release

- **What:** four ROOT files of binned, calibrated Skipper-CCD images (exposures 0, 2, 6 and 20
  hours) with per-pixel mask bits, and the release macro `plotRate.C`.
- **Source:** https://github.com/sensei-skipper/DataReleases, folder `Snolab_binned_1`,
  licence CC-BY-4.0; the data behind arXiv:2410.18716.
- **Obtained by:** `scripts/fetch_public_data.py` (written from scratch; standard library only).
- **Check:** SHA-256 of each file recorded on first download (2026-09-14) and verified on every
  run; the script exits non-zero on a mismatch.
- **Mask bits:** the public README of the release names the eight masks selected by `0x067d`
  (neighbour, bleeding, high-energy halo, crosstalk, noisy row, edge, bad pixel, bad column) but
  does not say which bit is which. The bit order is a hypothesis tested in R2 below; nothing
  about the bits is taken from any non-public document.
- **What the macro establishes about the files:** tree `calPixTree` with branches `x`, `y`,
  `ePix` (charge in electrons) and `mask`; active rows `0 < y <= 16`; mask selection `0x067d`
  for the single-electron analysis; images of 3200 columns by 20 rows. To be confirmed by
  reading the files directly, and recorded here when done.

### Simulator parameters

Defined in `src/skmask/presets.py`, each marked MEASURED (with its public source) or SCENARIO.

| Quantity | Value(s) | Source |
|---|---|---|
| Pixel size, thicknesses, formats | 15 um; 665 / 675 / 725 um | arXiv:2410.18716, arXiv:2004.11378, arXiv:2304.04401 |
| Readout noise | 0.14 / 0.14 / 0.17 e | same |
| Dark and spurious rates | 1.39e-5 e/pix/day, 6.94e-5 e/superpix; 1.594e-4 e/pix/day, 1.664e-4 e/pix; 0.03 e/pix/day, 8.4e-4 e/pix | same |
| High-energy background | 50 / 3370 / 3e4 events/(kg day keV) | same (surface: within the quoted O(1e4-1e5)) |
| Muon flux | 3.31e-10 cm^-2 s^-1; 0.80 m^-2 s^-1; 1 cm^-2 min^-1 | arXiv:0902.2776; Fermilab-thesis-2014-08; PDG cosmic-ray review 2019 |
| Muon zenith distribution | cos^2 intensity | PDG cosmic-ray review 2019 |
| Diffusion | sigma = sqrt(-A ln\|1-bz\|), A = 218.715 um^2, b = 1.015e-3 /um | arXiv:2004.11378 |
| Minimum ionisation | 1.664 MeV cm^2/g, 2.329 g/cm^3 | PDG atomic and nuclear properties |
| Pair energy | 3.75 eV | arXiv:2004.10709 |
| Halo yield | 1e-4 per deposited electron | SCENARIO, order of magnitude from arXiv:2011.13939 and arXiv:2304.04401 |
| Halo length, CTI, hot columns/pixels, serial hits, low-energy clusters | differ per preset | SCENARIO |

**Check:** the public numbers are re-derived in `tests/test_presets.py` (per-pixel mass of the
MINOS sensor, active mass of the SNOLAB sensor, surface 1e rate of order 1e-2 e/pix/day).

### Modelling choices with a defensible alternative (simulator)

All marked `CHOICE` in `src/skmask/simulate.py`:

1. Halo photons use a single exponential absorption length; photons leaving the silicon are lost.
   Alternative: wavelength-dependent absorption (arXiv:2011.13939).
2. Muon charge per length is Poisson around the mean minimum-ionising value. Alternative: Landau
   fluctuations with delta rays.
3. Muon entry points over the area enlarged by 20 thicknesses; projected lengths capped there.
   Alternative: exact geometry; the capped fraction is ~6e-6.
4. The zenith distribution at depth is kept cos^2, although it is steeper underground.
5. High-energy background is point-like (diffusion-limited). Alternative: extended electron tracks.
6. Deferred CTI charge is added without removing it from the source pixel (bias ~p << 1).
7. The MINOS 1e rate is used as dark rate, although part of it is halo.
8. Hot-column and hot-pixel positions are drawn from a per-sensor `defect_seed`, so they are the
   same in every image of a sensor (an earlier version redrew them per image, which would have
   made every stack-based calibration meaningless; caught by the hot-column injection test).

### Masks (`src/skmask/masks.py`) — written from scratch

Library calls: `scipy.ndimage` (labelling, distance transform, dilation), `scipy.signal.fftconvolve`
(neighbour counts in discs), `scipy.stats.poisson` and `binom` (tail probabilities).

Choices with a defensible alternative:

1. A pixel with >= 20 e (`TRIGGER_E`) triggers halo and CTI masks. Alternatives: 100 e (public
   SENSEI 2020 cuts) or a per-sensor trigger.
2. The single-electron threshold is the Bayes cut c = 1/2 + sigma^2 ln((1-mu)/mu) from the measured
   noise and density (`events.py`). Alternatives: a fixed 0.7 e, or maximising F1.
3. Every adaptive threshold is a Bonferroni-corrected tail probability at alpha = 0.01.
   Alternatives: false-discovery-rate control; a likelihood-ratio scan.
4. CTI length uses the downstream-over-upstream excess (binomial test), so symmetric sources
   (halo, dark current, signal) cancel without being modelled.
5. The halo radius compares each annulus with everything outside it pooled (annuli farther out
   plus any far field), and within the annuli that are significant it maximises a figure of merit
   estimated from the data, so the exposure a larger radius costs is weighed. With no outside pool
   of at least 10 000 valid pixels the procedure reports `calibrated = False` and returns 0, leaving
   the decision to the caller. (Both parts replaced earlier behaviour that failed: see the fixes
   recorded further down.)
6. The muon mask uses only physical properties of the sensor (thickness, pixel size, back-surface
   diffusion) and the minimum-ionising charge per length; tolerance factor 2 on the charge.

Errors found by the tests and fixed (kept here because they are part of how the masks were checked):

- **Hot columns flagged their neighbours.** In the first version every group of 2-3 adjacent
  columns containing one hot column was significant because of that column, so 17 innocent
  neighbours were flagged. Fix: single columns take precedence; wider groups are tested only when
  none of their columns is flagged, restarting from width 1 after every new flag.
- **The halo test treated the far-field rate as exact.** A 4-sigma fluctuation of the uniform
  field (seed 5 of the test configuration, present in the true charge maps) gave radius 45 with no
  halo. Fix: a conditional binomial test of the annulus count against the far-field count, which
  carries the far field's own Poisson uncertainty. The false-positive *rate* is measured over many
  realisations in R3, not by one unit test.
- **An underpowered injection test.** Hot columns adding ~4.5 counts per column over a stack
  cannot pass a Bonferroni cut that needs >= 8; the test now injects a defect it can detect, with
  the power calculation written next to it.
- **A halo null test that passed for the wrong reason** (see above, no far field).

Design errors found by the first cross-sensor run (R4 at one image per stack, 2026-09-14):

- **The halo calibration relied on a far field that does not exist at the surface.** With ~200
  muon tracks per image no pixel is 200 px from all of them; the calibration reported
  `calibrated = False` and fell back to the largest radius, which masked 100 % of the surface image.
  Redesign: each annulus is compared with everything outside it pooled; an uncalibrated result
  now returns radius 0 with the flag set, so the caller must decide, instead of silently masking
  everything. A unit test with triggers 60 px apart (no far field at all) checks the new behaviour.
- **The muon mask missed crossing tracks.** A diagnostic on a surface image found 24 % of muon
  pixels unmasked, 98 % of them in clusters that fail the straightness cut with 1.2-2.9 times a
  single track's charge: merged crossing muons. Fix: a pile-up clause (enough charge for the
  length and at least two vertical tracks' worth), still using only physical sensor properties.
- **Hot columns were calibrated before CTI.** On the surface preset the adaptive hot-column mask
  flagged 81 columns, 69 of them false, and found only 12 of 25 true ones. A diagnostic using the
  simulator's truth showed the low-charge occupancy of the false columns to be 64 per column from
  CTI against 8 from dark current: vertical CTI trails of ~200 muon tracks are a column-wise excess.
  With CTI and halo excluded first, 15 columns are flagged, 13 true and 2 false. Fix: calibration
  order CTI -> halo -> hot columns. There is no circularity: the CTI calibration compares
  downstream with upstream of the same trigger in the same row or column, which a hot column,
  uniform along its column, cannot bias.

Second cross-sensor run at one image per stack, after these fixes (statistics too small for
conclusions; used to find design faults):

- Surface, six masks combined: adaptive FoM 1.2 with 43 % of clean pixels kept; fixed masks tuned
  with truth on the same sensor 1.0 (21 %); transplanted from the other sensors 0.6 (6-12 %). The
  combined "oracle" is tuned mask by mask, not jointly, so an adaptive combination can beat it.
- Shallow, CTI: adaptive FoM 6.1 against 10.5 for the oracle; the adaptive trail lengths are
  limited by the number of triggers (seen before in the unit test), to be re-measured with more
  images.
- Surface, muons: adaptive removes 98.9 % of muon pixels, but the survivors (tracks cut by the
  image border, per the earlier diagnostic) weigh heavily in a pixel-counted FoM (4.8 against
  10.5). Fix: a clause for border-clipped tracks (straight, longer than two diffusion-blob lengths,
  charge at least a minimum-ionising charge for the projected length / tolerance), with a unit test
  that a clipped segment is masked and a clipped 2000 e blob is not.

**Check:** `tests/test_masks.py` — hand-computed geometry for each mask; a **null test** for every
adaptive mask (defect absent, mask does not fire) and an **injection test** (defect present,
most of its events removed, scored against the simulator's truth). An earlier halo null test
passed for the wrong reason (no far field, so no significance was possible); it now also asserts
that the calibration actually ran.

**Check:** `tests/test_simulator.py` compares the simulator with answers derived independently
of the code: the diffusion width at a depth, the minimum-ionising charge per length, the Poisson
mean of uniform rates, the muon count and mean cos(theta) = 4/5 for a cos^2 intensity through a
horizontal plane, the mean projected halo distance (pi/4 of the absorption length), CTI trails
only downstream with the set mean length, and serial hits confined to single rows.

### Figures

All figures share `src/skmask/plotstyle.py`: solid hairline grid and axes, markers with a
surface-coloured ring, text in ink colours only, and two series colours (blue `#2a78d6`, orange
`#eb6834`; dark-surface steps `#3987e5`, `#d95926`) from a reference categorical palette.

**Check:** the palette was run through the colour validator of the data-visualisation guidance
used for this project (`validate_palette.js`, Node 2026-09-14), light and dark:
all checks pass — adjacent colour-vision-deficiency separation Delta E 24.7 (light) and 26.8 (dark)
against a target of 8, normal-vision separation 33.6 and 31.8 against a floor of 15, contrast
>= 3:1 on both surfaces. Series identity never relies on colour alone: the R4 figure also uses
marker shape, and every figure has a legend or a single labelled series.

## Reproduction from a clean clone

### Short pass (2026-09-15)

- **What was done:** `git clone` of commit `ab04ab9` into an empty folder, `uv sync --group dev`,
  `scripts/fetch_public_data.py` (all checksums verified), `pytest` (51 passed), then R1, R2 and R4
  seed 20260920 re-run in the clone and compared with the committed files by
  `scripts/compare_results.py` (relative tolerance 1e-9).
- **Environment actually used by the clone:** Python 3.14.7, numpy 2.5.3, scipy 1.18.1, because that
  commit did not yet pin Python and uv chose the newest interpreter installed. The committed results
  were made with Python 3.11.16, numpy 2.4.6, scipy 1.17.1.
- **Result:** R4 seed 20260920 and R2 are identical to the committed files (largest relative
  difference 0). R1 differs only in the fit: the single-electron rate by 1.7e-5 relative, its
  uncertainty by 7.8e-5, the intercept by 6.6e-5; the largest relative difference, 6.6e-3, is on the
  pull, a number close to zero (-0.0231 against -0.0232). No physical statement changes.
- **Cause of the R1 difference, first attributed wrongly.** This entry first blamed scipy 1.18. A second
  clone with Python pinned to 3.11 (numpy 2.4.6, scipy 1.17.1, the same versions as the committed
  results) gave *exactly the same* R1 as the Python 3.14 clone, so the version was not the cause. The
  committed R1 had been made in the conda environment, whose numpy and scipy are built against a
  generic BLAS/LAPACK; uv installs the PyPI wheels built against OpenBLAS (checked with
  `numpy.show_config()`), and the minimiser of R1 follows a slightly different numerical path. R2 and
  R4 do not depend on that minimiser, which is why they matched. R1 and its figure were re-run under uv,
  the declared environment, and now match the pinned clone bit for bit (largest relative difference 0).
- **What it showed about the repository:** the unpinned Python version reported by the second
  independent review was real (the first clone ran Python 3.14); Python is now pinned to 3.11
  (`.python-version`, `requires-python`) and `uv.lock` resolves a single numpy. And a result made in
  another environment than the declared one does not reproduce to the last digit even at equal package
  versions, so every committed result must come from the uv environment.
- **Not yet done:** the complete pass (R3 and the other four R4 seeds, figures, page, PDF), about
  2.5 hours, which must run with the computer attended.

### Re-run of R4 with the wider grids (2026-09-15, evening)

- **Code frozen** at commit `121c688` (wider grids, declared physical bounds, and a run that aborts
  right after the oracle search if an optimum sits on any other edge). Every one of the five seeds
  records `"oracle_edge_check": "passed"`.
- **How the seeds were run:** two workers in parallel plus a third process for seed 20260921, which
  therefore carries commit `c4b59f2` in its sidecar; between `121c688` and `c4b59f2` only R1's result
  files and this document changed, not the code that produces R4.
- **Seed 20260917 crashed once and was re-run on unchanged code.** The first attempt stopped with
  `IndexError: index 32768 is out of bounds for axis 0 with size 1428` in `low_charge_occupancy`,
  where a label larger than the number of clusters appeared in the label map. A diagnostic that
  re-simulated the same seed and replayed the same grids in the same order (hot, then trails, then
  halo, then the serial register) found: the electron maps are bit-identical before and after every
  grid (SHA-1 per frame), the offending frame has 1427 clusters with `labels.max() == 1427`, and the
  serial-register grid runs without error. The failure therefore does not reproduce on the same data
  and the same code; 32768 is 2^15, a single bit set, which fits a transient memory error on a
  loaded machine better than a defect in the code. Nothing was changed and the seed was re-run.
  **This seed has since been re-run from a clean clone** of the frozen commit, in the pinned
  environment, and the result is bit-identical to the committed one (largest relative difference
  0.0e+00). So whatever the crash was, it left no trace in the numbers: the concern that it might
  have a silent counterpart, one that changed a value instead of raising, is answered for this seed
  by an independent re-run rather than by the diagnostic alone.
- **One commit message is wrong:** `d83ef7c` says "R4 re-run over the five seeds" but carries four;
  seed 20260917 was still re-running. The following commit adds it and says so.
- **Reproduction of R4 from a clean clone** of `121c688` with the pinned environment (Python 3.11.16,
  numpy 2.4.6, scipy 1.17.1), each seed re-run from scratch and compared with
  `scripts/compare_results.py` at a relative tolerance of 1e-9: all five seeds (20260915, 20260916, 20260917, 20260920 and 20260921) were re-run from scratch in the clone and every one is identical to the committed file, largest relative difference 0.0e+00. `relative_fom_summary.json`, which the figure script derives from the five, is identical too, and `compare_masks.png` came out byte for byte the same. Together with R1, R2 and R3, reproduced the same way, every result the deliverables rest on has now been produced twice: once here and once in a clone of the frozen commit that shares nothing with this working copy but the commit and the pinned environment (Python 3.11.16, numpy 2.4.6, scipy 1.17.1). The seeds were run in two passes, three on 2026-09-15 and two on 2026-09-16, because each takes about 70 minutes and the machine cannot be left running unattended.
- **How to read the commit in a sidecar.** A sidecar records the commit that was checked out *while the
  file was being written*, which is by construction the parent of the commit that stores the file: the
  output cannot be committed before it exists. So `report/report.pdf.provenance.json` names `d83ef7c`
  although the file is stored in `377dbea`, and seed 20260921, run in the third process, names
  `c4b59f2`. To
  check what code produced an output, read the commit named in its sidecar, not the commit that
  contains it; `changed_outputs` counts files under `results/` and `report/` that differed from the
  index at that moment, which is expected while a set of results is being regenerated.

### Complete pass at the current code (2026-09-16)

Every analysis except R5 has now been re-run from a clone that shares nothing with this working copy
but the commit and the pinned environment, and compared with `scripts/compare_results.py` at a
relative tolerance of 1e-9.

- **R1, R2, R3 and R6**, in a fresh clone of commit `e9c797b`: `git clone`, `uv sync --group dev`,
  `scripts/fetch_public_data.py` (all four checksums verified), then each analysis re-run. All four
  results identical, largest relative difference 0.0e+00; the figures of R1 and R3 came out byte for
  byte the same.
- **R7**, in a second fresh clone, of commit `3964385` (it needs no public data): identical.
- **R4**, the five seeds, was reproduced this way at commit `121c688` and that reproduction still
  stands: of every file R4 reads, `git diff 121c688 HEAD -- src/skmask/simulate.py src/skmask/masks.py
  src/skmask/events.py src/skmask/estimate.py src/skmask/presets.py
  analysis/compare_masks_across_sensors.py` reports one file changed, 6 insertions and 3 deletions,
  all of them inside a module docstring. (`src/skmask/stats.py` was added in between; R4 does not
  import it.)
- **R5 was not reproduced from a clone.** Twenty runs of fifteen configurations take about two hours
  and the machine cannot be left running unattended. What stands in its place is weaker and is worth
  naming as such: each run is seeded from the run index, so the numbers do not depend on how the run
  was split, and `tests/test_cross_defect_resume.py` checks that; the summary is recomputed from the
  counts the file carries; and the file records the commit of every stretch of running. None of that
  is an independent re-run.

## Results

### R1. Single-electron rate of the public SENSEI SNOLAB release

- **Output:** `results/release_rate/release_rate.json` and its `.provenance.json` sidecar.
- **Produced by:** `analysis/reproduce_release_rate.py`, reading through `src/skmask/sensei_public.py`.
- **Inputs:** the four release ROOT files (SHA-256 in the sidecar and in `scripts/fetch_public_data.py`).
- **Written from scratch:** the histogram, the two-Gaussian binned Poisson likelihood, the
  numerical Hessian for the density uncertainty, the exposure bookkeeping and the weighted
  straight-line fit. **Called from libraries:** `uproot` (reading ROOT), `scipy.optimize.minimize`
  (L-BFGS-B), `scipy.stats.norm`.
- **Choices with an alternative:** the Hessian-based uncertainty instead of MINOS/Minuit2
  (the macro's minimiser); fitting on bin centres rather than integrating each bin.
- **Result:** 1.386e-5 +- 0.11e-5 e-/pix/day.
- **Check:** against the published (1.39 +- 0.11) e-5 e-/pix/day (arXiv:2410.18716, Golden
  quadrant): pull -0.02. The fitted read noise, 0.141-0.142 e- in all four exposures, matches the
  published 0.14 e-. File contents confirmed directly: tree `calPixTree` with branches `x`, `y`,
  `ePix`, `mask`, `ohdu`, `RUNID`, `LTANAME`; 14, 14, 14 and 13 images (0, 2, 6, 20 h), each
  3200 x 20 superpixels, one quadrant (`ohdu` 2).
- **Second check:** the fitted exposure-independent density, about 7.2e-5 e-/superpix, agrees
  with the published (6.94 +- 0.85) e-5 e-/superpix/image.
- **Figure:** `results/release_rate/release_rate.png`, drawn by `analysis/figure_release_rate.py`
  from the JSON above (no recomputation); inspected by eye: four densities on the fitted line, the
  printed rate identical to the JSON.

### R2. Which bit of the public release's mask is which mask

- **Output:** `results/release_mask_bits/mask_bit_signatures.json` and sidecar.
- **Produced by:** `analysis/check_release_mask_bits.py` on the 6-hour file (the one with the most
  bits present).
- **Hypothesis:** the eight masks named in the public README appear in bit order: 0x1 neighbour,
  0x4 bleeding, 0x8 halo, 0x10 crosstalk, 0x20 noisy row, 0x40 edge, 0x200 bad pixel,
  0x400 bad column.
- **Written from scratch:** the per-image reshaping and the geometric signatures; `scipy.ndimage`
  for the distance transform and the 3x3 neighbourhood.
- **Choices with an alternative:** a column or row counts as fully masked at >= 90 % of its
  length; "bright" means > 100 e-; physical distance uses 32 rows per superpixel.
- **Check (the result is the check):** every bit with a distinctive geometry shows the geometry
  its name predicts. Neighbour: 96 % of masked superpixels touch charge > 0.5 e-. Bleeding: 100 %
  have a > 100 e- superpixel upstream and 0 % only downstream (one-sided, as charge-transfer
  trails are). Halo: masked superpixels at a median 51 px from > 100 e- charge, with upstream and
  downstream-only fractions 0.20 and 0.18 (symmetric). Edge: whole columns at a median 12 px from
  the border. Bad column: whole columns in the interior. Crosstalk and bad pixel: isolated, with
  no > 100 e- charge in the same quadrant, as expected for crosstalk from other quadrants; not
  distinguishable from each other by geometry. **Noisy row (0x20): never set inside the active
  area of this file, so its bit is not verified.**

### R3. How often the adaptive masks fire when their target is absent

- **Output:** `results/null_false_positive_rates/null_rates.json`, sidecar, and figure
  `null_rates.png` (`analysis/figure_null_rates.py`).
- **Produced by:** `analysis/null_false_positive_rates.py 50`: 50 independent runs of 2 images on each
  of the three presets with every target defect switched off (hot columns and pixels, charge-transfer
  trails, serial-register hits, low-energy clusters, halo), keeping dark current, spurious charge, the
  injected signal, muons and high-energy deposits. Seed 20260915. Code frozen at commit `1a0c649`.
- **Written from scratch:** the run loop, the two trial units and the Clopper-Pearson interval
  (`scipy.stats.beta`).
- **Choices with an alternative:** "fires" means any flag at all in the trial (the strictest reading);
  stack-calibrated masks (hot columns and pixels, trails, halo) are counted per run, per-image masks
  (serial register, low-energy clusters) per image; the hot-column trial fires on a column *or* a
  pixel, each tested at alpha, so its nominal per-run rate is about 2 alpha.
- **Result:** every one of the 15 sensor-mask intervals contains alpha = 0.01, with 0 or 1 firings each
  (numbers in the JSON and on the page). The deep-underground halo could be calibrated in only 14 of
  50 runs (too few triggers).
- **What this does and does not show (from the second independent review):** 0 of 50 gives an upper
  limit of 0.071 and 0 of 100 of 0.036, so the two-sided check cannot flag a mask as "too
  conservative" until several hundred trials; pooled, the masks fired less often than alpha. And
  because every defect is off at once, the test cannot see a mask firing on *another* defect: in R4 the
  adaptive low-energy-cluster mask on the surface sensor keeps only 0.85-0.92 of the signal even in a
  seed with no low-energy clusters, because it fires on trail electrons. The page states both.

### Redesign of R3 and R4 after the independent review (2026-09-15)

The review (recorded in full below) showed three problems that the old R3 and R4 could not answer,
so both were rewritten and re-run; the numbers produced before this date are not used anywhere.

- **R3 ran on toy sensors and tested one side.** It now runs on the three presets of the project
  with every target defect switched off (hot columns and pixels, charge-transfer trails,
  serial-register hits, low-energy clusters, halo) while keeping dark current, spurious charge, the
  injected signal, muon tracks and high-energy deposits, so one simulation serves all five masks.
  It reports the rate per run for the stack-calibrated masks and per image for the others, and the
  check is two-sided: the 95 % interval must contain alpha, and an interval entirely below alpha is
  reported as "more conservative than alpha". 50 runs of 2 images, chosen for run time; the
  intervals are wide and the page says so.
- **The oracle grids were too narrow**, so "the best a hand-tuned mask can do" was not established:
  the optimum sat on the edge of the grid in charge-transfer trails, muons and low-energy clusters.
  All six grids were widened (the sidecar of each run lists them).
- **The headline measured the wrong thing.** A transplanted mask scoring below the oracle may simply
  be doing nothing useful. R4 now also evaluates NO MASK for every sensor and mask, and the page
  distinguishes a transplant that *harms* (scores below no mask) from one that merely does not help.
  Masks with fewer than 20 target events are excluded from the counts, which removes the
  deep-underground sensor from most comparisons; its oracle was an arbitrary tie-break.
- **The summary statistics are now pre-registered** in `PLAN.md`, written before the corrected code
  was run, and the held-out seeds (20260920, 20260921) are run once afterwards and never used to
  change anything.

### What went wrong in the re-run itself (2026-09-15), and what was done

- **Two stale seeds reached the page.** The re-run script began with `git rm` of the old R3 and R4
  results, but git refused because a short trial of R3 had rewritten a tracked result file, and the
  script carried on. The folders of seeds 20260918 and 20260919, made with the old threshold at commit
  `8041c76`, stayed on disk; the figure and the page read every seed folder they find and reported
  "7 independent seeds". Caught when checking the page after the run; the two folders were removed and
  the figure, page and PDF regenerated over the five valid seeds (commit `11cc397`). The R3 file of that
  2-run trial was also committed in `1a0c649` and replaced by the real run in `d1cf187`.
- **Every R4 sidecar said `dirty: true`.** No tracked file outside `results/` and `report/` changed
  between the frozen commit `1a0c649` and the results commit `d1cf187` (checked with `git diff`); R3
  rewrote its tracked result file first, and each later run saw that rewrite. The flag now ignores
  changes under `results/` and `report/`, which are outputs, and reports how many there were.
- **The short answer contradicted itself** ("never worse than no mask (1 of 14 cases)"), because the
  sentence was fixed text; it now follows the count and names the case. And the no-mask marker in the
  R4 figure was invisible (a line marker given the surface-coloured ring); line markers now carry
  their colour on the edge.

### R4. The masks across sensors: oracle, transplant, adaptive

- **Outputs:** `results/compare_masks/seed_<seed>/compare_masks.json` for seeds 20260915, 20260916,
  20260917 (used while developing) and 20260920, 20260921 (run once after the code was frozen at
  `1a0c649`), each with a sidecar; `results/compare_masks/relative_fom_summary.json` and the figure
  `compare_masks.png` from `analysis/figure_compare_masks.py`, which, like the page, reads only
  `seed_*` folders.
- **Produced by:** `analysis/compare_masks_across_sensors.py 4 --seed <seed>`: for each preset, an
  independent calibration stack and test stack of 4 images; each mask evaluated on the test stack as
  oracle (fixed form tuned with truth on the same sensor's calibration stack), transplant (oracle
  constants of another sensor) and adaptive (calibrated on the calibration stack without truth).
- **Written from scratch:** the evaluation (target events removed, clean pixels kept, signal
  efficiency, figure of merit), the parameter grids and the oracle search.
- **Choices with an alternative:**
  1. Figure of merit S / sqrt(S + B), S = surviving injected signal events, B = surviving target
     events (target pixels for muons). Alternatives: background rejection at fixed exposure; the
     bias of the measured single-electron rate.
  2. Oracle tuned mask by mask on grids (listed in the sidecar), not jointly, so the combined
     oracle is not the best combined mask.
  3. 4 images per stack, chosen for run time (~19 min per seed); adaptive CTI lengths depend on it.
  4. Relative figure of merit = FoM / FoM(oracle); for "all six combined" the oracle is the union of
     the six oracle masks.
- **Checks:** (a) the oracle is the best fixed form *on the calibration stack*; on the independent test
  stack a transplant or the adaptive form can exceed it, and this is reported, not hidden; (b) every
  adaptive mask passed its null and injection tests and R3 before the run; (c) variation over seeds is
  shown as a range, and the headline counts are paired: a transplant counts as harmful only if it is
  more than 1 % below no mask in every seed with at least 20 target events (the pre-registered median
  rule is reported next to it); (d) the figure was inspected by eye for hidden markers; (e) the second
  independent review recomputed every headline number from the JSON files and found them to match.
- **Grid edges (second review, fixed):** oracle optima sat on the edge of their grids for trails, muons,
  hot columns and low-energy clusters, overstating the adaptive masks relative to the oracle. The grids
  were widened (trails to 1600 x 800 px, muon charge down to 2 e and length to 1 px, hot-column factor
  from 1.0 to 30, cluster radius to 1 px, serial-register count up to 64 in windows up to 200 px,
  muon dilation up to 16 px; a first trial stopped on the serial-register count at the old maximum of 12), the edges that are physical bounds were declared in
  `PHYSICAL_BOUNDS`, and `compare_masks_across_sensors.py` now aborts, right after the oracle search
  and before any evaluation, if an optimum lies on any other edge. Each result file records
  `"oracle_edge_check": "passed"`, and the page drops the "upper bound" caveat only when every seed
  passed. R4 was re-run over the five seeds with this check.

**History of R4 before the re-run (superseded; kept because it is part of how the result was checked):**

- **A failure found by the second seed, and what was done about it.** Seed 20260916 (4 images)
  showed the adaptive masks failing catastrophically on the surface sensor: combined, 0.1 % of
  clean pixels kept (figure of merit 0.2 against the oracle's 2.2). Cause, from the adaptive
  calibration recorded in the result file: the halo radius came out at 125 px, significant up to
  the last testable annulus, with muon tracks ~90 px apart, so the halo mask covered the image; the
  hot-column calibration, which excludes the halo mask, then had no valid pixels and found no hot
  columns. The significance radius ignores the exposure a mask costs, whereas the oracle is chosen
  by figure of merit. Fix: within the significant range the radius maximises a figure of merit
  estimated from the data alone (excess of each annulus over the outside density; kept pixels),
  so the null behaviour of R3 is unchanged (no significant annulus, radius 0). (At the time R3 was
  not re-run for this change; the new radius is never larger than the significance radius and is
  0 whenever that is 0, so the halo mask can fire only in runs where the previous version fired, and
  the R3 halo rate was an upper bound. R3 has since been redesigned and re-run.) A unit test with
  dense triggers and a faint extended excess checks that the chosen radius is below the
  significance radius and masks less than half the image.
- **Truth leaking into the "adaptive" masks, found by an independent review (2026-09-15).** The
  pixel threshold was placed with `image.total.mean()`, the simulator's true mean charge, in the two
  analysis scripts, in `events.score_mask` and in the tests. That is not measurable on real data, and
  where tracks carry most of the charge it is also wrong: on the surface preset the mean is ~14
  e/pix, the density saturates at its cap and the threshold collapses to 0.50 e, so read noise alone
  produced ~1750 spurious single-electron events per image, which then fed the halo, low-energy
  cluster, serial and hot-column calibrations. Fix: `src/skmask/estimate.py` fits the noise and the
  single-electron density to the charge histogram of the image itself (two Gaussians at 0 and 1 e,
  binned Poisson likelihood), the same measurement the public SENSEI macro makes; every call site now
  uses it. Checks in `tests/test_estimate.py`: known noise and density recovered from synthetic
  images, unaffected by a track carrying thousands of electrons, and on a surface image the fitted
  density agrees with the true fraction of single-electron pixels while the threshold lands between
  0.55 and 0.85 e instead of 0.50.
  **This invalidates the R3 and R4 numbers produced before it; both are re-run.**
- **Guarding against fitting the method to the seeds it was debugged on:** seeds 20260915-17 were
  seen during development. After the fix all three are re-run, and two further seeds, 20260918 and
  20260919, were run once (later superseded: after the threshold fix the held-out seeds are
  20260920 and 20260921).
- **After the halo fix (seeds 20260915-19, code at commit `8041c76`, still with the leaking threshold):**
  no catastrophic adaptive result. This entry also claimed that "the held-out seeds fall inside the
  ranges of the development seeds"; the first independent review found that false (48 held-out values
  fell outside, 8 by more than 0.03), and all these numbers were in any case invalidated by the
  threshold fix and replaced by the re-run.
- **A finding the fix exposed, reported and not tuned away:** on the surface sensor the adaptive halo
  radius is 0 in all five seeds, so the halo mask alone reaches ~0.7 of the oracle. Cause: the
  adaptive radius estimates "signal" as all uniform single electrons, because from data alone signal
  cannot be separated from dark current; at the surface dark current is ~20 times the injected
  signal, so masking looks unfavourable. The evaluation figure of merit counts only the injected
  signal as S and does not count dark current in B. Both are choices; the evaluation was fixed
  before this was seen and is deliberately not changed now, since changing it after seeing the
  result would tune the test to the method.

### R5. What a mask fires on when the defect present is not its own

- **Output:** `results/cross_defect_false_positives/cross_defect.json`, sidecar, and figure
  `cross_defect.png` (`analysis/figure_cross_defect.py`).
- **Produced by:** `analysis/cross_defect_false_positives.py 20`: 15 configurations (the three presets
  times the five switchable defects), 20 runs of 2 images each, the same chain of adaptive
  calibrations as R3 in the same order, seed 20260916. Code frozen at commit `c013d54`.
- **Why it exists:** R3 switches every defect off at once, so it cannot tell a well-behaved mask from
  one that fires on a defect that happens to be absent too. Its own entry above records that blind
  spot, and R4 showed the symptom: on the surface sensor the adaptive low-energy-cluster mask keeps
  only 0.85-0.92 of the signal even in a seed with no low-energy clusters. Here each defect is
  switched on alone, so a mask that fires when its own defect is absent is firing on the wrong thing.
- **Written from scratch:** the configuration builder, the per-run generator seeding, the resume from
  the counts stored in the output, and the verdict rule. From libraries: the masks' own use of
  `scipy.ndimage`, and `scipy.stats.beta` for the Clopper-Pearson interval (`src/skmask/stats.py`,
  shared with R3).
- **Choices with a defensible alternative:**
  - *Muons and high-energy deposits stay on in every configuration*, as in R3. The halo is generated
    from deposited charge, so with no tracks there is no halo to switch on and that configuration
    would be empty. The alternative, switching tracks off too, would make four configurations cleaner
    and the fifth impossible.
  - *A defect is present at its preset value*, not at an exaggerated one. A stronger defect would make
    cross-talk easier to see, but the presets are what the rest of the project claims about these
    sensors.
  - *The verdict is one-sided* (the 95 % interval lies entirely above alpha). Unlike R3, the question
    here is only whether a mask fires on the wrong defect, not whether it is too conservative.
  - *The exclusion order of R3 is kept* (trails, then halo, then hot columns, then serial rows and
    clusters), so a mask upstream can protect the one after it. That is how the masks are used
    together; testing each mask alone would measure a different thing.
  - *20 runs per configuration, chosen for run time.* 0 of 20 has an upper limit of 0.17, so this test
    finds gross cross-talk and cannot resolve a rate near alpha; the page says so.
- **Expectation registered before the numbers were looked at, though the repository cannot prove it:**
  the adaptive low-energy-cluster mask on the surface sensor is expected to fire when charge-transfer
  trails are the only defect present, since that is the explanation R4 gave for its loss of signal
  days earlier. Any other cell that fires is a finding this analysis was not built to expect. The
  independent review checked the commit order and found that PLAN.md's paragraph landed at 11:37,
  after the first attempt had completed one run at about 11:29; that file was deleted unread at 11:31
  and the first R5 numbers anyone read came from the second attempt at 12:21. So this is an account,
  not a proof, and it is weaker on this point than R6, where expectation and script share one commit
  made before any result existed.
- **Result:** 20 runs of each of the 15 configurations. Of the 60 cells in which the defect present is
  not the one the mask looks for, 7 fire above the rate that mask should not exceed and 53 are
  consistent with it. The registered expectation is confirmed and is the largest effect by far: on the
  surface sensor, with charge-transfer trails as the only defect present, the adaptive
  low-energy-cluster mask fires in 40 of 40 images and masks a median 0.121 of the image. The other
  six are, in order: the serial-register mask on the same trails, 34 of 40, masking 0.0009; the
  low-energy-cluster mask on trails at shallow depth, 19 of 40; the hot-column mask on trails at the
  surface, 9 of 20; the low-energy-cluster mask on the halo at shallow depth, 10 of 40; the
  serial-register mask on low-energy clusters at shallow depth, 6 of 40; and the serial-register mask
  on trails at shallow depth, 5 of 40, which sits at its detection floor. Charge-transfer trails are
  what most often make another mask fire, which is what R4 had suspected from one symptom.
- **What the result does not show:** six of the seven mask a median of 0.001 of the image or less, so
  they are detectable but not consequential; only the surface low-energy-cluster case removes enough
  to matter. With 20 runs a cell that fires in 0 of 20 still has an upper limit of 0.17, so this
  finds gross cross-talk and says nothing about rates near alpha. And the two images of a run share
  that run's stack calibrations, so for the two per-image masks the intervals are narrower than
  independent trials would give; the page states this beside the result.
- **How it was checked:** the summary is recomputed from the counts stored in the result file, and
  `tests/test_cross_defect_resume.py` requires a run interrupted and resumed to produce the same
  numbers as an uninterrupted one, to record the two stretches honestly, to start over when the
  configurations change order, and to discard a file written with other settings.

- **Incident: a segmentation fault (2026-09-16).** The first attempt, launched while two clean-clone
  reproductions were running and about 2.2 GB of memory was free, died after one completed run with
  exit 139 and no Python traceback. Nothing was diagnosed beyond that: with no traceback the evidence
  does not separate memory exhaustion in a C extension from a defect in one. It is recorded here
  because it is the second low-level failure on this machine, after the `IndexError` of seed
  20260917, and both happened with three heavy processes running. What was done: the analysis was
  made resumable (each run is seeded from the run index and the output carries its counts, so a
  resumed run lands on the same numbers), and it was re-run with the machine otherwise idle.

### R6. The adaptive masks on the public SENSEI release

- **Output:** `results/adaptive_on_public/adaptive_on_public.json` and sidecar. No figure: four
  exposures by five masks is a table, not a chart.
- **Produced by:** `analysis/adaptive_masks_on_public_data.py`, reading the four release files through
  `src/skmask/sensei_public.py`. The expectation was registered with the first version of the script,
  commit `bb8da49`; the script was then revised twice, to record what the trail test had to work with
  and the evidence behind each flagged column, and once more after the independent review, to measure
  distances with the right pixel geometry. The shipped result comes from that last version. An earlier
  version of this entry called `bb8da49` "the frozen commit", which the result's own sidecar
  contradicted; the sidecar names the commit that was checked out while the file was written, which is
  always the parent of the commit that stores it.
- **Inputs:** the four release ROOT files (SHA-256 in the sidecar), active area only: rows 1-16 of 20
  and columns 0-3071 of 3200, the selection the release macro defines.
- **Written from scratch:** the per-image reconstruction of the active area from the flat tree, the
  overlap statistics against the release mask, and the record of what the trail test had to work with.
  From libraries: `uproot` to read the files, and everything the masks themselves call.
- **Choices with a defensible alternative:**
  - *The masks are calibrated on the raw charge, ignoring the release's own mask.* That is the point:
    the procedure has to find the defects by itself. The alternative, excluding the release-masked
    pixels first, would measure something else, namely what is left after a human mask.
  - *One stack per exposure*, not all exposures pooled: the images of one exposure share conditions,
    and pooling would mix single-electron densities that differ by a factor of five.
  - *The pairing of one of our masks with a bit of the release mask is written down, not derived*, and
    it inherits whatever R2's bit hypothesis gets wrong.
  - *The muon mask is not run.* It predicts the charge and the length of a track from the sensor's
    thickness, pixel size and back-surface diffusion; an image that bins 32 physical rows into one
    superpixel does not have tracks of that shape.
  - *The low-energy-cluster mask is run but cannot be checked*: the release publishes no counterpart.
- **Result:** given only the images, the estimator measures a readout noise of 0.141-0.142 e, where
  R1's independent fit on the same files gives 0.1413-0.1420 e, and a single-electron density that
  grows with exposure. The hot-column and hot-pixel procedure flags 2, 4, 6 and 7 columns as the
  exposure grows, and the sets are nested: each keeps the previous columns and adds the next loudest.
  At the longest exposure the seven columns with the highest rate of charged pixels are exactly the
  seven it flags, all seven lie inside the release's own bad-column mask, and the next column down
  (1659) is in that mask and is not flagged, which is the conservative direction. Of what we flag the
  release also flags 0.99-1.00; of what the release flags we flag 0.027-0.095, and our masks together
  cover 0.0007-0.0094 of the image where the release's mask covers 0.042-0.130. Both directions are on
  the page now: the first is the easy one, since a procedure that flagged a single true column would
  score 1.00, and it rests on 2 to 7 column decisions rather than on the hundreds of pixels it is
  counted over. The halo radius is 0 at the three shorter exposures and 25 column widths at 20 hours.
- **The registered expectation was wrong about the trails, and the result file says why.** PLAN.md,
  written before the run, expected a non-zero horizontal trail length because the release masks
  "bleeding". The procedure chose zero on all four exposures. The cause is not a defect in the
  procedure but the data: the release blinds its hits, so a whole exposure holds 0, 70, 103 and 844
  pixels above the 20 e trigger. At 20 hours the test does see the excess on the correct side, 13
  single electrons within 5 superpixels downstream of a trigger against 2 upstream, p = 0.0037, but
  the Bonferroni-corrected threshold is 6.3e-05. The script now records the trigger counts and the
  smallest p-value per direction, so this statement is read from the result file rather than typed.
- **How it was checked:** (i) the noise measured here image by image scatters over 0.1392-0.1456 e and
  its median per exposure is 0.141-0.142 e, which agrees to 0.001 e with the 0.1413-0.1420 e that R1
  fits per exposure on the same files through a separate code path; (ii) the flagged sets are nested
  across the four exposures, 2 then 4 then 6 then 7 columns, each keeping the previous ones, which is
  what accumulating evidence should look like and what a spurious flag would not do (an earlier
  version of this entry claimed the sets were "stable" and named three columns that recur, which the
  independent review correctly called a cherry-pick: at the time one column was flagged at three
  exposures and lost at the fourth, for the reason recorded below); (iii) each flagged column carries
  102 to 1830 times the rate of charged low-charge pixels of the columns left alone, the same counts
  the calibration works with, so these are not marginal flags (recorded per column in the result
  file); (iv) the eight loudest columns of every exposure are recorded with the rate that makes them
  loud and whether the procedure flagged them, so a column it misses is visible rather than absent.
- **What it does not show:** the single-electron densities here are measured on unmasked data and are
  therefore larger than R1's, which are fitted after the release mask; they are not the same quantity
  and the two numbers must not be compared as if they were. And with four masks checkable on one real
  sensor, this is a demonstration that the procedure survives contact with real data, not a
  measurement of how well it does in general.

### R7. The muon mask with no tracks at all

- **Output:** `results/muon_mask_null_rate/muon_null.json` and sidecar. No figure: three numbers.
- **Produced by:** `analysis/muon_mask_null_rate.py 100`, 100 runs of 2 images on each of the three
  presets with every defect off and, unlike R3 and R5, the muon flux and the high-energy background
  set to zero as well. Seed 20260917, one generator per sensor seeded by (seed, index) rather than by
  the sensor's name, since `hash` of a string differs between processes.
- **Why it exists:** R3 and R5 keep tracks on, because a real sensor cannot switch them off and
  because the halo is generated from their charge, and on the public release of R6 the muon mask is
  not run at all. So the sixth mask had never been asked the question the other five were asked. This
  asks it in the only place tracks can be absent, the simulator.
- **Choices with a defensible alternative:** the high-energy background is switched off with the
  muons, because it produces the same kind of compact deposit the mask could mistake for a short
  track; leaving it on would measure something else. The mask's threshold is physical, so there is no
  alpha to compare the rate with, and the result is reported as a rate with its interval rather than
  as a verdict.
- **Result:** 0 firings in 200 images on each of the three sensors, 0 of 600 in all, an upper limit of
  0.018 at 95 % per sensor, and no pixel masked in any image. The expectation registered in PLAN.md
  before the run said zero, for a stated reason: the mask needs a cluster carrying about the charge a
  minimum-ionising particle leaves crossing the sensor, of order a hundred electrons, and what is left
  in these images is single electrons from dark current, spurious charge and the injected signal.
- **What it does not show:** 0 of 200 cannot exclude a rate below about 0.018, and nothing here says
  what the mask does to a real track it should catch, which is R4's business. It also does not test
  the mask against the one thing that could plausibly imitate a track and is not in these images: a
  long, thin defect such as a partially charged column, which the hot-column mask removes first in
  the chain but which would reach the muon mask if the hot-column calibration missed it.

## Second independent review (2026-09-16), of R5 and R6

A second agent reviewed the new work with no access to this conversation. Its findings are recorded
here with what was done about each, including the ones that were not acted on.

**Accepted and fixed:**

- **The halo mask was measuring distance in superpixels on an image whose superpixels bin 32 rows.**
  `adaptive_halo_radius` and `halo_mask` both take a `sampling` argument and R6 passed neither, so a
  radius of 15 meant 15 columns by 15 superpixel rows, which on a 16-row image is the whole height of
  the frame. The consequence was concrete and had gone unnoticed: at the longest exposure the halo
  mask covered every pixel of column 134, the loudest column in the image, so the hot-column
  calibration could not see it and did not flag it, while the earlier text called the flagged set
  "stable". R6 now passes `sampling = (32, 1)`; the halo radius becomes 0, 0, 0 and 25 column widths,
  and the flagged sets become nested, 2, 4, 6 and 7 columns, with column 134 present throughout. This
  is also the cross-defect failure R5 was built to find, found on real data by a reviewer instead.
- **"Fraction masked" on the page was the sum of five overlapping masks.** It is now their union,
  computed in the analysis, and the table also prints what the release's own mask covers, which the
  page never said.
- **The page told only the flattering direction of the agreement with the release.** It gave "of what
  we flag, the release also flags 0.99-1.00" and left the recall, 0.027-0.095, inside a collapsed
  table. Both are now in the prose, with the sizes of the two masks and a sentence saying why the
  first number is the easy one.
- **A number on the page came from a Python constant, not from a result file.** The measured noise was
  compared against `PRESETS["deep_underground"].noise_e` while the page attributed it to the release
  paper. It is now compared against R1's fitted noise, read from `release_rate.json`.
- **"Code frozen at commit `bb8da49`" was false**, as the result's own sidecar showed. Corrected above.
- **The correction in the trail test counts distance blocks that cannot exist** on a 16-row image: 80
  per direction, of which three can ever be filled vertically. The page now says so and gives the
  threshold a correct count would ask for, which this p-value still does not meet, so the decision
  stands.

**Accepted, recorded, not fixed:**

- **The R5 pre-registration cannot be demonstrated by the repository.** PLAN.md's R5 paragraph was
  committed at 11:37, after the first attempt at R5 had completed one run and written it to disk at
  about 11:29. The file was deleted unread at 11:31, when the analysis was made resumable, and the
  first R5 number that was looked at was at 12:21, from the second attempt; the expectation itself
  comes from R4, which had recorded the low-energy-cluster mask firing on trail electrons days
  earlier. That is the honest account, and it is an account: the repository can show only the commit
  order, so R5's expectation is weaker evidence than R6's, where the expectation and the script land
  in one commit before any result exists. The entry above has been reworded accordingly.
- **R5's per-image trials within a run are not independent** (the two images of a run share the
  stack-calibrated exclusions), so its intervals are narrower than they should be for the two
  per-image masks. Measuring it properly needs a re-run of about 100 minutes, which did not fit; the
  page and the entry state it.

## The deliverable pages

- **`report/index.html` is generated**, never edited: `analysis/build_report.py` reads the result
  JSONs and writes the page, so no number on it is typed by hand and every comparison in its prose
  ("worse than no mask in every seed", the counts, the medians) is computed from the same files. Its
  sidecar lists every input with a hash. A section whose result file is absent simply does not appear.
- **`report/report.pdf` is that page printed** by `scripts/make_pdf.py` with headless Chrome, so the
  two cannot disagree; its sidecar records the browser and the hash of the page it printed.
- **`index.html` at the root of the repository is written by hand** and carries no result: it is the
  front door GitHub Pages serves, with links to the report, the PDF and the two documents. `.nojekyll`
  next to it tells Pages to serve the files as they are instead of building them.
- **Authors:** D. Rodriguez Juiz and F. Perez.
