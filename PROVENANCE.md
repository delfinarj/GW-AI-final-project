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
- **Produced by:** `analysis/null_false_positive_rates.py`, 200 independent runs (3 images each)
  of two defect-free 800 x 400 sensors: one with only dark current and signal, one adding sparse
  2-7 keV deposits so that the halo and CTI calibrations have triggers. Seed 20260914.
- **Code state:** started at commit `26acd81` (clean), which contains the halo redesign; the
  sidecar also records the state when it was written (a later commit with uncommitted edits),
  which does not affect the run.
- **Written from scratch:** the run loop and the Clopper-Pearson interval (`scipy.stats.beta`).
- **Choice with an alternative:** "fires" is per run (any flag at all in the stack), the strictest
  reading; per-pixel or per-column false-positive fractions would be smaller.
- **Result (fired / runs, 95 % interval):** hot columns 1/200 [0.000, 0.028]; CTI 0/200
  [0.000, 0.018]; halo 2/200 [0.001, 0.036]; serial register 0/200 [0.000, 0.018]; low-energy
  clusters 1/200 [0.000, 0.028].
- **Check:** every interval contains alpha = 0.01, so each adaptive mask's false-positive rate is
  consistent with the rate it was built to have. The earlier far-field version of the halo test is
  the one that failed this kind of check (a single realisation gave radius 45 with no halo).
  Limitation: 200 runs cannot distinguish 0.01 from 0.03; the muon mask is not included because it
  has no false-positive parameter (it is tested against point-like deposits in `tests/test_masks.py`).

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

### R4. The masks across sensors: oracle, transplant, adaptive

- **Outputs:** `results/compare_masks/compare_masks.json` (seed 20260915) and
  `results/compare_masks/seed_<seed>/compare_masks.json` (further seeds), each with a sidecar;
  `results/compare_masks/relative_fom_summary.json` (median and range over seeds) and the figure
  `compare_masks.png`, both from `analysis/figure_compare_masks.py`.
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
- **Checks:** (a) the oracle is by construction the best fixed form on the calibration stack; on the
  independent test stack the adaptive form can exceed it (e.g. serial-register hits at the surface in
  seed 20260915), which is possible because the oracle is chosen on a different stack and on a grid,
  and is reported, not hidden, in the per-mask table of the page; (b) every adaptive mask passed its null and injection tests and R3 before
  this run; (c) the variation over seeds is reported as a range, and the statements on the page
  ("below half of the oracle in n of N cases") are computed from the per-seed minimum, the least
  favourable reading; (d) the figure was inspected by eye for hidden markers (an earlier version
  hid transplant markers under the oracle line; series are now offset within each row).
- **Every number on the page** is read from these files by `analysis/build_report.py`; none is typed.
- **A failure found by the second seed, and what was done about it.** Seed 20260916 (4 images)
  showed the adaptive masks failing catastrophically on the surface sensor: combined, 0.1 % of
  clean pixels kept (figure of merit 0.2 against the oracle's 2.2). Cause, from the adaptive
  calibration recorded in the result file: the halo radius came out at 125 px, significant up to
  the last testable annulus, with muon tracks ~90 px apart, so the halo mask covered the image; the
  hot-column calibration, which excludes the halo mask, then had no valid pixels and found no hot
  columns. The significance radius ignores the exposure a mask costs, whereas the oracle is chosen
  by figure of merit. Fix: within the significant range the radius maximises a figure of merit
  estimated from the data alone (excess of each annulus over the outside density; kept pixels),
  so the null behaviour of R3 is unchanged (no significant annulus, radius 0). R3 was run before
  this change and is not re-run: the new radius is never larger than the significance radius and is
  0 whenever that is 0, so the halo mask can fire only in runs where the previous version fired, and
  the R3 halo rate is an upper bound for the current code. A unit test with
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
  20260919, are run once, with no change to the code in between, and reported separately.
- **After the fix (all five seeds, code at commit `8041c76`):** no catastrophic adaptive result on any
  sensor or seed; on the surface sensor the six adaptive masks combined reach a figure of merit of
  2.0-2.3 against 1.8-2.7 for the oracle and 0.0-1.8 for transplanted constants. The held-out seeds
  fall inside the ranges of the development seeds.
- **A finding the fix exposed, reported and not tuned away:** on the surface sensor the adaptive halo
  radius is 0 in all five seeds, so the halo mask alone reaches ~0.7 of the oracle. Cause: the
  adaptive radius estimates "signal" as all uniform single electrons, because from data alone signal
  cannot be separated from dark current; at the surface dark current is ~20 times the injected
  signal, so masking looks unfavourable. The evaluation figure of merit counts only the injected
  signal as S and does not count dark current in B. Both are choices; the evaluation was fixed
  before this was seen and is deliberately not changed now, since changing it after seeing the
  result would tune the test to the method.
