# Plan: sensor-independent pixel masks for Skipper-CCD images

Deadline: 2026-09-30. Deliverables: an HTML page and a PDF (English), and this repository,
reproducible end to end by an agent that has never spoken to the author.

## The question

Skipper-CCD analyses discard pixels with masks. Each mask targets one problem:

| Family | Masks studied here |
|---|---|
| Deficient readout | charge-transfer inefficiency (CTI, "bleed"); serial-register hits |
| Material defects | hot columns and hot pixels |
| Unwanted physics | halo around high-energy events; low-energy clusters (LEC); muon tracks |

**Can each of these masks be made efficient at its own problem without depending on the
particular Skipper-CCD it runs on?**

## Working hypothesis

Published mask definitions share procedures but not numbers: a radius, a trail length or a
column-rate cut is tuned by hand for one sensor, readout configuration and site, and a value
that is right for one detector is wrong for another. What does transfer between sensors is the
*procedure* used to choose the number - for example "grow the mask until the single-electron
rate stops changing", or "cut at a p-value corrected for the number of tests".

So a sensor-independent mask is defined here as a **procedure that calibrates itself from the
images it is applied to**, using only quantities measurable in those images (read noise,
single-electron rate, diffusion, trail profiles), with thresholds set by a **controlled
false-positive rate** rather than by a hand-chosen constant.

The test: run the same procedure, with no retuning, on sensors that differ, and compare it with
(a) a mask tuned for each sensor with ground truth ("oracle") and (b) a mask tuned on one sensor
and transplanted unchanged to another (the failure mode the question is about).

## The six masks

For each: the **fixed** reference (hand-tuned constants, as in the public literature) and the
**adaptive** candidate.

1. **Halo.** Fixed: disc of radius R around pixels above a charge threshold. Adaptive: measure
   the 1e-event density as a function of distance from high-energy events, and set R where the
   excess over the far-field density is consistent with zero at level alpha.
2. **CTI / bleed.** Fixed: L pixels trailing along the row and column readout directions.
   Adaptive: the same profile measurement, separately along each readout direction, since CTI
   is directional and halo is not.
3. **Hot columns and hot pixels.** Fixed: column rate above k times the median, or k MADs.
   Adaptive: per-column count versus per-column unmasked exposure, Poisson survival function
   with the cut corrected for the number of columns, iterated until stable, with neighbour
   binning for warm pairs.
4. **Serial-register hits.** Fixed: window width and event count per row. Adaptive: row-window
   occupancy compared with the Poisson expectation from the measured single-electron rate,
   corrected for the number of windows; plus a shape cut derived from the diffusion model (a
   bulk event cannot be one row tall and many columns wide).
5. **Low-energy clusters.** Fixed: disc around every cluster. Adaptive: nearest-neighbour
   statistic D_j = sum 1/r^2 over the n nearest events, p-values from a Monte Carlo of uniform
   events over the unmasked area, regional combination with Fisher's method and a binomial
   excess criterion.
6. **Muon tracks.** Fixed: charge and length cuts on clusters. Adaptive: morphology that scales
   with the sensor rather than with pixel counts - straightness (principal-axis eigenvalue
   ratio), length consistent with the sensor thickness and incidence angle, uniform charge per
   unit length - with the threshold set on simulated tracks of the sensor at hand; plus a
   per-pixel temporal Poisson-outlier test for image stacks.

## Simulator (ground truth)

A NumPy image simulator in `src/`, parametrised by sensor properties so that "different
sensors" is a configuration change:

- geometry: columns, rows, pixel size, thickness, binning;
- readout: noise, single-electron rate (dark current plus spurious charge), CTI probabilities
  along both directions, serial-register hit rate;
- defects: hot columns and hot pixels with their rates;
- backgrounds: muons (rate and angular distribution), high-energy electrons, halo 1e events
  around high-energy deposits, clustered low-energy events;
- diffusion: sigma_xy(z) = sqrt(-A ln|1 - b z|), A = 218.715 um^2, b = 1.015e-3 /um
  (SENSEI, arXiv:2004.11378), with z measured from the back of the sensor;
- a uniform injected signal, to measure signal efficiency.

Every simulated charge carries a **label of its origin**, so each mask's efficiency against each
source is measured directly. Sensor presets are built from **public** parameters only, each with
its citation: a deep-underground dark-matter-like sensor, a shallow-site sensor, and a surface
laboratory sensor with a high muon rate (the regime of the author's own measurements).

## Metrics

Per mask, per preset:

- **background rejection**: fraction of 1e events from the targeted source that are masked;
- **exposure kept**: fraction of clean pixels left unmasked;
- **signal efficiency**: fraction of injected uniform signal events kept;
- **rate bias**: estimated single-electron rate after masking minus the true clean rate;
- **null false-positive rate**: with the targeted defect switched off, the fraction masked
  must be consistent with alpha.

## Check against real data

The public SENSEI SNOLAB release (`scripts/fetch_public_data.py`) ships binned images with
per-pixel mask bits. Two checks:

1. Reproduce the single-electron rate the release macro computes, which validates reading the
   data and the exposure bookkeeping.
2. Run the adaptive hot-column, halo, CTI and edge procedures on the same images and compare
   with the shipped bits: overlap, exposure kept, and resulting rate.

## Provenance and checks

- Every figure and number in the deliverables is produced by a script in `analysis/`, writes a
  sidecar record (script, git commit, inputs with hashes, parameters, random seed), and is listed
  in `PROVENANCE.md` together with **how it was checked**.
- `tests/` holds checks with known answers: simulator statistics (Poisson rates, diffusion
  widths, CTI trail lengths), null behaviour of every adaptive mask, and injected-defect
  recovery.

## Changes to this plan while doing the work

Recorded here so the plan and the code do not disagree; the reasons are in `PROVENANCE.md`.

- **Halo:** the reference is everything outside each annulus, not a far field (no far field exists
  at the surface), and within the significant range the radius maximises a figure of merit
  estimated from the data (the significance radius alone masked a whole surface image in one seed).
- **Calibration order:** charge-transfer trails, then halo, then hot columns and pixels (vertical
  trails of muon tracks are a column-wise excess).
- **Muons:** clauses for crossing tracks (pile-up) and for tracks cut by the image border.
- **Low-energy clusters:** a neighbour count against a Poisson expectation for the local valid area,
  at several radii, instead of the nearest-neighbour statistic with regional combination.
- **R4 statistics:** five seeds, three seen during development and two held out.

## Schedule

| Dates | Work |
|---|---|
| Sep 14-16 | environment, data fetch, public-data reader and rate reproduction, simulator core |
| Sep 17-20 | simulator backgrounds and presets; tests of the simulator |
| Sep 21-25 | six masks, fixed and adaptive; metrics across presets; real-data comparison |
| Sep 26-28 | figures, provenance review, independent re-check by a second agent |
| Sep 29-30 | HTML page, PDF, final reproduction from a clean clone |
