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

**Check:** `tests/test_simulator.py` compares the simulator with answers derived independently
of the code: the diffusion width at a depth, the minimum-ionising charge per length, the Poisson
mean of uniform rates, the muon count and mean cos(theta) = 4/5 for a cos^2 intensity through a
horizontal plane, the mean projected halo distance (pi/4 of the absorption length), CTI trails
only downstream with the set mean length, and serial hits confined to single rows.

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
