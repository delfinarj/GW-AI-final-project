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
- **What the macro establishes about the files:** tree `calPixTree` with branches `x`, `y`,
  `ePix` (charge in electrons) and `mask`; active rows `0 < y <= 16`; mask selection `0x067d`
  for the single-electron analysis; images of 3200 columns by 20 rows. To be confirmed by
  reading the files directly, and recorded here when done.

## Results

(none yet)
