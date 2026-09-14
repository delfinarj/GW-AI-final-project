"""Reader for the public SENSEI SNOLAB release (github.com/sensei-skipper/DataReleases).

What the release macro `plotRate.C` establishes, and this module relies on:
- tree `calPixTree`, one entry per superpixel, branches `x`, `y`, `ePix` (electrons), `mask`;
- images of 3200 columns by 20 rows, active area 0 < y <= 16 and x < 3072 (the rest is overscan);
- hardware binning of 32 physical rows per superpixel, so per-pixel rates are per-superpixel / 32;
- mask selection 0x067d for the single-electron analysis;
- readout of a whole image takes 965 s, and a superpixel read at (x, y) has accumulated an
  extra exposure of (965 / 20) * (y + x / 3200) seconds.

Confirmed by reading the files (2026-09-14): branches x, y, ePix, mask, ohdu, RUNID, LTANAME;
one quadrant (ohdu 2); 14, 14, 14 and 13 images for 0, 2, 6 and 20 h, 64000 entries each.
"""
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "public" / "sensei_snolab_binned_1"
EXPOSURES_S = {0: "RELEASE_hits_blinded_EXP0_13.root",
               7200: "RELEASE_hits_blinded_EXP7200_13.root",
               21600: "RELEASE_hits_blinded_EXP21600_13.root",
               72000: "RELEASE_hits_blinded_EXP72000_13.root"}

TREE = "calPixTree"
N_COLUMNS = 3200
N_ROWS = 20
ACTIVE_COLUMNS = 3072
ACTIVE_ROWS = 16
ROWS_PER_SUPERPIXEL = 32
MASK_1E_ANALYSIS = 0x067D
IMAGE_READOUT_S = 965.0


def open_tree(exposure_s):
    import uproot
    return uproot.open(DATA_DIR / EXPOSURES_S[exposure_s])[TREE]


def branches(exposure_s):
    tree = open_tree(exposure_s)
    return {name: str(tree[name].interpretation) for name in tree.keys()}, tree.num_entries


def load(exposure_s, extra=()):
    """All superpixels of one exposure as a dict of numpy arrays."""
    tree = open_tree(exposure_s)
    names = ["x", "y", "ePix", "mask", *extra]
    return tree.arrays(names, library="np")


def active_unmasked(data, mask_bits=MASK_1E_ANALYSIS):
    """The selection used by the release macro: active rows and no bit of `mask_bits` set."""
    return (data["y"] > 0) & (data["y"] <= ACTIVE_ROWS) & ((data["mask"].astype(np.int64) & mask_bits) == 0)


def extra_exposure_days(data):
    """Exposure accumulated during readout, per superpixel, in days (macro's formula)."""
    return IMAGE_READOUT_S / N_ROWS * (data["y"] + data["x"] / N_COLUMNS) / 86400.0
