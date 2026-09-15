"""From measured charge to electrons and events, and scoring of a mask against ground truth.

Everything here is shared by all masks and must itself be sensor-independent: the pixel
threshold is derived from quantities measured in the image (noise, single-electron density)
rather than fixed.
"""
import numpy as np
from scipy import ndimage

BACKGROUND_ORIGINS = ("hot_column", "hot_pixel", "muon", "highE", "halo", "serial",
                      "lowE_cluster", "cti")
CLEAN_ORIGINS = ("dark", "signal")


def one_electron_threshold(noise_e, density):
    """Bayes-optimal cut between 0 and 1 electron for Gaussian noise.

    With prior probabilities (1 - mu) for 0 e and mu for 1 e and equal-width Gaussians of sigma,
    the two posteriors are equal where
        (1 - mu) exp(-c^2 / 2 sigma^2) = mu exp(-(c - 1)^2 / 2 sigma^2),
    i.e.  c = 1/2 + sigma^2 ln((1 - mu) / mu).
    Derived here from scratch; for sigma = 0.14 e and mu = 1e-4 it gives c = 0.68 e.
    """
    density = float(np.clip(density, 1e-12, 0.5))
    return 0.5 + noise_e ** 2 * np.log((1.0 - density) / density)


def to_electrons(measured, noise_e, density):
    """Integer electrons per pixel: 1 e above the Bayes cut, n e above n - 1/2 for n >= 2."""
    c1 = one_electron_threshold(noise_e, density)
    electrons = np.where(measured >= 1.5, np.floor(measured + 0.5), 0.0)
    electrons = np.where((measured >= c1) & (measured < 1.5), 1.0, electrons)
    return electrons.astype(np.int64)


def clusters(electrons, connectivity=8):
    """Label contiguous non-empty pixels. Returns (labels, table) with one row per cluster:
    charge, n_pix, y_min, y_max, x_min, x_max, y_centroid, x_centroid (charge-weighted)."""
    structure = np.ones((3, 3)) if connectivity == 8 else ndimage.generate_binary_structure(2, 1)
    labels, n = ndimage.label(electrons > 0, structure=structure)
    if n == 0:
        return labels, np.zeros((0, 8))
    index = np.arange(1, n + 1)
    charge = ndimage.sum_labels(electrons, labels, index)
    n_pix = ndimage.sum_labels(electrons > 0, labels, index)
    boxes = ndimage.find_objects(labels)
    y_min = np.array([b[0].start for b in boxes]); y_max = np.array([b[0].stop - 1 for b in boxes])
    x_min = np.array([b[1].start for b in boxes]); x_max = np.array([b[1].stop - 1 for b in boxes])
    centroid = np.array(ndimage.center_of_mass(electrons, labels, index)).reshape(-1, 2)
    table = np.column_stack([charge, n_pix, y_min, y_max, x_min, x_max, centroid[:, 0], centroid[:, 1]])
    return labels, table


def single_electron_events(electrons):
    """Boolean map of isolated one-pixel, one-electron events (no charged 8-neighbour)."""
    occupied = electrons > 0
    neighbours = ndimage.convolve(occupied.astype(np.int64), np.ones((3, 3), dtype=np.int64),
                                  mode="constant") - occupied
    return (electrons == 1) & (neighbours == 0)


def pixel_origin(charge):
    """Dominant true origin of each pixel (name array), '' where no charge."""
    names = list(charge)
    stack = np.stack([charge[name] for name in names])
    winner = np.argmax(stack, axis=0)
    origin = np.array(names, dtype=object)[winner]
    origin[stack.sum(axis=0) == 0] = ""
    return origin


def score_mask(mask, image, electrons=None):
    """How well a boolean mask (True = discarded) does against the simulator's truth.

    Returns a dict with, per origin, the number of single-electron events and the fraction
    removed; the fraction of clean pixels kept (pixels whose true charge comes only from dark
    current and signal, including empty ones); and the signal efficiency.
    """
    if electrons is None:
        from .estimate import electrons_from_image
        electrons, _ = electrons_from_image(image.measured)
    events = single_electron_events(electrons)
    origin = pixel_origin(image.charge)
    result = {"per_origin": {}}
    for name in image.charge:
        at = events & (origin == name)
        n = int(at.sum())
        result["per_origin"][name] = {"events": n,
                                      "removed_fraction": float((at & mask).sum() / n) if n else None}
    background = sum(image.charge[name] for name in BACKGROUND_ORIGINS if name in image.charge)
    clean = background == 0
    result["clean_pixels_kept"] = float((clean & ~mask).sum() / clean.sum())
    signal = result["per_origin"].get("signal", {})
    result["signal_efficiency"] = (1.0 - signal["removed_fraction"]) if signal.get("removed_fraction") is not None else None
    result["masked_fraction"] = float(mask.mean())
    return result
