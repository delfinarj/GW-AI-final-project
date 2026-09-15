"""Read noise and single-electron density measured from an image, without any truth.

The pixel threshold that turns measured charge into electrons depends on the read noise and on the
single-electron density (`events.one_electron_threshold`). Both are properties of the image and are
measured here the way the public SENSEI release macro measures them: the charge histogram around
zero is two Gaussians of equal width, at 0 and 1 electron, whose amplitude ratio is the density.

An earlier version passed the simulator's true mean charge as the density. That is truth leaking
into a mask that claims to calibrate itself, and it is also wrong where tracks carry most of the
charge: on the surface preset the mean is ~14 e/pix, which drove the threshold to 0.5 e and turned
read-noise fluctuations into single-electron events.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

FIT_RANGE = (-1.0, 1.5)
HIST_RANGE = (-1.1, 1.6)
HIST_BINS = 400


def estimate_noise_and_density(measured, bins=HIST_BINS, hist_range=HIST_RANGE, fit_range=FIT_RANGE):
    """Fit (noise, density) to the 0e and 1e peaks of a measured charge image.

    Returns a dict with `noise_e`, `density`, the fitted normalisation and whether the fit converged.
    Pixels far above 1 e (tracks, defects) are outside the histogram range and do not enter the fit.
    """
    values = np.asarray(measured).ravel()
    counts, edges = np.histogram(values, bins=bins, range=hist_range)
    centres = 0.5 * (edges[1:] + edges[:-1])
    width = edges[1] - edges[0]
    use = (centres >= fit_range[0]) & (centres <= fit_range[1])
    k = counts[use].astype(float)
    x = centres[use]

    def expected(params):
        log_n, noise, log_density = params
        density = np.exp(log_density)
        shape = (1 - density) * norm.pdf(x, 0.0, noise) + density * norm.pdf(x, 1.0, noise)
        return np.exp(log_n) * width * shape

    def negative_log_likelihood(params):
        lam = np.maximum(expected(params), 1e-300)
        return float(np.sum(lam - k * np.log(lam)))

    start = (np.log(max(k.sum(), 1.0)), 0.15, np.log(1e-3))
    fit = minimize(negative_log_likelihood, start, method="Nelder-Mead",
                   options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 4000})
    log_n, noise, log_density = fit.x
    return {"noise_e": float(abs(noise)), "density": float(np.exp(log_density)),
            "normalisation": float(np.exp(log_n)), "converged": bool(fit.success)}


def electrons_from_image(measured, bins=HIST_BINS, hist_range=HIST_RANGE, fit_range=FIT_RANGE):
    """Integer electrons per pixel, with the threshold set from the image itself."""
    from .events import to_electrons

    fit = estimate_noise_and_density(measured, bins, hist_range, fit_range)
    return to_electrons(measured, fit["noise_e"], fit["density"]), fit
