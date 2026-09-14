"""Reproduce the single-electron rate of the public SENSEI SNOLAB release.

This checks that the data are read correctly and that the exposure bookkeeping is right before
anything else is built on them. It re-implements the release macro `plotRate.C` in Python:

1. per exposure, histogram the charge of active, unmasked superpixels (mask & 0x067d == 0) with
   400 bins on [-1.1, 1.6];
2. fit a two-Gaussian model on [-1.0, 1.5] by binned Poisson likelihood (the macro's "L" option):
   N [(1 - mu) G(x; 0, s) + mu G(x; 1, s)], with the gain fixed to 1; mu is the 1e density;
3. exposure = nominal exposure + mean readout exposure of the selected superpixels;
4. straight-line fit of density against exposure; slope / 32 is the rate per pixel per day.

Published value (arXiv:2410.18716, Golden quadrant): (1.39 +- 0.11) e-5 e-/pix/day.

Run:  python analysis/reproduce_release_rate.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

from skmask import sensei_public as sp
from skmask.provenance import write_sidecar

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "release_rate"
PUBLISHED = (1.39e-5, 0.11e-5)


def fit_density(charge):
    counts, edges = np.histogram(charge, bins=400, range=(-1.1, 1.6))
    centres = 0.5 * (edges[1:] + edges[:-1])
    width = edges[1] - edges[0]
    use = (centres >= -1.0) & (centres <= 1.5)
    k = counts[use]
    xc = centres[use]

    def expected(params):
        log_n, noise, mu = params
        return np.exp(log_n) * width * ((1 - mu) * norm.pdf(xc, 0.0, noise) + mu * norm.pdf(xc, 1.0, noise))

    def nll(params):
        lam = np.maximum(expected(params), 1e-300)
        return np.sum(lam - k * np.log(lam))

    start = (np.log(k.sum()), 0.2, 1e-3)
    fit = minimize(nll, start, method="L-BFGS-B",
                   bounds=[(None, None), (0.01, 1.0), (0.0, 1.0)])
    # uncertainty from the numerical Hessian of the negative log-likelihood
    hessian = _numerical_hessian(nll, fit.x)
    covariance = np.linalg.inv(hessian)
    return {"mu": fit.x[2], "mu_err": float(np.sqrt(covariance[2, 2])),
            "noise": fit.x[1], "n_pixels": int(len(charge)), "converged": bool(fit.success)}


def _numerical_hessian(f, x, rel_step=1e-4):
    x = np.asarray(x, dtype=float)
    n = len(x)
    h = np.maximum(np.abs(x), 1e-6) * rel_step
    hess = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            e_i = np.zeros(n); e_i[i] = h[i]
            e_j = np.zeros(n); e_j[j] = h[j]
            hess[i, j] = (f(x + e_i + e_j) - f(x + e_i - e_j) - f(x - e_i + e_j) + f(x - e_i - e_j)) / (4 * h[i] * h[j])
    return hess


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for exposure_s in sorted(sp.EXPOSURES_S):
        data = sp.load(exposure_s)
        selected = sp.active_unmasked(data)
        fit = fit_density(data["ePix"][selected])
        readout = sp.extra_exposure_days(data)[selected]
        rows.append({"exposure_s": exposure_s,
                     "exposure_days": exposure_s / 86400.0 + readout.mean(),
                     "exposure_days_sd": float(readout.std()),
                     "n_images": data["x"].size / (sp.N_COLUMNS * sp.N_ROWS),
                     **fit})
        print(f"{exposure_s:>6d} s  density {fit['mu']:.3e} +- {fit['mu_err']:.1e} e-/superpix  "
              f"noise {fit['noise']:.3f} e-  unmasked {fit['n_pixels']}")

    exposure = np.array([r["exposure_days"] for r in rows])
    density = np.array([r["mu"] for r in rows])
    error = np.array([r["mu_err"] for r in rows])
    weights = 1.0 / error ** 2
    design = np.column_stack([np.ones_like(exposure), exposure])
    covariance = np.linalg.inv(design.T @ (design * weights[:, None]))
    intercept, slope = covariance @ design.T @ (weights * density)
    rate = slope / sp.ROWS_PER_SUPERPIXEL
    rate_err = np.sqrt(covariance[1, 1]) / sp.ROWS_PER_SUPERPIXEL
    chi2 = float(np.sum(weights * (density - design @ [intercept, slope]) ** 2))

    pull = (rate - PUBLISHED[0]) / np.hypot(rate_err, PUBLISHED[1])
    results = {"rate_e_per_pix_day": rate, "rate_err": rate_err,
               "intercept_e_per_superpix": intercept, "chi2": chi2, "dof": len(rows) - 2,
               "published_rate": PUBLISHED[0], "published_err": PUBLISHED[1], "pull": pull,
               "per_exposure": rows}
    print(f"rate {rate:.3e} +- {rate_err:.1e} e-/pix/day   published {PUBLISHED[0]:.2e} +- {PUBLISHED[1]:.1e}"
          f"   pull {pull:+.2f}")

    output = OUT_DIR / "release_rate.json"
    output.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
    write_sidecar(output, __file__,
                  inputs=[sp.DATA_DIR / name for name in sp.EXPOSURES_S.values()],
                  parameters={"mask_bits": hex(sp.MASK_1E_ANALYSIS), "hist_bins": 400,
                              "hist_range": [-1.1, 1.6], "fit_range": [-1.0, 1.5]},
                  results={"rate": rate, "rate_err": rate_err, "pull": pull},
                  notes="Python re-implementation of the release macro plotRate.C")


if __name__ == "__main__":
    main()
