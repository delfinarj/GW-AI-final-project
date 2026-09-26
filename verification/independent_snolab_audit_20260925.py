"""Independent audit of the public SENSEI SNOLAB binned release.

This deliberately imports no project code.  It was written before reading
``analysis/`` or ``src/``.  The only methodological inputs from the release are
the public README and plotRate.C: the tree layout, active geometry, mask word,
32 physical pixels per released superpixel, and readout timing.

Run with:
    uv run python verification/independent_snolab_audit_20260925.py
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from scipy.stats import binom


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "public" / "sensei_snolab_binned_1"
OUTPUT = Path(__file__).with_suffix(".json")

FILES = [
    ("RELEASE_hits_blinded_EXP0_13.root", 0.0),
    ("RELEASE_hits_blinded_EXP7200_13.root", 7_200.0),
    ("RELEASE_hits_blinded_EXP21600_13.root", 21_600.0),
    ("RELEASE_hits_blinded_EXP72000_13.root", 72_000.0),
]
EXPECTED_SHA256 = {
    "RELEASE_hits_blinded_EXP0_13.root": "e38824424e019e2e5096786e8fa38834da7c2bf61a6087214fe0ec11df02a585",
    "RELEASE_hits_blinded_EXP7200_13.root": "21d7294f4238dd83450b9aec56178599ff4eabf62729aeb16a94621b6b58e712",
    "RELEASE_hits_blinded_EXP21600_13.root": "9864c9e342b34077759560ee4240faab42f57a8c960f8e832ee517cd71d0f7e6",
    "RELEASE_hits_blinded_EXP72000_13.root": "79a09c60ecf8e05ca2f742b7816985cf5a31b3ce51b9f925c96f0e47fe4624fb",
}

ACTIVE_COLUMNS = 3_072
ACTIVE_ROWS = (1, 16)
RELEASE_MASK_WORD = 0x067D
PHYSICAL_PIXELS_PER_SUPERPIXEL = 32
READOUT_SECONDS = 965.0
TOTAL_ROWS = 20
TOTAL_COLUMNS = 3_200
COLUMN_FAMILY_ALPHA = 0.01
CHARGED_THRESHOLD_E = 0.5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def load(path: Path) -> dict[str, np.ndarray]:
    with uproot.open(path) as root_file:
        tree = root_file["calPixTree"]
        return tree.arrays(["x", "y", "ePix", "mask", "RUNID"], library="np")


def mixture_fit(charge: np.ndarray) -> dict[str, float]:
    """Fine-binned two-Gaussian fit with fixed means 0e and 1e.

    Both components share sigma.  Fitting is restricted to [-0.5, 1.5], whose
    omitted Gaussian tail is below 0.05% at the observed sigma.  Parameters are
    optimized as log(sigma) and logit(single-electron fraction).
    """
    values = np.asarray(charge, dtype=np.float64)
    values = values[(values >= -0.5) & (values <= 1.5)]
    counts, edges = np.histogram(values, bins=2_000, range=(-0.5, 1.5))
    centers = 0.5 * (edges[:-1] + edges[1:])
    log_two_pi = math.log(2.0 * math.pi)

    def nll(theta: np.ndarray) -> float:
        sigma = math.exp(float(theta[0]))
        one_fraction = float(expit(theta[1]))
        log_zero = (
            math.log1p(-one_fraction)
            - math.log(sigma)
            - 0.5 * log_two_pi
            - 0.5 * (centers / sigma) ** 2
        )
        log_one = (
            math.log(one_fraction)
            - math.log(sigma)
            - 0.5 * log_two_pi
            - 0.5 * ((centers - 1.0) / sigma) ** 2
        )
        log_probability = logsumexp(np.vstack((log_zero, log_one)), axis=0)
        log_probability -= logsumexp(log_probability)
        return -float(np.sum(counts * log_probability))

    initial_fraction = max(float(np.mean(values > 0.5)), 1e-6)
    fit = minimize(
        nll,
        x0=np.array([math.log(0.14), math.log(initial_fraction / (1.0 - initial_fraction))]),
        method="Nelder-Mead",
        options={"xatol": 1e-10, "fatol": 1e-5, "maxiter": 1_000},
    )
    if not fit.success:
        raise RuntimeError(f"mixture fit failed: {fit.message}")
    sigma = math.exp(float(fit.x[0]))
    one_fraction = float(expit(fit.x[1]))

    # A transparent counting-scale uncertainty.  At sigma~0.14 the 0e/1e
    # overlap is tiny, so this agrees closely with the likelihood curvature.
    fraction_se = math.sqrt(one_fraction * (1.0 - one_fraction) / values.size)
    return {
        "n_fit": int(values.size),
        "sigma_e": sigma,
        "single_e_fraction": one_fraction,
        "single_e_fraction_se_counting": fraction_se,
        "nll": float(fit.fun),
    }


def select_release_active(arrays: dict[str, np.ndarray]) -> np.ndarray:
    y = arrays["y"]
    return (y >= ACTIVE_ROWS[0]) & (y <= ACTIVE_ROWS[1]) & (
        (arrays["mask"] & RELEASE_MASK_WORD) == 0
    )


def exposure_days(arrays: dict[str, np.ndarray], nominal_seconds: float) -> np.ndarray:
    return nominal_seconds / 86_400.0 + (READOUT_SECONDS / TOTAL_ROWS / 86_400.0) * (
        arrays["y"] + arrays["x"] / TOTAL_COLUMNS
    )


def weighted_line_fit(x: np.ndarray, y: np.ndarray, y_se: np.ndarray) -> dict[str, float]:
    design = np.column_stack((np.ones_like(x), x))
    weights = 1.0 / np.square(y_se)
    normal = design.T @ (weights[:, None] * design)
    covariance = np.linalg.inv(normal)
    beta = covariance @ (design.T @ (weights * y))
    residual = y - design @ beta
    chi2 = float(np.sum(np.square(residual / y_se)))
    return {
        "intercept_superpixel_per_image": float(beta[0]),
        "intercept_se": float(math.sqrt(covariance[0, 0])),
        "slope_superpixel_per_day": float(beta[1]),
        "slope_se": float(math.sqrt(covariance[1, 1])),
        "rate_pixel_per_day": float(beta[1] / PHYSICAL_PIXELS_PER_SUPERPIXEL),
        "rate_pixel_per_day_se": float(
            math.sqrt(covariance[1, 1]) / PHYSICAL_PIXELS_PER_SUPERPIXEL
        ),
        "chi2": chi2,
        "dof": int(len(x) - 2),
    }


def iterative_bad_columns(arrays: dict[str, np.ndarray]) -> dict[str, object]:
    """Flag columns with a Bonferroni-controlled binomial upper-tail test.

    This intentionally uses no release mask bits and no simulator truth.  A
    charged pixel is simply ePix>0.5 in the 3072x16 active area.  The common
    null probability is re-estimated after removing flagged columns, until the
    set is stable.  The family-wise type-I error is 1% over 3072 columns.
    """
    active = (
        (arrays["x"] >= 0)
        & (arrays["x"] < ACTIVE_COLUMNS)
        & (arrays["y"] >= ACTIVE_ROWS[0])
        & (arrays["y"] <= ACTIVE_ROWS[1])
    )
    x = arrays["x"][active]
    charged = arrays["ePix"][active] > CHARGED_THRESHOLD_E
    n = np.bincount(x, minlength=ACTIVE_COLUMNS).astype(np.int64)
    k = np.bincount(x, weights=charged.astype(np.int64), minlength=ACTIVE_COLUMNS).astype(
        np.int64
    )
    flagged = np.zeros(ACTIVE_COLUMNS, dtype=bool)
    history: list[dict[str, object]] = []
    threshold = COLUMN_FAMILY_ALPHA / ACTIVE_COLUMNS
    for iteration in range(20):
        keep = ~flagged
        probability = float(k[keep].sum() / n[keep].sum())
        pvalues = binom.sf(k - 1, n, probability)
        updated = pvalues < threshold
        history.append(
            {
                "iteration": iteration,
                "null_probability": probability,
                "columns": np.flatnonzero(updated).astype(int).tolist(),
            }
        )
        if np.array_equal(updated, flagged):
            break
        flagged = updated
    else:
        raise RuntimeError("bad-column selection did not stabilize")

    rate = np.divide(k, n, out=np.zeros_like(k, dtype=float), where=n > 0)
    rank = np.lexsort((np.arange(ACTIVE_COLUMNS), -rate))
    top = [
        {
            "rank": position + 1,
            "column": int(column),
            "charged": int(k[column]),
            "pixels": int(n[column]),
            "rate": float(rate[column]),
            "pvalue": float(pvalues[column]),
        }
        for position, column in enumerate(rank[:20])
    ]
    return {
        "criterion": (
            "iterated exact binomial upper-tail p-value < 0.01/3072; common null "
            "estimated from currently unflagged active columns; charged means ePix>0.5; "
            "release mask bits ignored"
        ),
        "flagged_columns": np.flatnonzero(flagged).astype(int).tolist(),
        "null_probability": float(probability),
        "bonferroni_threshold": threshold,
        "history": history,
        "top_20": top,
    }


def main() -> None:
    file_results = []
    arrays_by_file = []
    for filename, nominal_seconds in FILES:
        path = DATA / filename
        digest = sha256(path)
        if digest != EXPECTED_SHA256[filename]:
            raise RuntimeError(f"SHA-256 mismatch for {filename}: {digest}")
        arrays = load(path)
        arrays_by_file.append(arrays)
        selected = select_release_active(arrays)
        fit = mixture_fit(arrays["ePix"][selected])
        negative = arrays["ePix"][selected]
        negative = negative[(negative < 0.0) & (negative > -0.5)]
        negative_half_sigma = float(np.sqrt(np.mean(np.square(negative))))
        exp_days = exposure_days(arrays, nominal_seconds)[selected]
        columns = iterative_bad_columns(arrays)
        file_results.append(
            {
                "filename": filename,
                "sha256": digest,
                "nominal_exposure_seconds": nominal_seconds,
                "entries": int(len(arrays["x"])),
                "images": int(len(np.unique(arrays["RUNID"]))),
                "x_range": [int(arrays["x"].min()), int(arrays["x"].max())],
                "y_range": [int(arrays["y"].min()), int(arrays["y"].max())],
                "release_selected_pixels": int(selected.sum()),
                "release_selected_fraction": float(selected.mean()),
                "exposure_days_mean": float(exp_days.mean()),
                "exposure_days_std": float(exp_days.std()),
                "negative_half_peak_sigma_e": negative_half_sigma,
                "mixture_fit": fit,
                "independent_bad_columns": columns,
            }
        )

    x = np.array([item["exposure_days_mean"] for item in file_results])
    y = np.array([item["mixture_fit"]["single_e_fraction"] for item in file_results])
    y_se = np.array(
        [item["mixture_fit"]["single_e_fraction_se_counting"] for item in file_results]
    )
    rate_fit = weighted_line_fit(x, y, y_se)

    combined = {key: np.concatenate([arrays[key] for arrays in arrays_by_file]) for key in arrays_by_file[0]}
    combined_columns = iterative_bad_columns(combined)
    nested = all(
        set(file_results[index]["independent_bad_columns"]["flagged_columns"])
        <= set(file_results[index + 1]["independent_bad_columns"]["flagged_columns"])
        for index in range(len(file_results) - 1)
    )

    result = {
        "provenance": {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "script": str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
            "script_sha256": sha256(Path(__file__)),
            "repository_git_commit": git_commit(),
            "public_release_git_commit": "d30c4a7688b9ffbd4ee7a4ecc6a9957da691f98f",
            "imports_project_analysis_or_src": False,
        },
        "parameters": {
            "active_columns": ACTIVE_COLUMNS,
            "active_rows_inclusive": list(ACTIVE_ROWS),
            "release_mask_word": hex(RELEASE_MASK_WORD),
            "physical_pixels_per_superpixel": PHYSICAL_PIXELS_PER_SUPERPIXEL,
            "readout_seconds": READOUT_SECONDS,
            "column_family_alpha": COLUMN_FAMILY_ALPHA,
            "charged_threshold_e": CHARGED_THRESHOLD_E,
        },
        "files": file_results,
        "rate_fit": rate_fit,
        "bad_column_sets_nested_in_exposure_order": nested,
        "combined_bad_columns": combined_columns,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print("Independent SENSEI SNOLAB audit")
    for item in file_results:
        mix = item["mixture_fit"]
        cols = item["independent_bad_columns"]["flagged_columns"]
        print(
            f"{item['nominal_exposure_seconds']/3600:>4.0f} h: "
            f"sigma={mix['sigma_e']:.6f} e-, negative-half={item['negative_half_peak_sigma_e']:.6f} e-, "
            f"density={mix['single_e_fraction']:.8g} +/- {mix['single_e_fraction_se_counting']:.2g}, "
            f"bad columns={cols}"
        )
        print("       top 10:", [row["column"] for row in item["independent_bad_columns"]["top_20"][:10]])
    print("nested:", nested)
    print("combined bad columns:", combined_columns["flagged_columns"])
    print(
        "rate = "
        f"({rate_fit['rate_pixel_per_day']:.8g} +/- "
        f"{rate_fit['rate_pixel_per_day_se']:.2g}) e-/pixel/day; "
        f"chi2/dof={rate_fit['chi2']:.3f}/{rate_fit['dof']}"
    )
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
