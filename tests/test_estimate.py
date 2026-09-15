"""The noise and the single-electron density must be recoverable from the image alone."""
import numpy as np
import pytest

from skmask.estimate import electrons_from_image, estimate_noise_and_density
from skmask.events import one_electron_threshold
from skmask.presets import SURFACE_LAB
from skmask.simulate import simulate


def synthetic(noise, density, n=400_000, seed=0):
    rng = np.random.default_rng(seed)
    electrons = (rng.random(n) < density).astype(float)
    return electrons + rng.normal(0.0, noise, n)


@pytest.mark.parametrize("noise, density", [(0.14, 1e-3), (0.20, 5e-3), (0.12, 2e-2)])
def test_recovers_known_noise_and_density(noise, density):
    fit = estimate_noise_and_density(synthetic(noise, density))
    assert fit["noise_e"] == pytest.approx(noise, rel=0.05)
    assert fit["density"] == pytest.approx(density, rel=0.25)


def test_tracks_do_not_bias_the_fit():
    # a few pixels carrying thousands of electrons (a track) sit far outside the histogram range
    values = synthetic(0.15, 2e-3, seed=1)
    values[:200] = np.linspace(500, 5000, 200)
    fit = estimate_noise_and_density(values)
    assert fit["noise_e"] == pytest.approx(0.15, rel=0.05)
    assert fit["density"] == pytest.approx(2e-3, rel=0.25)


def test_surface_image_threshold_is_not_the_degenerate_half():
    # The true mean charge of a surface image is ~14 e/pix because of muon tracks; using it as the
    # density gave the degenerate threshold 0.5 e. The measured density must be far smaller.
    image = simulate(SURFACE_LAB, 4)
    electrons, fit = electrons_from_image(image.measured)
    true_single = float((image.total == 1).mean())
    assert fit["density"] == pytest.approx(true_single, rel=0.5)
    assert fit["noise_e"] == pytest.approx(SURFACE_LAB.noise_e, rel=0.1)
    threshold = one_electron_threshold(fit["noise_e"], fit["density"])
    assert 0.55 < threshold < 0.85
    # and far fewer 1e events than the degenerate threshold would produce
    degenerate = (image.measured >= 0.5) & (image.measured < 1.5)
    measured = (electrons == 1)
    assert measured.sum() < degenerate.sum()
