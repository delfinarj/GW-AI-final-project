"""Checks of thresholding, clustering and mask scoring against hand-computed answers."""
import numpy as np
import pytest
from scipy.stats import norm

from skmask.events import (
    clusters,
    one_electron_threshold,
    pixel_origin,
    score_mask,
    single_electron_events,
    to_electrons,
)
from skmask.simulate import Sensor, simulate


def test_threshold_hand_value():
    # 0.5 + 0.14^2 * ln(0.9999 / 1e-4) = 0.5 + 0.0196 * 9.2102 = 0.6805
    assert one_electron_threshold(0.14, 1e-4) == pytest.approx(0.6805, abs=1e-4)


def test_threshold_equalises_posteriors():
    sigma, mu = 0.17, 3e-3
    c = one_electron_threshold(sigma, mu)
    assert (1 - mu) * norm.pdf(c, 0, sigma) == pytest.approx(mu * norm.pdf(c, 1, sigma), rel=1e-9)


def test_threshold_minimises_misclassification():
    sigma, mu = 0.15, 1e-3
    def error(c):
        return (1 - mu) * norm.sf(c, 0, sigma) + mu * norm.cdf(c, 1, sigma)
    c = one_electron_threshold(sigma, mu)
    grid = np.linspace(0.3, 1.0, 7001)
    assert c == pytest.approx(grid[np.argmin(error(grid))], abs=2e-4)


def test_to_electrons_rounding():
    measured = np.array([-0.3, 0.5, 0.7, 1.2, 1.49, 1.51, 2.4, 2.6, 10.49])
    np.testing.assert_array_equal(to_electrons(measured, 0.14, 1e-4), [0, 0, 1, 1, 1, 2, 2, 3, 10])


def test_clusters_on_hand_made_image():
    e = np.zeros((6, 8), dtype=np.int64)
    e[0, 0] = 1                  # isolated
    e[2, 2] = 2; e[3, 3] = 1     # diagonal neighbours: one cluster with 8-connectivity
    e[5, 7] = 4
    labels, table = clusters(e)
    assert len(table) == 3
    assert sorted(table[:, 0]) == [1, 3, 4]
    _, table4 = clusters(e, connectivity=4)
    assert len(table4) == 4


def test_single_electron_events_require_isolation():
    e = np.zeros((5, 5), dtype=np.int64)
    e[1, 1] = 1
    e[3, 3] = 1; e[3, 4] = 1
    e[0, 4] = 2
    events = single_electron_events(e)
    assert events.sum() == 1 and events[1, 1]


def test_pixel_origin_picks_dominant_source():
    charge = {"dark": np.array([[1, 0], [0, 0]]), "muon": np.array([[3, 0], [0, 1]])}
    origin = pixel_origin(charge)
    assert origin[0, 0] == "muon" and origin[1, 1] == "muon" and origin[0, 1] == ""


def test_score_of_empty_and_full_masks():
    sensor = Sensor(nx=400, ny=400, noise_e=0.1, exposure_days=1.0, dark_e_per_pix_day=2e-3,
                    signal_e_per_pix=2e-3, n_hot_columns=5, hot_column_e_per_pix=0.05)
    image = simulate(sensor, 12)
    nothing = score_mask(np.zeros((400, 400), dtype=bool), image)
    everything = score_mask(np.ones((400, 400), dtype=bool), image)
    assert nothing["clean_pixels_kept"] == 1.0 and nothing["signal_efficiency"] == 1.0
    assert everything["clean_pixels_kept"] == 0.0 and everything["signal_efficiency"] == 0.0
    assert everything["per_origin"]["hot_column"]["removed_fraction"] == 1.0


def test_score_of_the_true_hot_columns():
    sensor = Sensor(nx=400, ny=400, noise_e=0.1, exposure_days=1.0, dark_e_per_pix_day=2e-3,
                    signal_e_per_pix=2e-3, n_hot_columns=5, hot_column_e_per_pix=0.05)
    image = simulate(sensor, 13)
    mask = np.zeros((400, 400), dtype=bool)
    mask[:, image.truth["hot_columns"]] = True
    s = score_mask(mask, image)
    assert s["per_origin"]["hot_column"]["removed_fraction"] == 1.0
    # 5 of 400 columns lost: clean pixels kept and signal efficiency close to 395/400
    assert s["clean_pixels_kept"] == pytest.approx(395 / 400, abs=0.002)
    assert s["signal_efficiency"] == pytest.approx(395 / 400, abs=0.05)
