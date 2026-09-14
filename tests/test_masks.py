"""Geometry, null and injection tests for every mask.

Null test: with the targeted defect absent, an adaptive mask must (almost) never fire.
Injection test: with the defect present, it must remove most of what the defect produced.
"""
import numpy as np

from skmask import masks as M
from skmask.events import score_mask, to_electrons
from skmask.simulate import Sensor, diffusion_sigma_um, simulate

BASE = dict(noise_e=0.12, exposure_days=1.0, dark_e_per_pix_day=5e-4, signal_e_per_pix=5e-4)


def electrons_of(image):
    return to_electrons(image.measured, image.sensor.noise_e, image.total.mean())


def stack(sensor, n, seed):
    rng = np.random.default_rng(seed)
    images = [simulate(sensor, rng) for _ in range(n)]
    return images, [electrons_of(im) for im in images]


def sparse_bright(**extra):
    """About 40 point deposits of 2-7 keV per image on 800 x 400 pixels: triggers for halo and CTI
    tests, sparse enough that most of the image is far from all of them."""
    # 7e4 dru * 1.13e-4 kg * 1 day * 5 keV = 40 deposits; the fraction of pixels farther than 60 px
    # from all of them is about exp(-40 * pi * 60^2 / 320000) = 24 %
    return Sensor(nx=800, ny=400, highE_dru=7e4, highE_keV=(2.0, 7.0), **BASE, **extra)


# ---------------------------------------------------------------- geometry

def test_halo_mask_is_a_disc():
    e = np.zeros((101, 101), dtype=np.int64)
    e[50, 50] = 100
    mask = M.halo_mask(e, radius=10)
    assert abs(mask.sum() - np.pi * 100) < 2 * np.pi * 10   # area within one circumference
    assert mask[50, 60] and not mask[50, 61]


def test_cti_mask_points_downstream():
    e = np.zeros((20, 60), dtype=np.int64)
    e[5, 10] = 50
    mask = M.cti_mask(e, length_h=20, length_v=3)
    assert mask[5, 11:31].all() and not mask[5, 31] and not mask[5, :11].any()
    assert mask[6:9, 10].all() and not mask[9, 10] and not mask[:5, 10].any()


def test_fixed_serial_rows_four_of_five():
    e = np.zeros((10, 50), dtype=np.int64)
    e[3, [10, 11, 13, 14]] = 1
    e[7, [10, 12, 14, 16]] = 1
    assert list(M.fixed_serial_rows(e, window=5, min_charged=4)) == [3]


def test_temporal_outlier_finds_injected_spike():
    rng = np.random.default_rng(0)
    data = rng.poisson(0.01, size=(50, 40, 40)).astype(float)
    data[17, 5, 9] = 40
    outliers = M.temporal_poisson_outliers(data)
    assert outliers[17, 5, 9] and outliers.sum() == 1


# ---------------------------------------------------------------- hot columns

def test_hot_columns_null():
    _, electrons = stack(Sensor(nx=600, ny=300, **BASE), 3, 1)
    assert len(M.adaptive_hot_columns(electrons)) == 0


def test_hot_columns_injection():
    # Power: 3 images x 300 rows = 900 pixels per column at ~1e-3 -> ~1 count expected; the
    # Bonferroni cut 0.01 / 1797 = 5.6e-6 needs >= 8 counts. 2e-2 per pixel adds ~18, so every hot
    # column is detectable. (5e-3 added only ~4.5 and was undetectable: a test design error.)
    sensor = Sensor(nx=600, ny=300, n_hot_columns=6, hot_column_e_per_pix=2e-2, defect_seed=3, **BASE)
    images, electrons = stack(sensor, 3, 2)
    truth = set(images[0].truth["hot_columns"])
    assert all(set(im.truth["hot_columns"]) == truth for im in images)   # same sensor, same defects
    found = set(M.adaptive_hot_columns(electrons))
    assert truth <= found
    assert len(found - truth) <= 2      # neighbours of a hot column may be flagged in a group


# ---------------------------------------------------------------- CTI

def test_cti_null():
    _, electrons = stack(sparse_bright(), 3, 3)
    h, v, _ = M.adaptive_cti_lengths(electrons)
    assert h == 0 and v == 0


def test_cti_injection():
    sensor = sparse_bright(cti_h_prob=3e-3, cti_h_length_pix=15.0, cti_v_prob=3e-3, cti_v_length_pix=6.0)
    images, electrons = stack(sensor, 3, 4)
    h, v, _ = M.adaptive_cti_lengths(electrons)
    assert h > v > 0

    def removed(length_h, length_v):
        return np.mean([score_mask(M.cti_mask(e, length_h, length_v), im, e)["per_origin"]["cti"]["removed_fraction"]
                        for im, e in zip(images, electrons)])

    # Oracle lengths hold > 99.9 % of both geometric trails: 1 - (14/15)^100, 1 - (5/6)^40.
    oracle = removed(100, 40)
    assert oracle > 0.9
    # The adaptive lengths are limited by statistics, so the check is that the mask removes what its
    # own lengths predict: equal deferral probabilities make the two directions equally populated,
    # and a geometric trail of mean m lies within L with probability 1 - (1 - 1/m)^L.
    predicted = oracle * 0.5 * ((1 - (1 - 1 / 15) ** h) + (1 - (1 - 1 / 6) ** v))
    assert abs(removed(h, v) - predicted) < 0.1


# ---------------------------------------------------------------- halo

def test_halo_null():
    # One realisation cannot measure a false-positive rate; that is done over many seeds in
    # analysis/null_false_positive_rates.py (R3). Seed 5 of this configuration is a documented
    # ~4 sigma fluctuation of the uniform dark/signal field itself (it is present in the true charge
    # maps, independent of thresholding) and is not used as the unit-test realisation.
    _, electrons = stack(sparse_bright(), 3, 0)
    radius, info = M.adaptive_halo_radius(electrons, max_radius=60)
    assert info["calibrated"] and radius == 0


def test_halo_injection():
    sensor = sparse_bright(halo_yield_per_e=2e-2, halo_length_um=90.0)
    images, electrons = stack(sensor, 3, 6)
    radius, info = M.adaptive_halo_radius(electrons, max_radius=60)
    assert info["calibrated"] and radius > 0
    removed = [score_mask(M.halo_mask(e, radius), im, e)["per_origin"]["halo"]["removed_fraction"]
               for im, e in zip(images, electrons)]
    assert np.mean(removed) > 0.7


def test_halo_radius_weighs_exposure_when_triggers_are_dense():
    # Triggers every 60 px; within 25 px of each, the single-electron density is 30 % above the clean
    # density (1.3e-3 against 1e-3). With 3 images of 1200 x 1200 the excess is highly significant
    # out to ~25 px, but masking removes 0.3 rho of excess per pixel while losing rho of exposure, so
    # rho P / sqrt(rho P + B) decreases as soon as anything is masked: the exposure-aware radius must
    # be below the significance radius and must not mask most of the image.
    rng = np.random.default_rng(16)
    size, spacing, reach = 1200, 60, 25
    trigger = np.zeros((size, size), dtype=bool)
    trigger[spacing // 2::spacing, spacing // 2::spacing] = True
    near = M.halo_mask(trigger.astype(np.int64) * 100, reach) & ~trigger
    stack_e = []
    for _ in range(3):
        u = rng.random((size, size))
        e = ((u < 1e-3) | (near & (u < 1.3e-3))).astype(np.int64)
        e[trigger] = 100
        stack_e.append(e)
    radius, info = M.adaptive_halo_radius(stack_e, max_radius=60)
    assert info["calibrated"] and info["radius_significance"] >= 20
    assert radius < info["radius_significance"]
    assert M.halo_mask(stack_e[0], radius).mean() < 0.5


def test_halo_without_any_reference_is_not_calibrated():
    # 2 500 pixels in total: no annulus has an outside pool of 10 000 pixels to compare with
    e = np.zeros((50, 50), dtype=np.int64)
    e[::10, ::10] = 50
    radius, info = M.adaptive_halo_radius([e], max_radius=60)
    assert not info["calibrated"] and radius == 0


def test_halo_calibrates_without_a_far_field():
    # triggers every 60 px: no pixel is farther than max_radius = 60 from all of them, yet the
    # outside pools of the inner annuli are large, so the calibration runs and finds the halo
    sensor = Sensor(nx=900, ny=900, **BASE)
    rng = np.random.default_rng(15)
    stack_e = []
    for _ in range(2):
        image = simulate(sensor, rng)
        e = to_electrons(image.measured, sensor.noise_e, image.total.mean())
        trigger = np.zeros(e.shape, dtype=bool)
        trigger[30::60, 30::60] = True
        e[trigger] = 100
        # inject 1e halo events within 15 px of each trigger, never on a trigger itself (overwriting
        # one would remove it and open a hole farther than 60 px from every trigger)
        ty, tx = np.nonzero(trigger)
        for y, x in zip(ty, tx):
            dy, dx = rng.integers(-15, 16, size=(2, 6))
            hy, hx = np.clip(y + dy, 0, 899), np.clip(x + dx, 0, 899)
            keep = ~trigger[hy, hx]
            e[hy[keep], hx[keep]] = 1
        stack_e.append(e)
    radius, info = M.adaptive_halo_radius(stack_e, max_radius=60)
    assert info["far_pixels"] == 0
    assert info["calibrated"] and 10 <= radius <= 30


# ---------------------------------------------------------------- serial register

def test_serial_null():
    _, electrons = stack(Sensor(nx=1000, ny=400, **BASE), 3, 7)
    assert sum(len(M.adaptive_serial_rows(e)) for e in electrons) == 0


def test_serial_injection():
    sensor = Sensor(nx=1000, ny=400, serial_hits_per_image=10, serial_hit_electrons=(10, 30),
                    serial_hit_span_pix=(10, 60), **BASE)
    images, electrons = stack(sensor, 3, 8)
    hit_rows = found = 0
    for image, e in zip(images, electrons):
        rows = set(image.truth["serial_hits"][:, 0])
        hit_rows += len(rows)
        found += len(rows & set(M.adaptive_serial_rows(e)))
    assert found / hit_rows > 0.8


# ---------------------------------------------------------------- low-energy clusters

def test_low_energy_cluster_null():
    _, electrons = stack(Sensor(nx=800, ny=400, **BASE), 3, 9)
    assert max(M.adaptive_low_energy_cluster_mask(e).mean() for e in electrons) < 0.01


def test_low_energy_cluster_injection():
    sensor = Sensor(nx=800, ny=400, lowE_clusters_per_image=3, lowE_cluster_mean_events=25,
                    lowE_cluster_sigma_pix=8, **BASE)
    images, electrons = stack(sensor, 3, 10)
    removed = [score_mask(M.adaptive_low_energy_cluster_mask(e), im, e)["per_origin"]["lowE_cluster"]["removed_fraction"]
               for im, e in zip(images, electrons)]
    assert np.mean([r for r in removed if r is not None]) > 0.6


# ---------------------------------------------------------------- muons

def test_muon_mask_takes_tracks_cut_by_the_border_but_not_blobs():
    e = np.zeros((100, 100), dtype=np.int64)
    # straight segment along row 20 entering from the left edge: 20 px = 300 um, and a minimum-ionising
    # charge for its projected length only (103 e/um * 15 um per pixel)
    e[20, 0:20] = 1550
    # a diffusion blob of 2000 e cut by the right edge
    e[60:63, 97:100] = 222
    sigma_back = 16.0
    mask = M.adaptive_muon_mask(e, 15.0, 675.0, sigma_back)
    assert mask[20, 0:20].all()
    assert not mask[60:63, 97:100].any()


def test_muon_mask_separates_tracks_from_point_deposits():
    sensor = Sensor(nx=600, ny=600, thickness_um=675.0, muon_flux_per_cm2_day=1440.0,
                    highE_dru=3e6, highE_keV=(0.5, 10.0), noise_e=0.12, exposure_days=20 / 1440,
                    dark_e_per_pix_day=5e-4, signal_e_per_pix=5e-4)
    images, electrons = stack(sensor, 2, 11)
    sigma_back = float(diffusion_sigma_um(sensor, sensor.thickness_um))
    muon_kept = highE_masked = muon_total = highE_total = 0
    for image, e in zip(images, electrons):
        mask = M.adaptive_muon_mask(e, sensor.pixel_um, sensor.thickness_um, sigma_back)
        muon_pixels = image.charge["muon"] > 0
        highE_pixels = (image.charge["highE"] > 0) & (image.charge["muon"] == 0)
        muon_total += muon_pixels.sum(); muon_kept += (muon_pixels & ~mask).sum()
        highE_total += highE_pixels.sum(); highE_masked += (highE_pixels & mask).sum()
    assert muon_total > 0 and highE_total > 0
    assert muon_kept / muon_total < 0.05
    assert highE_masked / highE_total < 0.10
