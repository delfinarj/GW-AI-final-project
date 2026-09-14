"""Checks of the simulator against answers known independently of the code."""
import numpy as np
import pytest

from skmask.simulate import (
    MIP_ELECTRONS_PER_UM,
    MUON_MARGIN_THICKNESSES,
    ORIGINS,
    Sensor,
    charge_transfer_inefficiency,
    diffuse,
    diffusion_sigma_um,
    halo_positions,
    simulate,
)


def test_mip_electrons_per_um_from_pdg_numbers():
    # 1.664 MeV cm^2/g * 2.329 g/cm^3 = 3.8755 MeV/cm = 387.55 eV/um; / 3.75 eV = 103.35 e/um
    assert MIP_ELECTRONS_PER_UM == pytest.approx(103.35, abs=0.01)


def test_diffusion_sigma_endpoints():
    sensor = Sensor(thickness_um=675.0)
    assert diffusion_sigma_um(sensor, 0.0) == 0.0
    expected = np.sqrt(-218.715 * np.log(1 - 1.015e-3 * 675.0))   # ~15.9 um, about one pixel
    assert diffusion_sigma_um(sensor, 675.0) == pytest.approx(expected)
    assert 15.0 < expected < 17.0


def test_diffuse_spread_matches_sigma():
    rng = np.random.default_rng(1)
    sensor = Sensor()
    n = 400_000
    x, y = diffuse(sensor, np.zeros(n), np.zeros(n), np.full(n, 500.0), rng)
    sigma = diffusion_sigma_um(sensor, 500.0)
    # the relative standard error of a sample std is 1/sqrt(2n) ~ 0.1%
    assert np.std(x) / sigma == pytest.approx(1.0, abs=0.006)
    assert np.std(y) / sigma == pytest.approx(1.0, abs=0.006)


def test_uniform_dark_rate_and_nothing_else():
    sensor = Sensor(nx=512, ny=512, exposure_days=1.0, dark_e_per_pix_day=0.01, spurious_e_per_pix=0.002)
    image = simulate(sensor, 2)
    expected = 0.012
    standard_error = np.sqrt(expected / image.charge["dark"].size)
    assert abs(image.charge["dark"].mean() - expected) < 5 * standard_error
    for origin in ORIGINS:
        if origin not in ("dark", "cti"):
            assert image.charge[origin].sum() == 0, origin


def test_total_is_sum_of_origins_and_noise_width():
    sensor = Sensor(nx=300, ny=300, noise_e=0.15, exposure_days=1.0, dark_e_per_pix_day=0.1)
    image = simulate(sensor, 3)
    assert set(image.charge) == set(ORIGINS)
    residual = image.measured - image.total
    assert np.std(residual) / 0.15 == pytest.approx(1.0, abs=0.01)


def test_same_seed_same_image():
    sensor = Sensor(nx=200, ny=100, dark_e_per_pix_day=1.0, muon_flux_per_cm2_day=1440.0,
                    highE_dru=1e5, halo_yield_per_e=1e-3, cti_h_prob=1e-3, serial_hits_per_image=2)
    a = simulate(sensor, 7)
    b = simulate(sensor, 7)
    np.testing.assert_array_equal(a.measured, b.measured)


def test_hot_columns_are_where_truth_says():
    sensor = Sensor(nx=200, ny=400, n_hot_columns=3, hot_column_e_per_pix=0.5)
    image = simulate(sensor, 4)
    columns_with_charge = np.nonzero(image.charge["hot_column"].sum(axis=0))[0]
    np.testing.assert_array_equal(columns_with_charge, image.truth["hot_columns"])
    per_pixel = image.charge["hot_column"][:, image.truth["hot_columns"]].mean()
    assert per_pixel == pytest.approx(0.5, abs=5 * np.sqrt(0.5 / (400 * 3)))


def test_defects_belong_to_the_sensor_not_the_image():
    sensor = Sensor(nx=300, ny=200, n_hot_columns=4, hot_column_e_per_pix=0.2, n_hot_pixels=10,
                    hot_pixel_e=2.0, defect_seed=42)
    rng = np.random.default_rng(1)
    a, b = simulate(sensor, rng), simulate(sensor, rng)
    np.testing.assert_array_equal(a.truth["hot_columns"], b.truth["hot_columns"])
    np.testing.assert_array_equal(a.truth["hot_pixels"], b.truth["hot_pixels"])
    other = simulate(sensor.with_(defect_seed=43), rng)
    assert not np.array_equal(a.truth["hot_columns"], other.truth["hot_columns"])


def test_high_energy_count_and_charge():
    sensor = Sensor(nx=400, ny=400, exposure_days=1.0, highE_dru=1e6, highE_keV=(1.0, 3.0))
    image = simulate(sensor, 14)
    expected = 1e6 * sensor.mass_kg * 1.0 * 2.0
    deposits = image.truth["highE"]
    assert abs(len(deposits) - expected) < 5 * np.sqrt(expected)
    # mean energy 2 keV -> 533 electrons; charge landing inside the sensor cannot exceed the total
    assert deposits[:, 3].mean() == pytest.approx(2000 / 3.75, rel=0.05)
    assert image.charge["highE"].sum() <= deposits[:, 3].sum()


def test_cti_trails_only_downstream_with_the_set_length():
    sensor = Sensor(nx=3000, ny=50, cti_h_prob=1e-2, cti_h_length_pix=20.0)
    total = np.zeros((50, 3000), dtype=np.int32)
    total[25, 100] = 1_000_000
    trail = charge_transfer_inefficiency(sensor, np.random.default_rng(5), total)
    assert trail[:, :101].sum() == 0            # nothing upstream, nothing in the source pixel
    assert trail.sum() == trail[25].sum()       # horizontal CTI stays in its row
    n = trail.sum()
    assert abs(n - 10_000) < 5 * np.sqrt(10_000)
    distance = np.arange(3000) - 100
    mean_distance = (trail[25] * distance).sum() / n
    # geometric with mean 20 has sd sqrt(1-p)/p ~ 19.5; standard error ~ 0.2
    assert mean_distance == pytest.approx(20.0, abs=1.0)


def test_vertical_cti_stays_in_its_column():
    sensor = Sensor(nx=50, ny=3000, cti_v_prob=1e-2, cti_v_length_pix=10.0)
    total = np.zeros((3000, 50), dtype=np.int32)
    total[100, 25] = 1_000_000
    trail = charge_transfer_inefficiency(sensor, np.random.default_rng(6), total)
    assert trail.sum() == trail[:, 25].sum()
    assert trail[:101, 25].sum() == 0


def test_muon_rate_angles_and_charge_per_length():
    sensor = Sensor(nx=500, ny=500, thickness_um=675.0, muon_flux_per_cm2_day=1440.0,
                    exposure_days=30.0 / 1440.0)
    image = simulate(sensor, 8)
    margin = MUON_MARGIN_THICKNESSES * 675.0
    width = 500 * 15.0 + 2 * margin
    expected = 1440.0 * width * width * 1e-8 * 30.0 / 1440.0     # 1 per cm^2 per minute for 30 min
    n = image.truth["muons_generated"]
    assert abs(n - expected) < 5 * np.sqrt(expected)
    # flux through a horizontal plane with cos^2 intensity: density 4u^3, mean cos = 4/5
    cos_theta = image.truth["muon_cos_theta_generated"]
    assert cos_theta.mean() == pytest.approx(0.8, abs=5 * np.sqrt(2 / 75 / len(cos_theta)))
    tracks = image.truth["muons"]
    assert len(tracks) > 0
    path = np.hypot(np.hypot(tracks[:, 2] - tracks[:, 0], tracks[:, 3] - tracks[:, 1]), 675.0)
    expected_e = MIP_ELECTRONS_PER_UM * path
    pulls = (tracks[:, 5] - expected_e) / np.sqrt(expected_e)
    assert abs(pulls.mean()) < 5 / np.sqrt(len(pulls))
    assert image.charge["muon"].sum() <= tracks[:, 5].sum()


def test_halo_projected_distance_matches_isotropic_exponential():
    # no escape: thick sensor with sources in the middle; no diffusion
    sensor = Sensor(thickness_um=1e7, diffusion_A_um2=0.0, halo_length_um=200.0)
    rng = np.random.default_rng(9)
    n = 400_000
    x, y, _ = halo_positions(sensor, rng, np.zeros(1), np.zeros(1), np.full(1, 5e6), n)
    assert len(x) == n
    # E[d sin(polar)] for d ~ Exp(L) and isotropic direction is L * pi / 4
    assert np.hypot(x, y).mean() == pytest.approx(200.0 * np.pi / 4, rel=0.005)


def test_halo_photons_escape_through_the_surfaces():
    sensor = Sensor(thickness_um=100.0, diffusion_A_um2=0.0, halo_length_um=1000.0)
    x, _, z = halo_positions(sensor, np.random.default_rng(10), np.zeros(1), np.zeros(1),
                             np.full(1, 50.0), 100_000)
    assert len(x) < 20_000
    assert z.min() >= 0.0 and z.max() <= 100.0


def test_serial_hits_are_single_rows():
    sensor = Sensor(nx=1000, ny=500, serial_hits_per_image=20)
    image = simulate(sensor, 11)
    hits = image.truth["serial_hits"]
    rows_with_charge = np.nonzero(image.charge["serial"].sum(axis=1))[0]
    assert set(rows_with_charge) <= set(hits[:, 0])
    assert image.charge["serial"].sum() == hits[:, 3].sum()
