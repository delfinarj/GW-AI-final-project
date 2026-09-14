"""The presets reproduce the public numbers they were calibrated against, to order of magnitude."""
import numpy as np

from skmask.presets import DEEP_UNDERGROUND, PRESETS, SHALLOW_UNDERGROUND, SURFACE_LAB
from skmask.simulate import simulate


def test_presets_differ_in_what_the_masks_must_adapt_to():
    names = ("halo_length_um", "cti_h_length_pix", "cti_v_length_pix", "hot_column_e_per_pix",
             "serial_hit_span_pix", "lowE_cluster_sigma_pix", "noise_e", "thickness_um")
    for name in names:
        values = {getattr(p, name) for p in PRESETS.values()}
        assert len(values) > 1, name


def test_active_masses_match_publications():
    # [S20] quotes the mass of one 15x15x675 um pixel as 3.537e-7 g.
    per_pixel_g = SHALLOW_UNDERGROUND.mass_kg * 1e3 / (SHALLOW_UNDERGROUND.nx * SHALLOW_UNDERGROUND.ny)
    assert abs(per_pixel_g / 3.537e-7 - 1) < 0.01
    # [S24]: 2.19 g per 6144x1024 CCD at 665 um; the preset is one quadrant of it
    assert abs(DEEP_UNDERGROUND.mass_kg * 1e3 * 4 / 2.19 - 1) < 0.02


def test_surface_single_electron_rate_is_order_1e_minus_2_per_pixel_per_day():
    # [O23]: at the surface the 1e rate is expected to be O(1e-2) e/pix/day, dominated by
    # low-energy radiation from high-energy events rather than dark current.
    image = simulate(SURFACE_LAB, 20260914)
    halo_rate = image.charge["halo"].mean() / SURFACE_LAB.exposure_days
    assert 3e-3 < halo_rate < 1e-1
