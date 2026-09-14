"""Three simulated sensors that differ the way real Skipper-CCD deployments differ.

Every number is either MEASURED (taken from the public source cited next to it) or SCENARIO
(chosen here to span a plausible range; defect populations are properties of an individual
sensor and are not published in a transferable form). The point of the presets is that the
masks must work on all three without retuning, so SCENARIO values deliberately differ between
them.

Sources
-------
[S20]   SENSEI, arXiv:2004.11378 (MINOS, 2020): geometry, noise, rates, diffusion, 3370 dru.
[S24]   SENSEI, arXiv:2410.18716 (SNOLAB, 2024): geometry, noise, 1e rate, ~50 events/kg/day/keV.
[O23]   Oscura sensors, arXiv:2304.04401: format, thickness, noise, surface 1e rate and
        spurious charge, surface high-energy rate O(1e4-1e5) dru.
[SNO09] SNO, arXiv:0902.2776: muon flux at SNOLAB (3.31 +- 0.09) e-10 cm^-2 s^-1.
[SB14]  Fermilab-thesis-2014-08 (SciBath, MINOS hall, 100 m): muon flux 0.80 +- 0.04 m^-2 s^-1.
[PDG]   PDG cosmic-ray review (2019): ~1 cm^-2 min^-1 through a horizontal detector at sea level.
[D22]   Du et al., arXiv:2011.13939: photons from high-energy tracks explain much of the SENSEI
        MINOS 1e rate, ~450/g-day beyond a 60-pixel halo against ~900/g-day at ~5 pixels.
"""
from .simulate import Sensor

SECONDS_PER_DAY = 86400.0

# SCENARIO, same for all presets: halo electrons per deposited electron. Order of magnitude set so
# that at the MINOS muon flux [SB14] and 3370 dru [S20] the halo supplies a few hundred 1e per
# g-day, the size of the halo-dependent part of the SENSEI MINOS rate quoted in [D22], and at the
# surface a 1e rate of order 1e-2 e/pix/day, as expected in [O23].
HALO_YIELD_PER_E = 1e-4

# Injected uniform signal, electrons per pixel per image, identical everywhere so signal
# efficiency is comparable across presets.
SIGNAL_E_PER_PIX = 1e-4

DEEP_UNDERGROUND = Sensor(
    name="deep_underground",
    nx=3072, ny=512,                              # one quadrant, 6144x1024 CCD [S24]
    pixel_um=15.0, thickness_um=665.0,            # [S24]
    noise_e=0.14,                                 # [S24]
    exposure_days=20.0 / 24.0,                    # longest exposure in the cycle [S24]
    dark_e_per_pix_day=1.39e-5,                   # exposure-dependent 1e rate [S24]
    spurious_e_per_pix=6.94e-5 / 32.0,            # per superpixel of 32 pixels [S24]
    signal_e_per_pix=SIGNAL_E_PER_PIX,
    muon_flux_per_cm2_day=3.31e-10 * SECONDS_PER_DAY,   # [SNO09]
    highE_dru=50.0,                               # [S24]
    halo_yield_per_e=HALO_YIELD_PER_E,
    halo_length_um=100.0,                         # SCENARIO
    cti_h_prob=2e-4, cti_h_length_pix=30.0,       # SCENARIO
    cti_v_prob=1e-4, cti_v_length_pix=10.0,       # SCENARIO
    serial_hits_per_image=0.5,                    # SCENARIO
    serial_hit_electrons=(4, 20), serial_hit_span_pix=(5, 40),
    n_hot_columns=4, hot_column_e_per_pix=2e-3,   # SCENARIO
    n_hot_pixels=30, hot_pixel_e=1.0,             # SCENARIO
    lowE_clusters_per_image=0.5, lowE_cluster_mean_events=15.0, lowE_cluster_sigma_pix=8.0,
)

SHALLOW_UNDERGROUND = Sensor(
    name="shallow_underground",
    nx=3072, ny=443,                              # one quadrant [S20]
    pixel_um=15.0, thickness_um=675.0,            # [S20]
    noise_e=0.14,                                 # [S20]
    exposure_days=20.0 / 24.0,                    # [S20]
    dark_e_per_pix_day=1.594e-4,                  # 1e rate with extra shielding [S20]; CHOICE:
                                                  # used as dark rate although part of it is halo
    spurious_e_per_pix=1.664e-4,                  # [S20]
    signal_e_per_pix=SIGNAL_E_PER_PIX,
    muon_flux_per_cm2_day=0.80e-4 * SECONDS_PER_DAY,    # [SB14]
    highE_dru=3370.0,                             # [S20]
    halo_yield_per_e=HALO_YIELD_PER_E,
    halo_length_um=150.0,                         # SCENARIO; 1e rate flattens by ~30-60 px [S20]
    cti_h_prob=5e-4, cti_h_length_pix=60.0,       # SCENARIO
    cti_v_prob=3e-4, cti_v_length_pix=25.0,       # SCENARIO
    serial_hits_per_image=2.0,                    # SCENARIO
    serial_hit_electrons=(4, 40), serial_hit_span_pix=(5, 60),
    n_hot_columns=10, hot_column_e_per_pix=5e-3,  # SCENARIO
    n_hot_pixels=80, hot_pixel_e=2.0,             # SCENARIO
    lowE_clusters_per_image=2.0, lowE_cluster_mean_events=25.0, lowE_cluster_sigma_pix=12.0,
)

SURFACE_LAB = Sensor(
    name="surface_lab",
    nx=1278, ny=1058,                             # Oscura-format sensor [O23]
    pixel_um=15.0, thickness_um=725.0,            # [O23]
    noise_e=0.17,                                 # between 0.15 and 0.19 e [O23]
    exposure_days=1.0 / 24.0,                     # SCENARIO: short exposures at the surface
    dark_e_per_pix_day=0.03,                      # lowest dark current measured, 140 K [O23]
    spurious_e_per_pix=8.4e-4,                    # [O23]
    signal_e_per_pix=SIGNAL_E_PER_PIX,
    muon_flux_per_cm2_day=1.0 * 1440.0,           # [PDG]
    highE_dru=3e4,                                # within O(1e4-1e5) [O23]
    halo_yield_per_e=HALO_YIELD_PER_E,
    halo_length_um=250.0,                         # SCENARIO
    cti_h_prob=1e-3, cti_h_length_pix=100.0,      # SCENARIO
    cti_v_prob=1e-3, cti_v_length_pix=40.0,       # SCENARIO
    serial_hits_per_image=5.0,                    # SCENARIO
    serial_hit_electrons=(10, 80), serial_hit_span_pix=(10, 100),
    n_hot_columns=25, hot_column_e_per_pix=2e-2,  # SCENARIO
    n_hot_pixels=200, hot_pixel_e=3.0,            # SCENARIO
    lowE_clusters_per_image=1.0, lowE_cluster_mean_events=30.0, lowE_cluster_sigma_pix=20.0,
)

PRESETS = {s.name: s for s in (DEEP_UNDERGROUND, SHALLOW_UNDERGROUND, SURFACE_LAB)}
