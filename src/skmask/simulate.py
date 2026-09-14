"""Skipper-CCD image simulator that keeps the origin of every electron.

Each source of charge fills its own integer map, so a mask applied to the summed image can be
scored against each source separately.

Conventions
-----------
- Images are arrays of shape (ny, nx), indexed [row, column] = [y, x].
- The amplifier sits at pixel (x=0, y=0) and charge is read in order of increasing y and x, so
  charge deferred by transfer inefficiency trails towards +x and +y.
- Depth z is measured from the collection (front) surface: a deposit at z = 0 does not diffuse,
  one at the back surface (z = thickness) diffuses the most.

Modelling choices with a defensible alternative are marked "CHOICE" and listed in PROVENANCE.md.
"""
from dataclasses import dataclass, field, replace

import numpy as np

# Effective energy per electron-hole pair in silicon, arXiv:2004.10709.
PAIR_ENERGY_EV = 3.75
# Silicon density and minimum-ionising stopping power, PDG atomic and nuclear properties.
SILICON_DENSITY_G_CM3 = 2.329
MIP_DEDX_MEV_CM2_G = 1.664
# MeV/cm -> eV/um is a factor 1e6 / 1e4 = 1e2.
MIP_ELECTRONS_PER_UM = MIP_DEDX_MEV_CM2_G * SILICON_DENSITY_G_CM3 * 1e2 / PAIR_ENERGY_EV

# CHOICE: muon entry points are drawn over the active area enlarged by this many sensor
# thicknesses on each side, and projected track lengths are capped at the same distance.
# For a cos^2 zenith distribution the capped fraction is P(tan theta > 20) ~ 6e-6.
MUON_MARGIN_THICKNESSES = 20.0

ORIGINS = (
    "dark",          # dark current and spurious charge, uniform
    "signal",        # injected uniform single electrons, the "physics" that must survive
    "hot_column",
    "hot_pixel",
    "muon",
    "highE",         # diffusion-limited high-energy deposits (radioactivity, Compton)
    "halo",          # single electrons from photons emitted by muon and high-energy deposits
    "serial",        # serial-register hits during readout
    "lowE_cluster",  # spatially clustered low-energy events of unknown origin
    "cti",           # charge deferred by transfer inefficiency
)


@dataclass(frozen=True)
class Sensor:
    """Everything that makes one simulated Skipper-CCD different from another.

    Rates are per pixel unless stated; `exposure_days` converts per-day rates to per-image.
    """

    name: str = "generic"
    nx: int = 1024
    ny: int = 512
    pixel_um: float = 15.0
    thickness_um: float = 675.0
    # Diffusion sigma_xy(z) = sqrt(-A ln|1 - b z|), SENSEI arXiv:2004.11378.
    diffusion_A_um2: float = 218.715
    diffusion_b_per_um: float = 1.015e-3
    noise_e: float = 0.15
    exposure_days: float = 1.0 / 24

    dark_e_per_pix_day: float = 0.0
    spurious_e_per_pix: float = 0.0
    signal_e_per_pix: float = 0.0

    muon_flux_per_cm2_day: float = 0.0
    highE_dru: float = 0.0                       # events / (kg day keV)
    highE_keV: tuple = (0.5, 10.0)

    halo_yield_per_e: float = 0.0                # halo electrons per deposited electron
    halo_length_um: float = 150.0                # photon absorption length

    cti_h_prob: float = 0.0                      # deferred fraction along the row (+x)
    cti_h_length_pix: float = 1.0                # mean trail length
    cti_v_prob: float = 0.0                      # deferred fraction along the column (+y)
    cti_v_length_pix: float = 1.0

    serial_hits_per_image: float = 0.0
    serial_hit_electrons: tuple = (4, 40)
    serial_hit_span_pix: tuple = (5, 60)

    # Defect positions belong to the sensor, not to an image: they are drawn from `defect_seed`
    # so that every image of the same sensor has its hot columns and pixels in the same place.
    defect_seed: int = 0
    n_hot_columns: int = 0
    hot_column_e_per_pix: float = 0.0            # extra electrons per pixel per image
    n_hot_pixels: int = 0
    hot_pixel_e: float = 0.0                     # extra electrons per image

    lowE_clusters_per_image: float = 0.0
    lowE_cluster_mean_events: float = 20.0
    lowE_cluster_sigma_pix: float = 10.0

    @property
    def area_cm2(self):
        return self.nx * self.ny * (self.pixel_um * 1e-4) ** 2

    @property
    def mass_kg(self):
        return self.area_cm2 * self.thickness_um * 1e-4 * SILICON_DENSITY_G_CM3 / 1e3

    def with_(self, **changes):
        return replace(self, **changes)


@dataclass
class SimulatedImage:
    sensor: Sensor
    charge: dict              # origin -> int32 array (ny, nx), electrons
    measured: np.ndarray      # total charge plus Gaussian readout noise, electrons
    truth: dict = field(default_factory=dict)

    @property
    def total(self):
        return sum(self.charge.values())


def diffusion_sigma_um(sensor, z_um):
    """Lateral diffusion sigma for charge created at depth z (from the collection surface)."""
    z = np.clip(np.asarray(z_um, dtype=float), 0.0, sensor.thickness_um)
    return np.sqrt(np.abs(-sensor.diffusion_A_um2 * np.log(np.abs(1.0 - sensor.diffusion_b_per_um * z))))


def diffuse(sensor, x_um, y_um, z_um, rng):
    """Move each electron laterally by a Gaussian step of width sigma(z)."""
    sigma = diffusion_sigma_um(sensor, z_um)
    x = np.asarray(x_um, dtype=float)
    y = np.asarray(y_um, dtype=float)
    return x + rng.normal(size=x.shape) * sigma, y + rng.normal(size=y.shape) * sigma


def _deposit(image, sensor, x_um, y_um):
    """Add one electron per (x, y) position to `image`, dropping positions off the sensor."""
    ix = np.floor(np.asarray(x_um) / sensor.pixel_um).astype(np.int64)
    iy = np.floor(np.asarray(y_um) / sensor.pixel_um).astype(np.int64)
    inside = (ix >= 0) & (ix < sensor.nx) & (iy >= 0) & (iy < sensor.ny)
    flat = iy[inside] * sensor.nx + ix[inside]
    image += np.bincount(flat, minlength=sensor.nx * sensor.ny).reshape(image.shape).astype(image.dtype)


def halo_positions(sensor, rng, x_um, y_um, z_um, n_photons):
    """Absorption points of photons emitted isotropically from randomly chosen deposit points.

    CHOICE: an exponential absorption distance with a single length, and photons that leave the
    silicon through either surface are lost. Real Cherenkov and recombination photons have a
    wavelength-dependent absorption length (arXiv:2011.13939); one length is the simplest model
    that gives a halo whose extent differs between sensors.
    """
    idx = rng.integers(0, len(x_um), n_photons)
    distance = rng.exponential(sensor.halo_length_um, n_photons)
    cos_polar = rng.uniform(-1.0, 1.0, n_photons)
    azimuth = rng.uniform(0.0, 2.0 * np.pi, n_photons)
    sin_polar = np.sqrt(1.0 - cos_polar ** 2)
    x = np.asarray(x_um)[idx] + distance * sin_polar * np.cos(azimuth)
    y = np.asarray(y_um)[idx] + distance * sin_polar * np.sin(azimuth)
    z = np.asarray(z_um)[idx] + distance * cos_polar
    kept = (z >= 0.0) & (z <= sensor.thickness_um)
    return x[kept], y[kept], z[kept]


def _emit_halo(sensor, rng, image, x_um, y_um, z_um, n_deposited):
    if sensor.halo_yield_per_e <= 0 or len(x_um) == 0:
        return
    n_photons = rng.poisson(sensor.halo_yield_per_e * n_deposited)
    x, y, z = halo_positions(sensor, rng, x_um, y_um, z_um, n_photons)
    xd, yd = diffuse(sensor, x, y, z, rng)
    _deposit(image, sensor, xd, yd)


def _muons(sensor, rng, charge, truth):
    """Straight minimum-ionising tracks through the full thickness.

    Zenith angles follow a cos^2 intensity (PDG cosmic-ray review); the flux through a horizontal
    plane adds one more cos, so cos(theta) has density 4 u^3 on [0, 1], sampled as U^(1/4).

    CHOICE: the charge per unit length is Poisson around the mean minimum-ionising value; the
    Landau tail from delta rays is not modelled.
    """
    t = sensor.thickness_um
    margin = MUON_MARGIN_THICKNESSES * t
    width = sensor.nx * sensor.pixel_um + 2 * margin
    height = sensor.ny * sensor.pixel_um + 2 * margin
    n = rng.poisson(sensor.muon_flux_per_cm2_day * width * height * 1e-8 * sensor.exposure_days)

    cos_theta = rng.uniform(size=n) ** 0.25
    tan_theta = np.sqrt(1.0 - cos_theta ** 2) / cos_theta
    projected = np.minimum(t * tan_theta, margin)
    azimuth = rng.uniform(0.0, 2.0 * np.pi, n)
    x0 = rng.uniform(-margin, width - margin, n)
    y0 = rng.uniform(-margin, height - margin, n)
    x1 = x0 + projected * np.cos(azimuth)
    y1 = y0 + projected * np.sin(azimuth)

    xmax = sensor.nx * sensor.pixel_um
    ymax = sensor.ny * sensor.pixel_um
    crosses = ((np.maximum(x0, x1) >= 0) & (np.minimum(x0, x1) < xmax)
               & (np.maximum(y0, y1) >= 0) & (np.minimum(y0, y1) < ymax))

    records = []
    for i in np.nonzero(crosses)[0]:
        n_e = rng.poisson(MIP_ELECTRONS_PER_UM * np.hypot(projected[i], t))
        s = rng.uniform(size=n_e)
        x = x0[i] + s * (x1[i] - x0[i])
        y = y0[i] + s * (y1[i] - y0[i])
        z = t * (1.0 - s)   # enters at the back surface, leaves at the front
        xd, yd = diffuse(sensor, x, y, z, rng)
        _deposit(charge["muon"], sensor, xd, yd)
        _emit_halo(sensor, rng, charge["halo"], x, y, z, n_e)
        records.append((x0[i], y0[i], x1[i], y1[i], cos_theta[i], n_e))

    truth["muons_generated"] = n
    truth["muon_cos_theta_generated"] = cos_theta
    truth["muons"] = np.array(records, dtype=float).reshape(-1, 6)   # x0 y0 x1 y1 cos n_e


def _high_energy(sensor, rng, charge, truth):
    """Point-like deposits with a flat energy spectrum at a uniform depth.

    CHOICE: all high-energy background is diffusion-limited (no extended electron tracks), which
    keeps the halo and CTI tests clean; muons provide the extended tracks.
    """
    emin, emax = sensor.highE_keV
    n = rng.poisson(sensor.highE_dru * sensor.mass_kg * sensor.exposure_days * (emax - emin))
    energy_keV = rng.uniform(emin, emax, n)
    n_e = np.rint(energy_keV * 1e3 / PAIR_ENERGY_EV).astype(np.int64)
    xc = rng.uniform(0, sensor.nx * sensor.pixel_um, n)
    yc = rng.uniform(0, sensor.ny * sensor.pixel_um, n)
    zc = rng.uniform(0, sensor.thickness_um, n)
    # all deposits at once; drawing the halo photons of all of them together, each photon from a
    # uniformly chosen electron, has the same distribution as one Poisson draw per deposit
    x = np.repeat(xc, n_e)
    y = np.repeat(yc, n_e)
    z = np.repeat(zc, n_e)
    xd, yd = diffuse(sensor, x, y, z, rng)
    _deposit(charge["highE"], sensor, xd, yd)
    _emit_halo(sensor, rng, charge["halo"], x, y, z, len(x))
    truth["highE"] = np.column_stack([xc, yc, zc, n_e]).reshape(-1, 4)


def _serial_hits(sensor, rng, charge, truth):
    """Charge scattered along a single row, as left by a hit in the serial register."""
    n = rng.poisson(sensor.serial_hits_per_image)
    records = []
    for _ in range(n):
        span = int(rng.integers(sensor.serial_hit_span_pix[0], sensor.serial_hit_span_pix[1] + 1))
        span = min(span, sensor.nx)
        row = int(rng.integers(0, sensor.ny))
        x0 = int(rng.integers(0, sensor.nx - span + 1))
        n_e = int(rng.integers(sensor.serial_hit_electrons[0], sensor.serial_hit_electrons[1] + 1))
        columns = x0 + rng.integers(0, span, n_e)
        np.add.at(charge["serial"], (np.full(n_e, row), columns), 1)
        records.append((row, x0, span, n_e))
    truth["serial_hits"] = np.array(records, dtype=np.int64).reshape(-1, 4)


def _low_energy_clusters(sensor, rng, charge, truth):
    n = rng.poisson(sensor.lowE_clusters_per_image)
    xc = rng.uniform(0, sensor.nx * sensor.pixel_um, n)
    yc = rng.uniform(0, sensor.ny * sensor.pixel_um, n)
    sigma_um = sensor.lowE_cluster_sigma_pix * sensor.pixel_um
    for i in range(n):
        m = rng.poisson(sensor.lowE_cluster_mean_events)
        _deposit(charge["lowE_cluster"], sensor,
                 xc[i] + rng.normal(size=m) * sigma_um, yc[i] + rng.normal(size=m) * sigma_um)
    truth["lowE_clusters"] = np.column_stack([xc, yc]).reshape(-1, 2)


def charge_transfer_inefficiency(sensor, rng, total):
    """Deferred charge trailing each pixel along +x (row transfer) and +y (column transfer).

    Each electron in a pixel is deferred with probability p and reappears d pixels downstream,
    d geometric with mean equal to the trail length.

    CHOICE: deferred charge is added without being removed from its source pixel; for p << 1 the
    bias on the source pixel is negligible and it keeps every other origin map exact.
    """
    out = np.zeros_like(total)
    iy, ix = np.nonzero(total)
    q = total[iy, ix]
    for prob, length, axis in ((sensor.cti_h_prob, sensor.cti_h_length_pix, "x"),
                               (sensor.cti_v_prob, sensor.cti_v_length_pix, "y")):
        if prob <= 0:
            continue
        n_deferred = rng.poisson(q * prob)
        has = n_deferred > 0
        sx = np.repeat(ix[has], n_deferred[has])
        sy = np.repeat(iy[has], n_deferred[has])
        step = rng.geometric(1.0 / max(length, 1.0), size=sx.size)
        if axis == "x":
            sx = sx + step
        else:
            sy = sy + step
        inside = (sx < sensor.nx) & (sy < sensor.ny)
        flat = sy[inside] * sensor.nx + sx[inside]
        out += np.bincount(flat, minlength=total.size).reshape(total.shape).astype(out.dtype)
    return out


def simulate(sensor, rng=None):
    """Simulate one image. `rng` is a numpy Generator or a seed."""
    rng = np.random.default_rng(rng)
    shape = (sensor.ny, sensor.nx)
    charge = {origin: np.zeros(shape, dtype=np.int32) for origin in ORIGINS}
    truth = {}

    charge["dark"] += rng.poisson(sensor.dark_e_per_pix_day * sensor.exposure_days
                                  + sensor.spurious_e_per_pix, shape).astype(np.int32)
    if sensor.signal_e_per_pix > 0:
        charge["signal"] += rng.poisson(sensor.signal_e_per_pix, shape).astype(np.int32)

    defects = np.random.default_rng(sensor.defect_seed)
    hot_columns = defects.choice(sensor.nx, size=sensor.n_hot_columns, replace=False)
    hot_flat = defects.choice(sensor.nx * sensor.ny, size=sensor.n_hot_pixels, replace=False)

    if len(hot_columns):
        charge["hot_column"][:, hot_columns] += rng.poisson(
            sensor.hot_column_e_per_pix, (sensor.ny, len(hot_columns))).astype(np.int32)
    truth["hot_columns"] = np.sort(hot_columns)

    if len(hot_flat):
        hy, hx = np.unravel_index(hot_flat, shape)
        charge["hot_pixel"][hy, hx] += rng.poisson(sensor.hot_pixel_e, len(hot_flat)).astype(np.int32)
    truth["hot_pixels"] = np.column_stack(np.unravel_index(hot_flat, shape)).reshape(-1, 2)  # y, x

    _muons(sensor, rng, charge, truth)
    _high_energy(sensor, rng, charge, truth)
    _serial_hits(sensor, rng, charge, truth)
    _low_energy_clusters(sensor, rng, charge, truth)

    charge["cti"] = charge_transfer_inefficiency(sensor, rng, sum(charge.values()))

    total = sum(charge.values())
    measured = total + rng.normal(0.0, sensor.noise_e, shape)
    return SimulatedImage(sensor=sensor, charge=charge, measured=measured, truth=truth)
