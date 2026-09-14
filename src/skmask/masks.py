"""Six pixel masks, each in a FIXED form (constants given by hand) and an ADAPTIVE form.

A fixed mask takes its sizes and thresholds as arguments; in the experiments they are tuned on one
sensor and transplanted to another. An adaptive mask obtains them from the images it is applied
to, with a threshold set by a false-positive rate `alpha` (corrected for the number of tests), or
from physical properties of the sensor that are measurable (thickness, pixel size, diffusion).

All functions work on integer electron maps (see `events.to_electrons`); masks are boolean arrays
with True = discarded. Geometry follows `simulate.py`: charge-transfer trails point to +x and +y.

`exclude` arguments take either one boolean array for every image or, for functions that work on
a stack, a list with one array per image.

Order of calibration of the adaptive masks, each excluding what the previous ones removed:
CTI -> halo -> hot columns/pixels -> serial register -> low-energy clusters; muon tracks are
independent of the others.

CTI goes first because its calibration compares downstream with upstream of the same trigger in
the same row or column, which a hot column (uniform along the column) or a halo (symmetric) cannot
bias; hot columns cannot go first because vertical CTI trails of tracks are themselves a
column-wise excess (on the surface preset they produced 69 false hot columns out of 81 flagged).
"""
import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve
from scipy.stats import binom, poisson

from .events import clusters, single_electron_events
from .simulate import MIP_ELECTRONS_PER_UM

# CHOICE: charge that makes a pixel a trigger for halo and CTI masks. Low enough to include the
# point-like high-energy deposits of every preset, high enough that pixel noise never reaches it.
TRIGGER_E = 20

# CHOICE: a far-field reference smaller than this many pixels is treated as no reference at all.
MIN_FAR_FIELD_PIXELS = 10_000


def _exclude_for(exclude, i):
    if exclude is None:
        return None
    if isinstance(exclude, (list, tuple)):
        return exclude[i]
    return exclude


def _valid(shape, exclude):
    return np.ones(shape, dtype=bool) if exclude is None else ~exclude


def _disc(radius):
    r = int(np.ceil(radius))
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return (x ** 2 + y ** 2 <= radius ** 2).astype(float)


def dilate_disc(seeds, radius):
    """Pixels within `radius` of any seed (Euclidean), via convolution."""
    if radius <= 0 or not seeds.any():
        return seeds.copy()
    return fftconvolve(seeds.astype(float), _disc(radius), mode="same") > 0.5


def neighbour_counts(events, radius):
    """For every pixel, the number of events within `radius`, not counting the pixel itself."""
    return np.rint(fftconvolve(events.astype(float), _disc(radius), mode="same")) - events


# --------------------------------------------------------------------------------------------
# 1. Hot columns and hot pixels
# --------------------------------------------------------------------------------------------

def low_charge_occupancy(electrons, trigger_e=TRIGGER_E):
    """Charged pixels that belong to clusters below the trigger charge (defects, not tracks)."""
    labels, table = clusters(electrons)
    small = np.zeros(len(table) + 1, dtype=bool)
    small[1:] = table[:, 0] < trigger_e
    return (electrons > 0) & small[labels]


def fixed_hot_columns(electron_stack, factor):
    """Columns whose count of charged low-charge pixels exceeds `factor` times the median column."""
    counts = sum(low_charge_occupancy(e).sum(axis=0) for e in electron_stack)
    return np.nonzero(counts > factor * max(np.median(counts), 1.0))[0]


def fixed_hot_pixels(electron_stack, min_images):
    """Pixels charged (in a low-charge cluster) in at least `min_images` images of the stack."""
    counts = sum(low_charge_occupancy(e).astype(int) for e in electron_stack)
    return counts >= min_images


def adaptive_hot_columns(electron_stack, alpha=0.01, widths=(1, 2, 3), exclude=None, max_iter=50):
    """Columns (and groups of adjacent columns) with more charge than a Poisson rate allows.

    For each column c: n_c charged low-charge pixels over the stack, e_c valid pixels. The common
    rate is estimated from columns not yet flagged; a group of w adjacent unflagged columns is
    flagged when P(N >= n | rate * e) < alpha / (number of groups tested).

    Single columns take precedence: whenever anything is flagged the search restarts from w = 1
    with the rate re-estimated, and a wider group is tested only if none of its columns is flagged.
    Without this, every group containing one hot column is significant because of that column and
    its innocent neighbours are flagged with it.
    """
    n_c = np.zeros(electron_stack[0].shape[1])
    e_c = np.zeros_like(n_c)
    for i, electrons in enumerate(electron_stack):
        valid = _valid(electrons.shape, _exclude_for(exclude, i))
        n_c += (low_charge_occupancy(electrons) & valid).sum(axis=0)
        e_c += valid.sum(axis=0)
    n_tests = sum(len(n_c) - w + 1 for w in widths)
    flagged = np.zeros(len(n_c), dtype=bool)
    for _ in range(max_iter):
        changed = False
        for w in widths:
            good = ~flagged
            rate = n_c[good].sum() / max(e_c[good].sum(), 1.0)
            kernel = np.ones(w)
            n_w = np.convolve(n_c * good, kernel, mode="valid")
            e_w = np.convolve(e_c * good, kernel, mode="valid")
            clean_group = np.convolve(flagged, kernel, mode="valid") == 0
            p = poisson.sf(n_w - 1, rate * e_w)
            hits = np.nonzero(clean_group & (p < alpha / n_tests))[0]
            if len(hits):
                for start in hits:
                    flagged[start:start + w] = True
                changed = True
                break
        if not changed:
            break
    return np.nonzero(flagged)[0]


def adaptive_hot_pixels(electron_stack, alpha=0.01, exclude=None):
    """Pixels charged more often across the stack than the common per-pixel rate allows.

    With n images a pixel's count is at most n, so a single image can never flag anything; the
    rate is re-estimated once without the pixels flagged in the first pass.
    """
    counts = np.zeros(electron_stack[0].shape)
    exposures = np.zeros(electron_stack[0].shape)
    for i, electrons in enumerate(electron_stack):
        valid = _valid(electrons.shape, _exclude_for(exclude, i))
        counts += low_charge_occupancy(electrons) & valid
        exposures += valid
    flagged = np.zeros(counts.shape, dtype=bool)
    for _ in range(2):
        rate = counts[~flagged].sum() / max(exposures[~flagged].sum(), 1)
        flagged = poisson.sf(counts - 1, rate * exposures) < alpha / counts.size
    return flagged


def column_mask(shape, columns):
    mask = np.zeros(shape, dtype=bool)
    mask[:, columns] = True
    return mask


# --------------------------------------------------------------------------------------------
# 2. Charge-transfer inefficiency (bleeding trails)
# --------------------------------------------------------------------------------------------

def _distance_after_trigger(trigger, axis):
    """Distance to the nearest trigger pixel before each pixel along `axis` (inf if none)."""
    n = trigger.shape[axis]
    shape = [1, 1]
    shape[axis] = n
    index = np.arange(n, dtype=float).reshape(shape)
    last = np.maximum.accumulate(np.where(trigger, index, -np.inf), axis=axis)
    return index - last


def _distance_before_trigger(trigger, axis):
    return np.flip(_distance_after_trigger(np.flip(trigger, axis=axis), axis), axis=axis)


def cti_mask(electrons, length_h, length_v, trigger_e=TRIGGER_E):
    """Mask `length_h` pixels after each trigger along its row and `length_v` along its column."""
    trigger = electrons >= trigger_e
    mask = np.zeros(electrons.shape, dtype=bool)
    if length_h > 0:
        d = _distance_after_trigger(trigger, axis=1)
        mask |= (d >= 1) & (d <= length_h)
    if length_v > 0:
        d = _distance_after_trigger(trigger, axis=0)
        mask |= (d >= 1) & (d <= length_v)
    return mask


def adaptive_cti_lengths(electron_stack, alpha=0.01, block=5, max_distance=400,
                         trigger_e=TRIGGER_E, exclude=None):
    """Trail lengths from the downstream-over-upstream excess of single-electron events.

    For each distance block, n_down counts 1e events that far after a trigger (along the readout
    direction) and n_up that far before one. Halo, dark current and signal are symmetric, so
    without CTI n_down | n_down + n_up ~ Binomial(1/2). The length is the far edge of the last
    block with a one-sided p-value below alpha / (number of blocks tested).
    """
    n_blocks = max_distance // block
    lengths = {}
    profiles = {}
    for name, axis in (("h", 1), ("v", 0)):
        n_down = np.zeros(n_blocks)
        n_up = np.zeros(n_blocks)
        for i, electrons in enumerate(electron_stack):
            trigger = electrons >= trigger_e
            events = single_electron_events(electrons) & _valid(electrons.shape, _exclude_for(exclude, i))
            for distance, counts in ((_distance_after_trigger(trigger, axis), n_down),
                                     (_distance_before_trigger(trigger, axis), n_up)):
                use = events & (distance >= 1) & (distance <= n_blocks * block)
                k = ((distance[use] - 1) // block).astype(int)
                counts += np.bincount(k, minlength=n_blocks)[:n_blocks]
        total = n_down + n_up
        p = np.where(total > 0, binom.sf(n_down - 1, total, 0.5), 1.0)
        significant = np.nonzero(p < alpha / (2 * n_blocks))[0]
        lengths[name] = int((significant.max() + 1) * block) if len(significant) else 0
        profiles[name] = {"n_down": n_down, "n_up": n_up, "p": p}
    return lengths["h"], lengths["v"], profiles


# --------------------------------------------------------------------------------------------
# 3. Halo around high-energy deposits
# --------------------------------------------------------------------------------------------

def halo_mask(electrons, radius, trigger_e=TRIGGER_E, sampling=(1.0, 1.0)):
    """Mask every pixel within `radius` pixels of a trigger pixel (physical distance)."""
    trigger = electrons >= trigger_e
    if radius <= 0 or not trigger.any():
        return np.zeros(electrons.shape, dtype=bool)
    return ndimage.distance_transform_edt(~trigger, sampling=sampling) <= radius


def adaptive_halo_radius(electron_stack, alpha=0.01, annulus=5, max_radius=200,
                         trigger_e=TRIGGER_E, exclude=None, sampling=(1.0, 1.0)):
    """Radius beyond which the 1e density around triggers stops decreasing.

    Counts of 1e events and valid pixels are accumulated in annuli of width `annulus` around the
    nearest trigger, up to `max_radius`, plus everything farther ("far field"). Annulus k is
    compared with everything outside it pooled (annuli k+1.. and the far field) by a conditional
    binomial test; the radius is the outer edge of the last annulus with p < alpha / (annuli
    tested). An annulus is tested only if its outside pool has at least MIN_FAR_FIELD_PIXELS
    valid pixels.

    CHOICE: the outside pool, not a far field, is the reference. A far field does not exist where
    triggers are dense (at the surface, muon tracks are ~90 px apart and no pixel is 200 px from all
    of them): an earlier version fell back to the largest radius there and masked the whole image.
    Halo in the pool makes the test conservative (shorter radius).

    Returns (radius, info); info["calibrated"] is False when no annulus could be tested (then the
    radius is 0 and the caller must decide), info["at_limit"] is True when the last significant
    annulus is the last testable one, i.e. the halo may extend beyond what could be measured.
    """
    n_annuli = max_radius // annulus
    counts = np.zeros(n_annuli)
    pixels = np.zeros(n_annuli)
    far_counts = 0.0
    far_pixels = 0.0
    for i, electrons in enumerate(electron_stack):
        valid = _valid(electrons.shape, _exclude_for(exclude, i))
        events = single_electron_events(electrons)
        trigger = electrons >= trigger_e
        if trigger.any():
            distance = ndimage.distance_transform_edt(~trigger, sampling=sampling)
        else:
            distance = np.full(electrons.shape, np.inf)
        near = valid & (distance > 0) & (distance <= n_annuli * annulus)
        k = np.minimum(((distance[near] - 1e-9) // annulus).astype(int), n_annuli - 1)
        counts += np.bincount(k, weights=events[near], minlength=n_annuli)
        pixels += np.bincount(k, minlength=n_annuli)
        far = valid & (distance > n_annuli * annulus)
        far_counts += events[far].sum()
        far_pixels += far.sum()
    # outside pool of annulus k: annuli k+1 .. n_annuli-1 plus the far field
    outside_counts = np.concatenate([np.cumsum(counts[::-1])[::-1][1:], [0.0]]) + far_counts
    outside_pixels = np.concatenate([np.cumsum(pixels[::-1])[::-1][1:], [0.0]]) + far_pixels
    testable = (pixels > 0) & (outside_pixels >= MIN_FAR_FIELD_PIXELS)
    # Conditional test of two Poisson rates, carrying the reference's own uncertainty:
    # n_k | n_k + n_out ~ Binomial(n_k + n_out, P_k / (P_k + P_out)) under equal densities.
    # (Treating a reference rate as exact overstated significance: caught by the null test.)
    share = np.divide(pixels, pixels + outside_pixels, out=np.zeros_like(pixels), where=testable)
    p = np.where(testable, binom.sf(counts - 1, counts + outside_counts, share), np.nan)
    n_tested = int(testable.sum())
    info = {"counts": counts, "pixels": pixels, "far_pixels": far_pixels, "p": p,
            "outside_pixels": outside_pixels, "n_tested": n_tested}
    if n_tested == 0:
        info.update(calibrated=False, at_limit=False)
        return 0, info
    significant = np.nonzero(testable & (p < alpha / n_tested))[0]
    radius = int((significant.max() + 1) * annulus) if len(significant) else 0
    last_testable = int(np.nonzero(testable)[0].max())
    info.update(calibrated=True, at_limit=bool(len(significant) and significant.max() == last_testable))
    return radius, info


# --------------------------------------------------------------------------------------------
# 4. Serial-register hits
# --------------------------------------------------------------------------------------------

def _row_window_sums(values, width):
    padded = np.concatenate([np.zeros((values.shape[0], 1)), np.cumsum(values, axis=1)], axis=1)
    return padded[:, width:] - padded[:, :-width]


def fixed_serial_rows(electrons, window, min_charged, exclude=None):
    """Rows with at least `min_charged` charged low-charge pixels in some `window` consecutive
    pixels (SENSEI 2020, arXiv:2004.11378: four of five)."""
    occupied = low_charge_occupancy(electrons) & _valid(electrons.shape, exclude)
    return np.nonzero((_row_window_sums(occupied.astype(float), window) >= min_charged).any(axis=1))[0]


def adaptive_serial_rows(electrons, alpha=0.01, widths=(5, 10, 20, 50, 100), exclude=None):
    """Rows holding a window with more charged pixels than the image's occupancy allows.

    Occupancy rho is measured on the valid pixels of the image. For each width w, the densest
    window of each row is compared with Poisson(rho * valid pixels in the window); the p-value is
    corrected for the windows in the row, and a row is flagged below alpha / (rows * widths).
    """
    valid = _valid(electrons.shape, exclude)
    occupied = (low_charge_occupancy(electrons) & valid).astype(float)
    rho = occupied.sum() / max(valid.sum(), 1)
    flagged = np.zeros(electrons.shape[0], dtype=bool)
    threshold = alpha / (electrons.shape[0] * len(widths))
    for w in widths:
        if w > electrons.shape[1]:
            continue
        n = _row_window_sums(occupied, w)
        v = _row_window_sums(valid.astype(float), w)
        p_window = poisson.sf(n - 1, rho * v)
        p_row = np.minimum(1.0, p_window.min(axis=1) * n.shape[1])
        flagged |= p_row < threshold
    return np.nonzero(flagged)[0]


def row_mask(shape, rows):
    mask = np.zeros(shape, dtype=bool)
    mask[rows, :] = True
    return mask


# --------------------------------------------------------------------------------------------
# 5. Low-energy clusters
# --------------------------------------------------------------------------------------------

def fixed_low_energy_cluster_mask(electrons, min_charge, half_width, trigger_e=TRIGGER_E):
    """Square of `half_width` pixels around clusters with charge in [min_charge, trigger)
    (SENSEI 2020, arXiv:2004.11378: 4 pixels around clusters of at least 5 e for 1e)."""
    labels, table = clusters(electrons)
    keep = np.zeros(len(table) + 1, dtype=bool)
    keep[1:] = (table[:, 0] >= min_charge) & (table[:, 0] < trigger_e)
    seeds = keep[labels]
    return ndimage.maximum_filter(seeds, size=2 * half_width + 1) if seeds.any() else seeds


def fixed_clustered_events_mask(electrons, radius, min_neighbours, exclude=None):
    """Disc of `radius` around every 1e event with at least `min_neighbours` others within it.

    The fixed counterpart of the adaptive mask below: same geometry, a hand-set count instead of a
    count derived from the measured density.
    """
    events = single_electron_events(electrons) & _valid(electrons.shape, exclude)
    seeds = events & (neighbour_counts(events, radius) >= min_neighbours)
    return dilate_disc(seeds, radius)


def adaptive_low_energy_cluster_mask(electrons, alpha=0.01, radii=(5, 10, 20, 40, 80), exclude=None):
    """Discs around single-electron events with too many neighbours for a uniform density.

    For each event and radius r: n = other events within r, a = valid pixels within r; under a
    uniform density rho (measured on the valid area) n ~ Poisson(rho * a). Events with
    p < alpha / (events * radii) are flagged, and a disc of the radius at which each was flagged
    is masked around it.
    """
    valid = _valid(electrons.shape, exclude)
    events = single_electron_events(electrons) & valid
    n_events = int(events.sum())
    mask = np.zeros(electrons.shape, dtype=bool)
    if n_events < 2:
        return mask
    rho = n_events / valid.sum()
    threshold = alpha / (n_events * len(radii))
    for r in radii:
        area = np.rint(fftconvolve(valid.astype(float), _disc(r), mode="same"))
        p = poisson.sf(neighbour_counts(events, r) - 1, rho * area)
        flagged = events & (p < threshold)
        if flagged.any():
            mask |= dilate_disc(flagged, r)
    return mask


# --------------------------------------------------------------------------------------------
# 6. Muon tracks
# --------------------------------------------------------------------------------------------

def fixed_muon_mask(electrons, min_charge, min_length, dilation, labels_table=None):
    """Clusters above `min_charge` whose bounding-box diagonal exceeds `min_length` pixels."""
    labels, table = labels_table if labels_table is not None else clusters(electrons)
    if len(table) == 0:
        return np.zeros(electrons.shape, dtype=bool)
    diagonal = np.hypot(table[:, 3] - table[:, 2] + 1, table[:, 5] - table[:, 4] + 1)
    keep = np.zeros(len(table) + 1, dtype=bool)
    keep[1:] = (table[:, 0] >= min_charge) & (diagonal >= min_length)
    selected = keep[labels]
    return dilate_disc(selected, dilation) if dilation > 0 else selected


def track_features(electrons, labels, table, pixel_um):
    """Charge-weighted principal axes of each cluster: projected length and width in um."""
    features = np.zeros((len(table), 2))
    for i, box in enumerate(ndimage.find_objects(labels)):
        if box is None:          # label absent (e.g. zeroed out by the caller): leave features at 0
            continue
        sub = labels[box] == i + 1
        q = electrons[box] * sub
        ys, xs = np.nonzero(sub)
        w = q[ys, xs].astype(float)
        if len(w) < 2:
            continue
        coords = np.column_stack([xs, ys]).astype(float)
        mean = np.average(coords, axis=0, weights=w)
        cov = np.cov((coords - mean).T, aweights=w)
        eig = np.sort(np.linalg.eigvalsh(np.atleast_2d(cov)))[::-1]
        # a uniform segment of length L has variance L^2 / 12 along its axis
        features[i] = (np.sqrt(12.0 * max(eig[0], 0.0)) * pixel_um, np.sqrt(max(eig[-1], 0.0)) * pixel_um)
    return features


def adaptive_muon_mask(electrons, pixel_um, thickness_um, diffusion_sigma_back_um,
                       charge_tolerance=2.0, trigger_e=TRIGGER_E, labels_table=None):
    """Clusters whose charge matches a minimum-ionising particle crossing this sensor.

    A straight track crossing thickness t with projected length L deposits on average
    Q_mip = (dE/dx)_min / E_pair * sqrt(L^2 + t^2). A cluster is masked as a track when
    - it is a single track: charge within a factor `charge_tolerance` of Q_mip for its measured
      length, and width compatible with diffusion (at most three back-surface sigmas plus a pixel); or
    - it is a pile-up of tracks: charge at least Q_mip / tolerance for its length and at least the
      charge of two vertical tracks, 2 (dE/dx)_min t / E_pair, whatever its shape.
    The mask covers the cluster dilated by three back-surface sigmas. Only physical properties of
    the sensor enter.

    The pile-up clause was added after the surface preset showed crossing muons (about 200 tracks
    per image) merging into wide clusters with 1.2-2.9 times a single track's charge, which the
    single-track clause rejected: 24 % of muon pixels were missed.
    """
    labels, table = labels_table if labels_table is not None else clusters(electrons)
    mask = np.zeros(electrons.shape, dtype=bool)
    if len(table) == 0 or not (table[:, 0] >= trigger_e).any():
        return mask
    big = np.nonzero(table[:, 0] >= trigger_e)[0]
    features = np.zeros((len(table), 2))
    sub_labels = np.where(np.isin(labels, big + 1), labels, 0)
    features[big] = track_features(electrons, sub_labels, table, pixel_um)[big]
    length_um, width_um = features[:, 0], features[:, 1]
    q_mip = MIP_ELECTRONS_PER_UM * np.hypot(length_um, thickness_um)
    ratio = table[:, 0] / q_mip
    straight = width_um <= 3.0 * diffusion_sigma_back_um + pixel_um
    big = table[:, 0] >= trigger_e
    single_track = big & (ratio >= 1 / charge_tolerance) & (ratio <= charge_tolerance) & straight
    pile_up = big & (ratio >= 1 / charge_tolerance) & (table[:, 0] >= 2 * MIP_ELECTRONS_PER_UM * thickness_um)
    # A track leaving the image keeps only part of its depth, so Q_mip overestimates its charge.
    # For a cluster touching the border the only safe bound is the projected length itself, and it
    # must be longer than a diffusion blob (a blob of width sigma has L = sqrt(12) sigma), so that a
    # point-like deposit cut by the border is not taken for a track.
    ny, nx = electrons.shape
    at_border = (table[:, 2] == 0) | (table[:, 4] == 0) | (table[:, 3] == ny - 1) | (table[:, 5] == nx - 1)
    clipped = (big & at_border & straight
               & (length_um >= 2.0 * np.sqrt(12.0) * diffusion_sigma_back_um)
               & (table[:, 0] >= MIP_ELECTRONS_PER_UM * length_um / charge_tolerance))
    is_muon = np.zeros(len(table) + 1, dtype=bool)
    is_muon[1:] = single_track | pile_up | clipped
    selected = is_muon[labels]
    return dilate_disc(selected, 3.0 * diffusion_sigma_back_um / pixel_um)


def temporal_poisson_outliers(stack, alpha=0.01, n_iter=2):
    """Per-pixel outliers in a stack of images of the same scene (shape: images, ny, nx).

    Each pixel's rate is its mean over the stack, re-estimated without values already flagged;
    a value n is an outlier when P(N >= n | rate) < alpha / (number of values).
    """
    stack = np.asarray(stack, dtype=float)
    threshold = alpha / stack.size
    data = stack.copy()
    outlier = np.zeros(stack.shape, dtype=bool)
    for _ in range(n_iter):
        rate = np.maximum(np.nanmean(data, axis=0), 1e-9)
        outlier = poisson.sf(stack - 1, rate[None]) < threshold
        data = np.where(outlier, np.nan, stack)
    return outlier
