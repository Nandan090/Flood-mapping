"""Sentinel-1 part of the impact label (spec section 5.3, with 14.11, 14.12, 14.17).

Each relative orbit and pass direction is processed on its own. Sentinel-1
GRD values on Earth Engine are already in dB, so no further conversion is done.

The layover and shadow mask ports the terrain geometry of gee_s1_ard
(Mullissa et al. 2021, Remote Sensing 13(10), 1954), which follows Vollrath,
Mullissa and Reiche 2020 (Remote Sensing 12(11), 1867):
slope steepness in range alpha_r = atan(tan(slope) * cos(phi_i - phi_s)),
layover where alpha_r >= theta_i, shadow where alpha_r <= -(90 deg - theta_i).
phi_i is the look direction (the aspect of the incidence-angle band) and
phi_s the terrain aspect, both clockwise from north as ee.Terrain.aspect
returns them. The gee_s1_ard Python code negates the terrain aspect before
the subtraction; with that sign a slope facing the sensor gets a negative
alpha_r, so here both angles use the same sign and a slope facing the sensor
has positive alpha_r.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy import ndimage

MAD_SCALE = 1.4826


def terrain_slope_aspect(dem: np.ndarray, pixel_m: float) -> tuple[np.ndarray, np.ndarray]:
    """Slope (deg) and downslope aspect (deg clockwise from north) from 4 neighbours.

    Rows run north to south and columns west to east, as in a north-up raster.
    """
    z = np.asarray(dem, dtype=float)
    d_row, d_col = np.gradient(z, pixel_m)
    east = d_col
    north = -d_row
    slope = np.degrees(np.arctan(np.hypot(east, north)))
    aspect = np.degrees(np.arctan2(-east, -north)) % 360.0
    return slope, aspect


def layover_shadow(slope_deg: np.ndarray, aspect_deg: np.ndarray,
                   incidence_deg: np.ndarray, heading_deg: float,
                   buffer_m: float = 0.0, pixel_m: float = 10.0
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Layover and shadow per pixel for one orbit. NaN inputs are neither."""
    alpha_s = np.radians(slope_deg)
    theta_i = np.radians(incidence_deg)
    phi_r = np.radians(heading_deg - np.asarray(aspect_deg, dtype=float))
    alpha_r = np.arctan(np.tan(alpha_s) * np.cos(phi_r))
    with np.errstate(invalid="ignore"):
        layover = alpha_r >= theta_i
        shadow = alpha_r <= -(np.pi / 2 - theta_i)
    if buffer_m > 0:
        layover = _grow(layover, buffer_m, pixel_m)
        shadow = _grow(shadow, buffer_m, pixel_m)
    return layover, shadow


def _grow(mask: np.ndarray, distance_m: float, pixel_m: float) -> np.ndarray:
    if not mask.any():
        return mask
    return ndimage.distance_transform_edt(~mask, sampling=pixel_m) <= distance_m


def focal_median(image: np.ndarray, size: int) -> np.ndarray:
    """size x size median. A window containing nodata gives nodata."""
    data = np.asarray(image, dtype=float)
    bad = ~np.isfinite(data)
    filtered = ndimage.median_filter(np.where(bad, 0.0, data), size=size,
                                     mode="nearest")
    touched = ndimage.maximum_filter(bad.astype(np.uint8), size=size,
                                     mode="nearest") > 0
    filtered[touched] = np.nan
    return filtered


def temporal_median(scenes: list[np.ndarray], min_scenes: int) -> np.ndarray:
    """Per-pixel median over scenes; fewer than min_scenes valid values is nodata."""
    stack = np.stack(scenes)
    count = np.isfinite(stack).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out = np.nanmedian(stack, axis=0)
    out[count < min_scenes] = np.nan
    return out


def robust_z(delta: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, dict]:
    """z = (delta − median_ref) / (1.4826 × MAD_ref)."""
    values = delta[np.asarray(reference, dtype=bool) & np.isfinite(delta)]
    if values.size == 0:
        raise RuntimeError("The reference zone has no valid radar change values.")
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    if mad == 0:
        raise RuntimeError(
            "The reference zone has zero spread in radar change, so z-scores "
            "are undefined. No substitute spread was used.")
    z = (delta - median) / (MAD_SCALE * mad)
    return z, {"median_db": median, "mad_db": mad, "pixels": int(values.size)}


def agreement(fired: list[np.ndarray], valid: list[np.ndarray],
              fraction: float, strict: bool = False) -> np.ndarray:
    """Orbit vote (14.12): 1 impacted, 0 not, NaN where no orbit is valid.

    Default: at least `fraction` of the valid orbits fired (1 of 1, 1 of 2,
    2 of 3, 2 of 4 with 0.5). strict: every valid orbit fired (sar_strict).
    """
    valid_stack = np.stack([np.asarray(v, dtype=bool) for v in valid])
    fired_stack = np.stack([np.asarray(f, dtype=bool) for f in fired]) & valid_stack
    n_valid = valid_stack.sum(axis=0)
    n_fired = fired_stack.sum(axis=0)
    if strict:
        hit = n_fired == n_valid
    else:
        hit = n_fired >= fraction * n_valid
    out = np.where(hit & (n_fired > 0), 1.0, 0.0)
    out[n_valid == 0] = np.nan
    return out
