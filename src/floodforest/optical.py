"""Sentinel-2 part of the impact label (spec section 5.2, with 14.18).

Reflectance is on a 0–1 scale. Labels are float arrays: 1 impacted,
0 not impacted, NaN nodata.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

_EIGHT = np.ones((3, 3), dtype=int)


def ndvi(nir: np.ndarray, red: np.ndarray) -> np.ndarray:
    """(B8 − B4) / (B8 + B4)."""
    return _ratio(np.asarray(nir, dtype=float) - red,
                  np.asarray(nir, dtype=float) + red)


def mndwi(green: np.ndarray, swir: np.ndarray) -> np.ndarray:
    """(B3 − B11) / (B3 + B11)."""
    return _ratio(np.asarray(green, dtype=float) - swir,
                  np.asarray(green, dtype=float) + swir)


def bsi(blue: np.ndarray, red: np.ndarray, nir: np.ndarray,
        swir: np.ndarray) -> np.ndarray:
    """((B11 + B4) − (B8 + B2)) / ((B11 + B4) + (B8 + B2))."""
    bare = np.asarray(swir, dtype=float) + red
    green_cover = np.asarray(nir, dtype=float) + blue
    return _ratio(bare - green_cover, bare + green_cover)


def bright(blue: np.ndarray, green: np.ndarray, red: np.ndarray) -> np.ndarray:
    """(B2 + B3 + B4) / 3."""
    return (np.asarray(blue, dtype=float) + green + red) / 3.0


def spectral_channel(ndvi_pre: np.ndarray, mndwi_pre: np.ndarray,
                     ndvi_max: float, mndwi_min: float) -> np.ndarray:
    """Water or bare gravel before the event: MNDWI_pre > min or NDVI_pre < max."""
    ok = np.isfinite(ndvi_pre) & np.isfinite(mndwi_pre)
    with np.errstate(invalid="ignore"):
        return ok & ((mndwi_pre > mndwi_min) | (ndvi_pre < ndvi_max))


def keep_components_touching(mask: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    """Pixels of `mask` whose 8-connected component contains an anchor pixel."""
    mask = np.asarray(mask, dtype=bool)
    labels, _ = ndimage.label(mask, structure=_EIGHT)
    ids = np.unique(labels[mask & np.asarray(anchor, dtype=bool)])
    ids = ids[ids != 0]
    if ids.size == 0:
        return np.zeros(mask.shape, dtype=bool)
    return np.isin(labels, ids)


def t_bsi_from_channel(bsi_post: np.ndarray, channel: np.ndarray,
                       percentile: float) -> float:
    """t_bsi = the given percentile of BSI_post inside the pre-event channel."""
    values = bsi_post[np.asarray(channel, dtype=bool) & np.isfinite(bsi_post)]
    if values.size == 0:
        raise RuntimeError(
            "The pre-event active channel has no valid post-event BSI value, "
            "so t_bsi cannot be computed. No default threshold was used.")
    return float(np.percentile(values, percentile))


def dark_mask(b8_pre: np.ndarray, b8_post: np.ndarray, channel: np.ndarray,
              threshold: float) -> np.ndarray:
    """B8 below threshold in the pre or post composite, outside the channel (14.19 a)."""
    with np.errstate(invalid="ignore"):
        dark = (np.asarray(b8_pre) < threshold) | (np.asarray(b8_post) < threshold)
    return dark & ~np.asarray(channel, dtype=bool)


def optical_impact(ndvi_pre: np.ndarray, ndvi_post: np.ndarray,
                   bsi_post: np.ndarray, channel: np.ndarray, t_bsi: float,
                   ndvi_unknown_max: float, ndvi_deposit_max: float,
                   dark: np.ndarray | None = None) -> np.ndarray:
    """deposit state AND not unknown AND not active channel (5.2 steps 4–6).

    Unknown (already bare or built outside the channel), dark pixels (14.19 a)
    and missing data are NaN. Channel pixels are 0: a deposit there is not a
    mark on the bank.
    """
    channel = np.asarray(channel, dtype=bool)
    valid = np.isfinite(ndvi_pre) & np.isfinite(ndvi_post) & np.isfinite(bsi_post)
    with np.errstate(invalid="ignore"):
        unknown = ~channel & (ndvi_pre < ndvi_unknown_max)
        deposit = (ndvi_post < ndvi_deposit_max) & (bsi_post > t_bsi)
    if dark is not None:
        unknown |= np.asarray(dark, dtype=bool)
    out = np.where(deposit & ~channel, 1.0, 0.0)
    out[unknown | ~valid] = np.nan
    return out


def _ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    num = np.asarray(numerator, dtype=float)
    den = np.asarray(denominator, dtype=float)
    out = np.full(num.shape, np.nan, dtype=float)
    ok = np.isfinite(num) & np.isfinite(den) & (den != 0)
    out[ok] = num[ok] / den[ok]
    return out
