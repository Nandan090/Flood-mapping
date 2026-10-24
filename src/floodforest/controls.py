"""Terrain controls per transect (spec section 7.4).

Curvature is positive where the river turns left looking downstream. The
outer bank of a left turn is the right bank.
"""

from __future__ import annotations

import math

import numpy as np


def gradient_per_km(river_km: np.ndarray, z0: np.ndarray, window_m: float,
                    min_points: int = 3) -> np.ndarray:
    """Least-squares slope of z0 against river_km over stations within
    ±window_m (metres of z0 per km downstream; negative where the bed falls).
    NaN when fewer than min_points stations have a z0."""
    km = np.asarray(river_km, dtype=float)
    z = np.asarray(z0, dtype=float)
    out = np.full(len(km), np.nan)
    for i in range(len(km)):
        near = (np.abs(km - km[i]) * 1000.0 <= window_m + 1e-6) & np.isfinite(z)
        if near.sum() >= min_points and np.ptp(km[near]) > 0:
            out[i] = float(np.polyfit(km[near], z[near], 1)[0])
    return out


def curvature_per_km(river, distance_m: np.ndarray, window_m: float) -> np.ndarray:
    """Signed turn of the verified line at each station, radians per km.

    The heading from window_m upstream to the station is compared with the
    heading from the station to window_m downstream, with both points clamped
    to the line ends. A station whose upstream or downstream arm has no length
    (the first or last vertex) gets NaN.
    """
    length = float(river.length)
    out = np.full(len(distance_m), np.nan)
    for i, d in enumerate(np.asarray(distance_m, dtype=float)):
        behind = max(0.0, d - window_m)
        ahead = min(length, d + window_m)
        here = np.array(river.interpolate(d).coords[0])
        u = here - np.array(river.interpolate(behind).coords[0])
        v = np.array(river.interpolate(ahead).coords[0]) - here
        lu, lv = float(np.linalg.norm(u)), float(np.linalg.norm(v))
        if lu < 1.0 or lv < 1.0:
            continue
        turn = math.atan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1])
        out[i] = turn / ((lu + lv) / 1000.0)
    return out


def outer_bend(curvature: np.ndarray, side: np.ndarray) -> np.ndarray:
    """1 if this bank is the outside of the bend, else 0. A straight reach
    (curvature 0) or unknown curvature gives 0 on both banks."""
    c = np.asarray(curvature, dtype=float)
    s = np.asarray(side)
    return (((c > 0) & (s == "right")) | ((c < 0) & (s == "left"))).astype(int)


def max_within(values: np.ndarray, transform, points: np.ndarray,
               radius_m: float) -> np.ndarray:
    """Largest finite value among pixels whose centres lie within radius_m of
    the centre of each point's pixel. NaN when there is none."""
    from rasterio.transform import rowcol
    pixel = abs(float(transform.a))
    reach = int(math.ceil(radius_m / pixel))
    dr, dc = np.mgrid[-reach:reach + 1, -reach:reach + 1]
    disc = np.hypot(dr, dc) * pixel <= radius_m + 1e-6
    rows, cols = rowcol(transform, points[:, 0], points[:, 1])
    height, width = values.shape
    out = np.full(len(points), np.nan)
    for i, (row, col) in enumerate(zip(np.asarray(rows), np.asarray(cols))):
        rr = row + dr[disc]
        cc = col + dc[disc]
        keep = (rr >= 0) & (cc >= 0) & (rr < height) & (cc < width)
        window = values[rr[keep], cc[keep]]
        if np.isfinite(window).any():
            out[i] = float(np.nanmax(window))
    return out
