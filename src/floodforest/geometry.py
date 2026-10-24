"""Centerline smoothing, stations and transect normals (spec section 6).

Left and right are defined looking downstream. For a downstream direction
(dx, dy) the left normal is (-dy, dx) and the right normal is (dy, -dx).
"""

from __future__ import annotations

import numpy as np
from shapely.geometry import LineString, Point


def _unit(vector: np.ndarray) -> np.ndarray:
    length = np.linalg.norm(vector)
    if length == 0:
        raise ValueError("Cannot normalise a zero-length direction vector.")
    return vector / length


def chaikin(points: np.ndarray, iterations: int = 3) -> np.ndarray:
    """Chaikin corner cutting. Endpoints are kept fixed."""
    pts = np.asarray(points, dtype=float)
    for _ in range(iterations):
        nxt = [pts[0]]
        for a, b in zip(pts[:-1], pts[1:]):
            nxt.append(0.75 * a + 0.25 * b)
            nxt.append(0.25 * a + 0.75 * b)
        nxt.append(pts[-1])
        pts = np.vstack(nxt)
    return pts


def densify(points: np.ndarray, step_m: float) -> np.ndarray:
    """Resample a polyline so consecutive vertices are exactly step_m apart.

    The final vertex is the original endpoint, so the last gap may be shorter.
    """
    pts = np.asarray(points, dtype=float)
    out = [pts[0]]
    carried = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        seg = b - a
        seg_len = np.linalg.norm(seg)
        if seg_len == 0:
            continue
        direction = seg / seg_len
        pos = step_m - carried
        while pos <= seg_len + 1e-9:
            out.append(a + direction * min(pos, seg_len))
            pos += step_m
        carried = seg_len - (pos - step_m)
    if np.linalg.norm(out[-1] - pts[-1]) > 1e-6:
        out.append(pts[-1])
    return np.vstack(out)


def _cumulative(points: np.ndarray) -> np.ndarray:
    d = np.zeros(len(points))
    d[1:] = np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))
    return d


def _point_at(points: np.ndarray, cum: np.ndarray, distance: float) -> np.ndarray:
    distance = float(np.clip(distance, 0, cum[-1]))
    i = np.searchsorted(cum, distance, side="right") - 1
    i = min(i, len(points) - 2)
    span = cum[i + 1] - cum[i]
    t = 0 if span == 0 else (distance - cum[i]) / span
    return points[i] + t * (points[i + 1] - points[i])


def stations(points: np.ndarray, spacing_m: float,
             direction_offset_m: float = 100.0) -> list[dict]:
    """Stations every spacing_m along the line, upstream to downstream.

    river_km is the distance downstream from the start in kilometres.
    segment_id is floor(river_km), giving 1 km blocks.
    """
    pts = np.asarray(points, dtype=float)
    cum = _cumulative(pts)
    out = []
    distance = 0.0
    while distance <= cum[-1] + 1e-9:
        here = _point_at(pts, cum, min(distance, cum[-1]))
        ahead = _point_at(pts, cum, min(distance + direction_offset_m, cum[-1]))
        behind = _point_at(pts, cum, max(distance - direction_offset_m, 0))
        direction = _unit(ahead - behind)
        river_km = distance / 1000.0
        out.append({
            "x": float(here[0]),
            "y": float(here[1]),
            "river_km": river_km,
            "segment_id": int(np.floor(river_km)),
            "dx": float(direction[0]),
            "dy": float(direction[1]),
        })
        if distance >= cum[-1]:
            break
        distance += spacing_m
    return out


def normals(dx: float, dy: float) -> dict[str, np.ndarray]:
    """Unit normals. Left and right are as seen looking downstream."""
    return {
        "left": _unit(np.array([-dy, dx])),
        "right": _unit(np.array([dy, -dx])),
    }


def transect(station: dict, side: str, length_m: float,
             step_m: float) -> np.ndarray:
    """Sample points from the station outward along one normal, every step_m."""
    n = normals(station["dx"], station["dy"])[side]
    origin = np.array([station["x"], station["y"]])
    count = int(round(length_m / step_m))
    return np.vstack([origin + n * step_m * i for i in range(count + 1)])


def inside_reach(river_km: np.ndarray, reach_start_km: float,
                 reach_end_km: float | None) -> np.ndarray:
    """Stations between the configured start and end (14.5, 14.16 rule 1).

    Stations outside stay in the table, flagged False. There is no
    open-ended reach, so a missing end is an error.
    """
    if reach_end_km is None:
        raise ValueError(
            "reach_end_km is not set for this event; the reach cannot be open-ended.")
    km = np.asarray(river_km, dtype=float)
    return (km >= reach_start_km) & (km <= reach_end_km)


def river_sections(river_km: np.ndarray,
                   section_breaks_km: list[float] | None) -> np.ndarray:
    """Section number per station: 0 above the first break, 1 below it, ...

    A station exactly at a break belongs to the section below it.
    """
    if section_breaks_km is None:
        raise ValueError(
            "section_breaks_km has not been confirmed for this event. Set it in "
            "config.yaml, using [] if the reach has no confluence.")
    breaks = np.sort(np.asarray(section_breaks_km, dtype=float))
    return np.searchsorted(breaks, np.asarray(river_km, dtype=float), side="right")


def propose_section_breaks(river_km: np.ndarray, upa: np.ndarray,
                           window_m: float, rise: float) -> list[float]:
    """Candidate confluences for the user to confirm (14.16 rule 4).

    A candidate is where upstream area grows by more than `rise` (a share)
    within `window_m` downstream. Each jump is reported once, at the first
    station past the rise.
    """
    km = np.asarray(river_km, dtype=float)
    area = np.asarray(upa, dtype=float)
    window_km = window_m / 1000.0
    candidates: list[float] = []
    for i in range(len(km)):
        ahead = np.flatnonzero((km > km[i]) & (km <= km[i] + window_km))
        jumped = ahead[area[ahead] > area[i] * (1 + rise)]
        if jumped.size:
            at = float(km[jumped[0]])
            if not candidates or at - candidates[-1] > window_km:
                candidates.append(at)
    return candidates


def landmark_distances(line_xy: np.ndarray,
                       landmarks_xy: dict[str, tuple[float, float]]) -> dict:
    """Along-river km from the line start to each landmark (14.16 rule 6).

    Coordinates must already be in the projected CRS. Also returns how far
    each landmark sits off the line, so a misplaced point is easy to spot.
    """
    line = LineString(np.asarray(line_xy, dtype=float))
    out = {}
    for name, (x, y) in landmarks_xy.items():
        point = Point(x, y)
        out[name] = {
            "river_km": line.project(point) / 1000.0,
            "offset_m": line.distance(point),
        }
    return out
