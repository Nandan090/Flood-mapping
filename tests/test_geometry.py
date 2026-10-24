"""Test 1: transect normals are perpendicular and left/right look downstream."""

import numpy as np
import pytest

from floodforest.geometry import (densify, inside_reach, landmark_distances,
                                  normals, propose_section_breaks,
                                  river_sections, stations, transect)


def _dot(station):
    direction = np.array([station["dx"], station["dy"]])
    return {side: float(np.dot(direction, normals(station["dx"], station["dy"])[side]))
            for side in ("left", "right")}


def test_normals_on_a_straight_river():
    # Flows west to east, so downstream is +x. Left must point north (+y).
    line = np.array([[0.0, 0.0], [1000.0, 0.0]])
    sampled = stations(densify(line, 10), spacing_m=100)
    for s in sampled[1:-1]:
        dots = _dot(s)
        assert abs(dots["left"]) < 1e-9
        assert abs(dots["right"]) < 1e-9
        assert normals(s["dx"], s["dy"])["left"][1] > 0.999
        assert normals(s["dx"], s["dy"])["right"][1] < -0.999


def test_normals_on_a_curved_river():
    # A quarter circle bending from eastward to northward.
    angles = np.linspace(0, np.pi / 2, 400)
    radius = 2000.0
    line = np.column_stack([radius * np.sin(angles), radius * (1 - np.cos(angles))])
    sampled = stations(densify(line, 10), spacing_m=100)
    for s in sampled[2:-2]:
        dots = _dot(s)
        assert abs(dots["left"]) < 1e-6
        assert abs(dots["right"]) < 1e-6
        # Looking downstream, left points to the inside of this left-hand
        # bend, toward the centre of curvature at (0, radius).
        n = normals(s["dx"], s["dy"])["left"]
        toward_centre = np.array([-s["x"], radius - s["y"]])
        assert np.dot(n, toward_centre) > 0


def test_transects_start_at_the_station_and_step_correctly():
    line = np.array([[0.0, 0.0], [500.0, 0.0]])
    s = stations(densify(line, 10), spacing_m=100)[2]
    pts = transect(s, "left", length_m=1500, step_m=10)
    assert len(pts) == 151
    assert np.allclose(pts[0], [s["x"], s["y"]])
    gaps = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    assert np.allclose(gaps, 10.0)


# ---------------------------------------------------------------------------
# Reach definitions and confluences (14.16)
# ---------------------------------------------------------------------------

def test_stations_past_the_reach_end_are_flagged_not_dropped():
    km = np.array([0.0, 10.0, 35.0, 35.1, 50.0])
    assert inside_reach(km, 0, 35).tolist() == [True, True, True, False, False]


def test_open_ended_reach_is_refused():
    with pytest.raises(ValueError):
        inside_reach(np.array([0.0, 1.0]), 0, None)


def test_river_sections_split_at_confirmed_breaks():
    km = np.array([0.0, 12.0, 18.0, 18.1, 30.0])
    assert river_sections(km, [18.0]).tolist() == [0, 0, 1, 1, 1]
    assert river_sections(km, []).tolist() == [0, 0, 0, 0, 0]


def test_unconfirmed_section_breaks_raise():
    with pytest.raises(ValueError):
        river_sections(np.array([0.0, 1.0]), None)


def test_confluence_candidates_where_upstream_area_jumps():
    km = np.round(np.arange(0, 10.01, 0.1), 1)
    upa = np.full(km.shape, 100.0)
    upa[km >= 4.0] = 180.0   # an 80% rise: a tributary joins at 4 km
    upa[km >= 7.0] = 200.0   # an 11% rise: too small to propose
    assert propose_section_breaks(km, upa, window_m=300, rise=0.5) == [4.0]


def test_landmark_distance_along_the_line():
    line = np.array([[0.0, 0.0], [10_000.0, 0.0], [10_000.0, 5_000.0]])
    out = landmark_distances(line, {"bridge": (10_050.0, 2_000.0)})
    assert out["bridge"]["river_km"] == pytest.approx(12.0)
    assert out["bridge"]["offset_m"] == pytest.approx(50.0)
