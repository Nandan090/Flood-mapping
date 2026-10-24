import numpy as np
import pytest
from shapely.geometry import LineString, box

from floodforest.centerline import (choose_status, crosses_river_again,
                                    excluded_as_neighbour, hits_polygon,
                                    manual_reason, station_river_km,
                                    stations_near_confluence, usable_lengths,
                                    verified_direction)


def test_usable_length_is_the_nearest_crossing_within_1_km():
    own = np.array([[0.0, 0.0], [0.0, 500.0]])
    crosses_at_80 = np.array([[100.0, 80.0], [-20.0, 80.0]])
    far_station = np.array([[0.0, 40.0], [20.0, 40.0]])  # would cross at 40 m
    lengths = usable_lengths(
        [own, crosses_at_80, far_station],
        river_m=np.array([0.0, 100.0, 1500.0]),
        search_m=1000, full_length_m=1500)
    assert lengths[0] == pytest.approx(80.0)
    assert lengths[1] == pytest.approx(100.0)
    assert lengths[2] == 1500.0
    assert not excluded_as_neighbour(100.0, 100.0)
    assert excluded_as_neighbour(lengths[0], 100.0)


def test_direction_uses_250_m_on_the_verified_line_and_clamps_at_the_ends():
    # East, then a short northward kink, then east again.
    river = LineString([(0, 0), (1000, 0), (1100, 100), (2500, 100)])
    dx, dy = verified_direction(river, 0.0, 250.0)
    assert (dx, dy) == pytest.approx((1.0, 0.0))
    dx, dy = verified_direction(river, 1000.0, 250.0)
    assert dx > 0.95 and dy > 0
    dx, dy = verified_direction(river, river.length, 250.0)
    assert (dx, dy) == pytest.approx((1.0, 0.0))


def test_a_second_meeting_with_the_river_is_a_meander():
    river = LineString([(0, 0), (0, 1000), (400, 1000), (400, 0)])
    across = np.column_stack([np.linspace(0, 400, 41), np.full(41, 100.0)])
    assert crosses_river_again(across, river, separation_m=100)
    # Two hits 40 m apart along the transect, though far apart along the river.
    narrow = LineString([(0, 0), (0, 1000), (40, 1000), (40, 0)])
    close = np.column_stack([np.linspace(0, 40, 5), np.full(5, 100.0)])
    assert not crosses_river_again(close, narrow, separation_m=100)
    once = np.column_stack([np.linspace(-20, 20, 5), np.full(5, 100.0)])
    assert not crosses_river_again(once, river, separation_m=100)


def test_confluence_flags_cover_300_m_and_a_small_rise_does_not():
    km = np.round(np.arange(0, 5.01, 0.1), 1)
    upa = np.full(km.shape, 100.0)
    upa[km >= 2.0] = 110.0
    near, found = stations_near_confluence(km, upa, 300, 0.05, 300)
    assert found == [2.0]
    assert near[np.isclose(km, 1.7)] and near[np.isclose(km, 2.3)]
    assert not near[np.isclose(km, 1.6)] and not near[np.isclose(km, 2.4)]

    small = np.full(km.shape, 100.0)
    small[km >= 2.0] = 104.0
    near, found = stations_near_confluence(km, small, 300, 0.05, 300)
    assert found == [] and not near.any()

    zero = np.zeros(km.shape)
    zero[km >= 2.0] = 50.0
    _, found = stations_near_confluence(km, zero, 300, 0.05, 300)
    assert found == []

    # Losing the river for one sample, then finding it again, is not a tributary.
    dipped = np.full(km.shape, 100.0)
    dipped[np.isclose(km, 1.0)] = 20.0
    dipped[km >= 2.0] = 160.0
    _, found = stations_near_confluence(km, dipped, 300, 0.05, 300)
    assert found == [2.0]


def test_manual_exclusion_includes_its_ends_and_only_its_bank():
    exclusions = [{"from_km": 3.7, "to_km": 4.7, "side": "right",
                   "reason": "hillslope_failure"}]
    assert manual_reason(3.7, "right", exclusions) == "hillslope_failure"
    assert manual_reason(4.7, "right", exclusions) == "hillslope_failure"
    assert manual_reason(4.8, "right", exclusions) is None
    assert manual_reason(4.0, "left", exclusions) is None


def test_reservoir_polygon_flags_only_a_transect_that_meets_it():
    transect = np.array([[0.0, 0.0], [0.0, 100.0]])
    assert hits_polygon(transect, [box(-10, 50, 10, 60)])
    assert not hits_polygon(transect, [box(50, 50, 60, 60)])
    assert not hits_polygon(transect, [])


def test_status_colour_prefers_the_manual_reason():
    assert choose_status([]) == "kept"
    assert choose_status(["intersects_neighbour", "tributary_fan"]) == "tributary_fan"
    assert choose_status(["outside_reach"]) == "outside_reach"


def test_station_km_follows_the_verified_line_and_must_not_go_backwards():
    river = LineString([(0, 0), (1000, 0)])
    points = np.array([[0.0, 5.0], [500.0, -5.0], [1000.0, 5.0]])
    assert station_river_km(points, river) == pytest.approx([0.0, 0.5, 1.0])
    with pytest.raises(RuntimeError):
        station_river_km(np.array([[500.0, 0.0], [100.0, 0.0]]), river)
