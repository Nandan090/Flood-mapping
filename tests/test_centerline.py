"""The OSM draft finds the flat from a profile and is never marked verified."""

import json

import numpy as np
import pytest

from floodforest.centerline import (CenterlineNotVerified,
                                    assert_centerline_verified,
                                    downstream_end_of_flat)


def test_flat_ends_where_the_river_steepens():
    # 3 km of steep ground, then 2 km dropping only 40 m, then a gorge.
    distance = np.arange(0, 6001, 100)
    elevation = np.empty_like(distance, dtype=float)
    elevation[distance <= 3000] = 4200 - 0.20 * distance[distance <= 3000]
    gentle = (distance > 3000) & (distance <= 5000)
    elevation[gentle] = 3600 - 0.02 * (distance[gentle] - 3000)
    steep = distance > 5000
    elevation[steep] = 3560 - 0.25 * (distance[steep] - 5000)

    end = downstream_end_of_flat(distance, elevation, max_slope=0.08,
                                 window_m=500, min_length_m=1500,
                                 min_elevation_m=3000)
    assert end == pytest.approx(5000, abs=100)


def test_no_flat_is_an_error_rather_than_a_guess():
    distance = np.arange(0, 5001, 100)
    elevation = 4000 - 0.2 * distance
    with pytest.raises(RuntimeError):
        downstream_end_of_flat(distance, elevation, max_slope=0.08,
                               window_m=500, min_length_m=1500,
                               min_elevation_m=3000)


def test_unverified_centerline_is_refused(tmp_path):
    path = tmp_path / "centerline.geojson"
    path.write_text(json.dumps({"type": "FeatureCollection", "features": [{
        "type": "Feature",
        "properties": {"role": "centerline", "verified_by_user": False},
        "geometry": {"type": "LineString", "coordinates": [[0, 0], [1, 1]]},
    }]}), encoding="utf-8")
    with pytest.raises(CenterlineNotVerified):
        assert_centerline_verified(path)
