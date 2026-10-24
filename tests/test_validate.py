"""The validation sampler draws stations and writes no forest or outcome fields."""

import numpy as np
import pandas as pd
import pytest

from floodforest.validate import (
    coverage_rows,
    draw_sample,
    eligible_station_forest,
    far_from_channel,
    group_sizes,
    sample_layers,
    station_group_labels,
    tercile_groups,
)


def _values(pairs):
    return pd.Series(dict(pairs), name="forest")


def test_equal_groups_and_extras_go_to_the_lower_ones():
    values = _values((i, i // 10) for i in range(30))
    groups = tercile_groups(values, 3)
    assert [len(g) for g in groups] == [10, 10, 10]
    assert set(groups[0]) == set(range(10))
    assert group_sizes(10, 3) == [4, 3, 3]
    drawn = draw_sample(groups, 10, seed=42)
    assert len(drawn) == len(set(drawn)) == 10
    assert [sum(i in set(g) for i in drawn) for g in groups] == [4, 3, 3]
    assert draw_sample(groups, 10, seed=42) == drawn


def test_tied_fractions_split_by_station_id():
    values = _values((i, 0.0) for i in range(6))
    groups = tercile_groups(values, 3)
    assert [g.tolist() for g in groups] == [[0, 1], [2, 3], [4, 5]]


def test_only_analysis_set_banks_count():
    df = pd.DataFrame({
        "station_id": [1, 1, 2, 2],
        "side": ["left", "right", "left", "right"],
        "analysis_set": [0, 0, 1, 0],
        "forest_frac_100": [1.0, 1.0, 0.2, 0.8],
    })
    values = eligible_station_forest(df, "forest_frac_100")
    assert list(values.index) == [2]
    assert values.loc[2] == pytest.approx(0.2)


def test_two_analysis_set_banks_are_averaged():
    df = pd.DataFrame({
        "station_id": [5, 5],
        "side": ["left", "right"],
        "analysis_set": [1, 1],
        "forest_frac_100": [0.0, 1.0],
    })
    values = eligible_station_forest(df, "forest_frac_100")
    assert values.loc[5] == pytest.approx(0.5)


def test_a_missing_forest_fraction_raises():
    df = pd.DataFrame({
        "station_id": [1, 1],
        "side": ["left", "right"],
        "analysis_set": [1, 1],
        "forest_frac_100": [0.4, np.nan],
    })
    with pytest.raises(RuntimeError, match="station 1 right"):
        eligible_station_forest(df, "forest_frac_100")


def test_a_group_smaller_than_its_quota_raises():
    groups = tercile_groups(_values([(1, 0.0)]), 1)
    with pytest.raises(RuntimeError, match="asks for 2"):
        draw_sample(groups, 2, seed=1)


def test_sample_layers_keep_only_place_fields():
    import geopandas as gpd
    from shapely.geometry import LineString

    rows = []
    for station_id, x in ((1, 0.0), (2, 100.0)):
        for side, y in (("left", 10.0), ("right", -10.0)):
            rows.append({
                "station_id": station_id,
                "side": side,
                "river_km": station_id / 10,
                "station_x": x,
                "station_y": 0.0,
                "forest_frac_100": 0.9,
                "reach_distance_m_fab": 40.0,
                "overtopped_fab": 1,
                "geometry": LineString([(x, 0.0), (x, y)]),
            })
    lines = gpd.GeoDataFrame(rows, crs="EPSG:32645")
    points, out = sample_layers(lines, [2, 1])
    assert list(points.columns) == ["station_id", "river_km", "geometry"]
    assert list(out.columns) == ["station_id", "side", "river_km", "geometry"]
    assert points["station_id"].tolist() == [1, 2]
    assert len(out) == 4
    assert points.geometry.iloc[0].coords[0] == (0.0, 0.0)


def test_coverage_is_checked_or_skipped_and_group_counts_hide_values():
    values = _values((i, i) for i in range(6))
    chosen = [0, 1, 2, 3, 4, 5]
    groups = station_group_labels(values, chosen, 3)
    stations = pd.DataFrame({"station_id": chosen, "river_km": [k / 10 for k in chosen]})
    coverage = coverage_rows(stations, {0, 5}, groups)
    assert list(coverage.status) == ["checked", "skipped", "skipped", "skipped", "skipped", "checked"]
    assert set(coverage.forest_group) == {"low", "mid", "high"}
    counts = coverage.groupby(["forest_group", "status"]).size()
    assert int(counts.loc[("low", "checked")]) == 1
    assert int(counts.loc[("high", "checked")]) == 1


def test_far_from_channel_lists_only_sampled_stations():
    table = pd.DataFrame({
        "station_id": [1, 1, 2, 2, 9, 9],
        "river_km": [0.1, 0.1, 0.2, 0.2, 0.9, 0.9],
        "centre_offset_m": [150.0, 150.0, 20.0, 20.0, 200.0, 200.0],
    })
    far = far_from_channel(table, [1, 2], 100)
    assert far.station_id.tolist() == [1]
