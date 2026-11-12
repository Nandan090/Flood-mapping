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
    prepare_validation_comparison,
    standardized_rank_slope,
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


def _comparison_fixture():
    import geopandas as gpd
    from shapely.geometry import LineString, Point

    sides = [(1, "left", 60.0), (2, "left", 60.0),
             (3, "right", 40.0), (328, "right", 60.0), (4, "left", 40.0)]
    edge_rows, line_rows, table_rows = [], [], []
    for station_id, side, along in sides:
        sign = 1 if side == "left" else -1
        line = LineString([(0.0, 0.0), (0.0, sign * 100.0)])
        line_rows.append({"station_id": station_id, "side": side, "geometry": line})
        edge_rows.append({
            "station_id": station_id, "side": side, "image_date": "2022-01-17",
            "confidence": 1, "note": "dot is slightly off" if station_id == 328 else "",
            "geometry": Point(0.0, sign * along),
        })
        table_rows.append({
            "station_id": station_id, "side": side, "block_id": station_id // 10,
            "station_x": 0.0, "station_y": 0.0,
            "bank_top_from_station_m": {1: 20.0, 2: 20.0, 3: 30.0,
                                         328: 20.0, 4: 10.0}[station_id],
            "reach_distance_m_fab": {1: 30.0, 2: 40.0, 3: np.nan,
                                     328: 40.0, 4: 25.0}[station_id],
            "valid_fraction_fab": {1: .9, 2: .9, 3: .8, 328: .9, 4: .5}[station_id],
            "near_bank_nodata_fab": int(station_id == 3),
            "censored_fab": np.nan if station_id == 3 else int(station_id == 2),
            "analysis_set": int(station_id in (1, 2, 328)),
            "exclude_reasons": {3: "near_bank_nodata,river_not_visible",
                                4: "low_valid_fraction"}.get(station_id, ""),
            "forest_frac_100": {1: .1, 2: .4, 3: .8, 328: .9, 4: .2}[station_id],
            "detached_reach_m_fab": {1: 50.0, 2: 40.0, 3: 30.0,
                                     328: 40.0, 4: 25.0}[station_id],
            "z0_fab": 100.0,
        })
    edges = gpd.GeoDataFrame(edge_rows, crs="EPSG:32645")
    lines = gpd.GeoDataFrame(line_rows, crs="EPSG:32645")
    table = pd.DataFrame(table_rows)
    coverage = pd.DataFrame({
        "station_id": [1, 2, 3, 328, 4],
        "forest_group": ["low", "mid", "high", "high", "low"],
    })
    return edges, lines, table, coverage


def test_comparison_sign_convention_exclusions_and_groups():
    edges, lines, table, coverage = _comparison_fixture()
    rows = prepare_validation_comparison(edges, lines, table, coverage)
    by_id = rows.set_index("station_id")
    assert by_id.loc[1, "validated_reach_m"] == pytest.approx(40.0)
    assert by_id.loc[1, "diff_m"] == pytest.approx(-10.0)
    assert bool(by_id.loc[1, "main_include"])
    assert not bool(by_id.loc[2, "main_include"])
    assert bool(by_id.loc[2, "censored"])
    assert not bool(by_id.loc[3, "main_include"])
    assert bool(by_id.loc[3, "near_bank_nodata"])
    assert bool(by_id.loc[3, "near_bank_zero_sensitivity_include"])
    assert by_id.loc[3, "diff_m_if_near_bank_zero"] == pytest.approx(-10.0)
    assert not bool(by_id.loc[4, "main_include"])
    assert bool(by_id.loc[4, "low_valid_fraction"])
    assert by_id.loc[1, "forest_group"] == "Low"
    assert by_id.loc[2, "forest_group"] == "Mid"
    assert by_id.loc[328, "forest_group"] == "High"
    assert by_id.loc[328, "station_flag"] == "braided reach, far from main channel"
    assert by_id.loc[328, "note"] == "dot is slightly off"


def test_exclusion_reasons_come_from_flags_and_are_reported_separately():
    from floodforest.validate import _exclusion_report_lines

    edges, lines, table, coverage = _comparison_fixture()
    rows = prepare_validation_comparison(edges, lines, table, coverage)
    reasons = rows.set_index("station_id")["main_exclusion_reasons"]
    assert reasons[1] == "" and reasons[328] == ""
    # Censored with no model exclude_reasons must still say why.
    assert reasons[2] == "censored"
    assert reasons[3] == ("near_bank_nodata;model_excluded(river_not_visible);"
                          "missing_computed_reach")
    assert reasons[4] == "low_valid_fraction"

    report = "\n".join(_exclusion_report_lines(rows[~rows.main_include]))
    assert "- Censored (reach edge not observed): 2-left." in report
    assert "- Near-bank nodata: 3-right." in report
    assert "- Low valid fraction: 4-left." in report
    assert "- Other model exclusion: 3-right (river_not_visible)." in report
    assert "- Computed reach unavailable: 3-right." in report
    assert "outcome unavailable" not in report


def test_detached_columns_keep_the_requested_signs():
    edges, lines, table, coverage = _comparison_fixture()
    by_id = prepare_validation_comparison(edges, lines, table, coverage).set_index("station_id")
    # Station 1: validated 40 m, reach 30 m, detached 50 m.
    assert by_id.loc[1, "detached_reach_m"] == pytest.approx(50.0)
    assert by_id.loc[1, "validated_minus_detached_m"] == pytest.approx(-10.0)
    assert by_id.loc[1, "diff_detached_m"] == pytest.approx(10.0)
    assert by_id.loc[1, "diff_m"] == pytest.approx(-10.0)


def test_compare_refuses_a_station_with_two_points_on_one_bank():
    import pandas.errors

    edges, lines, table, coverage = _comparison_fixture()
    twice = pd.concat([edges, edges[edges.station_id == 1]], ignore_index=True)
    with pytest.raises(pandas.errors.MergeError, match="not unique in left dataset"):
        prepare_validation_comparison(twice, lines, table, coverage)


def test_compare_checks_each_event_against_its_own_date_window():
    from floodforest.config import load_config
    from floodforest.validate import check_point_dates, comparison_date_window

    events = load_config().events
    assert comparison_date_window(events.melamchi2021) == ("2021-08-01", "2022-05-31")
    lende = comparison_date_window(events.lende2025)
    assert lende == ("2025-07-09", "2026-05-31")
    check_point_dates(pd.DataFrame({"image_date": ["2025-11-29", "2026-01-20"]}), lende)
    with pytest.raises(RuntimeError, match="capture-date window"):
        check_point_dates(pd.DataFrame({"image_date": ["2025-07-08"]}), lende)


def test_first_pass_review_file_marks_every_point_not_reviewed_and_is_never_overwritten(
        tmp_path):
    from floodforest.validate import load_validation_review, write_review_template

    points = pd.DataFrame({"station_id": [5, 2, 2], "side": ["left", "right", "left"]})
    path = tmp_path / "validation_review.csv"
    assert write_review_template(points, path)
    review = pd.read_csv(path, dtype={"note": str})
    assert review[["station_id", "side"]].values.tolist() == [[2, "left"], [2, "right"],
                                                               [5, "left"]]
    assert set(review.status) == {"not reviewed"}
    assert len(load_validation_review(review, points)) == 3
    path.write_text("station_id,side,status,note" + chr(10) + "2,left,visible change,x" + chr(10),
                    encoding="utf-8")
    assert not write_review_template(points, path)
    assert "visible change" in path.read_text(encoding="utf-8")


def test_no_reviewed_points_says_the_second_pass_is_still_to_do():
    from floodforest.validate import review_scope_limitation

    points = pd.DataFrame({"review_status": ["not reviewed"] * 2, "diff_m": [-30.0, 5.0]})
    assert "still to be done" in review_scope_limitation(points, 20.0)


def test_subset_coverage_counts_unassigned_stations_apart():
    from floodforest.validate import subset_coverage_note

    coverage = pd.DataFrame({"station_id": [1, 2, 3, 4, 5],
                             "status": ["checked", "skipped", "skipped", "skipped", "checked"],
                             "forest_group": ["low", "low", "low", "low", "high"]})
    subset = {"done_before": [1], "to_do": [2, 5], "drawn_on": "2026-11-12", "seed": 43}
    note = subset_coverage_note(coverage, subset)
    assert "2 of 4 sampled low-forest stations were assigned" in note
    assert "1 of those was skipped or removed" in note


def _cutoff_frame():
    return pd.DataFrame({
        "forest_group": ["Low"] * 4 + ["Mid"] * 2 + ["High"] * 4,
        "diff_m": [-5.0, -15.0, -25.0, 3.0, -35.0, 0.0, -12.0, -22.0, -31.0, -40.0],
        "station_id": range(10), "side": ["left"] * 10,
    })


def test_metrics_count_under_detected_at_the_given_cutoff():
    from floodforest.validate import _validation_metrics

    frame = _cutoff_frame()
    at_20 = _validation_metrics(frame, under_detect_m=20.0)
    at_10 = _validation_metrics(frame, under_detect_m=10.0)
    assert at_20["n_under_detected"] == 5
    assert at_20["share_under_detected"] == pytest.approx(0.5)
    assert at_10["n_under_detected"] == 7
    # Exactly the cut-off is not under-detected (strictly less than -cutoff).
    edge = _validation_metrics(pd.DataFrame({"diff_m": [-20.0, -20.1]}), under_detect_m=20.0)
    assert edge["n_under_detected"] == 1


def test_cutoff_table_shares_counts_gap_and_flag():
    from floodforest.validate import under_detection_cutoff_table

    table = under_detection_cutoff_table(_cutoff_frame(), [10, 20, 30], min_group_n=2)
    by_cut = table.set_index("cutoff_m")
    assert by_cut.loc[20, "low_under_detected"] == 1 and by_cut.loc[20, "low_n"] == 4
    assert by_cut.loc[20, "high_under_detected"] == 3 and by_cut.loc[20, "high_n"] == 4
    assert by_cut.loc[20, "high_minus_low"] == pytest.approx(0.75 - 0.25)
    assert by_cut.loc[30, "high_minus_low"] == pytest.approx(0.5 - 0.0)
    assert by_cut.loc[10, "high_minus_low"] == pytest.approx(1.0 - 0.5)
    assert by_cut.flagged.all()
    small = under_detection_cutoff_table(_cutoff_frame(), [20])  # default min n is 8
    assert not bool(small.flagged.iloc[0])


def test_detached_sensitivity_keeps_censored_rows_separate():
    from floodforest.validate import detached_sensitivity

    points = pd.DataFrame({
        "main_include": [True, True, False, False],
        "censored": [False, False, True, True],
        "diff_m": [-30.0, 0.0, -10.0, -50.0],
        "diff_detached_m": [-5.0, 0.0, 20.0, -50.0],
    })
    table = detached_sensitivity(points).set_index(["subset", "computed_outcome"])
    assert table.loc[("main", "reach_distance_m"), "n"] == 2
    assert table.loc[("main", "detached_reach_m"), "mean_diff_m"] == pytest.approx(-2.5)
    assert table.loc[("main", "reach_distance_m"), "n_under_detected"] == 1
    assert table.loc[("main", "detached_reach_m"), "n_under_detected"] == 0
    assert table.loc[("censored", "detached_reach_m"), "n"] == 2
    assert table.loc[("censored", "detached_reach_m"), "mean_diff_m"] == pytest.approx(-15.0)


def _review_points():
    return pd.DataFrame({"station_id": [346, 346, 1], "side": ["left", "right", "left"]})


def test_review_attaches_status_and_refuses_bad_files():
    from floodforest.validate import load_validation_review

    review = pd.DataFrame({
        "station_id": [346, 346, 1], "side": ["left", "right", "left"],
        "status": ["not reviewed", "no visible change", "not reviewed"],
        "note": [np.nan, "plain unchanged", np.nan],
    })
    out = load_validation_review(review, _review_points()).set_index(["station_id", "side"])
    assert out.loc[(346, "right"), "review_status"] == "no visible change"
    assert out.loc[(346, "right"), "review_note"] == "plain unchanged"
    assert out.loc[(1, "left"), "review_note"] == ""
    with pytest.raises(RuntimeError, match="unknown statuses"):
        load_validation_review(review.assign(status="maybe"), _review_points())
    with pytest.raises(RuntimeError, match="exactly one row per validation point"):
        load_validation_review(review.iloc[:2], _review_points())
    with pytest.raises(RuntimeError, match="duplicate"):
        load_validation_review(pd.concat([review, review.iloc[:1]]), _review_points())


def test_review_sensitivity_drops_only_no_visible_change_points():
    from floodforest.validate import REVIEW_SENSITIVITY, _validation_summary

    points = pd.DataFrame({
        "main_include": [True, True, True, False],
        "forest_group": ["Low", "High", "High", "Low"],
        "image_date": ["2022-01-17"] * 4, "confidence": [1, 1, 1, 1],
        "diff_m": [-5.0, -150.0, -25.0, -60.0],
        "diff_m_if_near_bank_zero": [np.nan] * 4,
        "near_bank_zero_sensitivity_include": [False] * 4,
        "is_station_328": [False] * 4,
        "review_status": ["not reviewed", "no visible change", "can't tell",
                          "no visible change"],
    })
    summary = _validation_summary(points)
    row = summary[summary.subset == REVIEW_SENSITIVITY].iloc[0]
    main = summary[(summary.subset == "main") & (summary.breakdown == "all")].iloc[0]
    assert main.n == 3 and row.n == 2
    assert row.mean_diff_m == pytest.approx(-15.0)


def test_review_scope_sentence_only_when_reviewed_errors_exceed_cutoff():
    from floodforest.validate import review_scope_limitation

    points = pd.DataFrame({
        "station_id": [346, 328, 1, 13], "side": ["right", "right", "left", "left"],
        "diff_m": [-150.0, -25.0, -5.0, np.nan],
        "review_status": ["no visible change", "visible change", "not reviewed",
                          "not reviewed"],
    })
    assert review_scope_limitation(points, 20.0) == (
        "Only points with an error above 20 m were reviewed against pre-event imagery; "
        "smaller errors were not.")
    with pytest.raises(RuntimeError, match="Update the text"):
        review_scope_limitation(points.assign(review_status="can't tell"), 20.0)


def test_traced_point_terrain_samples_nearest_pixel_and_subtracts_z0():
    import geopandas as gpd
    from rasterio.transform import from_origin
    from shapely.geometry import Point
    from floodforest.validate import traced_point_terrain

    grid = {"transform": from_origin(0.0, 30.0, 10.0, 10.0), "height": 3, "width": 3}
    ramp = np.arange(9, dtype=float).reshape(3, 3)
    layers = {"z_fab": 100.0 + ramp, "slope_fab": ramp * 5,
              "ndvi_pre": ramp / 10, "ndvi_post": ramp / 20}
    points = gpd.GeoDataFrame({"z0_fab": [100.0, 101.0]},
                              geometry=[Point(15.0, 15.0), Point(25.0, 5.0)],
                              crs="EPSG:32645")
    out = traced_point_terrain(points, layers, grid)
    # (15, 15) is row 1, col 1 -> value 4; (25, 5) is row 2, col 2 -> value 8.
    assert list(out.traced_slope_deg) == [20.0, 40.0]
    assert list(out.traced_height_above_z0_m) == [4.0, 7.0]
    assert list(out.traced_ndvi_pre) == pytest.approx([0.4, 0.8])
    assert list(out.traced_ndvi_post) == pytest.approx([0.2, 0.4])


def test_flood_edge_limitation_text_counts_from_the_data():
    from floodforest.validate import flood_edge_visibility_limitation

    main = pd.DataFrame({"diff_m": [-150.5, -29.9, -28.7, -26.3, -25.7, -21.0, -20.7,
                                    -5.0, 3.0]})
    text = flood_edge_visibility_limitation(main, 20.0, 10.0)
    assert text == (
        "Validation is weaker where the flood edge is not visible in imagery: on wide bare "
        "gravel plains, where pre-event bare ground is unknown to the satellite label and a "
        "traced edge may only mark the old plain boundary, and on steep forested slopes, "
        "where a flood that passed under trees may leave no bare deposit. Six of the seven "
        "largest main-set errors are 21 to 30 m, within a few 10 m pixels.")
    assert "outlier" not in text.lower()


def test_largest_errors_rank_by_absolute_diff_with_stable_ties():
    from floodforest.validate import largest_errors

    main = pd.DataFrame({
        "station_id": [5, 2, 2, 9, 1, 7],
        "side": ["left", "right", "left", "left", "left", "right"],
        "diff_m": [3.0, -40.0, 40.0, 10.0, -1.0, -12.0],
    })
    top = largest_errors(main, k=3)
    assert list(zip(top.station_id, top.side)) == [(2, "left"), (2, "right"), (7, "right")]
    assert list(top.abs_diff_m) == [40.0, 40.0, 12.0]


def test_correlation_interval_is_clipped_and_raw_bound_reported():
    from floodforest.validate import clip_correlation_interval

    (low, high), notes = clip_correlation_interval((-1.17, 0.188))
    assert (low, high) == (-1.0, 0.188)
    assert notes == ["lower bound was -1.170"]
    assert clip_correlation_interval((-0.5, 0.4)) == ((-0.5, 0.4), [])


def test_standardized_rank_slope_matches_spearman_with_ties():
    from scipy.stats import spearmanr

    x = np.array([1.0, 2.0, 2.0, 4.0, 5.0, 5.0])
    y = np.array([5.0, 1.0, 1.0, 3.0, 2.0, 4.0])
    assert standardized_rank_slope(x, y) == pytest.approx(spearmanr(x, y).statistic)
