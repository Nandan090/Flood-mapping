"""Test 3: a belt with a known forest share returns that share."""

import numpy as np
import pytest

from floodforest.exposure import forest_fraction, forest_mask, hansen_loss_window


def test_known_forest_share():
    trees = np.array([0.9, 0.9, 0.9, 0.1, 0.9, 0.2])
    canopy = np.array([12, 8, 3, 20, 6, 9])       # pixel 2 is too short
    loss = np.array([0, 0, 0, 0, 2022, 0])        # pixel 4 was cleared
    forest = forest_mask(trees, canopy, loss, pre_event_year=2024)
    # Only pixels 0 and 1 qualify, so 2 of 6.
    assert forest_fraction(forest) == np.float64(2 / 6)


def test_invalid_pixels_are_excluded_from_the_share():
    forest = np.array([True, True, False, False])
    valid = np.array([True, False, True, False])
    assert forest_fraction(forest, valid) == 0.5


def test_empty_belt_raises():
    with np.testing.assert_raises(ValueError):
        forest_fraction(np.array([]))


def test_hansen_condition_skipped_when_pre_event_year_is_before_2021():
    trees = np.array([0.8, 0.8])
    canopy = np.array([10, 10])
    loss = np.array([2020, 0])  # 2020 loss is already in the canopy layer
    forest = forest_mask(trees, canopy, loss, pre_event_year=2020)
    assert forest.tolist() == [True, True]


def test_forest_without_the_canopy_height_condition():
    trees = np.array([0.9, 0.9])
    canopy = np.array([2.0, 9.0])  # first pixel grew past 5 m only after 2020
    loss = np.array([0, 0])
    strict = forest_mask(trees, canopy, loss, pre_event_year=2024)
    relaxed = forest_mask(trees, canopy, loss, pre_event_year=2024,
                          use_canopy_height=False)
    assert strict.tolist() == [False, True]
    assert relaxed.tolist() == [True, True]


def test_hansen_window_is_skipped_for_a_2020_pre_composite():
    decision = hansen_loss_window(2020, newest_year=2025)
    assert decision["applied"] is False
    assert decision["uncovered_years"] == []


def test_belt_on_a_raster_returns_the_planted_forest_share():
    from rasterio.transform import from_origin

    from floodforest.exposure import belt_metrics, belt_pixel_indices, belt_polygon
    # 10 m pixels; the belt covers columns 2-11 and rows 0-9 (100 x 100 m).
    shape = (20, 20)
    transform = from_origin(0, 200, 10, 10)
    forest = np.zeros(shape)
    forest[0:10, 2:5] = 1.0          # 30 forest pixels in one patch
    forest[0:10, 9] = 1.0            # 10 more, a separate patch
    forest[0:5, 11] = np.nan         # 5 invalid pixels
    layers = {name: np.zeros(shape) for name in
              ("forest_dw", "canopy_height", "ndvi_pre", "prior_loss", "worldcover_forest")}
    layers["forest"] = forest
    source = np.zeros(shape)
    source[0:10, 9] = 1.0            # one column decided by WorldCover
    layers["forest_source"] = source
    flat = {name: values.ravel() for name, values in layers.items()}
    # Station at the left edge of the belt, normal pointing east.
    belt = belt_polygon((20, 150), (0, 1), (1, 0), 0, 100, 100)
    index = belt_pixel_indices([belt, None], shape, transform)
    assert index[0].size == 100 and index[1].size == 0
    out = belt_metrics(index[0], flat, shape[1])
    assert out["belt_pixels"] == 95
    assert out["forest_frac"] == pytest.approx(40 / 95)
    assert out["belt_continuity"] == pytest.approx(30 / 40)
    assert out["forest_source"] == "mixed"
    assert out["forest_source_worldcover_frac"] == pytest.approx(10 / 95)


def test_low_tree_probability_is_non_forest_without_a_canopy_height():
    from floodforest.exposure import forest_with_unknowns
    trees = np.array([0.1, 0.9, 0.9, np.nan])
    height = np.array([np.nan, np.nan, 12.0, 12.0])
    out, source = forest_with_unknowns(trees, height, np.zeros(4), 2020, 0.5, 5, 2021, 2025)
    assert out[0] == 0 and np.isnan(out[1]) and out[2] == 1 and np.isnan(out[3])
    assert source[0] == 0 and np.isnan(source[1]) and np.isnan(source[3])


def test_worldcover_stands_in_only_where_dynamic_world_is_missing():
    from floodforest.exposure import forest_with_unknowns
    trees = np.array([np.nan, np.nan, np.nan, np.nan, 0.1])
    height = np.array([12.0, 3.0, 12.0, np.nan, 12.0])
    wc = np.array([10.0, 10.0, 30.0, 10.0, 10.0])  # last pixel: Dynamic World decides
    loss = np.zeros(5)
    out, source = forest_with_unknowns(trees, height, loss, 2020, 0.5, 5, 2021, 2025,
                                       wc, 10)
    assert out[:3].tolist() == [1, 0, 0] and np.isnan(out[3]) and out[4] == 0
    assert source[:3].tolist() == [1, 1, 1] and np.isnan(source[3]) and source[4] == 0
    both_missing, _ = forest_with_unknowns(np.array([np.nan]), np.array([12.0]), np.zeros(1),
                                           2020, 0.5, 5, 2021, 2025, np.array([np.nan]), 10)
    assert np.isnan(both_missing[0])


def test_hansen_window_reports_years_past_the_newest_map():
    decision = hansen_loss_window(2026, newest_year=2024)
    assert decision["applied"] is True
    assert decision["loss_through_year"] == 2024
    assert decision["uncovered_years"] == [2025, 2026]
