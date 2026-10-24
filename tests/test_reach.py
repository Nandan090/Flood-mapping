"""Test 2: a planted impact extent is recovered within one pixel."""

import math

import numpy as np
import pytest

from types import SimpleNamespace

from floodforest.profiles import (bank_top_index, extract_reach, measure_side,
                                  outward_indices, river_visible, valley_wall_index)

STEP = 10.0
Z0 = 1000.0


def _valley(n=200, bank_top=5, reach_px=40, slope_deg=10.0, nodata_at=None,
            gap_at=None):
    """Flat channel, then a bank rising at a constant slope."""
    channel = np.zeros(n, dtype=bool)
    channel[:bank_top] = True
    label = np.zeros(n)
    label[bank_top:bank_top + reach_px] = 1
    if gap_at is not None:
        label[bank_top + gap_at] = 0
    if nodata_at is not None:
        label[nodata_at] = math.nan
    elevation = np.full(n, Z0)
    rise = STEP * math.tan(math.radians(slope_deg))
    for i in range(bank_top, n):
        elevation[i] = Z0 + (i - bank_top) * rise
    slope = np.zeros(n)
    slope[bank_top:] = slope_deg
    return label, elevation, slope, channel


def test_planted_extent_recovered_within_one_pixel():
    reach_px = 40
    label, z, slope, channel = _valley(reach_px=reach_px)
    top = bank_top_index(channel)
    out = extract_reach(label, z, Z0, top, STEP, slope)
    expected_height = reach_px * STEP * math.tan(math.radians(10))
    assert abs(out["reach_distance_m"] - reach_px * STEP) <= STEP
    assert abs(out["reach_height_m"] - expected_height) <= STEP * math.tan(
        math.radians(10)) + 1e-6
    assert out["overtopped"] == 1
    assert out["censored"] == 0


def test_short_non_impact_gap_is_crossed():
    label, z, slope, channel = _valley(reach_px=40, gap_at=10)  # one pixel = 10 m
    top = bank_top_index(channel)
    out = extract_reach(label, z, Z0, top, STEP, slope)
    assert out["reach_distance_m"] == 39 * STEP


def test_long_non_impact_gap_ends_the_run():
    label, z, slope, channel = _valley(reach_px=40)
    label[10:14] = 0  # 40 m of no impact, beyond the 20 m tolerance
    top = bank_top_index(channel)
    out = extract_reach(label, z, Z0, top, STEP, slope)
    assert out["reach_distance_m"] == 4 * STEP
    assert out["censored"] == 0


def test_nodata_gap_rules():
    label, z, slope, channel = _valley(reach_px=40, nodata_at=slice(10, 13))
    top = bank_top_index(channel)
    out = extract_reach(label, z, Z0, top, STEP, slope)
    assert out["reach_distance_m"] == 39 * STEP  # 30 m of nodata does not end it
    assert out["censored"] == 0

    # 50 m of nodata (four steps) is beyond the 30 m allowance, so the run
    # ends at the last impact before the gap and is flagged as censored.
    label, z, slope, channel = _valley(reach_px=40, nodata_at=slice(10, 15))
    out = extract_reach(label, z, Z0, top, STEP, slope)
    assert out["reach_distance_m"] == 4 * STEP
    assert out["censored"] == 1


def test_no_impact_near_the_bank_top_is_zero_reach():
    label, z, slope, channel = _valley(reach_px=0)
    label[20:60] = 1  # impact exists, but detached from the bank
    top = bank_top_index(channel)
    out = extract_reach(label, z, Z0, top, STEP, slope)
    assert out["reach_distance_m"] == 0
    assert out["reach_height_m"] == 0
    assert out["overtopped"] == 0
    assert out["detached_reach_m"] == 54 * STEP


def test_censoring_at_the_valley_wall_and_profile_end():
    label, z, slope, channel = _valley(n=80, reach_px=80)
    top = bank_top_index(channel)
    out = extract_reach(label, z, Z0, top, STEP, slope, valley_wall=30)
    assert out["censored"] == 1
    assert out["reach_edge_index"] == 30

    out = extract_reach(label, z, Z0, top, STEP, slope, valley_wall=79)
    assert out["censored"] == 1
    assert out["reach_edge_index"] == 79


def test_reaching_the_usable_limit_is_censored_at_the_neighbour():
    label, z, slope, _ = _valley(n=80, bank_top=0, reach_px=80)
    out = extract_reach(label, z, Z0, 0, STEP, slope, usable_index=12)
    assert out["censored"] == 1
    assert out["censor_reason"] == "neighbour_limit"
    assert out["reach_edge_index"] == 12
    assert out["reach_distance_m"] == 12 * STEP

    short, z, slope, _ = _valley(n=40, bank_top=0, reach_px=3)
    out = extract_reach(short, z, Z0, 0, STEP, slope, usable_index=20)
    assert out["censored"] == 0
    assert out["censor_reason"] is None
    assert out["reach_edge_index"] == 2


def test_censor_reason_names_the_cause():
    label, z, slope, _ = _valley(n=80, bank_top=0, reach_px=80)
    assert extract_reach(label, z, Z0, 0, STEP, slope)["censor_reason"] == "profile_end"
    assert extract_reach(label, z, Z0, 0, STEP, slope,
                         valley_wall=30)["censor_reason"] == "valley_wall"
    gap = label.copy()
    gap[20:25] = math.nan  # 5 nodata samples span 40 m, longer than 30 m
    out = extract_reach(gap, z, Z0, 0, STEP, slope)
    assert out["censored"] == 1 and out["censor_reason"] == "nodata_gap"


def test_valley_wall_censors_only_an_edge_within_one_sample_of_the_wall():
    wall = 30
    for reach_px, censored in ((31, 1), (30, 1), (28, 0)):
        label, z, slope, _ = _valley(n=80, bank_top=0, reach_px=reach_px)
        out = extract_reach(label, z, Z0, 0, STEP, slope, valley_wall=wall)
        edge = min(reach_px - 1, wall)
        assert out["reach_edge_index"] == edge
        assert out["censored"] == censored
        assert out["censor_reason"] == ("valley_wall" if censored else None)
    # The walk stops at the wall: impact beyond it is not part of the run.
    label = np.zeros(80)
    label[:10] = 1
    label[31:60] = 1
    out = extract_reach(label, z, Z0, 0, STEP, slope, valley_wall=11)
    assert out["reach_edge_index"] == 9 and out["censored"] == 0


def test_channel_wholly_on_one_side_walks_outward_from_each_end():
    from floodforest.profiles import centre_channel
    offsets = np.arange(-300, 301, 10.0)
    channel = (offsets >= -90) & (offsets <= -50)  # 50 m to the right of the station
    centred = centre_channel(channel, offsets, np.full(offsets.size, Z0), STEP, 20, 100)
    assert centred["channel_resolved"] == 1
    assert centred["bank_top_left_m"] == -50 and centred["bank_top_right_m"] == -90
    # The left bank top lies on the right of the station, so its distance from
    # the station along the left normal (bank_top_from_station_m) is negative.
    assert 1.0 * centred["bank_top_left_m"] < 0
    assert -1.0 * centred["bank_top_right_m"] == 90
    left, _ = outward_indices(offsets, "left", centred["bank_top_left_m"])
    right, _ = outward_indices(offsets, "right", centred["bank_top_right_m"])
    assert offsets[left][0] == -50 and np.all(np.diff(offsets[left]) > 0)
    assert offsets[right][0] == -90 and np.all(np.diff(offsets[right]) < 0)
    # Neither walk crosses the channel run: only its own end is a channel sample.
    assert channel[left].sum() == 1 and channel[right].sum() == 1


def test_the_samples_beside_the_bank_decide_how_the_row_starts():
    nan = math.nan
    z = np.full(12, Z0)
    slope = np.zeros(12)

    def start(values, offset=0.0):
        label = np.array(values + [0.0] * (12 - len(values)))
        return extract_reach(label, z, Z0, 0, STEP, slope, start_offset_m=offset)

    # Rule 1: impact within 20 m starts the run, even after a non-impacted
    # sample at 10 m.
    out = start([nan, 0, 1, 1])
    assert out["near_bank_nodata"] == 0 and out["reach_edge_index"] == 3
    # Rule 2: a valid non-impacted sample within 20 m and no impact there: reach 0.
    for values in ([nan, nan, 0, 1], [nan, 0, nan, 1]):
        out = start(values)
        assert out["reach_distance_m"] == 0 and out["near_bank_nodata"] == 0
    # Rule 3: all no-data within 20 m, then impacted at 30 m: the run starts there.
    out = start([nan, nan, nan, 1, 1])
    assert out["near_bank_nodata"] == 0 and out["reach_edge_index"] == 4
    # Otherwise impacted only at 40 m, or the first valid sample is a
    # non-impacted one at 30 m: outcomes missing.
    for values in ([nan, nan, nan, nan, 1], [nan, nan, nan, 0, 1], [nan] * 12):
        out = start(values)
        assert out["near_bank_nodata"] == 1
        assert math.isnan(out["reach_distance_m"]) and math.isnan(out["overtopped"])
    # A valid 0 on the bank-top sample itself (the channel edge) does not decide.
    out = start([0, nan, 1, 1])
    assert out["near_bank_nodata"] == 0 and out["reach_edge_index"] == 3
    # The 5 m offset of an unresolved bank top counts toward the distances.
    assert start([nan, nan, 0], offset=5.0)["near_bank_nodata"] == 1
    assert start([nan, 0], offset=5.0)["reach_distance_m"] == 0


def test_valley_wall_needs_50_m_of_steep_slope_beyond_the_bank_top():
    slope = np.array([40, 10, 35, 35, 35, 35, 35, 10, 31, 31, 31, 31, 31, 31, 5.0])
    # Index 0 is the bank top itself. 2-6 is five samples (40 m): too short.
    assert valley_wall_index(slope, STEP, 30, 50) == 8
    assert valley_wall_index(np.full(5, 45.0), STEP, 30, 50) is None


def test_river_visible_within_30_m_of_the_station():
    offsets = np.arange(-60, 61, 10.0)
    ndvi = np.full(offsets.size, 0.8)
    label = np.zeros(offsets.size)
    assert river_visible(offsets, ndvi, label, 30, 0.30) == 0
    ndvi[offsets == 40] = 0.1  # outside the window
    assert river_visible(offsets, ndvi, label, 30, 0.30) == 0
    label[offsets == -30] = 1  # impacted, at exactly 30 m
    assert river_visible(offsets, ndvi, label, 30, 0.30) == 1


def test_outward_profile_starts_at_the_bank_top_on_each_side():
    offsets = np.arange(-50, 51, 10.0)
    idx, extra = outward_indices(offsets, "left", 20.0)
    assert offsets[idx].tolist() == [20, 30, 40, 50] and extra == 0
    idx, extra = outward_indices(offsets, "right", -5.0)  # unresolved channel
    assert offsets[idx].tolist() == [-10, -20, -30, -40, -50] and extra == 5
    # A run entirely on the right: the left bank top is right of the station.
    idx, _ = outward_indices(offsets, "left", -20.0)
    assert offsets[idx][0] == -20


def test_unresolved_bank_offset_is_added_to_the_distance():
    cfg = SimpleNamespace(impact_gap_m=20, nodata_gap_m=30, overtop_m=20,
                          near_bank_impact_m=30, near_bank_zero_m=20)
    label = np.array([1, 1, 0, 0, 0, 0, 0, 0.0])
    z = np.full(8, Z0)
    out = measure_side(label, z, np.zeros(8), Z0, STEP, 5.0, None, False, cfg)
    assert out["reach_distance_m"] == 15.0
    assert out["overtopped"] == 0
    assert out["detached_reach_m"] == 15.0


def test_valid_fraction_leaves_out_the_channel_bank_top_sample():
    cfg = SimpleNamespace(impact_gap_m=20, nodata_gap_m=30, overtop_m=20,
                          near_bank_impact_m=30, near_bank_zero_m=20)
    label = np.array([math.nan, 1, 0, math.nan, 0, 0])
    z = np.full(6, Z0)
    resolved = measure_side(label, z, np.zeros(6), Z0, STEP, 0.0, 4, False, cfg)
    assert resolved["valid_fraction"] == 0.75  # samples 1-4
    unresolved = measure_side(label, z, np.zeros(6), Z0, STEP, 5.0, 4, False, cfg)
    assert unresolved["valid_fraction"] == 0.6  # samples 0-4
