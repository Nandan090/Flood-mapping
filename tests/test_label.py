import numpy as np
import pytest

from floodforest import optical, sar
from floodforest.label import combine, connectivity_filter, label_source
from floodforest.profiles import centre_channel


def test_bsi_matches_the_spec_formula():
    # ((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))
    value = optical.bsi(np.array(0.05), np.array(0.10), np.array(0.30), np.array(0.25))
    assert float(value) == pytest.approx((0.35 - 0.35) / 0.70)
    value = optical.bsi(np.array(0.05), np.array(0.20), np.array(0.10), np.array(0.30))
    assert float(value) == pytest.approx((0.50 - 0.15) / 0.65)


def test_optical_impact_rules():
    channel = np.array([True, False, False, False])
    ndvi_pre = np.array([0.1, 0.25, 0.6, 0.6])    # bare before -> unknown
    ndvi_post = np.array([0.1, 0.1, 0.1, 0.5])
    bsi_post = np.array([0.3, 0.3, 0.3, 0.3])
    out = optical.optical_impact(ndvi_pre, ndvi_post, bsi_post, channel, 0.1,
                                 0.30, 0.30)
    assert out[0] == 0          # channel itself is never impact
    assert np.isnan(out[1])     # already bare outside the channel
    assert out[2] == 1          # vegetated before, deposit after
    assert out[3] == 0          # still vegetated


def test_channel_component_must_come_near_the_line():
    spectral = np.zeros((6, 6), dtype=bool)
    spectral[0, 0:3] = True     # touches the near-line zone at (0, 2)
    spectral[5, 5] = True       # far away
    near = np.zeros_like(spectral)
    near[0:3, 2] = True
    kept = optical.keep_components_touching(spectral, near)
    assert kept[0, 0] and kept[0, 1] and kept[0, 2]
    assert not kept[5, 5]


def test_aspect_uses_the_earth_engine_convention():
    rows, cols = np.mgrid[0:5, 0:5].astype(float)
    rises_east = cols * 10.0
    rises_north = -rows * 10.0
    _, aspect = sar.terrain_slope_aspect(rises_east, 10.0)
    assert aspect[2, 2] == pytest.approx(270.0)
    slope, aspect = sar.terrain_slope_aspect(rises_north, 10.0)
    assert aspect[2, 2] == pytest.approx(180.0)
    assert slope[2, 2] == pytest.approx(45.0)


def test_layover_faces_the_sensor_and_shadow_faces_away():
    # Descending pass over Nepal: sensor to the east, look direction ~87 deg.
    heading, theta = 87.0, np.array([35.0])
    facing_sensor = np.array([90.0])
    facing_away = np.array([270.0])
    lay, sha = sar.layover_shadow(np.array([40.0]), facing_sensor, theta, heading)
    assert lay[0] and not sha[0]
    lay, sha = sar.layover_shadow(np.array([60.0]), facing_away, theta, heading)
    assert sha[0] and not lay[0]
    lay, sha = sar.layover_shadow(np.array([40.0]), facing_away, theta, heading)
    assert not lay[0] and not sha[0]
    lay, sha = sar.layover_shadow(np.array([20.0]), facing_sensor, theta, heading)
    assert not lay[0] and not sha[0]


def test_orbit_agreement_and_sar_strict():
    yes, no = np.array([True]), np.array([False])
    assert sar.agreement([yes, no], [yes, yes], 0.5)[0] == 1          # 1 of 2
    assert sar.agreement([yes, no, no], [yes, yes, yes], 0.5)[0] == 0  # 1 of 3
    assert sar.agreement([yes, yes, no], [yes, yes, yes], 0.5)[0] == 1  # 2 of 3
    assert sar.agreement([yes, no], [yes, no], 0.5)[0] == 1           # 1 of 1 valid
    assert sar.agreement([yes, yes, no], [yes, yes, yes], 0.5, strict=True)[0] == 0
    assert sar.agreement([yes, yes], [yes, yes], 0.5, strict=True)[0] == 1
    assert np.isnan(sar.agreement([yes], [no], 0.5)[0])


def test_robust_z_uses_the_reference_zone():
    delta = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 10.0])
    reference = np.array([True, True, True, True, True, False])
    z, stats = sar.robust_z(delta, reference)
    assert stats["median_db"] == 2.0 and stats["mad_db"] == 1.0
    assert z[-1] == pytest.approx(8.0 / 1.4826)


def test_focal_median_spreads_nodata_and_smooths_speckle():
    image = np.zeros((5, 5))
    image[2, 2] = 9.0
    assert sar.focal_median(image, 3)[2, 2] == 0.0
    image[0, 0] = np.nan
    out = sar.focal_median(image, 3)
    assert np.isnan(out[1, 1]) and np.isfinite(out[3, 3])


def test_combined_not_impacted_needs_a_valid_optical_zero():
    nan = np.nan
    o = np.array([1.0, 0.0, nan, nan, nan, 0.0])
    s = np.array([nan, nan, 1.0, 0.0, nan, 1.0])
    out = combine(o, s)
    assert out[0] == 1 and out[1] == 0 and out[2] == 1 and out[5] == 1
    assert np.isnan(out[3])     # optical nodata, radar 0: not a valid "no impact"
    assert np.isnan(out[4])


def test_dark_pixels_are_optical_nodata_except_in_the_channel():
    b8_pre = np.array([0.05, 0.30, 0.30, 0.05])
    b8_post = np.array([0.30, 0.05, 0.30, 0.05])
    channel = np.array([False, False, False, True])
    dark = optical.dark_mask(b8_pre, b8_post, channel, 0.08)
    assert list(dark) == [True, True, False, False]
    out = optical.optical_impact(np.full(4, 0.6), np.full(4, 0.1), np.full(4, 0.3),
                                 channel, 0.1, 0.30, 0.30, dark=dark)
    assert np.isnan(out[0]) and np.isnan(out[1])
    assert out[2] == 1 and out[3] == 0


def test_layover_buffer_grows_the_mask():
    slope = np.zeros((1, 9))
    slope[0, 4] = 60.0
    aspect = np.full((1, 9), 90.0)
    theta = np.full((1, 9), 35.0)
    lay, _ = sar.layover_shadow(slope, aspect, theta, 87.0)
    assert lay.sum() == 1
    lay, _ = sar.layover_shadow(slope, aspect, theta, 87.0, buffer_m=30.0, pixel_m=10.0)
    assert lay.sum() == 7


def test_manual_exclusions_are_validated():
    from types import SimpleNamespace as NS
    from floodforest.label import manual_exclusions
    good = NS(reach_start_km=0, reach_end_km=35.5, manual_exclusions=[
        NS(from_km=3.7, to_km=4.7, side="right", reason="hillslope_failure")])
    assert manual_exclusions(good, "x")[0]["side"] == "right"
    for bad in (NS(from_km=3.7, to_km=4.7, side="west", reason="tributary_fan"),
                NS(from_km=3.7, to_km=4.7, side="left", reason="landslide"),
                NS(from_km=4.7, to_km=3.7, side="left", reason="tributary_fan"),
                NS(from_km=35.0, to_km=36.0, side="left", reason="tributary_fan")):
        with pytest.raises(RuntimeError):
            manual_exclusions(NS(reach_start_km=0, reach_end_km=35.5,
                                 manual_exclusions=[bad]), "x")
    with pytest.raises(RuntimeError):
        manual_exclusions(NS(reach_start_km=0, reach_end_km=35.5), "x")


def test_exclusion_bar_sits_on_the_named_side():
    from floodforest.label import exclusion_bar
    samples = np.column_stack([np.zeros(101), -np.arange(101) * 10.0])  # flows south
    river_m = np.arange(101) * 10.0
    item = {"from_km": 0.2, "to_km": 0.5, "side": "right"}
    bar = exclusion_bar(samples, river_m, item, 40.0, 110.0)
    assert (bar[:, 0] < -39).all()       # flowing south, right bank is west
    item["side"] = "left"
    assert (exclusion_bar(samples, river_m, item, 40.0, 110.0)[:, 0] > 39).all()


def test_winsorise_caps_only_the_top():
    from floodforest.stats import winsorise_upper
    values = np.append(np.arange(1.0, 101.0), np.nan)
    out = winsorise_upper(values, 97.5)
    cap = np.percentile(np.arange(1.0, 101.0), 97.5)
    assert out[99] == pytest.approx(cap) and out[0] == 1.0 and np.isnan(out[-1])


def test_connectivity_keeps_patches_joined_to_channel_or_seed():
    label = np.zeros((7, 9))
    anchor = np.zeros((7, 9), dtype=bool)
    anchor[:, 0] = True
    label[3, 1:3] = 1          # touches the anchor
    label[3, 4:6] = 1          # one-pixel gap at column 3: bridged
    label[0, 8] = 1            # isolated
    out = connectivity_filter(label, anchor, 1)
    assert out[3, 1] == 1 and out[3, 5] == 1
    assert out[0, 8] == 0
    assert out[3, 3] == 0      # the bridge itself is not added


def test_label_source_codes():
    combined = np.array([1.0, 1.0, 1.0, 0.0, np.nan])
    o = np.array([1.0, 0.0, 1.0, 0.0, np.nan])
    s = np.array([0.0, 1.0, 1.0, 0.0, np.nan])
    assert list(label_source(combined, o, s)) == [1, 2, 3, 0, 255]


def _profile(channel_offsets):
    offsets = np.arange(-30, 31) * 10.0
    ch = np.isin(offsets, channel_offsets)
    z = 1000.0 + np.abs(offsets) / 10.0
    return ch, offsets, z


def test_centring_uses_the_run_at_the_station():
    ch, off, z = _profile([-20, -10, 0, 10])
    r = centre_channel(ch, off, z, 10.0, 20.0, 100.0)
    assert r["channel_resolved"] == 1
    assert r["bank_top_left_m"] == 10 and r["bank_top_right_m"] == -20
    assert r["channel_width_m"] == 40
    assert r["centre_offset_m"] == -5
    assert r["z0"] == pytest.approx(np.median([1002, 1001, 1000, 1001]))


def test_centring_recentres_within_100_m_and_merges_short_gaps():
    ch, off, z = _profile([60, 70, 110, 120])      # gap 80-100 m: 3 samples, 20 m
    r = centre_channel(ch, off, z, 10.0, 20.0, 100.0)
    assert r["channel_resolved"] == 1
    assert r["bank_top_right_m"] == 60 and r["bank_top_left_m"] == 120
    ch, off, z = _profile([60, 70, 120, 130])      # gap 80-110 m: 4 samples, 30 m
    r = centre_channel(ch, off, z, 10.0, 20.0, 100.0)
    assert r["bank_top_left_m"] == 70


def test_centring_falls_back_to_one_pixel_when_no_channel_is_close():
    ch, off, z = _profile([150, 160])
    r = centre_channel(ch, off, z, 10.0, 20.0, 100.0)
    assert r["channel_resolved"] == 0
    assert r["bank_top_left_m"] == 5 and r["bank_top_right_m"] == -5
    assert r["channel_width_m"] == 10
    assert np.isnan(r["centre_offset_m"])
    assert r["z0"] == 1000.0
