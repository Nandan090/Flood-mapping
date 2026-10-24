"""The config carries every window from section 4.1 and the section 14 values."""

from datetime import datetime, timedelta, timezone

from floodforest.config import load_config, parse_utc

UTC = timezone.utc


def test_windows_and_event_times():
    cfg = load_config()
    assert cfg.crs == "EPSG:32645"
    assert cfg.seed == 42

    m = cfg.events.melamchi2021
    assert m.s2_pre.start == "2020-10-01" and m.s2_post.end == "2021-12-15"
    assert m.pre_end_utc == "2021-06-15T00:00:00Z"
    assert m.post_start_utc == "2021-08-03T00:00:00Z"
    assert m.placebo.fake_event_date == "2019-06-15"
    assert m.placebo.s2_pre.start == "2018-10-01"

    l = cfg.events.lende2025
    assert l.post_start_utc == "2025-07-09T00:00:00Z"
    assert l.reach_end_km == 35
    assert l.s1_post.start == "2025-07-09T00:00:00Z"
    assert l.placebo.fake_event_date == "2023-07-08"

    b = cfg.events.bhotekoshi2026
    assert b.event_time_utc == "2026-08-26T02:52:00Z"
    assert b.reach_end_km == 72
    # Post imagery starts 12 h after the flood.
    flood = parse_utc(b.event_time_utc)
    assert parse_utc(b.post_start_utc) - flood == timedelta(hours=12)
    assert b.s1_post.start == "2026-08-26T14:52:00Z"
    assert b.placebo.s2_post.start == "2023-10-01"
    assert b.placebo.s1_post.start == "2023-08-26T14:52:00Z"


def test_asset_ids_match_section_4():
    assets = load_config().assets
    assert assets.s1 == "COPERNICUS/S1_GRD"
    assert assets.s2 == "COPERNICUS/S2_SR_HARMONIZED"
    assert assets.cloud_score == "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED"
    assert assets.fabdem == "projects/sat-io/open-datasets/FABDEM"
    assert assets.copernicus_dem == "COPERNICUS/DEM/GLO30_2024_1"
    assert assets.dynamic_world == "GOOGLE/DYNAMICWORLD/V1"
    assert assets.canopy_height == "users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1"
    assert assets.hansen_parent == "UMD/hansen"
    assert assets.merit == "MERIT/Hydro/v1_0_1"
    assert [(item.id, item.year) for item in assets.worldcover] == [
        ("ESA/WorldCover/v100", 2020),
        ("ESA/WorldCover/v200", 2021),
    ]
    assert load_config().events.melamchi2021.exclude_dates == []


def test_s1_baseline_is_the_same_season_one_year_earlier():
    events = load_config().events
    m = events.melamchi2021
    assert (m.s1_pre.start, m.s1_pre.end) == (
        "2020-08-03T00:00:00Z", "2020-09-15T00:00:00Z")
    assert (m.sar_near.start, m.sar_near.end) == (
        "2021-04-01T00:00:00Z", "2021-06-15T00:00:00Z")
    assert m.known_floods_between == []
    assert m.placebo.s1_pre.start == "2018-08-03T00:00:00Z"

    l = events.lende2025
    assert (l.s1_pre.start, l.s1_pre.end) == (
        "2024-07-09T00:00:00Z", "2024-08-15T00:00:00Z")
    assert l.sar_near.start == "2025-05-01T00:00:00Z"
    assert l.known_floods_between == []

    b = events.bhotekoshi2026
    assert (b.s1_pre.start, b.s1_pre.end) == (
        "2025-08-26T00:00:00Z", "2025-09-30T00:00:00Z")
    assert b.sar_near.end == "2026-08-26T00:00:00Z"
    assert b.known_floods_between == []
    assert b.placebo.s1_pre.start == "2022-08-26T14:52:00Z"


def test_decision_rule_thresholds():
    cfg = load_config()
    assert cfg.narrow_effect_sd == 0.25
    assert (cfg.narrow_or_low, cfg.narrow_or_high) == (0.67, 1.50)
    assert cfg.model_b_min_stations == 30
    assert cfg.model_b_min_forest_diff == 0.2
    assert cfg.reach_end_shift_km == 5
    assert cfg.reference_pixels_min == 5000
    assert cfg.hansen_loss_from_year == 2021
    assert cfg.validation_sample_size == 40


def test_uncertainty_settings_from_14_15():
    cfg = load_config()
    assert cfg.bootstrap_block_km == 3
    assert cfg.bootstrap_final_block_min_km == 1.5
    assert cfg.min_blocks == 8
    assert cfg.wild_bootstrap_draws == 1999
    assert cfg.ci_bisection_tolerance_se == 0.01
    assert cfg.bootstrap_reps == 2000   # sensitivity percentile bootstrap


def test_reaches_and_sections_from_14_16():
    cfg = load_config()
    events = cfg.events
    m = events.melamchi2021
    assert (m.reach_start_km, m.reach_end_km) == (0, 35.5486)
    assert m.section_breaks_km == [31.221]
    bazaar = next(mark for mark in m.landmarks if mark.name == "Melamchi Bazaar")
    assert (round(bazaar.lon, 4), round(bazaar.lat, 4)) == (85.5779, 27.8323)
    assert m.rivers == ["Melamchi River", "Indrawati"]
    assert (events.lende2025.reach_start_km,
            events.lende2025.reach_end_km) == (0, 35)
    assert (events.bhotekoshi2026.reach_start_km,
            events.bhotekoshi2026.reach_end_km) == (0, 72)

    for key in ("melamchi2021", "lende2025", "bhotekoshi2026"):
        event = getattr(events, key)
        assert event.reach_end_km is not None, f"{key}: no open-ended reach"
        assert hasattr(event, "section_breaks_km")
        assert isinstance(event.landmarks, list) and event.landmarks
        for landmark in event.landmarks:
            assert landmark.name

    assert cfg.section_break_upa_rise == 0.5
    assert cfg.section_break_window_m == 300


def test_lende_and_bhotekoshi_share_one_centerline_and_origin():
    cfg = load_config()
    lende = cfg.events.lende2025
    bhote = cfg.events.bhotekoshi2026
    assert lende.centerline_file == bhote.centerline_file
    assert lende.reach_start_km == bhote.reach_start_km == 0
    assert lende.start_point.name == bhote.start_point.name
    assert "Rasuwagadhi" in lende.start_point.name
    assert cfg.events.melamchi2021.centerline_file != lende.centerline_file
