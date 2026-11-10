"""The placebo (8.9, 14.2): a derived event run through the unchanged pipeline,
a confirmed flood-free period, a placebo-only WorldCover stand-in, and a report
whose diagnostic categories never decide what is impact."""

import copy
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from floodforest.config import derive_placebo_event, load_config, shift_years
from floodforest.placebo import (attach_visual_review, classify_detection, detection_table,
                                 forest_terciles, outcome_summary, resolve_placebo, standin_flag,
                                 tercile_table)

DIAG = SimpleNamespace(min_patch_px=5, hillslope_height_m=20, channel_margin_distance_m=30,
                       channel_margin_height_m=10, channel_search_m=300)


# --- the derived event -----------------------------------------------------

def test_melamchi_placebo_event_uses_placebo_windows_and_the_real_reach():
    cfg = load_config()
    real, p = cfg.events.melamchi2021, cfg.events.melamchi2021_placebo
    assert p.base_event == "melamchi2021" and p.is_placebo and p.year_shift == -2
    assert (p.s2_pre.start, p.s2_pre.end) == ("2018-10-01", "2018-12-15")
    assert (p.s2_post.start, p.s2_post.end) == ("2019-10-01", "2019-12-15")
    assert p.s1_pre.start == "2018-08-03T00:00:00Z" and p.s1_post.start == "2019-08-03T00:00:00Z"
    assert p.pre_end_utc == "2019-06-15T00:00:00Z"
    assert p.post_start_utc == "2019-08-03T00:00:00Z"
    assert (p.sar_near.start, p.sar_near.end) == ("2019-04-01T00:00:00Z", "2019-06-15T00:00:00Z")
    assert p.event_dates == ["2019-06-15"]
    assert p.centerline_file == real.centerline_file
    assert (p.reach_start_km, p.reach_end_km) == (real.reach_start_km, real.reach_end_km)
    assert p.manual_exclusions == real.manual_exclusions
    assert p.exclude_dates == [] and p.known_floods_between == []
    assert not hasattr(p, "placebo") and not hasattr(p, "dynamic_world_gap_check")
    assert p.flood_check.status == "confirmed"
    assert (p.worldcover_standin.year, p.worldcover_standin.for_year) == (2020, 2018)


def test_the_real_event_is_left_unchanged():
    real = load_config().events.melamchi2021
    assert real.s2_pre.start == "2020-10-01" and real.pre_end_utc == "2021-06-15T00:00:00Z"
    assert not getattr(real, "is_placebo", False)
    assert real.placebo.fake_event_date == "2019-06-15"


def test_shift_keeps_month_day_and_time():
    assert shift_years("2021-08-03T00:00:00Z", -2) == "2019-08-03T00:00:00Z"
    assert shift_years("2026-08-26T14:52:00Z", -3) == "2023-08-26T14:52:00Z"
    with pytest.raises(ValueError, match="29 February"):
        shift_years("2020-02-29", -1)


def test_fake_date_must_share_month_and_day():
    real = copy.deepcopy(load_config().events.melamchi2021)
    real.placebo.fake_event_date = "2019-06-20"
    with pytest.raises(ValueError, match="month and day"):
        derive_placebo_event(real, "melamchi2021")


# --- the confirmed flood-free period ---------------------------------------

@pytest.fixture
def unconfirmed_config(tmp_path):
    """A copy of config.yaml with the lende2025 placebo flood_check removed, so the
    refusal paths can be tested without the real (confirmed) config downloading anything."""
    import yaml
    from floodforest.config import PROJECT_ROOT
    raw = yaml.safe_load((PROJECT_ROOT / "config.yaml").read_text(encoding="utf-8"))
    del raw["events"]["lende2025"]["placebo"]["flood_check"]
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    return path


def test_lende_flood_check_is_confirmed_over_the_placebo_interval():
    from floodforest.ee_data import require_confirmed_flood_check
    from floodforest.placebo import flood_check_limitations
    p = load_config().events.lende2025_placebo
    require_confirmed_flood_check(p, "lende2025_placebo")
    fc = p.flood_check
    assert (fc.interval.start, fc.interval.end) == ("2022-07-09", "2023-12-15")
    assert fc.placebo_pair_review.date == "2026-11-10"
    lines = flood_check_limitations(fc, "lende2025_placebo")
    assert lines[0].startswith("- Reported inside the confirmed interval: sudden rise, "
                               "Bhote Koshi, originating in Tibet, 2023-07-12 to 2023-07-13 "
                               "(source: Radio Nepal, 2023-07-13): police warnings issued; "
                               "damaging flood not confirmed.")
    assert len(lines) == 1 + 5
    assert any("could not be seen in radar" in line for line in lines)
    assert any("Wayback imagery cannot separate" in line for line in lines)
    # Melamchi records neither field, so its report gains no lines.
    assert flood_check_limitations(load_config().events.melamchi2021_placebo.flood_check, "m") == []


def test_flood_check_limitations_refuse_bad_entries():
    from floodforest.placebo import flood_check_limitations
    with pytest.raises(RuntimeError, match="missing"):
        flood_check_limitations(SimpleNamespace(reported_floods=[
            SimpleNamespace(start="2023-07-12", end="2023-07-13")]), "k")
    with pytest.raises(RuntimeError, match="empty or non-text"):
        flood_check_limitations(SimpleNamespace(limitations=[" "]), "k")


def test_flood_check_must_be_confirmed_and_cover_the_windows(unconfirmed_config):
    from floodforest.ee_data import require_confirmed_flood_check
    cfg = load_config()
    require_confirmed_flood_check(cfg.events.melamchi2021_placebo, "m")
    with pytest.raises(RuntimeError, match="no confirmed flood_check"):
        require_confirmed_flood_check(
            load_config(unconfirmed_config).events.lende2025_placebo, "l")
    short = copy.deepcopy(cfg.events.melamchi2021_placebo)
    short.flood_check.interval.end = "2019-10-31"
    with pytest.raises(RuntimeError, match="does not cover"):
        require_confirmed_flood_check(short, "m")
    pending = copy.deepcopy(cfg.events.melamchi2021_placebo)
    pending.flood_check.status = "provisional"
    with pytest.raises(RuntimeError, match="no confirmed flood_check"):
        require_confirmed_flood_check(pending, "m")


def test_placebo_step_refuses_without_confirmation_before_any_download(unconfirmed_config):
    from floodforest.placebo import run_placebo
    with pytest.raises(RuntimeError, match="no confirmed flood_check"):
        run_placebo("lende2025", config_path=unconfirmed_config)


def test_placebo_key_must_be_the_real_event():
    cfg = load_config()
    assert resolve_placebo(cfg, "melamchi2021") == "melamchi2021_placebo"
    with pytest.raises(RuntimeError, match="real event key"):
        resolve_placebo(cfg, "melamchi2021_placebo")


def test_fetch_refuses_placebos_of_events_it_cannot_fetch(unconfirmed_config):
    from floodforest.ee_data import fetch_event
    with pytest.raises(RuntimeError, match="not listed"):
        fetch_event("bhotekoshi2026_placebo")
    # lende2025 is listed, but its placebo still needs a confirmed flood check first.
    with pytest.raises(RuntimeError, match="no confirmed flood_check"):
        fetch_event("lende2025_placebo", config_path=unconfirmed_config)


def test_fetch_events_come_from_config():
    from floodforest.ee_data import fetch_event, fetch_events
    cfg = load_config()
    assert fetch_events(cfg) == ("melamchi2021", "lende2025")
    with pytest.raises(RuntimeError, match="not listed"):
        fetch_event("bhotekoshi2026")


def test_worldcover_standin_for_lende_real_and_placebo():
    from floodforest.ee_data import AssetError, worldcover_for_event
    cfg = load_config()
    item, standin = worldcover_for_event(cfg, cfg.events.lende2025, 2024)
    assert item.year == 2021 and standin["for_year"] == 2024
    assert standin["note"].startswith("Stand-in (user decision")
    item, standin = worldcover_for_event(cfg, cfg.events.lende2025_placebo, 2022)
    assert item.year == 2021 and standin["for_year"] == 2022
    assert standin["note"].startswith("Placebo only")
    with pytest.raises(AssetError, match="standin is for 2024"):
        worldcover_for_event(cfg, cfg.events.lende2025, 2025)
    # An event without the setting still refuses a year with no map.
    with pytest.raises(AssetError, match="No other year"):
        worldcover_for_event(cfg, cfg.events.bhotekoshi2026, 2025)


def test_exposure_accepts_a_recorded_event_standin():
    from floodforest.exposure import check_worldcover_year
    hansen = {"pre_composite_year": 2024}
    check_worldcover_year({"hansen": hansen, "worldcover_exported": {
        "year": 2021, "standin_for_year": 2024, "standin_decided_on": "2026-11-03"}})
    with pytest.raises(RuntimeError, match="not the pre-event year"):
        check_worldcover_year({"hansen": hansen, "worldcover_exported": {
            "year": 2021, "standin_for_year": 2024}})


def test_cli_wires_the_placebo_step(unconfirmed_config, monkeypatch):
    from floodforest import placebo
    from floodforest.cli import main
    monkeypatch.setattr(placebo, "load_config", lambda path=None: load_config(unconfirmed_config))
    with pytest.raises(RuntimeError, match="no confirmed flood_check"):
        main(["placebo", "--event", "lende2025"])
    with pytest.raises(SystemExit, match="--refetch applies to placebo"):
        main(["label", "--event", "melamchi2021", "--refetch"])


def test_fetch_state_tells_full_and_partial_downloads_apart(tmp_path, monkeypatch):
    import json
    import floodforest.placebo as placebo
    monkeypatch.setattr(placebo, "PROJECT_ROOT", tmp_path)
    cfg = SimpleNamespace(paths=SimpleNamespace(data="data", outputs="outputs"))
    key = "melamchi2021_placebo"
    assert placebo.fetch_state(cfg, key) == (False, False)
    rasters = tmp_path / "data" / key / "rasters"
    rasters.mkdir(parents=True)
    (rasters / "s2_pre.tif").write_bytes(b"x")
    prov = {"placebo": {}, "files": {"s2_pre": f"data/{key}/rasters/s2_pre.tif"},
            "composites": {}}
    out = tmp_path / "outputs" / key
    out.mkdir(parents=True)
    (out / "provenance.json").write_text(json.dumps(prov))
    assert placebo.fetch_state(cfg, key) == (True, False)
    (rasters / "s1_near.tif").write_bytes(b"x")
    (rasters / "s1_angle.tif").write_bytes(b"x")
    prov["composites"]["s1_near"] = {}
    prov["s1_geometry"] = {}
    (out / "provenance.json").write_text(json.dumps(prov))
    assert placebo.fetch_state(cfg, key) == (True, True)
    (rasters / "s2_pre.tif").unlink()
    assert placebo.fetch_state(cfg, key) == (False, False)


# --- WorldCover stand-in: placebo only ------------------------------------

def test_worldcover_standin_is_used_only_for_a_placebo_that_declares_it():
    from floodforest.ee_data import AssetError, worldcover_for_event
    cfg = load_config()
    item, standin = worldcover_for_event(cfg, cfg.events.melamchi2021_placebo, 2018)
    assert item.year == 2020 and standin["for_year"] == 2018
    item, standin = worldcover_for_event(cfg, cfg.events.melamchi2021, 2020)
    assert item.year == 2020 and standin is None
    with pytest.raises(AssetError, match="No other year"):
        worldcover_for_event(cfg, cfg.events.melamchi2021, 2018)
    no_standin = copy.deepcopy(cfg.events.melamchi2021_placebo)
    no_standin.worldcover_standin = None
    with pytest.raises(AssetError, match="No other year"):
        worldcover_for_event(cfg, no_standin, 2018)
    with pytest.raises(AssetError, match="standin is for 2018"):
        worldcover_for_event(cfg, cfg.events.melamchi2021_placebo, 2019)


def test_exposure_accepts_the_standin_only_in_a_placebo_provenance():
    from floodforest.exposure import check_worldcover_year
    hansen = {"pre_composite_year": 2018}
    check_worldcover_year({"hansen": {"pre_composite_year": 2020},
                           "worldcover_exported": {"year": 2020}})
    with pytest.raises(RuntimeError, match="not the pre-event year"):
        check_worldcover_year({"hansen": hansen, "worldcover_exported": {"year": 2020}})
    stand = {"year": 2020, "placebo_standin_for_year": 2018}
    check_worldcover_year({"hansen": hansen, "worldcover_exported": stand, "placebo": {}})
    with pytest.raises(RuntimeError, match="not the pre-event year"):
        check_worldcover_year({"hansen": hansen, "worldcover_exported": stand})


# --- diagnostic categories --------------------------------------------------

def test_categories_follow_the_configured_rules():
    assert classify_detection(35.0, 5.0, DIAG) == "hillslope"
    assert classify_detection(2.0, 10.0, DIAG) == "channel_margin"
    assert classify_detection(-4.0, 30.0, DIAG) == "channel_margin"
    assert classify_detection(12.0, 10.0, DIAG) == "floodplain_terrace"
    assert classify_detection(3.0, 120.0, DIAG) == "floodplain_terrace"
    assert classify_detection(float("nan"), 10.0, DIAG) == "channel_margin"


def _synthetic_scene():
    from rasterio.transform import from_origin
    n = 40
    transform = from_origin(0.0, 400.0, 10.0, 10.0)
    channel = np.zeros((n, n), bool)
    channel[:, 20] = True                       # river at x = 205, flowing south
    dem = np.full((n, n), 100.0)
    dem[:, :6] = 150.0                          # hillside on the west
    dem[25:35, 28:34] = 105.0                   # low terrace on the east
    impact = np.zeros((n, n))
    source = np.zeros((n, n), np.uint8)
    impact[5:8, 21:23] = 1; source[5:8, 21:23] = 1      # A: channel margin, east
    impact[20:23, 2:4] = 1; source[20:23, 2:4] = 2      # B: hillslope, west
    impact[30:33, 30:32] = 1; source[30:33, 30:32] = 3  # C: terrace, 100 m east
    impact[10, 30:32] = 1                                # D: 2 px, not listed
    samples = np.column_stack([np.full(40, 205.0), 395.0 - 10.0 * np.arange(40)])
    river_m = 10.0 * np.arange(40)
    return impact, source, channel, dem, transform, samples, river_m


def test_detection_table_sorts_synthetic_patches():
    impact, source, channel, dem, transform, samples, river_m = _synthetic_scene()
    det, ids = detection_table(impact, source, channel, dem, transform, 10.0, samples, river_m, DIAG)
    assert len(det) == 3 and ids.max() == 4
    by_cat = det.set_index("category")
    a, b, c = by_cat.loc["channel_margin"], by_cat.loc["hillslope"], by_cat.loc["floodplain_terrace"]
    assert a.side == "left" and a.distance_to_channel_m == 10.0 and a.height_above_channel_m == 0.0
    assert a.optical_share == 1.0 and a.radar_share == 0.0
    assert b.side == "right" and b.height_above_channel_m == 50.0 and b.radar_share == 1.0
    assert c.distance_to_channel_m == 100.0 and c.height_above_channel_m == 5.0
    assert c.optical_share == 1.0 and c.radar_share == 1.0
    assert set(det.area_m2) == {600.0}
    # The categories describe the label; they do not change it.
    assert int((impact == 1).sum()) == 20


def test_visual_review_is_attached_only_within_the_match_distance():
    det = pd.DataFrame({"patch_id": [1, 2], "x": [100.0, 1000.0], "y": [100.0, 1000.0]})
    review = pd.DataFrame({"x": [150.0], "y": [100.0], "visual_class": ["agricultural"],
                           "reviewer": ["user"], "check": ["2018_2019"], "patch": [5],
                           "note": ["fields"]})
    out = attach_visual_review(det, review, 100.0)
    assert list(out.visual_class) == ["agricultural", ""]
    assert out.visual_note.iloc[0] == "2018_2019 patch 5: fields"


# --- tables -----------------------------------------------------------------

def _table():
    return pd.DataFrame({
        "station_id": [0, 0, 1, 1, 2, 2, 3],
        "side": ["left", "right"] * 3 + ["left"],
        "segment_id": [0, 0, 0, 0, 1, 1, 1],
        "analysis_set": [1, 1, 1, 1, 1, 1, 0],
        "exclude_reasons": ["", "", "", "", "", "", "near_bank_nodata"],
        "reach_distance_m_fab": [0.0, 0.0, 30.0, 0.0, 10.0, 0.0, np.nan],
        "reach_height_m_fab": [0.0, 0.0, 2.0, 0.0, 1.0, 0.0, np.nan],
        "overtopped_fab": [0, 0, 1, 0, 0, 0, np.nan],
        "censored_fab": [0, 0, 0, 0, 0, 0, np.nan],
        "forest_frac_100": [0.1, 0.2, 0.3, 0.6, 0.8, 0.9, 0.5],
        "forest_source_worldcover_frac_100": [0.0, 0.0, 0.4, 0.0, 1.0, 0.0, np.nan],
    })


def test_outcome_summary_counts_the_analysis_set():
    s = outcome_summary(_table())
    assert s["analysis_rows"] == 6 and s["analysis_stations"] == 3
    assert s["share_reach_above_0"] == pytest.approx(2 / 6)
    assert s["share_overtopped"] == pytest.approx(1 / 6)
    assert s["mean_reach_distance_m"] == pytest.approx(40 / 6)
    assert s["rows_near_bank_nodata"] == 1


def test_terciles_and_standin_flag():
    t = _table()
    assert list(standin_flag(t, 100)) == [0, 0, 1, 0, 1, 0, 0]
    groups = forest_terciles(t.loc[t.analysis_set == 1, "forest_frac_100"])
    assert list(groups) == ["Low", "Low", "Mid", "Mid", "High", "High"]
    table = tercile_table(t, "forest_frac_100")
    assert list(table.tercile) == ["Low", "Mid", "High"]
    assert list(table.rows) == [2, 2, 2]
    assert list(table.share_overtopped) == [0.0, 0.5, 0.0]
    assert list(table.overtopped_rows) == [0, 1, 0]


def test_river_not_visible_sensitivity_adds_only_those_rows_with_zero_reach():
    from floodforest.placebo import river_not_visible_as_zero
    t = pd.DataFrame({
        "in_reach": [1, 1, 1, 0], "analysis_set": [1, 0, 0, 0],
        "exclude_reasons": ["", "river_not_visible", "river_not_visible,low_valid_fraction",
                            "river_not_visible"],
        "reach_distance_m_fab": [30.0, np.nan, np.nan, np.nan],
        "reach_height_m_fab": [2.0, np.nan, np.nan, np.nan],
        "overtopped_fab": [1, np.nan, np.nan, np.nan], "censored_fab": [0, np.nan, np.nan, np.nan],
    })
    out = river_not_visible_as_zero(t)
    assert list(out.analysis_set) == [1, 1, 0, 0]
    assert out.reach_distance_m_fab.iloc[1] == 0 and out.overtopped_fab.iloc[1] == 0
    assert t.analysis_set.iloc[1] == 0  # the input table is not changed


def test_placebo_labels_come_from_each_events_windows():
    from floodforest.stats_run import placebo_labels
    cfg = load_config()
    assert placebo_labels(cfg.events.melamchi2021) == {
        "placebo_year": 2018, "real_year": 2020, "real_season": "Oct-Dec 2020"}
    assert placebo_labels(cfg.events.lende2025) == {
        "placebo_year": 2022, "real_year": 2024, "real_season": "Oct-Dec 2024"}


def test_placebo_report_years_come_from_each_events_windows():
    from floodforest.placebo import report_years
    cfg = load_config()
    m = report_years(cfg, "melamchi2021", "melamchi2021_placebo")
    assert (m["plac_pre"], m["plac_post"], m["real_pre"], m["real_event"]) == (2018, 2019, 2020, 2021)
    assert (m["plac_post_season"], m["real_pre_season"], m["standin_year"]) == (
        "Oct-Dec 2019", "Oct-Dec 2020", 2020)
    lende = report_years(cfg, "lende2025", "lende2025_placebo")
    assert (lende["plac_pre"], lende["plac_post"], lende["real_pre"], lende["real_event"]) == (
        2022, 2023, 2024, 2025)
    assert (lende["plac_pre_season"], lende["plac_post_season"], lende["standin_year"]) == (
        "Oct-Dec 2022", "Oct-Dec 2023", 2021)
    assert lende["centerline_traced_on"] is None and lende["validation_section"] == "14.13"
