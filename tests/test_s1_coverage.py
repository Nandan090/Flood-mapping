"""Sentinel-1 scenes must cover the corridor, and zero-filled scenes are never kept."""

from types import SimpleNamespace

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from floodforest import ee_data
from floodforest.ee_data import apply_s1_coverage, s1_band_check


def _scene(i, orbit, direction="DESCENDING"):
    return {"scene_id": f"S1/{i}", "time_utc": f"2025-07-{10 + i:02d}T00:00:00Z",
            "relativeOrbitNumber_start": orbit, "orbitProperties_pass": direction}


def _group(*orbits):
    scenes = [_scene(i, o) for i, o in enumerate(orbits)]
    return {"scenes": scenes, "date_range_used": {}}


def test_scenes_below_the_coverage_share_are_dropped_and_recorded():
    group = _group(19, 121, 19, 121, 19, 121)
    out = apply_s1_coverage(group, [1.0, 0.004, 1.0, 0.003, 0.97, 0.005], 0.95, 3, "s1_post")
    assert [s["relativeOrbitNumber_start"] for s in out["scenes"]] == [19, 19, 19]
    assert [s["corridor_coverage_share"] for s in out["scenes"]] == [1.0, 1.0, 0.97]
    dropped = out["dropped_scenes"]
    assert [d["scene_id"] for d in dropped] == ["S1/1", "S1/3", "S1/5"]
    assert dropped[0]["relativeOrbitNumber_start"] == 121
    assert dropped[0]["time_utc"] == "2025-07-11T00:00:00Z"
    assert dropped[0]["corridor_coverage_share"] == 0.004
    assert "below the 95% minimum" in dropped[0]["reason"]
    assert out["date_range_used"] == {"start": "2025-07-10T00:00:00Z",
                                      "end": "2025-07-14T00:00:00Z"}


def test_exactly_the_threshold_is_kept():
    out = apply_s1_coverage(_group(19, 19, 19), [0.95, 0.95, 1.0], 0.95, 3, "g")
    assert len(out["scenes"]) == 3 and out["dropped_scenes"] == []


def test_an_orbit_pushed_below_the_minimum_stops_the_run():
    with pytest.raises(RuntimeError, match="orbit 85 DESCENDING: 2 of 3 scenes kept"):
        apply_s1_coverage(_group(19, 19, 19, 85, 85, 85), [1, 1, 1, 1, 1, 0.5], 0.95, 3, "s1_pre")


def test_an_orbit_already_short_is_left_to_the_short_orbit_rule():
    # 2 scenes before the filter, 1 after: the existing 5.3 rule skips it later.
    out = apply_s1_coverage(_group(19, 19, 19, 85, 85), [1, 1, 1, 1, 0.1], 0.95, 3, "s1_pre")
    assert [s["relativeOrbitNumber_start"] for s in out["scenes"]] == [19, 19, 19, 85]


def test_nothing_left_stops_the_run():
    with pytest.raises(RuntimeError, match="no Sentinel-1 scene is left"):
        apply_s1_coverage(_group(121, 121), [0.0, 0.01], 0.95, 3, "s1_post")


def _write_stack(path, bands: dict, nodata=-np.inf):
    first = next(iter(bands.values()))
    with rasterio.open(path, "w", driver="GTiff", height=first.shape[0], width=first.shape[1],
                       count=len(bands), dtype="float32", crs="EPSG:32645",
                       transform=from_origin(0, 100, 10, 10), nodata=nodata) as dst:
        for i, (name, data) in enumerate(bands.items(), start=1):
            dst.write(data.astype("float32"), i)
            dst.set_band_description(i, name)


def test_band_check_finds_zero_filled_and_empty_scenes(tmp_path):
    good = np.full((10, 10), -12.0)
    zeros = np.zeros((10, 10))
    empty = np.full((10, 10), -np.inf)
    partial = good.copy()
    partial[:, :3] = 0.0
    path = tmp_path / "s1.tif"
    _write_stack(path, {"VV_0": good, "VH_0": good, "VV_1": zeros, "VH_1": zeros,
                        "VV_2": good, "VH_2": empty, "VV_3": partial, "VH_3": good})
    corridor = np.ones((10, 10), bool)
    bad, stats = s1_band_check(path, _group(19, 121, 121, 19), ("VV", "VH"), corridor)
    assert bad == [1, 2]
    assert stats[0] == {"corridor_nodata_share": 0.0, "corridor_zero_share": 0.0}
    assert stats[1]["corridor_zero_share"] == 1.0
    assert stats[2]["corridor_nodata_share"] == 1.0
    assert stats[3] == {"corridor_nodata_share": 0.0, "corridor_zero_share": 0.3}


def test_band_check_looks_only_inside_the_corridor(tmp_path):
    data = np.zeros((10, 10))
    data[0, 0] = -10.0
    path = tmp_path / "s1.tif"
    _write_stack(path, {"VV_0": data, "VH_0": data})
    corridor = np.zeros((10, 10), bool)
    corridor[5:, 5:] = True
    bad, _ = s1_band_check(path, _group(19), ("VV", "VH"), corridor)
    assert bad == [0]


def test_save_s1_drops_zero_filled_scenes_and_downloads_again(tmp_path, monkeypatch):
    """A zero-filled scene is never kept: it is recorded and the stack is rebuilt."""
    calls = []

    def fake_stack(group, pols, proj, region):
        return [s["scene_id"] for s in group["scenes"]]

    def fake_save(ids, path, grid, crs, dtype):
        calls.append(list(ids))
        bands = {}
        for i, sid in enumerate(ids):
            value = np.zeros((10, 10)) if sid == "S1/1" else np.full((10, 10), -11.0)
            bands[f"VV_{i}"] = value
            bands[f"VH_{i}"] = value
        _write_stack(path, bands)
        return path

    monkeypatch.setattr(ee_data, "_s1_stack", fake_stack)
    monkeypatch.setattr(ee_data, "_save", fake_save)
    cfg = SimpleNamespace(assets=SimpleNamespace(s1_polarisations=["VV", "VH"]),
                          crs="EPSG:32645", sar_scenes_min=3)
    group = _group(19, 19, 19, 19)
    path = tmp_path / "s1_post.tif"
    out = ee_data._save_s1(cfg, group, "s1_post", path, {}, None, None, np.ones((10, 10), bool))
    assert calls == [["S1/0", "S1/1", "S1/2", "S1/3"], ["S1/0", "S1/2", "S1/3"]]
    assert [s["scene_id"] for s in out["scenes"]] == ["S1/0", "S1/2", "S1/3"]
    assert out["dropped_scenes"][0]["scene_id"] == "S1/1"
    assert "all zero or nodata" in out["dropped_scenes"][0]["reason"]
    assert all(s["corridor_zero_share"] == 0.0 for s in out["scenes"])
    with rasterio.open(path) as src:
        assert src.descriptions == ("VV_0", "VH_0", "VV_1", "VH_1", "VV_2", "VH_2")


def test_save_s1_stops_when_a_drop_leaves_an_orbit_short(tmp_path, monkeypatch):
    monkeypatch.setattr(ee_data, "_s1_stack", lambda g, p, pr, r: [s["scene_id"] for s in g["scenes"]])

    def fake_save(ids, path, grid, crs, dtype):
        bands = {}
        for i, sid in enumerate(ids):
            value = np.zeros((4, 4)) if sid == "S1/1" else np.full((4, 4), -11.0)
            bands[f"VV_{i}"] = bands[f"VH_{i}"] = value
        _write_stack(path, bands)
        return path

    monkeypatch.setattr(ee_data, "_save", fake_save)
    cfg = SimpleNamespace(assets=SimpleNamespace(s1_polarisations=["VV", "VH"]),
                          crs="EPSG:32645", sar_scenes_min=3)
    path = tmp_path / "s1_pre.tif"
    with pytest.raises(RuntimeError, match="2 of 3 scenes kept"):
        ee_data._save_s1(cfg, _group(19, 19, 19), "s1_pre", path, {}, None, None,
                         np.ones((4, 4), bool))
    assert not path.exists()  # the zero-filled file is not left behind


def test_config_sets_the_coverage_share():
    from floodforest.config import load_config
    assert load_config().s1_corridor_coverage_min == 0.95
