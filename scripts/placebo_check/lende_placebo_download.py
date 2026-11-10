"""One-off visual check of the lende2025 placebo period, km 0-15 (not the placebo pipeline).

S2: the study's composite (ee_data._s2_composite: Cloud Score+ cs_cdf >= config,
per-band median, < clear_obs_min clear observations = nodata) for the placebo
S2 windows of config.yaml. S1: per-orbit median VV and VH (dB) for the placebo
S1 windows. Only orbits whose scenes have valid pixels over the region are used;
every scene is listed with its coverage. Nothing under data/ is written, and
config is only read.

Output: results/placebo_check_lende/ (rasters in rasters/, scene lists as CSV).
"""
import csv
import json
from datetime import datetime, timezone

import ee
import geopandas as gpd
import numpy as np
from shapely.ops import substring

from floodforest.config import PROJECT_ROOT, load_config
from floodforest.dates import window_bounds
from floodforest.ee_data import (_assert_grid, _grid, _initialize, _on_grid,
                                 _s1_collection, _s2_composite, _save)
from floodforest.qa import line_samples

import sys

# Usage: lende_placebo_download.py [km_from km_to subfolder]; default km 0-15 in the base folder.
KM_FROM, KM_TO = (float(sys.argv[1]), float(sys.argv[2])) if len(sys.argv) > 2 else (0.0, 15.0)
SUB = sys.argv[3] if len(sys.argv) > 3 else ""
REGION_PAD_M = 1000   # region = line km range buffered by this (chips use 300 m)
EXTRA_M = 500         # line extends this far past each end so the end chips are whole

cfg = load_config(PROJECT_ROOT / "config.yaml")
real = cfg.events.lende2025
plc = cfg.events.lende2025_placebo   # derived from config by the code; read only
_initialize(cfg.ee_project)

out = PROJECT_ROOT / "results" / "placebo_check_lende" / SUB
rast = out / "rasters"
rast.mkdir(parents=True, exist_ok=True)

line, _, _ = line_samples(real, cfg)   # refuses an unverified line
start_m = max(KM_FROM * 1000 - EXTRA_M, 0.0)
part = substring(line, start_m, KM_TO * 1000 + EXTRA_M)
area = part.buffer(REGION_PAD_M)
grid = _grid(area.bounds, cfg.grid_m)
region = ee.Geometry(gpd.GeoSeries([area], crs=cfg.crs).to_crs("EPSG:4326").iloc[0].__geo_interface__)
proj = ee.Projection(cfg.crs, transform=grid["transform"])
print("grid", grid)

S2 = {"s2_2022": plc.s2_pre, "s2_2023": plc.s2_post}
S1 = {"s1_2022": plc.s1_pre, "s1_2023": plc.s1_post}
manifest = {"purpose": f"visual check of lende2025 placebo period, km {KM_FROM:g}-{KM_TO:g}; "
                       "not a pipeline output",
            "km_from": KM_FROM, "km_to": KM_TO,
            "grid": grid, "crs": cfg.crs, "region": f"centerline km {start_m / 1000}-{KM_TO + EXTRA_M / 1000} "
            f"buffered {REGION_PAD_M} m", "s2": {}, "s1": {}}


def utc(ms):
    return datetime.fromtimestamp(ms / 1000, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------- Sentinel-2 ----------
score = cfg.assets.bands.cloud_score
s2_rows = []
for name, win in S2.items():
    start, end = window_bounds(win.start, win.end)
    col = ee.ImageCollection(cfg.assets.s2).filterBounds(region).filterDate(start, end)
    linked = col.linkCollection(ee.ImageCollection(cfg.assets.cloud_score), [score])

    def clear_share(image):
        clear = image.select(score).gte(cfg.cloud_score_min).unmask(0).clip(region)
        stats = clear.reduceRegion(ee.Reducer.mean(), region, 60, maxPixels=1e9)
        return image.set({"clear": stats.get(score)})

    info = (linked.map(clear_share).toList(500).map(lambda i: ee.Image(i).toDictionary(
        ["system:id", "system:time_start", "MGRS_TILE", "clear"])).getInfo())
    info.sort(key=lambda r: r["system:time_start"])
    ids = [r["system:id"] for r in info]
    if not ids:
        raise RuntimeError(f"No Sentinel-2 scenes for {name} {start} to {end}.")
    for r in info:
        s2_rows.append([name, utc(r["system:time_start"]), r["system:id"], r.get("MGRS_TILE"),
                        "" if r.get("clear") is None else round(r["clear"], 3)])
    composite, count = _s2_composite(cfg, {"scenes": [{"scene_id": i} for i in ids]}, proj)
    _save(composite, rast / f"{name}.tif", grid, cfg.crs, "float32")
    _save(count, rast / f"{name}_clear_count.tif", grid, cfg.crs, "uint16")
    manifest["s2"][name] = {"window_start_utc": start, "window_end_utc_exclusive": end,
                            "cloud_score_min": cfg.cloud_score_min,
                            "clear_obs_min": cfg.clear_obs_min, "scenes": ids}

with open(out / "scenes_s2.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["window", "time_utc", "scene", "tile", "region_clear_share"])
    w.writerows(s2_rows)

# ---------- Sentinel-1 ----------
s1 = _s1_collection(cfg)
pol0 = cfg.assets.s1_polarisations[0]


def coverage(image):
    valid = image.select(pol0).mask().gt(0).unmask(0).clip(region)
    stats = valid.reduceRegion(ee.Reducer.mean(), region, 60, maxPixels=1e9)
    return image.set({"valid_share": stats.get(pol0)})


s1_rows, by_window = [], {}
for name, win in S1.items():
    start, end = window_bounds(win.start, win.end)
    col = s1.filterBounds(region).filterDate(start, end).map(coverage)
    info = (col.toList(500).map(lambda i: ee.Image(i).toDictionary(
        ["system:id", "system:time_start", "relativeOrbitNumber_start",
         "orbitProperties_pass", "valid_share"])).getInfo())
    info.sort(key=lambda r: r["system:time_start"])
    by_window[name] = (start, end, info)
    for r in info:
        s1_rows.append([name, utc(r["system:time_start"]), r["system:id"],
                        r["relativeOrbitNumber_start"], r["orbitProperties_pass"],
                        round(r.get("valid_share") or 0.0, 3)])

with open(out / "scenes_s1.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["window", "time_utc", "scene", "relative_orbit", "pass", "region_valid_share"])
    w.writerows(s1_rows)

# An orbit is drawn if its scenes cover the region in both windows (share > 0.5).
orbits = {}
for name, (_, _, info) in by_window.items():
    for r in info:
        key = (r["relativeOrbitNumber_start"], r["orbitProperties_pass"])
        orbits.setdefault(key, {}).setdefault(name, []).append(r)
usable = sorted(k for k, v in orbits.items()
                if all(any((r.get("valid_share") or 0) > 0.5 for r in v.get(n, [])) for n in S1))
print("S1 orbits seen:", {f"{k[0]} {k[1]}": {n: len(v.get(n, [])) for n in S1}
                          for k, v in orbits.items()})
print("S1 orbits drawn:", usable)

for name, (start, end, info) in by_window.items():
    bands, used = [], {}
    for orbit, direction in usable:
        ids = [r["system:id"] for r in info
               if (r["relativeOrbitNumber_start"], r["orbitProperties_pass"]) == (orbit, direction)
               and (r.get("valid_share") or 0) > 0.5]
        used[f"{orbit}_{direction}"] = ids
        for pol in cfg.assets.s1_polarisations:
            stack = ee.ImageCollection([_on_grid(ee.Image(i).select(pol), proj, "bilinear", region)
                                        for i in ids])
            bands.append(stack.median().rename(f"{pol}_{orbit}{direction[0]}"))
    _save(ee.Image.cat(bands), rast / f"{name}.tif", grid, cfg.crs, "float32")
    manifest["s1"][name] = {"window_start_utc": start, "window_end_utc_exclusive": end,
                            "scenes_by_orbit": used}

for path in rast.glob("*.tif"):
    _assert_grid(path, cfg.crs, cfg.grid_m)
(out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print("done")
