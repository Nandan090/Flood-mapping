"""Earth Engine access for the fetch step.

Asset IDs come from config.yaml. Each one is opened at runtime. If an ID or a
required band is missing, the run stops; nothing is swapped in its place.
"""

from __future__ import annotations

import json
import re
import warnings
from datetime import datetime, timezone
from pathlib import Path

import ee
import geopandas as gpd
import rasterio
from geedim import Initialize, MaskedImage

from .centerline import assert_centerline_verified
from .config import PROJECT_ROOT, load_config, parse_utc
from .dates import (assert_exposure_is_pre_event, assert_outcome_is_post_event,
                    composite_year, in_window, window_bounds)
from .exposure import hansen_loss_window

FETCH_EVENTS = ("melamchi2021",)
_HANSEN_VERSION = re.compile(
    r"global_forest_change_(\d{4})(?:_v(\d+)_(\d+))?$"
)


class AssetError(RuntimeError):
    """An Earth Engine asset or band could not be used as configured."""


def fetch_event(event_key: str, config_path: str | Path | None = None,
                only: list[str] | None = None) -> Path:
    """Download the corridor rasters for one event and write provenance.json."""
    if event_key not in FETCH_EVENTS:
        raise RuntimeError(
            f"fetch is implemented for {', '.join(FETCH_EVENTS)} only. "
            f"{event_key} is a later step, and nothing was downloaded for it."
        )
    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    _require_exclude_dates(event, event_key)
    _initialize(cfg.ee_project)

    corridor, region, grid = _corridor(event, cfg)
    if only:
        return _fetch_only(cfg, event, event_key, region, grid, only)
    print(
        f"Corridor grid: {grid['width']} x {grid['height']} pixels "
        f"at {cfg.grid_m} m in {cfg.crs}."
    )

    verified, hansen, dem_notice = _verify_assets(cfg)
    pre_year = composite_year(event.s2_pre.start, event.s2_pre.end)
    hansen_rule = hansen_loss_window(
        pre_year, hansen["newest_year"], cfg.hansen_loss_from_year)
    print(hansen_rule["message"])
    print(f"Hansen asset: {hansen['id']} (newest year {hansen['newest_year']}).")
    if dem_notice:
        print(dem_notice)

    worldcover = _worldcover_for_year(cfg, pre_year)
    groups = _scene_groups(cfg, event, region)
    _print_scenes(groups)

    out_dir = PROJECT_ROOT / cfg.paths.data / event_key / "rasters"
    out_dir.mkdir(parents=True, exist_ok=True)
    files = _download_all(
        cfg, region, grid, groups, hansen, worldcover, out_dir)

    provenance = {
        "event": event_key,
        "crs": cfg.crs,
        "grid_m": cfg.grid_m,
        "corridor_buffer_m": cfg.corridor_buffer_m,
        "grid": {k: grid[k] for k in ("width", "height", "transform")},
        "assets_verified": verified,
        "copernicus_dem_notice": dem_notice,
        "hansen": {**hansen, **hansen_rule, "pre_composite_year": pre_year},
        "worldcover_exported": {
            "id": worldcover.id,
            "year": worldcover.year,
            "reason": f"map year matches the Sentinel-2 pre-composite year {pre_year}",
        },
        "worldcover_not_exported": [
            {"id": item.id, "year": item.year,
             "reason": "map year is not the pre-composite year, so it was not used"}
            for item in cfg.assets.worldcover if item.year != pre_year
        ],
        "exclude_dates": list(event.exclude_dates),
        "known_floods_between": _known_floods(event, event_key),
        "composites": groups,
        "files": {name: str(path.relative_to(PROJECT_ROOT)).replace("\\", "/")
                  for name, path in files.items()},
    }
    dest = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"Wrote {dest.relative_to(PROJECT_ROOT)}")
    pre = files["s2_pre"]
    print(f"Pre-flood Sentinel-2 composite: {pre.relative_to(PROJECT_ROOT)}")
    return dest


def dynamic_world_tile_dates(event_key: str, config_path: str | Path | None = None) -> list[dict]:
    """List Dynamic World images over the corridor for one Sentinel-2 tile (14.23 a).

    Tile and dates come from the event's dynamic_world_gap_check in config.
    The end date is included. The list goes into provenance.json; nothing
    is downloaded and the composite window is not changed.
    """
    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    check = getattr(event, "dynamic_world_gap_check", None)
    if check is None:
        raise RuntimeError(f"{event_key} has no dynamic_world_gap_check in config.yaml.")
    _initialize(cfg.ee_project)
    region = _corridor(event, cfg)[1]
    start, end = window_bounds(check.start, check.end)
    collection = (ee.ImageCollection(cfg.assets.dynamic_world)
                  .filterBounds(region).filterDate(start, end))
    ids = collection.aggregate_array("system:index").getInfo()
    times = collection.aggregate_array("system:time_start").getInfo()
    found = sorted((t, i) for t, i in zip(times, ids) if i.endswith(f"_{check.tile}"))
    images = [{"scene_id": f"{cfg.assets.dynamic_world}/{i}",
               "time_utc": datetime.fromtimestamp(t / 1000, timezone.utc)
               .strftime("%Y-%m-%dT%H:%M:%SZ")} for t, i in found]
    print(f"Dynamic World images over the corridor for tile {check.tile}, "
          f"{check.start} to {check.end} (inclusive): {len(images)}")
    for item in images:
        print(f"  {item['time_utc']}  {item['scene_id']}")
    prov_path = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    provenance = json.loads(prov_path.read_text(encoding="utf-8"))
    provenance["dynamic_world_gap_check"] = {
        "tile": check.tile, "start": check.start, "end_inclusive": check.end,
        "s2_pre_window": {"start": event.s2_pre.start, "end": event.s2_pre.end},
        "images": images,
    }
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"Updated {prov_path.relative_to(PROJECT_ROOT)}")
    return images


def _initialize(project: str) -> None:
    if not isinstance(project, str) or not project.strip():
        raise RuntimeError(
            "ee_project in config.yaml is empty. Earth Engine was not started."
        )
    try:
        Initialize(project=project)
    except Exception as exc:
        raise RuntimeError(
            "Earth Engine did not start. Run `earthengine authenticate` once "
            "in this environment, and check ee_project in config.yaml. "
            f"Earth Engine said: {exc}"
        ) from exc


def _require_exclude_dates(event, event_key: str) -> None:
    if not hasattr(event, "exclude_dates"):
        raise RuntimeError(
            f"{event_key} has no exclude_dates in config.yaml. Add a list "
            "(use [] when there is nothing to drop). No date was skipped silently."
        )


def _corridor(event, cfg):
    path = PROJECT_ROOT / event.centerline_file
    assert_centerline_verified(path)
    gdf = gpd.read_file(path)
    if "role" not in gdf.columns:
        raise RuntimeError(f"{path} has no 'role' property, so the centerline was not used.")
    lines = gdf[gdf["role"] == "centerline"]
    if len(lines) != 1:
        raise RuntimeError(
            f"{path} has {len(lines)} features with role 'centerline'; expected one."
        )
    projected = lines.to_crs(cfg.crs)
    if projected.crs is None or projected.crs.to_string() != cfg.crs:
        raise RuntimeError(
            f"Centerline did not reproject to {cfg.crs}. "
            f"Got {None if projected.crs is None else projected.crs.to_string()}."
        )
    polygon = projected.geometry.iloc[0].buffer(cfg.corridor_buffer_m)
    grid = _grid(polygon.bounds, cfg.grid_m)
    wgs = gpd.GeoSeries([polygon], crs=cfg.crs).to_crs("EPSG:4326").iloc[0]
    region = ee.Geometry(wgs.__geo_interface__)
    return polygon, region, grid


def _grid(bounds, scale: float) -> dict:
    import math
    minx, miny, maxx, maxy = bounds
    minx = math.floor(minx / scale) * scale
    miny = math.floor(miny / scale) * scale
    maxx = math.ceil(maxx / scale) * scale
    maxy = math.ceil(maxy / scale) * scale
    width = int(round((maxx - minx) / scale))
    height = int(round((maxy - miny) / scale))
    if width < 1 or height < 1 or width * height > 30_000_000:
        raise RuntimeError(
            f"Corridor grid is {width} x {height} pixels. That is not a "
            f"{scale:g} m corridor, so nothing was downloaded."
        )
    return {
        "width": width,
        "height": height,
        "transform": [scale, 0.0, minx, 0.0, -scale, maxy],
    }


def _verify_assets(cfg) -> tuple[list, dict, str | None]:
    assets = cfg.assets
    bands = assets.bands
    verified = []
    verified.append(_verify_collection(assets.s2, bands.s2))
    verified.append(_verify_collection(assets.cloud_score, [bands.cloud_score]))
    verified.append(_verify_collection(assets.dynamic_world, [bands.dynamic_world]))
    verified.append(_verify_collection(assets.fabdem, [bands.fabdem]))
    dem_notice = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        verified.append(_verify_collection(assets.copernicus_dem, [bands.copernicus_dem]))
    notes = [str(item.message).strip() for item in caught
             if issubclass(item.category, DeprecationWarning)]
    if notes:
        dem_notice = (
            f"{assets.copernicus_dem} opened, and Earth Engine reports it is "
            f"deprecated: {notes[0]} The configured asset was kept."
        )
    verified.append(_verify_image(assets.canopy_height, [bands.canopy_height]))
    verified.append(_verify_image(assets.merit, [bands.merit]))
    for item in assets.worldcover:
        verified.append(_verify_collection(item.id, [bands.worldcover]))
    hansen = _newest_hansen(assets.hansen_parent)
    hansen_bands = _verify_image(hansen["id"], [bands.hansen])
    verified.append(hansen_bands)
    # Sentinel-1's first global image can be an HH scene. VV and VH are
    # checked on the filtered window, where those bands are the ones used.
    verified.append(_verify_collection(assets.s1, None))
    print("  VV and VH are checked on the filtered window; "
          "the first scene in the global collection can be HH only.")
    return verified, hansen, dem_notice


def _verify_collection(asset_id: str, required: list[str] | None) -> dict:
    try:
        first = ee.ImageCollection(asset_id).limit(1).first()
        found = first.bandNames().getInfo()
    except Exception as exc:
        raise AssetError(
            f"Earth Engine asset '{asset_id}' could not be opened. "
            f"No replacement was substituted. Earth Engine said: {exc}"
        ) from exc
    return _checked(asset_id, "collection", found, required)


def _verify_image(asset_id: str, required: list[str] | None) -> dict:
    try:
        found = ee.Image(asset_id).bandNames().getInfo()
    except Exception as exc:
        raise AssetError(
            f"Earth Engine asset '{asset_id}' could not be opened. "
            f"No replacement was substituted. Earth Engine said: {exc}"
        ) from exc
    return _checked(asset_id, "image", found, required)


def _checked(asset_id, kind, found, required) -> dict:
    if not found:
        raise AssetError(f"'{asset_id}' opened but reported no bands.")
    if required:
        missing = [band for band in required if band not in found]
        if missing:
            raise AssetError(
                f"'{asset_id}' is missing {missing}. Bands present: {found}. "
                "No other band was substituted."
            )
    print(f"Verified {asset_id}")
    print(f"  bands: {found}")
    return {"id": asset_id, "kind": kind, "bands": found}


def _newest_hansen(parent: str) -> dict:
    try:
        assets = _list_assets(parent)
    except Exception as exc:
        raise AssetError(
            f"Could not list Hansen assets under '{parent}'. "
            f"No version was guessed. Earth Engine said: {exc}"
        ) from exc
    ranked = []
    for asset in assets:
        aid = asset.get("id") or asset.get("name") or ""
        name = aid.rstrip("/").split("/")[-1]
        match = _HANSEN_VERSION.search(name)
        if match:
            year = int(match.group(1))
            major = int(match.group(2) or 0)
            minor = int(match.group(3) or 0)
            ranked.append((year, major, minor, aid))
    if not ranked:
        raise AssetError(
            f"No global_forest_change asset was found under '{parent}'. "
            "No Hansen version was guessed."
        )
    year, _major, _minor, aid = max(ranked)
    return {"id": aid, "newest_year": year, "parent": parent}


def _list_assets(parent: str) -> list:
    found = []
    token = None
    while True:
        args = {"parent": parent}
        if token:
            args["pageToken"] = token
        page = ee.data.listAssets(args)
        found.extend(page.get("assets", []))
        token = page.get("nextPageToken")
        if not token:
            return found


def _worldcover_for_year(cfg, year: int):
    matches = [item for item in cfg.assets.worldcover if item.year == year]
    if len(matches) != 1:
        years = [item.year for item in cfg.assets.worldcover]
        raise AssetError(
            f"WorldCover has {len(matches)} entries for pre-composite year "
            f"{year} (configured years: {years}). No other year was used."
        )
    return matches[0]


def _scene_groups(cfg, event, region) -> dict:
    assets = cfg.assets
    groups = {}
    groups["s2_pre"] = _collect(
        ee.ImageCollection(assets.s2), event.s2_pre, region, event,
        role="exposure", extra=["MGRS_TILE"])
    groups["s2_post"] = _collect(
        ee.ImageCollection(assets.s2), event.s2_post, region, event,
        role="outcome", extra=["MGRS_TILE"])
    groups["dynamic_world_pre"] = _collect(
        ee.ImageCollection(assets.dynamic_world), event.s2_pre, region, event,
        role="exposure", extra=["MGRS_TILE"])
    s1 = _s1_collection(cfg)
    groups["s1_pre"] = _drop_short_orbits(_collect(
        s1, event.s1_pre, region, event, role="exposure",
        extra=["relativeOrbitNumber_start", "orbitProperties_pass"]),
        cfg.sar_scenes_min, "s1_pre")
    groups["s1_post"] = _collect(
        s1, event.s1_post, region, event, role="outcome",
        extra=["relativeOrbitNumber_start", "orbitProperties_pass"])
    _require_s1_bands(groups["s1_pre"], assets)
    _require_s1_bands(groups["s1_post"], assets)
    return groups


def _collect(collection, window, region, event, role: str, extra: list[str]) -> dict:
    start, end = window_bounds(window.start, window.end)
    scoped = collection.filterBounds(region).filterDate(start, end)
    before = scoped.size().getInfo()
    scoped = _exclude_dates(scoped, list(event.exclude_dates))
    after = scoped.size().getInfo()
    dropped = before - after
    rows = _scene_rows(scoped, extra)
    check = (assert_exposure_is_pre_event if role == "exposure"
             else assert_outcome_is_post_event)
    cutoff = event.pre_end_utc if role == "exposure" else event.post_start_utc
    scenes = []
    for row in sorted(rows, key=lambda item: item["time_ms"]):
        when = _utc(row["time_ms"])
        if not in_window(parse_utc(when), start, end):
            raise RuntimeError(
                f"Scene {row['scene_id']} at {when} is outside {start} to {end}."
            )
        check(when, cutoff, what=row["scene_id"])
        scene = {"scene_id": row["scene_id"], "time_utc": when}
        for key in extra:
            if row.get(key) is not None:
                scene[key] = row[key]
        scenes.append(scene)
    if not scenes:
        raise RuntimeError(
            f"No scenes in {start} to {end} after exclude_dates "
            f"{list(event.exclude_dates)}. Nothing was invented to fill the gap."
        )
    times = [scene["time_utc"] for scene in scenes]
    return {
        "role": role,
        "window_start_utc": start,
        "window_end_utc": end,
        "excluded_scene_count": dropped,
        "date_range_used": {"start": times[0], "end": times[-1]},
        "scenes": scenes,
    }


def _exclude_dates(collection, days: list):
    from datetime import timedelta
    for day in days:
        if not isinstance(day, str) or len(day) != 10:
            raise RuntimeError(
                f"exclude_dates entry {day!r} is not a YYYY-MM-DD date."
            )
        nxt = (datetime.fromisoformat(day) + timedelta(days=1)).date().isoformat()
        collection = collection.filter(ee.Filter.date(
            f"{day}T00:00:00", f"{nxt}T00:00:00").Not())
    return collection


def _scene_rows(collection, extra: list[str]) -> list:
    count = collection.size().getInfo()
    if count == 0:
        return []
    if count > 500:
        raise RuntimeError(
            f"The window contains {count} scenes. Refusing to truncate the list."
        )
    fields = ["scene_id", "time_ms", *extra]

    def pack(image):
        values = {
            "scene_id": image.get("system:id"),
            "time_ms": image.get("system:time_start"),
        }
        for name in extra:
            values[name] = image.get(name)
        return image.set(values)

    rows = (collection.map(pack).toList(count)
            .map(lambda image: ee.Image(image).toDictionary(fields))
            .getInfo())
    for row in rows:
        if not row.get("scene_id") or row.get("time_ms") is None:
            raise RuntimeError(f"A scene is missing its id or time: {row}")
    return rows


def _utc(time_ms) -> str:
    moment = datetime.fromtimestamp(float(time_ms) / 1000.0, timezone.utc)
    return moment.strftime("%Y-%m-%dT%H:%M:%SZ")


def _require_s1_bands(group: dict, assets) -> None:
    image = ee.Image(group["scenes"][0]["scene_id"])
    try:
        found = image.bandNames().getInfo()
    except Exception as exc:
        raise AssetError(
            f"Sentinel-1 scene {group['scenes'][0]['scene_id']} could not be "
            f"opened. Earth Engine said: {exc}"
        ) from exc
    missing = [band for band in assets.s1_polarisations if band not in found]
    if missing:
        raise AssetError(
            f"Sentinel-1 scene {group['scenes'][0]['scene_id']} is missing "
            f"{missing}. Bands present: {found}."
        )
    print(f"Sentinel-1 bands used, from {group['scenes'][0]['scene_id']}: {found}")
    print(
        f"  date range: {group['date_range_used']['start']} to "
        f"{group['date_range_used']['end']}"
    )


def _known_floods(event, event_key: str) -> list:
    if not hasattr(event, "known_floods_between"):
        raise RuntimeError(
            f"{event_key} has no known_floods_between in config.yaml. "
            "Add a list (use [] when none are known). None were assumed."
        )
    return list(event.known_floods_between)


def _fetch_only(cfg, event, event_key, region, grid, only) -> Path:
    """Refresh the named products and leave the other rasters as they are."""
    allowed = {"s1_pre", "copernicus_dem", "s1_near", "s1_geometry"}
    unknown = [name for name in only if name not in allowed]
    if unknown:
        raise RuntimeError(
            f"Partial fetch does not know {unknown}. Nothing was downloaded."
        )
    dest = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    if not dest.exists():
        raise RuntimeError(
            f"{dest} is missing, so a partial fetch would leave an incomplete "
            "record. No files were changed."
        )
    provenance = json.loads(dest.read_text(encoding="utf-8"))
    out_dir = PROJECT_ROOT / cfg.paths.data / event_key / "rasters"
    crs = cfg.crs
    proj = ee.Projection(crs, transform=grid["transform"])
    if "copernicus_dem" in only:
        checked, notice = _verify_copernicus(cfg)
        path = _save(
            _mosaic_on_grid(ee.ImageCollection(cfg.assets.copernicus_dem),
                            cfg.assets.bands.copernicus_dem, "elevation_m",
                            proj, "bilinear", region),
            out_dir / "copernicus_dem.tif", grid, crs, "float32")
        _assert_grid(path, crs, cfg.grid_m)
        _replace_verified_asset(provenance, checked)
        provenance["copernicus_dem_notice"] = notice
        print(f"Copernicus DEM: {cfg.assets.copernicus_dem} opened"
              + (f". {notice}" if notice else "."))
    if "s1_pre" in only:
        group = _drop_short_orbits(
            _collect(_s1_collection(cfg), event.s1_pre, region, event,
                     role="exposure",
                     extra=["relativeOrbitNumber_start", "orbitProperties_pass"]),
            cfg.sar_scenes_min, "s1_pre")
        _require_s1_bands(group, cfg.assets)
        path = _save(
            _s1_stack(group, cfg.assets.s1_polarisations, proj, region),
            out_dir / "s1_pre.tif", grid, crs, "float32")
        _assert_grid(path, crs, cfg.grid_m)
        provenance.setdefault("composites", {})["s1_pre"] = group
    if "s1_near" in only:
        if not hasattr(event, "sar_near"):
            raise RuntimeError(f"{event_key} has no sar_near window in config.yaml.")
        group = _drop_short_orbits(
            _collect(_s1_collection(cfg), event.sar_near, region, event,
                     role="exposure",
                     extra=["relativeOrbitNumber_start", "orbitProperties_pass"]),
            cfg.sar_scenes_min, "s1_near")
        _require_s1_bands(group, cfg.assets)
        path = _save(
            _s1_stack(group, cfg.assets.s1_polarisations, proj, region),
            out_dir / "s1_near.tif", grid, crs, "float32")
        _assert_grid(path, crs, cfg.grid_m)
        group["label"] = "sar_near: sensitivity baseline, never the primary"
        provenance.setdefault("composites", {})["s1_near"] = group
    if "s1_geometry" in only:
        post = provenance.get("composites", {}).get("s1_post")
        if not post:
            raise RuntimeError(
                "provenance.json has no s1_post scenes, so the radar geometry "
                "cannot be matched to them. Run the full fetch first.")
        geometry, image = _s1_geometry(post, proj, region)
        path = _save(image, out_dir / "s1_angle.tif", grid, crs, "float32")
        _assert_grid(path, crs, cfg.grid_m)
        provenance["s1_geometry"] = geometry
    provenance["known_floods_between"] = _known_floods(event, event_key)
    dest.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"Updated {dest.relative_to(PROJECT_ROOT)}")
    return dest


def _s1_collection(cfg):
    collection = (ee.ImageCollection(cfg.assets.s1)
                  .filter(ee.Filter.eq("instrumentMode", cfg.assets.s1_mode)))
    for polarisation in cfg.assets.s1_polarisations:
        collection = collection.filter(ee.Filter.listContains(
            "transmitterReceiverPolarisation", polarisation))
    return collection


def _verify_copernicus(cfg) -> tuple[dict, str | None]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        checked = _verify_collection(
            cfg.assets.copernicus_dem, [cfg.assets.bands.copernicus_dem])
    notes = [str(item.message).strip() for item in caught
             if issubclass(item.category, DeprecationWarning)]
    if not notes:
        return checked, None
    return checked, (
        f"{cfg.assets.copernicus_dem} opened, and Earth Engine reports it is "
        f"deprecated: {notes[0]} The configured asset was kept."
    )


def _replace_verified_asset(provenance: dict, checked: dict) -> None:
    verified = provenance.setdefault("assets_verified", [])
    verified[:] = [
        item for item in verified
        if not str(item.get("id", "")).startswith("COPERNICUS/DEM/GLO30")
    ]
    verified.append(checked)


def _drop_short_orbits(group: dict, minimum: int, label: str) -> dict:
    """Keep orbits with enough scenes. Shorter orbits are logged and left out."""
    by_orbit: dict[tuple, list] = {}
    for scene in group["scenes"]:
        key = (scene.get("relativeOrbitNumber_start"),
               scene.get("orbitProperties_pass"))
        by_orbit.setdefault(key, []).append(scene)
    print(f"{label} scenes by orbit:")
    kept = []
    skipped = []
    for key in sorted(by_orbit, key=lambda item: (item[0] is None, item[0] or 0)):
        rows = by_orbit[key]
        orbit, direction = key
        print(f"  orbit {orbit} {direction}: {len(rows)} scenes")
        for scene in rows:
            print(f"    {scene['time_utc']}  {scene['scene_id']}")
        if len(rows) < minimum:
            print(f"    skipped: fewer than {minimum} scenes")
            skipped.append({
                "relativeOrbitNumber_start": orbit,
                "orbitProperties_pass": direction,
                "scene_count": len(rows),
            })
        else:
            kept.extend(rows)
    if not kept:
        raise RuntimeError(
            f"{label} has no orbit with at least {minimum} scenes. "
            "No scenes were added to fill the gap."
        )
    kept = sorted(kept, key=lambda scene: scene["time_utc"])
    updated = dict(group)
    updated["scenes"] = kept
    updated["skipped_orbits"] = skipped
    updated["date_range_used"] = {"start": kept[0]["time_utc"], "end": kept[-1]["time_utc"]}
    return updated


def _print_scenes(groups: dict) -> None:
    print("Scenes used:")
    for name, group in groups.items():
        print(
            f"  {name}: {len(group['scenes'])} scenes, "
            f"{group['date_range_used']['start']} to {group['date_range_used']['end']}"
            + (f", {group['excluded_scene_count']} excluded by date"
               if group["excluded_scene_count"] else "")
        )
        if name.endswith("post"):
            print("  Post-window scenes (check these for a later flood):")
            for scene in group["scenes"]:
                extra = ""
                if "orbitProperties_pass" in scene:
                    extra = (f"  orbit {scene.get('relativeOrbitNumber_start')} "
                             f"{scene['orbitProperties_pass']}")
                elif "MGRS_TILE" in scene:
                    extra = f"  tile {scene['MGRS_TILE']}"
                print(f"    {scene['time_utc']}  {scene['scene_id']}{extra}")


def _download_all(cfg, region, grid, groups, hansen, worldcover, out_dir):
    files = {}
    bands = cfg.assets.bands
    crs = cfg.crs
    proj = ee.Projection(crs, transform=grid["transform"])
    files.update(_download_s2(cfg, region, grid, groups, out_dir, crs, proj))
    files["dynamic_world_trees_pre"] = _save(
        _dynamic_world(cfg, groups["dynamic_world_pre"], proj),
        out_dir / "dynamic_world_trees_pre.tif", grid, crs, "float32")
    files["canopy_height"] = _save(
        _on_grid(ee.Image(cfg.assets.canopy_height).select(bands.canopy_height)
                 .rename("canopy_height_m"), proj, "bilinear", region),
        out_dir / "canopy_height.tif", grid, crs, "float32")
    files["hansen_loss_year"] = _save(
        _on_grid(_hansen_calendar_year(hansen["id"], bands.hansen),
                 proj, "nearest", region),
        out_dir / "hansen_loss_year.tif", grid, crs, "uint16")
    files["worldcover"] = _save(
        _mosaic_on_grid(ee.ImageCollection(worldcover.id), bands.worldcover,
                        "class", proj, "nearest", region),
        out_dir / "worldcover.tif", grid, crs, "uint8")
    files["fabdem"] = _save(
        _mosaic_on_grid(ee.ImageCollection(cfg.assets.fabdem), bands.fabdem,
                        "elevation_m", proj, "bilinear", region),
        out_dir / "fabdem.tif", grid, crs, "float32")
    files["copernicus_dem"] = _save(
        _mosaic_on_grid(ee.ImageCollection(cfg.assets.copernicus_dem),
                        bands.copernicus_dem, "elevation_m", proj, "bilinear",
                        region),
        out_dir / "copernicus_dem.tif", grid, crs, "float32")
    files["merit_upa"] = _save(
        _on_grid(ee.Image(cfg.assets.merit).select(bands.merit).rename("upa"),
                 proj, "bilinear", region),
        out_dir / "merit_upa.tif", grid, crs, "float32")
    files["s1_pre"] = _save(
        _s1_stack(groups["s1_pre"], cfg.assets.s1_polarisations, proj, region),
        out_dir / "s1_pre.tif", grid, crs, "float32")
    files["s1_post"] = _save(
        _s1_stack(groups["s1_post"], cfg.assets.s1_polarisations, proj, region),
        out_dir / "s1_post.tif", grid, crs, "float32")
    for path in files.values():
        _assert_grid(path, cfg.crs, cfg.grid_m)
    _assert_images_have_spatial_detail(files)
    return files


def _download_s2(cfg, region, grid, groups, out_dir, crs: str, proj) -> dict:
    files = {}
    for name in ("s2_pre", "s2_post"):
        composite, count = _s2_composite(cfg, groups[name], proj)
        files[name] = _save(
            composite, out_dir / f"{name}.tif", grid, crs, "float32")
        files[f"{name}_clear_count"] = _save(
            count, out_dir / f"{name}_clear_count.tif", grid, crs, "uint16")
    return files


def _s2_composite(cfg, group: dict, proj):
    bands = list(cfg.assets.bands.s2)
    score = cfg.assets.bands.cloud_score
    images = [ee.Image(scene["scene_id"]) for scene in group["scenes"]]
    linked = ee.ImageCollection(images).linkCollection(
        ee.ImageCollection(cfg.assets.cloud_score), [score])
    sample = linked.limit(1).first().bandNames().getInfo()
    if score not in sample:
        raise AssetError(
            f"Cloud Score+ band '{score}' did not join onto "
            f"{group['scenes'][0]['scene_id']}. Bands present: {sample}."
        )

    def mask(image):
        clear = image.select(score).gte(cfg.cloud_score_min)
        scaled = (image.select(bands).multiply(cfg.s2_reflectance_scale)
                  .updateMask(clear))
        # Reproject each scene before the median. Reprojecting the median
        # only resamples an already-coarse composite.
        return scaled.resample("bilinear").reproject(proj)

    masked = linked.map(mask)
    count = masked.select(bands[0]).reduce(ee.Reducer.count()).rename("clear_count")
    composite = masked.median().updateMask(count.gte(cfg.clear_obs_min))
    return composite, count


def _dynamic_world(cfg, group: dict, proj):
    band = cfg.assets.bands.dynamic_world
    score = cfg.assets.bands.cloud_score
    images = [ee.Image(scene["scene_id"]) for scene in group["scenes"]]
    linked = ee.ImageCollection(images).linkCollection(
        ee.ImageCollection(cfg.assets.cloud_score), [score])
    sample = linked.limit(1).first().bandNames().getInfo()
    if score not in sample:
        raise AssetError(
            f"Cloud Score+ band '{score}' did not join onto Dynamic World "
            f"scene {group['scenes'][0]['scene_id']}. Bands present: {sample}."
        )

    def mask(image):
        trees = image.select(band).updateMask(
            image.select(score).gte(cfg.cloud_score_min))
        return trees.resample("bilinear").reproject(proj)

    return linked.map(mask).mean().rename("trees")


def _hansen_calendar_year(asset_id: str, band: str):
    """Hansen stores the loss year as years since 2000. Export the calendar year.

    0 stays 0 and means no loss.
    """
    raw = ee.Image(asset_id).select(band)
    calendar = raw.where(raw.gt(0), raw.add(2000))
    return calendar.rename("loss_year")


def _on_grid(image, proj, resampling: str, region):
    """Put one image on the analysis grid before it is combined or saved."""
    clipped = image.clip(region)
    if resampling == "bilinear":
        return clipped.resample("bilinear").reproject(proj)
    if resampling != "nearest":
        raise RuntimeError(f"Unknown resampling {resampling!r}.")
    return clipped.reproject(proj)


def _mosaic_on_grid(collection, source_band: str, name: str, proj,
                    resampling: str, region):
    def prepare(image):
        return _on_grid(image.select(source_band), proj, resampling, region)

    return (collection.filterBounds(region).map(prepare).mosaic().rename(name))


def _s1_stack(group: dict, polarisations, proj, region) -> ee.Image:
    bands = []
    for index, scene in enumerate(group["scenes"]):
        image = ee.Image(scene["scene_id"])
        for polarisation in polarisations:
            bands.append(_on_grid(image.select(polarisation), proj, "bilinear",
                                  region).rename(f"{polarisation}_{index}"))
    return ee.Image.cat(bands)


def _s1_geometry(post: dict, proj, region) -> tuple[dict, ee.Image]:
    """Incidence angle on the grid and look direction, per post-event orbit.

    The look direction follows gee_s1_ard (Mullissa et al. 2021; Vollrath et
    al. 2020): the mean ee.Terrain.aspect of the scene's angle band at 1000 m,
    wrapped to (-180, 180]. A scene without a value stops the run instead of
    being set to 0.
    """
    by_orbit: dict[tuple, list] = {}
    for scene in post["scenes"]:
        key = (scene["relativeOrbitNumber_start"], scene["orbitProperties_pass"])
        by_orbit.setdefault(key, []).append(scene)
    bands = []
    orbits = []
    for (orbit, direction), scenes in sorted(by_orbit.items()):
        images = [ee.Image(scene["scene_id"]) for scene in scenes]
        headings = [_scene_heading(image, scene["scene_id"])
                    for image, scene in zip(images, scenes)]
        name = f"angle_{orbit}_{direction}"
        angle = ee.ImageCollection(
            [_on_grid(image.select("angle"), proj, "bilinear", region)
             for image in images]).mean().rename(name)
        bands.append(angle)
        orbits.append({
            "relativeOrbitNumber_start": orbit,
            "orbitProperties_pass": direction,
            "band": name,
            "heading_deg": sum(headings) / len(headings),
            "scene_headings_deg": headings,
            "scenes": [scene["scene_id"] for scene in scenes],
        })
        print(f"  orbit {orbit} {direction}: look direction "
              f"{sum(headings) / len(headings):.2f} deg from {len(scenes)} scenes")
    geometry = {
        "source": "s1_post scenes",
        "method": "gee_s1_ard terrain geometry (Mullissa et al. 2021; "
                  "Vollrath et al. 2020); angle is the mean of the orbit's "
                  "post-event scenes on the analysis grid",
        "orbits": orbits,
    }
    return geometry, ee.Image.cat(bands)


def _scene_heading(image, scene_id: str) -> float:
    value = (ee.Terrain.aspect(image.select("angle"))
             .reduceRegion(ee.Reducer.mean(), image.geometry(), 1000)
             .get("aspect").getInfo())
    if value is None:
        raise RuntimeError(
            f"No look direction could be computed for {scene_id}. "
            "No default heading was used.")
    value = float(value)
    return value - 360.0 if value > 180.0 else value


def _save(image, path: Path, grid: dict, crs: str, dtype: str) -> Path:
    print(f"Downloading {path.name} ({dtype})")
    MaskedImage(image).download(
        path,
        overwrite=True,
        crs=crs,
        crs_transform=grid["transform"],
        shape=(grid["height"], grid["width"]),
        resampling="near",
        dtype=dtype,
    )
    return path


def _assert_grid(path: Path, crs: str, scale: float) -> None:
    with rasterio.open(path) as src:
        got = src.crs.to_string() if src.crs is not None else None
        if got != crs:
            raise RuntimeError(f"{path.name} is {got}, not {crs}.")
        if abs(src.transform.a - scale) > 1e-6 or abs(src.transform.e + scale) > 1e-6:
            raise RuntimeError(
                f"{path.name} pixel size is {src.transform.a} by {src.transform.e}, "
                f"not {scale} m."
            )


def _assert_images_have_spatial_detail(files: dict) -> None:
    """Fail if a composite was sampled so coarsely that the valley disappeared."""
    import numpy as np
    with rasterio.open(files["s2_pre"]) as src:
        nir = src.read(src.descriptions.index("B8") + 1, masked=True)
    if float(np.ma.std(nir)) < 0.03:
        raise RuntimeError(
            "The pre-flood Sentinel-2 image has almost no spatial detail "
            f"(near-infrared standard deviation {float(np.ma.std(nir)):.4f}). "
            "The 10 m grid was not applied. Nothing was saved as a finished fetch."
        )
    with rasterio.open(files["fabdem"]) as src:
        elev = src.read(1)
    middle = elev[elev.shape[0] // 2]
    if float(middle.max() - middle.min()) < 50:
        raise RuntimeError(
            "FABDEM shows less than 50 m of relief across the corridor, so the "
            "valley was not resolved. The 10 m grid was not applied."
        )
