"""OpenStreetMap draft of a river centerline (spec section 6).

The draft is a starting line for the user to correct in QGIS. It is saved with
verified_by_user set to false, and nothing downstream will use it until the
user sets that to true.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import requests
from pyproj import Transformer
from shapely import transform
from shapely.geometry import LineString, Point, mapping
from shapely.ops import linemerge, substring, unary_union

from .config import PROJECT_ROOT, load_config

WGS84 = "EPSG:4326"
ELEVATION_URL = "https://api.open-meteo.com/v1/elevation"


class CenterlineNotVerified(RuntimeError):
    """The centerline has not been checked by the user."""


def _centerline_is_verified(path: str | Path) -> bool:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return any(
        f.get("properties", {}).get("role") == "centerline"
        and f.get("properties", {}).get("verified_by_user") is True
        for f in data.get("features", [])
    )


def assert_centerline_verified(path: str | Path) -> None:
    """Refuse to analyse a centerline the user has not signed off."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Centerline not found: {path}")
    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)
    features = data.get("features", [])
    flags = [
        f.get("properties", {}).get("verified_by_user") is True
        for f in features
        if f.get("properties", {}).get("role") == "centerline"
    ]
    if not flags or not all(flags):
        raise CenterlineNotVerified(
            f"{path} is not verified. Check it in QGIS against the pre-event "
            f"imagery, correct it, and set verified_by_user to true."
        )


def downstream_end_of_flat(distance_m: np.ndarray, elevation_m: np.ndarray,
                           max_slope: float, window_m: float,
                           min_length_m: float, min_elevation_m: float) -> float:
    """Distance along the line where a high gentle reach gives way to a gorge.

    The line must run upstream to downstream. Slope is the drop over the next
    window_m, divided by window_m. The flat is the highest run that stays
    gentler than max_slope for at least min_length_m.
    """
    distance_m = np.asarray(distance_m, dtype=float)
    elevation_m = np.asarray(elevation_m, dtype=float)
    if len(distance_m) < 3:
        raise ValueError("Need at least three samples to find the flat.")
    step = distance_m[1] - distance_m[0]
    window = max(1, int(round(window_m / step)))
    # Trailing slope, so the reported end is the knickpoint itself rather
    # than one window short of it.
    slope = np.full(len(elevation_m), np.nan)
    for i in range(window, len(elevation_m)):
        slope[i] = (elevation_m[i - window] - elevation_m[i]) / (window * step)

    gentle = (slope < max_slope) & (elevation_m >= min_elevation_m)
    runs: list[tuple[int, int]] = []
    start = None
    for i, flag in enumerate(gentle):
        if flag and start is None:
            start = i
        if start is not None and (not flag or i == len(gentle) - 1):
            end = i - 1 if not flag else i
            runs.append((start, end))
            start = None
    long = [r for r in runs if distance_m[r[1]] - distance_m[r[0]] >= min_length_m]
    if not long:
        raise RuntimeError(
            "No gentle reach matched the Bhemathang flat "
            f"(slope < {max_slope}, longer than {min_length_m:.0f} m, "
            f"above {min_elevation_m:.0f} m). Not guessing a start point."
        )
    chosen = max(long, key=lambda r: elevation_m[r[0]])
    return float(distance_m[chosen[1]])


def _named(gdf, needle: str):
    """Features whose name contains the distinctive part of a river name."""
    key = needle.split()[0]
    names = gdf["name"].fillna("") if "name" in gdf.columns else ""
    found = gdf[names.str.contains(key, case=False)]
    if found.empty:
        raise RuntimeError(
            f"OpenStreetMap returned no waterway named like {needle!r} "
            f"in the search box."
        )
    return found


def _to_metres(gdf, to_utm: Transformer) -> list[LineString]:
    lines = []
    for geom in gdf.geometry:
        parts = [geom] if geom.geom_type == "LineString" else list(geom.geoms)
        for part in parts:
            if part.is_empty or len(part.coords) < 2:
                continue
            lines.append(LineString(
                [to_utm.transform(x, y) for x, y in part.coords]))
    return lines


def _as_one_line(lines: list[LineString], river: str) -> LineString:
    merged = linemerge(unary_union(lines))
    if merged.geom_type == "MultiLineString":
        gaps = ", ".join(f"{part.length:.0f} m" for part in merged.geoms)
        raise RuntimeError(
            f"The OSM {river} line has gaps and did not join into one line "
            f"(pieces of {gaps}). Not bridging them automatically."
        )
    return merged


def _elevations(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Ground elevation from the Open-Meteo Copernicus DEM service."""
    out: list[float] = []
    for i in range(0, len(lon), 100):
        response = requests.get(
            ELEVATION_URL,
            params={
                "latitude": ",".join(f"{v:.5f}" for v in lat[i:i + 100]),
                "longitude": ",".join(f"{v:.5f}" for v in lon[i:i + 100]),
            },
            timeout=60,
        )
        if response.status_code != 200:
            raise RuntimeError(
                f"Elevation request failed ({response.status_code}) from "
                f"{ELEVATION_URL}. The start point was not guessed."
            )
        block = response.json().get("elevation")
        if not block or any(v is None for v in block):
            raise RuntimeError(
                "The elevation service returned no value for part of the "
                "river. The start point was not guessed."
            )
        out.extend(block)
    return np.asarray(out, dtype=float)


def _sample(line: LineString, step_m: float, to_ll: Transformer):
    count = int(line.length // step_m)
    points = [line.interpolate(i * step_m) for i in range(count + 1)]
    lon, lat = [], []
    for point in points:
        x, y = to_ll.transform(point.x, point.y)
        lon.append(x)
        lat.append(y)
    return points, np.asarray(lon), np.asarray(lat)


def _orient_downstream(line: LineString, to_ll: Transformer) -> LineString:
    start = _elevations(*_ends_lonlat(line, to_ll, which="start"))
    end = _elevations(*_ends_lonlat(line, to_ll, which="end"))
    if start[0] < end[0]:
        return LineString(list(line.coords)[::-1])
    return line


def _ends_lonlat(line, to_ll, which):
    point = Point(line.coords[0] if which == "start" else line.coords[-1])
    lon, lat = to_ll.transform(point.x, point.y)
    return np.array([lon]), np.array([lat])


def _downstream_side(host: LineString, end: Point, snap_m: float,
                     to_ll: Transformer) -> LineString:
    """The part of host that leaves `end` heading downhill."""
    cut = host.project(end)
    candidates = []
    for piece in (substring(host, cut, host.length), substring(host, 0, cut)):
        if piece.length <= snap_m:
            continue
        if Point(piece.coords[0]).distance(end) > Point(piece.coords[-1]).distance(end):
            piece = LineString(list(piece.coords)[::-1])
        probe = piece.interpolate(min(1000, piece.length))
        lon, lat = to_ll.transform(probe.x, probe.y)
        height = _elevations(np.array([lon]), np.array([lat]))[0]
        candidates.append((height, piece))
    if not candidates:
        raise RuntimeError("The joining river has no downstream continuation.")
    return min(candidates, key=lambda item: item[0])[1]


def _extend_downstream(stem: LineString, others: list[LineString],
                       snap_m: float, to_ll: Transformer,
                       target_length_m: float) -> tuple[LineString, float]:
    """Continue stem downhill until it is long enough. Returns the junction distance."""
    junction_m = stem.length
    line = stem
    used: set[int] = set()
    while line.length < target_length_m:
        end = Point(line.coords[-1])
        remaining = [other for other in others if id(other) not in used]
        if not remaining:
            break
        host = min(remaining, key=lambda other: other.distance(end))
        gap = host.distance(end)
        if gap > snap_m:
            break
        used.add(id(host))
        onward = _downstream_side(host, end, snap_m, to_ll)
        line = LineString(list(line.coords) + list(onward.coords)[1:])
    return line, junction_m


def draft_centerline(event_key: str, config_path: str | Path | None = None) -> Path:
    """Fetch the OSM waterways and write an unverified centerline GeoJSON."""
    import osmnx as ox

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    out = PROJECT_ROOT / event.centerline_file
    if out.exists() and _centerline_is_verified(out):
        raise RuntimeError(
            f"{out} is already verified. The draft was not written, so the "
            "checked line is unchanged."
        )
    if not hasattr(event, "osm_bbox"):
        raise RuntimeError(
            f"{event_key} has no osm_bbox in config.yaml, so no draft was fetched."
        )
    start = event.start_point
    to_utm = Transformer.from_crs(WGS84, cfg.crs, always_xy=True)
    to_ll = Transformer.from_crs(cfg.crs, WGS84, always_xy=True)

    ox.settings.requests_timeout = 180
    waterways = ox.features.features_from_bbox(
        tuple(event.osm_bbox), tags={"waterway": ["river", "stream"]})
    print(f"OpenStreetMap waterways in the search box: {len(waterways)}")

    main = _as_one_line(_to_metres(_named(waterways, event.rivers[0]), to_utm),
                        event.rivers[0])
    main = _orient_downstream(main, to_ll)
    points, lon, lat = _sample(main, start.sample_m, to_ll)
    elevation = _elevations(lon, lat)
    distance = np.arange(len(elevation)) * start.sample_m
    flat_end_m = downstream_end_of_flat(
        distance, elevation, start.flat_max_slope, start.flat_slope_window_m,
        start.flat_min_length_m, start.flat_min_elevation_m)
    start_m = flat_end_m + start.downstream_of_flat_m
    if start_m >= main.length:
        raise RuntimeError(
            "The point 500 m downstream of the flat is past the end of the "
            "OSM river line."
        )

    length_m = (event.reach_end_km - event.reach_start_km) * 1000
    end_m = start_m + length_m
    joined, junction_m = _extend_downstream(
        main, _to_metres(_named(waterways, event.rivers[1]), to_utm),
        event.osm_snap_m, to_ll, end_m)
    if end_m > joined.length:
        raise RuntimeError(
            f"The joined OSM line is {joined.length / 1000:.1f} km, shorter "
            f"than the {length_m / 1000:.0f} km reach. Not inventing the rest."
        )
    reach = substring(joined, start_m, end_m)
    flat_pt = main.interpolate(flat_end_m)
    start_pt = joined.interpolate(start_m)
    junction_pt = joined.interpolate(junction_m)

    out.parent.mkdir(parents=True, exist_ok=True)
    _write_geojson(out, to_ll, reach, flat_pt, start_pt, junction_pt, {
        "event": event_key,
        "verified_by_user": False,
        "source": "OpenStreetMap, via osmnx; elevation from the Open-Meteo "
                  "Copernicus DEM service, used only to find the flat",
        "rivers": list(event.rivers),
        "start_rule": start.name,
        "flat_end_m_along_osm_line": round(flat_end_m, 1),
        "flat_end_elevation_m": round(float(np.interp(flat_end_m, distance, elevation)), 1),
        "start_elevation_m": round(float(np.interp(start_m, distance, elevation)), 1),
        "confluence_river_km": round((junction_m - start_m) / 1000, 3),
        "length_m": round(reach.length, 1),
        "note": "Draft only. Check against the pre-event Sentinel-2 image "
                "before use. verified_by_user is false.",
    })
    start_ll = to_ll.transform(start_pt.x, start_pt.y)
    flat_ll = to_ll.transform(flat_pt.x, flat_pt.y)
    join_ll = to_ll.transform(junction_pt.x, junction_pt.y)
    print(f"Draft written to {out}")
    print("verified_by_user: false")
    print(f"Downstream end of the flat: {flat_ll[1]:.5f} N, {flat_ll[0]:.5f} E, "
          f"{np.interp(flat_end_m, distance, elevation):.0f} m")
    print(f"Start, 500 m downstream:    {start_ll[1]:.5f} N, {start_ll[0]:.5f} E, "
          f"{np.interp(start_m, distance, elevation):.0f} m")
    print(f"Melamchi-Indrawati confluence at river_km "
          f"{(junction_m - start_m) / 1000:.2f} "
          f"({join_ll[1]:.5f} N, {join_ll[0]:.5f} E)")
    print(f"Line length: {reach.length / 1000:.2f} km")
    return out


def _write_geojson(path, to_ll, reach, flat_pt, start_pt, junction_pt, props):
    def wgs(geom):
        def project(coords):
            lon, lat = to_ll.transform(coords[:, 0], coords[:, 1])
            return np.column_stack([lon, lat])
        return transform(geom, project)

    features = []
    for role, geom in (("centerline", reach), ("flat_downstream_end", flat_pt),
                       ("start", start_pt), ("confluence", junction_pt)):
        feature_props = {"role": role, "verified_by_user": False}
        if role == "centerline":
            feature_props.update(props)
        features.append({
            "type": "Feature",
            "properties": feature_props,
            "geometry": mapping(wgs(geom)),
        })
    path.write_text(json.dumps(
        {"type": "FeatureCollection", "features": features}, indent=2),
        encoding="utf-8")


# ---------------------------------------------------------------------------
# Stations and transects (section 6, with 14.15, 14.16, 14.18, 14.19 e)
# ---------------------------------------------------------------------------

# One colour per status. Kept is green so it stays distinct from the yellow
# centerline. The first match is the colour when a transect has several reasons;
# the count table still counts every reason.
STATUS_PRIORITY = (
    "hillslope_failure",
    "tributary_fan",
    "reservoir_or_dam",
    "crosses_river",
    "intersects_neighbour",
    "tributary_confluence",
    "outside_reach",
)
STATUS_COLOUR = {
    "kept": "#2ca02c",
    "hillslope_failure": "#e31a1c",
    "tributary_fan": "#ff7f00",
    "reservoir_or_dam": "#6a3d9a",
    "crosses_river": "#1f78b4",
    "intersects_neighbour": "#33b5e5",
    "tributary_confluence": "#b15928",
    "outside_reach": "#888888",
}


def choose_status(reasons: list[str]) -> str:
    """The status drawn on the sheet. Every reason is still counted."""
    for name in STATUS_PRIORITY:
        if name in reasons:
            return name
    return "kept"


def verified_direction(river, distance_m: float, baseline_m: float) -> tuple[float, float]:
    """Unit direction along the verified line, baseline_m each side of distance_m.

    Near either end the point is clamped to the end of the line (14.21).
    """
    length = float(river.length)
    behind = max(0.0, distance_m - baseline_m)
    ahead = min(length, distance_m + baseline_m)
    if ahead - behind <= 0:
        raise RuntimeError("The verified line is too short to give a direction.")
    start = np.array(river.interpolate(behind).coords[0])
    end = np.array(river.interpolate(ahead).coords[0])
    vector = end - start
    norm = float(np.linalg.norm(vector))
    if norm == 0:
        raise RuntimeError("The direction baseline has no length.")
    vector = vector / norm
    return float(vector[0]), float(vector[1])


def usable_lengths(transects: list[np.ndarray], river_m: np.ndarray,
                   search_m: float, full_length_m: float) -> np.ndarray:
    """Distance from each station to the nearest same-bank transect crossing.

    Only stations within search_m along the river are considered. A transect
    that meets none of them keeps full_length_m (14.21).
    """
    from shapely.geometry import LineString
    lines = [LineString(np.asarray(pts)) for pts in transects]
    along = np.asarray(river_m, dtype=float)
    out = np.full(len(lines), full_length_m, dtype=float)
    for i, line in enumerate(lines):
        nearest = full_length_m
        for j, other in enumerate(lines):
            if i == j or abs(along[i] - along[j]) > search_m + 1e-6:
                continue
            if not line.intersects(other):
                continue
            for pt in _geom_points(line.intersection(other)):
                nearest = min(nearest, float(line.project(pt)))
        out[i] = nearest
    return out


def excluded_as_neighbour(usable_length_m: float, minimum_m: float) -> bool:
    """True only when the usable length is under the minimum (14.21)."""
    return float(usable_length_m) < float(minimum_m)


def _geom_points(geom) -> list:
    from shapely.geometry import Point
    if geom.is_empty:
        return []
    kind = geom.geom_type
    if kind == "Point":
        return [geom]
    if kind == "LineString":
        return [Point(geom.coords[0]), Point(geom.coords[-1])]
    if kind.startswith("Multi") or kind == "GeometryCollection":
        found = []
        for part in geom.geoms:
            found.extend(_geom_points(part))
        return found
    return []


def crosses_river_again(transect_xy: np.ndarray, river, separation_m: float) -> bool:
    """True when two meetings with the river are more than separation_m apart
    along the transect (14.21). A wiggle next to the station is one crossing.
    """
    from shapely.geometry import LineString
    line = LineString(np.asarray(transect_xy))
    hits = _geom_points(line.intersection(river))
    if len(hits) < 2:
        return False
    along = [line.project(pt) for pt in hits]
    return max(along) - min(along) > separation_m


def hits_polygon(transect_xy: np.ndarray, polygons: list) -> bool:
    from shapely.geometry import LineString
    if not polygons:
        return False
    line = LineString(np.asarray(transect_xy))
    return any(line.intersects(poly) for poly in polygons)


def manual_reason(river_km: float, side: str, exclusions: list[dict]) -> str | None:
    """The manual exclusion covering this station and bank, if there is one.

    The ends are included: a station at from_km or to_km is excluded.
    """
    for item in exclusions:
        if item["side"] == side and item["from_km"] <= river_km <= item["to_km"]:
            return item["reason"]
    return None


def stations_near_confluence(river_km: np.ndarray, upa: np.ndarray,
                             window_m: float, fraction: float,
                             distance_m: float) -> tuple[np.ndarray, list[float]]:
    """Stations within distance_m of a tributary joining with more than
    `fraction` of the mainstem's upstream area (section 6).

    Upstream area cannot shrink downstream, so a sample lower than the
    highest value so far is treated as a miss and does not create a second
    confluence when the river is found again. Non-positive values are ignored.
    """
    from .geometry import propose_section_breaks
    area = np.asarray(upa, dtype=float).copy()
    area[~np.isfinite(area) | (area <= 0)] = np.nan
    highest = np.nan
    for i in range(len(area)):
        if not np.isfinite(area[i]):
            area[i] = highest
        elif not np.isfinite(highest) or area[i] > highest:
            highest = area[i]
        else:
            area[i] = highest
    found = propose_section_breaks(river_km, area, window_m, fraction)
    km = np.asarray(river_km, dtype=float)
    near = np.zeros(len(km), dtype=bool)
    for brk in found:
        near |= np.abs(km - brk) * 1000.0 <= distance_m + 1e-6
    return near, found


def station_river_km(points_xy: np.ndarray, river) -> np.ndarray:
    """Distance in km along the verified line, for stations on the smoothed line.

    Smoothing shortens the line, and the reach end, section breaks and manual
    exclusions are all measured on the verified line.
    """
    from shapely.geometry import Point
    km = np.array([river.project(Point(xy)) / 1000.0 for xy in points_xy])
    if len(km) > 1 and np.any(np.diff(km) < -1e-6):
        at = int(np.flatnonzero(np.diff(km) < -1e-6)[0])
        raise RuntimeError(
            f"Station river_km goes backwards between stations {at} and {at + 1} "
            f"({km[at]:.4f} then {km[at + 1]:.4f}). No station was reordered."
        )
    return km


def run_transects(event_key: str, config_path: str | Path | None = None) -> Path:
    """Build every station and transect, flag the exclusions, and draw them."""
    import geopandas as gpd

    from .config import load_config
    from .geometry import chaikin, densify, inside_reach, river_sections, stations, transect
    from .label import manual_exclusions
    from .qa import line_samples
    from .stats import block_ids

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    river, samples, river_m = line_samples(event, cfg)
    smooth_pts = densify(chaikin(np.asarray(river.coords), 3), cfg.grid_m)
    smooth_length = float(np.linalg.norm(np.diff(smooth_pts, axis=0), axis=1).sum())
    station_list = stations(smooth_pts, cfg.station_spacing_m)
    xy = np.array([[s["x"], s["y"]] for s in station_list])
    km = station_river_km(xy, river)
    for station, dist_m in zip(station_list, km * 1000.0):
        station["dx"], station["dy"] = verified_direction(
            river, float(dist_m), cfg.normal_baseline_m)

    rasters = PROJECT_ROOT / cfg.paths.data / event_key / "rasters"
    post, grid = _read_bands(rasters / "s2_post.tif", ("B2", "B3", "B4"), cfg.crs)
    dem = _read_bands(rasters / "fabdem.tif", ("elevation_m",), cfg.crs, grid)[0]["elevation_m"]
    upa_grid = _read_bands(rasters / "merit_upa.tif", ("upa",), cfg.crs, grid)[0]["upa"]
    channel_path = (PROJECT_ROOT / cfg.paths.outputs / event_key / "label"
                    / "active_channel.tif")
    channel = _read_channel(channel_path, grid, cfg.crs)

    _assert_flows_downstream(dem, grid["transform"], river)
    # The station pixel is rarely the MERIT river pixel, so the on-pixel value
    # jumps up and down. The mainstem value is the largest upstream area
    # within one station spacing.
    upa = _sample_max(upa_grid, grid, xy, cfg.station_spacing_m)
    in_reach = inside_reach(km, event.reach_start_km, event.reach_end_km)
    section = river_sections(km, event.section_breaks_km)
    blocks = block_ids(km, event.reach_start_km, cfg.bootstrap_block_km,
                       cfg.bootstrap_final_block_min_km, event.reach_end_km)
    exclusions = manual_exclusions(event, event_key)
    polygons = _load_exclusion_polygons(event_key, cfg.crs)
    near_conf, confluences = stations_near_confluence(
        km, upa, cfg.confluence_distance_m, cfg.confluence_upa_fraction,
        cfg.confluence_distance_m)

    per_side = {}
    usable = {}
    for side in ("left", "right"):
        lines = [transect(s, side, cfg.max_half_length_m, cfg.grid_m) for s in station_list]
        per_side[side] = lines
        usable[side] = usable_lengths(
            lines, km * 1000.0, cfg.neighbour_search_m, cfg.max_half_length_m)

    rows = []
    resolved = []
    for i, station in enumerate(station_list):
        centred = _centre_one(station, channel, dem, grid, cfg,
                              usable["left"][i], usable["right"][i])
        resolved.append(centred["channel_resolved"])
        for side in ("left", "right"):
            lines = per_side[side]
            length_m = float(usable[side][i])
            reasons = []
            if not in_reach[i]:
                reasons.append("outside_reach")
            if excluded_as_neighbour(length_m, cfg.neighbour_min_usable_m):
                reasons.append("intersects_neighbour")
            if near_conf[i]:
                reasons.append("tributary_confluence")
            if crosses_river_again(lines[i], river, cfg.station_spacing_m):
                reasons.append("crosses_river")
            if hits_polygon(lines[i], polygons):
                reasons.append("reservoir_or_dam")
            named = manual_reason(float(km[i]), side, exclusions)
            if named:
                reasons.append(named)
            bank_top = (centred["bank_top_left_m"] if side == "left"
                        else centred["bank_top_right_m"])
            rows.append({
                "station_id": i,
                "river_km": round(float(km[i]), 4),
                "segment_id": int(np.floor(km[i])),
                "block_id": int(blocks[i]),
                "river_section": int(section[i]),
                "side": side,
                "in_reach": int(in_reach[i]),
                "channel_resolved": int(centred["channel_resolved"]),
                "centre_offset_m": centred["centre_offset_m"],
                "channel_width_m": centred["channel_width_m"],
                "bank_top_m": bank_top,
                "z0": centred["z0"],
                "usable_length_m": round(length_m, 1),
                "exclude_reasons": ",".join(reasons),
                "status": choose_status(reasons),
                "geometry": _clip_transect(lines[i], length_m, cfg.grid_m),
            })

    out_dir = PROJECT_ROOT / cfg.paths.outputs / event_key / "transects"
    out_dir.mkdir(parents=True, exist_ok=True)
    gpkg = out_dir / "transects.gpkg"
    if gpkg.exists():
        gpkg.unlink()
    gpd.GeoDataFrame(rows, crs=cfg.crs).to_file(gpkg, driver="GPKG")

    _print_transect_report(
        rows, km, confluences, upa, resolved, river.length, smooth_length,
        polygons, event_key, cfg.station_spacing_m)

    sheet_dir = out_dir / "sheets"
    sheet_dir.mkdir(parents=True, exist_ok=True)
    pngs = _write_transect_sheets(sheet_dir, grid, post, samples, river_m, rows)
    prov_path = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    if prov_path.exists():
        provenance = json.loads(prov_path.read_text(encoding="utf-8"))
        provenance["transects"] = {
            "smoothed_with": "Chaikin, 3 iterations, then densified to the grid",
            "verified_length_m": round(float(river.length), 1),
            "smoothed_length_m": round(smooth_length, 1),
            "river_km": "distance along the verified centerline",
            "direction": f"unit vector {cfg.normal_baseline_m:g} m up and downstream "
                         "on the verified line, clamped at the ends (14.21)",
            "usable_length_m": "distance to the first same-bank transect from a "
                               f"station within {cfg.neighbour_search_m:g} m; "
                               "excluded only when shorter than "
                               f"{cfg.neighbour_min_usable_m:g} m",
            "stations": int(len(station_list)),
            "channel_resolved": int(np.sum(resolved)),
            "confluence_km": confluences,
            "exclusion_polygons": len(polygons),
            "file": str(gpkg.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        }
        prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
        print(f"Updated {prov_path.relative_to(PROJECT_ROOT)}")
    for path in [gpkg, *pngs]:
        print(path.relative_to(PROJECT_ROOT))
    return gpkg


def _assert_flows_downstream(dem, transform, river) -> None:
    start = np.array(river.coords[0])
    end = np.array(river.coords[-1])
    z0 = _sample(dem, {"transform": transform, "height": dem.shape[0],
                       "width": dem.shape[1]}, np.vstack([start, end]))
    if not np.isfinite(z0).all():
        raise RuntimeError(
            "FABDEM has no elevation at an end of the centerline, so the "
            "upstream-to-downstream direction was not checked."
        )
    if z0[1] > z0[0]:
        raise RuntimeError(
            f"The verified line ends {z0[1] - z0[0]:.0f} m higher than it starts. "
            "It was not reversed, and no transects were built."
        )
    print(f"Centerline falls from {z0[0]:.0f} m to {z0[1]:.0f} m, "
          "upstream to downstream.")


def _centre_one(station, channel, dem, grid, cfg,
                usable_left_m: float, usable_right_m: float) -> dict:
    from .geometry import transect
    from .profiles import centre_channel
    step = cfg.grid_m
    left = _within_usable(transect(station, "left", cfg.max_half_length_m, step),
                          usable_left_m, step)
    right = _within_usable(transect(station, "right", cfg.max_half_length_m, step),
                           usable_right_m, step)
    points = np.vstack([right[::-1], left[1:]])
    offsets = np.concatenate([
        -np.arange(len(right))[::-1] * step,
        np.arange(1, len(left)) * step,
    ])
    sampled_channel = _sample(channel.astype(float), grid, points) == 1
    sampled_z = _sample(dem, grid, points)
    return centre_channel(sampled_channel, offsets, sampled_z, step,
                          cfg.channel_merge_gap_m, cfg.recentre_max_m)


def _within_usable(points: np.ndarray, usable_m: float, step_m: float) -> np.ndarray:
    keep = np.arange(len(points)) * step_m <= usable_m + 1e-6
    return points[keep] if keep.any() else points[:1]


def _sample_max(values: np.ndarray, grid: dict, points: np.ndarray,
                radius_m: float) -> np.ndarray:
    """Largest value within radius_m of each point."""
    from rasterio.transform import rowcol
    pixel = abs(float(grid["transform"].a))
    reach = int(round(radius_m / pixel))
    rows, cols = rowcol(grid["transform"], points[:, 0], points[:, 1])
    rows, cols = np.asarray(rows), np.asarray(cols)
    out = np.full(len(points), np.nan)
    height, width = values.shape
    for i, (row, col) in enumerate(zip(rows, cols)):
        if row < 0 or col < 0 or row >= height or col >= width:
            continue
        window = values[max(0, row - reach):min(height, row + reach + 1),
                        max(0, col - reach):min(width, col + reach + 1)]
        if np.isfinite(window).any():
            out[i] = np.nanmax(window)
    return out


def _sample(values: np.ndarray, grid: dict, points: np.ndarray) -> np.ndarray:
    from rasterio.transform import rowcol
    rows, cols = rowcol(grid["transform"], points[:, 0], points[:, 1])
    rows = np.asarray(rows)
    cols = np.asarray(cols)
    inside = ((rows >= 0) & (cols >= 0)
              & (rows < grid["height"]) & (cols < grid["width"]))
    out = np.full(len(points), np.nan)
    out[inside] = values[rows[inside], cols[inside]]
    return out


def _read_bands(path: Path, names, crs: str, grid: dict | None = None):
    import rasterio
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing.")
    with rasterio.open(path) as src:
        got = src.crs.to_string() if src.crs is not None else None
        if got != crs:
            raise RuntimeError(f"{path.name} is {got}, not {crs}.")
        here = {"transform": src.transform, "height": src.height,
                "width": src.width, "crs": crs}
        if grid is not None and (here["transform"] != grid["transform"]
                                 or (here["height"], here["width"]) != (grid["height"], grid["width"])):
            raise RuntimeError(f"{path.name} is not on the same grid as the post-event image.")
        found = list(src.descriptions)
        missing = [n for n in names if n not in found]
        if missing:
            raise RuntimeError(f"{path.name} is missing {missing}. Bands: {found}.")
        bands = {}
        for name in names:
            data = src.read(found.index(name) + 1, masked=True).astype(float)
            bands[name] = np.ma.filled(data, np.nan)
    return bands, here


def _read_channel(path: Path, grid: dict, crs: str) -> np.ndarray:
    import rasterio
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. The pre-event active channel comes from the label step."
        )
    with rasterio.open(path) as src:
        got = src.crs.to_string() if src.crs is not None else None
        if got != crs or src.transform != grid["transform"]:
            raise RuntimeError(f"{path.name} is not on the analysis grid.")
        band = src.read(1)
    return band == 1


def _load_exclusion_polygons(event_key: str, crs: str) -> list:
    """Reservoirs and dams from data/{event}/exclusions.geojson, or none."""
    import geopandas as gpd
    path = PROJECT_ROOT / "data" / event_key / "exclusions.geojson"
    if not path.exists():
        print(f"{path.relative_to(PROJECT_ROOT)} is not present, so no reservoir "
              "or dam was excluded.")
        return []
    gdf = gpd.read_file(path).to_crs(crs)
    polygons = [geom for geom in gdf.geometry if geom is not None and not geom.is_empty]
    print(f"Exclusion polygons: {len(polygons)} from {path.relative_to(PROJECT_ROOT)}")
    return polygons


def _print_transect_report(rows, km, confluences, upa, resolved, verified_length,
                           smooth_length, polygons, event_key, spacing_m) -> None:
    print(f"Verified line {verified_length:.1f} m; smoothed line {smooth_length:.1f} m "
          f"({smooth_length - verified_length:+.1f} m).")
    print("river_km is measured along the verified line, so it matches the "
          "contact sheets, the reach end and the manual exclusions.")
    print(f"Stations: {len(km)}, every {spacing_m:.0f} m. "
          f"Pre-event channel resolved at {int(np.sum(resolved))} of {len(km)}.")
    if confluences:
        print("Tributary confluences: the largest MERIT upstream area within "
              f"{spacing_m:.0f} m rises by more than 5% inside 300 m.")
        print("A later return to a river already reached is not counted.")
        for brk in confluences:
            i = int(np.argmin(np.abs(km - brk)))
            before = upa[max(0, i - 3)]
            print(f"  {brk:.3f} km, upstream area {before:.0f} -> {upa[i]:.0f}")
    else:
        print("No tributary confluence above 5% of the mainstem upstream area.")
    if not polygons:
        print(f"data/{event_key}/exclusions.geojson: none supplied.")

    reasons = ["kept", *STATUS_PRIORITY]
    extra = sorted({r["status"] for r in rows} - set(reasons) - {"kept"})
    print("Transect counts. A transect with more than one reason is counted in each.")
    print(f"{'reason':<24}{'left':>8}{'right':>8}{'total':>8}")
    for reason in reasons + extra:
        left = sum(1 for r in rows if r["side"] == "left" and _has(r, reason))
        right = sum(1 for r in rows if r["side"] == "right" and _has(r, reason))
        print(f"{reason:<24}{left:8d}{right:8d}{left + right:8d}")
    _print_usable(rows)


def _print_usable(rows: list[dict]) -> None:
    values = np.array([row["usable_length_m"] for row in rows], dtype=float)
    minimum, lower, median, upper = np.percentile(values, [0, 25, 50, 75])
    print("usable_length_m, both banks:")
    print(f"  minimum {minimum:.0f} m, lower quartile {lower:.0f} m, "
          f"median {median:.0f} m, upper quartile {upper:.0f} m")
    print(f"  share under 200 m: {100.0 * np.mean(values < 200):.1f}%")
    print(f"  share under 300 m: {100.0 * np.mean(values < 300):.1f}%")


def _clip_transect(points: np.ndarray, length_m: float, step_m: float):
    from shapely.geometry import LineString
    full = LineString(np.asarray(points))
    if length_m >= full.length - 1e-6:
        return full
    count = int(np.floor((length_m + 1e-6) / step_m))
    kept = [points[i] for i in range(count + 1)]
    end = full.interpolate(length_m)
    if np.linalg.norm(np.asarray(kept[-1]) - (end.x, end.y)) > 1e-3:
        kept.append((end.x, end.y))
    return LineString(kept)


def _has(row: dict, reason: str) -> bool:
    if reason == "kept":
        return row["status"] == "kept"
    return reason in row["exclude_reasons"].split(",") if row["exclude_reasons"] else False


def _write_transect_sheets(out_dir, grid, post, samples, river_m, rows) -> list:
    from matplotlib.lines import Line2D

    from . import sheets
    limits = sheets.fixed_limits([post], sheets.TRUE_COLOUR)
    present = []
    for name in ("kept", *STATUS_PRIORITY):
        if any(r["status"] == name for r in rows):
            present.append(name)
    legend = [Line2D([0], [0], color=STATUS_COLOUR[name], lw=2, label=name.replace("_", " "))
              for name in present]
    legend.append(Line2D([0], [0], color="yellow", lw=1.5, label="centerline"))

    def draw(ax, win):
        sheets.show(ax, sheets.rgb(post, sheets.TRUE_COLOUR, limits, win), win)
        for name in ("kept", *STATUS_PRIORITY):
            for row in rows:
                if row["status"] != name:
                    continue
                if not (win["start_m"] <= row["river_km"] * 1000 <= win["end_m"]):
                    continue
                coords = np.asarray(row["geometry"].coords)
                ax.plot(coords[:, 0], coords[:, 1], color=STATUS_COLOUR[name],
                        lw=1.1, solid_capstyle="round")
        sheets.finish_panel(ax, win, line_width=1.2)

    return sheets.write_sheets(
        out_dir, "transects",
        "Transects over the post-event image, each drawn only to its usable "
        "length. The frame is the usual 250 m around the river",
        grid, samples, river_m, draw, legend)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="python -m floodforest.centerline")
    parser.add_argument("command", choices=["draft"])
    parser.add_argument("--event", required=True)
    args = parser.parse_args(argv)
    draft_centerline(args.event)


if __name__ == "__main__":
    main()
