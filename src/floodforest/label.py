"""Impact label (spec section 5, with clarifications 14.11, 14.12, 14.17, 14.18).

Labels are 1 impacted, 0 not impacted, nodata. Five are written:
combined (primary), optical_only, sar_only, and the radar sensitivities
sar_strict (all valid orbits agree) and sar_near (near-event baseline).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.transform import rowcol, xy
from scipy import ndimage
from scipy.spatial import cKDTree

from . import optical, sar, sheets
from .config import PROJECT_ROOT, load_config
from .geometry import chaikin, densify, normals, stations
from .profiles import centre_channel
from .qa import line_samples

NODATA = 255
_EIGHT = np.ones((3, 3), dtype=int)
_S2_BANDS = ("B2", "B3", "B4", "B8", "B11")
LABELS = ("combined", "optical_only", "sar_only", "sar_strict", "sar_near")


def run_label(event_key: str, config_path: str | Path | None = None) -> list[Path]:
    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    rasters = PROJECT_ROOT / cfg.paths.data / event_key / "rasters"
    prov_path = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    if not prov_path.exists():
        raise FileNotFoundError(f"{prov_path} is missing. Run fetch first.")
    provenance = json.loads(prov_path.read_text(encoding="utf-8"))

    line, samples, river_m = line_samples(event, cfg)
    pre, grid = _read_named(rasters / "s2_pre.tif", _S2_BANDS, cfg.crs)
    post, _ = _read_named(rasters / "s2_post.tif", _S2_BANDS, cfg.crs, grid)
    dem = _read_named(rasters / "fabdem.tif", ("elevation_m",), cfg.crs, grid)[0]["elevation_m"]
    pixel_m = grid["pixel_m"]

    corridor = _burn(line.buffer(cfg.corridor_buffer_m), grid, all_touched=False)
    seed = _burn(line, grid, all_touched=True)

    # Optical (5.2 with 14.18).
    ndvi_pre = optical.ndvi(pre["B8"], pre["B4"])
    mndwi_pre = optical.mndwi(pre["B3"], pre["B11"])
    spectral = optical.spectral_channel(
        ndvi_pre, mndwi_pre, cfg.ndvi_channel_max, cfg.mndwi_channel_min) & corridor
    near_line = _within_line(line, grid, cfg.channel_seed_max_m)
    channel = optical.keep_components_touching(spectral, near_line)
    if not channel.any():
        raise RuntimeError(
            f"No pre-event active-channel component comes within "
            f"{cfg.channel_seed_max_m} m of the centerline. No other mask was used.")
    ndvi_post = optical.ndvi(post["B8"], post["B4"])
    bsi_post = optical.bsi(post["B2"], post["B4"], post["B8"], post["B11"])
    t_bsi = _frozen_t_bsi(provenance, optical.t_bsi_from_channel(
        bsi_post, channel, cfg.bsi_percentile))
    dark = optical.dark_mask(pre["B8"], post["B8"], channel, cfg.dark_b8_threshold) & corridor
    optical_raw = optical.optical_impact(
        ndvi_pre, ndvi_post, bsi_post, channel, t_bsi,
        cfg.ndvi_unknown_max, cfg.ndvi_deposit_max, dark=dark)
    optical_raw[~corridor] = np.nan
    exclusions = manual_exclusions(event, event_key, cfg)

    # Stations, z0 and the pixel -> nearest-station lookup (14.11, 14.18).
    slope, aspect = sar.terrain_slope_aspect(dem, pixel_m)
    station_list = _stations(line, cfg)
    z0 = np.array([_station_z0(s, channel, dem, grid, cfg) for s in station_list])
    resolved = np.array([r["channel_resolved"] for r in z0], dtype=int)
    z0 = np.array([r["z0"] for r in z0], dtype=float)
    nearest = _nearest_station(station_list, corridor, grid)
    pixel_z0 = np.full(dem.shape, np.nan)
    pixel_z0[corridor] = z0[nearest]
    pixel_km = np.full(dem.shape, -1, dtype=int)
    station_km = np.array([s["river_km"] for s in station_list])
    pixel_km[corridor] = np.floor(station_km[nearest]).astype(int)

    # Radar (5.3, 14.11, 14.12, 14.17).
    orbit_geometry = _orbit_geometry(provenance, rasters / "s1_angle.tif",
                                     slope, aspect, corridor, grid, cfg)
    post_medians = _orbit_medians(rasters / "s1_post.tif", provenance, "s1_post", cfg, grid)
    pre_medians = _orbit_medians(rasters / "s1_pre.tif", provenance, "s1_pre", cfg, grid)
    near_medians = _orbit_medians(rasters / "s1_near.tif", provenance, "s1_near", cfg, grid)
    stable = (corridor & np.isfinite(dem) & np.isfinite(pixel_z0)
              & (dem >= pixel_z0 + cfg.reference_height_above_channel_m)
              & (slope < cfg.reference_slope_max_deg)
              & (optical_raw != 1))
    primary = _sar_votes(post_medians, pre_medians, orbit_geometry, stable,
                         corridor, cfg, "s1_pre")
    near = _sar_votes(post_medians, near_medians, orbit_geometry, stable,
                      corridor, cfg, "s1_near")

    raw = {
        "optical_only": optical_raw,
        "sar_only": _vote(primary, cfg, strict=False),
        "sar_strict": _vote(primary, cfg, strict=True),
        "sar_near": _vote(near, cfg, strict=False),
    }
    for name in ("sar_only", "sar_strict", "sar_near"):
        raw[name][channel] = np.nan
    raw["combined"] = combine(raw["optical_only"], raw["sar_only"])
    raw["combined"][channel] = np.nan
    anchor = channel | seed
    labels = {name: connectivity_filter(raw[name], anchor, cfg.morph_close_radius_px)
              for name in LABELS}
    source = label_source(labels["combined"], raw["optical_only"], raw["sar_only"])

    out_dir = PROJECT_ROOT / cfg.paths.outputs / event_key / "label"
    out_dir.mkdir(parents=True, exist_ok=True)
    written = _write_rasters(out_dir, grid, labels, source, channel, seed,
                             orbit_geometry, corridor, dark)

    _print_summary(t_bsi, channel, primary, near, orbit_geometry, resolved)
    _print_km_table(pixel_km, corridor, labels, orbit_geometry, pixel_m, dark)
    forest = _pre_event_forest(rasters, provenance, cfg, grid)
    radar_only = _radar_only_breakdown(source, slope, forest, corridor, cfg, pixel_m)
    _print_exclusions(exclusions, event)

    provenance["label"] = _label_record(cfg, t_bsi, channel, resolved, primary,
                                        near, orbit_geometry, written)
    provenance["label"]["dark_pixels"] = int(dark.sum())
    provenance["label"]["dark_rule"] = (
        f"B8 < {cfg.dark_b8_threshold} in the pre or post composite, outside "
        "the active channel: optical nodata (14.19 a)")
    provenance["label"]["combined_rule"] = (
        "1 if optical or radar is 1; 0 only where optical is a valid 0; "
        "pre-event active channel is nodata in the combined and radar labels (14.19 b, c)")
    provenance["label"]["radar_only_breakdown"] = radar_only
    provenance["manual_exclusions"] = {
        "list": exclusions,
        "confirmed_on": getattr(event, "manual_exclusions_confirmed_on", None),
        "note": "Station x side rows in these stretches stay in the table with "
                "the reason code and are left out of the primary models (14.19 e).",
    }
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"Updated {prov_path.relative_to(PROJECT_ROOT)}")

    sheet_dir = out_dir / "sheets"
    sheet_dir.mkdir(parents=True, exist_ok=True)
    pngs = _write_label_sheets(sheet_dir, grid, pre, post, labels, source,
                               samples, river_m, exclusions)
    for path in written + pngs:
        print(path.relative_to(PROJECT_ROOT))
    return written + pngs


# ---------------------------------------------------------------------------
# Combining and the connectivity filter (5.4 with 14.18)
# ---------------------------------------------------------------------------

def combine(optical_label: np.ndarray, sar_label: np.ndarray) -> np.ndarray:
    """1 if optical or radar is 1; 0 only where optical is a valid 0 (14.19 b)."""
    out = np.full(optical_label.shape, np.nan)
    out[optical_label == 0] = 0.0
    out[(optical_label == 1) | (sar_label == 1)] = 1.0
    return out


def connectivity_filter(label: np.ndarray, anchor: np.ndarray,
                        close_radius_px: int) -> np.ndarray:
    """Keep impact pixels 8-connected to the anchor (channel ∪ channel_seed).

    A closing of the given radius bridges 1-pixel gaps for the connectivity
    test only; pixels it fills are not added to the label.
    """
    impact = label == 1
    closed = ndimage.binary_closing(impact, structure=_square(close_radius_px)) | impact
    network = closed | np.asarray(anchor, dtype=bool)
    components, _ = ndimage.label(network, structure=_EIGHT)
    ids = np.unique(components[np.asarray(anchor, dtype=bool)])
    ids = ids[ids != 0]
    kept = np.isin(components, ids) & impact
    out = label.copy()
    out[impact & ~kept] = 0.0
    return out


def label_source(combined: np.ndarray, optical_raw: np.ndarray,
                 sar_raw: np.ndarray) -> np.ndarray:
    """0 not impacted, 1 optical, 2 radar, 3 both, NODATA nodata."""
    out = np.full(combined.shape, NODATA, dtype=np.uint8)
    out[np.isfinite(combined)] = 0
    hit = combined == 1
    out[hit & (optical_raw == 1)] = 1
    out[hit & (sar_raw == 1)] = 2
    out[hit & (optical_raw == 1) & (sar_raw == 1)] = 3
    return out


def _square(radius: int) -> np.ndarray:
    # A cross of radius 1 cannot bridge a 1-pixel gap in a 1-pixel-wide strip.
    r = int(radius)
    return np.ones((2 * r + 1, 2 * r + 1), dtype=bool)


# ---------------------------------------------------------------------------
# Manual exclusions (14.19 e) and radar-only reporting (14.19 f)
# ---------------------------------------------------------------------------

EXCLUSION_SIDES = ("left", "right")
EXCLUSION_REASONS = ("hillslope_failure", "tributary_fan")


def manual_exclusions(event, event_key: str, cfg=None) -> list[dict]:
    """Validated manual_exclusions for one event. A missing key is an error."""
    if not hasattr(event, "manual_exclusions"):
        raise RuntimeError(
            f"{event_key} has no manual_exclusions in config.yaml. Add a list "
            "(use [] when there are none). None were assumed.")
    out = []
    for item in event.manual_exclusions or []:
        entry = {k: getattr(item, k, None) for k in ("from_km", "to_km", "side", "reason")}
        if entry["side"] not in EXCLUSION_SIDES:
            raise RuntimeError(f"manual_exclusions side must be left or right: {entry}")
        if entry["reason"] not in EXCLUSION_REASONS:
            raise RuntimeError(
                f"manual_exclusions reason must be one of {EXCLUSION_REASONS}: {entry}")
        if not (isinstance(entry["from_km"], (int, float))
                and isinstance(entry["to_km"], (int, float))
                and entry["from_km"] < entry["to_km"]):
            raise RuntimeError(f"manual_exclusions needs from_km < to_km: {entry}")
        if entry["from_km"] < event.reach_start_km or entry["to_km"] > event.reach_end_km:
            raise RuntimeError(
                f"manual_exclusions stretch {entry} lies outside the reach "
                f"{event.reach_start_km}–{event.reach_end_km} km.")
        out.append({**entry, "from_km": float(entry["from_km"]),
                    "to_km": float(entry["to_km"])})
    return out


def _print_exclusions(exclusions: list[dict], event) -> None:
    confirmed = getattr(event, "manual_exclusions_confirmed_on", None)
    print(f"Manual exclusions: {len(exclusions)} stretches; confirmed on: "
          f"{confirmed or 'not yet confirmed by the user'}")
    for item in exclusions:
        print(f"  {item['side']:5s} km {item['from_km']:g}–{item['to_km']:g}  "
              f"{item['reason']}")


def _pre_event_forest(rasters: Path, provenance: dict, cfg, grid) -> np.ndarray:
    """Forest pixel as in 7.3: 1 forest, 0 non-forest, NaN no data."""
    from .exposure import forest_mask
    trees = _read_named(rasters / "dynamic_world_trees_pre.tif", ("trees",),
                        cfg.crs, grid)[0]["trees"]
    height = _read_named(rasters / "canopy_height.tif", ("canopy_height_m",),
                         cfg.crs, grid)[0]["canopy_height_m"]
    loss = _read_named(rasters / "hansen_loss_year.tif", ("loss_year",),
                       cfg.crs, grid)[0]["loss_year"]
    hansen = provenance["hansen"]
    forest = forest_mask(trees, height, np.nan_to_num(loss, nan=0.0),
                         hansen["pre_composite_year"], cfg.trees_probability_min,
                         cfg.canopy_height_min_m, cfg.hansen_loss_from_year,
                         hansen["newest_year"])
    out = forest.astype(float)
    out[~(np.isfinite(trees) & np.isfinite(height))] = np.nan
    return out


def _radar_only_breakdown(source, slope, forest, corridor, cfg, pixel_m) -> dict:
    """Radar-only impact area by slope class and pre-event forest (14.19 f)."""
    ha = pixel_m * pixel_m / 10_000
    low, high = cfg.radar_only_slope_classes_deg
    radar_only = (source == 2) & corridor
    slopes = {f"under {low} deg": slope < low,
              f"{low}-{high} deg": (slope >= low) & (slope <= high),
              f"over {high} deg": slope > high}
    covers = {"forest": forest == 1, "non-forest": forest == 0,
              "no forest data": ~np.isfinite(forest)}
    table = {}
    print("Radar-only impact (combined label, radar fired and optical did not), ha:")
    print(f"  {'slope class':14s}" + "".join(f"{c:>16s}" for c in covers) + f"{'total':>10s}")
    for s_name, s_mask in slopes.items():
        row = {c: round(int((radar_only & s_mask & c_mask).sum()) * ha, 2)
               for c, c_mask in covers.items()}
        row["total"] = round(int((radar_only & s_mask).sum()) * ha, 2)
        table[s_name] = row
        print(f"  {s_name:14s}" + "".join(f"{row[c]:16.2f}" for c in covers)
              + f"{row['total']:10.2f}")
    total = {c: round(int((radar_only & c_mask).sum()) * ha, 2) for c, c_mask in covers.items()}
    total["total"] = round(int(radar_only.sum()) * ha, 2)
    table["all slopes"] = total
    print(f"  {'all slopes':14s}" + "".join(f"{total[c]:16.2f}" for c in covers)
          + f"{total['total']:10.2f}")
    return table


def exclusion_bar(samples: np.ndarray, river_m: np.ndarray, item: dict,
                  inner_m: float, outer_m: float) -> np.ndarray | None:
    """Polygon beside the river on the excluded side, between inner_m and outer_m."""
    keep = np.flatnonzero((river_m >= item["from_km"] * 1000)
                          & (river_m <= item["to_km"] * 1000))
    if keep.size < 2:
        return None
    reach = 10
    ahead = samples[np.minimum(keep + reach, len(samples) - 1)]
    behind = samples[np.maximum(keep - reach, 0)]
    d = ahead - behind
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    left = np.column_stack([-d[:, 1], d[:, 0]])
    normal = left if item["side"] == "left" else -left
    pts = samples[keep]
    return np.vstack([pts + normal * inner_m, (pts + normal * outer_m)[::-1]])


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def _read_named(path: Path, names, crs: str, grid: dict | None = None):
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing. Run fetch first.")
    with rasterio.open(path) as src:
        got = src.crs.to_string() if src.crs is not None else None
        if got != crs:
            raise RuntimeError(f"{path.name} is {got}, not {crs}.")
        here = {"transform": src.transform, "height": src.height,
                "width": src.width, "pixel_m": abs(float(src.transform.a)),
                "crs": crs}
        if grid is not None and (here["transform"] != grid["transform"]
                                 or here["height"] != grid["height"]
                                 or here["width"] != grid["width"]):
            raise RuntimeError(f"{path.name} is not on the s2_pre grid.")
        found = list(src.descriptions)
        missing = [n for n in names if n not in found]
        if missing:
            raise RuntimeError(f"{path.name} is missing {missing}. Bands: {found}.")
        bands = {}
        for name in names:
            data = src.read(found.index(name) + 1, masked=True).astype(float)
            data = np.ma.filled(data, np.nan)
            data[~np.isfinite(data)] = np.nan
            bands[name] = data
    return bands, here


def _burn(geometry, grid: dict, all_touched: bool) -> np.ndarray:
    return rasterize([(geometry, 1)], out_shape=(grid["height"], grid["width"]),
                     transform=grid["transform"], fill=0, dtype="uint8",
                     all_touched=all_touched).astype(bool)


def _within_line(line, grid: dict, distance_m: float) -> np.ndarray:
    """Pixels whose centre lies within distance_m of the line."""
    import shapely
    band = _burn(line.buffer(distance_m + grid["pixel_m"]), grid, all_touched=True)
    rows, cols = np.nonzero(band)
    xs, ys = xy(grid["transform"], rows, cols)
    dist = shapely.distance(shapely.points(np.asarray(xs), np.asarray(ys)), line)
    out = np.zeros(band.shape, dtype=bool)
    out[rows[dist <= distance_m], cols[dist <= distance_m]] = True
    return out


def _frozen_t_bsi(provenance: dict, computed: float) -> float:
    """t_bsi is computed once per event and then frozen in provenance (5.2 step 5)."""
    frozen = provenance.get("label", {}).get("t_bsi")
    if frozen is None:
        print(f"t_bsi computed and frozen: {computed:.4f}")
        return computed
    if not math.isclose(frozen, computed, rel_tol=0, abs_tol=1e-9):
        print(f"t_bsi frozen at {frozen:.4f}; the inputs now give {computed:.4f}. "
              "The frozen value is used. Delete label.t_bsi in provenance.json "
              "to recompute it.")
    else:
        print(f"t_bsi frozen value used: {frozen:.4f}")
    return float(frozen)


# ---------------------------------------------------------------------------
# Stations and z0 (14.18)
# ---------------------------------------------------------------------------

def _stations(line, cfg) -> list[dict]:
    smooth = densify(chaikin(np.asarray(line.coords), 3), cfg.grid_m)
    return stations(smooth, cfg.station_spacing_m)


def _station_z0(station: dict, channel, dem, grid, cfg) -> dict:
    step = cfg.grid_m
    count = int(round(cfg.max_half_length_m / step))
    offsets = np.arange(-count, count + 1) * step
    left = normals(station["dx"], station["dy"])["left"]
    points = np.array([station["x"], station["y"]]) + offsets[:, None] * left
    rows, cols = rowcol(grid["transform"], points[:, 0], points[:, 1])
    rows, cols = np.asarray(rows), np.asarray(cols)
    inside = (rows >= 0) & (cols >= 0) & (rows < grid["height"]) & (cols < grid["width"])
    ch = np.zeros(len(offsets), dtype=bool)
    z = np.full(len(offsets), np.nan)
    ch[inside] = channel[rows[inside], cols[inside]]
    z[inside] = dem[rows[inside], cols[inside]]
    return centre_channel(ch, offsets, z, step, cfg.channel_merge_gap_m,
                          cfg.recentre_max_m)


def _nearest_station(station_list, corridor, grid) -> np.ndarray:
    rows, cols = np.nonzero(corridor)
    xs, ys = xy(grid["transform"], rows, cols)
    tree = cKDTree(np.array([[s["x"], s["y"]] for s in station_list]))
    _, index = tree.query(np.column_stack([xs, ys]))
    return index


# ---------------------------------------------------------------------------
# Radar
# ---------------------------------------------------------------------------

def _orbit_key(orbit, direction) -> str:
    return f"{orbit} {direction}"


def _orbit_geometry(provenance, angle_path, slope, aspect, corridor, grid, cfg) -> dict:
    record = provenance.get("s1_geometry")
    if not record:
        raise RuntimeError(
            "provenance.json has no s1_geometry. Run "
            "`fetch --only s1_geometry` first.")
    names = tuple(o["band"] for o in record["orbits"])
    angles, _ = _read_named(angle_path, names, cfg.crs, grid)
    out = {}
    for item in record["orbits"]:
        layover, shadow = sar.layover_shadow(
            slope, aspect, angles[item["band"]], item["heading_deg"],
            cfg.layover_shadow_buffer_m, grid["pixel_m"])
        key = _orbit_key(item["relativeOrbitNumber_start"], item["orbitProperties_pass"])
        out[key] = {
            "layover": layover & corridor,
            "shadow": shadow & corridor,
            "heading_deg": item["heading_deg"],
            "angle_valid": np.isfinite(angles[item["band"]]),
        }
    return out


def _orbit_medians(path, provenance, group_name, cfg, grid) -> dict:
    group = provenance.get("composites", {}).get(group_name)
    if not group:
        raise RuntimeError(f"provenance.json has no {group_name} scenes. Run fetch first.")
    scenes = group["scenes"]
    names = [f"{pol}_{i}" for i in range(len(scenes)) for pol in cfg.assets.s1_polarisations]
    bands, _ = _read_named(path, names, cfg.crs, grid)
    with rasterio.open(path) as src:
        extra = [d for d in src.descriptions
                 if d.split("_")[0] in cfg.assets.s1_polarisations and d not in names]
    if extra:
        raise RuntimeError(
            f"{path.name} has bands {extra} that match no scene in provenance.")
    by_orbit: dict[str, list[int]] = {}
    for i, scene in enumerate(scenes):
        key = _orbit_key(scene["relativeOrbitNumber_start"], scene["orbitProperties_pass"])
        by_orbit.setdefault(key, []).append(i)
    out = {}
    for key, indices in by_orbit.items():
        if len(indices) < cfg.sar_scenes_min:
            out[key] = {"skipped": f"{len(indices)} scenes, fewer than {cfg.sar_scenes_min}"}
            continue
        medians = {}
        for pol in cfg.assets.s1_polarisations:
            filtered = [sar.focal_median(bands[f"{pol}_{i}"], cfg.sar_focal_median_px)
                        for i in indices]
            medians[pol] = sar.temporal_median(filtered, cfg.sar_scenes_min)
        out[key] = {"medians": medians, "scenes": len(indices)}
    return out


def _sar_votes(post: dict, baseline: dict, geometry: dict, stable: np.ndarray,
               corridor: np.ndarray, cfg, baseline_name: str) -> dict:
    """Per-orbit change, reference zone and z-scores. Returns fired/valid per orbit."""
    pols = cfg.assets.s1_polarisations
    orbits = {}
    for key in sorted(post):
        entry = {"baseline": baseline_name}
        reason = None
        if "skipped" in post[key]:
            reason = f"post window: {post[key]['skipped']}"
        elif key not in baseline:
            reason = f"no {baseline_name} scenes on this orbit"
        elif "skipped" in baseline[key]:
            reason = f"{baseline_name}: {baseline[key]['skipped']}"
        elif key not in geometry:
            reason = "no radar geometry for this orbit"
        if reason:
            entry["skipped"] = reason
            orbits[key] = entry
            print(f"  orbit {key} ({baseline_name}) skipped: {reason}")
            continue
        masked = geometry[key]["layover"] | geometry[key]["shadow"]
        deltas = {pol: post[key]["medians"][pol] - baseline[key]["medians"][pol]
                  for pol in pols}
        valid = corridor & ~masked
        for pol in pols:
            valid &= np.isfinite(deltas[pol])
        reference = stable & valid
        n_ref = int(reference.sum())
        entry["reference_pixels"] = n_ref
        if n_ref < cfg.reference_pixels_min:
            entry["skipped"] = (f"{n_ref} reference pixels, fewer than "
                                f"{cfg.reference_pixels_min}")
            orbits[key] = entry
            print(f"  orbit {key} ({baseline_name}) skipped: {entry['skipped']}")
            continue
        fired = np.zeros(corridor.shape, dtype=bool)
        entry["reference"] = {}
        for pol in pols:
            z, stats = sar.robust_z(deltas[pol], reference)
            entry["reference"][pol] = stats
            with np.errstate(invalid="ignore"):
                fired |= np.abs(z) > cfg.sar_z_threshold
        entry["fired"] = fired & valid
        entry["valid"] = valid
        orbits[key] = entry
    return orbits


def _vote(orbits: dict, cfg, strict: bool) -> np.ndarray:
    used = [o for o in orbits.values() if "fired" in o]
    if not used:
        raise RuntimeError("No radar orbit is usable, so no radar label was made.")
    return sar.agreement([o["fired"] for o in used], [o["valid"] for o in used],
                         cfg.sar_agreement_fraction, strict=strict)


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

def _profile(grid: dict, count: int) -> dict:
    return {"driver": "GTiff", "height": grid["height"], "width": grid["width"],
            "count": count, "dtype": "uint8", "crs": grid["crs"],
            "transform": grid["transform"], "nodata": NODATA,
            "compress": "deflate"}


def _to_uint8(label: np.ndarray, corridor: np.ndarray) -> np.ndarray:
    out = np.full(label.shape, NODATA, dtype=np.uint8)
    ok = np.isfinite(label) & corridor
    out[ok] = label[ok].astype(np.uint8)
    return out


def _write_rasters(out_dir, grid, labels, source, channel, seed, geometry,
                   corridor, dark) -> list[Path]:
    written = []
    for name in LABELS:
        path = out_dir / f"impact_{name}.tif"
        with rasterio.open(path, "w", **_profile(grid, 1)) as dst:
            dst.write(_to_uint8(labels[name], corridor), 1)
            dst.set_band_description(1, f"impact_{name}")
        written.append(path)
    path = out_dir / "label_source.tif"
    with rasterio.open(path, "w", **_profile(grid, 1)) as dst:
        dst.write(np.where(corridor, source, NODATA).astype(np.uint8), 1)
        dst.set_band_description(1, "0 none, 1 optical, 2 radar, 3 both")
    written.append(path)
    path = out_dir / "active_channel.tif"
    with rasterio.open(path, "w", **_profile(grid, 2)) as dst:
        dst.write(np.where(corridor, channel, NODATA).astype(np.uint8), 1)
        dst.write(np.where(corridor, seed, NODATA).astype(np.uint8), 2)
        dst.set_band_description(1, "active_channel")
        dst.set_band_description(2, "channel_seed")
    written.append(path)
    path = out_dir / "dark_mask.tif"
    with rasterio.open(path, "w", **_profile(grid, 1)) as dst:
        dst.write(np.where(corridor, dark, NODATA).astype(np.uint8), 1)
        dst.set_band_description(1, "dark: optical nodata (14.19 a)")
    written.append(path)
    path = out_dir / "layover_shadow.tif"
    keys = sorted(geometry)
    with rasterio.open(path, "w", **_profile(grid, len(keys))) as dst:
        for i, key in enumerate(keys, start=1):
            g = geometry[key]
            code = g["layover"].astype(np.uint8) + 2 * g["shadow"].astype(np.uint8)
            dst.write(np.where(corridor, code, NODATA).astype(np.uint8), i)
            dst.set_band_description(
                i, f"orbit {key}: 0 ok, 1 layover, 2 shadow, 3 both")
    written.append(path)
    return written


def _print_summary(t_bsi, channel, primary, near, geometry, resolved) -> None:
    print(f"Pre-event active channel: {int(channel.sum())} pixels.")
    print(f"t_bsi = {t_bsi:.4f}")
    print(f"Stations with a resolved pre-event channel: "
          f"{int(resolved.sum())} of {resolved.size}")
    for name, orbits in (("baseline s1_pre (primary)", primary),
                         ("baseline s1_near (sar_near)", near)):
        print(f"Radar reference zone, {name}:")
        for key, entry in sorted(orbits.items()):
            if "skipped" in entry:
                print(f"  orbit {key}: skipped, {entry['skipped']}")
                continue
            parts = ", ".join(
                f"{pol} median {s['median_db']:+.2f} dB, MAD {s['mad_db']:.2f} dB"
                for pol, s in entry["reference"].items())
            print(f"  orbit {key}: {entry['reference_pixels']} pixels; {parts}")
    for key, g in sorted(geometry.items()):
        print(f"Orbit {key}: look direction {g['heading_deg']:.2f} deg; "
              f"layover {int(g['layover'].sum())} px, shadow {int(g['shadow'].sum())} px")


def _print_km_table(pixel_km, corridor, labels, geometry, pixel_m, dark) -> None:
    ha = pixel_m * pixel_m / 10_000
    keys = sorted(geometry)
    header = ("river_km  corridor_ha  optical_ha  radar_ha  combined_ha  "
              "sar_strict_ha  sar_near_ha  dark_%  combined_nodata_%  "
              + "  ".join(f"lost_{k.replace(' ', '_')}_%" for k in keys))
    print("Impacted area per river km (label after the connectivity filter).")
    print("dark_% = share of the corridor masked as dark; combined_nodata_% = "
          "share of the corridor where the combined label is nodata;")
    print("lost_% = share of the corridor in layover or shadow for that orbit.")
    print(header)
    rows = list(range(int(pixel_km[corridor].min()), int(pixel_km[corridor].max()) + 1))
    for km in rows + [None]:
        cell = corridor & (pixel_km == km) if km is not None else corridor
        n = int(cell.sum())
        values = [n * ha] + [int((labels[name][cell] == 1).sum()) * ha for name in
                             ("optical_only", "sar_only", "combined", "sar_strict",
                              "sar_near")]
        lost = [100.0 * float((geometry[k]["layover"] | geometry[k]["shadow"])[cell].sum()) / n
                for k in keys]
        dark_pct = 100.0 * float(dark[cell].sum()) / n
        nodata_pct = 100.0 * float((~np.isfinite(labels["combined"][cell])).sum()) / n
        tag = f"{km:8d}" if km is not None else "   total"
        print(f"{tag}  {values[0]:11.1f}  {values[1]:10.2f}  {values[2]:8.2f}  "
              f"{values[3]:11.2f}  {values[4]:13.2f}  {values[5]:11.2f}  "
              f"{dark_pct:6.1f}  {nodata_pct:17.1f}  "
              + "  ".join(f"{v:{len(k) + 9}.1f}" for k, v in zip(keys, lost)))


def _label_record(cfg, t_bsi, channel, resolved, primary, near, geometry, written) -> dict:
    def orbits(entries):
        out = {}
        for key, entry in sorted(entries.items()):
            out[key] = {k: v for k, v in entry.items() if k not in ("fired", "valid")}
        return out

    return {
        "t_bsi": t_bsi,
        "t_bsi_rule": f"{cfg.bsi_percentile}th percentile of BSI_post inside "
                      "the pre-event active channel",
        "active_channel_pixels": int(channel.sum()),
        "channel_rule": f"MNDWI_pre > {cfg.mndwi_channel_min} or NDVI_pre < "
                        f"{cfg.ndvi_channel_max}; component within "
                        f"{cfg.channel_seed_max_m} m of the centerline (14.18)",
        "stations": int(resolved.size),
        "stations_channel_resolved": int(resolved.sum()),
        "dem": "FABDEM",
        "sar_db_note": "Sentinel-1 GRD on Earth Engine is already in dB; not converted again",
        "sar_primary": orbits(primary),
        "sar_near": orbits(near),
        "layover_shadow": {
            key: {"layover_pixels": int(g["layover"].sum()),
                  "shadow_pixels": int(g["shadow"].sum()),
                  "heading_deg": g["heading_deg"]}
            for key, g in sorted(geometry.items())
        },
        "files": [str(p.relative_to(PROJECT_ROOT)).replace("\\", "/") for p in written],
    }


def redraw_combined_sheets(event_key: str, config_path: str | Path | None = None) -> list[Path]:
    """Redraw the label_combined sheets from the saved label_source.tif.

    The label is not recomputed. The current manual_exclusions are drawn and
    copied into provenance with their confirmation date.
    """
    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    rasters = PROJECT_ROOT / cfg.paths.data / event_key / "rasters"
    out_dir = PROJECT_ROOT / cfg.paths.outputs / event_key / "label"
    _, samples, river_m = line_samples(event, cfg)
    pre, grid = _read_named(rasters / "s2_pre.tif", _S2_BANDS, cfg.crs)
    post, _ = _read_named(rasters / "s2_post.tif", _S2_BANDS, cfg.crs, grid)
    source_path = out_dir / "label_source.tif"
    if not source_path.exists():
        raise FileNotFoundError(f"{source_path} is missing. Run label first.")
    with rasterio.open(source_path) as src:
        if src.transform != grid["transform"] or src.shape != (grid["height"], grid["width"]):
            raise RuntimeError(f"{source_path.name} is not on the s2_pre grid.")
        source = src.read(1)
    exclusions = manual_exclusions(event, event_key, cfg)
    _print_exclusions(exclusions, event)
    prov_path = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    provenance = json.loads(prov_path.read_text(encoding="utf-8"))
    provenance["manual_exclusions"] = {
        "list": exclusions,
        "confirmed_on": getattr(event, "manual_exclusions_confirmed_on", None),
        "note": "Station x side rows in these stretches stay in the table with "
                "the reason code and are left out of the primary models (14.19 e).",
    }
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"Updated {prov_path.relative_to(PROJECT_ROOT)}")
    pngs = _write_label_sheets(out_dir / "sheets", grid, pre, post, None, source,
                               samples, river_m, exclusions, combined_only=True)
    for path in pngs:
        print(path.relative_to(PROJECT_ROOT))
    return pngs


def _write_label_sheets(out_dir, grid, pre, post, labels, source, samples,
                        river_m, exclusions, combined_only: bool = False) -> list[Path]:
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Polygon

    limits = sheets.fixed_limits([pre, post], sheets.TRUE_COLOUR)
    line = [Line2D([0], [0], color="yellow", lw=1.5, label="centerline")]
    nodata = Patch(color=(0.5, 0.5, 0.5, 0.55), label="nodata")
    excluded = Patch(facecolor=(0, 0, 0, 0.25), edgecolor="black", hatch="////",
                     label="manual exclusion (side shown)")
    bars = [(item, exclusion_bar(samples, river_m, item, 40.0, 110.0))
            for item in exclusions]
    written = []

    def background(ax, win, image):
        sheets.show(ax, sheets.rgb(image, sheets.TRUE_COLOUR, limits, win), win)

    def draw_exclusions(ax, win):
        for item, polygon in bars:
            if polygon is None:
                continue
            if item["to_km"] * 1000 < win["start_m"] or item["from_km"] * 1000 > win["end_m"]:
                continue
            ax.add_patch(Polygon(polygon, closed=True, facecolor=(0, 0, 0, 0.25),
                                 edgecolor="white", hatch="////", lw=0.8))

    def image_sheet(name, title, image):
        def draw(ax, win):
            background(ax, win, image)
            draw_exclusions(ax, win)
            sheets.finish_panel(ax, win)
        return sheets.write_sheets(out_dir, name, title, grid, samples, river_m,
                                   draw, line + ([excluded] if bars else []))

    def overlay(win, values, colours):
        rows, cols = win["rows"], win["cols"]
        part = values[rows, cols]
        rgba = np.zeros(part.shape + (4,), dtype=float)
        for value, colour in colours.items():
            rgba[part == value] = colour
        return rgba

    def label_sheet(name, title, values, colours, legend):
        def draw(ax, win):
            background(ax, win, post)
            sheets.show(ax, overlay(win, values, colours), win)
            draw_exclusions(ax, win)
            sheets.finish_panel(ax, win, line_width=0.8)
        return sheets.write_sheets(out_dir, name, title, grid, samples, river_m,
                                   draw, legend + [nodata] + line
                                   + ([excluded] if bars else []))

    grey = (0.5, 0.5, 0.5, 0.55)
    if not combined_only:
        written += image_sheet("pre_truecolour", "Pre-event Sentinel-2 true colour "
                               "(same stretch as post)", pre)
        written += image_sheet("post_truecolour", "Post-event Sentinel-2 true colour "
                               "(same stretch as pre)", post)
        for name, key, colour, text in (
                ("label_optical", "optical_only", (1.0, 0.45, 0.0, 0.85), "optical impact"),
                ("label_radar", "sar_only", (0.0, 0.85, 1.0, 0.85), "radar impact")):
            coded = np.where(np.isfinite(labels[key]), labels[key], NODATA)
            written += label_sheet(
                name, f"{text.capitalize()} label over the post-event image",
                coded, {1: colour, NODATA: grey},
                [Patch(color=colour, label=text)])
    colours = {1: (1.0, 0.45, 0.0, 0.85), 2: (0.0, 0.85, 1.0, 0.85),
               3: (1.0, 0.0, 1.0, 0.9), NODATA: grey}
    written += label_sheet(
        "label_combined", "Combined label over the post-event image",
        source.astype(float), colours,
        [Patch(color=colours[1], label="optical only"),
         Patch(color=colours[2], label="radar only"),
         Patch(color=colours[3], label="both")])
    return written
