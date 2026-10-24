"""Sampling along transects and reach extraction (spec section 7.1, 7.2,
clarifications 14.6, 14.7, 14.18, 14.20 and 14.21).

Samples are spaced `step_m` apart, starting at the station. Distances are
measured from `bank_top`.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

IMPACT = 1
NO_IMPACT = 0

# Primary label first. Every label is measured with both DEMs (7.5).
LABELS = ("combined", "optical_only", "sar_only", "sar_strict", "sar_near")
DEMS = {"fab": "fabdem.tif", "cop": "copernicus_dem.tif"}
OUTCOME_KEYS = ("reach_distance_m", "reach_height_m", "overtopped", "censored",
                "censor_reason", "reach_height_unc_m", "detached_reach_m",
                "valid_fraction", "reach_edge_from_station_m", "near_bank_nodata")


def bank_top_index(active_channel: np.ndarray, run: int = 3) -> int:
    """First sample after which `run` consecutive samples are not channel."""
    ch = np.asarray(active_channel, dtype=bool)
    for i in range(len(ch) - run + 1):
        if not ch[i:i + run].any():
            return i
    raise ValueError("No bank top found: the transect never leaves the channel.")


def centre_channel(channel: np.ndarray, offsets_m: np.ndarray,
                   elevation_m: np.ndarray, step_m: float, merge_gap_m: float,
                   recentre_max_m: float) -> dict:
    """Transect centring on the pre-event active channel (14.18).

    `offsets_m` is the signed distance of each sample from the station along
    the whole transect, ascending, with the station at 0 (positive = left
    bank looking downstream). Gaps are measured as in extract_reach: N
    non-channel samples in a row span (N - 1) steps, and a gap of exactly
    merge_gap_m still merges.
    """
    ch = np.asarray(channel, dtype=bool)
    off = np.asarray(offsets_m, dtype=float)
    z = np.asarray(elevation_m, dtype=float)
    station = int(np.argmin(np.abs(off)))
    chosen = None
    hits = np.flatnonzero(ch)
    if hits.size:
        max_gap_px = int(round(merge_gap_m / step_m)) + 1
        runs = []
        lo = hi = int(hits[0])
        for i in hits[1:]:
            if i - hi - 1 <= max_gap_px:
                hi = int(i)
            else:
                runs.append((lo, hi))
                lo = hi = int(i)
        runs.append((lo, hi))
        containing = [run for run in runs if run[0] <= station <= run[1]]
        if containing:
            chosen = containing[0]
        else:
            edge = [min(abs(off[a]), abs(off[b])) for a, b in runs]
            best = int(np.argmin(edge))
            if edge[best] <= recentre_max_m:
                chosen = runs[best]
    if chosen is None:
        return {
            "channel_resolved": 0,
            "bank_top_left_m": step_m / 2,
            "bank_top_right_m": -step_m / 2,
            "channel_width_m": step_m,
            "centre_offset_m": float("nan"),
            "z0": float(z[station]),
        }
    a, b = chosen
    return {
        "channel_resolved": 1,
        "bank_top_left_m": float(off[b]),
        "bank_top_right_m": float(off[a]),
        "channel_width_m": float(off[b] - off[a] + step_m),
        "centre_offset_m": float((off[a] + off[b]) / 2),
        "z0": float(np.nanmedian(z[a:b + 1])),
    }


def _is_impact(label: np.ndarray, i: int) -> bool:
    return label[i] == IMPACT


def extract_reach(label: np.ndarray, elevation_m: np.ndarray, z0: float,
                  bank_top: int, step_m: float, slope_deg: np.ndarray,
                  valley_wall: int | None = None,
                  impact_gap_m: float = 20.0, nodata_gap_m: float = 30.0,
                  overtop_m: float = 20.0,
                  usable_index: int | None = None,
                  near_bank_impact_m: float = 30.0, near_bank_zero_m: float = 20.0,
                  start_offset_m: float = 0.0) -> dict:
    """Walk outward from the bank top and measure the impacted run.

    The first sample with data beyond the bank top decides how the row
    starts (14.24, refining 14.6). With start_offset_m = 0 the bank-top
    sample itself is skipped. Its distance from the bank top is its index times step_m plus
    `start_offset_m`. If it is impacted and within `near_bank_impact_m`, the
    run starts there. If it is not impacted and within `near_bank_zero_m`, the
    reach is 0. Otherwise the outcomes are missing and near_bank_nodata = 1.

    The run continues across non-impact gaps up to `impact_gap_m`. Nodata gaps up to
    `nodata_gap_m` do not end it; a longer one does and sets censored = 1.
    Samples past `usable_index` are not used (14.21). The walk stops at
    `valley_wall`, and the row is censored there only when the reach edge is
    within one sample of it (14.23 b). censor_reason is nodata_gap,
    valley_wall, profile_end or neighbour_limit (the run reached
    `usable_index`), and None when the row is not censored.
    """
    lab = np.asarray(label, dtype=float)
    z = np.asarray(elevation_m, dtype=float)
    slope = np.asarray(slope_deg, dtype=float)
    end = len(lab) - 1
    limit = end if usable_index is None else min(end, int(usable_index))
    wall = end if valley_wall is None else valley_wall
    fraction_end = min(wall, limit)

    finite = np.isfinite(lab[bank_top:fraction_end + 1])
    valid_fraction = float(finite.mean()) if finite.size else 0.0

    # Furthest impact that passed the 2-D connectivity filter, gaps ignored.
    impacted = np.flatnonzero(lab[bank_top:limit + 1] == IMPACT)
    detached_reach_m = float(impacted.max() * step_m) if impacted.size else 0.0

    none = {
        "reach_distance_m": 0.0,
        "reach_height_m": 0.0,
        "overtopped": 0,
        "censored": 0,
        "censor_reason": None,
        "reach_edge_index": None,
        "reach_height_unc_m": None,
        "detached_reach_m": detached_reach_m,
        "valid_fraction": valid_fraction,
        "near_bank_nodata": 0,
    }

    # Only samples beyond the bank top decide. With no offset the bank-top
    # sample is the channel edge, which some labels code as a valid 0.
    skip = 1 if start_offset_m == 0 else 0
    beyond = lab[bank_top + skip:limit + 1]
    dist_m = (np.arange(beyond.size) + skip) * step_m + start_offset_m
    near = beyond[dist_m <= near_bank_zero_m + 1e-6]
    seen = np.flatnonzero(np.isfinite(beyond))
    first_impacted = bool(seen.size) and beyond[seen[0]] == IMPACT
    if np.any(near == IMPACT):
        pass
    elif np.any(np.isfinite(near)):
        return none
    elif first_impacted and dist_m[seen[0]] <= near_bank_impact_m + 1e-6:
        pass
    else:
        nan = float("nan")
        return {**none, "reach_distance_m": nan, "reach_height_m": nan,
                "overtopped": nan, "censored": nan, "near_bank_nodata": 1}

    # Distances are measured between samples, so N bare samples in a row span
    # (N - 1) steps. A gap of exactly the tolerance is still inside the run.
    gap_px = int(round(impact_gap_m / step_m)) + 1
    nodata_px = int(round(nodata_gap_m / step_m)) + 1
    last = None
    non_impact_run = 0
    nodata_run = 0
    censored = 0

    censor_reason = None
    for i in range(bank_top, limit + 1):
        if math.isnan(lab[i]):
            nodata_run += 1
            if nodata_run > nodata_px:
                censored = 1
                censor_reason = "nodata_gap"
                break
        else:
            nodata_run = 0
            if _is_impact(lab, i):
                last = i
                non_impact_run = 0
            else:
                non_impact_run += 1
                if non_impact_run > gap_px:
                    break
        if valley_wall is not None and i >= valley_wall:
            # The walk stops at the wall. Only a run whose edge is within one
            # sample of the wall is censored there (14.23 b).
            if last is not None and valley_wall - last <= 1:
                censored = 1
                censor_reason = "valley_wall"
                last = min(last, valley_wall)
            break

    if last is None:
        return none
    if usable_index is not None and last >= int(usable_index):
        censored = 1
        censor_reason = "neighbour_limit"
        last = min(last, int(usable_index))
    elif last == end:
        censored = 1
        censor_reason = "profile_end"

    distance = (last - bank_top) * step_m
    return {
        "reach_distance_m": float(distance),
        "reach_height_m": float(max(0.0, z[last] - z0)),
        "overtopped": int(distance >= overtop_m),
        "censored": censored,
        "censor_reason": censor_reason,
        "reach_edge_index": int(last),
        "reach_height_unc_m": float(step_m * math.tan(math.radians(slope[last]))),
        "detached_reach_m": detached_reach_m,
        "valid_fraction": valid_fraction,
        "near_bank_nodata": 0,
    }


def valley_wall_index(slope_deg: np.ndarray, step_m: float, slope_min_deg: float,
                      run_m: float) -> int | None:
    """First sample beyond the bank top (index 0) where the slope stays at or
    above slope_min_deg for at least run_m (7.1). None if there is no wall.

    N samples in a row span (N - 1) steps, as for the gaps in extract_reach.
    """
    steep = np.asarray(slope_deg, dtype=float) >= slope_min_deg
    need = int(round(run_m / step_m)) + 1
    count = 0
    for i in range(1, len(steep)):
        count = count + 1 if steep[i] else 0
        if count >= need:
            return i - need + 1
    return None


def river_visible(offsets_m: np.ndarray, ndvi_post: np.ndarray, label: np.ndarray,
                  window_m: float, ndvi_max: float) -> int:
    """1 if a sample within window_m of the station has NDVI_post below
    ndvi_max or is impacted in this label (14.20), else 0."""
    near = np.abs(np.asarray(offsets_m, dtype=float)) <= window_m + 1e-6
    seen = (np.asarray(ndvi_post, dtype=float) < ndvi_max) | (np.asarray(label) == IMPACT)
    return int(np.any(near & seen))


def outward_indices(offsets_m: np.ndarray, side: str,
                    bank_top_m: float) -> tuple[np.ndarray, float]:
    """Indices of the whole-transect profile from the bank top outward on one side.

    `offsets_m` is signed (left positive). The first index is the first sample
    at or beyond the bank top. Also returns the distance from the bank top to
    that sample, which is 5 m when the channel is unresolved (14.18) and 0
    otherwise.
    """
    off = np.asarray(offsets_m, dtype=float)
    if side == "left":
        found = np.flatnonzero(off >= bank_top_m - 1e-6)
        if not found.size:
            raise ValueError("The left bank top lies beyond the usable profile.")
        start = int(found[0])
        idx = np.arange(start, len(off))
    elif side == "right":
        found = np.flatnonzero(off <= bank_top_m + 1e-6)
        if not found.size:
            raise ValueError("The right bank top lies beyond the usable profile.")
        start = int(found[-1])
        idx = np.arange(start, -1, -1)
    else:
        raise ValueError(f"side must be left or right, got {side!r}")
    return idx, float(abs(off[start] - bank_top_m))


def measure_side(label: np.ndarray, elevation_m: np.ndarray, slope_deg: np.ndarray,
                 z0: float, step_m: float, start_offset_m: float,
                 valley_wall: int | None, neighbour_limited: bool, cfg) -> dict:
    """Outcomes for one bank along an outward profile that starts at the bank top.

    The profile already stops at usable_length_m. When that limit comes from a
    neighbouring transect, a run that reaches it is censored as
    neighbour_limit (14.21); otherwise it is the profile end.
    """
    lab = np.asarray(label, dtype=float)
    end = len(lab) - 1
    res = extract_reach(lab, elevation_m, z0, 0, step_m, slope_deg,
                        valley_wall=valley_wall, impact_gap_m=cfg.impact_gap_m,
                        nodata_gap_m=cfg.nodata_gap_m, overtop_m=cfg.overtop_m,
                        usable_index=end if neighbour_limited else None,
                        near_bank_impact_m=cfg.near_bank_impact_m,
                        near_bank_zero_m=cfg.near_bank_zero_m,
                        start_offset_m=start_offset_m)
    if res["reach_edge_index"] is not None:
        res["reach_distance_m"] += start_offset_m
        res["overtopped"] = int(res["reach_distance_m"] >= cfg.overtop_m)
    impacted = np.flatnonzero(lab == IMPACT)
    res["detached_reach_m"] = (float(impacted.max() * step_m + start_offset_m)
                               if impacted.size else 0.0)
    # valid_fraction counts samples beyond the bank top. A resolved bank top
    # is a channel sample, which the label always leaves as nodata (14.19 c).
    first = 0 if start_offset_m > 0 else 1
    stop = end if valley_wall is None else min(valley_wall, end)
    span = lab[first:stop + 1]
    res["valid_fraction"] = float(np.isfinite(span).mean()) if span.size else 0.0
    return res


def transect_points(geometry, usable_length_m: float, step_m: float) -> np.ndarray:
    """Sample points every step_m from the station to usable_length_m.

    The saved transect starts at the station and runs along the normal.
    """
    coords = np.asarray(geometry.coords, dtype=float)
    normal = coords[1] - coords[0]
    normal = normal / np.linalg.norm(normal)
    count = int(np.floor(usable_length_m / step_m + 1e-6))
    return coords[0] + normal * step_m * np.arange(count + 1)[:, None]


def station_profile(rows: dict, step_m: float) -> tuple[np.ndarray, np.ndarray]:
    """Points and signed offsets of the whole transect, right end to left end."""
    left = transect_points(rows["left"].geometry, rows["left"].usable_length_m, step_m)
    right = transect_points(rows["right"].geometry, rows["right"].usable_length_m, step_m)
    points = np.vstack([right[::-1], left[1:]])
    offsets = np.concatenate([-np.arange(len(right))[::-1] * step_m,
                              np.arange(1, len(left)) * step_m])
    return points, offsets


def assert_time_separation(provenance: dict, event) -> None:
    """Rule 0.4: exposure scenes before the event, outcome scenes after it."""
    from .dates import assert_exposure_is_pre_event, assert_outcome_is_post_event
    composites = provenance.get("composites", {})
    checks = (("s2_pre", "exposure"), ("dynamic_world_pre", "exposure"),
              ("s2_post", "outcome"), ("s1_post", "outcome"))
    for name, role in checks:
        if name not in composites or not composites[name].get("scenes"):
            raise RuntimeError(
                f"provenance.json has no scene list for {name}, so the "
                "pre/post dates cannot be checked. Run fetch first.")
        for scene in composites[name]["scenes"]:
            what = f"{name} {scene['scene_id']}"
            if role == "exposure":
                assert_exposure_is_pre_event(scene["time_utc"], event.pre_end_utc, what)
            else:
                assert_outcome_is_post_event(scene["time_utc"], event.post_start_utc, what)


def load_layers(event_key: str, cfg) -> tuple[dict, dict]:
    """Rasters sampled along the transects, all on the analysis grid."""
    from . import optical
    from .config import PROJECT_ROOT
    from .label import _read_named
    from .sar import terrain_slope_aspect

    rasters = PROJECT_ROOT / cfg.paths.data / event_key / "rasters"
    label_dir = PROJECT_ROOT / cfg.paths.outputs / event_key / "label"
    post, grid = _read_named(rasters / "s2_post.tif", ("B4", "B8"), cfg.crs)
    layers = {"ndvi_post": optical.ndvi(post["B8"], post["B4"])}
    for dem, name in DEMS.items():
        z = _read_named(rasters / name, ("elevation_m",), cfg.crs, grid)[0]["elevation_m"]
        layers[f"z_{dem}"] = z
        layers[f"slope_{dem}"] = terrain_slope_aspect(z, grid["pixel_m"])[0]
    channel = _read_named(label_dir / "active_channel.tif", ("active_channel",),
                          cfg.crs, grid)[0]["active_channel"]
    layers["channel"] = (channel == 1).astype(float)
    for label in LABELS:
        band = f"impact_{label}"
        layers[label] = _read_named(label_dir / f"{band}.tif", (band,), cfg.crs, grid)[0][band]
    return layers, grid


def run_profiles(event_key: str, config_path: str | Path | None = None) -> list[Path]:
    """Sample every label and both DEMs along each transect and measure the reach."""
    import geopandas as gpd
    import pandas as pd

    from .centerline import _sample
    from .config import PROJECT_ROOT, load_config

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    prov_path = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    if not prov_path.exists():
        raise FileNotFoundError(f"{prov_path} is missing. Run fetch first.")
    provenance = json.loads(prov_path.read_text(encoding="utf-8"))
    assert_time_separation(provenance, event)

    gpkg = PROJECT_ROOT / cfg.paths.outputs / event_key / "transects" / "transects.gpkg"
    if not gpkg.exists():
        raise FileNotFoundError(f"{gpkg} is missing. Run transects first.")
    transects = gpd.read_file(gpkg)
    if transects.crs is None or transects.crs.to_string() != cfg.crs:
        raise RuntimeError(f"{gpkg.name} is not in {cfg.crs}.")
    layers, grid = load_layers(event_key, cfg)
    step = cfg.grid_m

    geometry_rows, outcome_rows = [], []
    for station_id, pair in transects.groupby("station_id", sort=True):
        rows = {r.side: r for r in pair.itertuples(index=False)}
        if set(rows) != {"left", "right"}:
            raise RuntimeError(f"Station {station_id} does not have both banks.")
        points, offsets = station_profile(rows, step)
        sampled = {name: _sample(values, grid, points) for name, values in layers.items()}
        channel = sampled["channel"] == 1
        centred = {dem: centre_channel(channel, offsets, sampled[f"z_{dem}"], step,
                                       cfg.channel_merge_gap_m, cfg.recentre_max_m)
                   for dem in DEMS}
        for side in ("left", "right"):
            saved = float(rows[side].bank_top_m)
            if not math.isclose(centred["fab"][f"bank_top_{side}_m"], saved, abs_tol=1e-6):
                raise RuntimeError(
                    f"Station {station_id} {side}: bank top {centred['fab'][f'bank_top_{side}_m']} "
                    f"differs from the transect step ({saved}). Re-run transects.")
        visible = {label: river_visible(offsets, sampled["ndvi_post"], sampled[label],
                                        cfg.river_visible_window_m, cfg.river_visible_ndvi_max)
                   for label in LABELS}

        for side in ("left", "right"):
            row = rows[side]
            sign = 1.0 if side == "left" else -1.0
            bank = centred["fab"][f"bank_top_{side}_m"]
            idx, start_offset = outward_indices(offsets, side, bank)
            limited = float(row.usable_length_m) < cfg.max_half_length_m - 1e-6
            geo = {
                "station_id": int(station_id),
                "side": side,
                "bank_top_from_station_m": sign * bank,
                "bank_top_sample_offset_m": start_offset,
                "profile_end_from_station_m": float(sign * offsets[idx[-1]]),
                "usable_limited_by_neighbour": int(limited),
                **{f"river_visible_post_{label}": visible[label] for label in LABELS},
            }
            for dem in DEMS:
                z = sampled[f"z_{dem}"][idx]
                slope = sampled[f"slope_{dem}"][idx]
                wall = valley_wall_index(slope, step, cfg.wall_slope_deg, cfg.wall_run_m)
                run_px = int(round(cfg.wall_run_m / step)) + 1
                geo[f"z0_{dem}"] = centred[dem]["z0"]
                geo[f"valley_wall_found_{dem}"] = int(wall is not None)
                geo[f"valley_width_m_{dem}"] = ((wall if wall is not None else len(idx) - 1)
                                                * step + start_offset)
                geo[f"wall_slope_deg_{dem}"] = (float(np.nanmean(slope[wall:wall + run_px]))
                                                if wall is not None else float("nan"))
                geo[f"bank_elev_m_{dem}"] = float(z[0])
                for label in LABELS:
                    res = measure_side(sampled[label][idx], z, slope, centred[dem]["z0"],
                                       step, start_offset, wall, limited, cfg)
                    edge = res.pop("reach_edge_index")
                    res["reach_edge_from_station_m"] = (float(sign * offsets[idx[edge]])
                                                        if edge is not None else float("nan"))
                    if not visible[label]:
                        for key in OUTCOME_KEYS:
                            if key not in ("valid_fraction", "near_bank_nodata"):
                                res[key] = None if key == "censor_reason" else float("nan")
                    outcome_rows.append({"station_id": int(station_id), "side": side,
                                         "label": label, "dem": dem,
                                         "river_visible_post": visible[label],
                                         **{key: res[key] for key in OUTCOME_KEYS}})
            geometry_rows.append(geo)

    out_dir = PROJECT_ROOT / cfg.paths.data / event_key / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    geometry_csv = out_dir / "profile_geometry.csv"
    outcomes_csv = out_dir / "profile_outcomes.csv"
    pd.DataFrame(geometry_rows).to_csv(geometry_csv, index=False)
    pd.DataFrame(outcome_rows).to_csv(outcomes_csv, index=False)

    provenance["profiles"] = {
        "labels": list(LABELS),
        "dems": {"fab": "FABDEM (primary)", "cop": "Copernicus GLO-30 2024 (sensitivity)"},
        "samples": f"every {step:g} m, nearest pixel, only within usable_length_m (14.21)",
        "bank_top": "transect centring on the pre-event active channel (14.18)",
        "censor_reasons": ["nodata_gap", "valley_wall", "profile_end", "neighbour_limit"],
        "river_visible_post": (f"a sample within {cfg.river_visible_window_m:g} m of the station "
                               f"with NDVI_post < {cfg.river_visible_ndvi_max:g} or impacted "
                               "in the label used (14.20)"),
        "files": [str(p.relative_to(PROJECT_ROOT)).replace("\\", "/")
                  for p in (geometry_csv, outcomes_csv)],
    }
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"Profiles: {len(geometry_rows)} station x side rows, labels {', '.join(LABELS)}, "
          "DEMs FABDEM and Copernicus GLO-30.")
    print(f"Updated {prov_path.relative_to(PROJECT_ROOT)}")
    for path in (geometry_csv, outcomes_csv):
        print(path.relative_to(PROJECT_ROOT))
    return [geometry_csv, outcomes_csv]
