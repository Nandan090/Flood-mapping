"""The analysis table: one row per station x side (spec sections 7 and 10,
clarifications 14.5, 14.19, 14.20 and 14.21).

Built from the transects and the cached profiles. Exposure belts and controls
are added here. No model is fitted and no forest-reach relationship is shown.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .profiles import DEMS, LABELS, OUTCOME_KEYS

PRIMARY_LABEL = "combined"
OUTCOMES = ("reach_distance_m", "reach_height_m", "overtopped")
CENSOR_REASONS = ("nodata_gap", "valley_wall", "profile_end", "neighbour_limit")
TABLE_REASONS = ("hillslope_failure", "tributary_fan", "reservoir_or_dam",
                 "crosses_river", "intersects_neighbour", "tributary_confluence",
                 "outside_reach", "river_not_visible", "near_bank_nodata", "low_valid_fraction",
                 "belt_no_data")

COLUMNS = {
    "station_id": "station number, 0 at the reach start, every station_spacing_m",
    "side": "bank, left or right looking downstream",
    "river_km": "distance along the verified centerline from the start (km)",
    "segment_id": "floor(river_km), 1 km blocks for reporting",
    "block_id": "3 km bootstrap block (14.15)",
    "river_section": "0 above the first section break, 1 below it, ... (14.16)",
    "in_reach": "1 inside reach_start_km to reach_end_km (14.5)",
    "station_x": "station easting, EPSG:32645 (m)",
    "station_y": "station northing, EPSG:32645 (m)",
    "usable_length_m": "distance to the first same-bank transect crossing within 1 km; "
                       "1500 if none (14.21). Every measurement stops here",
    "usable_limited_by_neighbour": "1 if usable_length_m is set by a neighbouring transect",
    "channel_resolved": "1 if a pre-event active channel run was found (14.18)",
    "centre_offset_m": "signed distance from the station to the channel midpoint (m, left positive)",
    "channel_width_m": "length of the chosen channel run (m); 10 when unresolved",
    "bank_top_from_station_m": "bank top distance from the station along this side's normal (m)",
    "bank_top_sample_offset_m": "distance from the bank top to the first sample used (m); "
                                "5 when the channel is unresolved",
    "river_visible_post": "1 if within 30 m of the station a sample has NDVI_post < 0.30 "
                          "or is impacted in the combined label (14.20)",
    "exclude_reasons": "comma-separated reason codes; empty means the row is used",
    "analysis_set": "1 if in the reach with no reason code (primary models)",
    "reach_distance_m_{dem}": "bank top to the edge of the impacted run (m), combined label; "
                              "missing where river_not_visible",
    "reach_height_m_{dem}": "max(0, elevation at the reach edge - z0) (m), combined label",
    "overtopped_{dem}": "1 if reach_distance_m >= 20",
    "censored_{dem}": "1 if the true reach is at least the measured one",
    "censor_reason_{dem}": "nodata_gap, valley_wall (edge within 10 m of the wall, 14.23 b), "
                           "profile_end or neighbour_limit",
    "reach_height_unc_m_{dem}": "10 x tan(slope at the reach edge) (m)",
    "detached_reach_m_{dem}": "furthest connected impact sample, gaps ignored (m) (14.6)",
    "valid_fraction_{dem}": "share of samples beyond the bank top, up to the valley wall, "
                            "with label data",
    "reach_edge_from_station_m_{dem}": "reach edge distance from the station along this side (m)",
    "near_bank_nodata_{dem}": "1 if no valid sample decides the start (14.24): none within 30 m "
                              "of the bank top, or the first is a non-impacted sample beyond 20 m. "
                              "Outcomes are then missing",
    "reach_distance_m_wins_{dem}": "reach_distance_m capped at the analysis-set 97.5th percentile (14.19 f)",
    "reach_height_m_wins_{dem}": "reach_height_m capped at the analysis-set 97.5th percentile (14.19 f)",
    "z0_{dem}": "median DEM elevation over the channel run, or at the station when unresolved (m)",
    "valley_wall_found_{dem}": "1 if slope >= 30 deg for >= 50 m was found inside the usable length",
    "valley_width_m_{dem}": "bank top to the valley wall, or to the usable end when not found (m)",
    "wall_slope_deg_{dem}": "mean DEM slope over the 50 m run that defines the valley wall (deg)",
    "bank_elev_m_{dem}": "DEM elevation at the bank top sample (m)",
    "bank_slope_deg_{dem}": "mean DEM slope in the 100 m belt (deg)",
    "gradient_{dem}": "slope of z0 against river_km over +/-500 m (m per km; negative = falling)",
    "log10_upa": "log10 of the largest MERIT upstream area within 300 m of the station (km2)",
    "curvature": "signed turn of the verified line over +/-300 m (rad per km; positive = turns left)",
    "outer_bend": "1 if this bank is on the outside of the bend",
    "belt_pixels_{w}": "valid pixels in the W m belt",
    "belt_clipped_{w}": "1 if the belt was cut short at usable_length_m",
    "forest_frac_{w}": "forest share of valid belt pixels: DW trees >= 0.5, ETH height >= 5 m, "
                       "no Hansen loss in [2021, pre year] (exposure F). Where Dynamic World "
                       "has no observation, WorldCover class 10 replaces the trees test (14.23 a)",
    "forest_source_{w}": "dynamic_world, worldcover or mixed: which map decided the belt's "
                         "pixels; none when neither has data",
    "forest_source_worldcover_frac_{w}": "share of decided belt pixels taken from WorldCover",
    "forest_frac_dw_only_{w}": "forest share with Dynamic World alone (14.14)",
    "forest_definition_disagree_frac_{w}": "share of belt pixels where the two forest rules differ",
    "forest_frac_worldcover_{w}": "share of belt pixels in WorldCover class 10",
    "canopy_height_mean_{w}": "mean ETH canopy height in the belt (m)",
    "ndvi_pre_mean_{w}": "mean pre-event NDVI in the belt",
    "belt_continuity_{w}": "largest forest patch area / total forest area in the belt",
    "prior_loss_frac_{w}": "share of belt pixels with Hansen loss up to the pre-composite year",
}


def run_table(event_key: str, config_path: str | Path | None = None) -> list[Path]:
    import geopandas as gpd

    from . import controls, exposure, figures
    from .config import PROJECT_ROOT, load_config
    from .label import _read_named
    from .profiles import assert_time_separation, load_layers
    from .qa import line_samples
    from .stats import winsorise_upper

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    if not getattr(event, "manual_exclusions_confirmed_on", None):
        raise RuntimeError(
            f"manual_exclusions_confirmed_on is not set for {event_key}. The user "
            "fixes the manual exclusions before any outcome table is built (14.19 e).")
    prov_path = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    provenance = json.loads(prov_path.read_text(encoding="utf-8"))
    assert_time_separation(provenance, event)

    gpkg_in = PROJECT_ROOT / cfg.paths.outputs / event_key / "transects" / "transects.gpkg"
    cache = PROJECT_ROOT / cfg.paths.data / event_key / "tables"
    for path in (gpkg_in, cache / "profile_geometry.csv", cache / "profile_outcomes.csv"):
        if not path.exists():
            raise FileNotFoundError(f"{path} is missing. Run transects and profiles first.")
    transects = gpd.read_file(gpkg_in).sort_values(["station_id", "side"]).reset_index(drop=True)
    geometry = pd.read_csv(cache / "profile_geometry.csv")
    outcomes = pd.read_csv(cache / "profile_outcomes.csv")
    keys = ["station_id", "side"]
    df = transects.merge(geometry, on=keys, how="left", validate="one_to_one")
    if df["bank_top_from_station_m"].isna().any():
        raise RuntimeError("profile_geometry.csv does not match the transects. Re-run profiles.")
    for dem in DEMS:
        part = outcomes[(outcomes.label == PRIMARY_LABEL) & (outcomes.dem == dem)]
        part = part[keys + list(OUTCOME_KEYS)].rename(
            columns={k: f"{k}_{dem}" for k in OUTCOME_KEYS})
        df = df.merge(part, on=keys, how="left", validate="one_to_one")
    df["river_visible_post"] = df[f"river_visible_post_{PRIMARY_LABEL}"]

    coords = [np.asarray(g.coords) for g in df.geometry]
    df["station_x"] = [c[0, 0] for c in coords]
    df["station_y"] = [c[0, 1] for c in coords]
    normals = np.array([(c[1] - c[0]) / np.linalg.norm(c[1] - c[0]) for c in coords])

    # Exposure belts (7.3), cut at the usable length (14.21).
    layers, grid = load_layers(event_key, cfg)
    pre_grid = _read_named(PROJECT_ROOT / cfg.paths.data / event_key / "rasters" / "s2_pre.tif",
                           ("B4",), cfg.crs)[1]
    if pre_grid["transform"] != grid["transform"]:
        raise RuntimeError("s2_pre and s2_post are not on the same grid.")
    forest = exposure.forest_layers(event_key, cfg, provenance, grid)
    flat = {name: values.ravel() for name, values in forest.items()}
    shape = (grid["height"], grid["width"])
    start = df["bank_top_from_station_m"].to_numpy(float)
    usable = df["usable_length_m"].to_numpy(float)
    xy = df[["station_x", "station_y"]].to_numpy(float)
    along = np.column_stack([normals[:, 1], -normals[:, 0]])
    belt_index = {}
    for width in cfg.belt_widths_m:
        end = np.minimum(start + width, usable)
        polygons = [exposure.belt_polygon(xy[i], along[i], normals[i], start[i], end[i],
                                          cfg.station_spacing_m) for i in range(len(df))]
        index = exposure.belt_pixel_indices(polygons, shape, grid["transform"])
        belt_index[width] = index
        metrics = pd.DataFrame([exposure.belt_metrics(ix, flat, grid["width"]) for ix in index])
        for name in metrics.columns:
            df[f"{name}_{width}"] = metrics[name].to_numpy()
        df[f"belt_clipped_{width}"] = (start + width > usable + 1e-6).astype(int)
    for dem in DEMS:
        slope = layers[f"slope_{dem}"].ravel()
        df[f"bank_slope_deg_{dem}"] = [float(np.nanmean(slope[ix])) if ix.size and
                                       np.isfinite(slope[ix]).any() else np.nan
                                       for ix in belt_index[cfg.belt_width_m]]

    # Controls (7.4).
    river = line_samples(event, cfg)[0]
    stations = df.groupby("station_id").first().sort_index()
    for dem in DEMS:
        grad = controls.gradient_per_km(stations["river_km"].to_numpy(),
                                        stations[f"z0_{dem}"].to_numpy(), cfg.gradient_window_m)
        df[f"gradient_{dem}"] = df["station_id"].map(dict(zip(stations.index, grad)))
    curvature = controls.curvature_per_km(river, stations["river_km"].to_numpy() * 1000.0,
                                          cfg.curvature_window_m)
    df["curvature"] = df["station_id"].map(dict(zip(stations.index, curvature)))
    df["outer_bend"] = controls.outer_bend(df["curvature"].to_numpy(), df["side"].to_numpy())
    upa = _read_named(PROJECT_ROOT / cfg.paths.data / event_key / "rasters" / "merit_upa.tif",
                      ("upa",), cfg.crs, grid)[0]["upa"]
    upa_max = controls.max_within(upa, grid["transform"],
                                  stations[["station_x", "station_y"]].to_numpy(float),
                                  cfg.upa_window_m)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_upa = np.where(upa_max > 0, np.log10(upa_max), np.nan)
    df["log10_upa"] = df["station_id"].map(dict(zip(stations.index, log_upa)))

    # Reason codes (6, 7.2, 14.19 e, 14.20, 14.21) and the analysis set (8, 14.5).
    reasons = []
    for row in df.itertuples(index=False):
        found = [r for r in (row.exclude_reasons or "").split(",") if r]
        if row.river_visible_post == 0:
            found.append("river_not_visible")
        if row.near_bank_nodata_fab == 1:
            found.append("near_bank_nodata")
        if not (row.valid_fraction_fab >= cfg.valid_fraction_min):
            found.append("low_valid_fraction")
        if not np.isfinite(getattr(row, f"forest_frac_{cfg.belt_width_m}")):
            found.append("belt_no_data")
        reasons.append(",".join(found))
    df["exclude_reasons"] = reasons
    df["analysis_set"] = ((df["in_reach"] == 1) & (df["exclude_reasons"] == "")).astype(int)
    analysis = df["analysis_set"] == 1

    # Winsorised outcomes for the 14.19 (f) sensitivity run.
    # The cap is the percentile over the analysis set; other rows get the same cap.
    caps = {}
    mask = analysis.to_numpy()
    for dem in DEMS:
        for name in ("reach_distance_m", "reach_height_m"):
            values = df[f"{name}_{dem}"].to_numpy(float)
            out = values.copy()
            out[mask] = winsorise_upper(values[mask], cfg.winsorise_percentile)
            finite = values[mask & np.isfinite(values)]
            cap = float(np.percentile(finite, cfg.winsorise_percentile)) if finite.size else np.nan
            caps[f"{name}_{dem}"] = cap
            out[~mask & np.isfinite(out) & (out > cap)] = cap
            df[f"{name}_wins_{dem}"] = out

    written = _write_outputs(df, outcomes, cfg, event_key)
    _print_report(df, outcomes, cfg)

    pngs = figures.example_cross_sections(df, layers, forest, grid, cfg, event_key)
    provenance["table"] = {
        "rows": int(len(df)),
        "analysis_set_rows": int(analysis.sum()),
        "primary_label": PRIMARY_LABEL,
        "valid_fraction_min": cfg.valid_fraction_min,
        "winsorise_percentile": cfg.winsorise_percentile,
        "winsorise_caps_m": {k: (None if not np.isfinite(v) else round(v, 2))
                             for k, v in caps.items()},
        "manual_exclusions_confirmed_on": event.manual_exclusions_confirmed_on,
        "forest_fallback": (f"WorldCover {provenance['worldcover_exported']['id']} class "
                            f"{cfg.worldcover_forest_class} replaces the Dynamic World trees test "
                            "where Dynamic World has no observation in the S2 pre window, with "
                            "the same canopy height and Hansen conditions (14.23 a). See "
                            "dynamic_world_gap_check for the tile's image dates."),
        "valley_wall_censoring": "only when the reach edge is within one sample (10 m) of "
                                 "the valley wall; the walk stops at the wall (14.23 b)",
        "run_start": (f"beyond the bank top: any impacted sample within "
                      f"{cfg.near_bank_zero_m:g} m starts the run; otherwise a valid "
                      f"non-impacted sample within {cfg.near_bank_zero_m:g} m gives reach 0; "
                      f"otherwise a first valid sample that is impacted within "
                      f"{cfg.near_bank_impact_m:g} m starts the run, and any other row has "
                      "missing outcomes with near_bank_nodata; separately for every label "
                      "(14.24, revised)"),
        "example_cross_sections": "analysis-set rows nearest the 5th, 25th, 50th, 75th and "
                                  "95th percentile of reach_height_m_fab, plus the censored "
                                  "row nearest the median censored height; ties go to the "
                                  "lowest river_km, left bank first; no row is used twice",
        "files": [str(p.relative_to(PROJECT_ROOT)).replace("\\", "/") for p in written + pngs],
    }
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(f"Updated {prov_path.relative_to(PROJECT_ROOT)}")
    for path in written + pngs:
        print(path.relative_to(PROJECT_ROOT))
    return written + pngs


def _write_outputs(df, outcomes, cfg, event_key) -> list[Path]:
    import geopandas as gpd
    from shapely.geometry import Point

    from .config import PROJECT_ROOT
    out_dir = PROJECT_ROOT / cfg.paths.outputs / event_key
    # z0 and bank_top_m from the transect step are replaced by z0_fab and
    # bank_top_from_station_m.
    drop = ["status", "z0", "bank_top_m"] + [
        c for c in df.columns if c.startswith("river_visible_post_")]
    table = df.drop(columns=[c for c in drop if c in df.columns])
    first = ["station_id", "side", "river_km", "segment_id", "block_id", "river_section",
             "in_reach", "analysis_set", "exclude_reasons"]
    table = table[first + [c for c in table.columns if c not in first]]
    csv = out_dir / "transects.csv"
    pd.DataFrame(table.drop(columns="geometry")).to_csv(csv, index=False)

    labels_csv = out_dir / "transects_labels.csv"
    long = outcomes.merge(df[["station_id", "side", "river_km", "in_reach", "exclude_reasons"]],
                          on=["station_id", "side"], how="left")
    long.to_csv(labels_csv, index=False)

    gpkg = out_dir / "transects.gpkg"
    if gpkg.exists():
        gpkg.unlink()
    lines = gpd.GeoDataFrame(table, geometry="geometry", crs=cfg.crs)
    lines.to_file(gpkg, layer="transects", driver="GPKG")
    stations = df.groupby("station_id").first().reset_index()
    gpd.GeoDataFrame(
        stations[["station_id", "river_km", "segment_id", "block_id", "river_section",
                  "in_reach", "channel_resolved", "channel_width_m"]],
        geometry=[Point(x, y) for x, y in zip(stations.station_x, stations.station_y)],
        crs=cfg.crs).to_file(gpkg, layer="stations", driver="GPKG")
    edges = df[np.isfinite(df["reach_edge_from_station_m_fab"])]
    points = []
    for row in edges.itertuples(index=False):
        c = np.asarray(row.geometry.coords)
        n = (c[1] - c[0]) / np.linalg.norm(c[1] - c[0])
        points.append(Point(*(c[0] + n * row.reach_edge_from_station_m_fab)))
    gpd.GeoDataFrame(
        edges[["station_id", "side", "river_km", "analysis_set", "reach_distance_m_fab",
               "reach_height_m_fab", "censored_fab", "censor_reason_fab"]],
        geometry=points, crs=cfg.crs).to_file(gpkg, layer="reach_edges", driver="GPKG")

    columns_csv = out_dir / "transects_columns.csv"
    pd.DataFrame({"column": list(COLUMNS), "meaning": list(COLUMNS.values())}).to_csv(
        columns_csv, index=False)
    return [csv, labels_csv, gpkg, columns_csv]


def _km_ranges(km: np.ndarray, spacing_m: float) -> str:
    """Consecutive stations joined into 'a–b km' ranges."""
    if not km.size:
        return "none"
    parts, start, prev = [], km[0], km[0]
    for value in km[1:]:
        if (value - prev) * 1000.0 > spacing_m * 1.5:
            parts.append((start, prev))
            start = value
        prev = value
    parts.append((start, prev))
    return ", ".join(f"{a:.2f}" if a == b else f"{a:.2f}–{b:.2f}" for a, b in parts) + " km"


def _has(value: str, reason: str) -> bool:
    return reason in (value or "").split(",")


def _print_report(df, outcomes, cfg) -> None:
    analysis = df[df["analysis_set"] == 1]
    print("Rows by reason code. A row with more than one reason is counted in each.")
    print(f"{'reason':<24}{'left':>8}{'right':>8}{'total':>8}")
    for reason in TABLE_REASONS:
        left = int(sum(_has(v, reason) for v in df.loc[df.side == "left", "exclude_reasons"]))
        right = int(sum(_has(v, reason) for v in df.loc[df.side == "right", "exclude_reasons"]))
        print(f"{reason:<24}{left:8d}{right:8d}{left + right:8d}")
    for name, part in (("analysis set", analysis), ("all rows", df)):
        left = int((part.side == "left").sum())
        right = int((part.side == "right").sum())
        print(f"{name:<24}{left:8d}{right:8d}{left + right:8d}")

    source = f"forest_source_{cfg.belt_width_m}"
    print(f"Forest source, {cfg.belt_width_m} m belt (14.23 a):")
    print(f"{'source':<24}{'all rows':>10}{'analysis set':>14}")
    for kind in ("dynamic_world", "worldcover", "mixed", "none"):
        print(f"{kind:<24}{int((df[source] == kind).sum()):10d}"
              f"{int((analysis[source] == kind).sum()):14d}")
    uses_wc = df[df[source].isin(["worldcover", "mixed"])]
    print("River km using WorldCover (no Dynamic World image for Sentinel-2 tile T45RUM "
          "in the pre window):")
    for side in ("left", "right"):
        km = np.sort(uses_wc.loc[uses_wc.side == side, "river_km"].to_numpy())
        print(f"  {side}: {_km_ranges(km, cfg.station_spacing_m)}")

    not_visible = df.loc[df["river_visible_post"] == 0, "river_km"].drop_duplicates()
    print(f"Stations with no visible river after the event (14.20): {len(not_visible)}"
          + (f", at km {', '.join(f'{k:.1f}' for k in not_visible)}" if len(not_visible) else ""))

    near = df[df["exclude_reasons"].apply(lambda v: _has(v, "near_bank_nodata"))]
    print(f"near_bank_nodata rows (14.24), combined label, by river km: {len(near)}")
    if len(near):
        print(f"{'km':<10}{'left':>8}{'right':>8}{'total':>8}")
        km = np.floor(near["river_km"]).astype(int)
        for k in sorted(km.unique()):
            here = near[km == k]
            left = int((here.side == "left").sum())
            right = int((here.side == "right").sum())
            print(f"{f'{k}-{k + 1}':<10}{left:8d}{right:8d}{left + right:8d}")
    print("near_bank_nodata rows in the other labels (all rows):")
    for label in LABELS:
        part = outcomes[(outcomes.label == label) & (outcomes.dem == "fab")]
        print(f"  {label:<14}{int((part.near_bank_nodata == 1).sum()):6d}")

    print("Censored rows in the analysis set, combined label, by censor reason.")
    print(f"{'censor reason':<24}{'FABDEM':>8}{'GLO-30':>8}")
    for reason in CENSOR_REASONS:
        counts = [int((analysis[f"censor_reason_{dem}"] == reason).sum()) for dem in DEMS]
        print(f"{reason:<24}{counts[0]:8d}{counts[1]:8d}")
    totals = [int((analysis[f"censored_{dem}"] == 1).sum()) for dem in DEMS]
    print(f"{'all censored':<24}{totals[0]:8d}{totals[1]:8d}")

    print(f"Summary statistics, analysis set (N = {len(analysis)} rows, "
          f"{analysis.station_id.nunique()} stations), combined label.")
    print(f"{'variable':<28}{'n':>5}{'mean':>9}{'sd':>9}{'min':>9}{'q25':>9}"
          f"{'median':>9}{'q75':>9}{'max':>9}")
    columns = [f"{name}_{dem}" for dem in DEMS for name in OUTCOMES]
    columns.append(f"forest_frac_{cfg.belt_width_m}")
    for column in columns:
        values = analysis[column].to_numpy(float)
        values = values[np.isfinite(values)]
        if not values.size:
            print(f"{column:<28}{0:5d}")
            continue
        q = np.percentile(values, [0, 25, 50, 75, 100])
        print(f"{column:<28}{values.size:5d}{values.mean():9.3f}{values.std(ddof=1):9.3f}"
              f"{q[0]:9.3f}{q[1]:9.3f}{q[2]:9.3f}{q[3]:9.3f}{q[4]:9.3f}")
    print("overtopped is 0/1, so its mean is the share of rows overtopped.")
    for dem, name in (("fab", "FABDEM"), ("cop", "Copernicus GLO-30")):
        values = analysis[f"reach_distance_m_{dem}"].to_numpy(float)
        values = values[np.isfinite(values)]
        share = 100.0 * np.mean(values == 0) if values.size else float("nan")
        print(f"Share of analysis-set rows with zero reach ({name}): {share:.1f}% "
              f"of {values.size}")
