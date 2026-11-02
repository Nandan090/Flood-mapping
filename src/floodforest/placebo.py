"""Placebo (no-flood control) run and report (spec 8.9, 14.2).

The placebo is the real event's pipeline run unchanged on two pre-event years
with no notable flood: the derived event `<event>_placebo` (config.py) has the
placebo image windows and the real event's reach, exclusions and thresholds.
Its files go to data/<event>_placebo and outputs/<event>_placebo only.

The report compares the placebo with the real event (read only). It sorts the
impact the placebo label detects into diagnostic categories for review; those
categories never decide what is impact. Agricultural use and construction
cannot be told apart at 10 m except where a patch was reviewed by eye.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .config import PROJECT_ROOT, load_config, placebo_key

PRIMARY_DEM = "fab"
CATEGORIES = {
    "channel_margin": "Channel margin: ordinary seasonal river change or possible flood-related change",
    "floodplain_terrace": "Floodplain or terrace: agricultural, construction or other land use "
                          "(not separable at 10 m)",
    "hillslope": "Hillslope: landslide or gully (not river flooding)",
}
CATEGORY_COLOURS = {"channel_margin": "#e8412c", "floodplain_terrace": "#f2c12e",
                    "hillslope": "#9b59d0"}


# ---------------------------------------------------------------------------
# Running the pipeline
# ---------------------------------------------------------------------------

def run_placebo(event_key: str, config_path: str | Path | None = None,
                refetch: bool = False) -> list[Path]:
    """Run the placebo for one event and write its report."""
    from .centerline import run_transects
    from .ee_data import fetch_event, require_confirmed_flood_check
    from .label import run_label
    from .profiles import run_profiles
    from .table import run_table

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    key = resolve_placebo(cfg, event_key)
    event = getattr(cfg.events, key)
    require_confirmed_flood_check(event, key)
    print(f"Placebo {key}: fake event {event.event_dates[0]}, confirmed flood-free "
          f"{event.flood_check.interval.start} to {event.flood_check.interval.end}.")

    full, partial = fetch_state(cfg, key)
    if refetch or not full:
        fetch_event(key, config_path)
        full, partial = True, False
    if not partial:
        # The label also needs the sar_near baseline and the radar geometry, which
        # the full fetch leaves to a partial fetch (as for the real event).
        fetch_event(key, config_path, only=list(PARTIAL_PRODUCTS))
    if not refetch:
        print(f"Cached placebo rasters in {cfg.paths.data}/{key}/rasters are reused; "
              "pass --refetch to download them again.")
    for step in (run_label, run_transects, run_profiles, run_table):
        step(key, config_path)
    return write_placebo_report(event_key, config_path)


def resolve_placebo(cfg, event_key: str) -> str:
    if event_key.endswith("_placebo"):
        raise RuntimeError("Pass the real event key (for example melamchi2021); the placebo "
                           "key is derived from it.")
    if not hasattr(cfg.events, event_key):
        raise RuntimeError(f"Unknown event {event_key}.")
    key = placebo_key(event_key)
    if not hasattr(cfg.events, key):
        raise RuntimeError(f"{event_key} has no placebo block in config.yaml.")
    return key


PARTIAL_PRODUCTS = ("s1_near", "s1_geometry")


def fetch_state(cfg, key: str) -> tuple[bool, bool]:
    """(full fetch finished, partial fetch of PARTIAL_PRODUCTS finished) for this placebo."""
    prov = PROJECT_ROOT / cfg.paths.outputs / key / "provenance.json"
    if not prov.exists():
        return False, False
    provenance = json.loads(prov.read_text(encoding="utf-8"))
    files = provenance.get("files", {})
    rasters = PROJECT_ROOT / cfg.paths.data / key / "rasters"
    full = ("placebo" in provenance and bool(files)
            and all((PROJECT_ROOT / f).exists() for f in files.values()))
    partial = (full and "s1_near" in provenance.get("composites", {})
               and "s1_geometry" in provenance
               and (rasters / "s1_near.tif").exists() and (rasters / "s1_angle.tif").exists())
    return full, partial


# ---------------------------------------------------------------------------
# Pure helpers (tested on synthetic data)
# ---------------------------------------------------------------------------

def outcome_summary(table: pd.DataFrame, dem: str = PRIMARY_DEM) -> dict:
    """Outcomes over the analysis set of one transects table."""
    a = table[table["analysis_set"] == 1]
    reach = a[f"reach_distance_m_{dem}"].to_numpy(float)
    height = a[f"reach_height_m_{dem}"].to_numpy(float)
    ok = np.isfinite(reach)
    reasons = table["exclude_reasons"].fillna("")
    return {
        "rows": int(len(table)),
        "analysis_rows": int(len(a)),
        "analysis_stations": int(a["station_id"].nunique()),
        "share_reach_above_0": float(np.mean(reach[ok] > 0)) if ok.any() else np.nan,
        "share_overtopped": float(np.nanmean(a[f"overtopped_{dem}"].to_numpy(float)))
        if len(a) else np.nan,
        "mean_reach_distance_m": float(np.nanmean(reach)) if ok.any() else np.nan,
        "median_reach_distance_m": float(np.nanmedian(reach)) if ok.any() else np.nan,
        "p90_reach_distance_m": float(np.nanpercentile(reach, 90)) if ok.any() else np.nan,
        "mean_reach_height_m": float(np.nanmean(height)) if np.isfinite(height).any() else np.nan,
        "share_censored": float(np.nanmean(a[f"censored_{dem}"].to_numpy(float)))
        if len(a) else np.nan,
        "rows_near_bank_nodata": int(reasons.str.contains("near_bank_nodata").sum()),
        "rows_river_not_visible": int(reasons.str.contains("river_not_visible").sum()),
        "rows_low_valid_fraction": int(reasons.str.contains("low_valid_fraction").sum()),
    }


def per_km_summary(table: pd.DataFrame, dem: str = PRIMARY_DEM) -> pd.DataFrame:
    a = table[table["analysis_set"] == 1]
    g = a.groupby(["segment_id", "side"])
    return pd.DataFrame({
        "rows": g.size(),
        "share_overtopped": g[f"overtopped_{dem}"].mean(),
        "mean_reach_distance_m": g[f"reach_distance_m_{dem}"].mean(),
    }).reset_index()


def label_outcomes(labels: pd.DataFrame, table: pd.DataFrame, dem: str = PRIMARY_DEM) -> pd.DataFrame:
    """Share overtopped and mean reach per label variant, over the table's analysis set."""
    keep = table.loc[table["analysis_set"] == 1, ["station_id", "side"]]
    part = labels[labels["dem"] == dem].merge(keep, on=["station_id", "side"])
    g = part.groupby("label")
    return pd.DataFrame({"rows": g.size(),
                         "share_overtopped": g["overtopped"].mean(),
                         "mean_reach_distance_m": g["reach_distance_m"].mean()}).reset_index()


def classify_detection(height_above_channel_m: float, distance_to_channel_m: float,
                       diag) -> str:
    """Diagnostic category of one detected patch (never decides what is impact)."""
    if np.isfinite(height_above_channel_m) and height_above_channel_m > diag.hillslope_height_m:
        return "hillslope"
    if (distance_to_channel_m <= diag.channel_margin_distance_m
            and (not np.isfinite(height_above_channel_m)
                 or abs(height_above_channel_m) <= diag.channel_margin_height_m)):
        return "channel_margin"
    return "floodplain_terrace"


def detection_table(impact: np.ndarray, source: np.ndarray, channel: np.ndarray,
                    dem: np.ndarray, transform, pixel_m: float, samples: np.ndarray,
                    river_m: np.ndarray, diag) -> tuple[pd.DataFrame, np.ndarray]:
    """Connected patches of impact == 1, with position, terrain and category.

    Returns the table and a raster of patch ids (0 = none; small patches keep an id
    but are not listed)."""
    from rasterio.transform import xy
    from scipy import ndimage
    from scipy.spatial import cKDTree

    hit = np.asarray(impact) == 1
    ids, _ = ndimage.label(hit, structure=np.ones((3, 3), int))
    channel = np.asarray(channel, bool)
    dist = ndimage.distance_transform_edt(~channel) * pixel_m if channel.any() \
        else np.full(hit.shape, np.inf)
    ch_r, ch_c = np.nonzero(channel)
    ch_x, ch_y = (np.asarray(v, float) for v in xy(transform, ch_r, ch_c)) if ch_r.size \
        else (np.array([]), np.array([]))
    ch_z = dem[ch_r, ch_c]
    tree = cKDTree(samples)
    rows = []
    for idx, sl in enumerate(ndimage.find_objects(ids), start=1):
        if sl is None:
            continue
        mask = ids[sl] == idx
        px = int(mask.sum())
        if px < diag.min_patch_px:
            continue
        rr, cc = np.nonzero(mask)
        rr, cc = rr + sl[0].start, cc + sl[1].start
        xs, ys = (np.asarray(v, float) for v in xy(transform, rr, cc))
        cx, cy = float(xs.mean()), float(ys.mean())
        near = np.hypot(ch_x - cx, ch_y - cy) <= diag.channel_search_m
        height = (float(np.nanmedian(dem[rr, cc]) - np.nanmedian(ch_z[near]))
                  if near.any() else float("nan"))
        d_ch = float(dist[rr, cc].min())
        i = int(tree.query([cx, cy])[1])
        j = min(i + 1, len(samples) - 1)
        dx, dy = samples[j] - samples[max(j - 1, 0)]
        cross = dx * (cy - samples[i][1]) - dy * (cx - samples[i][0])
        src = source[rr, cc]
        rows.append({
            "patch_id": idx, "river_km": round(float(river_m[i]) / 1000, 2),
            "side": "left" if cross > 0 else "right", "area_m2": px * pixel_m * pixel_m,
            "distance_to_channel_m": round(d_ch, 1),
            "height_above_channel_m": round(height, 1),
            "optical_share": round(float(np.isin(src, (1, 3)).mean()), 2),
            "radar_share": round(float(np.isin(src, (2, 3)).mean()), 2),
            "x": round(cx, 1), "y": round(cy, 1),
            "category": classify_detection(height, d_ch, diag),
        })
    cols = ["patch_id", "river_km", "side", "area_m2", "distance_to_channel_m",
            "height_above_channel_m", "optical_share", "radar_share", "x", "y", "category"]
    return pd.DataFrame(rows, columns=cols), ids


def attach_visual_review(det: pd.DataFrame, review: pd.DataFrame, match_m: float) -> pd.DataFrame:
    """Carry the nearest visual review (within match_m) onto each detection."""
    out = det.copy()
    out["visual_class"] = ""
    out["visual_reviewer"] = ""
    out["visual_note"] = ""
    if review.empty or out.empty:
        return out
    rx, ry = review["x"].to_numpy(float), review["y"].to_numpy(float)
    for i, r in out.iterrows():
        d = np.hypot(rx - r.x, ry - r.y)
        j = int(np.argmin(d))
        if d[j] <= match_m:
            rv = review.iloc[j]
            out.at[i, "visual_class"] = rv["visual_class"]
            out.at[i, "visual_reviewer"] = rv["reviewer"]
            out.at[i, "visual_note"] = f"{rv['check']} patch {rv['patch']}: {rv['note']}"
    return out


def forest_terciles(values: pd.Series) -> pd.Series:
    """Low/Mid/High equal-count groups of forest fraction (ties by order of rank)."""
    ok = values.notna()
    out = pd.Series(pd.NA, index=values.index, dtype="object")
    if ok.sum() >= 3:
        ranks = values[ok].rank(method="first")
        out[ok] = pd.qcut(ranks, 3, labels=["Low", "Mid", "High"]).astype(str)
    return out


def tercile_table(table: pd.DataFrame, forest_col: str, dem: str = PRIMARY_DEM) -> pd.DataFrame:
    a = table[table["analysis_set"] == 1].copy()
    a["tercile"] = forest_terciles(a[forest_col])
    g = a.dropna(subset=["tercile"]).groupby("tercile")
    out = pd.DataFrame({
        "rows": g.size(),
        "forest_frac_min": g[forest_col].min(), "forest_frac_max": g[forest_col].max(),
        "overtopped_rows": g[f"overtopped_{dem}"].sum(),
        "share_overtopped": g[f"overtopped_{dem}"].mean(),
        "mean_reach_distance_m": g[f"reach_distance_m_{dem}"].mean(),
    })
    return out.reindex(["Low", "Mid", "High"]).reset_index()


def river_not_visible_as_zero(table: pd.DataFrame, dem: str = PRIMARY_DEM) -> pd.DataFrame:
    """The 14.20 sensitivity: in-reach rows left out only for river_not_visible join the
    analysis set with reach 0."""
    out = table.copy()
    only = (out["in_reach"] == 1) & (out["exclude_reasons"].fillna("") == "river_not_visible")
    out.loc[only, "analysis_set"] = 1
    for name in (f"reach_distance_m_{dem}", f"reach_height_m_{dem}", f"overtopped_{dem}",
                 f"censored_{dem}"):
        out.loc[only, name] = 0
    return out


def standin_flag(table: pd.DataFrame, width: int) -> pd.Series:
    """1 where any decided belt pixel came from the WorldCover stand-in."""
    share = pd.to_numeric(table[f"forest_source_worldcover_frac_{width}"], errors="coerce")
    return (share.fillna(0) > 0).astype(int)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_placebo_report(event_key: str, config_path: str | Path | None = None) -> list[Path]:
    import geopandas as gpd

    from .label import _read_named
    from .qa import line_samples
    from .table import _km_ranges

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    key = resolve_placebo(cfg, event_key)
    event = getattr(cfg.events, key)
    diag = cfg.placebo_diagnostics
    width = int(cfg.belt_width_m)
    forest_col = f"forest_frac_{width}"
    real_out = PROJECT_ROOT / cfg.paths.outputs / event_key
    plac_out = PROJECT_ROOT / cfg.paths.outputs / key
    out_dir = plac_out / "placebo_report"
    if real_out.resolve() in out_dir.resolve().parents:
        raise RuntimeError("The placebo report would be written inside the real event's outputs.")
    out_dir.mkdir(parents=True, exist_ok=True)

    plac = pd.read_csv(plac_out / "transects.csv")
    real = pd.read_csv(real_out / "transects.csv")
    plac_labels = pd.read_csv(plac_out / "transects_labels.csv")
    real_labels = pd.read_csv(real_out / "transects_labels.csv")
    provenance = json.loads((plac_out / "provenance.json").read_text(encoding="utf-8"))

    # Rows: flags, the 2020-forest sensitivity column (read from the real table).
    plac["worldcover_standin"] = standin_flag(plac, width)
    real_forest = real[["station_id", "side", forest_col]].rename(
        columns={forest_col: f"{forest_col}_2020"})
    plac = plac.merge(real_forest, on=["station_id", "side"], how="left", validate="one_to_one")

    # Detections in the placebo combined label.
    label_dir = plac_out / "label"
    impact, grid = _read_named(label_dir / "impact_combined.tif", ("impact_combined",), cfg.crs)
    impact = impact["impact_combined"]
    source = _read_first_band(label_dir / "label_source.tif", grid, cfg.crs)
    channel = _read_named(label_dir / "active_channel.tif", ("active_channel",), cfg.crs, grid)[0][
        "active_channel"] == 1
    rasters = PROJECT_ROOT / cfg.paths.data / key / "rasters"
    dem = _read_named(rasters / "fabdem.tif", ("elevation_m",), cfg.crs, grid)[0]["elevation_m"]
    line, samples, river_m = line_samples(event, cfg)
    det, ids = detection_table(impact, source, channel, dem, grid["transform"], grid["pixel_m"],
                               samples, river_m, diag)
    review_path = PROJECT_ROOT / diag.visual_review_file
    review = pd.read_csv(review_path) if review_path.exists() else pd.DataFrame()
    det = attach_visual_review(det, review, diag.review_match_m)
    from pyproj import Transformer
    to_ll = Transformer.from_crs(cfg.crs, "EPSG:4326", always_xy=True)
    if len(det):
        det["lon"], det["lat"] = to_ll.transform(det.x.to_numpy(), det.y.to_numpy())
    det["category_label"] = det["category"].map(CATEGORIES)

    # Which patch sets each analysis-row reach.
    edges = gpd.read_file(plac_out / "transects.gpkg", layer="reach_edges")
    from rasterio.transform import rowcol
    r, c = rowcol(grid["transform"], edges.geometry.x.to_numpy(), edges.geometry.y.to_numpy())
    r, c = np.asarray(r), np.asarray(c)
    inside = (r >= 0) & (c >= 0) & (r < grid["height"]) & (c < grid["width"])
    pid = np.zeros(len(edges), int)
    pid[inside] = ids[r[inside], c[inside]]
    edges["patch_id"] = pid
    cat = dict(zip(det.patch_id, det.category))
    edges["edge_category"] = [cat.get(p, "small_patch" if p else "") for p in pid]
    plac = plac.merge(edges[["station_id", "side", "patch_id", "edge_category"]],
                      on=["station_id", "side"], how="left")
    plac["edge_category"] = plac["edge_category"].fillna("")
    reach_rows = plac[(plac.analysis_set == 1) & (plac[f"reach_distance_m_{PRIMARY_DEM}"] > 0)]

    # Summaries.
    s_plac, s_real = outcome_summary(plac), outcome_summary(real)
    flagged = plac[plac.worldcover_standin == 1]
    summary = pd.DataFrame({"placebo": s_plac, "real_event": s_real})
    sens_rows = {
        "placebo, all analysis rows": plac,
        "placebo, without WorldCover stand-in rows": plac[plac.worldcover_standin == 0],
        "placebo, river-not-visible rows as reach 0 (14.20)": river_not_visible_as_zero(plac),
        "real event, river-not-visible rows as reach 0 (14.20)": river_not_visible_as_zero(real),
    }
    sens = pd.DataFrame({k: outcome_summary(v) for k, v in sens_rows.items()}).T
    terc = {
        "2018 forest (primary)": tercile_table(plac, forest_col),
        "2018 forest, without stand-in rows": tercile_table(plac[plac.worldcover_standin == 0],
                                                            forest_col),
        "2020 forest (sensitivity)": tercile_table(plac, f"{forest_col}_2020"),
        "real event, 2020 forest": tercile_table(real, forest_col),
    }
    a_plac = plac[plac.analysis_set == 1].copy()
    a_plac["tercile"] = forest_terciles(a_plac[forest_col])
    standin_by_tercile = (a_plac[a_plac.worldcover_standin == 1].groupby("tercile").size()
                          .reindex(["Low", "Mid", "High"]).fillna(0).astype(int))
    km_plac, km_real = per_km_summary(plac), per_km_summary(real)
    per_km = km_plac.merge(km_real, on=["segment_id", "side"], how="outer",
                           suffixes=("_placebo", "_real"))
    lab = label_outcomes(plac_labels, plac).merge(label_outcomes(real_labels, real), on="label",
                                                  how="outer", suffixes=("_placebo", "_real"))
    det_by_cat = (det.groupby("category").agg(patches=("patch_id", "size"),
                                              area_ha=("area_m2", lambda a: a.sum() / 1e4))
                  .reindex(list(CATEGORIES)).fillna(0))
    real_impact = _read_named(real_out / "label" / "impact_combined.tif", ("impact_combined",),
                              cfg.crs)[0]["impact_combined"]
    px_ha = grid["pixel_m"] ** 2 / 1e4
    area = {"placebo_impact_ha": float((impact == 1).sum() * px_ha),
            "real_event_impact_ha": float((real_impact == 1).sum() * px_ha)}

    # Files.
    keep = ["station_id", "side", "river_km", "segment_id", "in_reach", "analysis_set",
            "exclude_reasons", f"reach_distance_m_{PRIMARY_DEM}", f"reach_height_m_{PRIMARY_DEM}",
            f"overtopped_{PRIMARY_DEM}", f"censored_{PRIMARY_DEM}", forest_col,
            f"forest_source_{width}", f"forest_source_worldcover_frac_{width}",
            "worldcover_standin", f"{forest_col}_2020", "patch_id", "edge_category"]
    written = []

    def save(frame, name, index=False):
        path = out_dir / name
        frame.to_csv(path, index=index)
        written.append(path)

    save(plac[keep], "placebo_rows.csv")
    save(det.drop(columns=["x", "y"]), "placebo_detections.csv")
    save(summary, "placebo_vs_real_summary.csv", index=True)
    save(sens, "placebo_standin_sensitivity.csv", index=True)
    save(pd.concat({k: v for k, v in terc.items()}, names=["forest", "i"]).reset_index(level=0),
         "placebo_forest_terciles.csv")
    save(per_km, "placebo_vs_real_per_km.csv")
    save(lab, "placebo_vs_real_by_label.csv")
    written += _figures(out_dir, km_plac, km_real, plac, real, det, ids, grid, event, cfg, samples,
                        river_m, rasters)

    text = _report_text(event_key, key, event, cfg, provenance, summary, sens, terc,
                        standin_by_tercile, flagged, plac, real, det, det_by_cat, reach_rows, lab,
                        area, _km_ranges, written, out_dir)
    report = out_dir / "PLACEBO_REPORT.md"
    report.write_text(text, encoding="utf-8")
    written.append(report)
    provenance["placebo_report"] = {
        "files": [str(p.relative_to(PROJECT_ROOT)).replace("\\", "/") for p in written],
        "diagnostic_rules": vars(diag),
    }
    (plac_out / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    for path in written:
        print(path.relative_to(PROJECT_ROOT))
    return written


def _read_first_band(path: Path, grid: dict, crs: str) -> np.ndarray:
    """Band 1 of a raster on the given grid (label_source.tif names its band by its legend)."""
    import rasterio
    with rasterio.open(path) as src:
        if src.crs is None or src.crs.to_string() != crs:
            raise RuntimeError(f"{path.name} is not in {crs}.")
        if (src.transform != grid["transform"] or src.height != grid["height"]
                or src.width != grid["width"]):
            raise RuntimeError(f"{path.name} is not on the label grid.")
        return src.read(1)


def _figures(out_dir, km_plac, km_real, plac, real, det, ids, grid, event, cfg, samples, river_m,
             rasters) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    from . import sheets
    from .label import _read_named

    written = []
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, dpi=120)
    for ax, side in zip(axes, ("left", "right")):
        for frame, name, colour in ((km_real, "real event 2021", "#2b6cb0"),
                                    (km_plac, "placebo 2018-19", "#e8412c")):
            part = frame[frame.side == side].sort_values("segment_id")
            ax.step(part.segment_id + 0.5, part.share_overtopped, where="mid", label=name,
                    color=colour)
        ax.set_ylabel(f"{side} bank\nshare overtopped")
        ax.set_ylim(-0.02, 1.02)
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set_xlabel("river km")
    fig.suptitle("Share of analysis-set rows with reach >= 20 m, per river km (combined label, FABDEM)",
                 fontsize=9)
    fig.tight_layout()
    path = out_dir / "fig_along_river.png"
    fig.savefig(path)
    plt.close(fig)
    written.append(path)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
    for frame, name, colour in ((real, "real event 2021", "#2b6cb0"),
                                (plac, "placebo 2018-19", "#e8412c")):
        v = frame.loc[frame.analysis_set == 1, f"reach_distance_m_{PRIMARY_DEM}"].dropna()
        v = np.sort(v.to_numpy(float))
        if v.size:
            ax.step(v, np.arange(1, v.size + 1) / v.size, where="post", label=f"{name} (n={v.size})",
                    color=colour)
    ax.set_xlabel("reach distance (m)")
    ax.set_ylabel("share of analysis rows at or below")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    path = out_dir / "fig_reach_distribution.png"
    fig.savefig(path)
    plt.close(fig)
    written.append(path)

    fig, ax = plt.subplots(figsize=(10, 3.5), dpi=120)
    if len(det):
        km = np.floor(det.river_km).astype(int)
        table = (det.assign(km=km).pivot_table(index="km", columns="category", values="area_m2",
                                               aggfunc="sum", fill_value=0) / 1e4)
        bottom = np.zeros(len(table))
        for catg in CATEGORIES:
            if catg in table:
                ax.bar(table.index + 0.5, table[catg], bottom=bottom, width=0.9,
                       color=CATEGORY_COLOURS[catg], label=catg.replace("_", " "))
                bottom += table[catg].to_numpy()
    ax.set_xlabel("river km")
    ax.set_ylabel("detected impact (ha)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Placebo detections by diagnostic category (patches >= "
                 f"{cfg.placebo_diagnostics.min_patch_px} px)", fontsize=9)
    fig.tight_layout()
    path = out_dir / "fig_detections_by_km.png"
    fig.savefig(path)
    plt.close(fig)
    written.append(path)

    post, _ = _read_named(rasters / "s2_post.tif", sheets.TRUE_COLOUR, cfg.crs, grid)
    limits = sheets.fixed_limits([post], sheets.TRUE_COLOUR)
    code = np.zeros(ids.shape, int)
    for n, catg in enumerate(CATEGORIES, start=1):
        code[np.isin(ids, det.loc[det.category == catg, "patch_id"].to_numpy())] = n
    colours = [matplotlib.colors.to_rgba(CATEGORY_COLOURS[c]) for c in CATEGORIES]

    def draw(ax, win):
        sheets.show(ax, sheets.rgb(post, sheets.TRUE_COLOUR, limits, win), win)
        sub = code[win["rows"], win["cols"]]
        overlay = np.zeros(sub.shape + (4,))
        for n, colour in enumerate(colours, start=1):
            overlay[sub == n] = colour
        ax.imshow(overlay, extent=win["extent"], origin="upper", interpolation="nearest")
        sheets.finish_panel(ax, win, line_colour="cyan", line_width=0.6)

    legend = [Patch(color=CATEGORY_COLOURS[c], label=c.replace("_", " ")) for c in CATEGORIES]
    written += sheets.write_sheets(out_dir, "sheet_placebo_detections",
                                   "Placebo combined-label detections on Oct-Dec 2019 Sentinel-2",
                                   grid, samples, river_m, draw, legend)
    return written


def _fmt(v, d=2):
    return "NA" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{d}f}"


def _pct(v):
    return "NA" if v is None or not np.isfinite(v) else f"{100 * v:.1f}%"


def _report_text(event_key, key, event, cfg, provenance, summary, sens, terc, standin_by_tercile,
                 flagged, plac, real, det, det_by_cat, reach_rows, lab, area, km_ranges, written,
                 out_dir) -> str:
    sp, sr = summary["placebo"], summary["real_event"]
    wc = provenance.get("worldcover_exported", {})
    fc = event.flood_check
    a = plac[plac.analysis_set == 1]
    flagged_a = a[a.worldcover_standin == 1]
    lines = [
        f"# Placebo report: {event_key}",
        "",
        "Reproduce with `python -m floodforest.cli placebo --event " + event_key + "`.",
        "",
        "## What was run",
        "",
        f"The real event's pipeline (fetch, label, transects, profiles, table) was run unchanged as "
        f"`{key}`, a no-flood control (spec 8.9, 14.2). Same centerline, reach, manual exclusions, "
        "thresholds and label rules as the real event; only the image windows differ.",
        "",
        "| Input | Placebo window |",
        "|---|---|",
        f"| Sentinel-2 pre (forest, channel, controls) | {event.s2_pre.start} to {event.s2_pre.end} |",
        f"| Sentinel-2 post | {event.s2_post.start} to {event.s2_post.end} |",
        f"| Sentinel-1 pre / post | {event.s1_pre.start[:10]} to {event.s1_pre.end[:10]} / "
        f"{event.s1_post.start[:10]} to {event.s1_post.end[:10]} (end exclusive) |",
        f"| sar_near (sensitivity) | {event.sar_near.start[:10]} to {event.sar_near.end[:10]} |",
        f"| Fake event date | {event.event_dates[0]} (pre cutoff {event.pre_end_utc}, post cutoff "
        f"{event.post_start_utc}) |",
        "",
        f"The user confirmed no notable flood in {fc.interval.start} to {fc.interval.end} "
        f"(`flood_check`, {fc.placebo_pair_review.date}; record `{fc.record}`). That confirmation "
        "covers notable floods only; local floods, high flows and bank erosion are not ruled out.",
        "",
        "The real event's files were only read. The real event's Oct-Dec 2020 forest baseline is "
        "unchanged; it is used here only as a labelled sensitivity.",
        "",
        "## Expected outcome and result",
        "",
        "Expected with no flood (8.9): near-zero reach. The spec sets no numeric pass mark for "
        "'near-zero', so the placebo is shown beside the real event. The second part of 8.9 "
        "(beta1 not distinguishable from zero) needs the `stats` step, which is not built yet.",
        "",
        "| Analysis set, combined label, FABDEM | Placebo 2018-19 | Real event 2021 |",
        "|---|---:|---:|",
        f"| Rows (stations) | {int(sp.analysis_rows)} ({int(sp.analysis_stations)}) | "
        f"{int(sr.analysis_rows)} ({int(sr.analysis_stations)}) |",
        f"| Reach > 0 | {_pct(sp.share_reach_above_0)} | {_pct(sr.share_reach_above_0)} |",
        f"| Overtopped (reach >= {cfg.overtop_m:g} m) | {_pct(sp.share_overtopped)} | "
        f"{_pct(sr.share_overtopped)} |",
        f"| Mean / median / 90th pct reach distance (m) | {_fmt(sp.mean_reach_distance_m, 1)} / "
        f"{_fmt(sp.median_reach_distance_m, 1)} / {_fmt(sp.p90_reach_distance_m, 1)} | "
        f"{_fmt(sr.mean_reach_distance_m, 1)} / {_fmt(sr.median_reach_distance_m, 1)} / "
        f"{_fmt(sr.p90_reach_distance_m, 1)} |",
        f"| Mean reach height (m) | {_fmt(sp.mean_reach_height_m, 1)} | "
        f"{_fmt(sr.mean_reach_height_m, 1)} |",
        f"| Censored | {_pct(sp.share_censored)} | {_pct(sr.share_censored)} |",
        f"| Impact in the combined label, whole corridor (ha) | {_fmt(area['placebo_impact_ha'], 1)} | "
        f"{_fmt(area['real_event_impact_ha'], 1)} |",
        f"| Rows left out: near-bank nodata / river not visible / low valid fraction | "
        f"{int(sp.rows_near_bank_nodata)} / {int(sp.rows_river_not_visible)} / "
        f"{int(sp.rows_low_valid_fraction)} | {int(sr.rows_near_bank_nodata)} / "
        f"{int(sr.rows_river_not_visible)} / {int(sr.rows_low_valid_fraction)} |",
        "",
        "By label variant (analysis set of each run):",
        "",
        "| Label | Placebo overtopped | Placebo mean reach (m) | Real overtopped | Real mean reach (m) |",
        "|---|---:|---:|---:|---:|",
        *[f"| {r.label} | {_pct(r.share_overtopped_placebo)} | {_fmt(r.mean_reach_distance_m_placebo, 1)} "
          f"| {_pct(r.share_overtopped_real)} | {_fmt(r.mean_reach_distance_m_real, 1)} |"
          for r in lab.itertuples()],
        "",
        "Figures: `fig_along_river.png` (share overtopped per river km, both banks), "
        "`fig_reach_distribution.png`, `fig_detections_by_km.png`, `sheet_placebo_detections_*.png`.",
        "",
        "## What the placebo label detected (false positives)",
        "",
        "Every impact patch in the placebo combined label is a false positive in the sense that no "
        "notable flood occurred. Patches are sorted into diagnostic categories by the screening "
        "rules in `placebo_diagnostics` (config.yaml). The categories do not decide what is impact "
        "and do not say what caused a patch; agricultural use and construction cannot be separated "
        "at 10 m. Where a patch lies within "
        f"{cfg.placebo_diagnostics.review_match_m:g} m of a visually reviewed spot "
        f"(`{cfg.placebo_diagnostics.visual_review_file}`), the review is shown.",
        "",
        "| Category | Patches | Area (ha) |",
        "|---|---:|---:|",
        *[f"| {CATEGORIES[c]} | {int(r.patches)} | {r.area_ha:.2f} |"
          for c, r in det_by_cat.iterrows()],
        "",
    ]
    reviewed = det[det.visual_class != ""]
    if len(reviewed):
        lines += [
            "Detections that match a visually reviewed spot:",
            "",
            "| km | Bank | Area (m2) | Category | Visual review | Reviewer |",
            "|---:|---|---:|---|---|---|",
            *[f"| {r.river_km} | {r.side} | {int(r.area_m2)} | {r.category} | {r.visual_class} "
              f"({r.visual_note}) | {r.visual_reviewer} |" for r in reviewed.itertuples()],
            "",
        ]
    big = det.sort_values("area_m2", ascending=False).head(15)
    lines += [
        "Largest detections (all in `placebo_detections.csv` with coordinates):",
        "",
        "| km | Bank | Area (m2) | Distance to channel (m) | Height above channel (m) | Optical / radar share | Category |",
        "|---:|---|---:|---:|---:|---|---|",
        *[f"| {r.river_km} | {r.side} | {int(r.area_m2)} | {r.distance_to_channel_m} | "
          f"{r.height_above_channel_m} | {r.optical_share:.2f} / {r.radar_share:.2f} | {r.category} |"
          for r in big.itertuples()],
        "",
        f"Analysis rows with a reach above 0 m: {len(reach_rows)}. The patch at their reach edge, "
        "by category:",
        "",
        "| Category at the reach edge | Rows |",
        "|---|---:|",
        *[f"| {c or 'none'} | {n} |" for c, n in reach_rows.edge_category.value_counts().items()],
        "",
        "Affected river reaches (river km with any analysis row overtopped in the placebo):",
        "",
    ]
    for side in ("left", "right"):
        km = np.sort(a.loc[(a.side == side) & (a[f"overtopped_{PRIMARY_DEM}"] == 1), "river_km"]
                     .to_numpy())
        lines.append(f"- {side}: {km_ranges(km, cfg.station_spacing_m)}")
    lines += ["", "## Rows lost because the river is not visible (14.20)", ""]
    for name, t in (("Placebo (2018 forest)", plac), ("Real event (2020 forest)", real)):
        nv =t[t["exclude_reasons"].fillna("").str.contains("river_not_visible")]
        kept = t[t["analysis_set"] == 1]
        km = np.sort(nv.drop_duplicates("station_id")["river_km"].to_numpy())
        lines.append(
            f"- {name}: {nv.station_id.nunique()} stations ({len(nv)} rows); mean "
            f"forest fraction of those rows {_fmt(nv[f'forest_frac_{cfg.belt_width_m}'].mean())} "
            f"vs {_fmt(kept[f'forest_frac_{cfg.belt_width_m}'].mean())} in the analysis set. "
            f"River km: {km_ranges(km, cfg.station_spacing_m)}.")
    lines += [
        "",
        "Rule 14.20 keeps a row only if bare or impacted ground within 30 m of the station shows "
        "where the river is. After the 2021 flood the bed was bare almost everywhere; in an "
        "ordinary year a narrow river under canopy is often invisible at 10 m. The placebo "
        "analysis set therefore has fewer rows, mostly fewer forested gorge rows, than the real "
        "event. The 14.20 sensitivity (those rows with reach 0) is in the table at the end of "
        "the forest section.",
        "",
        "## Forest exposure",
        "",
        "Primary (14.2): forest from the placebo pre window (Oct-Dec 2018), with Dynamic World, "
        "ETH canopy height 2020 and the Hansen rule (skipped for 2018, 14.3). "
        f"{wc.get('reason', '')}",
        "",
        f"- Rows using the WorldCover stand-in (`worldcover_standin` = 1): {int(len(flagged))} of "
        f"{len(plac)} rows; {len(flagged_a)} of {len(a)} analysis rows "
        f"({_pct(len(flagged_a) / len(a)) if len(a) else 'NA'})."
        + (" The stand-in was not needed: every decided belt pixel had a Dynamic World value."
           if len(flagged) == 0 else ""),
        "- Forest source per row, 100 m belt: "
        + ", ".join(f"{k} {v}" for k, v in
                    plac[f"forest_source_{cfg.belt_width_m}"].value_counts().items()) + ".",
    ]
    for side in ("left", "right"):
        km = np.sort(flagged.loc[flagged.side == side, "river_km"].to_numpy())
        lines.append(f"- {side} bank, river km: {km_ranges(km, cfg.station_spacing_m)}")
    lines += [
        "- Stand-in analysis rows per 2018 forest tercile: "
        + ", ".join(f"{t} {n}" for t, n in standin_by_tercile.items()),
        "",
        "Descriptive only; no model is fitted (the `stats` step is not built). Terciles are "
        "equal-count groups of the analysis set; ties (many rows have no forest) are split by "
        "station order. With so few overtopped rows, differences between terciles are a few "
        "rows and should not be read as a trend.",
        "",
    ]
    for name, t in terc.items():
        lines += [f"**{name}**", "", "| Tercile | Rows | Forest fraction range | Overtopped rows | "
                  "Overtopped | Mean reach (m) |", "|---|---:|---|---:|---:|---:|"]
        lines += [f"| {r.tercile} | {_fmt(r.rows, 0)} | {_fmt(r.forest_frac_min)}-{_fmt(r.forest_frac_max)} "
                  f"| {_fmt(r.overtopped_rows, 0)} | {_pct(r.share_overtopped)} | "
                  f"{_fmt(r.mean_reach_distance_m, 1)} |" for r in t.itertuples()]
        lines.append("")
    lines += [
        "Stand-in sensitivity, beside the main placebo result:",
        "",
        "| Run | Analysis rows | Reach > 0 | Overtopped | Mean reach (m) |",
        "|---|---:|---:|---:|---:|",
        *[f"| {name} | {int(r.analysis_rows)} | {_pct(r.share_reach_above_0)} | "
          f"{_pct(r.share_overtopped)} | {_fmt(r.mean_reach_distance_m, 1)} |"
          for name, r in sens.iterrows()],
        "",
        "## Limitations",
        "",
        "- The placebo shows how much impact the method finds when no notable flood occurred. It does "
        "not prove that no local flood, high flow or bank erosion happened in the period.",
        "- The forest-effect half of the placebo test (beta1, 8.9) waits for the `stats` step.",
        "- Diagnostic categories use 10 m imagery and a 30 m DEM; they are review aids, not causes. "
        "Visual reviews exist only for the spots in the review file.",
        "- The placebo forest uses ETH canopy height 2020 (after the placebo period) and, where "
        "Dynamic World has no 2018 image, the WorldCover 2020 stand-in (flagged rows).",
        "- The 2021 validation edges (14.25) belong to the real event and are not applicable here.",
        "- The centerline was traced on post-2021 imagery; the placebo channel and bank tops come "
        "from the 2018 composite, so distances to the line are approximate in the lower valley.",
        "",
        "## Files",
        "",
        *[f"- `{p.relative_to(out_dir)}`" for p in written if p.parent == out_dir],
        "",
    ]
    return "\n".join(lines) + "\n"
