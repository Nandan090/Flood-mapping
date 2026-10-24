"""Check a verified centerline against the pre-event active channel.

The line file is only read. Nothing here writes a new centerline.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.features import rasterize
from rasterio.transform import rowcol
from scipy import ndimage

from . import sheets
from .centerline import assert_centerline_verified
from .config import PROJECT_ROOT, load_config
from .geometry import densify
from .optical import mndwi, ndvi

# Display choice for this check, not a study threshold.
_OFF_CHANNEL_M = 50.0


def qa_centerline(event_key: str, config_path: str | Path | None = None) -> list[Path]:
    """Sample the verified line every 10 m and write the contact sheets."""
    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    line, samples, river_m = line_samples(event, cfg)
    raster_path = (PROJECT_ROOT / cfg.paths.data / event_key / "rasters"
                   / "s2_pre.tif")
    if not raster_path.exists():
        raise FileNotFoundError(
            f"Pre-event Sentinel-2 composite not found: {raster_path}"
        )
    grid = _read_pre_event(raster_path, cfg.crs)
    channel = _active_channel(
        grid, line, cfg.mndwi_channel_min, cfg.ndvi_channel_max)
    rows, cols = _rows_cols(grid["transform"], samples, grid["height"],
                            grid["width"], river_m)
    inside = channel[rows, cols]
    distance_m = ndimage.distance_transform_edt(
        ~channel, sampling=grid["pixel_m"])[rows, cols]

    _print_km_table(river_m, inside, distance_m)
    _print_stretches(samples, river_m, inside, cfg.crs)

    out_dir = PROJECT_ROOT / cfg.paths.outputs / event_key / "qa"
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = _write_sheets(out_dir, grid, samples, river_m, inside)
    for path in paths:
        print(path.relative_to(PROJECT_ROOT))
    return paths


def off_channel_stretches(river_m: np.ndarray, inside: np.ndarray,
                          longer_than_m: float) -> list[tuple[int, int, float]]:
    """Runs of off-channel samples whose along-line length exceeds longer_than_m.

    Length is the distance from the first sample of the run to the last.
    """
    river_m = np.asarray(river_m, dtype=float)
    inside = np.asarray(inside, dtype=bool)
    stretches = []
    i = 0
    n = len(river_m)
    while i < n:
        if inside[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and not inside[j + 1]:
            j += 1
        length_m = float(river_m[j] - river_m[i])
        if length_m > longer_than_m:
            stretches.append((i, j, length_m))
        i = j + 1
    return stretches


def channel_touching_line(index: np.ndarray, on_line: np.ndarray) -> np.ndarray:
    """Index pixels kept only when 8-connected to the centerline (section 5.2)."""
    index = np.asarray(index, dtype=bool)
    on_line = np.asarray(on_line, dtype=bool)
    labels, _ = ndimage.label(index, structure=np.ones((3, 3), dtype=int))
    seeds = np.unique(labels[index & on_line])
    seeds = seeds[seeds != 0]
    if seeds.size == 0:
        return np.zeros(index.shape, dtype=bool)
    return np.isin(labels, seeds)


def line_samples(event, cfg):
    """The verified line in the analysis CRS, its 10 m samples and their river_m."""
    path = PROJECT_ROOT / event.centerline_file
    assert_centerline_verified(path)
    gdf = gpd.read_file(path)
    if "role" not in gdf.columns:
        raise RuntimeError(f"{path} has no 'role' property.")
    lines = gdf[gdf["role"] == "centerline"]
    if len(lines) != 1:
        raise RuntimeError(
            f"{path} has {len(lines)} centerline features; expected one."
        )
    projected = lines.to_crs(cfg.crs)
    if projected.crs is None or projected.crs.to_string() != cfg.crs:
        raise RuntimeError(f"Centerline did not reproject to {cfg.crs}.")
    line = projected.geometry.iloc[0]
    samples = densify(np.asarray(line.coords), cfg.grid_m)
    river_m = np.zeros(len(samples))
    river_m[1:] = np.cumsum(np.linalg.norm(np.diff(samples, axis=0), axis=1))
    return line, samples, river_m


def _read_pre_event(path: Path, crs: str) -> dict:
    needed = ("B2", "B3", "B4", "B8", "B11")
    with rasterio.open(path) as src:
        got = src.crs.to_string() if src.crs is not None else None
        if got != crs:
            raise RuntimeError(f"{path.name} is {got}, not {crs}.")
        names = list(src.descriptions)
        missing = [band for band in needed if band not in names]
        if missing:
            raise RuntimeError(
                f"{path.name} is missing {missing}. Bands present: {names}."
            )
        bands = {name: _read_band(src, names.index(name) + 1) for name in needed}
        pixel_m = abs(float(src.transform.a))
        return {
            "bands": bands,
            "transform": src.transform,
            "height": src.height,
            "width": src.width,
            "pixel_m": pixel_m,
            "crs": crs,
        }


def _read_band(src, index: int) -> np.ndarray:
    data = src.read(index, masked=True).astype(float)
    return np.ma.filled(data, np.nan)


def _active_channel(grid: dict, line, mndwi_min: float, ndvi_max: float) -> np.ndarray:
    bands = grid["bands"]
    green_red = ndvi(bands["B8"], bands["B4"])
    water = mndwi(bands["B3"], bands["B11"])
    index = ((water > mndwi_min) | (green_red < ndvi_max))
    index &= np.isfinite(water) & np.isfinite(green_red)
    on_line = rasterize(
        [(line, 1)],
        out_shape=(grid["height"], grid["width"]),
        transform=grid["transform"],
        fill=0,
        dtype="uint8",
        all_touched=True,
    ).astype(bool)
    channel = channel_touching_line(index, on_line)
    if not channel.any():
        raise RuntimeError(
            "No pre-event active-channel pixel is 8-connected to the "
            "centerline. The index rule was not replaced with another mask."
        )
    return channel


def _rows_cols(transform, samples, height, width, river_m):
    rows, cols = rowcol(transform, samples[:, 0], samples[:, 1])
    rows = np.asarray(rows)
    cols = np.asarray(cols)
    outside = (rows < 0) | (cols < 0) | (rows >= height) | (cols >= width)
    if outside.any():
        km = river_m[int(np.flatnonzero(outside)[0])] / 1000
        raise RuntimeError(
            f"The centerline leaves the pre-event image at river km {km:.3f}."
        )
    return rows, cols


def _print_km_table(river_m, inside, distance_m) -> None:
    km = np.floor(river_m / 1000).astype(int)
    print("river_km  points  share_inside  max_distance_m")
    for bin_id in range(int(km.min()), int(km.max()) + 1):
        keep = km == bin_id
        share = float(inside[keep].mean())
        furthest = float(distance_m[keep].max())
        print(f"{bin_id:8d}  {int(keep.sum()):6d}  {share:12.2f}  {furthest:14.1f}")


def _print_stretches(samples, river_m, inside, crs: str) -> None:
    stretches = off_channel_stretches(river_m, inside, _OFF_CHANNEL_M)
    to_ll = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    print(f"Off-channel stretches longer than {_OFF_CHANNEL_M:.0f} m: {len(stretches)}")
    for start, end, length_m in stretches:
        lon0, lat0 = to_ll.transform(samples[start, 0], samples[start, 1])
        lon1, lat1 = to_ll.transform(samples[end, 0], samples[end, 1])
        print(
            f"  {river_m[start] / 1000:.3f}–{river_m[end] / 1000:.3f} km "
            f"({length_m:.0f} m)"
        )
        print(f"    start {lat0:.6f} N, {lon0:.6f} E")
        print(f"    end   {lat1:.6f} N, {lon1:.6f} E")


def _write_sheets(out_dir, grid, samples, river_m, inside) -> list[Path]:
    from matplotlib.lines import Line2D

    legend = [
        Line2D([0], [0], color="yellow", lw=2, label="centerline"),
        Line2D([0], [0], color="red", lw=2, label="outside the channel"),
    ]
    written = []
    for mode, names in (("truecolour", sheets.TRUE_COLOUR),
                        ("falsecolour", sheets.FALSE_COLOUR)):
        limits = sheets.fixed_limits([grid["bands"]], names)

        def draw(ax, win, names=names, limits=limits):
            sheets.show(ax, sheets.rgb(grid["bands"], names, limits, win), win)
            sheets.finish_panel(ax, win)
            _plot_off_channel(ax, samples, river_m, inside,
                              win["start_m"], win["end_m"])

        written += sheets.write_sheets(
            out_dir, mode,
            f"{mode}. Same 2nd–98th percentile stretch on every panel",
            grid, samples, river_m, draw, legend)
    return written


def _plot_off_channel(ax, samples, river_m, inside, start_m, end_m) -> None:
    show = (river_m >= start_m) & (river_m <= end_m)
    off = show & ~inside
    if not off.any():
        return
    indices = np.flatnonzero(off)
    breaks = np.where(np.diff(indices) > 1)[0]
    for piece in np.split(indices, breaks + 1):
        # Include the neighbouring sample so a one-point gap still draws.
        lo = max(int(piece[0]) - 1, 0)
        hi = min(int(piece[-1]) + 1, len(samples) - 1)
        ax.plot(samples[lo:hi + 1, 0], samples[lo:hi + 1, 1],
                color="red", lw=2.2, solid_capstyle="round")


