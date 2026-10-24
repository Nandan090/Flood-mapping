"""Contact sheets: one north-up panel per 2 km of river, six panels per PNG.

Used by qa_centerline and label so every check is drawn on the same panels.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
import rasterio
from rasterio.windows import Window, from_bounds
from shapely.geometry import LineString

# Display choices, not study thresholds.
PANEL_KM = 2.0
PANELS_PER_SHEET = 6
CONTEXT_M = 250.0
STRETCH_PERCENTILES = (2.0, 98.0)

TRUE_COLOUR = ("B4", "B3", "B2")
FALSE_COLOUR = ("B8", "B4", "B3")


def panels(length_m: float) -> list[tuple[float, float]]:
    out = []
    start = 0.0
    while start < length_m - 1e-6:
        out.append((start, min(start + PANEL_KM * 1000, length_m)))
        start += PANEL_KM * 1000
    return out


def fixed_limits(images: list[dict], names: tuple[str, ...]) -> list[tuple[float, float]]:
    """2nd–98th percentile per displayed band, pooled over the given images."""
    low, high = STRETCH_PERCENTILES
    limits = []
    for name in names:
        values = np.concatenate([
            image[name][np.isfinite(image[name])] for image in images])
        if values.size == 0:
            raise RuntimeError(f"Band {name} has no valid pixels.")
        lo, hi = np.percentile(values, [low, high])
        if hi <= lo:
            hi = lo + 1e-6
        limits.append((float(lo), float(hi)))
    return limits


def panel_window(grid: dict, samples: np.ndarray, river_m: np.ndarray,
                 panel: tuple[float, float]) -> dict:
    start_m, end_m = panel
    keep = (river_m >= start_m) & (river_m <= end_m)
    part = samples[keep]
    segment = LineString(part)
    bounds = segment.buffer(CONTEXT_M).bounds
    window = from_bounds(*bounds, transform=grid["transform"])
    window = window.round_offsets().round_lengths()
    window = window.intersection(Window(0, 0, grid["width"], grid["height"]))
    row0, col0 = int(window.row_off), int(window.col_off)
    transform = rasterio.windows.transform(window, grid["transform"])
    left, top = transform.c, transform.f
    return {
        "rows": slice(row0, row0 + int(window.height)),
        "cols": slice(col0, col0 + int(window.width)),
        "extent": (left, left + transform.a * window.width,
                   top + transform.e * window.height, top),
        "bounds": bounds,
        "part": part,
        "segment": segment,
        "start_m": start_m,
        "end_m": end_m,
    }


def rgb(bands: dict, names: tuple[str, ...], limits, win: dict) -> np.ndarray:
    stack = np.stack([bands[name][win["rows"], win["cols"]] for name in names])
    out = np.zeros((stack.shape[1], stack.shape[2], 3), dtype=float)
    for i, (lo, hi) in enumerate(limits):
        out[:, :, i] = np.clip((stack[i] - lo) / (hi - lo), 0, 1)
    out[~np.isfinite(stack).all(axis=0)] = 0
    return out


def show(ax, image: np.ndarray, win: dict) -> None:
    ax.imshow(image, extent=win["extent"], origin="upper", interpolation="nearest")


def finish_panel(ax, win: dict, line_colour: str = "yellow",
                 line_width: float = 1.6) -> None:
    from matplotlib import patheffects
    part = win["part"]
    if line_colour:
        ax.plot(part[:, 0], part[:, 1], color=line_colour, lw=line_width,
                solid_capstyle="round")
    start_m, end_m = win["start_m"], win["end_m"]
    for km in range(int(np.ceil(start_m / 1000)), int(np.floor(end_m / 1000)) + 1):
        point = win["segment"].interpolate(km * 1000 - start_m)
        ax.annotate(
            str(km), (point.x, point.y), color="white", fontsize=8,
            ha="left", va="bottom",
            path_effects=[patheffects.withStroke(linewidth=2, foreground="black")],
        )
    minx, miny, maxx, maxy = win["bounds"]
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)
    ax.set_aspect("equal")
    ax.set_title(f"km {start_m / 1000:.1f}–{end_m / 1000:.1f}", fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])


def write_sheets(out_dir: Path, name: str, title: str, grid: dict,
                 samples: np.ndarray, river_m: np.ndarray,
                 draw: Callable, legend: list) -> list[Path]:
    """draw(ax, win) paints one panel; the km labels and frame are added here."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    all_panels = panels(float(river_m[-1]))
    written = []
    for sheet_i in range(0, len(all_panels), PANELS_PER_SHEET):
        chunk = all_panels[sheet_i:sheet_i + PANELS_PER_SHEET]
        fig, axes = plt.subplots(2, 3, figsize=(12, 10), dpi=120)
        for ax, panel in zip(axes.ravel(), chunk):
            draw(ax, panel_window(grid, samples, river_m, panel))
        for ax in axes.ravel()[len(chunk):]:
            ax.axis("off")
        first = chunk[0][0] / 1000
        last = chunk[-1][1] / 1000
        fig.suptitle(f"{title}, km {first:.1f}–{last:.1f}. North is up.",
                     fontsize=10)
        if legend:
            fig.legend(handles=legend, loc="lower center", ncol=len(legend),
                       frameon=False, fontsize=8)
        fig.tight_layout(rect=(0, 0.04, 1, 0.96))
        path = out_dir / f"{name}_km{first:04.1f}-{last:04.1f}.png"
        fig.savefig(path)
        plt.close(fig)
        written.append(path)
    return written
