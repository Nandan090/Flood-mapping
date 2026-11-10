"""Review chips and overviews for the lende2025 placebo-period check, km 0-15.

Reads results/placebo_check_lende/rasters/ (from lende_placebo_download.py).
Optical chip per km: Oct-Dec 2022 | Oct-Dec 2023 true colour | dNDVI | dBSI (2023 - 2022).
Radar chip per km: per orbit VV 2022 | VV 2023 | dVV | dVH (dB, 3x3 focal median as in 5.3).
North up, each chip = that km of river plus 300 m on each side. Display only;
no label rule is applied and nothing is written outside results/placebo_check_lende/.
"""
import json

import matplotlib
import numpy as np
import rasterio
from matplotlib import patheffects
from matplotlib.lines import Line2D
from shapely.ops import substring

from floodforest import optical, sar, sheets
from floodforest.config import PROJECT_ROOT, load_config
from floodforest.label import _read_named
from floodforest.qa import line_samples

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import sys

# Usage: lende_placebo_chips.py [km_from km_to subfolder]; default km 0-15 in the base folder.
KM_FROM, KM_TO = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, 15)
SUB = sys.argv[3] if len(sys.argv) > 3 else ""
SIDE_M = 300
D_NDVI, D_BSI, D_DB = 0.3, 0.15, 6.0      # colour-bar limits (display only)

cfg = load_config(PROJECT_ROOT / "config.yaml")
base = PROJECT_ROOT / "results" / "placebo_check_lende" / SUB
rast = base / "rasters"
(base / "chips_optical").mkdir(exist_ok=True)
(base / "chips_radar").mkdir(exist_ok=True)
manifest = json.loads((base / "manifest.json").read_text(encoding="utf-8"))

B = ("B2", "B3", "B4", "B8", "B11")
P, grid = _read_named(rast / "s2_2022.tif", B, cfg.crs)
Q = _read_named(rast / "s2_2023.tif", B, cfg.crs, grid)[0]
cnt = {n: _read_named(rast / f"{n}_clear_count.tif", ("clear_count",), cfg.crs, grid)[0]["clear_count"]
       for n in ("s2_2022", "s2_2023")}
d_ndvi = optical.ndvi(Q["B8"], Q["B4"]) - optical.ndvi(P["B8"], P["B4"])
d_bsi = (optical.bsi(Q["B2"], Q["B4"], Q["B8"], Q["B11"])
         - optical.bsi(P["B2"], P["B4"], P["B8"], P["B11"]))

with rasterio.open(rast / "s1_2022.tif") as src:
    s1_names = list(src.descriptions)
orbits = sorted({n.split("_")[1] for n in s1_names})          # e.g. 19D, 85A
S1 = {}
for year in ("2022", "2023"):
    bands = _read_named(rast / f"s1_{year}.tif", s1_names, cfg.crs, grid)[0]
    S1[year] = {n: sar.focal_median(v, cfg.sar_focal_median_px) for n, v in bands.items()}
n_scenes = {y: {k: len(v) for k, v in manifest["s1"][f"s1_{y}"]["scenes_by_orbit"].items()}
            for y in ("2022", "2023")}

line, _, _ = line_samples(cfg.events.lende2025, cfg)
t = grid["transform"]
x0, y0 = t.c, t.f
extent_all = [x0, x0 + grid["width"] * 10, y0 - grid["height"] * 10, y0]
TC = sheets.TRUE_COLOUR
lim = sheets.fixed_limits([P, Q], TC)
vv_lim = {o: tuple(np.nanpercentile(np.concatenate([S1[y][f"VV_{o}"].ravel() for y in S1]), (2, 98)))
          for o in orbits}
FX = [patheffects.withStroke(linewidth=2, foreground="black")]


def window(bounds):
    minx, miny, maxx, maxy = bounds
    c0, c1 = int((minx - x0) // 10), int(np.ceil((maxx - x0) / 10))
    r0, r1 = int((y0 - maxy) // 10), int(np.ceil((y0 - miny) / 10))
    c0, r0 = max(c0, 0), max(r0, 0)
    c1, r1 = min(c1, grid["width"]), min(r1, grid["height"])
    return (slice(r0, r1), slice(c0, c1),
            [x0 + c0 * 10, x0 + c1 * 10, y0 - r1 * 10, y0 - r0 * 10])


def tc(img, rs, cs):
    return sheets.rgb(img, TC, lim, {"rows": rs, "cols": cs})


def decorate(ax, ext, km_lo, km_hi, buffer_line=None, every=1, fs=8):
    xs, ys = line.xy
    ax.plot(xs, ys, color="yellow", lw=0.9)
    if buffer_line is not None:
        bx, by = buffer_line.exterior.xy
        ax.plot(bx, by, color="white", lw=0.6, ls="--")
    for km in range(km_lo, km_hi + 1, every):
        p = line.interpolate(km * 1000)
        ax.plot(p.x, p.y, "o", ms=3, color="yellow", mec="black", mew=0.5)
        ax.annotate(f"km {km}", (p.x, p.y), xytext=(4, 2), textcoords="offset points",
                    color="white", fontsize=fs, path_effects=FX)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def show(ax, img, ext, **kw):
    return ax.imshow(img, extent=ext, origin="upper", interpolation="nearest", **kw)


def bar(fig, ax, im, label):
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(label, fontsize=7)
    cb.ax.tick_params(labelsize=6)


written = []
for k in range(KM_FROM, KM_TO):
    seg = substring(line, k * 1000, (k + 1) * 1000)
    box = seg.buffer(SIDE_M)
    rs, cs, ext = window(box.bounds)
    cc = {n: (int(np.nanmin(c[rs, cs])), int(np.nanmax(c[rs, cs]))) for n, c in cnt.items()}

    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), dpi=120)
    show(axes[0], tc(P, rs, cs), ext)
    axes[0].set_title(f"Sentinel-2 Oct 1-Dec 15 2022\nclear obs {cc['s2_2022'][0]}-{cc['s2_2022'][1]}", fontsize=8)
    show(axes[1], tc(Q, rs, cs), ext)
    axes[1].set_title(f"Sentinel-2 Oct 1-Dec 15 2023\nclear obs {cc['s2_2023'][0]}-{cc['s2_2023'][1]}", fontsize=8)
    im = show(axes[2], d_ndvi[rs, cs], ext, cmap="BrBG", vmin=-D_NDVI, vmax=D_NDVI)
    axes[2].set_title("NDVI change 2023 - 2022\n(brown = less green)", fontsize=8)
    bar(fig, axes[2], im, "dNDVI")
    im = show(axes[3], d_bsi[rs, cs], ext, cmap="PuOr_r", vmin=-D_BSI, vmax=D_BSI)
    axes[3].set_title("Bare-soil index change 2023 - 2022\n(orange = more bare)", fontsize=8)
    bar(fig, axes[3], im, "dBSI")
    for ax in axes:
        decorate(ax, ext, k, k + 1, buffer_line=box)
    fig.suptitle(f"lende2025 placebo check, river km {k}-{k + 1}. Yellow = verified river line, "
                 f"dashed = {SIDE_M} m each side. North up.", fontsize=9)
    fig.tight_layout()
    p = base / "chips_optical" / f"optical_km{k:02d}-{k + 1:02d}.png"
    fig.savefig(p)
    plt.close(fig)
    written.append(p)

    fig, axes = plt.subplots(len(orbits), 4, figsize=(15, 4.4 * len(orbits)), dpi=120, squeeze=False)
    for row, o in zip(axes, orbits):
        key = next(kk for kk in n_scenes["2022"] if kk.startswith(o[:-1]))
        lo, hi = vv_lim[o]
        for ax, y in zip(row[:2], ("2022", "2023")):
            show(ax, S1[y][f"VV_{o}"][rs, cs], ext, cmap="gray", vmin=lo, vmax=hi)
            ax.set_title(f"VV orbit {o}, Jul 9-Aug 15 {y}\nmedian of {n_scenes[y][key]} scenes", fontsize=8)
        im = show(row[2], S1["2023"][f"VV_{o}"][rs, cs] - S1["2022"][f"VV_{o}"][rs, cs], ext,
                  cmap="RdBu_r", vmin=-D_DB, vmax=D_DB)
        row[2].set_title(f"VV change 2023 - 2022, orbit {o}", fontsize=8)
        bar(fig, row[2], im, "dB")
        im = show(row[3], S1["2023"][f"VH_{o}"][rs, cs] - S1["2022"][f"VH_{o}"][rs, cs], ext,
                  cmap="RdBu_r", vmin=-D_DB, vmax=D_DB)
        row[3].set_title(f"VH change 2023 - 2022, orbit {o}", fontsize=8)
        bar(fig, row[3], im, "dB")
        for ax in row:
            decorate(ax, ext, k, k + 1, buffer_line=box)
    fig.suptitle(f"lende2025 placebo check, Sentinel-1, river km {k}-{k + 1}. Per orbit, never "
                 f"across orbits. 3x3 focal median. No layover/shadow mask applied. North up.", fontsize=9)
    fig.tight_layout()
    p = base / "chips_radar" / f"radar_km{k:02d}-{k + 1:02d}.png"
    fig.savefig(p)
    plt.close(fig)
    written.append(p)

# Overviews: the whole region. A tall region gets one row of panels, a
# squarer one a grid, so the river stays legible.
km_all = substring(line, KM_FROM * 1000, KM_TO * 1000).buffer(SIDE_M)
rs, cs, ext = slice(0, grid["height"]), slice(0, grid["width"]), extent_all
square = grid["width"] / grid["height"] > 0.6
tag = f"km{KM_FROM:02d}-{KM_TO:02d}"
where = f"river km {KM_FROM}-{KM_TO}" + (" (Rasuwagadhi at km 0)" if KM_FROM == 0 else "")
if square:
    fig, axes = plt.subplots(2, 2, figsize=(16, 15.5), dpi=130)
    axes = axes.ravel()
else:
    fig, axes = plt.subplots(1, 4, figsize=(17, 10.5), dpi=130)
show(axes[0], tc(P, rs, cs), ext)
axes[0].set_title("Sentinel-2 Oct 1-Dec 15 2022", fontsize=10)
show(axes[1], tc(Q, rs, cs), ext)
axes[1].set_title("Sentinel-2 Oct 1-Dec 15 2023", fontsize=10)
im = show(axes[2], d_ndvi, ext, cmap="BrBG", vmin=-D_NDVI, vmax=D_NDVI)
axes[2].set_title("NDVI change 2023 - 2022", fontsize=10)
bar(fig, axes[2], im, "dNDVI")
im = show(axes[3], d_bsi, ext, cmap="PuOr_r", vmin=-D_BSI, vmax=D_BSI)
axes[3].set_title("Bare-soil index change 2023 - 2022", fontsize=10)
bar(fig, axes[3], im, "dBSI")
for ax in axes:
    decorate(ax, ext, KM_FROM, KM_TO, buffer_line=km_all, every=1, fs=7)
fig.suptitle(f"lende2025 placebo check, optical overview, {where}. "
             "Yellow = verified line, dashed = 300 m each side.", fontsize=11)
fig.tight_layout()
p = base / f"overview_optical_{tag}.png"
fig.savefig(p)
plt.close(fig)
written.append(p)

n_panels = 2 * len(orbits) + 1
if square:
    fig, axes = plt.subplots(2, 3, figsize=(22, 15.5), dpi=130)
    axes = axes.ravel()
    for ax in axes[n_panels:]:
        ax.axis("off")
    axes = axes[:n_panels]
else:
    fig, axes = plt.subplots(1, n_panels, figsize=(4.3 * n_panels, 10.5), dpi=130)
o0 = orbits[0]
show(axes[0], S1["2023"][f"VV_{o0}"], ext, cmap="gray", vmin=vv_lim[o0][0], vmax=vv_lim[o0][1])
axes[0].set_title(f"VV orbit {o0}, Jul 9-Aug 15 2023", fontsize=10)
i = 1
for o in orbits:
    for pol in ("VV", "VH"):
        im = show(axes[i], S1["2023"][f"{pol}_{o}"] - S1["2022"][f"{pol}_{o}"], ext,
                  cmap="RdBu_r", vmin=-D_DB, vmax=D_DB)
        axes[i].set_title(f"{pol} change 2023 - 2022, orbit {o}", fontsize=10)
        bar(fig, axes[i], im, "dB")
        i += 1
for ax in axes:
    decorate(ax, ext, KM_FROM, KM_TO, buffer_line=km_all, every=1, fs=7)
fig.suptitle(f"lende2025 placebo check, Sentinel-1 overview, river km {KM_FROM}-{KM_TO}. Jul 9-Aug 15 "
             "2023 minus Jul 9-Aug 15 2022, per orbit, 3x3 focal median, no layover/shadow mask.",
             fontsize=11)
fig.tight_layout()
p = base / f"overview_radar_{tag}.png"
fig.savefig(p)
plt.close(fig)
written.append(p)

for p in written:
    print(p.relative_to(PROJECT_ROOT))
