# Placebo-period visual check (lende2025)

These two scripts are **visual-check tools, not pipeline steps**. They are not part of
`python -m floodforest.cli`, they do not apply any label rule, and they change nothing in
`config.yaml`, `data/` or `outputs/`. They were used to review the lende2025 placebo period
(spec 14.2) before its `flood_check` was confirmed. The images they made are kept in
`results/placebo_check_lende/`.

## What they produce

`lende_placebo_download.py` writes to `results/placebo_check_lende/<subfolder>/`:

- `rasters/`: Sentinel-2 composites (`s2_2022.tif`, `s2_2023.tif`), their clear-observation
  counts, and Sentinel-1 per-orbit median VV and VH (`s1_2022.tif`, `s1_2023.tif`). These are
  git-ignored.
- `scenes_s2.csv`: every Sentinel-2 scene in each window, with the share of the region that is
  clear.
- `scenes_s1.csv`: every Sentinel-1 scene in each window, with orbit, pass, and the share of the
  region with valid pixels.
- `manifest.json`: the grid and the scene IDs used in each composite.

`lende_placebo_chips.py` reads those rasters and writes:

- `chips_optical/optical_kmAA-BB.png`, one per river km: 2022 and 2023 true colour, NDVI change
  and bare-soil index change (2023 minus 2022).
- `chips_radar/radar_kmAA-BB.png`, one per river km: for each orbit, VV 2022, VV 2023, VV change
  and VH change.
- `overview_optical_kmAA-BB.png` and `overview_radar_kmAA-BB.png`: the whole km range.

Every image shows the verified river line, km labels and a dashed band 300 m each side of the
line, north up. Each chip covers its km of river plus 300 m each side.

## Arguments

Both scripts take the same optional arguments: `km_from km_to subfolder`. Without arguments they
cover km 0-15 and write to the base folder. Run them from the repository root:

```
PYTHONPATH=src python scripts/placebo_check/lende_placebo_download.py 15 35 km15-35
PYTHONPATH=src python scripts/placebo_check/lende_placebo_chips.py 15 35 km15-35
```

The download script fetches the line from km_from - 0.5 to km_to + 0.5, with 1 km each side, on
the study's 10 m grid (EPSG:32645). Run the download script first; the chip script reads its
output.

## Windows and composite rules

The windows are read from `config.yaml` (`lende2025_placebo`, derived by `floodforest.config`).
None are written into the scripts.

| Input | Window (UTC, end exclusive) |
|---|---|
| Sentinel-2 placebo "pre" | 2022-10-01 to 2022-12-16 (1 Oct to 15 Dec 2022) |
| Sentinel-2 placebo "post" | 2023-10-01 to 2023-12-16 (1 Oct to 15 Dec 2023) |
| Sentinel-1 placebo "pre" | 2022-07-09 to 2022-08-15 |
| Sentinel-1 placebo "post" | 2023-07-09 to 2023-08-15 |

- **Sentinel-2.** The study's own composite (`ee_data._s2_composite`): Cloud Score+
  `cs_cdf >= cloud_score_min`, per-band median, and pixels with fewer than `clear_obs_min` clear
  observations set to nodata.
- **Sentinel-1.** IW, VV and VH (`ee_data._s1_collection`). One median per relative orbit and
  pass, never across orbits. An orbit is used only if its scenes have valid pixels over more
  than half of the region in both windows. For display, the chip script applies the study's
  3x3 focal median (`sar_focal_median_px`). No layover or shadow mask is applied.
- **Change maps.** For display only. NDVI and BSI follow spec 5.2. The colour limits are fixed
  in the chip script (NDVI +/-0.3, BSI +/-0.15, radar +/-6 dB).

## Requirements

- **Earth Engine:** the download script needs it. Run `earthengine authenticate` once and set
  `ee_project` in `config.yaml`.
- **Verified centerline:** both scripts use the verified lende2025 centerline and refuse a line
  that is not verified.
- **Python environment:** the project's own (`requirements.txt`).
