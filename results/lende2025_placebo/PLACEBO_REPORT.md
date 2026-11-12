# Placebo report: lende2025

Reproduce with `python -m floodforest.cli placebo --event lende2025`.

## What was run

The real event's pipeline (fetch, label, transects, profiles, table) was run unchanged as `lende2025_placebo`, a no-flood control (spec 8.9, 14.2). Same centerline, reach, manual exclusions, thresholds and label rules as the real event; only the image windows differ.

| Input | Placebo window |
|---|---|
| Sentinel-2 pre (forest, channel, controls) | 2022-10-01 to 2022-12-15 |
| Sentinel-2 post | 2023-10-01 to 2023-12-15 |
| Sentinel-1 pre / post | 2022-07-09 to 2022-08-15 / 2023-07-09 to 2023-08-15 (end exclusive) |
| sar_near (sensitivity) | 2023-05-01 to 2023-07-08 |
| Fake event date | 2023-07-08 (pre cutoff 2023-07-08T00:00:00Z, post cutoff 2023-07-09T00:00:00Z) |

The user confirmed no notable flood in 2022-07-09 to 2023-12-15 (`flood_check`, 2026-11-10; record `results/placebo_check_lende/`). That confirmation covers notable floods only; local floods, high flows and bank erosion are not ruled out.

The real event's files were only read. The real event's Oct-Dec 2024 forest baseline is unchanged; it is used here only as a labelled sensitivity.

## Expected outcome and result

Expected with no flood (8.9): near-zero reach. The spec sets no numeric pass mark for 'near-zero', so the placebo is shown beside the real event. The second part of 8.9 (beta1 not distinguishable from zero) is fitted by the `stats` step through the same Model A and Model B code as the real event: see `outputs/lende2025/stats/STATS_REPORT.md`.

| Analysis set, combined label, FABDEM | Placebo 2022-23 | Real event 2025 |
|---|---:|---:|
| Rows (stations) | 432 (280) | 463 (297) |
| Reach > 0 | 13.0% | 40.2% |
| Overtopped (reach >= 20 m) | 5.6% | 18.1% |
| Mean / median / 90th pct reach distance (m) | 2.5 / 0.0 / 10.0 | 7.7 / 0.0 / 20.0 |
| Mean reach height (m) | 1.3 | 3.4 |
| Censored | 3.5% | 14.3% |
| Impact in the combined label, whole corridor (ha) | 63.2 | 119.0 |
| Rows left out: near-bank nodata / river not visible / low valid fraction | 103 / 36 / 258 | 89 / 28 / 221 |

By label variant (analysis set of each run):

| Label | Placebo overtopped | Placebo mean reach (m) | Real overtopped | Real mean reach (m) |
|---|---:|---:|---:|---:|
| combined | 5.6% | 2.5 | 18.1% | 7.7 |
| optical_only | 3.0% | 1.4 | 14.2% | 6.2 |
| sar_near | 2.0% | 0.8 | 9.2% | 3.6 |
| sar_only | 2.7% | 1.2 | 5.3% | 1.9 |
| sar_strict | 0.2% | 0.1 | 0.0% | 0.1 |

Figures: `fig_along_river.png` (share overtopped per river km, both banks), `fig_reach_distribution.png`, `fig_detections_by_km.png`, `sheet_placebo_detections_*.png`.

## What the placebo label detected (false positives)

Every impact patch in the placebo combined label is a false positive in the sense that no notable flood occurred. Patches are sorted into diagnostic categories by the screening rules in `placebo_diagnostics` (config.yaml). The categories do not decide what is impact and do not say what caused a patch; agricultural use and construction cannot be separated at 10 m. Where a patch lies within 100 m of a visually reviewed spot (`results/placebo_check/visual_review.csv`), the review is shown.

| Category | Patches | Area (ha) |
|---|---:|---:|
| Channel margin: ordinary seasonal river change or possible flood-related change | 142 | 19.29 |
| Floodplain or terrace: agricultural, construction or other land use (not separable at 10 m) | 80 | 16.78 |
| Hillslope: landslide or gully (not river flooding) | 43 | 8.32 |

Largest detections (all in `placebo_detections.csv` with coordinates):

| km | Bank | Area (m2) | Distance to channel (m) | Height above channel (m) | Optical / radar share | Category |
|---:|---|---:|---:|---:|---|---|
| 3.59 | left | 21700 | 10.0 | 14.3 | 0.56 / 0.60 | floodplain_terrace |
| 56.26 | right | 18300 | 36.1 | 16.3 | 0.07 / 0.95 | floodplain_terrace |
| 51.71 | left | 12800 | 22.4 | 5.5 | 0.98 / 0.30 | channel_margin |
| 68.02 | left | 11000 | 28.3 | 14.8 | 0.65 / 0.45 | floodplain_terrace |
| 23.71 | right | 10800 | 10.0 | 42.1 | 0.74 / 0.30 | hillslope |
| 29.95 | right | 10500 | 10.0 | 17.1 | 0.99 / 0.01 | floodplain_terrace |
| 67.65 | right | 10400 | 14.1 | 4.9 | 0.94 / 0.08 | channel_margin |
| 54.55 | right | 8400 | 10.0 | 11.5 | 0.99 / 0.01 | floodplain_terrace |
| 70.4 | left | 7900 | 10.0 | 1.5 | 0.70 / 0.33 | channel_margin |
| 61.13 | right | 7000 | 10.0 | 6.6 | 0.61 / 0.39 | channel_margin |
| 36.67 | right | 6700 | 10.0 | 22.1 | 0.93 / 0.12 | hillslope |
| 52.23 | left | 6400 | 10.0 | 16.7 | 0.75 / 0.52 | floodplain_terrace |
| 34.5 | right | 6300 | 36.1 | 4.2 | 1.00 / 0.00 | floodplain_terrace |
| 12.25 | right | 6300 | 10.0 | 92.5 | 1.00 / 0.00 | hillslope |
| 14.19 | left | 6100 | 10.0 | 13.7 | 0.79 / 0.57 | floodplain_terrace |

Analysis rows with a reach above 0 m: 56. The patch at their reach edge, by category:

| Category at the reach edge | Rows |
|---|---:|
| small_patch | 37 |
| channel_margin | 11 |
| floodplain_terrace | 4 |
| hillslope | 4 |

Affected river reaches (river km with any analysis row overtopped in the placebo):

- left: 0.70–0.80, 4.63, 10.64, 13.76, 17.27, 21.48, 23.50, 24.80, 32.11 km
- right: 0.40, 1.21, 2.62, 3.72, 4.42, 4.63, 13.66, 14.56, 15.26, 15.76, 17.67, 18.67, 30.01, 31.81 km

## Rows lost because the river is not visible (14.20)

- Placebo (2022 forest): 18 stations (36 rows); mean forest fraction of those rows 0.72 vs 0.40 in the analysis set. River km: 17.07, 18.27–18.57, 19.47, 22.29, 22.49–22.59, 22.79, 23.40, 24.40, 24.60, 25.61–25.71, 26.11, 27.41, 27.61 km.
- Real event (2024 forest): 14 stations (28 rows); mean forest fraction of those rows 0.56 vs 0.45 in the analysis set. River km: 7.63, 15.06, 18.37, 19.17, 19.47, 20.58–20.68, 22.39, 22.59, 22.79, 25.61–25.71, 27.61–27.71 km.

The placebo drops 11% of High-forest rows for river_not_visible (25 of 233), against 1% of Low-forest rows (2 of 233); Mid 4% (terciles of the 2022 forest over all in-reach rows). The placebo therefore cannot speak for the densest forest.

Rule 14.20 keeps a row only if bare or impacted ground within 30 m of the station shows where the river is. After the 2025 flood the bed was bare almost everywhere; in an ordinary year a narrow river under canopy is often invisible at 10 m. The placebo analysis set therefore has fewer rows, mostly fewer forested gorge rows, than the real event. The 14.20 sensitivity (those rows with reach 0) is in the table at the end of the forest section.

## Forest exposure

Primary (14.2): forest from the placebo pre window (Oct-Dec 2022), with Dynamic World, ETH canopy height 2020 and the Hansen rule (loss years 2021-2022, 14.3). Placebo only: WorldCover 2021 (ESA/WorldCover/v200) stands in for 2022, which has no map, where Dynamic World has no observation (user decision 2026-11-03). Rows using it are flagged.

- Rows using the WorldCover stand-in (`worldcover_standin` = 1): 18 of 1432 rows; 9 of 432 analysis rows (2.1%).
- Forest source per row, 100 m belt: dynamic_world 1414, mixed 14, worldcover 4.
- left bank, river km: 23.90–24.30, 26.11–26.51, 27.21–27.31 km
- right bank, river km: 26.11–26.51, 27.21 km
- Stand-in analysis rows per 2022 forest tercile: Low 1, Mid 1, High 7

Descriptive only; the models are in the `stats` report. Terciles are equal-count groups of the analysis set; ties (many rows have no forest) are split by station order. With so few overtopped rows, differences between terciles are a few rows and should not be read as a trend.

**2022 forest (primary)**

| Tercile | Rows | Forest fraction range | Overtopped rows | Overtopped | Mean reach (m) |
|---|---:|---|---:|---:|---:|
| Low | 144 | 0.00-0.04 | 11 | 7.6% | 3.7 |
| Mid | 144 | 0.04-0.66 | 8 | 5.6% | 2.4 |
| High | 144 | 0.66-1.00 | 5 | 3.5% | 1.4 |

**2022 forest, without stand-in rows**

| Tercile | Rows | Forest fraction range | Overtopped rows | Overtopped | Mean reach (m) |
|---|---:|---|---:|---:|---:|
| Low | 141 | 0.00-0.03 | 11 | 7.8% | 3.7 |
| Mid | 141 | 0.04-0.65 | 7 | 5.0% | 2.3 |
| High | 141 | 0.65-1.00 | 6 | 4.3% | 1.5 |

**2024 forest (sensitivity)**

| Tercile | Rows | Forest fraction range | Overtopped rows | Overtopped | Mean reach (m) |
|---|---:|---|---:|---:|---:|
| Low | 144 | 0.00-0.08 | 11 | 7.6% | 3.5 |
| Mid | 144 | 0.08-0.65 | 7 | 4.9% | 2.5 |
| High | 144 | 0.65-1.00 | 6 | 4.2% | 1.4 |

**real event, 2024 forest**

| Tercile | Rows | Forest fraction range | Overtopped rows | Overtopped | Mean reach (m) |
|---|---:|---|---:|---:|---:|
| Low | 154 | 0.00-0.18 | 37 | 24.0% | 10.1 |
| Mid | 154 | 0.19-0.69 | 32 | 20.8% | 8.3 |
| High | 155 | 0.69-1.00 | 15 | 9.7% | 4.7 |

Stand-in sensitivity, beside the main placebo result:

| Run | Analysis rows | Reach > 0 | Overtopped | Mean reach (m) |
|---|---:|---:|---:|---:|
| placebo, all analysis rows | 432 | 13.0% | 5.6% | 2.5 |
| placebo, without WorldCover stand-in rows | 423 | 13.2% | 5.7% | 2.5 |
| placebo, river-not-visible rows as reach 0 (14.20) | 459 | 12.2% | 5.2% | 2.3 |
| real event, river-not-visible rows as reach 0 (14.20) | 481 | 38.7% | 17.5% | 7.4 |

## Limitations

- The placebo shows how much impact the method finds when no notable flood occurred. It does not prove that no local flood, high flow or bank erosion happened in the period.
- Reported inside the confirmed interval: sudden rise, Bhote Koshi, originating in Tibet, 2023-07-12 to 2023-07-13 (source: Radio Nepal, 2023-07-13): police warnings issued; damaging flood not confirmed. The user's confirmation covers notable floods; this report is listed so a reader can judge it.
- No corridor-specific flood record was found for Oct-Dec 2022 or the rest of 2023.
- No Sentinel-1 scene exists for 2023-07-12 to 2023-07-15 in this corridor, so a short surge could not be seen in radar.
- Small isolated optical changes near the river (km 2.2-2.6, 3.3-3.8, 4.2-4.5, 7.3-7.8, 14.2-14.6, 23.6-23.9, 24.3-24.5, 26.9-27.0, 28.3, 29.8-30.0, 30.0, 33.6-33.8, 34.3-34.5) are consistent with bar, bank, slope or construction change and are not matched in radar.
- In the gorge (km 21-25) the 2023 channel is darker and less bare than 2022 (less bare-soil index), the opposite of flood scour.
- Wayback imagery cannot separate 2023 change from 2018-2024.
- The forest-effect half of the placebo test (beta1, 8.9) is in the `stats` report; it covers only rows whose river was visible, which leaves out much of the densest forest.
- Diagnostic categories use 10 m imagery and a 30 m DEM; they are review aids, not causes. Visual reviews exist only for the spots in the review file.
- The placebo forest uses ETH canopy height 2020 (before the placebo period) and, where Dynamic World has no 2022 image, the WorldCover 2021 stand-in (flagged rows).
- The 2025 validation edges (14.13) belong to the real event and are not applicable here.
- The centerline was traced on OpenStreetMap waterway (osmnx), checked by the user in QGIS against Esri World Imagery Wayback; the placebo channel and bank tops come from the 2022 composite, so distances to the line are approximate in the lower valley.

## Files

- `placebo_rows.csv`
- `placebo_detections.csv`
- `placebo_vs_real_summary.csv`
- `placebo_standin_sensitivity.csv`
- `placebo_forest_terciles.csv`
- `placebo_vs_real_per_km.csv`
- `placebo_vs_real_by_label.csv`
- `fig_along_river.png`
- `fig_reach_distribution.png`
- `fig_detections_by_km.png`
- `sheet_placebo_detections_km00.0-12.0.png`
- `sheet_placebo_detections_km12.0-24.0.png`
- `sheet_placebo_detections_km24.0-36.0.png`
- `sheet_placebo_detections_km36.0-48.0.png`
- `sheet_placebo_detections_km48.0-60.0.png`
- `sheet_placebo_detections_km60.0-71.8.png`

