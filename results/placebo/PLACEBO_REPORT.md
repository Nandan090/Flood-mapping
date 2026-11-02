# Placebo report: melamchi2021

Reproduce with `python -m floodforest.cli placebo --event melamchi2021`.

## What was run

The real event's pipeline (fetch, label, transects, profiles, table) was run unchanged as `melamchi2021_placebo`, a no-flood control (spec 8.9, 14.2). Same centerline, reach, manual exclusions, thresholds and label rules as the real event; only the image windows differ.

| Input | Placebo window |
|---|---|
| Sentinel-2 pre (forest, channel, controls) | 2018-10-01 to 2018-12-15 |
| Sentinel-2 post | 2019-10-01 to 2019-12-15 |
| Sentinel-1 pre / post | 2018-08-03 to 2018-09-15 / 2019-08-03 to 2019-09-15 (end exclusive) |
| sar_near (sensitivity) | 2019-04-01 to 2019-06-15 |
| Fake event date | 2019-06-15 (pre cutoff 2019-06-15T00:00:00Z, post cutoff 2019-08-03T00:00:00Z) |

The user confirmed no notable flood in 2018-08-03 to 2019-12-15 (`flood_check`, 2026-11-02; record `results/placebo_period_check.md`). That confirmation covers notable floods only; local floods, high flows and bank erosion are not ruled out.

The real event's files were only read. The real event's Oct-Dec 2020 forest baseline is unchanged; it is used here only as a labelled sensitivity.

## Expected outcome and result

Expected with no flood (8.9): near-zero reach. The spec sets no numeric pass mark for 'near-zero', so the placebo is shown beside the real event. The second part of 8.9 (beta1 not distinguishable from zero) needs the `stats` step, which is not built yet.

| Analysis set, combined label, FABDEM | Placebo 2018-19 | Real event 2021 |
|---|---:|---:|
| Rows (stations) | 339 (199) | 543 (304) |
| Reach > 0 | 5.3% | 82.5% |
| Overtopped (reach >= 20 m) | 1.8% | 70.3% |
| Mean / median / 90th pct reach distance (m) | 1.0 / 0.0 / 0.0 | 53.4 / 40.0 / 130.0 |
| Mean reach height (m) | 0.3 | 6.7 |
| Censored | 1.2% | 17.3% |
| Impact in the combined label, whole corridor (ha) | 7.0 | 370.0 |
| Rows left out: near-bank nodata / river not visible / low valid fraction | 51 / 226 / 86 | 31 / 16 / 56 |

By label variant (analysis set of each run):

| Label | Placebo overtopped | Placebo mean reach (m) | Real overtopped | Real mean reach (m) |
|---|---:|---:|---:|---:|
| combined | 1.8% | 1.0 | 70.3% | 53.4 |
| optical_only | 1.5% | 0.8 | 68.7% | 50.6 |
| sar_near | 0.6% | 0.2 | 13.6% | 8.0 |
| sar_only | 0.3% | 0.2 | 20.7% | 13.5 |
| sar_strict | 0.0% | 0.0 | 6.7% | 4.9 |

Figures: `fig_along_river.png` (share overtopped per river km, both banks), `fig_reach_distribution.png`, `fig_detections_by_km.png`, `sheet_placebo_detections_*.png`.

## What the placebo label detected (false positives)

Every impact patch in the placebo combined label is a false positive in the sense that no notable flood occurred. Patches are sorted into diagnostic categories by the screening rules in `placebo_diagnostics` (config.yaml). The categories do not decide what is impact and do not say what caused a patch; agricultural use and construction cannot be separated at 10 m. Where a patch lies within 100 m of a visually reviewed spot (`results/placebo_check/visual_review.csv`), the review is shown.

| Category | Patches | Area (ha) |
|---|---:|---:|
| Channel margin: ordinary seasonal river change or possible flood-related change | 11 | 1.49 |
| Floodplain or terrace: agricultural, construction or other land use (not separable at 10 m) | 10 | 1.15 |
| Hillslope: landslide or gully (not river flooding) | 10 | 1.19 |

Detections that match a visually reviewed spot:

| km | Bank | Area (m2) | Category | Visual review | Reviewer |
|---:|---|---:|---|---|---|
| 0.24 | right | 600 | floodplain_terrace | ordinary_seasonal_river_change (2018_2019 patch 1: edge of a bare fan that was already there; mostly before the 2019 monsoon) | claude |
| 0.29 | right | 2600 | channel_margin | ordinary_seasonal_river_change (2018_2019 patch 1: edge of a bare fan that was already there; mostly before the 2019 monsoon) | claude |
| 0.29 | left | 600 | channel_margin | ordinary_seasonal_river_change (2018_2019 patch 1: edge of a bare fan that was already there; mostly before the 2019 monsoon) | claude |
| 1.53 | right | 600 | hillslope | landslide_or_gully (monsoon2018 patch 21: hillside gully about 160 m above the river) | claude+independent |
| 5.51 | right | 1700 | channel_margin | uncertain (2018_2019 patch 4: user: minor ambiguity, no convincing flood deposit) | user |
| 16.89 | left | 2300 | channel_margin | agricultural (2018_2019 patch 2: fields bare in the dry season) | claude |
| 17.37 | left | 2100 | floodplain_terrace | uncertain (2018_2019 patch 3: user: no clear flood signal) | user |
| 20.74 | right | 1000 | channel_margin | uncertain (2018_2019 patch 6: user: minor ambiguity, no convincing flood deposit) | user |
| 20.9 | right | 1100 | channel_margin | ordinary_seasonal_river_change (2018_2019 patch 5: user: no clear flood signal; narrow channel-edge change) | user |

Largest detections (all in `placebo_detections.csv` with coordinates):

| km | Bank | Area (m2) | Distance to channel (m) | Height above channel (m) | Optical / radar share | Category |
|---:|---|---:|---:|---:|---|---|
| 1.77 | right | 2700 | 40.0 | 9.1 | 0.00 / 1.00 | floodplain_terrace |
| 0.29 | right | 2600 | 10.0 | -5.5 | 1.00 / 0.04 | channel_margin |
| 16.89 | left | 2300 | 10.0 | -0.3 | 1.00 / 0.00 | channel_margin |
| 33.07 | right | 2300 | 22.4 | 2.1 | 0.96 / 0.04 | channel_margin |
| 33.09 | right | 2300 | 89.4 | 10.3 | 0.78 / 0.22 | floodplain_terrace |
| 33.29 | right | 2200 | 70.0 | 30.3 | 0.68 / 0.45 | hillslope |
| 17.37 | left | 2100 | 36.1 | 12.0 | 1.00 / 0.00 | floodplain_terrace |
| 5.51 | right | 1700 | 10.0 | 3.4 | 0.94 / 0.06 | channel_margin |
| 18.74 | left | 1700 | 143.2 | 26.2 | 0.06 / 0.94 | hillslope |
| 33.14 | right | 1700 | 102.0 | 30.8 | 0.88 / 0.12 | hillslope |
| 33.34 | left | 1600 | 30.0 | 0.2 | 0.12 / 0.88 | channel_margin |
| 11.83 | right | 1400 | 10.0 | 53.2 | 0.57 / 0.50 | hillslope |
| 0.91 | right | 1200 | 10.0 | 101.7 | 0.00 / 1.00 | hillslope |
| 20.9 | right | 1100 | 10.0 | -0.7 | 1.00 / 0.00 | channel_margin |
| 2.9 | right | 1100 | 40.0 | 29.0 | 0.00 / 1.00 | hillslope |

Analysis rows with a reach above 0 m: 18. The patch at their reach edge, by category:

| Category at the reach edge | Rows |
|---|---:|
| small_patch | 11 |
| channel_margin | 4 |
| floodplain_terrace | 2 |
| hillslope | 1 |

Affected river reaches (river km with any analysis row overtopped in the placebo):

- left: 11.20, 16.97, 20.61 km
- right: 5.51, 10.39, 20.81 km

## Rows lost because the river is not visible (14.20)

- Placebo (2018 forest): 113 stations (226 rows); mean forest fraction of those rows 0.50 vs 0.36 in the analysis set. River km: 0.00–0.10, 1.13–1.54, 1.74, 1.95–2.05, 2.26, 2.46, 2.67–2.77, 2.97–3.28, 3.48, 3.79–3.89, 4.10, 4.30–4.60, 4.80, 5.01–5.41, 5.61, 5.82–6.43, 6.63–7.14, 7.65, 7.85–7.95, 8.16–8.26, 8.86–8.97, 9.27–9.38, 9.59, 10.09–10.19, 10.49–11.00, 11.40–11.50, 12.42, 13.94, 14.75, 15.76, 16.37–16.47, 16.87, 17.48, 17.68, 17.98, 18.69, 18.89, 19.10–19.50, 19.70, 20.10, 21.82–22.22, 23.12–23.22, 23.43–23.83, 24.03–24.33, 25.75–25.85, 27.26, 27.77–27.87, 28.17–28.27, 29.48, 33.31, 34.41, 34.61–34.72 km.
- Real event (2020 forest): 8 stations (16 rows); mean forest fraction of those rows 0.82 vs 0.39 in the analysis set. River km: 1.44, 2.26, 2.77, 5.11, 5.82, 6.03–6.13, 10.59 km.

Rule 14.20 keeps a row only if bare or impacted ground within 30 m of the station shows where the river is. After the 2021 flood the bed was bare almost everywhere; in an ordinary year a narrow river under canopy is often invisible at 10 m. The placebo analysis set therefore has fewer rows, mostly fewer forested gorge rows, than the real event. The 14.20 sensitivity (those rows with reach 0) is in the table at the end of the forest section.

## Forest exposure

Primary (14.2): forest from the placebo pre window (Oct-Dec 2018), with Dynamic World, ETH canopy height 2020 and the Hansen rule (skipped for 2018, 14.3). Placebo only: WorldCover 2020 (ESA/WorldCover/v100) stands in for 2018, which has no map, where Dynamic World has no observation (user decision 2026-11-02). Rows using it are flagged.

- Rows using the WorldCover stand-in (`worldcover_standin` = 1): 0 of 704 rows; 0 of 339 analysis rows (0.0%). The stand-in was not needed: every decided belt pixel had a Dynamic World value.
- Forest source per row, 100 m belt: dynamic_world 704.
- left bank, river km: none
- right bank, river km: none
- Stand-in analysis rows per 2018 forest tercile: Low 0, Mid 0, High 0

Descriptive only; no model is fitted (the `stats` step is not built). Terciles are equal-count groups of the analysis set; ties (many rows have no forest) are split by station order. With so few overtopped rows, differences between terciles are a few rows and should not be read as a trend.

**2018 forest (primary)**

| Tercile | Rows | Forest fraction range | Overtopped rows | Overtopped | Mean reach (m) |
|---|---:|---|---:|---:|---:|
| Low | 113 | 0.00-0.00 | 3 | 2.7% | 1.4 |
| Mid | 113 | 0.00-0.63 | 2 | 1.8% | 1.1 |
| High | 113 | 0.63-1.00 | 1 | 0.9% | 0.6 |

**2018 forest, without stand-in rows**

| Tercile | Rows | Forest fraction range | Overtopped rows | Overtopped | Mean reach (m) |
|---|---:|---|---:|---:|---:|
| Low | 113 | 0.00-0.00 | 3 | 2.7% | 1.4 |
| Mid | 113 | 0.00-0.63 | 2 | 1.8% | 1.1 |
| High | 113 | 0.63-1.00 | 1 | 0.9% | 0.6 |

**2020 forest (sensitivity)**

| Tercile | Rows | Forest fraction range | Overtopped rows | Overtopped | Mean reach (m) |
|---|---:|---|---:|---:|---:|
| Low | 113 | 0.00-0.00 | 4 | 3.5% | 1.9 |
| Mid | 113 | 0.00-0.65 | 1 | 0.9% | 0.7 |
| High | 113 | 0.65-1.00 | 1 | 0.9% | 0.5 |

**real event, 2020 forest**

| Tercile | Rows | Forest fraction range | Overtopped rows | Overtopped | Mean reach (m) |
|---|---:|---|---:|---:|---:|
| Low | 181 | 0.00-0.00 | 165 | 91.2% | 85.4 |
| Mid | 181 | 0.00-0.70 | 145 | 80.1% | 57.3 |
| High | 181 | 0.70-1.00 | 72 | 39.8% | 17.6 |

Stand-in sensitivity, beside the main placebo result:

| Run | Analysis rows | Reach > 0 | Overtopped | Mean reach (m) |
|---|---:|---:|---:|---:|
| placebo, all analysis rows | 339 | 5.3% | 1.8% | 1.0 |
| placebo, without WorldCover stand-in rows | 339 | 5.3% | 1.8% | 1.0 |
| placebo, river-not-visible rows as reach 0 (14.20) | 524 | 3.4% | 1.1% | 0.7 |
| real event, river-not-visible rows as reach 0 (14.20) | 554 | 80.9% | 69.0% | 52.4 |

## Limitations

- The placebo shows how much impact the method finds when no notable flood occurred. It does not prove that no local flood, high flow or bank erosion happened in the period.
- The forest-effect half of the placebo test (beta1, 8.9) waits for the `stats` step.
- Diagnostic categories use 10 m imagery and a 30 m DEM; they are review aids, not causes. Visual reviews exist only for the spots in the review file.
- The placebo forest uses ETH canopy height 2020 (after the placebo period) and, where Dynamic World has no 2018 image, the WorldCover 2020 stand-in (flagged rows).
- The 2021 validation edges (14.25) belong to the real event and are not applicable here.
- The centerline was traced on post-2021 imagery; the placebo channel and bank tops come from the 2018 composite, so distances to the line are approximate in the lower valley.

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
- `sheet_placebo_detections_km24.0-35.5.png`

