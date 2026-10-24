# Build spec: Do riverbank forests limit the lateral reach of Himalayan debris floods?

This file is the complete instruction set for the AI coding assistant. Read all of it before writing code. Build it in the phase order given at the end. Where this file and `research-premise.md` disagree, **this file wins** (see section 1.3).

---

## 0. Rules for the coding assistant (non-negotiable)

1. **Never fabricate.** No synthetic numbers, placeholder results, mock imagery or hardcoded statistics anywhere outside `tests/`. If data is missing, raise an error that names what is missing.
2. **No silent fallbacks.** If a dataset, band or asset is unavailable, stop and report. Do not substitute another dataset without a config change made by the user.
3. **Verify every Earth Engine asset ID at runtime** (`ee.data.getInfo` / `.first().bandNames().getInfo()`), and print the ID, band names and date range actually used. Asset IDs in this file were correct to the author's knowledge but must be checked; if one fails, stop and tell the user rather than guessing a replacement.
4. **Strict time separation.** Exposure (forest) and controls use only pre-event imagery or terrain. Outcomes use post-event imagery. Add assertions on acquisition dates that fail the run if violated.
5. **Vegetation loss is never the outcome.** Do not use NDVI drop as the impact label (section 5 explains why).
6. **Everything is config-driven** (`config.yaml`). No dates, thresholds, paths or coordinates hardcoded in modules.
7. **Deterministic.** One global random seed in config, passed to every bootstrap and sampler.
8. **Cache every expensive step** to `data/` and make each step re-runnable alone via the CLI.
9. **Log provenance:** every scene ID and acquisition date used in every composite goes into `outputs/{event}/provenance.json` and into the report.
10. **Do not tune thresholds to get a significant result.** Thresholds are set by the rules in section 5 and then frozen. All sensitivity runs are reported whether or not they agree.
11. Ask the user before changing any analytical decision in this file.

---

## 1. What the study is

### 1.1 Question

When a very large debris flood comes down a Himalayan valley, is pre-event forest on the riverbank associated with **less lateral reach** of the flow on that bank, after accounting for terrain?

### 1.2 Roles of the variables

| Role                                    | Variable                                                                                                                                | Source timing       |
| --------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- | ------------------- |
| **Exposure** (the "cause" being tested) | Forest belt fraction on the bank                                                                                                        | Pre-event only      |
| **Outcome** (what the flood did)        | Reach distance (m), reach height (m), overtopped (yes/no)                                                                               | Post-event          |
| **Controls**                            | Channel width, gradient, valley width, wall slope, bank slope, upstream area, bend side, curvature, bank elevation, position downstream | Terrain / pre-event |

### 1.3 Warning about `research-premise.md`

The first half of that file (goal, hypothesis, transect method, four traps) is correct. The second half (the tutorial beginning "Absolutely. These concepts fit together…") explains the individual concepts correctly but **frames the study backwards** in its statistics sections and final diagram: it treats _flood exposure_ as the cause and _vegetation change_ as the outcome, and suggests sampling NDVI change along transects as the result. Do **not** implement that. Here forest is the exposure, flood reach is the outcome, and vegetation loss is explicitly banned as an outcome.

### 1.4 Events

| Key              | Event                                           | Event date(s)                                            | Reach                                                                                | Status                                                                                                                                                                                     |
| ---------------- | ----------------------------------------------- | -------------------------------------------------------- | ------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `melamchi2021`   | Melamchi landslide-dam flood                    | 2021-06-15, second surge 2021-07-31 (treat as one event) | ~40 km, Melamchi River → Indrawati corridor                                          | All data exists. **Build and debug on this first.**                                                                                                                                        |
| `lende2025`      | Lende Khola lake-outburst flood                 | 2025-07-08                                               | ~35 km, Lende/Lhende Khola → Bhotekoshi → Trishuli corridor, starting at Rasuwagadhi | All data exists                                                                                                                                                                            |
| `bhotekoshi2026` | Ice–rock-avalanche-triggered debris/flash flood | 2026-08-26                                               | ~72 km, Bhotekoshi → Trishuli corridor, Rasuwagadhi → Bidur                          | **Main test.** Post-event Oct–Dec Sentinel-2 is still accumulating; the pipeline must be re-runnable as scenes arrive and must print how many clear post-event observations each pixel has |

Exact start and end points, and how river confluences inside a reach are handled, are fixed in section 14.16.

Events are analysed **separately with identical code**. Never pool them into one model.

---

## 2. Environment

- Python 3.11 or newer, virtual environment, `pip install -r requirements.txt`, then `pip freeze > requirements.lock.txt`.
- Projected CRS for all geometry and rasters: **EPSG:32645** (UTM 45N). All distances in metres. Assert CRS on every load.
- Analysis grid: **10 m**, snapped to the Sentinel-2 grid.

### 2.1 Libraries

```
earthengine-api   # Earth Engine Python client
geemap            # EE helpers, image download
geedim            # tiled GeoTIFF export from EE
numpy
pandas
scipy
xarray
rioxarray
rasterio
geopandas
shapely
pyproj
rasterstats       # zonal statistics for belts
scikit-image      # connected components, morphology
statsmodels       # regression models
patsy             # formulas and splines
lifelines         # censored (AFT) regression
matplotlib
jinja2            # HTML report
pyyaml
tqdm
osmnx             # draft river line from OpenStreetMap
pytest
# optional
earthaccess       # only if using the NASA HMA 8 m DEM
```

### 2.2 Accounts (all free; no paid API keys)

| Service              | Needed for                | What the user must do                                                                                                                                                                                                                                             |
| -------------------- | ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Google Earth Engine  | Almost all data           | Google account, a Google Cloud project registered for **noncommercial** Earth Engine use (Community tier). A billing account must be attached but is not charged for noncommercial use. Run `earthengine authenticate` once; put the project ID in `config.yaml`. |
| NASA Earthdata Login | Optional HMA 8 m DEM only | Free account                                                                                                                                                                                                                                                      |
| OpenStreetMap        | Draft centerline          | None                                                                                                                                                                                                                                                              |

The code must call `ee.Initialize(project=cfg.ee_project)` and fail with a clear message if not authenticated. No credentials in the repo.

---

## 3. Project layout

```
config.yaml
requirements.txt
src/floodforest/
  cli.py            # python -m floodforest.cli <step> --event <key>
  config.py
  ee_data.py        # all Earth Engine access
  sar.py            # Sentinel-1 preprocessing, layover/shadow mask
  optical.py        # Sentinel-2 masking, composites, indices
  label.py          # impact label
  centerline.py     # river line, stations, transects
  profiles.py       # sampling rasters along transects, reach extraction
  exposure.py       # forest belt metrics
  controls.py
  stats.py          # models, bootstrap, overlap, E-value
  validate.py       # comparison with published extents
  figures.py
  report.py
tests/
data/{event}/       # cached rasters, centerline, tables
outputs/{event}/    # transects.csv, figures, report.html, provenance.json
outputs/summary_report.html
```

CLI steps, in order: `fetch`, `label`, `transects`, `profiles`, `table`, `stats`, `validate`, `placebo`, `report`. Plus `all`.

---

## 4. Data

All fetched through Earth Engine, clipped to a corridor polygon (centerline buffered by 2,000 m), exported as 10 m GeoTIFFs in EPSG:32645 using `geedim`/`geemap.download_ee_image`.

| Purpose                | Dataset                                               | Asset ID (verify at runtime)                                                   |
| ---------------------- | ----------------------------------------------------- | ------------------------------------------------------------------------------ |
| Radar                  | Sentinel-1 GRD, IW, VV+VH                             | `COPERNICUS/S1_GRD`                                                            |
| Optical                | Sentinel-2 L2A harmonized                             | `COPERNICUS/S2_SR_HARMONIZED`                                                  |
| Cloud mask             | Cloud Score+                                          | `GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED`                                     |
| Main DEM (bare earth)  | FABDEM                                                | `projects/sat-io/open-datasets/FABDEM`                                         |
| Sensitivity DEM        | Copernicus GLO-30                                     | `COPERNICUS/DEM/GLO30_2024_1` (replaces the deprecated `COPERNICUS/DEM/GLO30`) |
| Pre-event land cover   | Dynamic World                                         | `GOOGLE/DYNAMICWORLD/V1`                                                       |
| Land cover cross-check | ESA WorldCover                                        | `ESA/WorldCover/v100` (2020), `ESA/WorldCover/v200` (2021)                     |
| Canopy height          | ETH Global Canopy Height 2020                         | `users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1`                               |
| Forest loss history    | Hansen Global Forest Change                           | latest `UMD/hansen/global_forest_change_*` (look up newest version)            |
| Upstream area          | MERIT Hydro, band `upa`                               | `MERIT/Hydro/v1_0_1`                                                           |
| Validation (2026)      | Vantor open data                                      | `projects/sat-io/open-datasets/VANTOR-DISASTER-DATA/NEPAL_FLOOD_2026`          |
| Validation (2026)      | Planet disaster data                                  | `projects/sat-io/open-datasets/PLANET-DISASTER-DATA/NEPAL_FLOOD_2026`          |
| Validation (2026)      | Copernicus EMS activation **EMSR927** vector products | Downloaded manually by the user into `data/bhotekoshi2026/validation/`         |

FABDEM and both validation collections are licensed for non-commercial use (CC BY-NC); record licences and citations in the report.

DEMs are 30 m. Resample to the 10 m grid with **bilinear** interpolation for elevation; never nearest-neighbour. Categorical rasters use nearest; probabilities use bilinear.

### 4.1 Time windows (in `config.yaml`, editable)

| Event            | S2 pre composite         | S2 post composite        | S1 pre                   | S1 post                                   |
| ---------------- | ------------------------ | ------------------------ | ------------------------ | ----------------------------------------- |
| `melamchi2021`   | 2020-10-01 to 2020-12-15 | 2021-10-01 to 2021-12-15 | 2020-08-03 to 2020-09-15 | 2021-08-03 to 2021-09-15                  |
| `lende2025`      | 2024-10-01 to 2024-12-15 | 2025-10-01 to 2025-12-15 | 2024-07-09 to 2024-08-15 | 2025-07-09 to 2025-08-15                  |
| `bhotekoshi2026` | 2025-10-01 to 2025-12-15 | 2026-10-01 to 2026-12-15 | 2025-08-26 to 2025-09-30 | after 2026-08-26 14:52 UTC, to 2026-09-30 |

Window boundaries are enforced on scene **timestamps**, not dates; see section 14.1.

Reasons: pre and post optical composites are the **same season one year apart**, which cancels seasonal greenness and crop-cycle differences. The 2026 pre-composite is taken after the 2025 flood on purpose, so 2025 scars count as non-forest in 2026; state this in the report.

The `fetch` step must list every post-window scene by date so the user can exclude windows contaminated by a later flood (`exclude_dates` in config).

---

## 5. The impact label (most error-prone part)

Goal: a 10 m raster per event, `impact ∈ {1 impacted, 0 not impacted, nodata}`, marking where the debris flow left its mark, built so that detection does not depend on whether trees stood there.

### 5.1 Why not vegetation loss

If impact = NDVI drop, only vegetated pixels can register impact, and forest has the most NDVI to lose. That bakes the answer into the label.

### 5.2 Optical component (Sentinel-2)

1. Mask each scene with Cloud Score+ (`cs_cdf >= 0.60`), then take the per-band **median** over the window. Export also the per-pixel count of clear observations; pixels with count < 2 are nodata.
2. Indices on pre and post composites (reflectance scaled 0–1):
   - `NDVI = (B8 − B4) / (B8 + B4)`
   - `BSI = ((B11 + B4) − (B8 + B2)) / ((B11 + B4) + (B8 + B2))`
   - `MNDWI = (B3 − B11) / (B3 + B11)`
   - `BRIGHT = (B2 + B3 + B4) / 3`
3. **Pre-event active channel** = pixels with `MNDWI_pre > 0` or `NDVI_pre < 0.20`, kept only if 8-connected to the centerline. This is water plus bare gravel bars.
4. **Unknown mask:** pixels outside the active channel that were already bare or built before the event (`NDVI_pre < 0.30`) are set to nodata. A deposit there cannot be distinguished from what was already bare.
5. **Post-event deposit state:** `NDVI_post < 0.30 AND BSI_post > t_bsi`, where `t_bsi` is the 5th percentile of `BSI_post` inside the pre-event active channel (a surface known to be bare sediment). Compute per event, write to provenance, then freeze.
6. `optical_impact = deposit state AND not unknown AND not active channel`.

### 5.3 Radar component (Sentinel-1)

1. Process each **relative orbit and pass direction separately**. Never difference images from different orbits.
2. Convert to dB, apply a 3×3 focal median, take the temporal median of the pre window and of the post window (require ≥ 3 scenes each; otherwise that orbit is skipped and logged).
3. Change: `ΔVV = VV_post − VV_pre`, `ΔVH = VH_post − VH_pre` (dB).
4. **Layover and shadow mask** from the DEM and the radar geometry, per orbit. Port the terrain-geometry method from the open-source `gee_s1_ard` project (Mullissa et al. 2021; based on Vollrath et al. 2020) and credit it. Masked pixels are nodata for that orbit.
5. Standardise change against a **stable reference zone**: corridor pixels more than 150 m above the channel, with no optical impact. Robust z-score:
   `z = (Δ − median_ref) / (1.4826 × MAD_ref)`
6. `sar_impact = |z_VV| > 3 OR |z_VH| > 3` in at least one unmasked orbit. If more than one orbit is valid, require agreement in at least half of the valid orbits.

### 5.4 Combining

- `impact_combined = optical_impact OR sar_impact`; nodata only where both are nodata.
- Keep three labels throughout: `combined` (primary), `optical_only`, `sar_only`.
- **Connectivity filter:** keep only impact pixels in components 8-connected to the pre-event active channel (allow bridging gaps of 1 pixel with a morphological closing, radius 1). This removes harvested fields, landslides and other change unrelated to the river.
- Record per pixel which source fired (`label_source`).

### 5.5 Known label biases that must be reported

- **Canopy concealment:** a flow that passed under surviving trees may be invisible to both sensors, which would make forested banks look less impacted than they were. This biases _toward_ a protective finding. Section 9 tests for it; the report must state it.
- **Cropland noise:** fields that are bare after harvest can mimic deposits. The same-season composites and the connectivity filter reduce but do not remove this.

---

## 6. Centerline, stations and transects

1. **Centerline.** Load `data/{event}/centerline.geojson`. A helper (`centerline.py draft`) fetches the OSM waterway with `osmnx` as a starting draft, but the user must check and correct it in QGIS against the pre-event Sentinel-2 composite before it is used. The code refuses to run if the file lacks the property `"verified_by_user": true`.
2. Reproject to EPSG:32645, smooth (Chaikin, 3 iterations, or a 200 m moving average), densify to 10 m vertices, orient upstream to downstream.
3. **Stations** every `station_spacing_m` (default **100**) along the line. `river_km` = distance downstream from the start. `segment_id = floor(river_km)` (1 km blocks).
4. **Local direction** at each station = unit vector between the points 100 m upstream and 100 m downstream on the smoothed line. **Normal** = that vector rotated 90°: for direction `(dx, dy)`, left normal is `(−dy, dx)`, right normal is `(dy, −dx)`. "Left" and "right" are defined looking downstream.
5. **Transect** = line from the station outward along each normal to `max_half_length_m` (default 1,500), sampled every 10 m.
6. Flag and exclude from analysis (but keep in the table with a reason code):
   - transects that intersect a neighbouring transect before their reach edge,
   - transects within 300 m of a tributary confluence where MERIT `upa` of the tributary exceeds 5% of the mainstem's (tributary fans carry their own debris),
   - transects that cross the river a second time (meander loops),
   - reservoirs and dam sites supplied by the user in `exclusions.geojson`.

---

## 7. Per-transect measurements

Each row of the output table is one **station × side**. Sample every raster along the transect at 10 m.

### 7.1 Geometry

- `bank_top`: walking outward from the station, the first sample after which 3 consecutive samples are **not** pre-event active channel.
- `channel_width_m`: distance between left and right `bank_top`.
- `z0` (channel elevation): median DEM elevation of samples inside the active channel on that transect.
- `valley_wall`: first sample beyond `bank_top` where DEM slope ≥ 30° for ≥ 50 m continuously. `valley_width_m` = distance from `bank_top` to `valley_wall` (or to profile end, flagged).

### 7.2 Outcomes (post-event)

Walking outward from `bank_top` along the chosen label:

- The impacted run continues across gaps of non-impact ≤ 20 m. Nodata gaps ≤ 30 m do not end the run; a longer nodata gap ends it and sets `censored = 1`.
- `reach_edge` = last impacted sample of that run.
- `reach_distance_m` = distance from `bank_top` to `reach_edge` (0 if none).
- `reach_height_m` = `max(0, z(reach_edge) − z0)`.
- `overtopped` = 1 if `reach_distance_m ≥ 20`.
- `censored` = 1 if the run reaches the valley wall, the profile end, or a long nodata gap. The true reach is then "at least this much".
- `reach_height_unc_m` = `10 × tan(local slope at reach_edge)`; one pixel of horizontal error converted to height.
- `valid_fraction` = share of samples between `bank_top` and `valley_wall` with usable data. Rows with `valid_fraction < 0.7` are excluded from models.

### 7.3 Exposure (pre-event only)

Belt = rectangle centred on the transect, `station_spacing_m` long (along river) × `W` wide (from `bank_top` outward). Default **W = 100 m**; also 50 and 200.

- **Forest pixel** = mean Dynamic World `trees` probability over the S2 pre window ≥ 0.5 **and** ETH canopy height ≥ 5 m **and** no Hansen loss between 2021 and the pre-event year.
- `forest_frac_W` = forest share of valid belt pixels (via `rasterstats`). This is the exposure `F`.
- Also record: `canopy_height_mean`, `ndvi_pre_mean`, `belt_continuity` (largest forest patch area ÷ total forest area in belt), `prior_loss_frac`.
- Cross-check column: `forest_frac_worldcover` (class 10).

### 7.4 Controls

`channel_width_m`; `gradient` (slope of `z0` vs `river_km`, linear fit over ±500 m); `valley_width_m`; `wall_slope_deg`; `bank_slope_deg` (mean DEM slope in the belt); `log10_upa` (max MERIT `upa` within 300 m of the station); `curvature` (signed, over ±300 m); `outer_bend` (1 if this side is the outside of the bend); `bank_elev_m`; `river_km`; `side`.

### 7.5 Run everything for both DEMs

All DEM-derived columns are computed twice, suffix `_fab` (primary) and `_cop` (sensitivity).

---

## 8. Statistics

Analysis set: rows that pass exclusions, `valid_fraction ≥ 0.7`, inside the mapped impact reach.

### 8.1 Common support (do this before any model)

1. `forest_high = F ≥ 0.5`.
2. Logistic regression of `forest_high` on all controls gives a propensity score `p`.
3. Keep rows with `p` inside `[max(min_p_forest, min_p_nonforest), min(max_p_forest, max_p_nonforest)]`.
4. Report how many rows were dropped and plot both `p` distributions. If fewer than 30% of rows or fewer than 15 segments remain, the report must say the comparison is too thin to interpret.

### 8.2 Model A: between-transect

For outcome `y`:

`y = β0 + β1·F + β2·log(channel_width) + β3·gradient + β4·log(valley_width) + β5·wall_slope + β6·bank_slope + β7·log10_upa + β8·outer_bend + β9·curvature + β10·bank_elev + side + cr(river_km, df=5) + ε`

- `cr(river_km, df=5)` is a natural cubic spline in downstream distance. It is essential: the flood loses energy downstream and forest cover also changes downstream.
- Fit by OLS for `reach_height_m` and `log1p(reach_distance_m)`; logistic regression for `overtopped`.

### 8.3 Model B: within-transect paired banks (co-primary)

Both banks of one transect met the same flood at the same moment. Using only stations where both sides are valid:

`(y_L − y_R) = α + β1·(F_L − F_R) + γ·Δ(valley_width, wall_slope, bank_slope, outer_bend) + ε`

This removes everything shared by the two banks (discharge, sediment load, position downstream). If A and B disagree in sign, say so prominently.

### 8.4 Uncertainty: block bootstrap

Neighbouring transects are not independent. Resample **1 km segments** with replacement (all rows of a segment move together), refit, repeat `B = 2000` times, seed from config. Report the percentile 95% interval for `β1`. Never report plain OLS standard errors as the headline.

### 8.5 Censoring sensitivity

Refit Model A for reach distance and height with `lifelines.WeibullAFTFitter`, treating `censored = 1` rows as right-censored. Report alongside OLS.

### 8.6 Effect size reporting

For each model give `β1` as: change in outcome going from 0% to 100% forest belt, and per interquartile range of `F`, with the bootstrap interval, the number of rows, stations and segments.

### 8.7 E-value (how strong a hidden confounder would need to be)

1. Standardised effect `d = β1 / SD(y)`.
2. Approximate risk ratio `RR = exp(0.91 × d)`; if `RR < 1` use `RR* = 1/RR`, else `RR* = RR`.
3. `E = RR* + sqrt(RR* × (RR* − 1))`.
4. Repeat for the confidence limit nearest the null (if the interval includes 0, that E-value is 1).
   Explain in the report: an E-value of 2 means an unmeasured factor would need to be associated with both forest and flood reach by a risk ratio of 2 each to erase the result.

### 8.8 Required sensitivity runs (all reported in one table)

| Dimension            | Variants                                      |
| -------------------- | --------------------------------------------- |
| Label                | combined, optical_only, sar_only              |
| DEM                  | FABDEM, Copernicus GLO-30                     |
| Belt width           | 50, 100, 200 m                                |
| Station spacing      | 50, 100 m                                     |
| Forest definition    | Dynamic World rule, WorldCover                |
| Outcome              | height, distance, overtopped                  |
| Steep-terrain filter | all rows; rows with `reach_height_unc_m ≤ 10` |
| Terraces only        | rows with `bank_slope_deg ≤ 15`               |
| Label thresholds     | `t_bsi` and SAR z-threshold ±20%              |

### 8.9 Placebo test

Run the full label and model pipeline on two **pre-event** years with no flood (same season, e.g. Oct–Dec two years before vs one year before the event). Expected: near-zero reach, and `β1` not distinguishable from zero. If the placebo shows a "forest effect", the method is producing artefacts and the report must say the main result cannot be trusted.

### 8.10 Decision rule, fixed in advance

Report the verdict for each event as exactly one of:

- **Consistent with a protective association:** `β1 < 0` with the bootstrap interval excluding 0 in **both** Model A and Model B for reach height, same sign for reach distance, same sign under optical-only and radar-only labels and under both DEMs, placebo null, and validation omission error not materially higher on forested banks (section 9).
- **No evidence of association:** intervals include 0 and are reasonably narrow.
- **Inconclusive:** anything else, including thin overlap, sign flips across sensitivity runs, or a failed placebo.

Across events: tabulate sign, size and interval of `β1` side by side and state whether direction agrees. Do not pool. Wording must always be "associated with", never "protects" or "causes".

---

## 9. Validation of the impact label

Compare `impact_combined` with independent extents (published maps, EMSR927 products, or polygons the user digitises from the Vantor/Planet imagery):

- Intersection over union, precision, recall, overall.
- **The same metrics separately for pixels that were forest and non-forest before the event.** If recall is clearly lower on forest (suggested flag: more than 10 percentage points), the label is missing impact under canopy and the report must warn that the main result is biased toward finding protection.
- Per transect: difference between validated reach distance and our reach distance, plotted against `F`. A trend with `F` is the same warning.

For the Vantor collection: mask nodata (value 0) and use the `cloud_probability` band with a loose threshold (20–30); post-event scenes are very cloudy. Do not difference Planet pre/post analytic bands (different processing levels); use them for visual digitising only.

---

## 10. Output: the report

`outputs/{event}/report.html`, self-contained (figures embedded as base64 PNG), built with Jinja2. Equations via MathJax with a plain-text fallback. Every figure has a caption with three parts: **what it shows, how to read it, what it means here**. Every number has units and the N behind it.

Sections and required figures:

1. **Plain-language summary:** question, verdict from the decision rule, three key numbers, main caveats.
2. **Study area and event:** corridor map with centerline, stations, exclusions.
3. **Data used:** table of every scene ID and date per composite; clear-observation count map; licences.
4. **Imagery:** pre and post Sentinel-2 true-colour and false-colour (B8-B4-B3) panels; pre/post Sentinel-1 VV per orbit; layover/shadow mask map with % of corridor lost.
5. **Building the impact label:** index maps, histograms showing how `t_bsi` and the SAR reference statistics were derived, the three labels side by side, agreement matrix between optical and radar.
6. **Transect method:** diagram of one transect; 6 example cross-section plots (elevation profile with `bank_top`, `z0`, `reach_edge`, valley wall, forest samples and impact samples marked), chosen by fixed rule (stations at the 5th, 25th, 50th, 75th, 95th percentile of reach height plus one censored), not hand-picked.
7. **The data table:** column dictionary, summary statistics, missingness, counts of exclusions by reason, link to `transects.csv`.
8. **Along-river plot:** reach height and forest fraction vs `river_km`, both banks.
9. **Overlap check:** propensity histograms, rows kept and dropped.
10. **Results:** formulas written out; scatter of reach height vs `F` with fitted line (partial residual plot); paired-bank difference plot; coefficient tables with bootstrap intervals; bootstrap distribution histogram for `β1`.
11. **Sensitivity:** forest plot of `β1` across all runs in 8.8; censored-model results.
12. **Placebo and validation:** results of sections 8.9 and 9, including the forest/non-forest recall comparison.
13. **E-value** with the worked calculation.
14. **Limitations,** always including: observational design; canopy concealment; 30 m DEM in a narrow gorge; pre-event DEM vs changed riverbed; unmeasured soil and disturbance history; forest correlating with less settlement; small number of events.
15. **Conclusion:** restates the verdict in two or three sentences, using the fixed wording.

`outputs/summary_report.html` puts the three events' headline rows side by side.

Also write `outputs/{event}/transects.csv`, `transects.gpkg` (stations, transects, reach edges, for QGIS), and all label rasters.

---

## 11. Tests (write these before running real data)

Use small synthetic rasters in `tests/`:

1. **Transect geometry:** on a straight and a curved synthetic river, normals are perpendicular (dot product ≈ 0) and left/right are correct looking downstream.
2. **Reach extraction:** a synthetic valley with a planted impact extent of known distance and height is recovered within one pixel; gap, nodata and censoring rules behave as specified.
3. **Exposure:** a belt with a known forest share returns that share.
4. **Statistics, planted effect:** simulate a table with a known `β1` and spatially correlated noise; Model A and B recover it within the bootstrap interval.
5. **Statistics, null:** with `β1 = 0`, the 95% interval covers 0 in roughly 95% of 200 simulations (block bootstrap coverage check).
6. **Confounding check:** simulate forest placed on steeper banks with no true effect; confirm the model with `bank_slope` returns null and the model without it does not (documents why the control is there).
7. **Date assertions:** exposure built from a post-event date raises an error.
8. **E-value:** matches a hand-computed example.

---

## 12. Build order

1. Scaffolding, config, tests 1–8 passing on synthetic data.
2. `fetch` for `melamchi2021`; user verifies centerline in QGIS.
3. `label` → user inspects label maps before continuing. **Stop here and show the maps.**
4. `transects`, `profiles`, `table` → user inspects the six example cross-sections. **Stop and show.**
5. `stats`, `placebo`, `validate`, `report` for Melamchi.
6. Repeat for `lende2025`, then `bhotekoshi2026`.
7. Summary report.

## 13. Things the user must supply or confirm

- Earth Engine project ID.
- Verified centerline per event, and start/end points of each reach (see 14.16).
- `exclusions.geojson` (dams, reservoirs).
- Validation extents (EMSR927 download; any ICIMOD maps; digitised polygons).
- Any post-event dates to exclude because of later floods.

---

## 14. Clarifications (these override anything above that conflicts)

### 14.1 Event-day scenes

Each event has `event_time_utc` in config. Filter on the scene acquisition timestamp:

- **Pre** scenes: acquired before 00:00 UTC on the event date.
- **Post** scenes: acquired later than `event_time_utc + 12 h`.
- If the event hour is unknown, post starts 00:00 UTC on the day after.

Values: `bhotekoshi2026` = 2026-08-26T02:52Z. `lende2025` = hour unknown, so post starts 2025-07-09T00:00Z. `melamchi2021` = two surges (2021-06-15 and 2021-07-31 to 08-01); pre ends before 2021-06-15, post starts 2021-08-03T00:00Z. The date assertions in rule 0.4 use these timestamps.

### 14.2 Placebo periods

The placebo pair must contain no known flood in the corridor. Fixed choices:
| Event | Placebo S2 "pre" | Placebo S2 "post" | Fake event date for S1 windows |
|---|---|---|---|
| `melamchi2021` | Oct–Dec 2018 | Oct–Dec 2019 | 2019-06-15 |
| `lende2025` | Oct–Dec 2022 | Oct–Dec 2023 | 2023-07-08 |
| `bhotekoshi2026` | Oct–Dec 2022 | Oct–Dec 2023 | 2023-08-26 |

Placebo S1 windows are the real event's windows shifted to the fake event date (same month-day offsets). Forest, active channel and controls for the placebo are built from the placebo "pre" period. The user must confirm no notable flood hit the corridor in the placebo interval; the report lists the interval so a reader can check. (2024 is avoided because of the late-September 2024 floods in central Nepal and the 2025 Lende flood.)

### 14.3 Hansen loss window

Condition = no Hansen loss with loss year in `[2021, Y]`, where `Y` = calendar year of the S2 pre composite. The canopy height layer is from 2020, so loss up to 2020 is already reflected in it.

- If `Y < 2021` (Melamchi: `Y = 2020`), the range is empty: skip the Hansen condition and log that it was skipped.
- If Hansen's newest year is earlier than `Y`, use loss through the newest available year, log the uncovered years, and rely on Dynamic World for them. This is not an error.

### 14.4 2026 data still arriving

Expected. Runs of `bhotekoshi2026` before the post window has closed are allowed but the report must carry a **PROVISIONAL** banner, show the clear-observation count map, and give no verdict (state "verdict withheld until the post-event window is complete"). The final run happens after 2026-12-15.

### 14.5 "Inside the mapped impact reach"

This is a **geographic** filter, never an outcome filter. It means stations with `river_km` between `reach_start_km` and `reach_end_km` in config. Rows with zero reach stay in the analysis. Never select rows by whether they show impact.

- `lende2025`: `reach_end_km = 35` (from the draft paper) until the user supplies a better-documented end point.
- Sensitivity run for every event: reach end moved 5 km upstream and 5 km downstream.

### 14.6 Where an impacted run starts

The run must begin within 20 m of `bank_top` (the same gap tolerance used later). If the first 20 m beyond `bank_top` contains no impact, then `reach_distance_m = 0`, `reach_height_m = 0`, `overtopped = 0`.
Also record `detached_reach_m`: distance to the furthest impact sample on the transect that passed the 2-D connectivity filter, regardless of gaps. Add one sensitivity run that uses `detached_reach_m` as the outcome.

### 14.7 Zero reach in the censored model

Do not add constants to zeros. Use a two-part approach:

1. The logistic model for `overtopped` (already specified) handles zero vs non-zero.
2. The Weibull AFT model is fitted **only on rows with reach > 0**, with right-censoring.
   Label part 2 in the report as a sensitivity analysis conditional on overtopping, not as a headline result. The primary estimates remain the OLS models in 8.2 and 8.3, which include zeros.

### 14.8 Overlap check vs models

- The 50% split is used only to define common support. Models always use `F` as a continuous share.
- If either group at the 50% split has fewer than 15% of rows, split at the median of `F` instead and log the value used.
- Trimming applies to **Model A only**.
- Model B is not trimmed (each station is its own match). Report the distribution of `|F_L − F_R|`. If fewer than 30 stations have `|F_L − F_R| ≥ 0.2`, mark Model B as **uninformative** for that event; the verdict can then be at most "inconclusive".
- Sensitivity run: Model B restricted to stations where both banks survived Model A trimming.

### 14.9 E-value for the yes/no outcome

Use the odds ratio for `F` going 0 → 1 from the logistic model. If the prevalence of `overtopped` in the analysis set exceeds 15%, convert with `RR ≈ sqrt(OR)`; otherwise `RR ≈ OR`. Then apply the same E-value formula.

### 14.10 "Reasonably narrow" (for the "no evidence" verdict)

Fixed thresholds:

- Continuous outcomes: the whole 95% bootstrap interval for the 0 → 100% forest effect, divided by `SD(y)`, lies within **±0.25**.
- `overtopped`: the whole 95% interval for the odds ratio lies within **[0.67, 1.50]**.
  If the interval includes 0 (or OR = 1) but is wider than this, the verdict is "inconclusive", not "no evidence".

### 14.11 Radar reference zone

A pixel is in the stable reference zone if all hold: inside the corridor; elevation at least 150 m above `z0` of the **nearest station** (FABDEM); DEM slope below 35°; not layover or shadow in that orbit; no optical impact; valid data. Require at least 5,000 reference pixels per orbit, otherwise skip that orbit and log it.

### 14.12 SAR agreement rule

Intended as written: at least half of the valid orbits (1 of 1, 1 of 2, 2 of 3, 2 of 4). Layover leaves many gorge pixels with a single usable orbit, so a stricter rule would discard too much. Add a sensitivity run, `sar_strict`, where **all** valid orbits must agree.

### 14.13 Validation for 2021 and 2025

No ready-made open extent is assumed. Procedure for every event:

1. The code draws a random sample of 40 stations, stratified into forest-fraction terciles, and writes `data/{event}/validation/sample_stations.geojson` **without** the forest values.
2. The user digitises the flood edge on both banks at those stations from independent high-resolution post-event imagery (for example Google Earth Pro historical imagery), and saves `validation_edges.geojson`. Any published extent maps the user obtains are added alongside.
3. `validate` compares digitised vs computed reach, overall and by forest tercile.
   If no validation data is supplied, `validate` reports "not validated" and the verdict for that event can be at most **inconclusive**.

### 14.14 2020 canopy height used for later events

Accepted. Clearing after 2020 is caught by Dynamic World (pre-event) and Hansen. The known weakness is the reverse: forest that grew past 5 m after 2020 is counted as non-forest. Add a sensitivity run with forest defined by Dynamic World alone (no canopy-height condition), and report the share of belt pixels where the two definitions disagree.

### 14.15 Uncertainty method (replaces 8.4 as the primary interval)

The 1 km pairs block bootstrap with a percentile interval in 8.4 gives intervals that are too narrow when there are few segments and the data are correlated along the river. In a simulation with a true effect of zero, 24 to 72 km of river and 20 rows per km, its nominal 95% interval covered zero only 89–93% of the time. The method below covered 94–96% in the same simulation.

**Blocks.** `bootstrap_block_km: 3` in config. `block_id = floor((river_km − reach_start_km) / 3)`; a final block shorter than 1.5 km is merged into the previous one. If an event has fewer than 8 blocks, its verdict can be at most "inconclusive". `segment_id` (1 km) is kept in the table for reporting only.

**Primary interval for the OLS models (Model A and Model B): wild cluster bootstrap-t with the null imposed.**
To test `β1 = b0`:

1. Fit the full model; compute `t_obs = (β1_hat − b0) / SE_cluster`, where `SE_cluster` is the cluster-robust (CR1) standard error with clusters = 3 km blocks.
2. Fit the restricted model with `β1` fixed at `b0` (regress `y − b0·F` on all other terms); keep its fitted values `f_r` and residuals `u_r`.
3. For each of `B = 1999` draws: draw one Rademacher weight (+1 or −1) per block, set `y* = b0·F + f_r + w·u_r`, refit the full model, compute `t* = (β1* − b0) / SE_cluster*`.
4. `p(b0)` = share of draws with `|t*| ≥ |t_obs|`.

Report `p(0)`. The 95% interval is the set of `b0` with `p(b0) > 0.05`, found by bisection outward from `β1_hat` on each side (tolerance 1% of `SE_cluster`).

- "Interval excludes 0" in 8.10 means `p(0) ≤ 0.05`.
- The thresholds in 14.10 are applied to this interval.
- The E-value confidence limit in 8.7 uses this interval.

**Logistic model (`overtopped`):** cluster-robust standard error with 3 km blocks and a t critical value with `(number of blocks − 1)` degrees of freedom. Label it in the report as approximate.

**Sensitivity, always reported beside the primary interval:** the original 1 km pairs block bootstrap percentile interval from 8.4, labelled "known to be too narrow with few segments".

**Test 5 is redefined and its settings are frozen.** Simulate with a true effect of zero at reach lengths of 24, 36 and 72 km, 20 rows per km, and along-river correlation of 0.8 between neighbouring stations in the exposure, a control and the noise. 400 simulations each, 399 bootstrap draws, each simulation with its own seed. The primary interval must cover zero in **92% to 98%** of simulations at every reach length. Do not change the number of segments, the correlation, or the band to make this test pass. If it fails, stop and report the coverage numbers to the user.

### 14.16 Reach definitions and confluences

The lengths in the table in 1.4 are approximate. The analysed reach is always the config values `reach_start_km` and `reach_end_km`, measured along the user-verified centerline.

| Event            | Start (`river_km = 0`)                                                             | End                                                            | Rivers in order                                       |
| ---------------- | ---------------------------------------------------------------------------------- | -------------------------------------------------------------- | ----------------------------------------------------- |
| `melamchi2021`   | User-chosen point on the upper Melamchi River, recorded in config with coordinates | 40 km downstream                                               | Melamchi River, then Indrawati below their confluence |
| `lende2025`      | Rasuwagadhi border crossing                                                        | 35 km downstream                                               | Bhotekoshi, then Trishuli                             |
| `bhotekoshi2026` | Rasuwagadhi border crossing                                                        | 72 km downstream, or Bidur if the user confirms that is nearer | Bhotekoshi, then Trishuli                             |

Rules:

1. **No open-ended reach.** "And beyond" is not analysed. Stations past `reach_end_km` are kept in the table but flagged outside the reach.
2. **`lende2025` and `bhotekoshi2026` share one centerline file and the same `river_km = 0`,** so a station has the same ID in both events. Each event still uses its own pre-event active channel, forest and controls.
3. **Above Rasuwagadhi is excluded** for both events: that stretch is the source and runout zone of the flow, not a river corridor with banks to compare.
4. **Confluences inside a reach.** The user lists them in config as `section_breaks_km` (the code proposes candidates where MERIT `upa` rises by more than 50% within 300 m, and the user confirms). Each station gets a categorical `river_section`. Model A includes `river_section` as a fixed effect in addition to the spline. Model B is unchanged.
5. **Sensitivity run for any event with a section break:** Model A and Model B on the first section only (for `melamchi2021`, the Melamchi River above the Indrawati confluence). The report states that other rivers joining below a confluence carried their own flood water in the same storm.
6. After the centerline is verified, the code prints the along-river distance from start to each named place in config (`landmarks`), so the user can check the approximate lengths above against the real line.

### 14.17 Radar baseline is the same season one year earlier

The original S1 pre windows sat just before each event, in a different season from the post windows (dry pre-monsoon versus monsoon). Irrigated rice terraces beside the river are dry in spring and flooded in the monsoon, which changes radar backscatter strongly. That would mark cropland banks as "impacted" and never forested banks, a bias toward finding a protective association. The reference-zone standardisation in 5.3 does not remove it, because it only corrects the corridor-wide shift.

- **Primary S1 pre window = the post window's calendar dates one year earlier** (table in 4.1, already updated). Same relative orbit and pass direction as the post scenes, as before.
- **Sensitivity label `sar_near`:** the previous near-event windows (`melamchi2021` 2021-04-01 to 2021-06-14; `lende2025` 2025-05-01 to 2025-07-07; `bhotekoshi2026` 2026-06-01 to 2026-08-25). Report it beside `sar_only`; it is never the primary.
- The report lists any known flood between the baseline window and the event for each corridor, because its footprint would be attributed to the event. The user supplies these in config (`known_floods_between`).
- If an orbit has fewer than 3 scenes in the baseline year, skip that orbit and log it (rule in 5.3).
- Placebo S1 windows follow the same rule: the fake post window, and the same dates one year before it.
- Report, per event, the share of non-forest belt pixels flagged by SAR but not by optical, split by Dynamic World class (crops versus other). A high crop share is a warning that seasonal cropland change is still leaking into the label.

### 14.18 Narrow channels and transect centring (overrides 5.2 step 3, 5.4 connectivity, 7.1)

The centerline QA for `melamchi2021` showed that the spectral rule for the pre-event active channel (`MNDWI_pre > 0` or `NDVI_pre < 0.20`) does not detect the river where it is narrower than about one 10 m pixel or overhung by canopy, which is most of the upper 12 km. It also showed that a centerline traced on post-flood imagery can sit a few pixels from the pre-flood channel in the wide lower reach. Neither may be allowed to break the label or the measurements.

**Channel seed.** Every pixel touched by the verified centerline is a `channel_seed` pixel.

- The connectivity filter in 5.4 keeps impact pixels 8-connected to (pre-event active channel ∪ `channel_seed`).
- The active-channel rule in 5.2 step 3 keeps spectral channel pixels whose component comes within 30 m of the centerline.
- `channel_seed` is used for connectivity only. A seed pixel is **not** treated as active channel when excluding pixels from the impact label.

**Transect centring (replaces the bank-top rule in 7.1).** For each station, sample the pre-event active channel along the whole transect line (both sides, as one profile with the station at 0).

1. Find runs of active-channel samples; merge runs separated by 20 m or less.
2. Choose the run that contains the station; if none, the run whose nearest edge is closest to the station, provided that edge is within `recentre_max_m` (default 100).
3. If a run is chosen: `channel_resolved = 1`; left and right `bank_top` are the two ends of the run; `channel_width_m` is its length; `centre_offset_m` is the signed distance from the station to the run's midpoint; `z0` is the median DEM elevation over the run.
4. If no run is chosen: `channel_resolved = 0`; `bank_top` is 5 m either side of the station; `channel_width_m = 10` (one pixel, an upper bound); `centre_offset_m` is missing; `z0` is the DEM elevation at the station.

**Models.** Model A adds the indicator `channel_resolved`. Sensitivity run: rows with `channel_resolved = 1` only.

**QA output (extend `qa_centerline` once transects exist).** Per river km: share of stations with `channel_resolved = 1` and the median absolute `centre_offset_m`. List every stretch of 200 m or more where `|centre_offset_m|` exceeds 50 m, and every stretch of 500 m or more with `channel_resolved = 0` where the pre-event composite shows a channel wider than 30 m (use `MNDWI_pre > 0` within 150 m of the line as the test). The user reviews only those stretches in QGIS. Do not move the centerline automatically.

**Report.** State the share of stations with an unresolved pre-event channel, and that reach distance at those stations is measured from an assumed bank position with an uncertainty of about one pixel.

### 14.19 Label corrections after the first visual review (melamchi2021)

The first label maps matched the visible post-flood debris band well, but showed four problems. These rules apply to every event.

**(a) Dark and shadowed pixels are "unknown", not "no impact".** In deep terrain shadow the optical indices cannot detect a deposit, and several gorge stretches read as zero impact although the flood passed through them.

- Optical label = nodata where `B8 < 0.08` in the pre **or** the post composite (`dark_b8_threshold` in config), except pixels already classed as pre-event active channel.
- Report the share of corridor pixels masked this way, per river km.

**(b) A combined "not impacted" needs a valid optical observation.**

- `impact_combined = 1` if optical = 1 or radar = 1.
- `impact_combined = 0` only if optical = 0 (valid).
- Otherwise nodata (optical nodata and radar not 1).
  `sar_only` keeps its own 0/1/nodata as before.

**(c) Radar label excludes the pre-event active channel.** Radar change inside the pre-event channel and its gravel bars is set to nodata in `sar_only`, `sar_strict`, `sar_near` and the combined label, as the optical label already does.

**(d) Layover/shadow buffer.** `layover_shadow_buffer_m: 30`.

**(e) Manual exclusions for hillslope failures and tributary fans.** Landslides and tributary debris fans that touch the river pass the connectivity filter and would be measured as lateral flood reach. They are a different process.

- Config list per event: `manual_exclusions: [{from_km, to_km, side, reason}]`, `side` is `left` or `right` looking downstream, `reason` is `hillslope_failure` or `tributary_fan`.
- The user fills this list from the post-event contact sheets **before** any outcome table or model has been produced for that event, and the date it was fixed is written to provenance.
- Affected station × side rows stay in the table with the reason code and are left out of the primary models.
- Sensitivity run: models with these rows included.
- Initial list for `melamchi2021` (west bank = right, east bank = left; the user confirms the km values against the sheets):
  - right, km 3.7–4.7, `hillslope_failure`
  - left, km 10.1–10.55, `tributary_fan`
  - left, km 13.6–13.85, `tributary_fan`
  - right, km 19.7–20.0, `tributary_fan`
  - left, km 22.8–23.3, `tributary_fan`
    (Revised once after the hatched sheets showed three bars sitting beside, not on, the fans. No outcome table existed at the time.)

**(f) Extra reporting and one robustness run.**

- Area of radar-only impact split by slope class (under 15°, 15–30°, over 30°) and by pre-event forest versus non-forest.
- Sensitivity run with `reach_height_m` and `reach_distance_m` winsorised at the 97.5th percentile.

### 14.20 A zero needs a visible river

In a few forested gorge stretches (for `melamchi2021`, about km 2.0–2.3 and 5.6–6.1) the post-event image shows no trace of the river at all: it is narrower than a pixel and under canopy or partial shadow, and the dark-pixel mask in 14.19(a) does not catch it. The flood passed through, so "no impact detected" there is not evidence of zero reach.

- For each station compute `river_visible_post = 1` if, within 30 m of the station along the transect (either side), at least one sample has `NDVI_post < 0.30` or is impacted in the label being used; otherwise 0.
- Where `river_visible_post = 0`, both rows of that station get missing outcomes and the reason code `river_not_visible`, and are left out of the primary models.
- Sensitivity run: those rows included with reach set to 0.
- Report the number of such stations, their river km, and their mean forest fraction beside the same figures for the stations kept, so a reader can see what kind of bank was lost.

### 14.21 Transect direction and neighbour crossings (overrides section 6 items 4 and 6, first bullet)

Testing for crossings over the full 1,500 m flagged 408 of 704 transects for `melamchi2021`, almost all for crossings far up the hillside that no measurement would ever reach. The crossing test cannot use the reach edge either, because excluding rows by how far the flood reached would select on the outcome.

**Direction.** The local river direction at a station is the unit vector between the points `normal_baseline_m` (default **250**) upstream and downstream along the verified line, clamped at the ends of the line. This replaces the 100 m baseline, which follows small kinks in a hand-traced line and makes neighbouring transects fan out. Chaikin smoothing is no longer needed for direction.

**Usable length instead of exclusion.**

- `usable_length_m` = distance from the station to the first point where the transect meets another transect on the same bank from any station within 1 km along the river; 1,500 if it meets none.
- Every measurement on the transect (bank top, impacted run, valley wall, `valid_fraction`, forest belt) uses only samples within `usable_length_m`.
- If the impacted run reaches the usable limit, set `censored = 1` with `censor_reason = neighbour_limit`.
- A transect is excluded as `intersects_neighbour` only if `usable_length_m < 100`.
- The forest belt needs the full belt width: for the 200 m belt sensitivity run, rows with `usable_length_m < 200` are left out of that run.
- Sensitivity run: rows with `usable_length_m ≥ 300` only.
- Report the distribution of `usable_length_m` and the number of rows censored at the neighbour limit.

**Confirmed implementation choices.** `river_km` is measured on the verified line, not the smoothed one. The tributary test uses the largest MERIT upstream area within 100 m of the station and ignores decreases. A second river crossing counts only if the two crossings are more than 100 m apart along the transect.

### 14.22 Station position sensitivity

Transect results can depend on where the stations happen to fall, not only on their spacing (Scherelis, Doering & Laube 2023, _Transactions in GIS_ 27:1503–1525, test transect location separately from spacing).

- Sensitivity run `station_offset_50`: all stations shifted 50 m downstream (half the default spacing), with transects, exclusions, usable length and all measurements rebuilt from the shifted stations. Manual exclusions and section breaks stay at their river km.
- Report the primary estimates from the shifted stations beside the primary run in the sensitivity table.

### 14.23 Decisions after the first outcome table (melamchi2021)

**(a) Missing Dynamic World coverage.** Where a belt pixel has no Dynamic World observation in the S2 pre window, the forest test uses ESA WorldCover for the pre-event year instead (class 10 tree cover) together with the same ETH canopy height ≥ 5 m and Hansen conditions. WorldCover 2020 (`ESA/WorldCover/v100`) is the map for `melamchi2021`.

- Record `forest_source` per row: `dynamic_world`, `worldcover`, or `mixed` (share of belt pixels from each).
- `belt_no_data` applies only where both sources are missing.
- Sensitivity run: rows with `forest_source = dynamic_world` only.
- The report states which river km use WorldCover and why (no Dynamic World image for Sentinel-2 tile T45RUM in Oct–Dec 2020).
- First, list Dynamic World image dates for tile T45RUM between 2020-09-01 and 2021-03-31, so the report can say whether the gap is specific to the window. Do not widen the window.

**(b) Valley-wall censoring.** `censor_reason = valley_wall` only if the reach edge is within 10 m (one sample) of `valley_wall`. A run that ends further from the wall is not censored. The walk still stops at the valley wall.

**(c) Bank top on the far side of the station.** Where the re-centred channel lies wholly to one side of the station, one bank top is on the opposite side of the station and `bank_top_from_station_m` is negative. This is expected. The outward walk always starts at that bank top and moves away from the channel. Add a test with a synthetic channel lying 50 m to one side of the station that checks both walks go outward from the correct ends of the channel run.

**(d) Confirmed choices.** Pixels with tree probability < 0.5 are non-forest even without canopy height. `valid_fraction` counts only samples beyond the bank top. The definitions of `wall_slope_deg`, `bank_elev_m`, `gradient`, `curvature` and `belt_clipped` reported in `transects_columns.csv` are accepted.

### 14.24 A zero needs a valid observation beside the bank (refines 14.6)

In the second set of example cross-sections (station 6, km 0.62, left bank) the first 50 m beyond the bank top are all no-data, there is no impacted sample, and the row was recorded as reach 0. No data beside the bank is not evidence of no flood, and it is most common in shaded, forested gorge, so counting it as zero would bias toward a protective finding.

Walking outward from `bank_top` (excluding the bank-top sample itself), apply these in order:

1. If any impacted sample lies within 20 m of the bank top, the run starts at the first impacted sample, as in 14.6. Non-impacted or no-data samples before it within those 20 m do not matter.
2. Otherwise, if at least one valid, non-impacted sample lies within 20 m, then `reach_distance_m = 0`, `reach_height_m = 0`, `overtopped = 0`.
3. Otherwise (every sample within 20 m is no-data): if the first valid sample lies within 30 m and is impacted, the run starts there; in all other cases the outcomes are missing, with reason code `near_bank_nodata`, and the row is left out of the primary models.

(Revised: an earlier wording let a non-impacted sample at 10 m override an impacted sample at 20 m. That was not intended; 14.6 still decides when impact is present within 20 m.)

- Apply this to every label variant separately (combined, optical-only, radar-only, `sar_strict`, `sar_near`).
- Sensitivity run: `near_bank_nodata` rows included with reach 0.
- Report the number of such rows, their river km, and their mean forest fraction beside the rows kept, as in 14.20.

  14.25 Independent validation edges (melamchi2021)

Imagery: Esri Wayback, traced by hand in QGIS. Sentinel-2 and all flood-map
layers are not used. A point is valid only if its image_date (capture date,
not release date) is between 2021-08-01 and 2022-05-31. Stations with no
photo in that window are skipped (no points for either bank); skipping is
for imagery reasons only.

Input: data/melamchi2021/validation/validation_edges.gpkg with fields
station_id, side, image_date, wayback_release, confidence (1 sure, 2 fairly
sure, 3 unsure), note. "flood stayed in channel" means the edge was placed at
the channel edge; it is a real observation, not a skip.

Conversion to validation_edges.geojson (EPSG:32645):

- Snap each point to its own station's transect. Refuse the run (write
  nothing) if any point is more than 20 m from that transect.
- Ignore the typed side. Derive side (looking downstream) from which side of
  the re-centred channel the point lies on. List every point where the
  derived side differs from the typed side.
- List every sampled station whose station point is more than 100 m from the
  re-centred channel.

Coverage: for all 40 sampled stations, checked or skipped, with the
forest-fraction group of each (counts only). No forest or reach values are
printed for any station, and the comparison is not run until the coverage
table has been reviewed.
