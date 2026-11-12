# Review guide

This guide is for reviewing the results without running the code. Every number below is
copied from the reports in `results/`; the reports are the source.

## Current deliverables

The 15-section `report.html` described in spec section 10 is **not built**. The current
deliverables are three reports per event:

| Report | Melamchi 2021 | Lende 2025 |
|---|---|---|
| Statistics and verdict | `results/stats/STATS_REPORT.md` | `results/lende2025/stats/STATS_REPORT.md` |
| Placebo (no-flood control) | `results/placebo/PLACEBO_REPORT.md` | `results/lende2025_placebo/PLACEBO_REPORT.md` |
| Validation against hand-digitised edges | `results/REPORT.md` | `results/lende2025/REPORT.md` |

## Research question

Do riverbank forests limit how far a Himalayan debris flood reaches across the bank? The
outcome is the reach of mapped flood impact along transects perpendicular to the river
(height above the channel, and distance). The exposure is the forest fraction of a 100 m
belt along the bank. Model A compares transects with each other and controls for valley and
channel shape, slope, upstream area, bend, curvature, bank elevation, side and position
along the river. Model B compares the two banks of the same transect. The verdict follows
the decision rule fixed in advance (spec 8.10). The wording is that forest is *associated
with* reach; nothing here shows that forest protects or causes anything.

## Melamchi 2021

**Verdict: Inconclusive.** Model A: 380 rows, 258 stations, 12 blocks. Model B: 239
stations. Forest fraction from 0 to 1, primary interval (wild cluster bootstrap-t):

- Model A reach height: −0.52 m [−1.70, 0.74], p(0) = 0.309
- Model A log1p reach distance: −0.31 [−0.69, 0.12], p(0) = 0.114
- Model B reach height: 2.14 m [−0.31, 5.16], p(0) = 0.078
- Model B log1p reach distance: 0.09 [−0.78, 1.10], p(0) = 0.825
- Overtopped: odds ratio 0.75 [0.29, 1.96] (approximate)

Model A and Model B disagree in sign, and 25 sensitivity estimates flip sign relative to the
primary run. **Placebo:** not null. Model A log1p reach distance −0.38 [−0.86, −0.00],
p(0) = 0.050; Model B log1p reach distance −0.39 [−0.80, −0.08], p(0) = 0.020. Under 8.9
the method is producing artefacts and the main result cannot be trusted. **Validation:** 52
digitised points at 26 stations, 35 in the main comparison; mean difference −13.3 m
(computed minus digitised), 80% within 20 m, 20% under-detected. Under-detected share High
26.7% (4 of 15) vs Low 10.0% (1 of 10), a gap of 16.7 points, which is flagged (above 10):
the main result may be biased toward finding protection.

## Lende 2025

**Verdict: Inconclusive.** Model A: 413 rows, 285 stations, 12 blocks. Model B: 166
stations. Forest fraction from 0 to 1, primary interval:

- Model A reach height: −3.51 m [−7.57, −0.31], p(0) = 0.029
- Model A log1p reach distance: −0.75 [−1.44, −0.29], p(0) < 0.001
- Model B reach height: −3.32 m [−7.34, 0.34], p(0) = 0.075
- Model B log1p reach distance: −0.85 [−2.13, 0.26], p(0) = 0.128
- Overtopped: odds ratio 0.17 [0.06, 0.51] (approximate)

Two sensitivity estimates flip sign relative to the primary run. **Placebo:** not null.
Model A reach height −1.81 m [−3.83, −0.33], p(0) = 0.011. Under 8.9 the method is
producing artefacts and the main result cannot be trusted. Both models point the same way
(negative), unlike Melamchi's, but the flood-free placebo shows a smaller association in the
same direction, so the real effect cannot be separated from the method's artefact. **Validation:** 61 digitised
points at 31 stations, 37 in the main comparison (27 stations, 6 blocks); mean difference
−12.5 m, 75.7% within 20 m, 18.9% under-detected. Under-detected share High 8.3% (1 of 12)
vs Low 25.0% (3 of 12), a gap of −16.7 points, not flagged. With about 12 points per
group this is weak evidence, and most traced banks were "flood stayed in channel", so it
tests mainly where the map says nothing happened. The main set spans fewer than 8
blocks, so the validation intervals are approximate and the rank correlation has none.

## Limitations

- **Timure flood.** The 30 July 2025 Timure/Bhote Koshi flash flood falls inside the Lende
  post window and affects km 0–5. Imagery cannot separate it from the GLOF; no exclusion was
  made for it.
- **km 0–5 in validation.** The Lende validation sample was drawn from km 5–10 and 25–35
  only: km 0–5 was left out because of the Timure flood, and km 10–25 is the gorge.
- **km 15 block.** The block at km 15.0 was kept as ambiguous (possible village-edge or road
  ground).
- **Valid-fraction filter.** On Lende, the filter (below 0.7) removed 96 of 559 rows, and
  removed non-forest banks more often than forested ones (share with no forest 37.5% dropped
  vs 22.2% kept), so the forest comparison rests on fewer non-forest banks. It also removed
  15 of the 24 Lende validation points left out of the main comparison.
- **Bank-side rule at validation stations 66 and 338.** The side of a digitised point comes
  from the re-centred channel. At both stations the rule moved a point to the other bank.
  Station 66 was removed (the visible channel was about 20 m off the channel centre the rule
  uses); station 338 was removed (headworks reach). Stations 263 (shadow), 267 (trees) and
  335 (no points) were also removed. All removals are recorded in `config.yaml`
  (`edge_corrections`) and the provenance.
- **Six "no visible change" review points.** Of the nine Lende points that differ from the
  computed reach by more than 20 m, three were reviewed as visible change and six as no
  visible change (69-right, 249-right, 251-right, 252-right, 253-right, 283-left). Stations
  249–253 lie within about 340 m of each other in one 3 km block. Without the six, 31
  points remain, mean difference −6.0 m. Four points outside the main set that differ by
  more than 20 m were not reviewed (52-right, 260-right, 272-right, 315-right).
- **Collinearity.** Bank elevation, upstream area and river km are strongly correlated with
  each other in both events (|r| up to 0.995), and their VIFs are high in the fitted design.
  Forest fraction is in no pair above 0.7 and its VIF is below 5 in both events
  (`results/diagnostics/`). Model A refitted without bank elevation, upstream area and
  valley_wall_found is in `stats_collinear_sensitivity.csv` for each event; it is
  descriptive and not part of the verdict.
- **Lende placebo period** (`flood_check` in `config.yaml`, user check of 2026-11-10):
  - A sudden rise of the Bhote Koshi was reported on 2023-07-12/13 (Radio Nepal); police
    warnings were issued, and a damaging flood was not confirmed.
  - No Sentinel-1 scene exists for 2023-07-12 to 2023-07-15 in this corridor, so a short
    surge could not be seen in radar.
  - Small isolated optical changes near the river (km 2.2–2.6, 3.3–3.8, 4.2–4.5, 7.3–7.8,
    14.2–14.6, 23.6–23.9, 24.3–24.5, 26.9–27.0, 28.3, 29.8–30.0, 30.0, 33.6–33.8,
    34.3–34.5) are consistent with bar, bank, slope or construction change and are not
    matched in radar.
  - In the gorge (km 21–25) the 2023 channel is darker and less bare than in 2022 (lower
    bare-soil index), the opposite of flood scour.
- **Lende reach end.** km 35 is the draft-paper value (14.5), not a documented end point.
- **Lende km 0–5 exclusions.** 20 of the 100 rows in km 0–5 were excluded as man-made
  (settlement or construction polygons, coded `reservoir_or_dam`).
- **Lende WorldCover stand-in.** WorldCover 2021 stands in for the missing 2024 map in the
  pre-event forest where Dynamic World has no observation: 10 of 463 analysis rows.
- **Lende reach height baseline.** Reach height is measured from the median channel
  elevation (z0), which was more than 10 m above the lowest channel point in 24 of 463
  analysis rows.
- The general limitations of the design (observational data, canopy concealment, a 30 m DEM
  in a narrow gorge, a pre-event DEM against a changed riverbed, unmeasured soil and
  disturbance history, forest correlating with less settlement, few events) apply to both
  events.

## Not yet done

- The final report (`report.html`, spec section 10) for Lende and Melamchi, and
  `summary_report.html`.
- bhotekoshi2026: configured, not run.

## File map

| Path | Contents |
|---|---|
| `results/stats/STATS_REPORT.md` | Melamchi statistics report and verdict |
| `results/stats/stats_results.csv`, `results/stats/stats_summary.json` | Melamchi: every run and outcome; summary and checks |
| `results/stats/stats_collinear_sensitivity.csv` | Melamchi: Model A without the collinear controls |
| `results/stats/fig_forest_plot_A_height.png`, `results/stats/fig_forest_plot_B_height.png`, `results/stats/fig_model_b_forest_difference.png`, `results/stats/fig_propensity.png` | Melamchi statistics figures |
| `results/placebo/PLACEBO_REPORT.md` | Melamchi placebo report |
| `results/REPORT.md`, `results/limitations.md` | Melamchi validation report and its limitations |
| `results/validation_per_point.csv`, `results/validation_summary.csv`, `results/validation_cutoff_sensitivity.csv`, `results/validation_detached_sensitivity.csv`, `results/validation_plot.png` | Melamchi validation tables and plot |
| `data/melamchi2021/validation/validation_edges.geojson`, `data/melamchi2021/validation/validation_review.csv` | Melamchi digitised points and review statuses |
| `results/placebo_check/`, `results/placebo_period_check.md` | Melamchi check that the placebo period had no notable flood |
| `results/lende2025/stats/STATS_REPORT.md` | Lende statistics report and verdict |
| `results/lende2025/stats/stats_results.csv`, `results/lende2025/stats/stats_sign_flips.csv`, `results/lende2025/stats/stats_summary.json` | Lende: every run and outcome; sign flips; summary and checks |
| `results/lende2025/stats/stats_collinear_sensitivity.csv` | Lende: Model A without the collinear controls |
| `results/lende2025/stats/fig_forest_plot_A_height.png`, `results/lende2025/stats/fig_forest_plot_B_height.png`, `results/lende2025/stats/fig_model_b_forest_difference.png`, `results/lende2025/stats/fig_propensity.png` | Lende statistics figures |
| `results/lende2025/REPORT.md`, `results/lende2025/limitations.md` | Lende validation report and its limitations |
| `results/lende2025/validation_per_point.csv`, `results/lende2025/validation_summary.csv`, `results/lende2025/validation_cutoff_sensitivity.csv`, `results/lende2025/validation_detached_sensitivity.csv`, `results/lende2025/validation_plot.png` | Lende validation tables and plot |
| `results/lende2025/validation/sample_stations.geojson` | Lende validation sample (50 stations; ids, river km, location) |
| `results/lende2025/validation/validation_edges.csv`, `results/lende2025/validation/validation_edges.geojson` | Lende digitised points (as entered; snapped) |
| `results/lende2025/validation/validation_review.csv` | Lende second-pass review statuses |
| `results/lende2025/validation/sample_coverage.csv` | Lende checked or skipped per sample station (no forest values) |
| `results/lende2025/table/transects.csv`, `results/lende2025/table/transects_columns.csv` | Lende analysis table and its column dictionary |
| `results/lende2025_placebo/PLACEBO_REPORT.md` | Lende placebo report; its figures are in the same folder |
| `results/placebo_check_lende/` | Lende check that the placebo period had no notable flood |
| `results/diagnostics/README.md` | Collinearity diagnostic (correlations and VIFs, both events); tables and script in the same folder |
| `results/melamchi2021_checksums.sha1` | SHA-1 checksums of the committed Melamchi results |
| `config.yaml` | Every event setting and user decision, with dates |
