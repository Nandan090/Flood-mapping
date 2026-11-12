# Statistics: lende2025

Reproduce with `python -m floodforest.cli stats --event lende2025`. Seed 42.

## Verdict (8.10): Inconclusive

Wording follows the spec: forest is *associated with* reach; nothing here shows that forest protects or causes anything.

**The placebo shows a forest association (Model A reach height (m): -1.81 [-3.83, -0.33], p(0) = 0.011). Under 8.9 the method is producing artefacts and the main result cannot be trusted.**

**2 sensitivity estimates flip sign relative to the primary run (8.10); see the sign-flip table.**

- placebo not null (8.9): the main result cannot be trusted
- sign flips across sensitivity runs (8.10)
- not met: Model B reach height: beta < 0, interval excludes 0
- not met: placebo null

| Check | Result |
|---|---|
| a_height_neg_sig | True |
| b_height_neg_sig | False |
| dist_same_sign | True |
| labels_same_sign | True |
| dems_same_sign | True |
| placebo_null | False |
| placebo_effects | Model A reach height (m): -1.81 [-3.83, -0.33], p(0) = 0.011 |
| validation_ok | True |
| validation_note | under-detected High 8.3% vs Low 25.0% (gap -16.7 points, flag above 10) |
| validated | True |
| enough_blocks | True |
| overlap_ok | True |
| model_b_informative | True |
| no_sign_flips | False |
| sign_flips | 2 (run x model x outcome); see the sign-flip table |
| all_include_null | False |
| all_narrow | False |

## Known floods around the event (14.17)

Supplied by the user in config.yaml. The effect column follows from the event's own image windows.

**Between the Sentinel-1 baseline and the event**

| Dates | Place | Type | Note | Effect on the labels | Handling |
|---|---|---|---|---|---|
| 2024-08-25 | Kerung rivulet near Rasuwagadhi | flash_flood | local tributary | radar label only: after the Sentinel-1 baseline (2024-07-09 to 2024-08-15) and before the Sentinel-2 pre composite (2024-10-01 to 2024-12-15), so its scars are part of the optical pre-event state |  |
| 2024-09-27 to 2024-09-29 | Trishuli basin | river_flood | basin-wide; not shown to inundate the whole reach | radar label only: after the Sentinel-1 baseline (2024-07-09 to 2024-08-15) and before the Sentinel-2 pre composite (2024-10-01 to 2024-12-15), so its scars are part of the optical pre-event state |  |

**After the event, inside the post window**

| Dates | Place | Type | Note | Effect on the labels | Handling |
|---|---|---|---|---|---|
| 2025-07-30 | Timure and Bhotekoshi | flash_flood | road erosion and downstream flooding | inside or before the end of: Sentinel-1 post (2025-07-09 to 2025-08-15); Sentinel-2 post (2025-10-01 to 2025-12-15), whose scenes all come after it | none; see limitations |

The floods between the baseline and the event (2024-08-25, 2024-09-27 to 2024-09-29) affect only the radar label: their footprint can appear as radar change and be attributed to the event. The optical label is not affected, because its pre-event composite already shows their scars.
The 2025-07-30 flood at Timure and Bhotekoshi falls inside the post window, so the post-event labels can contain its footprint as well as the event's. Handling: none; see limitations.

## Models

Model A (8.2, 14.16, 14.18): `y ~ forest_frac + np.log(channel_width_m) + gradient + np.log(valley_width_m) + wall_slope_deg + bank_slope_deg + log10_upa + outer_bend + curvature + bank_elev_m + C(side) + cr(river_km, df=5) + channel_resolved`. OLS for reach height and log1p(reach distance); logistic for overtopped. `valley_wall_found` is in the formula before trimming but left out of the fit because it is constant in the trimmed Model A rows.

Model B (8.3): `(y_L - y_R) ~ d_forest + d_valley_width_m + d_wall_slope_deg + d_valley_wall_found + d_bank_slope_deg + d_outer_bend`, stations with both banks in the analysis set; not trimmed (14.8).

Wall slope is 0 where no valley wall was found, with valley_wall_found in both models (user decision); a complete-case run is in the sensitivity table. Model B pairs are built from Model B's own variables, so a bank missing only a Model A control (for example curvature at the end of the line) still pairs.

### Methods notes

- **Gap tolerance (7.2).** Samples are 10 m apart and a run of N samples spans (N - 1) x 10 m. A non-impact gap of up to 3 consecutive samples (span <= 20 m) stays inside the run; a nodata gap of up to 4 samples (span <= 30 m) does not end it, and a longer one ends it with censored = 1. A nodata sample neither resets nor adds to the non-impact count. The same convention sets the valley-wall run (50 m = 6 samples) and the channel-run merge (20 m, 14.18).
- **Channel width (14.18).** channel_width_m is the length of the chosen channel run counted in whole pixels: last sample minus first sample plus one pixel (10 m). Where no channel run is found it is 10 m (one pixel, an upper bound).
- **Censored model (14.7).** The Weibull AFT uses lifelines' robust (sandwich) standard errors, not clustered on blocks; it is a sensitivity conditional on reach > 0.
- **Percentile intervals (8.4, 14.15).** The 1 km pairs percentile interval is computed for the primary run and the placebo's primary variant only; sensitivity runs report the 14.15 interval alone.
- **E-values (8.7, 14.9).** Each model uses its own SD(y): Model A the SD of y in its rows, Model B the SD of y_L - y_R in its pairs. The 14.9 prevalence of overtopped is taken over the run's full analysis set. 14.10 narrowness uses Model A's SD(y) for both models, as before.
- **WorldCover stand-in (14.23 a, user decision 2026-11-03).** WorldCover 2021 stands in for 2024, which has no map, where Dynamic World has no observation in the S2 pre window. 10 of 463 analysis-set rows use it (forest_source worldcover or mixed); the 'forest_source dynamic_world only' run leaves them out.

Primary interval (14.15): wild cluster bootstrap-t, null imposed, CR1, Rademacher, 1999 draws, 3 km blocks; p(0) reported; interval = b0 with p(b0) > 0.05. Logistic: cluster-robust SE with a t critical value (blocks - 1 df), approximate. Beside it: the 1 km pairs block bootstrap percentile interval (2000 reps), known to be too narrow with few segments.

## Common support (8.1, 14.8)

Split: F >= 0.5. Kept 413 of 461 rows (dropped 48), 35 segments; p bounds (0.052, 0.903). Too thin to interpret: False.
Rows without every Model A value (curvature at the line end): 2; these rows can still enter Model B.

## Primary results

Model A: 413 rows, 285 stations, 35 segments, 12 blocks. Model B: 166 stations, 35 segments, 12 blocks; 125 with |F_L - F_R| >= 0.2 (minimum 30); informative: True.

| Model | Outcome | beta1 for F 0 -> 1 [95% interval], p(0) | Percentile interval (too narrow with few segments) | Per IQR of F | E-value (point / CI) | Reasonably narrow (14.10) |
|---|---|---|---|---|---|---|
| A | reach height (m) | -3.51 [-7.57, -0.31], p(0) = 0.029 | [-6.60, -0.70] | -2.44 [-5.27, -0.21] (IQR 0.70) | 2.34 / 1.23 | False |
| A | log1p reach distance | -0.75 [-1.44, -0.29], p(0) = 0.000 | [-1.35, -0.31] | -0.52 [-1.00, -0.20] (IQR 0.70) | 2.61 / 1.70 | False |
| B | reach height (m) | -3.32 [-7.34, 0.34], p(0) = 0.075 | [-7.42, 0.29] | -2.31 [-5.11, 0.24] (IQR 0.70) | 1.90 / 1.00 | False |
| B | log1p reach distance | -0.85 [-2.13, 0.26], p(0) = 0.128 | [-1.77, -0.07] | -0.59 [-1.48, 0.18] (IQR 0.70) | 2.29 / 1.00 | False |
| A | overtopped | OR 0.17 [0.06, 0.51] (approximate, t with 11 df) | | | 4.32 / 2.16 (prevalence 18%) | False |

log1p distance: exp(beta1) - 1 is the proportional change in (1 + distance).

E-value (8.7): an E-value of 2 means an unmeasured factor would need to be associated with both forest and flood reach by a risk ratio of 2 each to erase the result. The CI E-value is 1 when the interval includes the null.

### Censored model, conditional on reach > 0 (14.7, sensitivity)

- y_dist: time ratio for F 0 -> 1 0.48 [0.31, 0.75], p 0.001; 163 rows, 57 censored; 0 rows with a zero value left out. Weibull AFT, right-censored, robust (not cluster) SE; conditional on reach > 0.
- y_height: time ratio for F 0 -> 1 0.54 [0.25, 1.15], p 0.109; 129 rows, 56 censored; 34 rows with a zero value left out. Weibull AFT, right-censored, robust (not cluster) SE; conditional on reach > 0.

## Sensitivity runs (8.8 and section 14)

Model A reach height and Model B reach height, 14.15 interval. All runs, all outcomes: `stats_results.csv`. Forest plots: `fig_forest_plot_A_height.png`, `fig_forest_plot_B_height.png`.

| Run | Group | Note | A rows | A height beta1 [interval], p(0) | A log1p distance beta1 [interval], p(0) | B stations / segments | B height beta1 [interval], p(0) |
|---|---|---|---:|---|---|---:|---|
| primary | primary | combined label, FABDEM, 100 m belt, trimmed A | 413 | -3.51 [-7.57, -0.31], p(0) = 0.029 | -0.75 [-1.44, -0.29], p(0) = 0.000 | 166 / 35 | -3.32 [-7.34, 0.34], p(0) = 0.075 |
| label optical_only | label | 8.8 | 411 | -2.82 [-6.11, -0.17], p(0) = 0.035 | -0.67 [-1.34, -0.15], p(0) = 0.012 | 155 / 34 | -2.88 [-6.81, 0.92], p(0) = 0.131 |
| label sar_only | label | 8.8 | 485 | -0.77 [-2.20, 0.45], p(0) = 0.267 | -0.25 [-0.60, 0.04], p(0) = 0.099 | 232 / 34 | -0.85 [-2.11, 0.16], p(0) = 0.125 |
| label sar_strict | label | 14.12 | 485 | -0.17 [-0.46, 0.02], p(0) = 0.107 | -0.06 [-0.18, 0.03], p(0) = 0.192 | 232 / 34 | -0.32 [-0.93, 0.03], p(0) = 0.118 |
| label sar_near | label | 14.17 | 524 | -0.34 [-2.33, 1.61], p(0) = 0.712 | -0.05 [-0.58, 0.50], p(0) = 0.827 | 235 / 34 | 0.55 [-0.93, 2.16], p(0) = 0.483 |
| DEM Copernicus GLO-30 | DEM | 8.8; outcomes and DEM controls from GLO-30 | 387 | -3.20 [-7.86, 0.19], p(0) = 0.065 | -0.65 [-1.43, -0.12], p(0) = 0.013 | 149 / 35 | -2.95 [-7.04, 0.73], p(0) = 0.140 |
| belt 50 m | belt width | 8.8 | 446 | -3.63 [-7.41, -0.51], p(0) = 0.019 | -1.09 [-1.86, -0.47], p(0) = 0.003 | 166 / 35 | -3.34 [-7.96, 0.56], p(0) = 0.093 |
| belt 200 m | belt width | 8.8; rows with usable_length_m >= 200 (14.21) | 416 | -3.58 [-6.32, -0.98], p(0) = 0.011 | -0.66 [-1.20, -0.23], p(0) = 0.008 | 163 / 35 | -4.01 [-6.57, -0.47], p(0) = 0.035 |
| forest WorldCover | forest definition | 8.8: WorldCover class 10 share | 435 | -4.64 [-8.72, -1.10], p(0) = 0.012 | -0.73 [-1.42, -0.08], p(0) = 0.022 | 166 / 35 | -4.82 [-10.64, -0.01], p(0) = 0.049 |
| forest Dynamic World only | forest definition | 14.14: no canopy-height condition; the two rules disagree on 0.1% of belt pixels, so this run nearly equals the dynamic_world-only run | 408 | -3.35 [-7.27, -0.19], p(0) = 0.036 | -0.72 [-1.34, -0.29], p(0) = 0.002 | 165 / 35 | -3.39 [-7.40, 0.20], p(0) = 0.073 |
| forest_source dynamic_world only | forest definition | 14.23 a | 410 | -3.52 [-7.76, -0.21], p(0) = 0.032 | -0.73 [-1.41, -0.25], p(0) = 0.001 | 162 / 35 | -3.37 [-7.44, 0.27], p(0) = 0.072 |
| steep: bank-top slope <= 45 deg | steep terrain | 8.8 (user decision): pre-event terrain, 10 x tan(bank-top slope) <= 10 m | 415 | -3.68 [-7.54, -0.75], p(0) = 0.016 | -0.83 [-1.55, -0.36], p(0) = 0.002 | 156 / 35 | -4.94 [-9.75, -0.07], p(0) = 0.042 |
| steep: unc <= 10 m plus all zero-reach rows | steep terrain | 8.8 variant (user decision): keeps every zero-reach row | 397 | -3.32 [-6.63, -0.84], p(0) = 0.003 | -0.78 [-1.46, -0.35], p(0) = 0.000 | 163 / 35 | -3.87 [-7.98, 0.51], p(0) = 0.069 |
| steep: literal unc <= 10 m | steep terrain | selects rows by outcome (conflicts with 14.5); not a headline result | 146 | -4.40 [-10.26, 1.54], p(0) = 0.144 | -0.42 [-0.73, -0.14], p(0) = 0.013 | 31 / 18 | -6.97 [-15.35, 5.69], p(0) = 0.105 |
| terraces: bank slope <= 15 deg | terrain | 8.8 | 15 | 0.28 [NA, NA], p(0) = 0.000 | -175.70 [NA, NA], p(0) = 0.000 | 0 / 0 | not estimable (ValueError: zero-size array to reduction operation maximum which has no identity) |
| reach end 30.0 km (5 km upstream) | reach | 14.5 | 386 | -3.10 [-7.83, 0.32], p(0) = 0.074 | -0.75 [-1.61, -0.28], p(0) = 0.000 | 152 / 30 | -2.99 [-8.07, 0.64], p(0) = 0.139 |
| detached reach | outcome | 14.6: detached_reach_m as the distance outcome | 413 | -3.51 [-7.57, -0.31], p(0) = 0.029 | -0.78 [-1.59, -0.22], p(0) = 0.002 | 166 / 35 | -3.32 [-7.34, 0.34], p(0) = 0.075 |
| winsorised outcomes | outcome | 14.19 f: 97.5th percentile caps | 413 | -2.45 [-5.48, -0.14], p(0) = 0.034 | -0.74 [-1.41, -0.28], p(0) = 0.000 | 166 / 35 | -3.45 [-7.54, 0.64], p(0) = 0.107 |
| Model A untrimmed | overlap | 14.8: no common-support trimming | 461 | -3.67 [-7.65, -0.53], p(0) = 0.024 | -0.81 [-1.46, -0.32], p(0) = 0.000 | 166 / 35 | -3.32 [-7.34, 0.34], p(0) = 0.075 |
| Model B: both banks survived trimming | overlap | 14.8 | 413 | -3.51 [-7.57, -0.31], p(0) = 0.029 | -0.75 [-1.44, -0.29], p(0) = 0.000 | 128 / 33 | -4.25 [-8.36, -0.60], p(0) = 0.026 |
| first section only | river section | 14.16: above the first section break | 413 | -3.51 [-7.57, -0.31], p(0) = 0.029 | -0.75 [-1.44, -0.29], p(0) = 0.000 | 166 / 35 | -3.32 [-7.34, 0.34], p(0) = 0.075 |
| channel_resolved = 1 | channel | 14.18 | 402 | -3.57 [-7.62, -0.40], p(0) = 0.027 | -0.79 [-1.44, -0.35], p(0) = 0.000 | 161 / 35 | -3.33 [-7.42, 0.39], p(0) = 0.074 |
| manual exclusions included | exclusions | 14.19 e | 417 | -3.77 [-7.97, -0.46], p(0) = 0.024 | -0.78 [-1.44, -0.32], p(0) = 0.001 | 171 / 35 | -2.69 [-6.82, 1.00], p(0) = 0.167 |
| river_not_visible rows as reach 0 | zeros | 14.20: 18 rows added | 427 | -3.32 [-7.19, -0.32], p(0) = 0.027 | -0.69 [-1.35, -0.24], p(0) = 0.000 | 170 / 35 | -3.20 [-7.31, 0.53], p(0) = 0.079 |
| near_bank_nodata rows as reach 0 | zeros | 14.24: 10 rows added | 422 | -3.43 [-7.41, -0.34], p(0) = 0.025 | -0.70 [-1.30, -0.29], p(0) = 0.000 | 172 / 35 | -3.03 [-6.86, 0.41], p(0) = 0.080 |
| usable length >= 300 m | neighbours | 14.21 | 417 | -3.18 [-7.15, -0.06], p(0) = 0.047 | -0.72 [-1.41, -0.23], p(0) = 0.001 | 160 / 33 | -2.83 [-7.09, 0.94], p(0) = 0.175 |
| wall slope: complete cases | wall slope | user decision: rows with a valley wall only | 413 | -3.51 [-7.57, -0.31], p(0) = 0.029 | -0.75 [-1.44, -0.29], p(0) = 0.000 | 164 / 35 | -3.34 [-7.37, 0.35], p(0) = 0.073 |
| stations with a visible river in the placebo | placebo comparability | (b): real event restricted to stations whose river was visible in the placebo | 415 | -3.38 [-7.52, -0.13], p(0) = 0.039 | -0.73 [-1.42, -0.27], p(0) = 0.002 | 161 / 35 | -3.35 [-7.47, 0.33], p(0) = 0.074 |
| station spacing 50 m (8.8) | rebuilt | rebuilt from label to table (lende2025__spacing50) | 882 | -2.52 [-4.51, -0.90], p(0) = 0.000 | -0.67 [-1.05, -0.34], p(0) = 0.000 | 327 / 35 | -1.71 [-5.15, 1.49], p(0) = 0.290 |
| stations shifted 50 m downstream (14.22) | rebuilt | rebuilt from label to table (lende2025__offset50) | 441 | -2.04 [-3.85, -0.19], p(0) = 0.038 | -0.73 [-1.11, -0.36], p(0) = 0.002 | 166 / 35 | -0.42 [-5.13, 3.55], p(0) = 0.827 |
| t_bsi x 0.8 (8.8) | rebuilt | rebuilt from label to table (lende2025__tbsi_minus20) | 413 | -3.50 [-7.56, -0.29], p(0) = 0.031 | -0.75 [-1.42, -0.29], p(0) = 0.001 | 166 / 35 | -3.35 [-7.36, 0.29], p(0) = 0.070 |
| t_bsi x 1.2 (8.8) | rebuilt | rebuilt from label to table (lende2025__tbsi_plus20) | 413 | -3.32 [-7.17, -0.35], p(0) = 0.031 | -0.63 [-1.36, -0.11], p(0) = 0.005 | 166 / 35 | -3.05 [-6.73, 0.31], p(0) = 0.076 |
| SAR z threshold 2.4 (8.8) | rebuilt | rebuilt from label to table (lende2025__sarz_minus20) | 432 | -4.32 [-8.90, -0.45], p(0) = 0.030 | -0.76 [-1.41, -0.27], p(0) = 0.001 | 178 / 35 | -2.78 [-7.56, 2.11], p(0) = 0.264 |
| SAR z threshold 3.6 (8.8) | rebuilt | rebuilt from label to table (lende2025__sarz_plus20) | 415 | -3.20 [-6.90, -0.30], p(0) = 0.033 | -0.78 [-1.46, -0.31], p(0) = 0.002 | 159 / 34 | -2.89 [-7.33, 1.38], p(0) = 0.199 |
| reach end +5 km downstream | reach | 14.5 | | not run: the verified centerline ends at reach_end_km (35 km) | | | |

### Sign flips across sensitivity runs (8.10)

A run flips when its beta1 has the opposite sign to the primary run's for the same model and outcome (Model A and B reach height and log1p distance, Model A overtopped). Any flip makes the verdict inconclusive. All runs are compared, including the literal steep filter. Full list: `stats_sign_flips.csv`.

| Run | Model | Outcome | Primary beta1 | Run beta1 | Run p(0) |
|---|---|---|---:|---:|---:|
| label sar_near | B | reach height (m) | -3.32 | 0.55 | 0.483 |
| terraces: bank slope <= 15 deg | A | reach height (m) | -3.51 | 0.28 | 0.000 |

(d) The real event's 14.20 sensitivity (river_not_visible rows included as reach 0) is reported: yes, row above.

## Placebo through the same models (8.9)

Same Model A and Model B code and the same 14.15 interval as the real event.

| Variant | Note | A rows | A reach height | A log1p distance | A overtopped | B stations / segments | B reach height | B log1p distance |
|---|---|---:|---|---|---|---:|---|---|
| placebo: kept rows (primary) |  | 413 | -1.81 [-3.83, -0.33], p(0) = 0.011 | -0.22 [-0.66, 0.19], p(0) = 0.256 | OR 0.88 [0.17, 4.59] (approximate, t with 11 df) | 152 / 34 | -1.10 [-3.89, 1.99], p(0) = 0.413 | -0.03 [-0.56, 0.57], p(0) = 0.896 |
| placebo: + 14.20 assumed-zero rows | pushes toward null by construction | 418 | -1.76 [-3.58, -0.37], p(0) = 0.009 | -0.22 [-0.63, 0.16], p(0) = 0.229 | OR 0.85 [0.00, 2087.93] (approximate, t with 11 df) | 163 / 34 | -0.97 [-3.67, 2.02], p(0) = 0.449 | -0.01 [-0.52, 0.54], p(0) = 0.947 |
| placebo: 2024 forest |  | 404 | -1.63 [-3.69, 0.02], p(0) = 0.052 | -0.22 [-0.62, 0.12], p(0) = 0.193 | OR 0.68 [0.11, 4.24] (approximate, t with 11 df) | 152 / 34 | -0.24 [-2.65, 2.98], p(0) = 0.833 | -0.25 [-0.87, 0.42], p(0) = 0.398 |

(c) The placebo drops 11% of High-forest rows for river_not_visible (25 of 233), against 1% of Low-forest rows (2 of 233) (terciles of the 2022 forest over all in-reach rows). The placebo therefore cannot speak for the densest forest.

The 14.20 assumed-zero variant adds those rows with reach 0; it pushes toward null by construction, so it is not evidence that the placebo is null.

## Limitations

- The 30 July 2025 Timure/Bhote Koshi flash flood falls inside the post window and affects km 0-5; imagery cannot separate it from the GLOF; no exclusion was made for it.
- The block at km 15.0 was kept as ambiguous (possible village-edge or road ground).
- Exclusion areas a, b, c, e_NW, h (settlement/construction) were identified by sheet review; b and h also checked in Google Earth; the others were not. Rows inside these polygons are coded reservoir_or_dam in the outcome table although they are settlement/construction areas, and, like the dam and bridge polygons, they are not added back in the exclusions-included sensitivity run.
- The reach end, km 35, is the draft-paper value (14.5), not a documented end point; it was kept by the user.
- The valid-fraction filter (below 0.7) removed 96 of 559 rows and removed non-forest banks more often than forested ones (share with no forest 37.5% dropped vs 22.2% kept; mean forest fraction 0.32 vs 0.45), so the forest comparison rests on fewer non-forest banks and on those with less pre-event bare or built ground.

## Files

- `stats_results.csv`, `stats_summary.json`, `stats_sign_flips.csv`
- `fig_propensity.png`, `fig_model_b_forest_difference.png`, `fig_forest_plot_A_height.png`, `fig_forest_plot_B_height.png`

