# Statistics: melamchi2021

Reproduce with `python -m floodforest.cli stats --event melamchi2021`. Seed 42.

## Verdict (8.10): Inconclusive

Wording follows the spec: forest is *associated with* reach; nothing here shows that forest protects or causes anything.

**The placebo shows a forest association (Model A log1p reach distance: -0.38 [-0.86, -0.00], p(0) = 0.050; Model B log1p reach distance: -0.39 [-0.80, -0.08], p(0) = 0.020). Under 8.9 the method is producing artefacts and the main result cannot be trusted.**

- placebo not null (8.9): the main result cannot be trusted
- not met: Model A reach height: beta < 0, interval excludes 0
- not met: Model B reach height: beta < 0, interval excludes 0
- not met: reach distance: same sign in Model A and B
- not met: same sign under optical-only and radar-only labels
- not met: same sign under both DEMs
- not met: placebo null
- not met: validation omission not materially higher on forested banks
- intervals include the null but are not reasonably narrow (14.10)

| Check | Result |
|---|---|
| a_height_neg_sig | False |
| b_height_neg_sig | False |
| dist_same_sign | False |
| labels_same_sign | False |
| dems_same_sign | False |
| placebo_null | False |
| placebo_effects | Model A log1p reach distance: -0.38 [-0.86, -0.00], p(0) = 0.050; Model B log1p reach distance: -0.39 [-0.80, -0.08], p(0) = 0.020 |
| validation_ok | False |
| validation_note | under-detected High 26.7% vs Low 10.0% (gap 16.7 points, flag above 10) |
| validated | True |
| enough_blocks | True |
| overlap_ok | True |
| model_b_informative | True |
| all_include_null | True |
| all_narrow | False |

## Models

Model A (8.2, 14.16, 14.18): `y ~ forest_frac + np.log(channel_width_m) + gradient + np.log(valley_width_m) + wall_slope_deg + bank_slope_deg + log10_upa + outer_bend + curvature + bank_elev_m + C(side) + cr(river_km, df=5) + valley_wall_found + channel_resolved + C(river_section)`. OLS for reach height and log1p(reach distance); logistic for overtopped.

Model B (8.3): `(y_L - y_R) ~ d_forest + d_valley_width_m + d_wall_slope_deg + d_valley_wall_found + d_bank_slope_deg + d_outer_bend`, stations with both banks in the analysis set; not trimmed (14.8).

Wall slope is 0 where no valley wall was found, with valley_wall_found in both models (user decision); a complete-case run is in the sensitivity table.

Primary interval (14.15): wild cluster bootstrap-t, null imposed, CR1, Rademacher, 1999 draws, 3 km blocks; p(0) reported; interval = b0 with p(b0) > 0.05. Logistic: cluster-robust SE with a t critical value (blocks - 1 df), approximate. Beside it: the 1 km pairs block bootstrap percentile interval (2000 reps), known to be too narrow with few segments.

## Common support (8.1, 14.8)

Split: F >= 0.5. Kept 380 of 542 rows (dropped 162), 36 segments; p bounds (0.036, 0.996). Too thin to interpret: False.
Rows dropped for a missing model value (curvature at the line end): 1.

## Primary results

Model A: 380 rows, 258 stations, 36 segments, 12 blocks. Model B: 239 stations, 12 blocks; 113 with |F_L - F_R| >= 0.2 (minimum 30); informative: True.

| Model | Outcome | beta1 for F 0 -> 1 [95% interval], p(0) | Percentile interval (too narrow with few segments) | Per IQR of F | E-value (point / CI) | Reasonably narrow (14.10) |
|---|---|---|---|---|---|---|
| A | reach height (m) | -0.52 [-1.70, 0.74], p(0) = 0.309 | [-2.64, 1.66] | -0.41 [-1.36, 0.60] (IQR 0.80) | 1.33 / 1.00 | True |
| A | log1p reach distance | -0.31 [-0.69, 0.12], p(0) = 0.114 | [-0.87, 0.36] | -0.25 [-0.56, 0.09] (IQR 0.80) | 1.67 / 1.00 | False |
| B | reach height (m) | 2.14 [-0.31, 5.16], p(0) = 0.078 | [-0.37, 4.70] | 1.72 [-0.25, 4.14] (IQR 0.80) | 1.93 / 1.00 | False |
| B | log1p reach distance | 0.09 [-0.78, 1.10], p(0) = 0.825 | [-0.54, 0.74] | 0.07 [-0.62, 0.88] (IQR 0.80) | 1.28 / 1.00 | False |
| A | overtopped | OR 0.75 [0.29, 1.96] (approximate, t with 11 df) | | | 1.57 / 1.00 (prevalence 64%) | False |

log1p distance: exp(beta1) - 1 is the proportional change in (1 + distance).

### Censored model, conditional on reach > 0 (14.7, sensitivity)

- y_dist: time ratio for F 0 -> 1 0.98 [0.63, 1.51], p 0.918; 302 rows, 81 censored; 0 rows with a zero value left out. Weibull AFT, right-censored, robust (not cluster) SE; conditional on reach > 0.
- y_height: time ratio for F 0 -> 1 0.88 [0.47, 1.65], p 0.694; 271 rows, 79 censored; 31 rows with a zero value left out. Weibull AFT, right-censored, robust (not cluster) SE; conditional on reach > 0.

## Sensitivity runs (8.8 and section 14)

Model A reach height and Model B reach height, 14.15 interval. All runs, all outcomes: `stats_results.csv`. Forest plots: `fig_forest_plot_A_height.png`, `fig_forest_plot_B_height.png`.

| Run | Group | Note | A rows | A height beta1 [interval], p(0) | A log1p distance beta1 [interval], p(0) | B stations | B height beta1 [interval], p(0) |
|---|---|---|---:|---|---|---:|---|
| primary | primary | combined label, FABDEM, 100 m belt, trimmed A | 380 | -0.52 [-1.70, 0.74], p(0) = 0.309 | -0.31 [-0.69, 0.12], p(0) = 0.114 | 239 | 2.14 [-0.31, 5.16], p(0) = 0.078 |
| label optical_only | label | 8.8 | 366 | -0.14 [-1.24, 1.07], p(0) = 0.797 | -0.25 [-0.68, 0.24], p(0) = 0.238 | 227 | 1.27 [-2.62, 6.04], p(0) = 0.480 |
| label sar_only | label | 8.8 | 407 | 1.67 [-0.38, 3.95], p(0) = 0.111 | 0.36 [-0.20, 1.01], p(0) = 0.200 | 278 | 3.32 [-0.79, 9.27], p(0) = 0.141 |
| label sar_strict | label | 14.12 | 407 | 0.62 [-0.48, 2.00], p(0) = 0.292 | 0.03 [-0.30, 0.42], p(0) = 0.830 | 278 | 1.20 [-0.32, 3.50], p(0) = 0.156 |
| label sar_near | label | 14.17 | 407 | 2.25 [-0.01, 4.93], p(0) = 0.051 | 0.57 [0.11, 1.09], p(0) = 0.019 | 278 | 3.70 [-0.17, 9.25], p(0) = 0.063 |
| DEM Copernicus GLO-30 | DEM | 8.8; outcomes and DEM controls from GLO-30 | 344 | 0.15 [-2.57, 2.80], p(0) = 0.897 | -0.13 [-0.61, 0.42], p(0) = 0.613 | 232 | 2.15 [-3.39, 6.87], p(0) = 0.348 |
| belt 50 m | belt width | 8.8 | 366 | -1.75 [-4.95, 1.78], p(0) = 0.234 | -0.79 [-1.34, -0.07], p(0) = 0.032 | 239 | 0.03 [-3.36, 2.22], p(0) = 0.979 |
| belt 200 m | belt width | 8.8; rows with usable_length_m >= 200 (14.21) | 376 | 0.88 [-0.53, 2.53], p(0) = 0.218 | 0.19 [-0.12, 0.53], p(0) = 0.201 | 235 | 2.54 [-1.12, 6.80], p(0) = 0.163 |
| forest WorldCover | forest definition | 8.8: WorldCover class 10 share | 456 | -2.43 [-5.56, 0.56], p(0) = 0.097 | -0.84 [-1.21, -0.45], p(0) = 0.003 | 239 | 1.59 [-3.29, 6.27], p(0) = 0.466 |
| forest Dynamic World only | forest definition | 14.14: no canopy-height condition; the two rules disagree on 0.0% of belt pixels, so this run nearly equals the dynamic_world-only run | 317 | -0.53 [-2.22, 1.15], p(0) = 0.416 | -0.25 [-0.60, 0.24], p(0) = 0.233 | 225 | 1.80 [-1.06, 4.93], p(0) = 0.167 |
| forest_source dynamic_world only | forest definition | 14.23 a | 317 | -0.53 [-2.22, 1.16], p(0) = 0.416 | -0.25 [-0.60, 0.24], p(0) = 0.232 | 225 | 1.79 [-1.06, 4.93], p(0) = 0.170 |
| steep: bank-top slope <= 45 deg | steep terrain | 8.8 (user decision): pre-event terrain, 10 x tan(bank-top slope) <= 10 m | 377 | -0.44 [-1.61, 0.79], p(0) = 0.366 | -0.31 [-0.71, 0.12], p(0) = 0.114 | 237 | 2.04 [-0.31, 5.00], p(0) = 0.082 |
| steep: unc <= 10 m plus all zero-reach rows | steep terrain | 8.8 variant (user decision): keeps every zero-reach row | 379 | -0.44 [-1.63, 0.81], p(0) = 0.362 | -0.30 [-0.68, 0.13], p(0) = 0.125 | 238 | 2.07 [-0.29, 5.00], p(0) = 0.080 |
| steep: literal unc <= 10 m | steep terrain | selects rows by outcome (conflicts with 14.5); not a headline result | 261 | -0.07 [-4.69, 2.77], p(0) = 0.958 | -0.04 [-0.27, 0.30], p(0) = 0.763 | 174 | 0.01 [-2.84, 3.10], p(0) = 0.991 |
| terraces: bank slope <= 15 deg | terrain | 8.8 | 240 | -0.64 [-5.56, 6.31], p(0) = 0.780 | -0.33 [-0.79, 0.11], p(0) = 0.110 | 46 | 13.44 [-9.69, 53.72], p(0) = 0.125 |
| reach end 30.5 km (5 km upstream) | reach | 14.5 | 317 | -0.15 [-1.69, 1.68], p(0) = 0.802 | -0.17 [-0.63, 0.45], p(0) = 0.410 | 206 | 2.69 [0.52, 5.91], p(0) = 0.021 |
| detached reach | outcome | 14.6: detached_reach_m as the distance outcome | 380 | -0.52 [-1.70, 0.74], p(0) = 0.309 | -0.36 [-1.02, 0.31], p(0) = 0.239 | 239 | 2.14 [-0.31, 5.16], p(0) = 0.078 |
| winsorised outcomes | outcome | 14.19 f: 97.5th percentile caps | 380 | -0.46 [-1.62, 0.88], p(0) = 0.346 | -0.31 [-0.69, 0.12], p(0) = 0.115 | 239 | 2.07 [-0.20, 4.86], p(0) = 0.070 |
| Model A untrimmed | overlap | 14.8: no common-support trimming | 542 | -0.99 [-2.79, 0.79], p(0) = 0.186 | -0.35 [-0.76, 0.12], p(0) = 0.113 | 239 | 2.14 [-0.31, 5.16], p(0) = 0.078 |
| Model B: both banks survived trimming | overlap | 14.8 | 380 | -0.52 [-1.70, 0.74], p(0) = 0.309 | -0.31 [-0.69, 0.12], p(0) = 0.114 | 122 | 3.47 [-0.01, 7.07], p(0) = 0.051 |
| first section only | river section | 14.16: above the first section break | 322 | -0.22 [-1.69, 1.49], p(0) = 0.689 | -0.20 [-0.64, 0.39], p(0) = 0.319 | 209 | 2.69 [0.49, 6.00], p(0) = 0.022 |
| channel_resolved = 1 | channel | 14.18 | 213 | 1.79 [-2.45, 6.71], p(0) = 0.372 | 0.33 [-0.30, 1.24], p(0) = 0.351 | 175 | 3.07 [-1.33, 8.34], p(0) = 0.143 |
| manual exclusions included | exclusions | 14.19 e | 389 | -0.70 [-2.56, 1.15], p(0) = 0.408 | -0.39 [-0.74, 0.08], p(0) = 0.075 | 255 | 2.94 [-0.17, 6.67], p(0) = 0.066 |
| river_not_visible rows as reach 0 | zeros | 14.20: 11 rows added | 398 | -0.48 [-1.66, 0.76], p(0) = 0.334 | -0.29 [-0.65, 0.15], p(0) = 0.149 | 243 | 2.11 [-0.20, 4.80], p(0) = 0.073 |
| near_bank_nodata rows as reach 0 | zeros | 14.24: 11 rows added | 380 | -0.45 [-1.73, 1.10], p(0) = 0.409 | -0.33 [-0.80, 0.15], p(0) = 0.143 | 248 | 2.74 [0.44, 5.55], p(0) = 0.025 |
| usable length >= 300 m | neighbours | 14.21 | 330 | -0.64 [-2.32, 1.05], p(0) = 0.346 | -0.37 [-0.82, 0.22], p(0) = 0.152 | 214 | 2.81 [-0.13, 6.78], p(0) = 0.066 |
| wall slope: complete cases | wall slope | user decision: rows with a valley wall only | 263 | -0.33 [-3.51, 3.19], p(0) = 0.807 | -0.31 [-0.73, 0.45], p(0) = 0.233 | 157 | 2.80 [-0.69, 6.11], p(0) = 0.099 |
| stations with a visible river in the placebo | placebo comparability | (b): real event restricted to stations whose river was visible in the placebo | 229 | -1.45 [-5.94, 1.77], p(0) = 0.353 | -0.53 [-0.98, 0.00], p(0) = 0.050 | 161 | -1.75 [-8.87, 5.04], p(0) = 0.529 |
| station spacing 50 m (8.8) | rebuilt | rebuilt from label to table (melamchi2021__spacing50) | 746 | -1.59 [-3.48, -0.13], p(0) = 0.039 | -0.52 [-1.12, -0.06], p(0) = 0.031 | 455 | 0.96 [-0.95, 3.38], p(0) = 0.274 |
| stations shifted 50 m downstream (14.22) | rebuilt | rebuilt from label to table (melamchi2021__offset50) | 343 | -2.17 [-3.80, -0.53], p(0) = 0.019 | -0.66 [-1.43, -0.01], p(0) = 0.047 | 239 | -0.95 [-4.02, 4.17], p(0) = 0.618 |
| t_bsi x 0.8 (8.8) | rebuilt | rebuilt from label to table (melamchi2021__tbsi_minus20) | 380 | -0.53 [-1.73, 0.71], p(0) = 0.296 | -0.24 [-0.64, 0.19], p(0) = 0.208 | 239 | 2.13 [-0.31, 5.10], p(0) = 0.078 |
| t_bsi x 1.2 (8.8) | rebuilt | rebuilt from label to table (melamchi2021__tbsi_plus20) | 380 | -0.59 [-1.61, 0.48], p(0) = 0.198 | -0.33 [-0.74, 0.13], p(0) = 0.120 | 239 | 1.82 [-1.07, 5.44], p(0) = 0.200 |
| SAR z threshold 2.4 (8.8) | rebuilt | rebuilt from label to table (melamchi2021__sarz_minus20) | 393 | -1.16 [-3.24, 1.03], p(0) = 0.207 | -0.38 [-0.74, 0.04], p(0) = 0.064 | 246 | 0.36 [-3.82, 5.56], p(0) = 0.843 |
| SAR z threshold 3.6 (8.8) | rebuilt | rebuilt from label to table (melamchi2021__sarz_plus20) | 376 | -0.07 [-0.88, 0.80], p(0) = 0.825 | -0.23 [-0.63, 0.22], p(0) = 0.249 | 235 | 1.98 [-0.32, 4.97], p(0) = 0.087 |
| reach end +5 km downstream | reach | 14.5 | | not run: the verified centerline ends at reach_end_km (35.5486 km) | | | |

(d) The real event's 14.20 sensitivity (river_not_visible rows included as reach 0) is reported: yes, row above.

## Placebo through the same models (8.9)

Same Model A and Model B code and the same 14.15 interval as the real event.

| Variant | Note | A rows | A reach height | A log1p distance | A overtopped | B stations | B reach height | B log1p distance |
|---|---|---:|---|---|---|---:|---|---|
| placebo: kept rows (primary) |  | 217 | -0.72 [-1.90, 0.17], p(0) = 0.130 | -0.38 [-0.86, -0.00], p(0) = 0.050 | not estimable: 4 rows in the smaller outcome class, fewer than the 20 model terms | 140 | -0.85 [-2.08, 0.16], p(0) = 0.065 | -0.39 [-0.80, -0.08], p(0) = 0.020 |
| placebo: + 14.20 assumed-zero rows | pushes toward null by construction | 340 | -0.26 [-0.87, 0.18], p(0) = 0.288 | -0.22 [-0.51, 0.03], p(0) = 0.083 | not estimable: 4 rows in the smaller outcome class, fewer than the 20 model terms | 220 | -0.55 [-1.18, 0.16], p(0) = 0.089 | -0.25 [-0.49, -0.06], p(0) = 0.016 |
| placebo: 2020 forest |  | 207 | -0.57 [-2.39, 0.65], p(0) = 0.407 | -0.44 [-1.45, 0.13], p(0) = 0.166 | not estimable: 3 rows in the smaller outcome class, fewer than the 20 model terms | 140 | -0.85 [-3.81, 1.17], p(0) = 0.478 | -0.39 [-1.14, 0.07], p(0) = 0.131 |

(c) The placebo drops 42% of High-forest rows for river_not_visible (99 of 235), against 24% of Low-forest rows (56 of 235) (terciles of the 2018 forest over all in-reach rows). The placebo therefore cannot speak for the densest forest.

The 14.20 assumed-zero variant adds those rows with reach 0; it pushes toward null by construction, so it is not evidence that the placebo is null.

## Files

- `stats_results.csv`, `stats_summary.json`
- `fig_propensity.png`, `fig_model_b_forest_difference.png`, `fig_forest_plot_A_height.png`, `fig_forest_plot_B_height.png`

