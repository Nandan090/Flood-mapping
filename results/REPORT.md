# Validation comparison report

## What was run

Compared 52 hand-digitised edge points from 26 stations with the primary combined-label/FABDEM reach-distance table output. The signed error is computed reach minus validated reach; negative values mean under-detection. The main summary includes 35 points in 9 3 km blocks.

Files produced: `validation_per_point.csv`, `validation_summary.csv`, `validation_cutoff_sensitivity.csv`, `validation_detached_sensitivity.csv`, `validation_plot.png`, `limitations.md`, and this report.

## Main results

| Group | n | Mean diff (m) | Median diff (m) | MAE (m) | Within 20 m | Under-detected |
|---|---:|---:|---:|---:|---:|---:|
| Low | 10 | -3.17 | -2.06 | 6.44 | 90.0% | 10.0% |
| Mid | 10 | -25.03 | -12.95 | 26.63 | 80.0% | 20.0% |
| High | 15 | -12.20 | -13.83 | 14.29 | 73.3% | 26.7% |

### Main result beside the imagery-review sensitivity

Under-detected = computed reach more than 20 m shorter than the traced edge.

| Set | n | Mean diff (m) | Low under-detected | Mid under-detected | High under-detected | High minus Low (points) | Flag (>10 points) | OLS slope (m per unit F) | Spearman rho |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| Main | 35 | -13.29 | 1 of 10 | 2 of 10 | 4 of 15 | 16.7 | on | 1.159 | -0.325 |
| Without "no visible change" | 29 | -6.50 | 1 of 10 | 0 of 8 | 0 of 11 | -10.0 | off | -4.317 | -0.313 |

No visible change can mean the traced edge is not a flood edge, or that the flood left no visible trace (for example under trees). Imagery cannot separate these, so this row does not show that the bias is absent. The main row stays primary.

Points reviewed as "no visible change" and left out of the second row: 30-left, 105-right, 115-left, 118-left, 281-left, 346-right. Review status counts over all 52 points: visible change 1, no visible change 6, can't tell 0, not reviewed 45. The slope and rho in this table are point estimates; the bootstrap intervals are for the main set only.

### Trend with forest fraction

- OLS slope of diff on forest fraction: **1.159 m per unit F**; 95% interval -30.794 to 31.793 m per unit F (approximate, few validation points).
- Spearman rho = **-0.325**; 95% interval -1.000 to 0.188 (approximate, few validation points; clipped to [-1, 1], the raw bootstrap lower bound was -1.170).

Both intervals use the 14.15 null-imposed wild cluster bootstrap-t (CR1, Rademacher, 1999 draws, 3 km block IDs).

The OLS slope and Spearman rho intervals are wide and include 0, so these data cannot rule out either no trend or a strong negative trend in the error with forest fraction. This is not evidence that there is no forest bias.

### Forest-bias check

Under-detected share (computed reach more than 20 m shorter than the traced edge): High 26.7% (4 of 15) vs Low 10.0% (1 of 10). High minus Low: 16.7 percentage points. Flagged (>10 percentage points); the main result may be biased toward finding protection.

### Under-detection cut-off sensitivity

The 20 m cut-off was fixed before the results were seen and stays primary. The flag is sensitive to the cut-off: it is raised at 10 m and 20 m, and not at 30 m.

| Cut-off (m) | Low | Mid | High | High minus Low (points) | Flag (>10 points) |
|---:|---:|---:|---:|---:|---|
| 10 | 10.0% (1 of 10) | 60.0% (6 of 10) | 60.0% (9 of 15) | 50.0 | yes |
| 20 (primary) | 10.0% (1 of 10) | 20.0% (2 of 10) | 26.7% (4 of 15) | 16.7 | yes |
| 30 | 0.0% (0 of 10) | 10.0% (1 of 10) | 0.0% (0 of 15) | 0.0 | no |

### Date and confidence breakdowns

- 2022-01-17: n=17, mean diff -11.09 m, MAE 12.93 m, under-detected 23.5%.
- 2022-03-26: n=18, mean diff -15.36 m, MAE 18.06 m, under-detected 16.7%.
- Confidence 1: n=31, mean diff -12.43 m, MAE 15.01 m.
- Confidence 2: n=4, mean diff -19.90 m, MAE 19.90 m.
- Confidence 3: n=0, mean diff NA m, MAE NA m.

### Largest errors

The 5 main-set points with the largest absolute error. All are kept in every summary.

| Station | Side | Forest group | Confidence | diff (m) |
|---:|---|---|---:|---:|
| 346 | right | Mid | 1 | -150.53 |
| 281 | left | Mid | 1 | -29.86 |
| 115 | left | High | 2 | -28.65 |
| 328 | right | Low | 1 | -26.26 |
| 118 | left | High | 1 | -25.73 |

Station 328 is kept in every main summary and flagged descriptively as `braided reach, far from main channel`; the sensitivity run below repeats the summary without it.

### Exclusions

Excluded from main numbers (17 point rows). A row with several reasons is listed under each:

- Censored (reach edge not observed): 13-right, 34-left, 38-left, 38-right, 265-left, 266-left, 279-right, 288-left, 322-left, 323-right.
- Near-bank nodata: 13-left, 20-left, 30-right, 105-left, 322-right.
- Low valid fraction: 13-left, 105-left, 118-right, 288-left, 317-right, 322-right, 323-right.
- Other model exclusion: 38-right (hillslope_failure).
- Validated reach unavailable: none.
- Computed reach unavailable: 13-left, 20-left, 30-right, 105-left, 322-right.
- Forest fraction unavailable: none.

### Sensitivity checks

| Sensitivity | n | Mean diff (m) | MAE (m) | Under-detected | OLS slope | Spearman rho |
|---|---:|---:|---:|---:|---:|---:|
| Without confidence 3 | 35 | -13.29 | 15.57 | 20.0% | 1.159 | -0.325 |
| Near-bank nodata set to zero | 37 | -13.62 | 15.78 | 18.9% | 0.310 | -0.356 |
| Without station 328 | 33 | -13.25 | 15.67 | 18.2% | 1.149 | -0.369 |

### Detached reach sensitivity (not part of the main summary)

Errors when the computed reach is `detached_reach_m` (furthest connected impact on the transect, gaps ignored, 14.6) instead of `reach_distance_m`. diff is still computed minus validated. Censored rows are shown separately; they are not in the main set.

| Set and computed outcome | n | Mean diff (m) | Median diff (m) | MAE (m) | Under-detected |
|---|---:|---:|---:|---:|---:|
| Main set, reach_distance_m (primary) | 35 | -13.29 | -8.24 | 15.57 | 20.0% (7 of 35) |
| Main set, detached_reach_m | 35 | -13.29 | -8.24 | 15.57 | 20.0% (7 of 35) |
| Censored rows, reach_distance_m | 10 | -25.39 | -7.43 | 28.88 | 20.0% (2 of 10) |
| Censored rows, detached_reach_m | 10 | -6.39 | -5.54 | 11.37 | 10.0% (1 of 10) |

| Station | Side | Forest group | Validated (m) | reach_distance_m | detached_reach_m |
|---:|---|---|---:|---:|---:|
| 13 | right | High | 7.7 | 5.0 | 5.0 |
| 34 | left | High | 40.4 | 45.0 | 45.0 |
| 38 | left | High | 46.9 | 35.0 | 35.0 |
| 38 | right | High | 112.3 | 75.0 | 75.0 |
| 265 | left | Mid | 56.9 | 40.0 | 40.0 |
| 266 | left | Mid | 75.1 | 70.0 | 70.0 |
| 279 | right | Mid | 37.2 | 50.0 | 50.0 |
| 288 | left | Low | 26.0 | 20.0 | 20.0 |
| 322 | left | Mid | 58.9 | 50.0 | 50.0 |
| 323 | right | High | 202.6 | 20.0 | 210.0 |

Negative diff means the computed reach is shorter than the digitised reach. These descriptive validation results do not establish that forest causes a difference.

## Still needed

Review the station-level comparison and limitations, then complete the study's separate statistical, placebo, and reporting steps. Validation is one check of the impact label; it does not by itself establish the forest-reach association.

See `limitations.md` for additional data and coverage limitations.
