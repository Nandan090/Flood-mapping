# Validation comparison report

## What was run

Compared 61 hand-digitised edge points from 31 stations with the primary combined-label/FABDEM reach-distance table output. The signed error is computed reach minus validated reach; negative values mean under-detection. The main summary includes 37 points in 6 3 km blocks.

Files produced: `validation_per_point.csv`, `validation_summary.csv`, `validation_cutoff_sensitivity.csv`, `validation_detached_sensitivity.csv`, `validation_plot.png`, `limitations.md`, and this report.

## Main results

| Group | n | Mean diff (m) | Median diff (m) | MAE (m) | Within 20 m | Under-detected |
|---|---:|---:|---:|---:|---:|---:|
| Low | 12 | -11.70 | -11.31 | 21.42 | 58.3% | 25.0% |
| Mid | 13 | -14.83 | -9.46 | 17.75 | 76.9% | 23.1% |
| High | 12 | -10.68 | -13.21 | 12.99 | 91.7% | 8.3% |

### Main result beside the imagery-review sensitivity

Under-detected = computed reach more than 20 m shorter than the traced edge.

| Set | n | Mean diff (m) | Low under-detected | Mid under-detected | High under-detected | High minus Low (points) | Flag (>10 points) | OLS slope (m per unit F) | Spearman rho |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| Main | 37 | -12.47 | 3 of 12 | 3 of 13 | 1 of 12 | -16.7 | off | 18.135 | 0.265 |
| Without "no visible change" | 31 | -6.04 | 1 of 10 | 0 of 10 | 0 of 11 | -10.0 | off | 1.258 | 0.170 |

No visible change can mean the traced edge is not a flood edge, or that the flood left no visible trace (for example under trees). Imagery cannot separate these, so this row does not show that the bias is absent. The main row stays primary.

Points reviewed as "no visible change" and left out of the second row: 69-right, 249-right, 251-right, 252-right, 253-right, 283-left. Review status counts over all 61 points: visible change 3, no visible change 6, can't tell 0, not reviewed 52. The slope and rho in this table are point estimates; the bootstrap intervals are for the main set only.

### Trend with forest fraction

- OLS slope of diff on forest fraction: **18.135 m per unit F**; 95% interval -9.837 to 70.374 m per unit F (approximate, few validation points).
- Spearman rho = 0.265; too few clusters (6); no interval.

Both intervals use the 14.15 null-imposed wild cluster bootstrap-t (CR1, Rademacher, 1999 draws, 3 km block IDs).

The OLS slope interval is wide and includes 0, so these data cannot rule out either no trend or a strong negative trend in the error with forest fraction. This is not evidence that there is no forest bias.

### Forest-bias check

Under-detected share (computed reach more than 20 m shorter than the traced edge): High 8.3% (1 of 12) vs Low 25.0% (3 of 12). High minus Low: -16.7 percentage points. Not flagged (>10 percentage points).

### Under-detection cut-off sensitivity

The 20 m cut-off was fixed before the results were seen and stays primary. The flag is the same at every cut-off shown (not raised), but the gap itself changes with the cut-off.

| Cut-off (m) | Low | Mid | High | High minus Low (points) | Flag (>10 points) |
|---:|---:|---:|---:|---:|---|
| 10 | 50.0% (6 of 12) | 46.2% (6 of 13) | 58.3% (7 of 12) | 8.3 | no |
| 20 (primary) | 25.0% (3 of 12) | 23.1% (3 of 13) | 8.3% (1 of 12) | -16.7 | no |
| 30 | 16.7% (2 of 12) | 7.7% (1 of 13) | 0.0% (0 of 12) | -16.7 | no |

### Date and confidence breakdowns

- 2025-11-29: n=13, mean diff -22.50 m, MAE 23.90 m, under-detected 30.8%.
- 2026-01-02: n=13, mean diff -3.78 m, MAE 16.12 m, under-detected 15.4%.
- 2026-01-20: n=11, mean diff -10.88 m, MAE 11.22 m, under-detected 9.1%.
- Confidence 1: n=37, mean diff -12.47 m, MAE 17.40 m.
- Confidence 2: n=0, mean diff NA m, MAE NA m.
- Confidence 3: n=0, mean diff NA m, MAE NA m.

### Largest errors

The 5 main-set points with the largest absolute error. All are kept in every summary.

| Station | Side | Forest group | Confidence | diff (m) |
|---:|---|---|---:|---:|
| 251 | right | Mid | 1 | -90.94 |
| 252 | right | Low | 1 | -63.48 |
| 253 | right | Low | 1 | -46.57 |
| 68 | right | Low | 1 | 32.58 |
| 69 | right | Mid | 1 | -26.98 |

### Exclusions

Excluded from main numbers (24 point rows). A row with several reasons is listed under each:

- Censored (reach edge not observed): 52-left, 53-left, 55-right, 59-left, 62-left, 331-left, 343-left.
- Near-bank nodata: 75-left, 249-left, 252-left, 296-right, 311-right, 317-right.
- Low valid fraction: 52-right, 53-right, 60-left, 75-left, 89-left, 249-left, 252-left, 253-left, 260-right, 272-right, 283-right, 295-right, 296-right, 311-right, 315-right.
- Other model exclusion: 52-right (reservoir_or_dam), 311-right (tributary_fan), 331-right (reservoir_or_dam).
- Validated reach unavailable: none.
- Computed reach unavailable: 75-left, 249-left, 252-left, 296-right, 311-right, 317-right.
- Forest fraction unavailable: none.

### Sensitivity checks

| Sensitivity | n | Mean diff (m) | MAE (m) | Under-detected | OLS slope | Spearman rho |
|---|---:|---:|---:|---:|---:|---:|
| Without confidence 3 | 37 | -12.47 | 17.40 | 18.9% | 18.135 | 0.265 |
| Near-bank nodata set to zero | 38 | -12.51 | 17.31 | 18.4% | 17.446 | 0.275 |

### Detached reach sensitivity (not part of the main summary)

Errors when the computed reach is `detached_reach_m` (furthest connected impact on the transect, gaps ignored, 14.6) instead of `reach_distance_m`. diff is still computed minus validated. Censored rows are shown separately; they are not in the main set.

| Set and computed outcome | n | Mean diff (m) | Median diff (m) | MAE (m) | Under-detected |
|---|---:|---:|---:|---:|---:|
| Main set, reach_distance_m (primary) | 37 | -12.47 | -12.45 | 17.40 | 18.9% (7 of 37) |
| Main set, detached_reach_m | 37 | -12.47 | -12.45 | 17.40 | 18.9% (7 of 37) |
| Censored rows, reach_distance_m | 7 | 6.27 | 11.24 | 12.72 | 0.0% (0 of 7) |
| Censored rows, detached_reach_m | 7 | 7.69 | 11.24 | 11.29 | 0.0% (0 of 7) |

| Station | Side | Forest group | Validated (m) | reach_distance_m | detached_reach_m |
|---:|---|---|---:|---:|---:|
| 52 | left | High | -3.7 | 10.0 | 10.0 |
| 53 | left | Mid | 35.1 | 20.0 | 30.0 |
| 55 | right | Low | 30.5 | 50.0 | 50.0 |
| 59 | left | Mid | 5.2 | 10.0 | 10.0 |
| 62 | left | Mid | -1.2 | 10.0 | 10.0 |
| 331 | left | Low | 22.8 | 40.0 | 40.0 |
| 343 | left | Mid | 17.5 | 10.0 | 10.0 |

Negative diff means the computed reach is shorter than the digitised reach. These descriptive validation results do not establish that forest causes a difference.

### Notes

- Stations 249, 251, 252 and 253 (river km 25.0 to 25.4) lie within about 340 m of each other and count as one cluster (they are in the same 3 km block).
- Points outside the main set whose computed and validated reach differ by more than 20 m were not reviewed: 52-right, 260-right, 272-right, 315-right.

## Still needed

Review the station-level comparison and limitations, then complete the study's separate statistical, placebo, and reporting steps. Validation is one check of the impact label; it does not by itself establish the forest-reach association.

See `limitations.md` for additional data and coverage limitations.
