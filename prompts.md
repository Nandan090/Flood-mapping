No, don't paste the spec into the chat as a prompt. Save it over the `CURSOR_BUILD_SPEC.md` in your project folder, then send the short prompt in row 2b, which points Cursor to the file. The attachment you sent matches my latest version, so it's the right one to save.

Each row is a new chat. Start every prompt with "Follow @CURSOR_BUILD_SPEC.md, including section 14, which overrides earlier sections" and end it with "Stop when done and tell me what you built and what you need from me."

| # | Model | Prompt | Your job afterwards |
|---|---|---|---|
| 1 | — | Done: plan and `info.md`. | — |
| 2 | — | Done: scaffolding, config, tests. | — |
| **2b** | Opus 5.5 | "Sections 14.15 and 14.16 are new. (1) Implement 14.15 in `stats.py`; it replaces 8.4 as the primary uncertainty method. (2) Revert the coverage test to the settings frozen in 14.15: reach lengths of 24, 36 and 72 km, not 80 segments. (3) Keep the old 1 km percentile bootstrap as the labelled sensitivity output. (4) Add `river_section`, `section_breaks_km`, `landmarks`, `bootstrap_block_km` and the shared-centerline rule from 14.16 to the config and tests, and update the reach descriptions in `info.md`. (5) Run `pytest` and `pytest -m slow` and show me the coverage percentage at each reach length. If coverage is outside 92–98%, do not change the test: stop and show me the numbers." | **Bring me the three coverage numbers.** Commit to git if they pass. |
| 3 | Grok 4.7 | "Implement `ee_data.py` and the `fetch` step for `melamchi2021` only. Verify every asset ID at runtime (rule 0.3), enforce timestamps per 14.1, apply the Hansen rule in 14.3, and write `provenance.json`. Also implement the OSM centerline draft helper from section 6, covering the Melamchi River and the Indrawati down to the reach end in 14.16." | Set up Earth Engine first. Correct the draft river line in QGIS, choose the start point, and mark it verified. |
| 4 | Opus 5.5 | "Implement `optical.py`, `sar.py` and `label.py` exactly as section 5, with the reference zone from 14.11 and the agreement rule plus `sar_strict` variant from 14.12. Run `label` for `melamchi2021` and show me the pre/post images and the three label maps." | **Pause.** Compare the label maps with the satellite images. Bring screenshots here if unsure. |
| 5 | Grok 4.7 | "Implement `centerline.py` per section 6 and 14.16: stations, normals, transects, all exclusion flags, reach start and end, section breaks, and the printed distances to landmarks." | Check transects in QGIS. Confirm the Indrawati junction distance and set `section_breaks_km`. |
| 6 | Opus 5.5 | "Implement `profiles.py`, `exposure.py`, `controls.py` and the `table` step per section 7, with the run-start rule and `detached_reach_m` from 14.6, for both DEMs. Produce the six example cross-section plots chosen by the fixed rule in section 10." | **Pause.** Check the six cross-sections. |
| 7 | Grok 4.7 | "Implement the validation sampler from 14.13: write `sample_stations.geojson` with 40 stations stratified by forest tercile, without forest values." | Trace the flood edge at those 40 stations before seeing any results. |
| 8 | Opus 5.5 | "Implement the rest of `stats.py` as section 8, with 14.7, 14.8, 14.9, 14.10, the intervals from 14.15, and `river_section` plus the first-section-only run from 14.16. Include every sensitivity run in 8.8 and those added in section 14." | Run `pytest -m slow` again. |
| 9 | Grok 4.7 | "Implement `placebo` using the fixed periods in 14.2, and `validate` per section 9 and 14.13." | Confirm no flood hit the Melamchi valley between late 2018 and late 2019. |
| 10 | Grok 4.7 | "Implement `figures.py` and `report.py` per section 10, including the PROVISIONAL banner logic from 14.4. Build the Melamchi report." | Read the report start to finish. |
| 11 | Sonnet 5.5 | "Review `label.py`, `profiles.py`, `exposure.py` and `stats.py` against the spec. List every place the code differs from it, with file and line. Change nothing." | Decide which differences to fix, then fix them in a new Opus chat. |
| 12 | Grok 4.7 | "Run the full pipeline for `lende2025`. Change only config, not code. The centerline is shared with `bhotekoshi2026` and starts at Rasuwagadhi (14.16). Report anything that fails." | Verify the shared river line once, trace 40 edges, check the pause-point outputs. |
| 13 | Grok 4.7 | "Run the full pipeline for `bhotekoshi2026` as a provisional run per 14.4, using the same centerline as `lende2025`, then build the summary report." | Re-run after 15 December for the final verdict. |

Three changes from the last table:

- **Row 2b is new** and comes next. It fixes the uncertainty method and adds the reach settings.
- **One river line now serves both the 2025 and 2026 floods**, so you only verify it once, at row 12.
- **Use standard Grok 4.7, not the Fast variant**, which costs twice as much.