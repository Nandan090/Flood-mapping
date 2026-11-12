# Validation limitations

- Low-forest coverage: 12 of 17 sampled low-forest stations were assigned in the subset drawn on 2026-11-12 (seed 43); 1 of those was skipped or removed. Reasons for stations removed after digitising are recorded in provenance (validation_edges).
- The 61 validation points were hand-marked from Google Earth Pro historical imagery. The recorded capture dates are 2025-11-29, 2026-01-02, 2026-01-20.
- The points cover 31 stations from the 50-station sample. The wild-bootstrap analysis has 6 blocks among its eligible points; intervals are labelled approximate because the validation sample is small.
- This is a transect-distance comparison against the primary combined-label/FABDEM table outcome, not a full-raster intersection-over-union assessment. Hand-digitised edges and the computed label both have measurement uncertainty.
- Skipped low-forest stations weaken coverage of that group; comparisons with any group having fewer than 8 eligible points are descriptive only.
- Validation is weaker where the flood edge is not visible in imagery: on wide bare gravel plains, where pre-event bare ground is unknown to the satellite label and a traced edge may only mark the old plain boundary, and on steep forested slopes, where a flood that passed under trees may leave no bare deposit. Five of the nine largest main-set errors are 23 to 27 m, within a few 10 m pixels.
- Only points with an error above 20 m were reviewed against pre-event imagery; smaller errors were not.
