# Validation limitations

- Low-forest coverage: 8 of 14 sampled stations were skipped. The coverage file does not record station-specific reasons. The digitising instructions say to skip when cloud, shadow, or trees hide the edge, so the imagery protocol is clear but the cause for each skipped station cannot be determined from the saved data.
- Station 328 is retained and flagged as `braided reach, far from main channel`; its original note about the imagery dot being a couple of metres off is preserved.
- The 52 validation points were hand-marked from Esri Wayback imagery. The recorded capture dates are 2022-01-17, 2022-03-26; the Wayback release is 2024-11-18.
- The points cover 26 stations from the 40-station sample. The wild-bootstrap analysis has 9 blocks among its eligible points; intervals are labelled approximate because the validation sample is small.
- This is a transect-distance comparison against the primary combined-label/FABDEM table outcome, not a full-raster intersection-over-union assessment. Hand-digitised edges and the computed label both have measurement uncertainty.
- Skipped low-forest stations weaken coverage of that group; comparisons with any group having fewer than 8 eligible points are descriptive only.
- Validation is weaker where the flood edge is not visible in imagery: on wide bare gravel plains, where pre-event bare ground is unknown to the satellite label and a traced edge may only mark the old plain boundary, and on steep forested slopes, where a flood that passed under trees may leave no bare deposit. Six of the seven largest main-set errors are 21 to 30 m, within a few 10 m pixels.
- Only points with an error above 20 m were reviewed against pre-event imagery; smaller errors were not.
