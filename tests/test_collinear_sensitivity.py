"""Model A without the collinear controls: a descriptive sensitivity (not in the verdict)."""

import pandas as pd

from floodforest.stats_run import COLLINEAR_CONTROLS, COLLINEAR_RUN, collinear_rows, without_terms

RHS = ("forest_frac + np.log(channel_width_m) + gradient + log10_upa + bank_elev_m"
       " + C(side) + cr(river_km, df=5) + valley_wall_found + channel_resolved")


def test_collinear_controls_are_removed_and_nothing_else():
    assert COLLINEAR_CONTROLS == ("bank_elev_m", "log10_upa", "valley_wall_found")
    assert without_terms(RHS, COLLINEAR_CONTROLS) == (
        "forest_frac + np.log(channel_width_m) + gradient + C(side) + cr(river_km, df=5)"
        " + channel_resolved")


def test_collinear_rows_put_the_primary_beside_the_sensitivity():
    primary = pd.DataFrame([
        {"run": "primary", "model": "A", "outcome": "y_height", "rows": 10, "beta": -1.0,
         "lo": -2.0, "hi": 0.5, "p0": 0.2, "interval_scale": "beta"},
        {"run": "primary", "model": "A", "outcome": "y_logdist", "rows": 10, "beta": -0.1,
         "lo": -0.3, "hi": 0.1, "p0": 0.4, "interval_scale": "beta"},
        {"run": "primary", "model": "A", "outcome": "y_over", "rows": 10, "error": "not estimable"},
    ])
    sens = {"rhs": "forest_frac + gradient", "rows": 10, "A": {
        "y_height": {"beta": -0.8, "interval": (-1.9, 0.4), "p0": 0.3},
        "y_logdist": {"error": "boom"},
        "y_over": {"beta": 0.1, "odds_ratio": 1.1, "or_interval": (0.5, 2.0)}}}
    out = collinear_rows(primary, RHS, sens)
    assert out.run.tolist() == ["primary", COLLINEAR_RUN] * 3
    height = out[out.outcome == "y_height"].set_index("run")
    assert height.loc["primary", "beta"] == -1.0 and height.loc[COLLINEAR_RUN, "hi"] == 0.4
    assert out[(out.outcome == "y_logdist") & (out.run == COLLINEAR_RUN)].error.iloc[0] == "boom"
    assert out[(out.outcome == "y_over") & (out.run == COLLINEAR_RUN)].odds_ratio.iloc[0] == 1.1
