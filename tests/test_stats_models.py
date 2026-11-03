"""Section 8 pieces added for the stats step: model frame, common support (8.1,
14.8), Model B pairs, censored model (14.7), effect sizes and E-values (8.6, 8.7,
14.9), narrowness (14.10), the decision rule (8.10), the sensitivity rebuild
variants (8.8, 14.22) and the runner's row builders."""

import copy

import numpy as np
import pandas as pd
import pytest

from floodforest import stats as S
from floodforest.config import effective_config, load_config, rasters_dir, variant_key

CFG = load_config()


def _table(n=60, seed=0):
    rng = np.random.default_rng(seed)
    st = np.repeat(np.arange(n // 2), 2)
    side = np.tile(["left", "right"], n // 2)
    km = st * 0.5
    wall = (rng.random(n) > 0.3).astype(int)
    return pd.DataFrame({
        "station_id": st, "side": side, "river_km": km, "segment_id": np.floor(km).astype(int),
        "block_id": np.floor(km / 3).astype(int), "river_section": (km > 10).astype(int),
        "channel_resolved": (rng.random(n) > 0.2).astype(int),
        "forest_frac_100": rng.random(n), "channel_width_m": 20 + rng.random(n) * 30,
        "gradient_fab": rng.normal(size=n), "valley_width_m_fab": 50 + rng.random(n) * 300,
        "wall_slope_deg_fab": np.where(wall == 1, 30 + rng.random(n) * 20, np.nan),
        "valley_wall_found_fab": wall, "bank_slope_deg_fab": rng.random(n) * 30,
        "log10_upa": rng.normal(2, 0.1, n), "outer_bend": rng.integers(0, 2, n),
        "curvature": rng.normal(size=n), "bank_elev_m_fab": 800 + rng.random(n) * 50,
        "reach_height_m_fab": rng.random(n) * 5, "reach_distance_m_fab": rng.random(n) * 80,
        "overtopped_fab": rng.integers(0, 2, n), "censored_fab": rng.integers(0, 2, n),
    })


def test_model_frame_fills_wall_slope_and_flags_it():
    t = _table()
    t.loc[3, "curvature"] = np.nan          # a Model A control only
    t.loc[5, "forest_frac_100"] = np.nan    # needed by both models
    df = S.model_frame(t, "forest_frac_100")
    assert df.attrs["dropped_missing"] == 2 and df.attrs["dropped_missing_b"] == 1
    assert len(df) == len(t) - 1                       # row 5 is usable by neither
    assert not df.loc[3, "a_ok"] and df.loc[3, "b_ok"]
    no_wall = df.valley_wall_found == 0
    assert (df.loc[no_wall, "wall_slope_deg"] == 0).all()
    assert (df.loc[~no_wall, "wall_slope_deg"] >= 30).all()
    assert np.allclose(df.y_logdist, np.log1p(df.y_dist))


def test_common_support_keeps_rows_inside_the_overlap():
    t = _table(200, seed=1)
    df = S.model_frame(t, "forest_frac_100")
    cs = S.common_support(df, "bank_slope_deg + gradient", 0.5, 0.15, 0.3, 1)
    lo, hi = cs["bounds"]
    assert ((cs["p"][cs["keep"]] >= lo) & (cs["p"][cs["keep"]] <= hi)).all()
    assert cs["kept"] + cs["dropped"] == len(df)
    assert lo == max(cs["p"][cs["high"]].min(), cs["p"][~cs["high"]].min())


def test_common_support_switches_to_the_median_when_a_group_is_small():
    t = _table(200, seed=2)
    t["forest_frac_100"] = t["forest_frac_100"] * 0.55     # few rows reach 0.5
    df = S.model_frame(t, "forest_frac_100")
    cs = S.common_support(df, "bank_slope_deg", 0.5, 0.15, 0.3, 1)
    assert cs["rule"].startswith("F >= median")
    assert abs(cs["high"].mean() - 0.5) < 0.05


def test_too_thin_overlap_is_flagged():
    df = S.model_frame(_table(200, seed=3), "forest_frac_100")
    cs = S.common_support(df, "bank_slope_deg", 0.5, 0.15, 0.3, 10_000)
    assert cs["too_thin"]


def test_paired_frame_and_informative_rule():
    df = S.model_frame(_table(), "forest_frac_100").iloc[:-1]   # last station loses a bank
    pairs = S.paired_frame(df)
    assert len(pairs) == df.station_id.nunique() - 1
    row = pairs.iloc[0]
    left = df[(df.station_id == row.station_id) & (df.side == "left")].iloc[0]
    right = df[(df.station_id == row.station_id) & (df.side == "right")].iloc[0]
    assert row.d_forest == pytest.approx(left.forest_frac - right.forest_frac)
    assert row.d_y_height == pytest.approx(left.y_height - right.y_height)
    info = S.model_b_informative(pairs, 30, 0.2)
    assert info["stations_large_diff"] == int((pairs.d_forest.abs() >= 0.2).sum())
    assert info["informative"] == (info["stations_large_diff"] >= 30)
    assert S.model_b_rhs().startswith("d_forest + d_valley_width_m")


def test_aft_recovers_a_planted_time_ratio_and_counts_zero_rows():
    rng = np.random.default_rng(4)
    n = 600
    forest = rng.random(n)
    slope = rng.random(n)
    scale = np.exp(3.0 - 1.0 * forest + 0.3 * slope)
    t = scale * rng.weibull(2.0, n)
    cap = np.quantile(t, 0.85)
    df = pd.DataFrame({"forest_frac": forest, "slope": slope, "y_dist": np.minimum(t, cap),
                       "censored": (t >= cap).astype(int)})
    df["y_height"] = df["y_dist"]
    df.loc[:9, "y_height"] = 0.0
    res = S.aft_conditional(df, "y_dist", "forest_frac + slope")
    assert res["time_ratio"] == pytest.approx(np.exp(-1.0), rel=0.2)
    assert res["censored"] == int(df.censored.sum())
    assert S.aft_conditional(df, "y_height", "forest_frac + slope")["zero_outcome_rows_left_out"] == 10


def test_effect_sizes_and_evalues():
    df = pd.DataFrame({"forest_frac": [0.0, 0.25, 0.5, 0.75, 1.0]})
    r = {"beta": -2.0, "interval": (-3.0, -1.0)}
    e = S.effect_sizes(df, r)
    assert e["iqr_F"] == pytest.approx(0.5) and e["effect_per_iqr"] == pytest.approx(-1.0)
    assert e["interval_per_iqr"] == pytest.approx((-1.5, -0.5))
    ev = S.evalue_for(r, sd_y=4.0)
    assert ev["evalue_ci"] == pytest.approx(S.evalue_continuous(-1.0, 4.0, -1.0)["evalue"])
    assert S.evalue_for({"beta": -2.0, "interval": (-3.0, 0.5)}, 4.0)["evalue_ci"] == 1.0
    lg = S.evalue_logistic({"odds_ratio": 0.5, "or_interval": (0.3, 0.8)}, 0.4, 0.15)
    assert lg["evalue_ci"] == pytest.approx(S.evalue_binary(0.8, 0.4, 0.15))
    assert S.evalue_logistic({"odds_ratio": 0.5, "or_interval": (0.3, 1.2)}, 0.4, 0.15)["evalue_ci"] == 1.0
    assert S.or_narrow((0.7, 1.4), 0.67, 1.5) and not S.or_narrow((0.6, 1.4), 0.67, 1.5)


def _checks(**over):
    base = {k: True for k in ("a_height_neg_sig", "b_height_neg_sig", "dist_same_sign",
                              "labels_same_sign", "dems_same_sign", "placebo_null",
                              "validation_ok", "validated", "enough_blocks", "overlap_ok",
                              "model_b_informative", "no_sign_flips")}
    base.update(all_include_null=False, all_narrow=False)
    base.update(over)
    return base


def test_decision_rule():
    assert S.decide_verdict(_checks())[0] == "Consistent with a protective association"
    verdict, reasons = S.decide_verdict(_checks(validation_ok=False))
    assert verdict == "Inconclusive"
    assert any("validation omission" in r for r in reasons)
    nul = dict(a_height_neg_sig=False, b_height_neg_sig=False, dist_same_sign=False,
               labels_same_sign=False, dems_same_sign=False, all_include_null=True)
    assert S.decide_verdict(_checks(**nul, all_narrow=True))[0] == "No evidence of association"
    verdict, reasons = S.decide_verdict(_checks(**nul, all_narrow=False))
    assert verdict == "Inconclusive" and any("reasonably narrow" in r for r in reasons)
    verdict, reasons = S.decide_verdict(_checks(**nul, all_narrow=True, placebo_null=False))
    assert verdict == "Inconclusive" and any("cannot be trusted" in r for r in reasons)
    for cap in ("enough_blocks", "overlap_ok", "model_b_informative", "validated",
                "no_sign_flips"):
        assert S.decide_verdict(_checks(**{cap: False}))[0] == "Inconclusive"
    verdict, reasons = S.decide_verdict(_checks(**nul, all_narrow=True, no_sign_flips=False))
    assert verdict == "Inconclusive" and any("sign flips" in r for r in reasons)


def test_sign_flips_compare_each_run_with_the_primary():
    primary = {("A", "y_height"): {"beta": -0.5}, ("B", "y_height"): {"beta": 2.0},
               ("A", "y_over"): {"error": "not estimable"}}
    runs = {"r1": {("A", "y_height"): {"beta": -1.0, "p0": 0.2},
                   ("B", "y_height"): {"beta": -0.3, "p0": 0.6}},
            "r2": {("A", "y_height"): {"beta": 0.4, "p0": 0.01},
                   ("A", "y_over"): {"beta": 1.0}},
            "r3": {("A", "y_height"): {"error": "x"}}}
    flips = S.sign_flips(primary, runs, [("A", "y_height"), ("B", "y_height"), ("A", "y_over")])
    assert {(f["run"], f["model"], f["outcome"]) for f in flips} == {
        ("r1", "B", "y_height"), ("r2", "A", "y_height")}
    assert next(f for f in flips if f["run"] == "r2")["run_p0"] == 0.01


# --- rebuild variants -------------------------------------------------------

def test_variant_events_change_only_their_overrides():
    cfg = load_config()
    names = {"spacing50", "offset50", "tbsi_minus20", "tbsi_plus20", "sarz_minus20", "sarz_plus20"}
    for name in names:
        assert hasattr(cfg.events, variant_key("melamchi2021", name))
    eff = effective_config(cfg, "melamchi2021__sarz_minus20")
    assert eff.sar_z_threshold == pytest.approx(2.4) and cfg.sar_z_threshold == 3.0
    assert effective_config(cfg, "melamchi2021__spacing50").station_spacing_m == 50
    assert effective_config(cfg, "melamchi2021__offset50").station_offset_m == 50
    assert effective_config(cfg, "melamchi2021__tbsi_plus20").t_bsi_factor == pytest.approx(1.2)
    assert effective_config(cfg, "melamchi2021") is cfg
    assert rasters_dir(cfg, "melamchi2021__spacing50") == rasters_dir(cfg, "melamchi2021")
    v = cfg.events.melamchi2021__offset50
    assert v.manual_exclusions == cfg.events.melamchi2021.manual_exclusions
    assert cfg.station_offset_m == 0 and cfg.t_bsi_factor == 1.0


def test_unknown_override_is_refused():
    cfg = load_config()
    bad = copy.deepcopy(cfg)
    bad.events.melamchi2021__spacing50.cfg_overrides = {"no_such_setting": 1}
    with pytest.raises(KeyError, match="not a config.yaml setting"):
        effective_config(bad, "melamchi2021__spacing50")


def test_station_offset_moves_every_station_downstream():
    from floodforest.geometry import stations
    line = np.column_stack([np.zeros(101), np.arange(101) * -10.0])     # 1 km straight
    base = stations(line, 100)
    shifted = stations(line, 100, start_m=50)
    assert base[0]["river_km"] == 0 and shifted[0]["river_km"] == pytest.approx(0.05)
    assert [round(s["river_km"] - 0.05, 6) for s in shifted] == [s["river_km"] for s in base[:-1]]


def test_variant_provenance_keeps_fetch_and_frozen_tbsi_only():
    from floodforest.sensitivity import variant_provenance
    cfg = load_config()
    event = cfg.events.melamchi2021__tbsi_plus20
    base = {"event": "melamchi2021", "composites": {"s2_pre": {}}, "files": {"a": "b"},
            "label": {"t_bsi": -0.01, "other": 1}, "table": {"rows": 1}, "transects": {}}
    out = variant_provenance(base, "melamchi2021__tbsi_plus20", event)
    assert out["label"] == {"t_bsi": -0.01}
    assert "table" not in out and "transects" not in out
    assert out["variant"]["config_overrides"] == {"t_bsi_factor": pytest.approx(1.2)}
    with pytest.raises(RuntimeError, match="no frozen t_bsi"):
        variant_provenance({"composites": {}}, "x", event)


# --- runner helpers -----------------------------------------------------------

def _planted(beta: float, seed: int = 11) -> pd.DataFrame:
    """Synthetic table in transects.csv column names with a planted forest effect.

    Reach height and log1p distance both carry beta * F per bank; noise is shared
    within 1 km segments, so neighbouring rows are correlated.
    """
    rng = np.random.default_rng(seed)
    n_st = 240                                   # 24 km at 100 m, 8 blocks of 3 km
    st = np.repeat(np.arange(n_st), 2)
    side = np.tile(["left", "right"], n_st)
    km = st * 0.1
    seg = np.floor(km).astype(int)
    n = len(st)
    wall = (rng.random(n) > 0.25).astype(int)
    forest = rng.random(n)
    bank_slope = rng.random(n) * 30
    seg_noise = rng.normal(0, 0.3, seg.max() + 1)[seg]
    height = 3 + beta * forest + 0.02 * bank_slope + seg_noise + rng.normal(0, 0.5, n)
    logdist = 3 + 0.3 * beta * forest + 0.01 * bank_slope + seg_noise + rng.normal(0, 0.3, n)
    dist = np.expm1(np.clip(logdist, 0, None))
    return pd.DataFrame({
        "station_id": st, "side": side, "river_km": km, "segment_id": seg,
        "block_id": np.floor(km / 3).astype(int), "river_section": 0,
        "channel_resolved": (rng.random(n) > 0.2).astype(int),
        "forest_frac_100": forest, "channel_width_m": 20 + rng.random(n) * 30,
        "gradient_fab": rng.normal(size=n), "valley_width_m_fab": 50 + rng.random(n) * 300,
        "wall_slope_deg_fab": np.where(wall == 1, 30 + rng.random(n) * 20, np.nan),
        "valley_wall_found_fab": wall, "bank_slope_deg_fab": bank_slope,
        "log10_upa": rng.normal(2, 0.1, n), "outer_bend": rng.integers(0, 2, n),
        "curvature": rng.normal(size=n), "bank_elev_m_fab": 800 + rng.random(n) * 50,
        "reach_height_m_fab": np.clip(height, 0, None), "reach_distance_m_fab": dist,
        "overtopped_fab": (dist >= 20).astype(int), "censored_fab": 0,
    })


@pytest.mark.parametrize("beta", [-3.0, 3.0])
def test_planted_effect_recovered_in_sign_through_the_production_path(beta):
    """Test 4 through model_frame and fit_run (the code the stats step runs)."""
    from types import SimpleNamespace
    from floodforest.stats_run import fit_run
    cfg = SimpleNamespace(**vars(CFG))
    cfg.wild_bootstrap_draws = 399          # fewer draws for speed; same code path
    df = S.model_frame(_planted(beta), "forest_frac_100")
    run = fit_run(df, cfg, with_sections=False)
    for model in ("A", "B"):
        for y in ("y_height", "y_logdist"):
            r = run[model][y]
            assert "error" not in r, r
            assert np.sign(r["beta"]) == np.sign(beta), (model, y, r["beta"])
            assert r["p0"] <= cfg.interval_alpha, (model, y, r["p0"])
    lo, hi = run["A"]["y_height"]["interval"]
    assert lo < beta < hi
    assert run["b_info"]["stations"] == 240 and run["a_blocks"] == 8


def test_model_b_keeps_a_pair_whose_bank_lacks_only_a_model_a_control():
    from types import SimpleNamespace
    from floodforest.stats_run import fit_run
    cfg = SimpleNamespace(**vars(CFG))
    cfg.wild_bootstrap_draws = 99
    t = _planted(-3.0)
    t.loc[0, "curvature"] = np.nan          # station 0, left bank
    run = fit_run(S.model_frame(t, "forest_frac_100"), cfg, with_sections=False, trim=False)
    assert run["a_rows"] == len(t) - 1
    assert run["b_info"]["stations"] == 240
    assert 0 in set(run["_pairs"].station_id)


def test_logistic_is_not_estimable_with_fewer_events_than_terms():
    from floodforest.stats_run import _logit
    df = S.model_frame(_table(60, seed=5), "forest_frac_100")
    df["y_over"] = 0
    df.iloc[:3, df.columns.get_loc("y_over")] = 1
    res = _logit(df, "forest_frac + bank_slope_deg + gradient + curvature", CFG)
    assert res["error"].startswith("not estimable: 3 rows in the smaller outcome class")


def test_include_as_zero_and_rhs_without_constant_terms():
    from floodforest.stats_run import include_as_zero, rhs_for
    t = pd.DataFrame({"in_reach": [1, 1, 1], "analysis_set": [1, 0, 0],
                      "exclude_reasons": ["", "near_bank_nodata", "near_bank_nodata,low_valid_fraction"],
                      "reach_distance_m_fab": [10.0, np.nan, np.nan],
                      "reach_height_m_fab": [1.0, np.nan, np.nan],
                      "overtopped_fab": [0, np.nan, np.nan], "censored_fab": [0, np.nan, np.nan]})
    z = include_as_zero(t, "near_bank_nodata")
    assert list(z.analysis_set) == [1, 1, 0] and z.attrs["added_rows"] == 1
    assert z.reach_distance_m_fab.iloc[1] == 0 and t.analysis_set.iloc[1] == 0
    df = S.model_frame(_table(), "forest_frac_100").assign(channel_resolved=1)
    rhs = rhs_for(df, CFG, with_sections=True)
    assert "channel_resolved" not in rhs and "valley_wall_found" in rhs
    assert "C(river_section)" in rhs and rhs.startswith("forest_frac + ")
    assert "C(river_section)" not in rhs_for(df[df.river_section == 0], CFG, with_sections=True)
