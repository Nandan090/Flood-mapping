"""Tests 4 to 6 (section 11, with test 5 as redefined in 14.15), plus checks
of the pieces the primary interval is built from."""

import numpy as np
import pandas as pd
import pytest
import statsmodels.formula.api as smf

from floodforest.config import load_config
from floodforest.stats import (PERCENTILE_LABEL, WildClusterBootstrap,
                               block_ids, enough_blocks,
                               logistic_cluster_interval, model_a_rhs,
                               pairs_percentile_interval, wild_cluster_result)

CFG = load_config()
SEED = CFG.seed
RHS = "forest + bank_slope"

# Frozen by 14.15. Do not change these to make test 5 pass.
COVERAGE_REACH_KM = (24, 36, 72)
COVERAGE_ROWS_PER_KM = 20          # 10 stations per km, two banks each
COVERAGE_RHO = 0.8                 # between neighbouring stations
COVERAGE_SIMULATIONS = 400
COVERAGE_DRAWS = 399
COVERAGE_BAND = (0.92, 0.98)


def _ar1(n, rho, rng):
    """Stationary AR(1) series with unit variance."""
    out = np.empty(n)
    out[0] = rng.normal()
    innovation_sd = np.sqrt(1 - rho ** 2)
    for i in range(1, n):
        out[i] = rho * out[i - 1] + innovation_sd * rng.normal()
    return out


def simulate_reach(length_km, seed, beta=0.0, rho=COVERAGE_RHO,
                   confound=False):
    """Synthetic transect table along one reach.

    Exposure, one control (bank_slope) and the noise are each correlated
    between neighbouring stations with coefficient rho, separately on each
    bank. With confound=True, forest follows bank slope and only bank slope
    moves the outcome.
    """
    rng = np.random.default_rng(seed)
    stations = length_km * COVERAGE_ROWS_PER_KM // 2
    river_km = np.arange(stations) * 0.1
    blocks = block_ids(river_km, 0.0, CFG.bootstrap_block_km,
                       CFG.bootstrap_final_block_min_km, reach_end_km=length_km)
    frames = []
    for side in ("left", "right"):
        slope = _ar1(stations, rho, rng)
        noise = _ar1(stations, rho, rng)
        if confound:
            forest = 0.8 * slope + 0.6 * _ar1(stations, rho, rng)
            y = 1.0 * slope + noise
        else:
            forest = _ar1(stations, rho, rng)
            y = beta * forest + 1.0 * slope + noise
        frames.append(pd.DataFrame({
            "station": np.arange(stations), "side": side,
            "river_km": river_km, "block_id": blocks,
            "segment_id": np.floor(river_km).astype(int),
            "forest": forest, "bank_slope": slope, "reach_height": y,
        }))
    return pd.concat(frames, ignore_index=True)


def paired_banks(table):
    """Model B rows: left minus right at each station."""
    left = table[table.side == "left"].set_index("station")
    right = table[table.side == "right"].set_index("station")
    return pd.DataFrame({
        "block_id": left.block_id,
        "segment_id": left.segment_id,
        "d_forest": left.forest - right.forest,
        "d_bank_slope": left.bank_slope - right.bank_slope,
        "d_reach_height": left.reach_height - right.reach_height,
    }).reset_index()


# ---------------------------------------------------------------------------
# Fast checks
# ---------------------------------------------------------------------------

def test_blocks_are_three_km_and_a_short_final_block_merges():
    km = np.arange(0, 10.01, 0.1)
    ids = block_ids(km, 0.0, 3, 1.5, reach_end_km=10)
    # 9 to 10 km is 1 km, shorter than 1.5 km, so it joins the 6 to 9 km block.
    assert sorted(set(ids)) == [0, 1, 2]
    assert ids[km.searchsorted(9.5)] == 2

    km = np.arange(0, 10.51, 0.1)
    ids = block_ids(km, 0.0, 3, 1.5, reach_end_km=10.5)
    assert sorted(set(ids)) == [0, 1, 2, 3]   # 1.5 km final block is kept


def test_blocks_count_from_the_reach_start():
    km = np.array([5.0, 7.9, 8.0, 10.9, 11.0])
    assert block_ids(km, 5.0, 3, 1.5, reach_end_km=14).tolist() == [0, 0, 1, 1, 2]


def test_fewer_than_eight_blocks_is_flagged():
    assert not enough_blocks(7, CFG.min_blocks)
    assert enough_blocks(8, CFG.min_blocks)


def test_cluster_standard_error_matches_statsmodels_cr1():
    table = simulate_reach(24, seed=SEED)
    wb = WildClusterBootstrap(table, "reach_height", RHS, "forest",
                              "block_id", draws=9, seed=SEED)
    ref = smf.ols(f"reach_height ~ {RHS}", data=table).fit(
        cov_type="cluster", cov_kwds={"groups": table.block_id})
    assert wb.beta_hat == pytest.approx(ref.params["forest"], rel=1e-10)
    assert wb.se_cluster == pytest.approx(ref.bse["forest"], rel=1e-8)


def test_p_value_is_one_at_the_estimate_and_falls_away_from_it():
    table = simulate_reach(36, seed=SEED)
    wb = WildClusterBootstrap(table, "reach_height", RHS, "forest",
                              "block_id", draws=399, seed=SEED)
    assert wb.p_value(wb.beta_hat) == 1.0
    assert wb.p_value(wb.beta_hat + 10 * wb.se_cluster) <= 0.05


def test_interval_edges_sit_where_p_crosses_alpha():
    table = simulate_reach(36, seed=SEED)
    wb = WildClusterBootstrap(table, "reach_height", RHS, "forest",
                              "block_id", draws=399, seed=SEED)
    lo, hi = wb.interval(tolerance_se=CFG.ci_bisection_tolerance_se)
    step = 2 * CFG.ci_bisection_tolerance_se * wb.se_cluster
    assert lo < wb.beta_hat < hi
    assert wb.p_value(lo + step) > 0.05 and wb.p_value(lo - step) <= 0.05
    assert wb.p_value(hi - step) > 0.05 and wb.p_value(hi + step) <= 0.05


def test_sensitivity_interval_is_the_labelled_one_km_bootstrap():
    table = simulate_reach(24, seed=SEED)
    out = pairs_percentile_interval(table, "reach_height", RHS, "forest",
                                    reps=50, seed=SEED)
    assert out["label"] == PERCENTILE_LABEL
    assert "1 km" in out["method"]
    lo, hi = out["interval"]
    assert lo < hi


def test_logistic_interval_uses_t_with_blocks_minus_one_df():
    table = simulate_reach(36, seed=SEED)
    table["overtopped"] = (table.reach_height > 0).astype(int)
    out = logistic_cluster_interval(table, "overtopped", RHS, "forest",
                                    "block_id")
    assert out["df"] == table.block_id.nunique() - 1
    lo, hi = out["or_interval"]
    assert lo < out["odds_ratio"] < hi
    assert "approximate" in out["method"]


def test_model_a_adds_river_section_only_when_there_are_sections():
    assert "C(river_section)" not in model_a_rhs(spline_df=CFG.spline_df)
    rhs = model_a_rhs(spline_df=CFG.spline_df, with_sections=True)
    assert "C(river_section)" in rhs
    assert f"cr(river_km, df={CFG.spline_df})" in rhs


# ---------------------------------------------------------------------------
# Section 11 tests 4 to 6
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_planted_effect_is_recovered_by_models_a_and_b():
    beta = -0.5
    table = simulate_reach(36, seed=SEED, beta=beta)

    a = wild_cluster_result(table, "reach_height", RHS, "forest", "block_id",
                            draws=CFG.wild_bootstrap_draws, seed=SEED,
                            tolerance_se=CFG.ci_bisection_tolerance_se)
    lo, hi = a["interval"]
    assert lo < beta < hi
    assert a["excludes_zero"]

    b = wild_cluster_result(paired_banks(table), "d_reach_height",
                            "d_forest + d_bank_slope", "d_forest", "block_id",
                            draws=CFG.wild_bootstrap_draws, seed=SEED,
                            tolerance_se=CFG.ci_bisection_tolerance_se)
    lo, hi = b["interval"]
    assert lo < beta < hi
    assert b["excludes_zero"]


@pytest.mark.slow
def test_null_interval_covers_zero_at_nominal_rate():
    """Test 5 as frozen in 14.15.

    The primary interval is the set of b0 with p(b0) > 0.05, so it covers
    zero exactly when p(0) > 0.05.
    """
    coverage = {}
    for length_km in COVERAGE_REACH_KM:
        covered = 0
        for i in range(COVERAGE_SIMULATIONS):
            seed = SEED + 100_000 * length_km + i   # own seed per simulation
            table = simulate_reach(length_km, seed=seed)
            wb = WildClusterBootstrap(table, "reach_height", RHS, "forest",
                                      "block_id", draws=COVERAGE_DRAWS,
                                      seed=seed)
            covered += wb.p_value(0.0) > 0.05
        coverage[length_km] = covered / COVERAGE_SIMULATIONS

    report = ", ".join(f"{km} km: {rate:.1%}" for km, rate in coverage.items())
    print(f"\nCoverage of zero by the primary 95% interval: {report}")
    low, high = COVERAGE_BAND
    assert all(low <= rate <= high for rate in coverage.values()), report


@pytest.mark.slow
def test_bank_slope_control_removes_confounding():
    # Forest sits on steeper banks but has no effect. Without the slope
    # control the model invents one; with it, zero stays inside the interval.
    table = simulate_reach(36, seed=SEED, confound=True)
    controlled = WildClusterBootstrap(table, "reach_height", RHS, "forest",
                                      "block_id", draws=CFG.wild_bootstrap_draws,
                                      seed=SEED)
    assert controlled.p_value(0.0) > 0.05

    naive = WildClusterBootstrap(table, "reach_height", "forest", "forest",
                                 "block_id", draws=CFG.wild_bootstrap_draws,
                                 seed=SEED)
    assert naive.p_value(0.0) <= 0.05
