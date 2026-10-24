"""Test 8: the E-value matches a hand-computed example."""

import math

from floodforest.stats import evalue_binary, evalue_continuous

# Worked by hand:
#   beta = -4, SD(y) = 8, so d = -0.5 and |d| = 0.5
#   RR* = exp(0.91 * 0.5) = exp(0.455) = 1.5762...
#   E = RR* + sqrt(RR* * (RR* - 1)) = 2.5291...
HAND = {
    "d": -0.5,
    "rr_star": math.exp(0.455),
    "evalue": 2.529142,
}


def test_evalue_matches_hand_computation():
    out = evalue_continuous(beta=-4.0, sd_y=8.0, ci_limit_nearest_null=-1.0)
    assert out["evalue"] == pytest_approx(HAND["evalue"])


def test_evalue_is_one_when_the_interval_includes_zero():
    out = evalue_continuous(beta=-4.0, sd_y=8.0, ci_limit_nearest_null=0.5)
    assert out["evalue_ci"] == 1.0


def test_binary_evalue_converts_odds_ratio_when_overtopping_is_common():
    # prevalence 20% > 15%, so RR = sqrt(4) = 2 and E = 2 + sqrt(2) = 3.4142
    assert evalue_binary(odds_ratio=4.0, prevalence=0.20) == pytest_approx(
        2 + math.sqrt(2))
    # prevalence 10% <= 15%, so RR = OR = 4 and E = 4 + sqrt(12) = 7.4641
    assert evalue_binary(odds_ratio=4.0, prevalence=0.10) == pytest_approx(
        4 + math.sqrt(12))


def pytest_approx(value, tol=1e-4):
    import pytest
    return pytest.approx(value, abs=tol)
