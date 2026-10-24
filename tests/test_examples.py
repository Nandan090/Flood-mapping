"""The six example cross-sections are chosen by a fixed rule (section 10, item 6)."""

import numpy as np
import pandas as pd

from floodforest.figures import choose_examples


def test_examples_follow_the_percentile_rule_and_never_repeat():
    heights = [0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, np.nan]
    n = len(heights)
    table = pd.DataFrame({
        "analysis_set": [1] * (n - 2) + [0, 1],
        "reach_height_m_fab": heights,
        "censored_fab": [0] * 12 + [1, 1, 0, 1, 0],
        "river_km": np.arange(n) * 0.1,
        "side": ["left", "right"] * 8 + ["left"],
    })
    chosen = choose_examples(table)
    index = [i for _, i in chosen]
    assert len(index) == 6 and len(set(index)) == 6
    # Row 15 (50 m) is outside the analysis set and row 16 has no height.
    assert 15 not in index and 16 not in index
    # The 5th and 25th percentiles are 0 m: the two lowest-km zero rows.
    assert index[:2] == [0, 1]
    # The sixth is a censored row that was not already used.
    assert table.loc[index[5], "censored_fab"] == 1
