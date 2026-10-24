"""The centerline QA rules: channel connectivity and off-channel runs."""

import numpy as np
import pytest

from floodforest.optical import mndwi, ndvi
from floodforest.qa import channel_touching_line, off_channel_stretches


def test_indices_match_the_spec_formulas():
    assert float(ndvi(np.array(0.5), np.array(0.1))) == pytest.approx(0.4 / 0.6)
    assert float(mndwi(np.array(0.2), np.array(0.05))) == pytest.approx(0.15 / 0.25)


def test_bare_pixels_away_from_the_line_are_not_the_channel():
    index = np.zeros((5, 5), dtype=bool)
    index[0, 0] = True  # bare ground with no connection to the line
    index[3, 2] = True  # water the line touches
    index[3, 3] = True
    on_line = np.zeros((5, 5), dtype=bool)
    on_line[3, 2] = True
    channel = channel_touching_line(index, on_line)
    assert channel[3, 2] and channel[3, 3]
    assert not channel[0, 0]


def test_a_run_of_exactly_50_m_is_not_listed():
    river_m = np.arange(0, 101, 10)  # 0, 10, ... 100
    inside = np.ones(len(river_m), dtype=bool)
    inside[1:6] = False  # 10 m through 50 m: length 40 m
    assert off_channel_stretches(river_m, inside, 50) == []
    inside[1:7] = False  # 10 m through 60 m: length 50 m, not longer
    assert off_channel_stretches(river_m, inside, 50) == []
    inside[1:8] = False  # 10 m through 70 m: length 60 m
    stretches = off_channel_stretches(river_m, inside, 50)
    assert stretches == [(1, 7, 60.0)]
