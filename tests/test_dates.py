"""Test 7: exposure built from a post-event scene fails the run."""

import pytest

from floodforest.dates import (TimeSeparationError, assert_exposure_is_pre_event,
                               assert_outcome_is_post_event, composite_year,
                               in_window, window_bounds)
from floodforest.config import parse_utc


def test_post_event_exposure_raises():
    with pytest.raises(TimeSeparationError):
        assert_exposure_is_pre_event("2026-10-05T04:30:00Z",
                                     "2026-08-26T00:00:00Z",
                                     what="forest belt")


def test_pre_event_exposure_passes():
    assert_exposure_is_pre_event("2025-11-02T05:00:00Z", "2026-08-26T00:00:00Z")


def test_scene_before_the_post_window_raises():
    # 10:00 UTC is only 7 h after the 02:52 flood, inside the 12 h buffer.
    with pytest.raises(TimeSeparationError):
        assert_outcome_is_post_event("2026-08-26T10:00:00Z",
                                     "2026-08-26T14:52:00Z")


def test_scene_after_the_buffer_passes():
    assert_outcome_is_post_event("2026-08-26T15:10:00Z", "2026-08-26T14:52:00Z")


def test_bare_end_date_includes_that_calendar_day():
    start, end = window_bounds("2020-10-01", "2020-12-15")
    assert (start, end) == ("2020-10-01T00:00:00Z", "2020-12-16T00:00:00Z")
    assert in_window(parse_utc("2020-12-15T18:00:00Z"), start, end)
    assert not in_window(parse_utc("2020-12-16T00:00:00Z"), start, end)


def test_timestamp_window_is_used_as_written():
    assert window_bounds("2021-04-01T00:00:00Z", "2021-06-15T00:00:00Z") == (
        "2021-04-01T00:00:00Z", "2021-06-15T00:00:00Z")


def test_composite_year_rejects_a_window_that_spans_two_years():
    with pytest.raises(ValueError):
        composite_year("2020-11-01", "2021-01-15")
