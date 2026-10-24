"""Time-separation checks (spec rule 0.4 and clarification 14.1).

Exposure and controls may only use scenes acquired before the event; outcomes
only scenes acquired after it. Boundaries are enforced on timestamps.
"""

from __future__ import annotations

from datetime import datetime, timedelta

from .config import parse_utc


class TimeSeparationError(RuntimeError):
    """A scene was used on the wrong side of the event."""


def assert_exposure_is_pre_event(scene_time_utc: str, pre_end_utc: str,
                                 what: str = "exposure") -> None:
    """Fail the run if an exposure or control scene is not strictly pre-event.

    Pre scenes must be acquired before 00:00 UTC on the event date.
    """
    scene = parse_utc(scene_time_utc)
    cutoff = parse_utc(pre_end_utc)
    if scene >= cutoff:
        raise TimeSeparationError(
            f"{what} scene acquired {scene_time_utc} is not before the "
            f"pre-event cutoff {pre_end_utc}. Exposure and controls must use "
            f"pre-event imagery only."
        )


def assert_outcome_is_post_event(scene_time_utc: str, post_start_utc: str,
                                 what: str = "outcome") -> None:
    """Fail the run if an outcome scene is earlier than the post window.

    Post scenes must be acquired later than event_time_utc + 12 h (or 00:00 UTC
    the day after, when the event hour is unknown).
    """
    scene = parse_utc(scene_time_utc)
    cutoff = parse_utc(post_start_utc)
    if scene < cutoff:
        raise TimeSeparationError(
            f"{what} scene acquired {scene_time_utc} is before the post-event "
            f"start {post_start_utc}. Outcomes must use post-event imagery only."
        )


def in_window(scene_time: datetime, start_utc: str, end_utc: str) -> bool:
    """True when start_utc <= scene_time < end_utc."""
    return parse_utc(start_utc) <= scene_time < parse_utc(end_utc)


def window_bounds(start: str, end: str) -> tuple[str, str]:
    """Turn a config window into an inclusive-start, exclusive-end UTC pair.

    A bare date (YYYY-MM-DD) is that calendar day. A bare end date is included,
    so the exclusive end is the next midnight. A full timestamp is used as
    written, which is how the Sentinel-1 windows are already stored.
    """
    return _bound(start, is_end=False), _bound(end, is_end=True)


def composite_year(start: str, end: str) -> int:
    """Calendar year of a composite window. Fails if the window spans two years."""
    years = {int(start[:4]), int(end[:4])}
    if len(years) != 1:
        raise ValueError(
            f"Window {start} to {end} spans more than one calendar year, so "
            "the Hansen year is ambiguous. No year was assumed."
        )
    return years.pop()


def _bound(value: str, is_end: bool) -> str:
    if isinstance(value, str) and value.endswith("Z"):
        parse_utc(value)
        return value
    if isinstance(value, str) and len(value) == 10:
        day = datetime.fromisoformat(value)
        if is_end:
            day = day + timedelta(days=1)
        return day.strftime("%Y-%m-%dT00:00:00Z")
    raise ValueError(
        f"Expected a UTC timestamp or a YYYY-MM-DD date, got {value!r}"
    )
