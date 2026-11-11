"""Optional event limitations (printed in the stats report) and the lende2025 exclusions."""

from types import SimpleNamespace

import pytest

from floodforest.config import load_config
from floodforest.label import manual_exclusions
from floodforest.stats_run import event_limitations


def test_events_without_limitations_add_nothing():
    events = load_config().events
    for key in ("melamchi2021", "melamchi2021_placebo", "bhotekoshi2026"):
        assert event_limitations(getattr(events, key), key) == []
    assert event_limitations(SimpleNamespace(), "e") == []
    assert event_limitations(SimpleNamespace(limitations=[]), "e") == []


def test_limitations_become_a_section():
    lines = event_limitations(SimpleNamespace(limitations=[" First caveat. ", "Second."]), "e")
    assert lines == ["", "## Limitations", "", "- First caveat.", "- Second."]


def test_bad_limitations_are_refused():
    with pytest.raises(RuntimeError, match="empty or non-text"):
        event_limitations(SimpleNamespace(limitations=["ok", ""]), "e")
    with pytest.raises(RuntimeError, match="empty or non-text"):
        event_limitations(SimpleNamespace(limitations=[3]), "e")


def test_lende_limitations_stay_with_the_real_event():
    events = load_config().events
    lines = event_limitations(events.lende2025, "lende2025")
    assert lines[:3] == ["", "## Limitations", ""]
    assert len(lines) == 3 + 3
    assert "30 July 2025 Timure/Bhote Koshi flash flood" in lines[3]
    assert "km 15.0" in lines[4]
    assert "coded reservoir_or_dam in the outcome table" in lines[5]
    # The placebo is a different period; it does not inherit the real event's caveats.
    assert event_limitations(events.lende2025_placebo, "lende2025_placebo") == []


def test_lende_manual_exclusions_as_decided():
    cfg = load_config()
    event = cfg.events.lende2025
    assert manual_exclusions(event, "lende2025") == [
        {"from_km": 20.17, "to_km": 20.40, "side": "left", "reason": "tributary_fan"},
        {"from_km": 24.51, "to_km": 24.72, "side": "right", "reason": "hillslope_failure"},
        {"from_km": 30.81, "to_km": 31.26, "side": "right", "reason": "tributary_fan"},
    ]
    assert event.manual_exclusions_confirmed_on == "2026-11-11"
    assert all("contact-sheet review" in item.note for item in event.manual_exclusions)
    # The placebo uses the real event's exclusions (same reach and stations).
    assert manual_exclusions(cfg.events.lende2025_placebo, "p") == manual_exclusions(
        event, "lende2025")


def test_lende_exclusion_polygons_are_all_loaded():
    """The loader takes every verified polygon, whatever its reason text (14.19 e decision)."""
    import json

    from floodforest.centerline import _load_exclusion_polygons
    from floodforest.config import PROJECT_ROOT
    path = PROJECT_ROOT / "data" / "lende2025" / "exclusions.geojson"
    if not path.exists():
        pytest.skip("data/lende2025/exclusions.geojson is not present")
    features = json.loads(path.read_text(encoding="utf-8"))["features"]
    assert all(f["properties"]["verified_by_user"] is True for f in features)
    names = [f["properties"]["name"] or "" for f in features]
    for area in ("Area a:", "Area b:", "Area c:", "Area e_NW:", "Area h:"):
        assert sum(n.startswith(area) for n in names) == 1
    assert len(_load_exclusion_polygons("lende2025", "EPSG:32645")) == len(features) == 12
