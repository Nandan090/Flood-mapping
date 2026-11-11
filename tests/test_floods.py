"""Known floods around an event (14.17): validation, label effects, provenance, report."""

import json
from types import SimpleNamespace

import pytest

from floodforest.config import load_config
from floodforest.floods import flood_records, report_lines


def test_lende_floods_and_their_label_effects():
    cfg = load_config()
    rec = flood_records(cfg.events.lende2025, "lende2025")
    json.dumps(rec)  # provenance must be able to store them
    between, after = rec["known_floods_between"], rec["known_floods_after"]
    assert [r["start"] for r in between] == ["2024-08-25", "2024-09-27"]
    # Both are after the S1 baseline (ends 2024-08-15) and before the S2 pre
    # composite (starts 2024-10-01): radar label only.
    assert all(r["label_effect"].startswith("radar label only") for r in between)
    assert after[0]["label_effect"].startswith("inside")
    assert "Sentinel-1 post (2025-07-09 to 2025-08-15)" in after[0]["label_effect"]
    assert "Sentinel-2 post (2025-10-01 to 2025-12-15), whose scenes all come after it" in (
        after[0]["label_effect"])
    text = "\n".join(report_lines(rec))
    assert "affect only the radar label" in text
    assert "The 2025-07-30 flood at Timure and Bhotekoshi falls inside the post window" in text
    assert "Handling: none; see limitations." in text


def test_events_without_floods_add_nothing_to_the_report():
    cfg = load_config()
    for key in ("melamchi2021", "bhotekoshi2026", "lende2025_placebo"):
        rec = flood_records(getattr(cfg.events, key), key)
        assert rec == {"known_floods_between": [], "known_floods_after": []}
        assert report_lines(rec) == []


def _event(**over):
    w = lambda a, b: SimpleNamespace(start=a, end=b)  # noqa: E731
    ev = SimpleNamespace(
        s1_pre=w("2024-07-09T00:00:00Z", "2024-08-15T00:00:00Z"),
        s2_pre=w("2024-10-01", "2024-12-15"),
        s1_post=w("2025-07-09T00:00:00Z", "2025-08-15T00:00:00Z"),
        s2_post=w("2025-10-01", "2025-12-15"),
        known_floods_between=[])
    vars(ev).update(over)
    return ev


def _flood(start, end=None, **kw):
    return SimpleNamespace(start=start, end=end or start, place="p", type="t", note="n", **kw)


def test_effects_follow_the_windows():
    rec = flood_records(_event(known_floods_between=[
        _flood("2024-08-01"), _flood("2024-10-05")], known_floods_after=[_flood("2025-09-01")]), "e")
    b = [r["label_effect"] for r in rec["known_floods_between"]]
    assert b[0].startswith("overlaps the Sentinel-1 baseline")
    assert b[1].startswith("on or after the start of the Sentinel-2 pre composite")
    a = rec["known_floods_after"][0]["label_effect"]
    assert "Sentinel-1 post" not in a and "Sentinel-2 post" in a


def test_bad_flood_entries_are_refused():
    with pytest.raises(RuntimeError, match="no known_floods_between"):
        flood_records(SimpleNamespace(), "e")
    with pytest.raises(RuntimeError, match="ends before it starts"):
        flood_records(_event(known_floods_between=[_flood("2024-09-02", "2024-09-01")]), "e")
    with pytest.raises(RuntimeError, match="YYYY-MM-DD"):
        flood_records(_event(known_floods_between=[_flood("Sept 2024")]), "e")
    with pytest.raises(RuntimeError, match="missing"):
        flood_records(_event(known_floods_between=[SimpleNamespace(start="2024-09-01")]), "e")
