"""Known floods around an event (14.17): recorded in provenance and the report.

known_floods_between: floods between the Sentinel-1 baseline and the event, whose
footprint the radar label would attribute to the event.
known_floods_after: floods inside the post window, whose footprint any post-event
label can contain.

Which labels a flood can reach follows from the event's own windows, never from
typed text.
"""

from __future__ import annotations

from datetime import date
from types import SimpleNamespace

from .dates import window_bounds

FIELDS = ("start", "end", "place", "type", "note")


def _records(event, key: str, event_key: str, required: bool) -> list[dict]:
    if not hasattr(event, key):
        if required:
            raise RuntimeError(
                f"{event_key} has no {key} in config.yaml. Add a list (use [] when none "
                "are known). None were assumed.")
        return []
    out = []
    for item in getattr(event, key) or []:
        rec = vars(item) if isinstance(item, SimpleNamespace) else dict(item)
        missing = [f for f in FIELDS if not rec.get(f)]
        if missing:
            raise RuntimeError(f"{event_key} {key} entry {rec} is missing {missing}.")
        try:
            start, end = date.fromisoformat(rec["start"]), date.fromisoformat(rec["end"])
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"{event_key} {key} dates must be YYYY-MM-DD: {rec}") from exc
        if end < start:
            raise RuntimeError(f"{event_key} {key} entry ends before it starts: {rec}")
        out.append(dict(rec))
    return out


def _day(stamp: str) -> date:
    return date.fromisoformat(stamp[:10])


def label_effect_between(flood: dict, event) -> str:
    """Which labels a flood between the S1 baseline and the event can reach."""
    s1_end = _day(window_bounds(event.s1_pre.start, event.s1_pre.end)[1])
    s2_start = _day(event.s2_pre.start)
    start, end = date.fromisoformat(flood["start"]), date.fromisoformat(flood["end"])
    if start >= s1_end and end < s2_start:
        return (f"radar label only: after the Sentinel-1 baseline ({event.s1_pre.start[:10]} to "
                f"{event.s1_pre.end[:10]}) and before the Sentinel-2 pre composite "
                f"({event.s2_pre.start} to {event.s2_pre.end}), so its scars are part of the "
                "optical pre-event state")
    if start < s1_end:
        return ("overlaps the Sentinel-1 baseline window; it can affect the radar baseline "
                "itself")
    return ("on or after the start of the Sentinel-2 pre composite; it can affect the "
            "optical pre-event state and the radar label")


def label_effect_after(flood: dict, event) -> str:
    """Which post windows a flood after the event falls in."""
    start = date.fromisoformat(flood["start"])
    hits = []
    for name, win in (("Sentinel-1 post", event.s1_post), ("Sentinel-2 post", event.s2_post)):
        lo, hi = (_day(v) for v in window_bounds(win.start, win.end))
        if start < hi:
            hits.append(f"{name} ({win.start[:10]} to {win.end[:10]})"
                        + ("" if start >= lo else ", whose scenes all come after it"))
    if not hits:
        return "after every post window; not in the label"
    return "inside or before the end of: " + "; ".join(hits)


def flood_records(event, event_key: str) -> dict:
    """Validated flood lists with their label effects, for provenance and the report."""
    between = _records(event, "known_floods_between", event_key, required=True)
    after = _records(event, "known_floods_after", event_key, required=False)
    for rec in between:
        rec["label_effect"] = label_effect_between(rec, event)
    for rec in after:
        rec["label_effect"] = label_effect_after(rec, event)
    return {"known_floods_between": between, "known_floods_after": after}


def report_lines(records: dict) -> list[str]:
    """Report section; empty when no flood is recorded."""
    between, after = records["known_floods_between"], records["known_floods_after"]
    if not between and not after:
        return []
    lines = ["## Known floods around the event (14.17)", "",
             "Supplied by the user in config.yaml. The effect column follows from the "
             "event's own image windows.", ""]
    for title, recs in (("Between the Sentinel-1 baseline and the event", between),
                        ("After the event, inside the post window", after)):
        if not recs:
            continue
        lines += [f"**{title}**", "", "| Dates | Place | Type | Note | Effect on the labels | Handling |",
                  "|---|---|---|---|---|---|"]
        for r in recs:
            dates = r["start"] if r["start"] == r["end"] else f"{r['start']} to {r['end']}"
            lines.append(f"| {dates} | {r['place']} | {r['type']} | {r['note']} | "
                         f"{r['label_effect']} | {r.get('handling', '')} |")
        lines.append("")
    radar_only = [r for r in between if r["label_effect"].startswith("radar label only")]
    if radar_only:
        lines.append("The floods between the baseline and the event ("
                     + ", ".join(r["start"] if r["start"] == r["end"] else
                                 f"{r['start']} to {r['end']}" for r in radar_only)
                     + ") affect only the radar label: their footprint can appear as radar "
                       "change and be attributed to the event. The optical label is not "
                       "affected, because its pre-event composite already shows their scars.")
    for r in after:
        if r["label_effect"].startswith("inside"):
            lines.append(f"The {r['start']} flood at {r['place']} falls inside the post window, "
                         "so the post-event labels can contain its footprint as well as the "
                         f"event's. Handling: {r.get('handling') or 'not set in config'}.")
    lines.append("")
    return lines
