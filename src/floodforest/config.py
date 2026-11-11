"""Load config.yaml and expose it as nested objects.

Nothing in this module hardcodes an analytical value; every number comes from
the file the user can edit.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config.yaml"


def _namespace(value):
    if isinstance(value, dict):
        return SimpleNamespace(**{k: _namespace(v) for k, v in value.items()})
    if isinstance(value, list):
        return [_namespace(v) for v in value]
    return value


def load_config(path: str | Path = DEFAULT_CONFIG_PATH) -> SimpleNamespace:
    """Read the YAML config. Raises FileNotFoundError naming the missing file."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open(encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    cfg = _namespace(raw)
    cfg.config_path = str(path)
    _add_placebo_events(cfg)
    return cfg


PLACEBO_SUFFIX = "_placebo"


def placebo_key(event_key: str) -> str:
    return f"{event_key}{PLACEBO_SUFFIX}"


def shift_years(stamp: str, years: int) -> str:
    """Move a 'YYYY-MM-DD' or 'YYYY-MM-DDTHH:MM:SSZ' stamp by whole years,
    keeping month, day and time (spec 14.2: same month-day offsets)."""
    if not isinstance(stamp, str) or len(stamp) < 10:
        raise ValueError(f"Expected a date or UTC timestamp, got {stamp!r}")
    year, rest = int(stamp[:4]), stamp[4:]
    if rest.startswith("-02-29"):
        raise ValueError(f"{stamp} is 29 February; it has no same-day shift.")
    return f"{year + years:04d}{rest}"


def _window(start: str, end: str) -> SimpleNamespace:
    return SimpleNamespace(start=start, end=end)


def derive_placebo_event(event: SimpleNamespace, event_key: str) -> SimpleNamespace:
    """The no-flood control of one event, as an event the pipeline can run (8.9, 14.2).

    Everything (centerline, reach, manual exclusions, section breaks) is the real
    event's. The image windows are the placebo block's, as written. The pre/post
    cutoffs and the sar_near sensitivity window are the real ones moved to the
    fake event date by whole years. Nothing here changes the real event.
    """
    block = event.placebo
    first = event.event_dates[0]
    if block.fake_event_date[4:] != first[4:]:
        raise ValueError(
            f"{event_key} placebo fake_event_date {block.fake_event_date} does not share "
            f"the month and day of the event date {first}; 14.2 shifts by whole years.")
    years = int(block.fake_event_date[:4]) - int(first[:4])
    out = SimpleNamespace(**{k: v for k, v in vars(event).items() if k != "placebo"})
    out.name = f"{event.name}: placebo (no-flood control, 14.2)"
    out.base_event = event_key
    out.is_placebo = True
    out.year_shift = years
    out.event_dates = [block.fake_event_date]
    out.event_time_utc = (shift_years(event.event_time_utc, years)
                          if getattr(event, "event_time_utc", None) else None)
    out.pre_end_utc = shift_years(event.pre_end_utc, years)
    out.post_start_utc = shift_years(event.post_start_utc, years)
    for name in ("s2_pre", "s2_post", "s1_pre", "s1_post"):
        win = getattr(block, name)
        setattr(out, name, _window(win.start, win.end))
    if hasattr(event, "sar_near"):
        out.sar_near = _window(shift_years(event.sar_near.start, years),
                               shift_years(event.sar_near.end, years))
    out.exclude_dates = []
    out.known_floods_between = []
    out.known_floods_after = []
    out.limitations = []  # the real event's caveats are not the placebo's
    out.flood_check = getattr(block, "flood_check", None)
    out.worldcover_standin = getattr(block, "worldcover_standin", None)
    if hasattr(out, "dynamic_world_gap_check"):
        del out.dynamic_world_gap_check
    return out


def _add_placebo_events(cfg) -> None:
    events = getattr(cfg, "events", None)
    if events is None:
        return
    for key, event in list(vars(events).items()):
        if hasattr(event, "placebo"):
            setattr(events, placebo_key(key), derive_placebo_event(event, key))
    for key, event in list(vars(events).items()):
        if not getattr(event, "is_placebo", False):
            for name, (label, overrides) in rebuild_variants(cfg).items():
                setattr(events, variant_key(key, name),
                        derive_variant_event(event, key, name, label, overrides))


VARIANT_SEP = "__"


def variant_key(event_key: str, name: str) -> str:
    return f"{event_key}{VARIANT_SEP}{name}"


def rebuild_variants(cfg) -> dict:
    """Sensitivity runs that need the label-to-table steps rebuilt (8.8, 14.22).

    name -> (description, config overrides). Values come from config.yaml.
    """
    f = float(cfg.label_threshold_sensitivity)
    z = float(cfg.sar_z_threshold)
    spacing = cfg.station_spacing_sensitivity_m
    offset = cfg.station_offset_sensitivity_m
    return {
        f"spacing{spacing:g}": (f"station spacing {spacing:g} m (8.8)",
                                {"station_spacing_m": spacing}),
        f"offset{offset:g}": (f"stations shifted {offset:g} m downstream (14.22)",
                              {"station_offset_m": offset}),
        "tbsi_minus20": (f"t_bsi x {1 - f:g} (8.8)", {"t_bsi_factor": 1 - f}),
        "tbsi_plus20": (f"t_bsi x {1 + f:g} (8.8)", {"t_bsi_factor": 1 + f}),
        "sarz_minus20": (f"SAR z threshold {round(z * (1 - f), 6):g} (8.8)",
                         {"sar_z_threshold": round(z * (1 - f), 6)}),
        "sarz_plus20": (f"SAR z threshold {round(z * (1 + f), 6):g} (8.8)",
                        {"sar_z_threshold": round(z * (1 + f), 6)}),
    }


def derive_variant_event(event: SimpleNamespace, event_key: str, name: str, label: str,
                         overrides: dict) -> SimpleNamespace:
    """A sensitivity rebuild of one event: same event, a few config values changed,
    the real event's downloaded rasters, its own outputs folder."""
    out = SimpleNamespace(**{k: v for k, v in vars(event).items() if k != "placebo"})
    out.base_event = event_key
    out.data_from = event_key
    out.is_variant = True
    out.variant_name = name
    out.variant_label = label
    out.cfg_overrides = dict(overrides)
    return out


def effective_config(cfg, event_key: str):
    """cfg with the event's sensitivity overrides applied (a copy; cfg is unchanged)."""
    event = getattr(cfg.events, event_key)
    overrides = getattr(event, "cfg_overrides", None)
    if not overrides:
        return cfg
    out = SimpleNamespace(**vars(cfg))
    for name, value in overrides.items():
        if not hasattr(cfg, name):
            raise KeyError(f"Override {name} is not a config.yaml setting.")
        setattr(out, name, value)
    return out


def rasters_dir(cfg, event_key: str) -> Path:
    """Downloaded rasters of an event; a rebuild variant reads its real event's."""
    event = getattr(cfg.events, event_key)
    return PROJECT_ROOT / cfg.paths.data / getattr(event, "data_from", event_key) / "rasters"


def parse_utc(timestamp: str) -> datetime:
    """Parse an ISO-8601 timestamp ending in Z into an aware UTC datetime."""
    if not isinstance(timestamp, str) or not timestamp.endswith("Z"):
        raise ValueError(
            f"Expected a UTC timestamp like '2026-08-26T14:52:00Z', got {timestamp!r}"
        )
    return datetime.fromisoformat(timestamp.replace("Z", "+00:00")).astimezone(
        timezone.utc
    )
