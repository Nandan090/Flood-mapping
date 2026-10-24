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
    return cfg


def parse_utc(timestamp: str) -> datetime:
    """Parse an ISO-8601 timestamp ending in Z into an aware UTC datetime."""
    if not isinstance(timestamp, str) or not timestamp.endswith("Z"):
        raise ValueError(
            f"Expected a UTC timestamp like '2026-08-26T14:52:00Z', got {timestamp!r}"
        )
    return datetime.fromisoformat(timestamp.replace("Z", "+00:00")).astimezone(
        timezone.utc
    )
