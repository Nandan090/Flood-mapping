"""Sensitivity rebuilds (8.8 station spacing and label thresholds, 14.22 station offset).

Each variant is an event derived in config.py: the real event with a few config
values changed. It reads the real event's downloaded rasters, keeps the real
event's frozen t_bsi, and writes only to outputs/<event>__<variant>.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

from .config import PROJECT_ROOT, load_config, rebuild_variants, variant_key

KEEP_FROM_FETCH = ("event", "crs", "grid_m", "corridor_buffer_m", "grid", "assets_verified",
                   "copernicus_dem_notice", "hansen", "worldcover_exported",
                   "worldcover_not_exported", "exclude_dates", "known_floods_between",
                   "composites", "files", "s1_geometry")


def variant_provenance(base: dict, key: str, event) -> dict:
    """The real event's fetch record plus its frozen t_bsi; nothing downstream of fetch."""
    if "label" not in base or "t_bsi" not in base["label"]:
        raise RuntimeError("The real event has no frozen t_bsi. Run its label step first.")
    out = {k: copy.deepcopy(base[k]) for k in KEEP_FROM_FETCH if k in base}
    out["event"] = key
    out["label"] = {"t_bsi": base["label"]["t_bsi"]}
    out["variant"] = {"base_event": event.base_event, "name": event.variant_name,
                      "description": event.variant_label,
                      "config_overrides": dict(event.cfg_overrides),
                      "rasters": f"data/{event.data_from}/rasters (read only)"}
    return out


def variant_complete(cfg, key: str) -> bool:
    return (PROJECT_ROOT / cfg.paths.outputs / key / "transects.csv").exists()


def build_variants(event_key: str, config_path: str | Path | None = None,
                   force: bool = False, only: list[str] | None = None) -> list[str]:
    """Run label, transects, profiles and table for each missing (or every) variant."""
    from .centerline import run_transects
    from .label import run_label
    from .profiles import run_profiles
    from .table import run_table

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    base_prov = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    if not base_prov.exists():
        raise FileNotFoundError(f"{base_prov} is missing. Run the real event first.")
    base = json.loads(base_prov.read_text(encoding="utf-8"))
    built = []
    for name in rebuild_variants(cfg):
        if only and name not in only:
            continue
        key = variant_key(event_key, name)
        if variant_complete(cfg, key) and not force:
            print(f"{key}: cached.")
            continue
        event = getattr(cfg.events, key)
        out_dir = PROJECT_ROOT / cfg.paths.outputs / key
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "provenance.json").write_text(
            json.dumps(variant_provenance(base, key, event), indent=2), encoding="utf-8")
        print(f"Building {key}: {event.variant_label}")
        for step in (run_label, run_transects, run_profiles, run_table):
            step(key, config_path)
        built.append(key)
    return built
