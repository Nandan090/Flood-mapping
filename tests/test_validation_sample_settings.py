"""Per-event validation sample settings (lende2025), the no-redraw guard and the instructions."""

from types import SimpleNamespace

import pandas as pd
import pytest

from floodforest.config import load_config
from floodforest.validate import (dated_digitising_instructions, digitising_instructions,
                                  in_reaches, refuse_redraw, sample_settings)


def test_lende_sample_settings_as_decided():
    cfg = load_config()
    s = sample_settings(cfg.events.lende2025, cfg)
    assert (s["n"], s["seed"], s["reaches"]) == (50, 42, [(5.0, 10.0), (25.0, 35.0)])


def test_melamchi_keeps_the_global_settings():
    cfg = load_config()
    s = sample_settings(cfg.events.melamchi2021, cfg)
    assert (s["n"], s["seed"], s["reaches"], s["block"]) == (40, int(cfg.seed), None, None)


def test_an_event_sample_needs_a_recorded_seed():
    cfg = load_config()
    event = SimpleNamespace(validation=SimpleNamespace(sample_size=50, sample_seed=None))
    with pytest.raises(RuntimeError, match="sample_seed"):
        sample_settings(event, cfg)


def test_reaches_include_both_ends_and_nothing_between():
    km = pd.Series([4.99, 5.0, 7.5, 10.0, 10.01, 24.9, 25.0, 35.0, 35.1])
    assert in_reaches(km, [(5, 10), (25, 35)]).tolist() == [
        False, True, True, True, False, False, True, True, False]


def test_the_sample_is_never_redrawn(tmp_path):
    refuse_redraw("e", tmp_path, {}, None)          # nothing drawn yet: passes
    with pytest.raises(RuntimeError, match="never redrawn"):
        refuse_redraw("e", tmp_path, {"validation_sample": {}}, None)
    with pytest.raises(RuntimeError, match="sample_drawn_on"):
        refuse_redraw("e", tmp_path, {}, "2026-11-12")
    (tmp_path / "sample_stations.geojson").write_text("{}", encoding="utf-8")
    with pytest.raises(RuntimeError, match="sample_stations.geojson"):
        refuse_redraw("e", tmp_path, {}, None)


def test_lende_instructions_say_what_the_user_decided():
    cfg = load_config()
    block = cfg.events.lende2025.validation
    text = dated_digitising_instructions(
        "Lende", "8 July 2025", cfg.crs, 50, block.imagery_source, block.image_date_window,
        block.digitising_notes)
    assert "at the 50 places" in text
    assert ("Use Google Earth Pro historical imagery captured between 9 July 2025 and "
            "31 May 2026.") in text
    assert "Prefer the earliest clear image after the flood" in text
    assert "Record the capture date for every point." in text
    assert "If shadow, cloud or trees hide the edge, skip that bank." in text
    assert "Winter shadow will hide many right banks at km 25-35" in text
    for field in ("`station_id`", "`side`", "`image_date`", "`confidence`", "`note`"):
        assert field in text
    assert "forest" not in text.lower()


def test_melamchi_instructions_unchanged():
    text = digitising_instructions("X", "1 June 2021", "EPSG:32645")
    assert "at the 40 places" in text and "image_date" not in text
