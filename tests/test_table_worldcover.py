"""The WorldCover heading in the table log (14.23 a): Melamchi wording unchanged, lende gaps."""

from types import SimpleNamespace

import pandas as pd

from floodforest.config import load_config
from floodforest.table import worldcover_heading

OLD_MELAMCHI = ("River km using WorldCover (no Dynamic World image for Sentinel-2 tile T45RUM "
                "in the pre window):")


def _df(sources, km):
    return pd.DataFrame({"forest_source_100": sources, "river_km": km})


def test_melamchi_keeps_its_wording_byte_for_byte():
    event = load_config().events.melamchi2021
    df = _df(["worldcover", "dynamic_world"], [1.0, 2.0])
    assert worldcover_heading(event, df, "forest_source_100", 2020) == OLD_MELAMCHI


def test_an_event_without_a_tile_gap_reports_counts_and_km():
    event = SimpleNamespace()
    df = _df(["dynamic_world", "worldcover", "mixed", "mixed", "dynamic_world"],
             [1.0, 26.21, 26.81, 27.21, 30.0])
    text = worldcover_heading(event, df, "forest_source_100", 2021)
    assert text == ("River km using WorldCover: Dynamic World had small gaps near the line in "
                    "the pre window, and WorldCover 2021 stood in there (1 rows fully and 2 rows "
                    "partly, of 5; km 26.21–27.21):")
    assert "T45RUM" not in text


def test_no_worldcover_rows():
    df = _df(["dynamic_world"], [1.0])
    assert worldcover_heading(SimpleNamespace(), df, "forest_source_100", 2021).endswith(
        "none (Dynamic World covered every belt pixel):")


def test_lende_uses_the_gap_wording():
    event = load_config().events.lende2025
    assert not hasattr(event, "dynamic_world_gap_check")
    df = _df(["worldcover", "mixed"], [26.21, 27.21])
    assert worldcover_heading(event, df, "forest_source_100", 2021).startswith(
        "River km using WorldCover: Dynamic World had small gaps")
