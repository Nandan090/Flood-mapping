"""The sample KML for Google Earth: stations named by station_id, lines by station_id_side."""

from pathlib import Path

import geopandas as gpd
from shapely.geometry import LineString, Point

from floodforest.validation_kml import export_sample_kml

STATION = (500000.0, 3100000.0)
LEFT = LineString([STATION, (500000.0, 3100200.0)])
RIGHT = LineString([STATION, (500000.0, 3099800.0)])


def _sample(folder: Path) -> Path:
    gpkg = folder / "sample_transects.gpkg"
    gpd.GeoDataFrame({"station_id": [7, 7], "side": ["left", "right"], "river_km": [5.2, 5.2]},
                     geometry=[LEFT, RIGHT], crs="EPSG:32645").to_file(
        gpkg, layer="transects", driver="GPKG")
    gpd.GeoDataFrame({"station_id": [7], "river_km": [5.2]}, geometry=[Point(STATION)],
                     crs="EPSG:32645").to_file(gpkg, layer="stations", driver="GPKG", mode="a")
    return gpkg


def test_kml_names_stations_by_id_and_lines_by_id_and_side(tmp_path):
    text = export_sample_kml(_sample(tmp_path), tmp_path / "s.kml", station_name="{}").read_text(
        encoding="utf-8")
    assert "<name>7</name>" in text and "station_7" not in text
    assert "<name>7_left</name>" in text and "<name>7_right</name>" in text
    # Blind: names and coordinates only.
    assert "ExtendedData" not in text and "forest" not in text.lower()
    assert "river_km" not in text and "5.2" not in text


def test_the_default_station_names_are_unchanged(tmp_path):
    text = export_sample_kml(_sample(tmp_path), tmp_path / "s.kml").read_text(encoding="utf-8")
    assert "<name>station_7</name>" in text
