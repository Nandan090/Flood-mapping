"""KML export names each line station_id_side, and the converter snaps points."""

import json
from pathlib import Path

import geopandas as gpd
import pytest
from pyproj import Transformer
from shapely.geometry import LineString, Point

from floodforest.validation_kml import (
    convert_edge_gpkg,
    convert_edge_kml,
    export_sample_kml,
    side_from_channel,
    snap_to_station_transects,
)

LINE = LineString([(500000.0, 3100000.0), (500000.0, 3100200.0)])
TO_WGS = Transformer.from_crs(32645, 4326, always_xy=True)


def _sample(folder: Path) -> Path:
    gpkg = folder / "sample_transects.gpkg"
    lines = gpd.GeoDataFrame(
        {"station_id": [7], "side": ["left"], "river_km": [1.0]},
        geometry=[LINE], crs="EPSG:32645")
    points = gpd.GeoDataFrame(
        {"station_id": [7], "river_km": [1.0]},
        geometry=[Point(500000.0, 3100000.0)], crs="EPSG:32645")
    lines.to_file(gpkg, layer="transects", driver="GPKG")
    points.to_file(gpkg, layer="stations", driver="GPKG", mode="a")
    return gpkg


def _kml(folder: Path, easting: float, northing: float, description: str,
         name: str = "7_left") -> Path:
    lon, lat = TO_WGS.transform(easting, northing)
    path = folder / "marks.kml"
    path.write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2"><Document><Placemark>
<name>{name}</name>
<description><![CDATA[{description}]]></description>
<Point><coordinates>{lon},{lat},0</coordinates></Point>
</Placemark></Document></kml>
""",
        encoding="utf-8",
    )
    return path


def test_exported_lines_are_named_station_and_side(tmp_path):
    gpkg = _sample(tmp_path)
    kml = export_sample_kml(gpkg, tmp_path / "sample.kml")
    text = kml.read_text(encoding="utf-8")
    assert "<name>7_left</name>" in text
    assert "<name>station_7</name>" in text
    assert "<name>stations</name>" in text
    assert "<name>transects</name>" in text


def test_a_nearby_point_snaps_onto_the_line(tmp_path):
    gpkg = _sample(tmp_path)
    kml = _kml(tmp_path, 500010.0, 3100050.0,
               "image_date: 2021-10-15<br>confidence: 2<br>note: terrace edge")
    out = convert_edge_kml(kml, gpkg, tmp_path / "validation_edges.geojson")
    edges = gpd.read_file(out)
    props = json.loads(Path(out).read_text(encoding="utf-8"))["features"][0]["properties"]
    assert list(edges.columns) == ["station_id", "side", "image_date", "confidence", "note", "geometry"]
    assert edges.crs.to_string() == "EPSG:32645"
    assert props["station_id"] == 7
    assert props["side"] == "left"
    assert props["image_date"] == "2021-10-15"
    assert props["confidence"] == 2
    assert props["note"] == "terrace edge"
    snapped = edges.geometry.iloc[0]
    assert snapped.x == pytest.approx(500000.0, abs=0.01)
    assert snapped.y == pytest.approx(3100050.0, abs=0.01)


def test_a_point_more_than_20_m_away_is_refused(tmp_path):
    gpkg = _sample(tmp_path)
    kml = _kml(tmp_path, 500030.0, 3100050.0, "image_date: 2021-10-15\nconfidence: 1\nnote:")
    out = tmp_path / "validation_edges.geojson"
    with pytest.raises(RuntimeError, match="20 m"):
        convert_edge_kml(kml, gpkg, out)
    assert not out.exists()


def test_exactly_20_m_is_accepted(tmp_path):
    gpkg = _sample(tmp_path)
    kml = _kml(tmp_path, 500020.0, 3100050.0, "image_date: 2021-10-15\nconfidence: 3")
    out = convert_edge_kml(kml, gpkg, tmp_path / "validation_edges.geojson")
    props = json.loads(Path(out).read_text(encoding="utf-8"))["features"][0]["properties"]
    assert props["confidence"] == 3
    assert props["note"] == ""


def test_side_comes_from_the_re_centred_channel():
    assert side_from_channel(10.0, 50.0, recorded="left") == "right"
    assert side_from_channel(80.0, 50.0, recorded="right") == "left"
    assert side_from_channel(-10.0, float("nan"), recorded="left") == "right"


def test_a_gpkg_point_far_from_both_lines_is_refused(tmp_path):
    station = Point(500000.0, 3100000.0)
    left = LineString([(500000.0, 3100000.0), (500000.0, 3100200.0)])
    right = LineString([(500000.0, 3100000.0), (500000.0, 3099800.0)])
    lines = tmp_path / "sample_transects.gpkg"
    gpd.GeoDataFrame(
        {"station_id": [7, 7], "side": ["left", "right"], "river_km": [1.0, 1.0]},
        geometry=[left, right], crs="EPSG:32645",
    ).to_file(lines, layer="transects", driver="GPKG")
    gpd.GeoDataFrame(
        {"station_id": [7], "river_km": [1.0]},
        geometry=[station], crs="EPSG:32645",
    ).to_file(lines, layer="stations", driver="GPKG", mode="a")
    edges = tmp_path / "validation_edges.gpkg"
    gpd.GeoDataFrame(
        {"station_id": [7, 7], "side": ["left", "left"], "image_date": ["2022-01-17", "2022-01-17"],
         "wayback_release": ["2024-11-18", "2024-11-18"], "confidence": [1, 1], "note": ["", ""]},
        geometry=[Point(500000.0, 3100080.0), Point(500030.0, 3100080.0)],
        crs="EPSG:32645",
    ).to_file(edges, layer="validation_edges", driver="GPKG")
    out = tmp_path / "validation_edges.geojson"
    path, report = convert_edge_gpkg(edges, lines, out, {7: 120.0}, max_distance_m=20)
    props = json.loads(Path(path).read_text(encoding="utf-8"))["features"][0]["properties"]
    assert props["side"] == "right"
    assert props["station_id"] == 7
    assert report["n_written"] == 1
    assert len(report["refused"]) == 1
    assert report["mismatches"][0]["recorded"] == "left"
    assert report["mismatches"][0]["derived"] == "right"
    snapped, offset, gap, line = snap_to_station_transects(
        Point(500000.0, 3100080.0), left, right)
    assert line == "left" and offset == pytest.approx(80.0) and gap == pytest.approx(0.0)


def test_a_bad_name_or_confidence_is_refused(tmp_path):
    gpkg = _sample(tmp_path)
    out = tmp_path / "validation_edges.geojson"
    bad_name = _kml(tmp_path, 500005.0, 3100050.0, "image_date: 2021-10-15\nconfidence: 1",
                    name="7-left")
    with pytest.raises(RuntimeError, match="13_left"):
        convert_edge_kml(bad_name, gpkg, out)
    bad_confidence = _kml(tmp_path, 500005.0, 3100050.0, "image_date: 2021-10-15\nconfidence: sure")
    with pytest.raises(RuntimeError, match="confidence"):
        convert_edge_kml(bad_confidence, gpkg, out)
    assert not out.exists()
