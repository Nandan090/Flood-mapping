"""The CSV of edge points: read, converted to validation_edges.gpkg, snapped as before."""

from pathlib import Path

import geopandas as gpd
import pytest
from pyproj import Transformer
from shapely.geometry import LineString, Point

from floodforest.config import load_config
from floodforest.validate import edge_date_window
from floodforest.validation_kml import (EDGE_CSV_COLUMNS, convert_edge_gpkg, csv_to_edge_gpkg,
                                        read_edge_csv)

TO_WGS = Transformer.from_crs(32645, 4326, always_xy=True)
STATION = (500000.0, 3100000.0)
LEFT = LineString([STATION, (500000.0, 3100200.0)])     # left bank runs north
RIGHT = LineString([STATION, (500000.0, 3099800.0)])    # right bank runs south
HEADER = ",".join(EDGE_CSV_COLUMNS)


def _sample(folder: Path) -> Path:
    gpkg = folder / "sample_transects.gpkg"
    gpd.GeoDataFrame({"station_id": [7, 7], "side": ["left", "right"], "river_km": [5.2, 5.2]},
                     geometry=[LEFT, RIGHT], crs="EPSG:32645").to_file(
        gpkg, layer="transects", driver="GPKG")
    gpd.GeoDataFrame({"station_id": [7], "river_km": [5.2]}, geometry=[Point(STATION)],
                     crs="EPSG:32645").to_file(gpkg, layer="stations", driver="GPKG", mode="a")
    return gpkg


def _latlon(easting, northing):
    lon, lat = TO_WGS.transform(easting, northing)
    return f"{lat:.8f},{lon:.8f}"


def _csv(folder: Path, *lines: str) -> Path:
    path = folder / "validation_edges.csv"
    path.write_text("\n".join([HEADER, *lines]) + "\n", encoding="utf-8")
    return path


def test_a_good_csv_is_read_as_wgs84_points(tmp_path):
    path = _csv(tmp_path,
                f"7,left,2025-10-03,1,,{_latlon(500000.0, 3100050.0)}",
                f"7,Right,2025-10-03,2,flood stayed in channel,{_latlon(500000.0, 3099990.0)}")
    points = read_edge_csv(path)
    assert points.crs.to_string() == "EPSG:4326"
    assert points.side.tolist() == ["left", "right"]
    assert points.confidence.tolist() == [1, 2]
    assert points.note.tolist() == ["", "flood stayed in channel"]


@pytest.mark.parametrize("line, message", [
    ("7,up,2025-10-03,1,,27.9,85.3", "must be left or right"),
    ("7,left,2025-10-03,4,,27.9,85.3", "confidence"),
    ("7,left,,1,,27.9,85.3", "missing image_date"),
    ("7,left,15/10/2025,1,,27.9,85.3", "Use a date like"),   # day/month is not month/day
    ("7,left,Oct 3 2025,1,,27.9,85.3", "Use a date like"),
    ("7,left,2025-10-03,1,,north,85.3", "decimal degrees"),
    ("seven,left,2025-10-03,1,,27.9,85.3", "not a whole number"),
])
def test_a_bad_csv_line_refuses_the_file(tmp_path, line, message):
    with pytest.raises(RuntimeError, match=message):
        read_edge_csv(_csv(tmp_path, line))


def test_wrong_columns_and_duplicates_are_refused(tmp_path):
    bad = tmp_path / "validation_edges.csv"
    bad.write_text("station_id,side,image_date,confidence,note,lat,lon\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="exactly the columns"):
        read_edge_csv(bad)
    dup = _csv(tmp_path, "7,left,2025-10-03,1,,27.9,85.3", "7,left,2025-10-04,1,,27.9,85.3")
    with pytest.raises(RuntimeError, match="more than once"):
        read_edge_csv(dup)


def _convert(tmp_path, csv):
    gpkg = csv_to_edge_gpkg(csv, tmp_path / "validation_edges.gpkg", "EPSG:32645")
    return convert_edge_gpkg(gpkg, _sample(tmp_path), tmp_path / "validation_edges.geojson",
                             {7: 0.0}, date_window=("2025-07-09", "2026-05-31"))


def test_csv_points_snap_to_their_transect(tmp_path):
    csv = _csv(tmp_path,
               f"7,left,2025-10-03,1,,{_latlon(500005.0, 3100050.0)}",      # 5 m off the line
               f"7,right,2025-10-03,1,,{_latlon(499992.0, 3099960.0)}")
    out, report = _convert(tmp_path, csv)
    edges = gpd.read_file(out)
    assert report["n_written"] == 2 and report["refused"] == []
    left = edges[edges.side == "left"].geometry.iloc[0]
    assert left.x == pytest.approx(500000.0, abs=0.01)
    assert left.y == pytest.approx(3100050.0, abs=0.05)


def test_a_csv_point_more_than_20_m_away_refuses_everything(tmp_path):
    csv = _csv(tmp_path,
               f"7,left,2025-10-03,1,,{_latlon(500005.0, 3100050.0)}",
               f"7,right,2025-10-03,1,,{_latlon(500025.0, 3099960.0)}")   # 25 m off
    with pytest.raises(RuntimeError, match="no GeoJSON was written"):
        _convert(tmp_path, csv)
    assert not (tmp_path / "validation_edges.geojson").exists()


def test_points_outside_the_event_window_are_left_out_and_counted(tmp_path):
    csv = _csv(tmp_path,
               f"7,left,2025-07-08,1,,{_latlon(500000.0, 3100050.0)}",     # the flood day
               f"7,right,2025-07-09,1,,{_latlon(500000.0, 3099960.0)}")
    _, report = _convert(tmp_path, csv)
    assert report["n_written"] == 1 and report["outside_date_window"] == 1


def test_a_qgis_gpkg_is_never_overwritten_by_the_csv(tmp_path):
    gpd.GeoDataFrame({"station_id": [7], "side": ["left"]}, geometry=[Point(STATION)],
                     crs="EPSG:32645").to_file(tmp_path / "validation_edges.gpkg",
                                               layer="validation_edges", driver="GPKG")
    csv = _csv(tmp_path, f"7,left,2025-10-03,1,,{_latlon(500000.0, 3100050.0)}")
    with pytest.raises(RuntimeError, match="keep only one"):
        csv_to_edge_gpkg(csv, tmp_path / "validation_edges.gpkg", "EPSG:32645")


def test_the_csv_gpkg_is_rewritten_from_an_edited_csv(tmp_path):
    gpkg = tmp_path / "validation_edges.gpkg"
    csv_to_edge_gpkg(_csv(tmp_path, f"7,left,2025-10-03,1,,{_latlon(500000.0, 3100050.0)}"),
                     gpkg, "EPSG:32645")
    csv_to_edge_gpkg(_csv(tmp_path, f"7,left,2025-10-03,1,,{_latlon(500000.0, 3100050.0)}",
                          f"7,right,2025-10-03,1,,{_latlon(500000.0, 3099960.0)}"),
                     gpkg, "EPSG:32645")
    assert len(gpd.read_file(gpkg, layer="validation_edges")) == 2


def test_each_event_uses_its_own_capture_date_window():
    events = load_config().events
    assert edge_date_window(events.lende2025) == {"date_window": ("2025-07-09", "2026-05-31")}
    assert edge_date_window(events.melamchi2021) == {}    # converter default, 2021-08-01..2022-05-31

def test_month_day_year_dates_are_read_and_written_as_iso(tmp_path):
    path = _csv(tmp_path,
                f"7,left,1/2/2026,1,,{_latlon(500000.0, 3100050.0)}",
                f"7,right,11/29/2025,1,,{_latlon(500000.0, 3099960.0)}")
    assert read_edge_csv(path).image_date.tolist() == ["2026-01-02", "2025-11-29"]
    out, _ = _convert(tmp_path, path)
    assert sorted(gpd.read_file(out).image_date.astype(str)) == ["2025-11-29", "2026-01-02"]


def test_each_snap_distance_is_reported(tmp_path):
    csv = _csv(tmp_path,
               f"7,left,2025-10-03,1,,{_latlon(500005.0, 3100050.0)}",
               f"7,right,2025-10-03,1,,{_latlon(499988.0, 3099960.0)}")
    _, report = _convert(tmp_path, csv)
    snaps = {s["side"]: s["snap_m"] for s in report["snap_distances"]}
    assert snaps["left"] == pytest.approx(5.0, abs=0.05)
    assert snaps["right"] == pytest.approx(12.0, abs=0.05)

def test_a_rerun_records_the_corrections_and_keeps_the_earlier_run():
    from types import SimpleNamespace

    from floodforest.validate import edge_rerun_record

    fix = [SimpleNamespace(date="2026-11-13", station_id=66, change="re-placed", banks=["left"],
                           reason="r"),
           SimpleNamespace(date="2026-11-13", station_id=263, change="removed", reason="shadow")]
    event = SimpleNamespace(validation=SimpleNamespace(edge_corrections=fix))
    first = {"n_written": 68}
    record = edge_rerun_record(event, first)
    assert record["rerun_after_corrections"]["corrected_on"] == "2026-11-13"
    assert [c["station_id"] for c in record["rerun_after_corrections"]["corrections"]] == [66, 263]
    assert record["earlier_runs"] == [first]
    again = edge_rerun_record(event, {**record, "n_written": 63})
    assert [r["n_written"] for r in again["earlier_runs"]] == [68, 63]
    plain = edge_rerun_record(SimpleNamespace(validation=SimpleNamespace()), None)
    assert set(plain) == {"run_on"}



def test_lende_edge_corrections_as_given():
    fixes = load_config().events.lende2025.validation.edge_corrections
    assert {(c.station_id, c.change) for c in fixes} == {
        (66, "re-placed"), (69, "re-placed"), (263, "removed"), (267, "removed"),
        (338, "removed"), (335, "removed")}
