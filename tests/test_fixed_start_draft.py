"""Fixed-start OSM draft, landmark distances and confluence candidates (14.16)."""

import json
from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pytest
import yaml
from pyproj import Transformer
from shapely.geometry import LineString

from floodforest import centerline as C
from floodforest.config import PROJECT_ROOT

CRS = "EPSG:32645"
TO_LL = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)
TO_UTM = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
X0, Y0 = 640_000.0, 3_130_000.0


def _ll(line: LineString) -> LineString:
    return LineString([TO_LL.transform(x, y) for x, y in line.coords])


def _waterways():
    # Main river: 20 km north to south, drawn upstream-up so orientation is tested.
    main = LineString([(X0, Y0 - d) for d in range(0, 20_001, 500)])
    # Second river starts at the main river's end and runs 80 km south-west, in two
    # OSM pieces that join, one of them named in the local script.
    lower = [(X0 - 0.6 * d, Y0 - 20_000 - 0.8 * d) for d in range(0, 80_001, 500)]
    piece_a, piece_b = LineString(lower[:81]), LineString(lower[80:])
    trib = LineString([(X0 + 5_000, Y0 - 10_000), (X0, Y0 - 10_000)])
    return gpd.GeoDataFrame(
        {"name": ["Bhote Koshi", "Trishuli Ganga River", "त्रिशुली नदी", "Other Khola"],
         "name:en": [None, None, "Trisuli", None]},
        geometry=[_ll(g) for g in (LineString(list(main.coords)[::-1]), piece_a, piece_b, trib)],
        crs="EPSG:4326")


@pytest.fixture
def elevations_from_latitude(monkeypatch):
    # Ground falls to the south, so downstream is south.
    monkeypatch.setattr(C, "_elevations", lambda lon, lat: np.asarray(lat) * 1e4)


def _event(**over):
    start_lon, start_lat = TO_LL.transform(X0 + 40, Y0 - 2_000)
    mark_lon, mark_lat = TO_LL.transform(X0, Y0 - 16_000)
    event = SimpleNamespace(
        centerline_file="data/shared/centerline.geojson",
        start_point=SimpleNamespace(name="border", lon=start_lon, lat=start_lat,
                                    rule="fixed_point"),
        osm_river_names=["Bhote ?Koshi|भोटे कोशी", "Tri[s]h?uli|त्रिशुली|त्रिशूली"],
        rivers=["Bhotekoshi", "Trishuli"], osm_snap_m=150, draft_length_km=72,
        reach_start_km=0, reach_end_km=35,
        landmarks=[SimpleNamespace(name="town", lon=mark_lon, lat=mark_lat),
                   SimpleNamespace(name="unknown", lon=None, lat=None)])
    vars(event).update(over)
    return event


def _cfg(event):
    other = SimpleNamespace(centerline_file=event.centerline_file)
    return SimpleNamespace(crs=CRS, events=SimpleNamespace(ev=event, ev2026=other))


def _draft(tmp_path, event):
    out = tmp_path / "centerline.geojson"
    C._draft_from_fixed_start(_cfg(event), event, "ev", _waterways(), TO_UTM, TO_LL, out)
    return json.loads(out.read_text(encoding="utf-8"))


def test_fixed_start_draft_starts_at_the_point_and_runs_72_km(tmp_path, elevations_from_latitude,
                                                              capsys):
    data = _draft(tmp_path, _event())
    feats = {f["properties"]["role"]: f for f in data["features"]}
    assert set(feats) == {"centerline", "start", "confluence"}
    props = feats["centerline"]["properties"]
    assert props["verified_by_user"] is False and props["start_rule"] == "fixed_point"
    assert props["shared_by"] == ["ev", "ev2026"]
    assert props["length_m"] == pytest.approx(72_000, abs=1)
    # The line starts at the config point, 40 m off the OSM line, linked straight to
    # the OSM point 2 km below its upstream end.
    assert props["start_link_to_osm_m"] == pytest.approx(40, abs=1)
    line = LineString([TO_UTM.transform(*c) for c in feats["centerline"]["geometry"]["coordinates"]])
    assert line.coords[0] == pytest.approx((X0 + 40, Y0 - 2_000), abs=0.5)
    assert line.coords[1] == pytest.approx((X0, Y0 - 2_000), abs=0.5)
    assert line.coords[-1][1] < Y0 - 20_000  # runs downstream, onto the second river
    assert props["confluence_river_km"] == pytest.approx(18.04, abs=0.01)
    assert props["landmarks_on_draft"]["town"]["river_km"] == pytest.approx(14.04, abs=0.01)
    printed = capsys.readouterr().out
    assert "town: river_km 14.04" in printed and "DRAFT" in printed
    assert "Landmark unknown: no lon/lat" in printed


def test_fixed_start_refuses_a_point_off_the_line_and_a_short_line(tmp_path,
                                                                    elevations_from_latitude):
    lon, lat = TO_LL.transform(X0 + 400, Y0 - 2_000)
    far = SimpleNamespace(name="border", lon=lon, lat=lat, rule="fixed_point")
    with pytest.raises(RuntimeError, match="not joined"):
        _draft(tmp_path, _event(start_point=far, start_link_max_m=250))
    with pytest.raises(RuntimeError, match="Not inventing"):
        _draft(tmp_path, _event(draft_length_km=120))
    missing = SimpleNamespace(name="border", lon=None, lat=None, rule="fixed_point")
    with pytest.raises(RuntimeError, match="no lon/lat"):
        _draft(tmp_path, _event(start_point=missing))
    assert not (tmp_path / "centerline.geojson").exists()


def _write_merit(path, line: LineString, upa_along):
    import rasterio
    from rasterio.transform import from_origin
    xs, ys = zip(*line.coords)
    west, north = min(xs) - 500, max(ys) + 500
    width = int((max(xs) - west + 500) // 10) + 1
    height = int((north - min(ys) + 500) // 10) + 1
    data = np.zeros((height, width), dtype="float32")
    transform = from_origin(west, north, 10, 10)
    for d in np.arange(0, line.length + 1, 10):
        p = line.interpolate(d)
        r, c = rasterio.transform.rowcol(transform, p.x, p.y)
        data[r, c] = upa_along(d)
    with rasterio.open(path, "w", driver="GTiff", width=width, height=height, count=1,
                       dtype="float32", crs=CRS, transform=transform) as dst:
        dst.write(data, 1)
        dst.set_band_description(1, "upa")


def test_river_checks_prints_landmarks_and_the_50_percent_confluence(tmp_path, capsys):
    raw = yaml.safe_load((PROJECT_ROOT / "config.yaml").read_text(encoding="utf-8"))
    line = LineString([(X0, Y0 - d) for d in range(0, 10_001, 100)])
    lon, lat = TO_LL.transform(X0, Y0 - 3_000)
    centre = tmp_path / "centerline.geojson"
    centre.write_text(json.dumps({"type": "FeatureCollection", "features": [{
        "type": "Feature", "properties": {"role": "centerline", "verified_by_user": True},
        "geometry": {"type": "LineString", "coordinates": list(_ll(line).coords)}}]}),
        encoding="utf-8")
    ev = raw["events"]["lende2025"]
    ev["centerline_file"] = str(centre)
    ev["landmarks"] = [{"name": "bridge", "lon": lon, "lat": lat}]
    raw["paths"]["data"] = str(tmp_path / "data")
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(raw, allow_unicode=True), encoding="utf-8")

    out = C.river_checks("lende2025", path)
    printed = capsys.readouterr().out
    assert out["verified"] and out["landmarks"]["bridge"]["river_km"] == pytest.approx(3.0, abs=0.01)
    assert "on the verified centerline" in printed and "bridge: river_km 3.00" in printed
    assert "need MERIT upa" in printed and out["candidates"] == []

    rasters = tmp_path / "data" / "lende2025" / "rasters"
    rasters.mkdir(parents=True)
    # +30 % at 4 km (a small tributary, not a candidate); +80 % at 7 km (a candidate).
    _write_merit(rasters / "merit_upa.tif", line,
                 lambda d: 1000.0 * (1.3 if d >= 4000 else 1.0) * (1.8 if d >= 7000 else 1.0))
    out = C.river_checks("lende2025", path)
    printed = capsys.readouterr().out
    assert out["candidates"] == [pytest.approx(7.0, abs=0.11)]
    assert "Candidate confluences" in printed


def test_transects_refuse_a_draft_exclusions_file(tmp_path, monkeypatch):
    monkeypatch.setattr(C, "PROJECT_ROOT", tmp_path)
    folder = tmp_path / "data" / "ev"
    folder.mkdir(parents=True)
    gpd.GeoDataFrame({"verified_by_user": [False]},
                     geometry=[_ll(LineString([(X0, Y0), (X0 + 10, Y0)]))],
                     crs="EPSG:4326").to_file(folder / "exclusions.geojson", driver="GeoJSON")
    with pytest.raises(RuntimeError, match="is a draft"):
        C._load_exclusion_polygons("ev", CRS)
