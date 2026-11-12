"""Google Earth KML for the validation sample, and the path back.

The sample lines are named station_id_side. A placemark KML with those
names comes back as validation_edges.geojson. A point more than 20 m from
its transect line is refused, and nothing is written.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from datetime import datetime
from html import unescape
from pathlib import Path

SNAP_LIMIT_M = 20.0
VALIDATION_DATE_START = "2021-08-01"
VALIDATION_DATE_END = "2022-05-31"
_EDGE_NAME = re.compile(r"^(\d+)_(left|right)$", re.IGNORECASE)
_STATION_NAME = re.compile(r"^station_(\d+)$")
_FIELD_KEYS = ("image_date", "confidence", "note")
_EDGE_COLUMNS = ("station_id", "side", "image_date", "wayback_release", "confidence", "note")


def export_sample_kml(gpkg_path: str | Path, kml_path: str | Path,
                      document_name: str = "sample transects",
                      station_name: str = "station_{}") -> Path:
    """Write both sample layers to one KML. Each line is named station_id_side.

    station_name formats the station points' names from the station_id."""
    import geopandas as gpd

    gpkg_path = Path(gpkg_path)
    kml_path = Path(kml_path)
    transects = gpd.read_file(gpkg_path, layer="transects")
    stations = gpd.read_file(gpkg_path, layer="stations")
    if transects.empty or stations.empty:
        raise RuntimeError(f"{gpkg_path} is missing the stations or transects layer.")
    transects = transects.to_crs("EPSG:4326")
    stations = stations.to_crs("EPSG:4326")

    line_marks = []
    for row in transects.sort_values(["station_id", "side"]).itertuples(index=False):
        name = f"{int(row.station_id)}_{row.side}"
        style = "left" if row.side == "left" else "right"
        line_marks.append(_placemark(name, style, _line_body(row.geometry)))
    point_marks = []
    for row in stations.sort_values("station_id").itertuples(index=False):
        name = station_name.format(int(row.station_id))
        point_marks.append(_placemark(name, "station", _point_body(row.geometry)))

    kml = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
  <name>{_xml(document_name)}</name>
  <Style id="left"><LineStyle><color>ffff0000</color><width>3</width></LineStyle><LabelStyle><scale>0.9</scale></LabelStyle></Style>
  <Style id="right"><LineStyle><color>ff0000ff</color><width>3</width></LineStyle><LabelStyle><scale>0.9</scale></LabelStyle></Style>
  <Style id="station"><IconStyle><color>ff00ffff</color><scale>0.7</scale></IconStyle><LabelStyle><scale>0.8</scale></LabelStyle></Style>
  <Folder>
    <name>stations</name>
    {"".join(point_marks)}
  </Folder>
  <Folder>
    <name>transects</name>
    <description>Each line is named station_id_side, for example 13_left. Add a point placemark with that same name where the flood edge crosses the line.</description>
    {"".join(line_marks)}
  </Folder>
</Document>
</kml>
"""
    kml_path.write_text(kml, encoding="utf-8")
    return kml_path


def convert_edge_kml(kml_path: str | Path, transects_gpkg: str | Path, out_path: str | Path,
                     *, crs: str = "EPSG:32645", max_distance_m: float = SNAP_LIMIT_M) -> Path:
    """Snap named placemarks onto their transect lines and write the GeoJSON.

    Raises, and does not write the file, if a point is more than max_distance_m
    from its line or a required field is missing.
    """
    import geopandas as gpd
    from shapely.geometry import Point

    edges = read_edge_placemarks(kml_path)
    lines = gpd.read_file(transects_gpkg, layer="transects")
    if lines.crs is None or lines.crs.to_string() != crs:
        raise RuntimeError(f"Transect CRS is {lines.crs}, expected {crs}.")
    lookup = {(int(r.station_id), r.side): r.geometry for r in lines.itertuples(index=False)}

    points = gpd.GeoDataFrame(
        {"name": [e["name"] for e in edges]},
        geometry=[Point(e["lon"], e["lat"]) for e in edges],
        crs="EPSG:4326",
    ).to_crs(crs)

    rows = []
    problems = []
    for edge, point in zip(edges, points.geometry):
        key = (edge["station_id"], edge["side"])
        line = lookup.get(key)
        if line is None:
            problems.append(f"{edge['name']} is not one of the sample transects.")
            continue
        snapped = line.interpolate(line.project(point))
        gap_m = float(point.distance(snapped))
        if gap_m > max_distance_m:
            problems.append(
                f"{edge['name']} is {gap_m:.1f} m from its transect line "
                f"(limit {max_distance_m:.0f} m).")
            continue
        rows.append({
            "station_id": edge["station_id"],
            "side": edge["side"],
            "image_date": edge["image_date"],
            "confidence": edge["confidence"],
            "note": edge["note"],
            "geometry": snapped,
        })
    if problems:
        raise RuntimeError("No validation edges were written.\n" + "\n".join(problems))

    out = gpd.GeoDataFrame(rows, crs=crs).sort_values(["station_id", "side"]).reset_index(drop=True)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_file(out_path, driver="GeoJSON")
    return out_path


def side_from_channel(offset_m: float, centre_offset_m: float | None,
                      recorded: str | None = None) -> str:
    """Bank of a point from the re-centred channel (14.18, 14.25).

    Offsets are metres from the station, left-positive. A point on the
    left of the channel centre is left-bank. An unresolved channel is
    treated as centred on the station. A point sitting on the centre
    keeps the recorded side when that is left or right.
    """
    import math

    centre = 0.0 if centre_offset_m is None or not math.isfinite(float(centre_offset_m)) \
        else float(centre_offset_m)
    offset = float(offset_m)
    if offset > centre + 1e-9:
        return "left"
    if offset < centre - 1e-9:
        return "right"
    if recorded in ("left", "right"):
        return recorded
    return "left"


def snap_to_station_transects(point, left_line, right_line):
    """Snap to the nearer of the two lines. Offset is metres from the station, left-positive.

    Both lines start at the station. The left line runs outward on the left
    bank looking downstream; the right line runs outward on the right.
    """
    d_left = float(point.distance(left_line))
    d_right = float(point.distance(right_line))
    if d_left <= d_right:
        along = float(left_line.project(point))
        return left_line.interpolate(along), along, d_left, "left"
    along = float(right_line.project(point))
    return right_line.interpolate(along), -along, d_right, "right"


def find_edge_gpkg(folder: str | Path):
    """validation_edges.gpkg, or the QGIS file saved as `.gpkg` with that layer."""
    import pyogrio

    folder = Path(folder)
    named = folder / "validation_edges.gpkg"
    if named.exists():
        return named
    hidden = folder / ".gpkg"
    if hidden.exists():
        layers = [name for name, *_ in pyogrio.list_layers(hidden)]
        if "validation_edges" in layers:
            return hidden
    raise FileNotFoundError(
        f"{named} is missing, and {hidden} has no validation_edges layer.")


EDGE_CSV_NAME = "validation_edges.csv"
EDGE_CSV_COLUMNS = ("station_id", "side", "image_date", "confidence", "note", "edge_lat",
                    "edge_lon")
_FROM_CSV = {"written_from": EDGE_CSV_NAME}


def read_edge_csv(csv_path: str | Path):
    """The user's edge CSV as WGS84 points. Refuses a malformed file; writes nothing.

    Columns: station_id, side (left/right), image_date (YYYY-MM-DD, or month/day/year
    such as 1/2/2026 for 2 January 2026; written out as YYYY-MM-DD), confidence (1-3),
    note (may be empty), edge_lat, edge_lon (decimal degrees, WGS84).
    """
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import Point

    path = Path(csv_path)
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    frame.columns = [c.strip() for c in frame.columns]
    missing = [c for c in EDGE_CSV_COLUMNS if c not in frame.columns]
    extra = [c for c in frame.columns if c not in EDGE_CSV_COLUMNS]
    if missing or extra:
        raise RuntimeError(f"{path.name} must have exactly the columns "
                           f"{', '.join(EDGE_CSV_COLUMNS)}; missing {missing}, extra {extra}.")
    frame = frame[[any(str(v).strip() for v in row) for row in frame.itertuples(index=False)]]
    rows, problems = [], []
    for i, row in enumerate(frame.itertuples(index=False), start=2):   # line 1 is the header
        label = f"line {i}"
        try:
            station_id = int(str(row.station_id).strip())
        except ValueError:
            problems.append(f"{label}: station_id {row.station_id!r} is not a whole number.")
            continue
        label = f"line {i} (station {station_id})"
        side = str(row.side).strip().lower()
        if side not in ("left", "right"):
            problems.append(f"{label}: side {row.side!r} must be left or right.")
            continue
        try:
            lat, lon = float(row.edge_lat), float(row.edge_lon)
        except ValueError:
            problems.append(f"{label}: edge_lat/edge_lon must be decimal degrees.")
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            problems.append(f"{label}: edge_lat {lat} / edge_lon {lon} are not valid degrees.")
            continue
        try:
            image_date = _csv_date(str(row.image_date).strip(), label)
            confidence = _as_confidence(str(row.confidence).strip() or None, label)
        except RuntimeError as exc:
            problems.append(str(exc))
            continue
        rows.append({"station_id": station_id, "side": side, "image_date": image_date,
                     "confidence": confidence, "note": str(row.note).strip(),
                     "geometry": Point(lon, lat)})
    keys = [(r["station_id"], r["side"]) for r in rows]
    dupes = sorted({k for k in keys if keys.count(k) > 1})
    problems += [f"station {s} {side} appears more than once." for s, side in dupes]
    if problems:
        raise RuntimeError(f"{path.name} was refused; nothing was written.\n" + "\n".join(problems))
    if not rows:
        raise RuntimeError(f"{path.name} has no points.")
    return gpd.GeoDataFrame(rows, crs="EPSG:4326")


def csv_to_edge_gpkg(csv_path: str | Path, gpkg_path: str | Path, crs: str) -> Path:
    """Write the edge CSV as validation_edges.gpkg (layer validation_edges, in crs).

    A validation_edges.gpkg that was not written from the CSV (a QGIS layer) is never
    overwritten."""
    import pyogrio

    gpkg_path = Path(gpkg_path)
    if gpkg_path.exists():
        meta = pyogrio.read_info(gpkg_path, layer="validation_edges").get("layer_metadata") or {}
        if meta.get("written_from") != EDGE_CSV_NAME:
            raise RuntimeError(f"Both {EDGE_CSV_NAME} and a {gpkg_path.name} not written from "
                               "it are present; keep only one. Nothing was written.")
    points = read_edge_csv(csv_path).to_crs(crs)
    if gpkg_path.exists():
        gpkg_path.unlink()
    points.to_file(gpkg_path, layer="validation_edges", driver="GPKG",
                   layer_metadata=_FROM_CSV)
    return gpkg_path


def export_event_sample_kml(event_key: str, config_path: str | Path | None = None) -> Path:
    """CLI entry: sample_transects.kml beside the drawn sample. Never overwrites a KML.

    Station points are named by station_id, lines station_id_side. Only names and
    geometry are written (no forest or outcome values)."""
    from .config import PROJECT_ROOT, load_config

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    folder = PROJECT_ROOT / cfg.paths.data / event_key / "validation"
    gpkg = folder / "sample_transects.gpkg"
    if not gpkg.exists():
        raise FileNotFoundError(f"{gpkg} is missing. Draw the validation sample first.")
    kml = folder / "sample_transects.kml"
    if kml.exists():
        raise RuntimeError(f"{kml} already exists; it is not overwritten.")
    export_sample_kml(gpkg, kml, document_name=f"{event_key} validation sample",
                      station_name="{}")
    print(kml.relative_to(PROJECT_ROOT))
    return kml


def convert_edge_gpkg(gpkg_path: str | Path, transects_gpkg: str | Path, out_path: str | Path,
                      centres: dict, *, crs: str = "EPSG:32645",
                      max_distance_m: float = SNAP_LIMIT_M,
                      date_window: tuple[str, str] = (VALIDATION_DATE_START,
                                                      VALIDATION_DATE_END)) -> tuple[Path, dict]:
    """Snap digitised points, assign side from the re-centred channel, write GeoJSON.

    A point more than max_distance_m from both of its station's transect
    lines refuses the whole conversion: nothing is written. Points whose
    image_date is outside date_window (both ends included) are left out and counted.
    """
    date_start, date_end = date_window
    import geopandas as gpd
    import pandas as pd

    gpkg_path = Path(gpkg_path)
    raw = gpd.read_file(gpkg_path, layer="validation_edges")
    if raw.crs is None:
        raise RuntimeError(f"{gpkg_path} has no CRS.")
    if raw.crs.to_string() != crs:
        raw = raw.to_crs(crs)
    lines = gpd.read_file(transects_gpkg, layer="transects")
    if lines.crs is None or lines.crs.to_string() != crs:
        raise RuntimeError(f"Transect CRS is {lines.crs}, expected {crs}.")
    lookup = {}
    for row in lines.itertuples(index=False):
        lookup.setdefault(int(row.station_id), {})[row.side] = row.geometry

    rows = []
    refused = []
    mismatches = []
    snaps = []
    outside_window = 0
    for i, row in enumerate(raw.itertuples(index=False)):
        geom = getattr(row, "geometry", None)
        if geom is None or geom.is_empty:
            continue
        image_date = _as_iso_date(getattr(row, "image_date", None),
                                  f"row {i + 1}", required=False)
        if image_date is None or not date_start <= image_date <= date_end:
            outside_window += 1
            continue
        try:
            station_id = int(row.station_id)
        except (TypeError, ValueError):
            refused.append(f"row {i + 1} has no station_id.")
            continue
        recorded = str(getattr(row, "side", "") or "").strip().lower() or None
        pair = lookup.get(station_id)
        if pair is None or set(pair) != {"left", "right"}:
            refused.append(f"station {station_id} is not one of the sample transects.")
            continue
        snapped, offset, gap_m, line_side = snap_to_station_transects(
            geom, pair["left"], pair["right"])
        if gap_m > max_distance_m + 1e-9:
            refused.append(
                f"station {station_id} is {gap_m:.1f} m from its transect "
                f"(limit {max_distance_m:.0f} m).")
            continue
        derived = side_from_channel(offset, centres.get(station_id), recorded)
        if recorded in ("left", "right") and recorded != derived:
            mismatches.append({
                "station_id": station_id,
                "recorded": recorded,
                "derived": derived,
                "line": line_side,
            })
        wayback = _as_iso_date(getattr(row, "wayback_release", None), f"station {station_id}",
                               required=False)
        confidence = _as_confidence(getattr(row, "confidence", None), f"station {station_id}")
        note = getattr(row, "note", None)
        note = "" if note is None or (isinstance(note, float) and pd.isna(note)) else str(note)
        rows.append({
            "station_id": station_id,
            "side": derived,
            "image_date": image_date,
            "wayback_release": wayback,
            "confidence": confidence,
            "note": note,
            "geometry": snapped,
        })
        snaps.append({"station_id": station_id, "side": derived, "snap_m": float(gap_m)})
    if refused:
        raise RuntimeError("Conversion refused; no GeoJSON was written.\n" + "\n".join(refused))
    if not rows:
        raise RuntimeError(
            "No validation edges were written.\n" + "\n".join(refused or ["No point features."]))
    out = gpd.GeoDataFrame(rows, crs=crs)
    out = out[list(_EDGE_COLUMNS) + ["geometry"]]
    out = out.sort_values(["station_id", "side"]).reset_index(drop=True)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_file(out_path, driver="GeoJSON")
    report = {
        "n_read": int(raw.geometry.notna().sum()) if "geometry" in raw else len(raw),
        "n_written": len(out),
        "refused": refused,
        "outside_date_window": outside_window,
        "mismatches": mismatches,
        "snap_distances": sorted(snaps, key=lambda s: (s["station_id"], s["side"])),
        "points_per_image_date": out["image_date"].value_counts().sort_index().to_dict(),
        "stations_written": sorted(int(s) for s in out["station_id"].unique()),
    }
    return out_path, report


def _as_iso_date(value, label: str, required: bool = True) -> str | None:
    import pandas as pd

    if value is None or (isinstance(value, float) and pd.isna(value)) or value == "":
        if required:
            raise RuntimeError(f"{label} is missing image_date.")
        return None
    if hasattr(value, "strftime"):
        return value.strftime("%Y-%m-%d")
    text = str(value).strip()
    try:
        return datetime.strptime(text[:10], "%Y-%m-%d").date().isoformat()
    except ValueError:
        raise RuntimeError(f"{label} has date {value!r}. Use a date like 2021-10-15.") from None


def _csv_date(text: str, label: str) -> str:
    """A CSV image_date as YYYY-MM-DD. Accepts YYYY-MM-DD or month/day/year (M/D/YYYY)."""
    if not text:
        raise RuntimeError(f"{label} is missing image_date.")
    for fmt in ("%Y-%m-%d", "%m/%d/%Y"):
        try:
            return datetime.strptime(text, fmt).date().isoformat()
        except ValueError:
            pass
    raise RuntimeError(f"{label} has date {text!r}. Use a date like 2021-10-15 or "
                       "month/day/year like 10/15/2021.")


def _as_confidence(value, label: str) -> int:
    import pandas as pd

    if value is None or (isinstance(value, float) and pd.isna(value)):
        raise RuntimeError(
            f"{label} is missing confidence. Use 1 (sure), 2 (fairly sure) or 3 (unsure).")
    try:
        number = int(value)
    except (TypeError, ValueError):
        raise RuntimeError(
            f"{label} has confidence {value!r}. Use 1 (sure), 2 (fairly sure) or 3 (unsure)."
        ) from None
    if number not in (1, 2, 3):
        raise RuntimeError(
            f"{label} has confidence {value!r}. Use 1 (sure), 2 (fairly sure) or 3 (unsure).")
    return number


def import_edge_kml(event_key: str, kml_path: str | Path,
                    config_path: str | Path | None = None) -> Path:
    """CLI entry: read the placemark KML and write validation_edges.geojson."""
    from .config import PROJECT_ROOT, load_config

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    if cfg.crs != "EPSG:32645":
        raise RuntimeError(f"Validation edges are written in EPSG:32645, and config crs is {cfg.crs}.")
    folder = PROJECT_ROOT / cfg.paths.data / event_key / "validation"
    gpkg = folder / "sample_transects.gpkg"
    if not gpkg.exists():
        raise FileNotFoundError(f"{gpkg} is missing. Draw the validation sample first.")
    out = folder / "validation_edges.geojson"
    written = convert_edge_kml(kml_path, gpkg, out, crs=cfg.crs, max_distance_m=SNAP_LIMIT_M)
    print(f"Wrote {written.relative_to(PROJECT_ROOT)}")
    return written


def read_edge_placemarks(kml_path: str | Path) -> list[dict]:
    """Point placemarks named station_id_side, with image_date, confidence and note."""
    path = Path(kml_path)
    if not path.exists():
        raise FileNotFoundError(f"KML file not found: {path}")
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise RuntimeError(f"{path} is not valid KML: {exc}") from exc

    edges = []
    seen = set()
    unexpected = []
    for placemark in _named(root, "Placemark"):
        lonlat = _point_lonlat(placemark)
        if lonlat is None:
            continue
        name = _direct_text(placemark, "name").strip()
        if _STATION_NAME.match(name):
            continue
        match = _EDGE_NAME.match(name)
        if match is None:
            unexpected.append(name or "(no name)")
            continue
        station_id = int(match.group(1))
        side = match.group(2).lower()
        label = f"{station_id}_{side}"
        if label in seen:
            raise RuntimeError(f"More than one point placemark is named {label}.")
        seen.add(label)
        fields = _description_fields(_direct_text(placemark, "description"))
        fields.update(_extended_fields(placemark))
        edges.append({
            "name": label,
            "station_id": station_id,
            "side": side,
            "lon": lonlat[0],
            "lat": lonlat[1],
            **_checked_fields(label, fields),
        })
    if unexpected:
        raise RuntimeError(
            "Point placemarks must be named like 13_left. "
            f"These are not: {', '.join(unexpected)}.")
    if not edges:
        raise RuntimeError(
            f"{path} has no point placemarks named station_id_side. "
            "In the description, one field per line: "
            "image_date: 2021-10-15, confidence: 1, note: text. "
            "confidence is 1 (sure), 2 (fairly sure) or 3 (unsure).")
    return edges


def _checked_fields(label: str, fields: dict) -> dict:
    missing = [key for key in ("image_date", "confidence") if key not in fields or fields[key] == ""]
    if missing:
        raise RuntimeError(
            f"{label} is missing {', '.join(missing)}. In the placemark description, "
            "one field per line: image_date: 2021-10-15, confidence: 1, note: text. "
            "confidence is 1 (sure), 2 (fairly sure) or 3 (unsure).")
    try:
        image_date = datetime.strptime(fields["image_date"], "%Y-%m-%d").date().isoformat()
    except ValueError:
        raise RuntimeError(
            f"{label} has image_date {fields['image_date']!r}. Use a date like 2021-10-15."
        ) from None
    if fields["confidence"] not in {"1", "2", "3"}:
        raise RuntimeError(
            f"{label} has confidence {fields['confidence']!r}. "
            "Use 1 (sure), 2 (fairly sure) or 3 (unsure).")
    return {
        "image_date": image_date,
        "confidence": int(fields["confidence"]),
        "note": fields.get("note", ""),
    }


def _description_fields(text: str) -> dict:
    plain = unescape(text or "")
    plain = re.sub(r"(?i)<br\s*/?>", "\n", plain)
    plain = re.sub(r"(?i)</div>|</p>|</tr>", "\n", plain)
    plain = re.sub(r"<[^>]+>", "", plain)
    fields = {}
    for line in plain.splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip().lower().replace(" ", "_")
        if key in _FIELD_KEYS:
            fields[key] = value.strip()
    return fields


def _extended_fields(placemark: ET.Element) -> dict:
    fields = {}
    for el in placemark.iter():
        if _local(el.tag) not in ("Data", "SimpleData"):
            continue
        key = el.attrib.get("name", "").strip().lower().replace(" ", "_")
        if key not in _FIELD_KEYS:
            continue
        if _local(el.tag) == "SimpleData":
            fields[key] = (el.text or "").strip()
            continue
        for child in el:
            if _local(child.tag) == "value":
                fields[key] = (child.text or "").strip()
    return fields


def _point_lonlat(placemark: ET.Element) -> tuple[float, float] | None:
    for el in placemark.iter():
        if _local(el.tag) != "Point":
            continue
        for child in el.iter():
            if _local(child.tag) == "coordinates" and (child.text or "").strip():
                token = child.text.strip().split()[0]
                lon, lat, *_ = token.split(",")
                return float(lon), float(lat)
        name = _direct_text(placemark, "name") or "(no name)"
        raise RuntimeError(f"Point placemark {name} has no coordinates.")
    return None


def _direct_text(el: ET.Element, name: str) -> str:
    for child in el:
        if _local(child.tag) == name:
            return child.text or ""
    return ""


def _named(el: ET.Element, name: str):
    return [child for child in el.iter() if _local(child.tag) == name]


def _local(tag: str) -> str:
    return tag.split("}", 1)[-1]


def _xml(text: str) -> str:
    return (text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _placemark(name: str, style: str, body: str) -> str:
    return (
        "<Placemark>"
        f"<name>{_xml(name)}</name>"
        f"<styleUrl>#{style}</styleUrl>"
        f"{body}"
        "</Placemark>\n"
    )


def _coords(geom) -> list[str]:
    if geom.geom_type == "LineString":
        sequences = [geom.coords]
    elif geom.geom_type == "MultiLineString":
        sequences = [part.coords for part in geom.geoms]
    elif geom.geom_type == "Point":
        sequences = [geom.coords]
    else:
        raise RuntimeError(f"Cannot write a {geom.geom_type} to the validation KML.")
    tuples = []
    for sequence in sequences:
        tuples.append(" ".join(f"{x:.7f},{y:.7f},0" for x, y, *rest in sequence))
    return tuples


def _line_body(geom) -> str:
    parts = []
    for coordinates in _coords(geom):
        parts.append(
            "<LineString><tessellate>1</tessellate>"
            "<altitudeMode>clampToGround</altitudeMode>"
            f"<coordinates>{coordinates}</coordinates></LineString>"
        )
    body = "".join(parts)
    if len(parts) > 1:
        body = f"<MultiGeometry>{body}</MultiGeometry>"
    return body


def _point_body(geom) -> str:
    coordinates = _coords(geom)[0]
    return (
        "<Point><altitudeMode>clampToGround</altitudeMode>"
        f"<coordinates>{coordinates}</coordinates></Point>"
    )
