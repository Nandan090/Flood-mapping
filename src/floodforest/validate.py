"""Validation sample for digitising (spec section 9, clarification 14.13).

Draws the station sample and writes it without forest or outcome fields.
The comparison with digitised edges is not built yet: with no edges file,
the step reports that the event is not validated.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

_MONTHS = ("January", "February", "March", "April", "May", "June", "July",
           "August", "September", "October", "November", "December")


def eligible_station_forest(df: pd.DataFrame, forest_column: str) -> pd.Series:
    """Mean forest fraction per station, using only analysis-set banks.

    A station with neither bank in the analysis set is left out. Where one
    bank is in the set, that bank's value is the station's value.
    """
    missing = [c for c in ("analysis_set", forest_column) if c not in df.columns]
    if missing:
        raise RuntimeError(
            f"The transect table is missing {missing}. Run the table step first.")
    use = df[df["analysis_set"] == 1]
    if use.empty:
        raise RuntimeError("The analysis set is empty, so there is nothing to sample.")
    gap = use[forest_column].isna()
    if gap.any():
        row = use.loc[gap].iloc[0]
        raise RuntimeError(
            f"{forest_column} is missing for {int(gap.sum())} analysis-set rows, "
            f"starting with station {int(row.station_id)} {row.side}.")
    return use.groupby("station_id")[forest_column].mean()


def tercile_groups(values: pd.Series, n_bins: int) -> list[np.ndarray]:
    """Equal-count groups, lowest forest fraction first.

    Stations with the same fraction are ordered by station id, so a station
    on a boundary always falls the same way.
    """
    if n_bins < 1:
        raise ValueError(f"validation_terciles must be at least 1, got {n_bins}.")
    frame = values.rename("value").rename_axis("station_id").reset_index()
    frame = frame.sort_values(["value", "station_id"], kind="mergesort")
    if len(frame) < n_bins:
        raise RuntimeError(
            f"Only {len(frame)} stations are eligible, which is fewer than "
            f"{n_bins} forest-fraction groups.")
    ids = frame["station_id"].to_numpy()
    return [np.asarray(part) for part in np.array_split(ids, n_bins)]


def group_sizes(n: int, n_bins: int) -> list[int]:
    """How many stations to draw from each group, lowest group first.

    The counts differ by at most one. Any extras go to the lower groups.
    """
    if n < n_bins:
        raise ValueError(
            f"validation_sample_size ({n}) is smaller than the number of groups ({n_bins}).")
    base, remainder = divmod(n, n_bins)
    return [base + (1 if i < remainder else 0) for i in range(n_bins)]


def draw_sample(groups: list[np.ndarray], n: int, seed: int) -> list[int]:
    """Station ids, stratified across groups. The result is sorted by id."""
    sizes = group_sizes(n, len(groups))
    rng = np.random.default_rng(seed)
    picked: list[int] = []
    for index, (group, size) in enumerate(zip(groups, sizes), start=1):
        if size > len(group):
            raise RuntimeError(
                f"Forest-fraction group {index} has {len(group)} stations, "
                f"and the sample asks for {size}.")
        if size:
            picked.extend(int(v) for v in rng.choice(group, size=size, replace=False))
    return sorted(picked)


def sample_layers(lines, station_ids: list[int]):
    """Point per station and both of its transect lines, with place fields only."""
    import geopandas as gpd
    from shapely.geometry import Point

    needed = {"station_id", "side", "river_km", "station_x", "station_y"}
    missing = sorted(needed - set(lines.columns))
    if missing:
        raise RuntimeError(f"The transect table is missing {missing}.")
    chosen = lines[lines["station_id"].isin(set(station_ids))].copy()
    found = set(chosen["station_id"].unique())
    absent = [int(s) for s in station_ids if s not in found]
    if absent:
        raise RuntimeError(
            f"Sampled stations are missing from the transect table: {absent}.")
    incomplete = []
    for station_id, part in chosen.groupby("station_id"):
        if set(part["side"]) != {"left", "right"}:
            incomplete.append(int(station_id))
    if incomplete:
        raise RuntimeError(
            "These sampled stations do not have both a left and a right "
            f"transect: {incomplete}.")

    out_lines = gpd.GeoDataFrame(
        {
            "station_id": chosen["station_id"].astype(int),
            "side": chosen["side"].astype(str),
            "river_km": chosen["river_km"].astype(float),
        },
        geometry=list(chosen.geometry),
        crs=lines.crs,
    ).sort_values(["station_id", "side"]).reset_index(drop=True)

    plain = chosen[["station_id", "side", "river_km", "station_x", "station_y"]]
    first = plain.sort_values(["station_id", "side"]).groupby("station_id", as_index=False).first()
    points = gpd.GeoDataFrame(
        {
            "station_id": first["station_id"].astype(int),
            "river_km": first["river_km"].astype(float),
        },
        geometry=[Point(float(x), float(y)) for x, y in zip(first.station_x, first.station_y)],
        crs=lines.crs,
    ).sort_values("station_id").reset_index(drop=True)
    _assert_place_only(points, ("station_id", "river_km"))
    _assert_place_only(out_lines, ("station_id", "side", "river_km"))
    return points, out_lines


def _assert_place_only(frame, allowed: tuple[str, ...]) -> None:
    extra = [c for c in frame.columns if c != frame.geometry.name and c not in allowed]
    if extra:
        raise RuntimeError(f"The validation sample would include {extra}.")


def digitising_instructions(event_name: str, when: str, crs: str) -> str:
    """Plain-language steps for marking the flood edge. No measured results."""
    return f"""# Marking the flood edge

This is for the {event_name} ({when}).

Mark where the debris flood stopped, on each bank, at the 40 places in this folder. Mark only what you can see in imagery taken after the flood. Use the earliest clear view you can find: imagery from years later can show repairs and newer deposits.

Do not open the results table, the cross-section figures, or any file outside this folder.

## Files

- `sample_stations.geojson` — one point on the river at each place.
- `sample_transects.gpkg` — the same points (layer `stations`) and the two lines at each place (layer `transects`).

The layers use {crs}. Left and right are as you look downstream. Each line starts at the station and runs out across one bank. The line's `side` field says which bank.

## In QGIS

1. Add both layers from `sample_transects.gpkg`, and a satellite image from after the flood.
2. Turn on snapping to the transect lines, so a point you add sits on the line.
3. Create a new point layer in the same CRS, with two fields: `station_id` (a whole number) and `side` (text). Save it in this folder as `validation_edges.geojson`.

## At each place

1. Read the station's `station_id`.
2. On the left line and on the right line, add one point at the outer limit of the flood: the furthest point from the river, along that line, where debris, erosion, or damaged ground from this flood is still visible.
3. Snap the point to the line. Set `station_id` from the station, and set `side` to `left` or `right`.
4. If the flood did not leave the river on that bank, put the point on the line at the edge of the channel.
5. If cloud, shadow, or trees hide the edge, skip that bank. Leave it out. Do not guess.

A place can have two points, one, or none. Save the layer when you finish.

Any published flood-extent map goes in this same folder, beside `validation_edges.geojson`.
"""


def _long_date(iso: str) -> str:
    year, month, day = (int(part) for part in str(iso).split("-"))
    return f"{day} {_MONTHS[month - 1]} {year}"


def _when(event) -> str:
    dates = [_long_date(d) for d in event.event_dates]
    if len(dates) == 1:
        return dates[0]
    return f"{dates[0]} to {dates[-1]}"


def run_validation_sample(event_key: str, config_path: str | Path | None = None) -> list[Path]:
    """Write the blind station sample, the QGIS layers, and the digitising notes."""
    import geopandas as gpd

    from .config import PROJECT_ROOT, load_config

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    gpkg_in = PROJECT_ROOT / cfg.paths.outputs / event_key / "transects.gpkg"
    if not gpkg_in.exists():
        raise FileNotFoundError(f"{gpkg_in} is missing. Run the table step first.")
    lines = gpd.read_file(gpkg_in, layer="transects")
    if lines.crs is None or lines.crs.to_string() != cfg.crs:
        raise RuntimeError(f"Transects CRS is {lines.crs}, expected {cfg.crs}.")

    forest_column = f"forest_frac_{int(cfg.belt_width_m)}"
    values = eligible_station_forest(lines, forest_column)
    n_bins = int(cfg.validation_terciles)
    groups = tercile_groups(values, n_bins)
    chosen = draw_sample(groups, int(cfg.validation_sample_size), int(cfg.seed))
    points, sample_lines = sample_layers(lines, chosen)
    counts = [int(np.isin(group, chosen).sum()) for group in groups]

    folder = PROJECT_ROOT / cfg.paths.data / event_key / "validation"
    folder.mkdir(parents=True, exist_ok=True)
    stations_path = folder / "sample_stations.geojson"
    transects_path = folder / "sample_transects.gpkg"
    notes_path = folder / "how_to_digitise.md"
    if transects_path.exists():
        transects_path.unlink()
    points.to_file(stations_path, driver="GeoJSON")
    sample_lines.to_file(transects_path, layer="transects", driver="GPKG")
    points.to_file(transects_path, layer="stations", driver="GPKG", mode="a")
    notes_path.write_text(digitising_instructions(event.name, _when(event), cfg.crs), encoding="utf-8")

    prov_path = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    if not prov_path.exists():
        raise FileNotFoundError(f"{prov_path} is missing. Run the earlier steps first.")
    provenance = json.loads(prov_path.read_text(encoding="utf-8"))
    written = [stations_path, transects_path, notes_path]
    provenance["validation_sample"] = {
        "n": len(chosen),
        "seed": int(cfg.seed),
        "eligible_stations": int(values.size),
        "per_group_low_to_high": counts,
        "station_ids": chosen,
        "rule": (
            f"Stations with at least one analysis-set side. Equal-count groups of "
            f"the mean {forest_column} on those sides, lowest first. Extra stations "
            f"come from the lower groups. Files omit forest and outcome fields (14.13)."
        ),
        "files": [str(p.relative_to(PROJECT_ROOT)).replace("\\", "/") for p in written],
    }
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")

    print(f"Validation sample for {event_key}: {len(chosen)} stations drawn from "
          f"{int(values.size)} eligible stations, {counts} from the lowest "
          f"forest-fraction group to the highest.")
    print(f"{'station_id':>10}  {'river_km':>8}")
    for row in points.itertuples(index=False):
        print(f"{int(row.station_id):10d}  {row.river_km:8.4f}")
    for path in written:
        print(path.relative_to(PROJECT_ROOT))
    edges = folder / "validation_edges.geojson"
    if edges.exists():
        print(f"{edges.name} is present. The comparison with it is not built yet.")
    else:
        print("not validated")
    return written


GROUP_NAMES = ("low", "mid", "high")


def station_group_labels(values: pd.Series, chosen: list[int], n_bins: int) -> pd.Series:
    """Forest-fraction group of each sampled station, lowest group first."""
    groups = tercile_groups(values, n_bins)
    names = GROUP_NAMES[:n_bins] if n_bins <= len(GROUP_NAMES) \
        else tuple(f"group_{i + 1}" for i in range(n_bins))
    label = {}
    for name, group in zip(names, groups):
        for station_id in group:
            if int(station_id) in set(chosen):
                label[int(station_id)] = name
    missing = [s for s in chosen if s not in label]
    if missing:
        raise RuntimeError(f"Sampled stations have no forest-fraction group: {missing}.")
    return pd.Series(label, name="forest_group")


def coverage_rows(stations, checked: set[int], groups: pd.Series) -> pd.DataFrame:
    """One row per sampled station: checked or skipped. No forest or reach values."""
    rows = []
    for row in stations.sort_values("station_id").itertuples(index=False):
        station_id = int(row.station_id)
        rows.append({
            "station_id": station_id,
            "river_km": float(row.river_km),
            "status": "checked" if station_id in checked else "skipped",
            "forest_group": groups.loc[station_id],
        })
    return pd.DataFrame(rows)


def far_from_channel(table: pd.DataFrame, station_ids: list[int], limit_m: float) -> pd.DataFrame:
    """Sampled stations whose re-centred channel is more than limit_m from the station."""
    use = table[table["station_id"].isin(station_ids)]
    first = use.groupby("station_id", as_index=False).first()
    offset = pd.to_numeric(first["centre_offset_m"], errors="coerce")
    far = first[offset.abs() > float(limit_m)].copy()
    return far[["station_id", "river_km"]].sort_values("station_id").reset_index(drop=True)


def import_validation_edges(event_key: str, config_path: str | Path | None = None) -> list[Path]:
    """Convert the digitised GPKG, write GeoJSON, and print coverage. No comparison."""
    import geopandas as gpd

    from .config import PROJECT_ROOT, load_config
    from .validation_kml import convert_edge_gpkg, find_edge_gpkg

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    folder = PROJECT_ROOT / cfg.paths.data / event_key / "validation"
    sample_path = folder / "sample_stations.geojson"
    lines_path = folder / "sample_transects.gpkg"
    table_path = PROJECT_ROOT / cfg.paths.outputs / event_key / "transects.csv"
    for path in (sample_path, lines_path, table_path):
        if not path.exists():
            raise FileNotFoundError(f"{path} is missing. Run the table and sample steps first.")
    source = find_edge_gpkg(folder)
    stations = gpd.read_file(sample_path)
    chosen = [int(s) for s in stations["station_id"].tolist()]
    table = pd.read_csv(table_path, usecols=[
        "station_id", "side", "river_km", "analysis_set",
        f"forest_frac_{int(cfg.belt_width_m)}", "centre_offset_m",
    ])
    centres = (table.groupby("station_id")["centre_offset_m"].first()
               .apply(lambda v: float(v) if pd.notna(v) else float("nan"))
               .to_dict())
    out_path = folder / "validation_edges.geojson"
    written, report = convert_edge_gpkg(
        source, lines_path, out_path, centres,
        crs=cfg.crs, max_distance_m=float(cfg.validation_snap_max_m))

    forest_column = f"forest_frac_{int(cfg.belt_width_m)}"
    values = eligible_station_forest(table, forest_column)
    groups = station_group_labels(values, chosen, int(cfg.validation_terciles))
    coverage = coverage_rows(stations, set(report["stations_written"]), groups)
    coverage_path = folder / "sample_coverage.csv"
    coverage[["station_id", "river_km", "status"]].to_csv(coverage_path, index=False)
    far = far_from_channel(table, chosen, float(cfg.validation_channel_offset_warn_m))

    counts = (coverage.groupby(["forest_group", "status"]).size()
              .unstack(fill_value=0).reindex(columns=["checked", "skipped"], fill_value=0))
    names = GROUP_NAMES[:int(cfg.validation_terciles)]
    counts = counts.reindex(names, fill_value=0)
    counts["total"] = counts["checked"] + counts["skipped"]

    print(f"Read {source.name}, layer validation_edges.")
    print(f"Points written: {report['n_written']}")
    print("Points per image_date:")
    for date, n in report["points_per_image_date"].items():
        print(f"  {date}  {n}")
    if report["refused"]:
        print("Refused points:")
        for item in report["refused"]:
            print(f"  {item}")
    else:
        print("Refused points: 0")
    print("Side mismatches (recorded vs re-centred channel): "
          f"{len(report['mismatches'])}")
    if report["mismatches"]:
        print(f"{'station_id':>10}  {'recorded':<8}  {'derived':<8}  snapped_line")
        for item in report["mismatches"]:
            print(f"{item['station_id']:10d}  {item['recorded']:<8}  {item['derived']:<8}  "
                  f"{item['line']}")
    print(f"Sampled stations more than {cfg.validation_channel_offset_warn_m:g} m "
          f"from the re-centred channel: {len(far)}")
    if len(far):
        print(f"{'station_id':>10}  {'river_km':>8}")
        for row in far.itertuples(index=False):
            print(f"{int(row.station_id):10d}  {row.river_km:8.4f}")
    print("Coverage of the 40 sampled stations:")
    print(f"{'station_id':>10}  {'river_km':>8}  status")
    for row in coverage.itertuples(index=False):
        print(f"{int(row.station_id):10d}  {row.river_km:8.4f}  {row.status}")
    print(f"checked  {int((coverage.status == 'checked').sum())}")
    print(f"skipped  {int((coverage.status == 'skipped').sum())}")
    print("Forest-fraction group (counts only; no forest or reach values):")
    print(f"{'group':<6}  {'checked':>8}  {'skipped':>8}  {'total':>8}")
    for name in names:
        row = counts.loc[name]
        print(f"{name:<6}  {int(row.checked):8d}  {int(row.skipped):8d}  {int(row.total):8d}")
    print("Comparison with digitised reach is not run.")
    print(written.relative_to(PROJECT_ROOT))
    print(coverage_path.relative_to(PROJECT_ROOT))

    prov_path = PROJECT_ROOT / cfg.paths.outputs / event_key / "provenance.json"
    provenance = json.loads(prov_path.read_text(encoding="utf-8"))
    provenance["validation_edges"] = {
        "source": str(source.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        "n_written": report["n_written"],
        "points_per_image_date": {str(k): int(v) for k, v in report["points_per_image_date"].items()},
        "n_refused": len(report["refused"]),
        "n_side_mismatches": len(report["mismatches"]),
        "side_mismatches": report["mismatches"],
        "n_far_from_channel": int(len(far)),
        "far_from_channel_station_ids": [int(s) for s in far["station_id"]],
        "checked": int((coverage.status == "checked").sum()),
        "skipped": int((coverage.status == "skipped").sum()),
        "group_counts": {name: {"checked": int(counts.loc[name, "checked"]),
                                "skipped": int(counts.loc[name, "skipped"])} for name in names},
        "rule": (
            f"Snap to the nearer sample transect; refuse if further than "
            f"{cfg.validation_snap_max_m:g} m. Side is left if the snapped point "
            "lies to the left of the re-centred channel centre (14.25). "
            "Coverage lists checked/skipped only; forest-fraction groups as counts."
        ),
        "files": [
            str(written.relative_to(PROJECT_ROOT)).replace("\\", "/"),
            str(coverage_path.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        ],
    }
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    return [written, coverage_path]
