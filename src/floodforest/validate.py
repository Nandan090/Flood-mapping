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
    coverage[["station_id", "status", "forest_group"]].to_csv(coverage_path, index=False)
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
        print("station_id")
        for row in far.itertuples(index=False):
            print(int(row.station_id))
    print("Coverage of the 40 sampled stations:")
    print(f"{'station_id':>10}  {'status':<8}  forest_group")
    for row in coverage.itertuples(index=False):
        print(f"{int(row.station_id):10d}  {row.status:<8}  {row.forest_group}")
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


def standardized_rank_slope(x, y) -> float:
    """OLS slope of average-ranked, z-scored y on x; equals Spearman's rho."""
    from scipy.stats import rankdata

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 2:
        raise ValueError("Rank correlation needs two equally sized one-dimensional arrays.")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Rank correlation inputs must be finite.")
    x_rank = rankdata(x, method="average")
    y_rank = rankdata(y, method="average")
    x_sd = float(x_rank.std(ddof=0))
    y_sd = float(y_rank.std(ddof=0))
    if x_sd == 0 or y_sd == 0:
        raise ValueError("Rank correlation is undefined for a constant input.")
    x_z = (x_rank - x_rank.mean()) / x_sd
    y_z = (y_rank - y_rank.mean()) / y_sd
    return float(np.mean(x_z * y_z))


def prepare_validation_comparison(edges, lines, table: pd.DataFrame,
                                  coverage: pd.DataFrame,
                                  forest_column: str = "forest_frac_100",
                                  valid_fraction_min: float = 0.7,
                                  snap_max_m: float = 20.0) -> pd.DataFrame:
    """Join digitised edges to table outcomes and derive bank-top distances."""
    edge_required = {"station_id", "side", "image_date", "confidence", "note", "geometry"}
    line_required = {"station_id", "side", "geometry"}
    table_required = {
        "station_id", "side", "block_id", "station_x", "station_y",
        "bank_top_from_station_m", "reach_distance_m_fab", "valid_fraction_fab",
        "near_bank_nodata_fab", "censored_fab", "analysis_set", "exclude_reasons",
        "detached_reach_m_fab", "z0_fab", forest_column,
    }
    coverage_required = {"station_id", "forest_group"}
    for label, frame, required in (
        ("validation edges", edges, edge_required),
        ("validation transects", lines, line_required),
        ("transect table", table, table_required),
        ("sample coverage", coverage, coverage_required),
    ):
        missing = sorted(required - set(frame.columns))
        if missing:
            raise RuntimeError(f"{label} is missing required columns: {missing}.")
    if edges.crs is None or edges.crs.to_epsg() != 32645:
        raise RuntimeError(f"Validation edges CRS is {edges.crs}, expected EPSG:32645.")
    if lines.crs is None or lines.crs.to_epsg() != 32645:
        raise RuntimeError(f"Validation transect CRS is {lines.crs}, expected EPSG:32645.")
    if table.duplicated(["station_id", "side"]).any():
        raise RuntimeError("Transect table has duplicate station_id/side rows.")
    if coverage["station_id"].duplicated().any():
        raise RuntimeError("Sample coverage has duplicate station_id rows.")

    fields = [
        "station_id", "side", "block_id", "station_x", "station_y",
        "bank_top_from_station_m", "reach_distance_m_fab", "valid_fraction_fab",
        "near_bank_nodata_fab", "censored_fab", "analysis_set", "exclude_reasons",
        "detached_reach_m_fab", "z0_fab", forest_column,
    ]
    joined =edges.merge(table[fields], on=["station_id", "side"], how="left",
                         validate="one_to_one")
    joined = joined.merge(coverage[["station_id", "forest_group"]], on="station_id",
                          how="left", validate="many_to_one")
    missing_rows = joined[["block_id", "bank_top_from_station_m", forest_column,
                           "forest_group"]].isna().any(axis=1)
    if missing_rows.any():
        bad = joined.loc[missing_rows, ["station_id", "side"]].to_dict("records")
        raise RuntimeError(f"Validation points lack table or forest coverage joins: {bad}.")
    groups = joined["forest_group"].astype(str).str.strip().str.lower()
    group_names = {"low": "Low", "mid": "Mid", "high": "High"}
    if not groups.isin(group_names).all():
        raise RuntimeError(f"Unexpected forest-fraction groups: {sorted(groups.unique())}.")
    joined["forest_group"] = groups.map(group_names)
    if joined[forest_column].isna().any():
        raise RuntimeError("Forest fraction is missing for one or more validation points.")
    if joined[["image_date", "confidence"]].isna().any(axis=None):
        raise RuntimeError("Validation image_date or confidence is missing.")

    line_lookup = {}
    for row in lines.itertuples(index=False):
        key = (int(row.station_id), str(row.side))
        if key in line_lookup:
            raise RuntimeError(f"Validation transects contain duplicate key {key}.")
        line_lookup[key] = row.geometry

    validated = []
    for row in joined.itertuples(index=False):
        key = (int(row.station_id), str(row.side))
        line = line_lookup.get(key)
        if line is None or line.is_empty:
            raise RuntimeError(f"No usable transect geometry for station/side {key}.")
        point = row.geometry
        if point is None or point.is_empty or point.geom_type != "Point":
            raise RuntimeError(f"Invalid validation point for station/side {key}.")
        gap_m = float(line.distance(point))
        if gap_m > snap_max_m + 1e-9:
            raise RuntimeError(
                f"Validation point {key} is {gap_m:.2f} m from its transect; "
                f"the existing {snap_max_m:.0f} m snap rule is not met.")
        station = (float(row.station_x), float(row.station_y))
        start = line.interpolate(0.0)
        if float(start.distance(point.__class__(*station))) > 0.01:
            raise RuntimeError(f"Transect {key} does not start at its station point.")
        along_m = float(line.project(point))
        validated.append(along_m - float(row.bank_top_from_station_m))

    joined["validated_reach_m"] = validated
    joined["computed_reach_m"] = pd.to_numeric(joined["reach_distance_m_fab"],
                                                errors="coerce")
    joined["diff_m"] = joined["computed_reach_m"] - joined["validated_reach_m"]
    joined["detached_reach_m"] = pd.to_numeric(joined["detached_reach_m_fab"], errors="coerce")
    joined["validated_minus_detached_m"] = joined["validated_reach_m"] - joined["detached_reach_m"]
    # Same sign as diff_m (computed minus validated), for the detached sensitivity table.
    joined["diff_detached_m"] = joined["detached_reach_m"] - joined["validated_reach_m"]
    joined["forest_fraction"] = pd.to_numeric(joined[forest_column], errors="coerce")
    joined["valid_fraction"] = pd.to_numeric(joined["valid_fraction_fab"], errors="coerce")
    joined["near_bank_nodata"] = pd.to_numeric(
        joined["near_bank_nodata_fab"], errors="coerce").fillna(0).astype(int).astype(bool)
    joined["censored"] = pd.to_numeric(
        joined["censored_fab"], errors="coerce").fillna(0).astype(int).astype(bool)
    joined["model_excluded"] = pd.to_numeric(joined["analysis_set"], errors="coerce").ne(1)
    joined["low_valid_fraction"] = (joined["valid_fraction"].isna()
                                    | (joined["valid_fraction"] < valid_fraction_min))
    joined["is_station_328"] = joined["station_id"].astype(int).eq(328)
    joined["station_flag"] = np.where(
        joined["is_station_328"], "braided reach, far from main channel", "")
    reasons = joined["exclude_reasons"].fillna("").astype(str)
    joined["missing_validated_reach"] = ~np.isfinite(joined["validated_reach_m"])
    joined["missing_computed_reach"] = ~np.isfinite(joined["computed_reach_m"])
    joined["missing_forest_fraction"] = ~np.isfinite(joined["forest_fraction"])
    joined["main_exclusion_reasons"] = [
        ";".join(_main_exclusion_categories(row, model_reasons))
        for row, model_reasons in zip(joined.itertuples(), reasons)
    ]
    joined["main_include"] = joined["main_exclusion_reasons"].eq("")
    allowed_nodata_reasons = {"near_bank_nodata", "river_not_visible"}
    other_reason = reasons.apply(
        lambda value: any(token and token not in allowed_nodata_reasons
                          for token in value.split(",")))
    joined["near_bank_zero_sensitivity_include"] = (
        joined["near_bank_nodata"] & ~joined["censored"]
        & ~joined["low_valid_fraction"] & ~other_reason
        & np.isfinite(joined["validated_reach_m"])
        & np.isfinite(joined["forest_fraction"])
    )
    joined["computed_reach_m_if_near_bank_zero"] = joined["computed_reach_m"].where(
        ~joined["near_bank_nodata"], 0.0)
    joined["diff_m_if_near_bank_zero"] = (
        joined["computed_reach_m_if_near_bank_zero"] - joined["validated_reach_m"])
    return joined


def largest_errors(main: pd.DataFrame, k: int = 5) -> pd.DataFrame:
    """The k main-set rows with the largest |diff_m|; ties go to lower station_id, left first."""
    ranked = main.assign(abs_diff_m=main["diff_m"].abs(),
                         _side_order=(main["side"] != "left").astype(int))
    ranked = ranked.sort_values(["abs_diff_m", "station_id", "_side_order"],
                                ascending=[False, True, True], kind="mergesort")
    return ranked.head(k).drop(columns="_side_order")


def clip_correlation_interval(interval) -> tuple[tuple[float, float], list[str]]:
    """Clip a correlation interval to [-1, 1] and describe any bound that was moved."""
    lower, upper = (float(interval[0]), float(interval[1]))
    notes = []
    if lower < -1:
        notes.append(f"lower bound was {lower:.3f}")
    if upper > 1:
        notes.append(f"upper bound was {upper:.3f}")
    return (max(lower, -1.0), min(upper, 1.0)), notes


EXCLUSION_CATEGORIES = {
    "censored": "Censored (reach edge not observed)",
    "near_bank_nodata": "Near-bank nodata",
    "low_valid_fraction": "Low valid fraction",
    "model_excluded": "Other model exclusion",
    "missing_validated_reach": "Validated reach unavailable",
    "missing_computed_reach": "Computed reach unavailable",
    "missing_forest_fraction": "Forest fraction unavailable",
}
_FLAG_REASON_TOKENS = {"censored", "near_bank_nodata", "low_valid_fraction"}


def _main_exclusion_categories(row, model_reasons: str) -> list[str]:
    """Every reason a point row is left out of the main set, one token per reason."""
    categories = [flag for flag in ("censored", "near_bank_nodata", "low_valid_fraction")
                  if bool(getattr(row, flag))]
    if bool(row.model_excluded):
        extra = [token for token in model_reasons.split(",")
                 if token and token not in _FLAG_REASON_TOKENS]
        if extra or not categories:
            categories.append("model_excluded" + (f"({','.join(extra)})" if extra else ""))
    categories += [flag for flag in ("missing_validated_reach", "missing_computed_reach",
                                     "missing_forest_fraction") if bool(getattr(row, flag))]
    return categories


def _exclusion_report_lines(excluded: pd.DataFrame) -> list[str]:
    """One bullet per exclusion category; a row with several reasons appears under each."""
    if not len(excluded):
        return ["None."]
    lines = []
    for key, label in EXCLUSION_CATEGORIES.items():
        rows = []
        for row in excluded.itertuples():
            tokens = str(row.main_exclusion_reasons).split(";")
            match = next((t for t in tokens if t == key or t.startswith(key + "(")), None)
            if match is None:
                continue
            detail = match[len(key):]
            rows.append(f"{int(row.station_id)}-{row.side}{' ' + detail if detail else ''}")
        lines.append(f"- {label}: {', '.join(rows) if rows else 'none'}.")
    return lines


def _validation_metrics(frame: pd.DataFrame, metric: str = "diff_m",
                        under_detect_m: float = 20.0) -> dict:
    values = pd.to_numeric(frame[metric], errors="coerce").to_numpy(float)
    values = values[np.isfinite(values)]
    result = {"n": int(values.size), "n_under_detected": int((values < -under_detect_m).sum())}
    if not values.size:
        return {**result, **{name: float("nan") for name in (
            "mean_diff_m", "median_diff_m", "sd_diff_m", "min_diff_m", "max_diff_m",
            "mean_absolute_error_m", "share_within_20m", "share_under_detected")}}
    result.update({
        "mean_diff_m": float(np.mean(values)),
        "median_diff_m": float(np.median(values)),
        "sd_diff_m": float(np.std(values, ddof=1)) if values.size > 1 else float("nan"),
        "min_diff_m": float(np.min(values)),
        "max_diff_m": float(np.max(values)),
        "mean_absolute_error_m": float(np.mean(np.abs(values))),
        "share_within_20m": float(np.mean(np.abs(values) <= 20.0)),
        "share_under_detected": float(np.mean(values < -under_detect_m)),
    })
    return result


def _validation_summary(points: pd.DataFrame, under_detect_m: float = 20.0) -> pd.DataFrame:
    main = points[points["main_include"]].copy()
    rows = []

    def add(subset: str, breakdown: str, category: str, frame: pd.DataFrame,
            metric: str = "diff_m") -> None:
        rows.append({"subset": subset, "breakdown": breakdown, "category": category,
                     **_validation_metrics(frame, metric, under_detect_m)})

    add("main", "all", "All", main)
    for group in ("Low", "Mid", "High"):
        add("main", "forest_group", group, main[main.forest_group == group])
    for date in sorted(points.image_date.astype(str).unique()):
        add("main", "image_date", date, main[main.image_date.astype(str) == date])
    for confidence in (1, 2, 3):
        add("main", "confidence", str(confidence),
            main[pd.to_numeric(main.confidence, errors="coerce") == confidence])

    sensitivities = {
        "without_confidence_3": main[pd.to_numeric(main.confidence, errors="coerce") != 3],
        "near_bank_nodata_as_zero": pd.concat([
            main,
            points[points.near_bank_zero_sensitivity_include]
            .assign(diff_m=lambda f: f["diff_m_if_near_bank_zero"]),
        ], ignore_index=True),
        "without_station_328": main[~main.is_station_328],
    }
    if "review_status" in main.columns:
        sensitivities[REVIEW_SENSITIVITY] = main[main.review_status != NO_VISIBLE_CHANGE]
    for name, frame in sensitivities.items():
        add(name, "all", "All", frame)
    return pd.DataFrame(rows)


REVIEW_STATUSES = ("visible change", "no visible change", "can't tell", "not reviewed")
NO_VISIBLE_CHANGE = "no visible change"
REVIEW_SENSITIVITY = "without_no_visible_change"


def load_validation_review(review: pd.DataFrame, points: pd.DataFrame) -> pd.DataFrame:
    """Attach the imagery-review status to every point; one review row per point."""
    required = {"station_id", "side", "status", "note"}
    missing = sorted(required - set(review.columns))
    if missing:
        raise RuntimeError(f"validation_review.csv is missing columns: {missing}.")
    review = review.copy()
    review["station_id"] = pd.to_numeric(review["station_id"], errors="raise").astype(int)
    review["side"] = review["side"].astype(str).str.strip()
    review["status"] = review["status"].astype(str).str.strip()
    if review.duplicated(["station_id", "side"]).any():
        raise RuntimeError("validation_review.csv has duplicate station_id/side rows.")
    bad = sorted(set(review["status"]) - set(REVIEW_STATUSES))
    if bad:
        raise RuntimeError(f"validation_review.csv has unknown statuses {bad}; "
                           f"allowed: {list(REVIEW_STATUSES)}.")
    keys = points[["station_id", "side"]].astype({"station_id": int})
    merged = keys.merge(review, on=["station_id", "side"], how="outer", indicator=True)
    no_review = merged[merged._merge == "left_only"][["station_id", "side"]]
    no_point = merged[merged._merge == "right_only"][["station_id", "side"]]
    if len(no_review) or len(no_point):
        raise RuntimeError(
            "validation_review.csv must have exactly one row per validation point. "
            f"Points without a review row: {no_review.to_dict('records')}. "
            f"Review rows without a point: {no_point.to_dict('records')}.")
    out = points.merge(
        review.rename(columns={"status": "review_status", "note": "review_note"})
        [["station_id", "side", "review_status", "review_note"]],
        on=["station_id", "side"], how="left", validate="one_to_one")
    out["review_note"] = out["review_note"].fillna("")
    return out


def review_scope_limitation(points: pd.DataFrame, under_detect_m: float) -> str:
    """The limitation sentence on which points were reviewed; refuses if it is untrue."""
    reviewed = points[points.review_status != "not reviewed"]
    errors = reviewed.diff_m.abs()
    if errors.isna().any() or (errors <= under_detect_m).any():
        raise RuntimeError(
            "limitations.md says only points with an error above "
            f"{under_detect_m:g} m were reviewed, but validation_review.csv reviews: "
            f"{reviewed[['station_id', 'side']].to_dict('records')}. Update the text.")
    return (f"Only points with an error above {under_detect_m:g} m were reviewed "
            "against pre-event imagery; smaller errors were not.")


def traced_point_terrain(points: pd.DataFrame, layers: dict, grid: dict) -> pd.DataFrame:
    """Slope, height above z0 and NDVI pre/post at each traced point (nearest pixel)."""
    from .centerline import _sample

    xy =np.column_stack([points.geometry.x.to_numpy(float),
                          points.geometry.y.to_numpy(float)])
    out = points.copy()
    out["traced_slope_deg"] = _sample(layers["slope_fab"], grid, xy)
    out["traced_height_above_z0_m"] = (_sample(layers["z_fab"], grid, xy)
                                       - pd.to_numeric(out["z0_fab"], errors="coerce"))
    out["traced_ndvi_pre"] = _sample(layers["ndvi_pre"], grid, xy)
    out["traced_ndvi_post"] = _sample(layers["ndvi_post"], grid, xy)
    return out


def under_detection_cutoff_table(main: pd.DataFrame, cutoffs, gap_flag: float = 0.10,
                                 min_group_n: int = 8) -> pd.DataFrame:
    """Under-detected share by forest group and the High-minus-Low gap at each cut-off."""
    rows = []
    for cutoff in cutoffs:
        row = {"cutoff_m": float(cutoff)}
        for group in ("Low", "Mid", "High"):
            values = main.loc[main.forest_group == group, "diff_m"].to_numpy(float)
            values = values[np.isfinite(values)]
            under = int((values < -float(cutoff)).sum())
            row[f"{group.lower()}_n"] = int(values.size)
            row[f"{group.lower()}_under_detected"] = under
            row[f"{group.lower()}_share"] = under / values.size if values.size else float("nan")
        row["high_minus_low"] = row["high_share"] - row["low_share"]
        row["flagged"] = bool(row["low_n"] >= min_group_n and row["high_n"] >= min_group_n
                              and row["high_minus_low"] > gap_flag)
        rows.append(row)
    return pd.DataFrame(rows)


def detached_sensitivity(points: pd.DataFrame, under_detect_m: float = 20.0) -> pd.DataFrame:
    """Errors using detached_reach_m, for the main set and the censored rows separately."""
    rows = []
    for subset, frame in (("main", points[points.main_include]),
                          ("censored", points[points.censored])):
        for metric, outcome in (("diff_m", "reach_distance_m"),
                                ("diff_detached_m", "detached_reach_m")):
            usable = frame[np.isfinite(frame[metric]) & np.isfinite(frame["diff_detached_m"])]
            rows.append({"subset": subset, "computed_outcome": outcome,
                         **_validation_metrics(usable, metric, under_detect_m)})
    return pd.DataFrame(rows)


_NUMBER_WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight",
                 "nine", "ten", "eleven", "twelve")


def _number_word(value: int, capital: bool = False) -> str:
    word = _NUMBER_WORDS[value] if 0 <= value < len(_NUMBER_WORDS) else str(value)
    return word[0].upper() + word[1:] if capital else word


def flood_edge_visibility_limitation(main: pd.DataFrame, under_detect_m: float,
                                     pixel_m: float) -> str:
    """The limitation text on where the flood edge is not visible, counts from the data."""
    text = ("Validation is weaker where the flood edge is not visible in imagery: on wide bare "
            "gravel plains, where pre-event bare ground is unknown to the satellite label and a "
            "traced edge may only mark the old plain boundary, and on steep forested slopes, "
            "where a flood that passed under trees may leave no bare deposit.")
    errors = main["diff_m"].abs().to_numpy(float)
    large = np.sort(errors[np.isfinite(errors) & (errors > under_detect_m)])
    if not large.size:
        return text
    near = large[large <= 3 * pixel_m]
    if not near.size:
        return text
    low, high = int(round(near.min())), int(round(near.max()))
    span = f"{low} m" if low == high else f"{low} to {high} m"
    return (f"{text} {_number_word(int(near.size), capital=True)} of the "
            f"{_number_word(int(large.size))} largest main-set errors are {span}, within a few "
            f"{pixel_m:g} m pixels.")


def _wild_validation_association(frame: pd.DataFrame, outcome: str, predictor: str,
                                 cfg) -> dict:
    from .stats import WildClusterBootstrap

    clusters = int(frame.block_id.nunique())
    if clusters < 2:
        raise RuntimeError("The 14.15 wild cluster bootstrap needs at least 2 blocks.")
    bootstrap = WildClusterBootstrap(
        frame, outcome, predictor, predictor, "block_id",
        draws=int(cfg.wild_bootstrap_draws), seed=int(cfg.seed))
    return {
        "estimate": bootstrap.beta_hat,
        "interval": bootstrap.interval(
            alpha=0.05, tolerance_se=float(cfg.ci_bisection_tolerance_se)),
        "p0": bootstrap.p_value(0.0),
        "n": int(len(frame)),
        "n_blocks": clusters,
        "draws": int(cfg.wild_bootstrap_draws),
        "method": "approximate (few validation points); 14.15 wild cluster bootstrap-t",
    }


def _fit_validation_associations(main: pd.DataFrame, cfg) -> dict:
    from scipy.stats import rankdata

    if len(main) < 3 or main.forest_fraction.nunique() < 2:
        raise RuntimeError("Too few eligible validation points or forest values for associations.")
    raw = _wild_validation_association(main, "diff_m", "forest_fraction", cfg)
    ranked = main.copy()
    for source, target in (("forest_fraction", "forest_rank_z"), ("diff_m", "diff_rank_z")):
        ranks = rankdata(ranked[source].to_numpy(float), method="average")
        sd = float(ranks.std(ddof=0))
        if sd == 0:
            raise RuntimeError(f"Spearman correlation is undefined: {source} is constant.")
        ranked[target] = (ranks - ranks.mean()) / sd
    spearman = _wild_validation_association(
        ranked, "diff_rank_z", "forest_rank_z", cfg)
    return {"ols_slope": raw, "spearman_rho": spearman}


def _association_estimates(frame: pd.DataFrame) -> dict:
    if len(frame) < 2 or frame.forest_fraction.nunique() < 2:
        return {"ols_slope": float("nan"), "spearman_rho": float("nan")}
    slope = float(np.polyfit(frame.forest_fraction, frame.diff_m, 1)[0])
    return {"ols_slope": slope,
            "spearman_rho": standardized_rank_slope(frame.forest_fraction, frame.diff_m)}


def _fmt(value, digits: int = 2) -> str:
    return "NA" if pd.isna(value) else f"{float(value):.{digits}f}"


def _markdown_summary(frame: pd.DataFrame, subset: str = "main") -> str:
    selected = frame[(frame.subset == subset) & (frame.breakdown == "forest_group")]
    lines = [
        "| Group | n | Mean diff (m) | Median diff (m) | MAE (m) | Within 20 m | Under-detected |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in selected.itertuples(index=False):
        lines.append(
            f"| {row.category} | {row.n} | {_fmt(row.mean_diff_m)} | "
            f"{_fmt(row.median_diff_m)} | {_fmt(row.mean_absolute_error_m)} | "
            f"{_fmt(row.share_within_20m * 100, 1)}% | "
            f"{_fmt(row.share_under_detected * 100, 1)}% |")
    return "\n".join(lines)


def _write_validation_plot(points: pd.DataFrame, path: Path) -> None:
    import matplotlib.pyplot as plt

    main = points[points.main_include]
    colors = {"Low": "#356b8c", "Mid": "#d28b32", "High": "#4e8962"}
    fig, axis = plt.subplots(figsize=(7, 5))
    for group in ("Low", "Mid", "High"):
        data = main[main.forest_group == group]
        axis.scatter(data.forest_fraction, data.diff_m, label=group,
                     color=colors[group], alpha=0.85, edgecolor="white", linewidth=0.5)
    axis.axhline(0, color="#333333", linewidth=1)
    axis.axhspan(-20, 20, color="#777777", alpha=0.08)
    axis.set_xlabel("Forest fraction F (100 m bank belt)")
    axis.set_ylabel("Computed minus validated reach (m)")
    axis.set_title("Validation reach error by forest fraction")
    axis.legend(title="Forest group", frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def run_validation_comparison(event_key: str = "melamchi2021",
                              config_path: str | Path | None = None) -> list[Path]:
    """Compare hand-digitised reach points with the primary table outcome."""
    import geopandas as gpd

    from .config import PROJECT_ROOT, load_config

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    if event_key != "melamchi2021":
        raise RuntimeError("Independent validation edges are configured for melamchi2021 only.")
    validation_dir = PROJECT_ROOT / cfg.paths.data / event_key / "validation"
    edges_path = validation_dir / "validation_edges.geojson"
    lines_path = validation_dir / "sample_transects.gpkg"
    coverage_path = validation_dir / "sample_coverage.csv"
    table_path = PROJECT_ROOT / cfg.paths.outputs / event_key / "transects.csv"
    review_path = validation_dir / "validation_review.csv"
    required = (edges_path, lines_path, coverage_path, table_path, review_path)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Validation comparison inputs are missing: {missing}.")

    edges = gpd.read_file(edges_path)
    lines = gpd.read_file(lines_path, layer="transects")
    coverage = pd.read_csv(coverage_path)
    table = pd.read_csv(table_path)
    forest_column = f"forest_frac_{int(cfg.belt_width_m)}"
    points = prepare_validation_comparison(
        edges, lines, table, coverage, forest_column=forest_column,
        valid_fraction_min=float(cfg.valid_fraction_min),
        snap_max_m=float(cfg.validation_snap_max_m))
    if len(points) != len(edges):
        raise RuntimeError("Comparison did not retain exactly one row per validation point.")
    point_dates = pd.to_datetime(points.image_date, errors="coerce")
    if point_dates.isna().any() or (point_dates < "2021-08-01").any() \
            or (point_dates > "2022-05-31").any():
        raise RuntimeError("Validation points violate the existing 14.25 capture-date window.")
    points = load_validation_review(pd.read_csv(review_path, dtype={"note": str}), points)

    from . import optical
    from .label import _read_named
    from .profiles import load_layers

    layers, grid = load_layers(event_key, cfg)
    pre = _read_named(PROJECT_ROOT / cfg.paths.data / event_key / "rasters" / "s2_pre.tif",
                      ("B4", "B8"), cfg.crs, grid)[0]
    layers["ndvi_pre"] = optical.ndvi(pre["B8"], pre["B4"])
    points = traced_point_terrain(points, layers, grid)
    del layers, pre

    under_detect_m = float(cfg.validation_under_detect_m)
    cutoffs = [float(c) for c in cfg.validation_under_detect_sensitivity_m]
    all_blocks = int(points.block_id.nunique())
    main = points[points.main_include].copy()
    if main.empty:
        raise RuntimeError("No validation points meet the main comparison inclusion rules.")
    associations = _fit_validation_associations(main, cfg)
    summary = _validation_summary(points, under_detect_m)
    cutoff_table = under_detection_cutoff_table(main, cutoffs)
    detached_table = detached_sensitivity(points, under_detect_m)
    review_frame = main[main.review_status != NO_VISIBLE_CHANGE]
    review_associations = _association_estimates(review_frame)
    confidence_sensitivity = summary[
        (summary.subset == "without_confidence_3") & (summary.breakdown == "all")].iloc[0]
    nodata_sensitivity = summary[
        (summary.subset == "near_bank_nodata_as_zero") & (summary.breakdown == "all")].iloc[0]
    station_sensitivity = summary[
        (summary.subset == "without_station_328") & (summary.breakdown == "all")].iloc[0]
    sensitivity_frames = {
        "without_confidence_3": main[pd.to_numeric(main.confidence, errors="coerce") != 3],
        "near_bank_nodata_as_zero": pd.concat([
            main,
            points[points.near_bank_zero_sensitivity_include]
            .assign(diff_m=lambda frame: frame["diff_m_if_near_bank_zero"]),
        ], ignore_index=True),
        "without_station_328": main[~main.is_station_328],
    }
    sensitivity_associations = {
        name: _association_estimates(frame) for name, frame in sensitivity_frames.items()
    }

    groups = summary[(summary.subset == "main") & (summary.breakdown == "forest_group")]
    low = groups[groups.category == "Low"].iloc[0]
    high = groups[groups.category == "High"].iloc[0]
    high_low_gap = (float(high.share_under_detected - low.share_under_detected)
                    if low.n and high.n else float("nan"))
    forest_bias_flag = bool(low.n >= 8 and high.n >= 8 and high_low_gap > 0.10)
    excluded = points[~points.main_include]
    top_errors = largest_errors(main, k=5)
    low_coverage = coverage[coverage.forest_group.astype(str).str.lower() == "low"]
    low_skipped = int((low_coverage.status.astype(str).str.lower() == "skipped").sum())
    low_total = int(len(low_coverage))
    unique_stations = int(points.station_id.nunique())

    from .config import validation_results_dir
    results_dir = validation_results_dir(cfg, event_key)
    results_dir.mkdir(parents=True, exist_ok=True)
    per_point_path = results_dir / "validation_per_point.csv"
    summary_path = results_dir / "validation_summary.csv"
    plot_path = results_dir / "validation_plot.png"
    limitations_path = results_dir / "limitations.md"
    report_path = results_dir / "REPORT.md"

    output_columns = [
        "station_id", "side", "image_date", "confidence", "forest_fraction",
        "forest_group", "block_id", "validated_reach_m", "computed_reach_m", "diff_m",
        "valid_fraction", "near_bank_nodata", "censored", "analysis_set",
        "exclude_reasons", "model_excluded", "low_valid_fraction",
        "missing_validated_reach", "missing_computed_reach", "missing_forest_fraction",
        "main_include",
        "main_exclusion_reasons", "computed_reach_m_if_near_bank_zero",
        "diff_m_if_near_bank_zero", "near_bank_zero_sensitivity_include",
        "is_station_328", "station_flag", "note",
        "traced_slope_deg", "traced_height_above_z0_m", "traced_ndvi_pre", "traced_ndvi_post",
        "detached_reach_m", "validated_minus_detached_m", "diff_detached_m",
        "review_status", "review_note",
    ]
    cutoff_path = results_dir / "validation_cutoff_sensitivity.csv"
    detached_path = results_dir / "validation_detached_sensitivity.csv"
    points[output_columns].to_csv(per_point_path, index=False)
    summary.to_csv(summary_path, index=False)
    cutoff_table.to_csv(cutoff_path, index=False)
    detached_table.to_csv(detached_path, index=False)
    _write_validation_plot(points, plot_path)

    ols = associations["ols_slope"]
    rho = associations["spearman_rho"]
    ols_line = (
        f"- OLS slope of diff on forest fraction: **{_fmt(ols['estimate'], 3)} m per unit F**; "
        f"95% interval {_fmt(ols['interval'][0], 3)} to {_fmt(ols['interval'][1], 3)} "
        "m per unit F (approximate, few validation points).")
    intervals_with_zero = []
    if ols["interval"][0] <= 0 <= ols["interval"][1]:
        intervals_with_zero.append("OLS slope")
    if all_blocks < 8:
        rho_line = (f"- Spearman rho = {_fmt(rho['estimate'], 3)}; too few clusters "
                    f"({all_blocks}); no interval.")
    else:
        (rho_low, rho_high), clip_notes = clip_correlation_interval(rho["interval"])
        clip_text = ("; clipped to [-1, 1], the raw bootstrap " + " and ".join(clip_notes)
                     if clip_notes else "")
        rho_line = (f"- Spearman rho = **{_fmt(rho['estimate'], 3)}**; 95% interval "
                    f"{_fmt(rho_low, 3)} to {_fmt(rho_high, 3)} (approximate, few validation "
                    f"points{clip_text}).")
        if rho_low <= 0 <= rho_high:
            intervals_with_zero.append("Spearman rho")
    if intervals_with_zero:
        subject = (" and ".join(intervals_with_zero) + " intervals are"
                   if len(intervals_with_zero) > 1 else intervals_with_zero[0] + " interval is")
        trend_text = (
            f"The {subject} wide and include{'s' if len(intervals_with_zero) == 1 else ''} 0, "
            "so these data cannot rule out either no trend or a strong negative trend in the "
            "error with forest fraction. This is not evidence that there is no forest bias.")
    else:
        trend_text = ("Neither interval includes 0. These are approximate intervals from few "
                      "validation points and describe the validation error only.")
    group_counts = (f"High {_fmt(high.share_under_detected * 100, 1)}% "
                    f"({int(high.n_under_detected)} of {int(high.n)}) vs Low "
                    f"{_fmt(low.share_under_detected * 100, 1)}% "
                    f"({int(low.n_under_detected)} of {int(low.n)})")
    bias_text = (
        f"Under-detected share (computed reach more than {under_detect_m:g} m shorter than the "
        f"traced edge): {group_counts}. "
        + ("Cannot assess the High-versus-Low difference because at least one group has n < 8."
           if low.n < 8 or high.n < 8 else
           f"High minus Low: {high_low_gap * 100:.1f} percentage points. "
           + ("Flagged (>10 percentage points); the main result may be biased toward finding "
              "protection." if forest_bias_flag else "Not flagged (>10 percentage points).")))
    cutoff_lines = [
        "| Cut-off (m) | Low | Mid | High | High minus Low (points) | Flag (>10 points) |",
        "|---:|---:|---:|---:|---:|---|",
        *[f"| {row.cutoff_m:g}{' (primary)' if row.cutoff_m == under_detect_m else ''} | "
          + " | ".join(f"{_fmt(getattr(row, g + '_share') * 100, 1)}% "
                       f"({getattr(row, g + '_under_detected')} of {getattr(row, g + '_n')})"
                       for g in ("low", "mid", "high"))
          + f" | {_fmt(row.high_minus_low * 100, 1)} | {'yes' if row.flagged else 'no'} |"
          for row in cutoff_table.itertuples()],
    ]
    flag_values = set(cutoff_table.flagged)

    def _and_list(values) -> str:
        items = [f"{c:g} m" for c in values]
        return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]

    if len(flag_values) > 1:
        raised = _and_list(cutoff_table.cutoff_m[cutoff_table.flagged])
        not_raised = _and_list(cutoff_table.cutoff_m[~cutoff_table.flagged])
        cutoff_note = (f"The flag is sensitive to the cut-off: it is raised at {raised}, "
                       f"and not at {not_raised}.")
    else:
        cutoff_note = ("The flag is the same at every cut-off shown "
                       f"({'raised' if flag_values == {True} else 'not raised'}), but the gap "
                       "itself changes with the cut-off.")
    cutoff_text = (f"The {under_detect_m:g} m cut-off was fixed before the results were seen "
                   f"and stays primary. {cutoff_note}")

    def detached_row(subset: str, outcome: str, label: str) -> str:
        row = detached_table[(detached_table.subset == subset)
                             & (detached_table.computed_outcome == outcome)].iloc[0]
        return (f"| {label} | {row.n} | {_fmt(row.mean_diff_m)} | {_fmt(row.median_diff_m)} | "
                f"{_fmt(row.mean_absolute_error_m)} | {_fmt(row.share_under_detected * 100, 1)}% "
                f"({row.n_under_detected} of {row.n}) |")

    detached_lines = [
        "| Set and computed outcome | n | Mean diff (m) | Median diff (m) | MAE (m) | "
        "Under-detected |",
        "|---|---:|---:|---:|---:|---:|",
        detached_row("main", "reach_distance_m", "Main set, reach_distance_m (primary)"),
        detached_row("main", "detached_reach_m", "Main set, detached_reach_m"),
        detached_row("censored", "reach_distance_m", "Censored rows, reach_distance_m"),
        detached_row("censored", "detached_reach_m", "Censored rows, detached_reach_m"),
    ]
    censored_rows = points[points.censored].sort_values(["station_id", "side"])
    censored_lines = [
        "| Station | Side | Forest group | Validated (m) | reach_distance_m | detached_reach_m |",
        "|---:|---|---|---:|---:|---:|",
        *[f"| {int(r.station_id)} | {r.side} | {r.forest_group} | {_fmt(r.validated_reach_m, 1)} "
          f"| {_fmt(r.computed_reach_m, 1)} | {_fmt(r.detached_reach_m, 1)} |"
          for r in censored_rows.itertuples()],
    ]
    main_row = summary[(summary.subset == "main") & (summary.breakdown == "all")].iloc[0]
    review_row = summary[(summary.subset == REVIEW_SENSITIVITY)
                         & (summary.breakdown == "all")].iloc[0]
    review_cut = under_detection_cutoff_table(review_frame, [under_detect_m]).iloc[0]
    main_cut = cutoff_table[cutoff_table.cutoff_m == under_detect_m].iloc[0]
    main_estimates = _association_estimates(main)
    review_excluded = main[main.review_status == NO_VISIBLE_CHANGE]

    def review_line(label, row, cut, estimates) -> str:
        groups = " | ".join(
            f"{getattr(cut, g + '_under_detected')} of {getattr(cut, g + '_n')}"
            for g in ("low", "mid", "high"))
        return (f"| {label} | {row.n} | {_fmt(row.mean_diff_m)} | {groups} | "
                f"{_fmt(cut.high_minus_low * 100, 1)} | {'on' if cut.flagged else 'off'} | "
                f"{_fmt(estimates['ols_slope'], 3)} | {_fmt(estimates['spearman_rho'], 3)} |")

    review_lines = [
        f"Under-detected = computed reach more than {under_detect_m:g} m shorter than the "
        "traced edge.",
        "",
        "| Set | n | Mean diff (m) | Low under-detected | Mid under-detected | "
        "High under-detected | High minus Low (points) | Flag (>10 points) | "
        "OLS slope (m per unit F) | Spearman rho |",
        "|---|---:|---:|---:|---:|---:|---:|---|---:|---:|",
        review_line("Main", main_row, main_cut, main_estimates),
        review_line("Without \"no visible change\"", review_row, review_cut, review_associations),
    ]
    review_note = (
        "Points reviewed as \"no visible change\" and left out of the second row: "
        + (", ".join(f"{int(r.station_id)}-{r.side}" for r in review_excluded.itertuples())
           or "none")
        + f". Review status counts over all {len(points)} points: "
        + ", ".join(f"{status} {int((points.review_status == status).sum())}"
                    for status in REVIEW_STATUSES)
        + ". The slope and rho in this table are point estimates; the bootstrap intervals "
        "are for the main set only.")
    review_caveat = (
        "No visible change can mean the traced edge is not a flood edge, or that the flood left "
        "no visible trace (for example under trees). Imagery cannot separate these, so this row "
        "does not show that the bias is absent. The main row stays primary.")
    review_scope = review_scope_limitation(points, under_detect_m)
    largest_error_lines = [
        "| Station | Side | Forest group | Confidence | diff (m) |",
        "|---:|---|---|---:|---:|",
        *[f"| {int(row.station_id)} | {row.side} | {row.forest_group} | "
          f"{_fmt(row.confidence, 0)} | {_fmt(row.diff_m)} |"
          for row in top_errors.itertuples()],
    ]
    excluded_lines = _exclusion_report_lines(excluded)
    sensitivity_text = [
        f"| Without confidence 3 | {confidence_sensitivity.n} | "
        f"{_fmt(confidence_sensitivity.mean_diff_m)} | {_fmt(confidence_sensitivity.mean_absolute_error_m)} | "
        f"{_fmt(confidence_sensitivity.share_under_detected * 100, 1)}% | "
        f"{_fmt(sensitivity_associations['without_confidence_3']['ols_slope'], 3)} | "
        f"{_fmt(sensitivity_associations['without_confidence_3']['spearman_rho'], 3)} |",
        f"| Near-bank nodata set to zero | {nodata_sensitivity.n} | "
        f"{_fmt(nodata_sensitivity.mean_diff_m)} | {_fmt(nodata_sensitivity.mean_absolute_error_m)} | "
        f"{_fmt(nodata_sensitivity.share_under_detected * 100, 1)}% | "
        f"{_fmt(sensitivity_associations['near_bank_nodata_as_zero']['ols_slope'], 3)} | "
        f"{_fmt(sensitivity_associations['near_bank_nodata_as_zero']['spearman_rho'], 3)} |",
        f"| Without station 328 | {station_sensitivity.n} | "
        f"{_fmt(station_sensitivity.mean_diff_m)} | {_fmt(station_sensitivity.mean_absolute_error_m)} | "
        f"{_fmt(station_sensitivity.share_under_detected * 100, 1)}% | "
        f"{_fmt(sensitivity_associations['without_station_328']['ols_slope'], 3)} | "
        f"{_fmt(sensitivity_associations['without_station_328']['spearman_rho'], 3)} |",
    ]
    date_table = summary[(summary.subset == "main") & (summary.breakdown == "image_date")]
    date_lines = [f"- {row.category}: n={row.n}, mean diff {_fmt(row.mean_diff_m)} m, "
                  f"MAE {_fmt(row.mean_absolute_error_m)} m, under-detected "
                  f"{_fmt(row.share_under_detected * 100, 1)}%."
                  for row in date_table.itertuples()]
    confidence_table = summary[(summary.subset == "main") & (summary.breakdown == "confidence")]
    confidence_lines = [f"- Confidence {row.category}: n={row.n}, mean diff "
                        f"{_fmt(row.mean_diff_m)} m, MAE {_fmt(row.mean_absolute_error_m)} m."
                        for row in confidence_table.itertuples()]
    low_skip_note = (
        f"- Low-forest coverage: {low_skipped} of {low_total} sampled stations were skipped. "
        "The coverage file does not record station-specific reasons. The digitising instructions "
        "say to skip when cloud, shadow, or trees hide the edge, so the imagery protocol is clear "
        "but the cause for each skipped station cannot be determined from the saved data."
    )
    limitations_text = "\n".join([
        "# Validation limitations",
        "",
        low_skip_note,
        "- Station 328 is retained and flagged as `braided reach, far from main channel`; its "
        "original note about the imagery dot being a couple of metres off is preserved.",
        f"- The {len(points)} validation points were hand-marked from Esri Wayback imagery. "
        f"The recorded capture dates are {', '.join(sorted(points.image_date.astype(str).unique()))}; "
        "the Wayback release is 2024-11-18.",
        f"- The points cover {unique_stations} stations from the 40-station sample. "
        f"The wild-bootstrap analysis has {int(main.block_id.nunique())} blocks among its "
        "eligible points; intervals are labelled approximate because the validation sample is small.",
        "- This is a transect-distance comparison against the primary combined-label/FABDEM "
        "table outcome, not a full-raster intersection-over-union assessment. Hand-digitised "
        "edges and the computed label both have measurement uncertainty.",
        "- Skipped low-forest stations weaken coverage of that group; comparisons with any group "
        "having fewer than 8 eligible points are descriptive only.",
        "- " + flood_edge_visibility_limitation(main, under_detect_m, float(cfg.grid_m)),
        "- " + review_scope,
    ]) + "\n"
    report_text = "\n".join([
        "# Validation comparison report",
        "",
        "## What was run",
        "",
        f"Compared {len(points)} hand-digitised edge points from {unique_stations} stations "
        "with the primary combined-label/FABDEM reach-distance table output. The signed error is "
        "computed reach minus validated reach; negative values mean under-detection. The main "
        f"summary includes {len(main)} points in {int(main.block_id.nunique())} 3 km blocks.",
        "",
        "Files produced: `validation_per_point.csv`, `validation_summary.csv`, "
        "`validation_cutoff_sensitivity.csv`, `validation_detached_sensitivity.csv`, "
        "`validation_plot.png`, `limitations.md`, and this report.",
        "",
        "## Main results",
        "",
        _markdown_summary(summary),
        "",
        "### Main result beside the imagery-review sensitivity",
        "",
        *review_lines,
        "",
        review_caveat,
        "",
        review_note,
        "",
        "### Trend with forest fraction",
        "",
        ols_line,
        rho_line,
        "",
        "Both intervals use the 14.15 null-imposed wild cluster bootstrap-t (CR1, Rademacher, "
        f"{int(cfg.wild_bootstrap_draws)} draws, 3 km block IDs).",
        "",
        trend_text,
        "",
        "### Forest-bias check",
        "",
        bias_text,
        "",
        "### Under-detection cut-off sensitivity",
        "",
        cutoff_text,
        "",
        *cutoff_lines,
        "",
        "### Date and confidence breakdowns",
        "",
        *date_lines,
        *confidence_lines,
        "",
        "### Largest errors",
        "",
        f"The {len(top_errors)} main-set points with the largest absolute error. All are kept "
        "in every summary.",
        "",
        *largest_error_lines,
        "",
        "Station 328 is kept in every main summary and flagged descriptively as `braided reach, "
        "far from main channel`; the sensitivity run below repeats the summary without it.",
        "",
        "### Exclusions",
        "",
        f"Excluded from main numbers ({len(excluded)} point rows). A row with several "
        "reasons is listed under each:",
        "",
        *excluded_lines,
        "",
        "### Sensitivity checks",
        "",
        "| Sensitivity | n | Mean diff (m) | MAE (m) | Under-detected | OLS slope | Spearman rho |",
        "|---|---:|---:|---:|---:|---:|---:|",
        *sensitivity_text,
        "",
        "### Detached reach sensitivity (not part of the main summary)",
        "",
        "Errors when the computed reach is `detached_reach_m` (furthest connected impact on the "
        "transect, gaps ignored, 14.6) instead of `reach_distance_m`. diff is still computed "
        "minus validated. Censored rows are shown separately; they are not in the main set.",
        "",
        *detached_lines,
        "",
        *censored_lines,
        "",
        "Negative diff means the computed reach is shorter than the digitised reach. These "
        "descriptive validation results do not establish that forest causes a difference.",
        "",
        "## Still needed",
        "",
        "Review the station-level comparison and limitations, then complete the study's separate "
        "statistical, placebo, and reporting steps. Validation is one check of the impact label; "
        "it does not by itself establish the forest-reach association.",
        "",
        "See `limitations.md` for additional data and coverage limitations.",
    ]) + "\n"
    limitations_path.write_text(limitations_text, encoding="utf-8")
    report_path.write_text(report_text, encoding="utf-8")
    return [per_point_path, summary_path, cutoff_path, detached_path, plot_path,
            limitations_path, report_path]
