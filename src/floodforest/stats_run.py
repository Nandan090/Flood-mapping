"""The stats step: Model A and B with the 14.15 intervals, every sensitivity run in
8.8 and section 14, the placebo through the same code, and the verdict (8.10).

Run: python -m floodforest.cli stats --event melamchi2021
Writes outputs/<event>/stats/ only.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from . import stats as S
from .config import (PROJECT_ROOT, load_config, placebo_key, rebuild_variants,
                     validation_results_dir, variant_key)
from .floods import flood_records, report_lines

LABEL_DYNAMIC = {"river_not_visible", "near_bank_nodata", "low_valid_fraction", "belt_no_data"}
NULL_PUSH = "pushes toward null by construction"


# ---------------------------------------------------------------------------
# Building the rows of one run
# ---------------------------------------------------------------------------

def analysis_rows(table: pd.DataFrame, labels: pd.DataFrame, cfg, label: str = "combined",
                  dem: str = "fab", belt: int | None = None,
                  ignore_reasons: frozenset = frozenset()) -> pd.DataFrame:
    """The table with one label's outcomes and that label's analysis set (in_set).

    Label- and DEM-specific reasons (14.20, 14.24, 7.2 valid fraction) are rebuilt from
    transects_labels.csv; belt_no_data follows the belt width; the other reason codes
    are the table's. For the combined label and FABDEM this equals analysis_set.
    """
    belt = belt or int(cfg.belt_width_m)
    part = labels[(labels.label == label) & (labels.dem == dem)]
    keep = ["reach_distance_m", "reach_height_m", "overtopped", "censored", "reach_height_unc_m",
            "detached_reach_m", "near_bank_nodata", "valid_fraction", "river_visible_post"]
    part = part[["station_id", "side"] + keep].rename(
        columns={k: f"L_{k}" for k in keep})
    t = table.merge(part, on=["station_id", "side"], how="left", validate="one_to_one")
    for k in ("reach_distance_m", "reach_height_m", "overtopped", "censored",
              "reach_height_unc_m", "detached_reach_m"):
        t[f"{k}_{dem}"] = t[f"L_{k}"]
    static = t["exclude_reasons"].fillna("").map(
        lambda v: {r for r in v.split(",") if r} - LABEL_DYNAMIC - set(ignore_reasons))
    t["in_set"] = ((t["in_reach"] == 1) & (static.map(len) == 0)
                   & (t["L_river_visible_post"] == 1) & (t["L_near_bank_nodata"] != 1)
                   & (t["L_valid_fraction"] >= cfg.valid_fraction_min)
                   & t[f"forest_frac_{belt}"].notna())
    return t


def include_as_zero(table: pd.DataFrame, reason: str, dem: str = "fab") -> pd.DataFrame:
    """In-reach rows left out only for `reason` join the set with reach 0 (14.20, 14.24)."""
    t = table.copy()
    only = (t["in_reach"] == 1) & (t["exclude_reasons"].fillna("") == reason)
    t.loc[only, "analysis_set"] = 1
    for k in ("reach_distance_m", "reach_height_m", "overtopped", "censored"):
        t.loc[only, f"{k}_{dem}"] = 0
    t.attrs["added_rows"] = int(only.sum())
    return t


def bank_top_slope(cfg, event_key: str, table: pd.DataFrame) -> pd.Series:
    """FABDEM slope (deg) at each row's bank top: pre-event terrain, not the outcome."""
    import geopandas as gpd
    from rasterio.transform import rowcol

    from .config import rasters_dir
    from .label import _read_named
    from .sar import terrain_slope_aspect

    dem, grid = _read_named(rasters_dir(cfg, event_key) / "fabdem.tif", ("elevation_m",), cfg.crs)
    slope = terrain_slope_aspect(dem["elevation_m"], grid["pixel_m"])[0]
    lines = gpd.read_file(PROJECT_ROOT / cfg.paths.outputs / event_key / "transects.gpkg",
                          layer="transects")
    lines = lines.set_index(["station_id", "side"])
    out = []
    for row in table.itertuples(index=False):
        c = np.asarray(lines.loc[(row.station_id, row.side)].geometry.coords)
        n = (c[1] - c[0]) / np.linalg.norm(c[1] - c[0])
        x, y = c[0] + n * row.bank_top_from_station_m
        r, cc = rowcol(grid["transform"], x, y)
        ok = 0 <= r < grid["height"] and 0 <= cc < grid["width"]
        out.append(float(slope[r, cc]) if ok else np.nan)
    return pd.Series(out, index=table.index)


def _const(df: pd.DataFrame, col: str) -> bool:
    return df[col].nunique(dropna=True) < 2


def rhs_for(df: pd.DataFrame, cfg, with_sections: bool) -> str:
    """Model A right-hand side; terms that are constant in these rows are left out
    (a constant column cannot be estimated beside the intercept)."""
    extra = tuple(t for t in S.MODEL_A_EXTRA if not _const(df, t))
    rhs = S.model_a_rhs("forest_frac", int(cfg.spline_df),
                        with_sections and not _const(df, "river_section"), extra)
    if _const(df, "side"):
        rhs = rhs.replace(" + C(side)", "")
    return rhs


def b_rhs_for(pairs: pd.DataFrame) -> str:
    terms = ["d_forest"] + [f"d_{n}" for n in S.MODEL_B_DIFFS if not _const(pairs, f"d_{n}")]
    return " + ".join(terms)


# ---------------------------------------------------------------------------
# Fitting one run
# ---------------------------------------------------------------------------

def _wild(df, outcome, rhs, term, cfg, percentile=False) -> dict:
    try:
        res = S.wild_cluster_result(df, outcome, rhs, term, "block_id",
                                    int(cfg.wild_bootstrap_draws), int(cfg.seed),
                                    float(cfg.ci_bisection_tolerance_se),
                                    alpha=float(cfg.interval_alpha))
    except Exception as exc:  # reported in the table, never hidden
        return {"error": f"{type(exc).__name__}: {exc}"}
    res["interval"] = tuple(float(v) for v in res["interval"])
    if percentile:
        res["percentile"] = S.pairs_percentile_interval(
            df, outcome, rhs, term, int(cfg.bootstrap_reps), int(cfg.seed))
    return res


def _logit(df, rhs, cfg) -> dict:
    import patsy
    k = patsy.dmatrix(rhs, df, return_type="dataframe").shape[1]
    events = int(df["y_over"].sum())
    fewer = min(events, len(df) - events)
    if fewer < k:
        return {"error": f"not estimable: {fewer} rows in the smaller outcome class, fewer than "
                         f"the {k} model terms"}
    try:
        res = S.logistic_cluster_interval(df, "y_over", rhs, "forest_frac", "block_id",
                                          level=1 - float(cfg.interval_alpha))
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}
    if not all(np.isfinite([res["odds_ratio"], *res["or_interval"]])):
        return {"error": "not estimable (odds ratio or interval not finite)"}
    return res


def fit_run(df: pd.DataFrame, cfg, with_sections: bool, trim: bool = True,
            b_df: pd.DataFrame | None = None, percentile: bool = False,
            aft: bool = False) -> dict:
    """Model A (trimmed by common support unless trim=False) and Model B (never trimmed).

    Model A uses rows with every Model A value (a_ok); Model B pairs rows with Model B's
    own variables (b_ok), so a missing Model A control does not cost a Model B pair.
    """
    a_all = df[df["a_ok"]]
    rhs = rhs_for(a_all, cfg, with_sections)
    controls = rhs.split(" + ", 1)[1]
    out = {"rows": int(len(df)), "rhs_a": rhs, "dropped_missing": df.attrs.get("dropped_missing", 0)}
    a_df = a_all
    if trim:
        try:
            cs = S.common_support(a_all, controls, float(cfg.forest_high_threshold),
                                  float(cfg.forest_high_min_share),
                                  float(cfg.overlap_min_row_fraction), int(cfg.overlap_min_segments))
            a_df = a_all[cs["keep"]]
            out["support"] = {k: cs[k] for k in ("rule", "bounds", "rows", "kept", "dropped",
                                                  "kept_segments", "too_thin")}
            out["_cs"] = cs
        except Exception as exc:
            out["support"] = {"error": f"{type(exc).__name__}: {exc}; Model A fitted untrimmed"}
    out["a_rows"] = int(len(a_df))
    out["_rhs_a_fitted"] = rhs_for(a_df, cfg, with_sections)   # after trimming: what is fitted
    out["a_stations"] = int(a_df.station_id.nunique())
    out["a_segments"] = int(a_df.segment_id.nunique())
    out["a_blocks"] = int(a_df.block_id.nunique())
    out["A"] = {y: _wild(a_df, y, rhs_for(a_df, cfg, with_sections), "forest_frac", cfg, percentile)
                for y in ("y_height", "y_logdist")}
    out["A"]["y_over"] = _logit(a_df, rhs_for(a_df, cfg, with_sections), cfg)
    # SD(y) of Model A's own rows: used for Model A E-values and, unchanged, for 14.10.
    out["sd"] = {y: float(a_df[y].std(ddof=1)) for y in ("y_height", "y_logdist")}
    # 14.9: prevalence of overtopped in the full analysis set of this run.
    out["prevalence"] = float(df["y_over"].dropna().mean())
    out["iqr_F"] = float(a_df.forest_frac.quantile(0.75) - a_df.forest_frac.quantile(0.25))
    b_src = df if b_df is None else b_df
    pairs = S.paired_frame(b_src[b_src["b_ok"]])
    out["b_info"] = S.model_b_informative(pairs, int(cfg.model_b_min_stations),
                                          float(cfg.model_b_min_forest_diff))
    out["b_blocks"] = int(pairs.block_id.nunique())
    out["b_segments"] = int(pairs.segment_id.nunique())
    # SD of Model B's own outcome (y_L - y_R), for Model B E-values (8.7).
    out["sd_b"] = {y: float(pairs[f"d_{y}"].std(ddof=1)) for y in ("y_height", "y_logdist")}
    out["_pairs"] = pairs
    out["B"] = {y: _wild(pairs, f"d_{y}", b_rhs_for(pairs), "d_forest", cfg, percentile)
                for y in ("y_height", "y_logdist")}
    if aft:
        out["aft"] = {}
        for y in ("y_dist", "y_height"):
            try:
                out["aft"][y] = S.aft_conditional(a_df, y, rhs_for(a_df, cfg, with_sections))
            except Exception as exc:
                out["aft"][y] = {"error": f"{type(exc).__name__}: {exc}"}
    out["_a_df"] = a_df
    return out


def flat_rows(name: str, group: str, note: str, run: dict) -> list[dict]:
    rows = []
    for model, part in (("A", run.get("A", {})), ("B", run.get("B", {}))):
        for outcome, r in part.items():
            base = {"run": name, "group": group, "note": note, "model": model, "outcome": outcome,
                    "rows": run["a_rows"] if model == "A" else run["b_info"]["stations"],
                    "segments": run["a_segments"] if model == "A" else run["b_segments"],
                    "blocks": run["a_blocks"] if model == "A" else run["b_blocks"]}
            if "error" in r:
                rows.append({**base, "error": r["error"]})
            elif outcome == "y_over":
                rows.append({**base, "beta": r["beta"], "se": r["se_cluster"],
                             "odds_ratio": r["odds_ratio"], "lo": r["or_interval"][0],
                             "hi": r["or_interval"][1], "interval_scale": "odds ratio",
                             "method": r["method"]})
            else:
                rows.append({**base, "beta": r["beta"], "se": r["se_cluster"], "p0": r["p0"],
                             "lo": r["interval"][0], "hi": r["interval"][1],
                             "interval_scale": "beta", "method": r["method"]})
    return rows


# ---------------------------------------------------------------------------
# The step
# ---------------------------------------------------------------------------

def run_stats(event_key: str, config_path: str | Path | None = None,
              rebuild: bool = False) -> list[Path]:
    from .placebo import forest_terciles, river_not_visible_as_zero
    from .sensitivity import build_variants

    cfg = load_config(config_path or PROJECT_ROOT / "config.yaml")
    event = getattr(cfg.events, event_key)
    out_root = PROJECT_ROOT / cfg.paths.outputs
    out_dir = out_root / event_key / "stats"
    out_dir.mkdir(parents=True, exist_ok=True)
    sections = bool(event.section_breaks_km)

    def load(key):
        return (pd.read_csv(out_root / key / "transects.csv"),
                pd.read_csv(out_root / key / "transects_labels.csv"))

    T, L = load(event_key)
    base = analysis_rows(T, L, cfg)
    if not (base["in_set"] == (T["analysis_set"] == 1)).all():
        raise RuntimeError("Rebuilt analysis set differs from analysis_set in transects.csv.")
    forest = f"forest_frac_{int(cfg.belt_width_m)}"

    def frame(t, mask=None, forest_col=forest, dem="fab", outcomes=None):
        rows = t[t["in_set"]] if mask is None else t[mask]
        return S.model_frame(rows, forest_col, dem, outcomes)

    print("Primary run ...")
    primary = fit_run(frame(base), cfg, sections, percentile=True, aft=True)
    results = flat_rows("primary", "primary", "combined label, FABDEM, 100 m belt, trimmed A",
                        primary)

    specs = []

    def spec(name, group, note, df, **kw):
        specs.append((name, group, note, df, kw))

    for lab, ref in (("optical_only", "8.8"), ("sar_only", "8.8"), ("sar_strict", "14.12"),
                     ("sar_near", "14.17")):
        t = analysis_rows(T, L, cfg, label=lab)
        spec(f"label {lab}", "label", ref, frame(t))
    t = analysis_rows(T, L, cfg, dem="cop")
    spec("DEM Copernicus GLO-30", "DEM", "8.8; outcomes and DEM controls from GLO-30",
         frame(t, dem="cop"))
    for w in (int(v) for v in cfg.belt_widths_m):
        if w == int(cfg.belt_width_m):
            continue
        t = analysis_rows(T, L, cfg, belt=w)
        mask = t["in_set"] & ((t.usable_length_m >= w) if w > cfg.belt_width_m else True)
        spec(f"belt {w} m", "belt width",
             "8.8" + (f"; rows with usable_length_m >= {w} (14.21)" if w > cfg.belt_width_m else ""),
             frame(t, mask, forest_col=f"forest_frac_{w}"))
    spec("forest WorldCover", "forest definition", "8.8: WorldCover class 10 share",
         frame(base, base["in_set"] & base[f"forest_frac_worldcover_{int(cfg.belt_width_m)}"].notna(),
               forest_col=f"forest_frac_worldcover_{int(cfg.belt_width_m)}"))
    disagree = base.loc[base["in_set"], f"forest_definition_disagree_frac_{int(cfg.belt_width_m)}"].mean()
    spec("forest Dynamic World only", "forest definition",
         f"14.14: no canopy-height condition; the two rules disagree on {disagree:.1%} of belt "
         "pixels, so this run nearly equals the dynamic_world-only run",
         frame(base, base["in_set"] & base[f"forest_frac_dw_only_{int(cfg.belt_width_m)}"].notna(),
               forest_col=f"forest_frac_dw_only_{int(cfg.belt_width_m)}"))
    spec("forest_source dynamic_world only", "forest definition", "14.23 a",
         frame(base, base["in_set"] & (base[f"forest_source_{int(cfg.belt_width_m)}"] == "dynamic_world")))
    slope_top = bank_top_slope(cfg, event_key, base)
    smax = float(cfg.steep_bank_top_slope_max_deg)
    umax = float(cfg.steep_unc_max_m)
    spec(f"steep: bank-top slope <= {smax:g} deg", "steep terrain",
         f"8.8 (user decision): pre-event terrain, {cfg.grid_m:g} x tan(bank-top slope) <= {umax:g} m",
         frame(base, base["in_set"] & (slope_top <= smax)))
    unc = base["reach_height_unc_m_fab"]
    spec(f"steep: unc <= {umax:g} m plus all zero-reach rows", "steep terrain",
         "8.8 variant (user decision): keeps every zero-reach row",
         frame(base, base["in_set"] & ((base["reach_distance_m_fab"] == 0) | (unc <= umax))))
    spec(f"steep: literal unc <= {umax:g} m", "steep terrain",
         "selects rows by outcome (conflicts with 14.5); not a headline result",
         frame(base, base["in_set"] & (unc <= umax)))
    tmax = float(cfg.terrace_bank_slope_max_deg)
    spec(f"terraces: bank slope <= {tmax:g} deg", "terrain", "8.8",
         frame(base, base["in_set"] & (base["bank_slope_deg_fab"] <= tmax)))
    end = float(event.reach_end_km) - float(cfg.reach_end_shift_km)
    spec(f"reach end {end:.1f} km (5 km upstream)", "reach", "14.5",
         frame(base, base["in_set"] & (base.river_km <= end)))
    spec("detached reach", "outcome", "14.6: detached_reach_m as the distance outcome",
         frame(base, outcomes={"y_height": "reach_height_m_fab", "y_dist": "detached_reach_m_fab",
                               "y_over": "overtopped_fab", "censored": "censored_fab"}))
    spec("winsorised outcomes", "outcome", "14.19 f: 97.5th percentile caps",
         frame(base, outcomes={"y_height": "reach_height_m_wins_fab",
                               "y_dist": "reach_distance_m_wins_fab",
                               "y_over": "overtopped_fab", "censored": "censored_fab"}))
    spec("Model A untrimmed", "overlap", "14.8: no common-support trimming", frame(base),
         trim=False)
    keep = primary["_cs"]["keep"] if "_cs" in primary else None
    if keep is not None:
        pdf = frame(base)
        spec("Model B: both banks survived trimming", "overlap", "14.8", pdf,
             b_df=pdf[keep.reindex(pdf.index).fillna(False)])
    spec("first section only", "river section", "14.16: above the first section break",
         frame(base, base["in_set"] & (base.river_section == 0)), with_sections=False)
    spec("channel_resolved = 1", "channel", "14.18",
         frame(base, base["in_set"] & (base.channel_resolved == 1)))
    t = analysis_rows(T, L, cfg, ignore_reasons=frozenset({"hillslope_failure", "tributary_fan"}))
    spec("manual exclusions included", "exclusions", "14.19 e", frame(t))
    for reason, ref in (("river_not_visible", "14.20"), ("near_bank_nodata", "14.24")):
        z = include_as_zero(T, reason)
        zt = analysis_rows(z, L, cfg).assign(in_set=z["analysis_set"] == 1)
        for k in ("reach_distance_m", "reach_height_m", "overtopped", "censored"):
            zt[f"{k}_fab"] = z[f"{k}_fab"]
        spec(f"{reason} rows as reach 0", "zeros", f"{ref}: {z.attrs['added_rows']} rows added",
             frame(zt))
    umin = float(cfg.usable_length_sensitivity_m)
    spec(f"usable length >= {umin:g} m", "neighbours", "14.21",
         frame(base, base["in_set"] & (base.usable_length_m >= umin)))
    spec("wall slope: complete cases", "wall slope", "user decision: rows with a valley wall only",
         frame(base, base["in_set"] & (base.valley_wall_found_fab == 1)))

    pk = placebo_key(event_key)
    placebo_ok = (out_root / pk / "transects.csv").exists()
    if placebo_ok:
        P, PL = load(pk)
        vis = P.drop_duplicates("station_id").set_index("station_id")["river_visible_post"]
        seen = base.station_id.map(vis) == 1
        spec("stations with a visible river in the placebo", "placebo comparability",
             "(b): real event restricted to stations whose river was visible in the placebo",
             frame(base, base["in_set"] & seen))

    for name, (desc, _) in rebuild_variants(cfg).items():
        key = variant_key(event_key, name)
        if rebuild or not (out_root / key / "transects.csv").exists():
            build_variants(event_key, config_path, force=rebuild, only=[name])
        VT, VL = load(key)
        vt = analysis_rows(VT, VL, cfg).assign(in_set=VT["analysis_set"] == 1)
        spec(desc, "rebuilt", f"rebuilt from label to table ({key})", frame(vt))

    runs = {"primary": primary}
    notes = {"primary": ("primary", "combined label, FABDEM, 100 m belt, trimmed A")}
    for name, group, note, df, kw in specs:
        print(f"Sensitivity: {name} ({len(df)} rows) ...")
        kw.setdefault("with_sections", sections)
        runs[name] = fit_run(df, cfg, **kw)
        notes[name] = (group, note)
        results += flat_rows(name, group, note, runs[name])
    results.append({"run": "reach end +5 km downstream", "group": "reach", "model": "A",
                    "outcome": "all", "note": "14.5",
                    "error": "not run: the verified centerline ends at reach_end_km "
                             f"({event.reach_end_km} km)"})

    # Placebo through the same code (a).
    placebo_runs = {}
    if placebo_ok:
        pbase = analysis_rows(P, PL, cfg)
        if not (pbase["in_set"] == (P["analysis_set"] == 1)).all():
            raise RuntimeError("Placebo analysis set could not be rebuilt.")
        z = river_not_visible_as_zero(P)
        zt = pbase.assign(in_set=z["analysis_set"] == 1)
        for k in ("reach_distance_m", "reach_height_m", "overtopped", "censored"):
            zt[f"{k}_fab"] = z[f"{k}_fab"]
        lab = placebo_labels(event)
        real_col = f"forest_{lab['real_year']}"
        p20 = pbase.merge(T[["station_id", "side", forest]].rename(columns={forest: real_col}),
                          on=["station_id", "side"], how="left")
        p20.index = pbase.index
        variants = {
            "placebo: kept rows (primary)": (frame(pbase),
                                             f"{lab['placebo_year']} forest, analysis set as built"),
            "placebo: + 14.20 assumed-zero rows": (frame(zt), NULL_PUSH),
            f"placebo: {lab['real_year']} forest": (
                frame(p20, p20["in_set"] & p20[real_col].notna(), forest_col=real_col),
                f"real event's {lab['real_season']} forest (sensitivity)"),
        }
        for name, (df, note) in variants.items():
            print(f"Placebo: {name} ({len(df)} rows) ...")
            placebo_runs[name] = fit_run(df, cfg, sections, percentile=name.endswith("(primary)"))
            results += flat_rows(name, "placebo", note, placebo_runs[name])

    table = pd.DataFrame(results)
    written = [_save(table, out_dir / "stats_results.csv")]
    flips = S.sign_flips(_by_key(primary), {n: _by_key(r) for n, r in runs.items() if n != "primary"},
                         FLIP_KEYS)
    written.append(_save(pd.DataFrame(flips, columns=["run", "model", "outcome", "primary_beta",
                                                      "run_beta", "run_p0"]),
                         out_dir / "stats_sign_flips.csv"))
    checks, verdict, reasons = _verdict(primary, runs, placebo_runs, cfg, flips, event_key)
    placebo_drop = _placebo_dropping(P, forest, forest_terciles) if placebo_ok else None
    written += _figures(out_dir, primary, table, cfg)
    text = _report(event_key, event, cfg, primary, runs, notes, table, placebo_runs,
                   placebo_drop, checks, verdict, reasons, flips,
                   standin_note=_worldcover_standin_note(event, cfg, T, base))
    report = out_dir / "STATS_REPORT.md"
    report.write_text(text, encoding="utf-8")
    written.append(report)
    summary = {"verdict": verdict, "reasons": reasons, "checks": checks,
               "primary": _plain(primary), "placebo": {k: _plain(v) for k, v in placebo_runs.items()}}
    path = out_dir / "stats_summary.json"
    path.write_text(json.dumps(summary, indent=2, default=_json), encoding="utf-8")
    written.append(path)
    for p in written:
        print(p.relative_to(PROJECT_ROOT))
    return written


# Outcomes compared for 8.10 sign flips: both models' OLS outcomes and Model A's logistic.
FLIP_KEYS = (("A", "y_height"), ("A", "y_logdist"), ("A", "y_over"),
             ("B", "y_height"), ("B", "y_logdist"))


def _by_key(run: dict) -> dict:
    return {(m, y): r for m in ("A", "B") for y, r in run.get(m, {}).items()}


def _save(frame: pd.DataFrame, path: Path) -> Path:
    frame.to_csv(path, index=False)
    return path


def _plain(run: dict) -> dict:
    return {k: v for k, v in run.items() if not k.startswith("_")}


def _json(v):
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, (pd.Series, pd.DataFrame)):
        return None
    return str(v)


def _neg_sig(r: dict, alpha: float) -> bool:
    return "error" not in r and r["beta"] < 0 and r["p0"] <= alpha


def _neg(r: dict) -> bool:
    return "error" not in r and r["beta"] < 0


def _includes_null(r: dict, alpha: float) -> bool:
    if "error" in r:
        return False
    if "or_interval" in r:
        return r["or_interval"][0] <= 1 <= r["or_interval"][1]
    return r["p0"] > alpha


def _validation_inputs(cfg, event_key: str) -> tuple:
    """(validation_ok, validation_note, validated) from this event's own validation files.

    An event with neither file has not been validated yet: all three are NOT_RUN, never
    another event's numbers and never a failure."""
    folder = validation_results_dir(cfg, event_key)
    summary, per_point = folder / "validation_summary.csv", folder / "validation_per_point.csv"
    if not summary.exists() and not per_point.exists():
        return S.NOT_RUN, S.NOT_RUN, S.NOT_RUN
    return (*_validation_flag(cfg, summary), per_point.exists())


def _validation_flag(cfg, path) -> tuple[bool, str]:
    if not path.exists():
        return False, "no validation summary"
    v = pd.read_csv(path)
    g = v[(v.subset == "main") & (v.breakdown == "forest_group")].set_index("category")
    low, high = g.loc["Low"], g.loc["High"]
    n_min = int(cfg.validation_min_group_n)
    if low.n < n_min or high.n < n_min:
        return False, f"a forest group has fewer than {n_min} validation points"
    gap = float(high.share_under_detected - low.share_under_detected)
    return gap <= cfg.validation_recall_gap, (
        f"under-detected High {high.share_under_detected:.1%} vs Low "
        f"{low.share_under_detected:.1%} (gap {gap * 100:.1f} points, flag above "
        f"{cfg.validation_recall_gap * 100:.0f})")


def pending_lines(checks: dict) -> list[str]:
    """A bold note naming the verdict inputs this event has not run yet (none: no lines)."""
    names = [name for name, keys in S.PENDING_INPUTS
             if any(S.not_run(checks.get(k)) for k in keys)]
    if not names:
        return []
    one = len(names) == 1
    return [f"**Not yet run for this event: {' and '.join(names)}. The verdict is pending "
            f"until {'it is' if one else 'they are'} run; {'its' if one else 'their'} checks "
            "below read 'not yet run'.**", ""]


def _verdict(primary, runs, placebo_runs, cfg, flips, event_key):
    A, B = primary["A"], primary["B"]
    alpha = float(cfg.interval_alpha)
    val_ok, val_note, validated = _validation_inputs(cfg, event_key)
    narrow = []
    for model in (A, B):
        for y in ("y_height", "y_logdist"):
            r = model[y]
            narrow.append("error" not in r and S.is_reasonably_narrow(
                r["interval"], primary["sd"][y], float(cfg.narrow_effect_sd)))
    over = A["y_over"]
    narrow.append("error" not in over and S.or_narrow(over["or_interval"], cfg.narrow_or_low,
                                                       cfg.narrow_or_high))
    pl = placebo_runs.get("placebo: kept rows (primary)")
    placebo_effects = []
    if pl:
        for model in ("A", "B"):
            for y in ("y_height", "y_logdist", "y_over"):
                r = pl[model].get(y)
                if r is None or "error" in r:      # not estimated: no evidence either way
                    continue
                if not _includes_null(r, alpha):
                    placebo_effects.append(f"Model {model} {OUT_NAMES[y]}: {_ci(r)}")
    # No placebo run yet is not a failed placebo (8.9): the check waits for it.
    placebo_null = S.NOT_RUN if not placebo_runs else bool(pl) and not placebo_effects
    placebo_note = S.NOT_RUN if not placebo_runs else "; ".join(placebo_effects) or "none"
    checks = {
        "a_height_neg_sig": _neg_sig(A["y_height"], alpha),
        "b_height_neg_sig": _neg_sig(B["y_height"], alpha),
        "dist_same_sign": _neg(A["y_logdist"]) and _neg(B["y_logdist"]),
        "labels_same_sign": all(_neg(runs[f"label {l}"]["A"]["y_height"])
                                for l in ("optical_only", "sar_only")),
        "dems_same_sign": _neg(runs["DEM Copernicus GLO-30"]["A"]["y_height"]),
        "placebo_null": placebo_null,
        "placebo_effects": placebo_note,
        "validation_ok": val_ok, "validation_note": val_note,
        "validated": validated,
        "enough_blocks": S.enough_blocks(primary["a_blocks"], int(cfg.min_blocks)),
        "overlap_ok": not primary.get("support", {}).get("too_thin", True),
        "model_b_informative": primary["b_info"]["informative"],
        "no_sign_flips": not flips,
        "sign_flips": f"{len(flips)} (run x model x outcome); see the sign-flip table",
        "all_include_null": all(_includes_null(r, alpha) for r in (A["y_height"], A["y_logdist"],
                                                                     B["y_height"], B["y_logdist"], over)),
        "all_narrow": all(narrow),
    }
    verdict, reasons = S.decide_verdict(checks)
    return checks, verdict, reasons


def placebo_labels(event) -> dict:
    """Year and season labels for the placebo runs, from the event's own windows:
    the placebo forest is from the placebo S2 pre window, the sensitivity forest from
    the real event's S2 pre window (14.2)."""
    import calendar

    from .dates import composite_year
    real, plac = event.s2_pre, event.placebo.s2_pre
    real_year = composite_year(real.start, real.end)
    months = "-".join(calendar.month_abbr[int(d[5:7])] for d in (real.start, real.end))
    return {"placebo_year": composite_year(plac.start, plac.end), "real_year": real_year,
            "real_season": f"{months} {real_year}"}


def _placebo_dropping(P, forest, forest_terciles) -> dict:
    t = P[P.in_reach == 1].copy()
    t["tercile"] = forest_terciles(t[forest])
    nv = t.exclude_reasons.fillna("").str.contains("river_not_visible")
    share = nv.groupby(t.tercile).mean()
    count = nv.groupby(t.tercile).sum()
    total = t.groupby("tercile").size()
    return {g: (float(share[g]), int(count[g]), int(total[g])) for g in ("Low", "Mid", "High")}


# ---------------------------------------------------------------------------
# Figures and report
# ---------------------------------------------------------------------------

def _figures(out_dir: Path, primary: dict, table: pd.DataFrame, cfg) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    written = []
    if "_cs" in primary:
        cs = primary["_cs"]
        fig, ax = plt.subplots(figsize=(6, 3.6), dpi=120)
        bins = np.linspace(0, 1, 26)
        ax.hist(cs["p"][cs["high"]], bins=bins, alpha=0.6, label="forest_high", color="#2f855a")
        ax.hist(cs["p"][~cs["high"]], bins=bins, alpha=0.6, label="not forest_high", color="#b7791f")
        for b in cs["bounds"]:
            ax.axvline(b, color="black", lw=0.8, ls="--")
        ax.set_xlabel("propensity score p")
        ax.set_ylabel("rows")
        ax.set_title(f"Common support (8.1): kept {cs['kept']} of {cs['rows']} rows", fontsize=9)
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        written.append(out_dir / "fig_propensity.png")
        fig.savefig(written[-1])
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 3.6), dpi=120)
    ax.hist(primary["_pairs"]["d_forest"].abs(), bins=20, color="#4a5568")
    ax.axvline(float(cfg.model_b_min_forest_diff), color="red", lw=1)
    ax.set_xlabel("|F_L - F_R|")
    ax.set_ylabel("stations")
    ax.set_title("Model B: forest difference between banks (14.8)", fontsize=9)
    fig.tight_layout()
    written.append(out_dir / "fig_model_b_forest_difference.png")
    fig.savefig(written[-1])
    plt.close(fig)

    for model in ("A", "B"):
        part = table[(table.model == model) & (table.outcome == "y_height") & table["beta"].notna()]
        part = part.iloc[::-1]
        fig, ax = plt.subplots(figsize=(7.5, 0.28 * len(part) + 1.2), dpi=120)
        y = np.arange(len(part))
        colours = ["#c53030" if g == "primary" else "#805ad5" if g == "placebo" else "#2b6cb0"
                   for g in part.group]
        ax.errorbar(part.beta, y, xerr=[part.beta - part.lo, part.hi - part.beta], fmt="o",
                    ecolor=colours, color="none", lw=1)
        ax.scatter(part.beta, y, c=colours, s=14, zorder=3)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(part.run, fontsize=7)
        ax.set_xlabel("beta1, reach height (m) for forest 0 -> 100 % (14.15 interval)")
        ax.set_title(f"Model {model}: reach height across runs", fontsize=9)
        fig.tight_layout()
        written.append(out_dir / f"fig_forest_plot_{model}_height.png")
        fig.savefig(written[-1])
        plt.close(fig)
    return written


def _f(v, d=2):
    return "NA" if v is None or (isinstance(v, float) and not math.isfinite(v)) else f"{v:.{d}f}"


def _ci(r, d=2):
    if "error" in r:
        msg = r["error"]
        return msg if msg.startswith("not estimable") else f"not estimable ({msg[:120]})"
    if "or_interval" in r:
        return (f"OR {_f(r['odds_ratio'])} [{_f(r['or_interval'][0])}, {_f(r['or_interval'][1])}] "
                f"(approximate, t with {r['df']} df)")
    return f"{_f(r['beta'], d)} [{_f(r['interval'][0], d)}, {_f(r['interval'][1], d)}], p(0) = {_f(r['p0'], 3)}"


OUT_NAMES = {"y_height": "reach height (m)", "y_logdist": "log1p reach distance",
             "y_over": "overtopped (logistic)"}


def _sign_disagreements(primary: dict) -> list[str]:
    """8.3: outcomes where Model A and Model B point in opposite directions."""
    out = []
    for y in ("y_height", "y_logdist"):
        a, b = primary["A"][y], primary["B"][y]
        if "beta" in a and "beta" in b and np.sign(a["beta"]) != np.sign(b["beta"]):
            out.append(f"{OUT_NAMES[y]}: Model A {_f(a['beta'])}, Model B {_f(b['beta'])}")
    return out


def _worldcover_standin_note(event, cfg, T, base) -> str | None:
    """Report line for an event with a WorldCover stand-in decision (14.23 a); None otherwise."""
    standin = getattr(event, "worldcover_standin", None)
    if standin is None:
        return None
    share = pd.to_numeric(T[f"forest_source_worldcover_frac_{int(cfg.belt_width_m)}"],
                          errors="coerce").fillna(0)
    used = int((base["in_set"] & (share > 0)).sum())
    return (f"- **WorldCover stand-in (14.23 a, user decision {standin.decided_on}).** WorldCover "
            f"{standin.year} stands in for {standin.for_year}, which has no map, where Dynamic "
            f"World has no observation in the S2 pre window. {used} of {int(base['in_set'].sum())} "
            "analysis-set rows use it (forest_source worldcover or mixed); the "
            "'forest_source dynamic_world only' run leaves them out.")


def event_limitations(event, event_key: str) -> list[str]:
    """A Limitations section from the event's optional `limitations` list in config.

    Empty when the event has none, so reports of events without the list are unchanged.
    """
    items = getattr(event, "limitations", None) or []
    lines = []
    for text in items:
        if not isinstance(text, str) or not text.strip():
            raise RuntimeError(f"{event_key} limitations has an empty or non-text entry.")
        lines.append(f"- {text.strip()}")
    return ["", "## Limitations", "", *lines] if lines else []


def model_a_formula_line(primary: dict) -> str:
    """The Model A formula as fitted (after common-support trimming), naming any term left
    out because it is constant in the trimmed rows. Unchanged text when nothing is left out."""
    fitted = primary.get("_rhs_a_fitted", primary["rhs_a"])
    line = (f"Model A (8.2, 14.16, 14.18): `y ~ {fitted}`. OLS for reach height and "
            "log1p(reach distance); logistic for overtopped.")
    dropped = [t for t in primary["rhs_a"].split(" + ") if t not in fitted.split(" + ")]
    if dropped:
        names = ", ".join(f"`{t}`" for t in dropped)
        line += (f" {names} {'is' if len(dropped) == 1 else 'are'} in the formula before "
                 "trimming but left out of the fit because "
                 f"{'it is' if len(dropped) == 1 else 'they are'} constant in the trimmed "
                 "Model A rows.")
    return line


def _report(event_key, event, cfg, primary, runs, notes, table, placebo_runs, placebo_drop,
            checks, verdict, reasons, flips, standin_note=None) -> str:
    A, B = primary["A"], primary["B"]
    sup = primary.get("support", {})
    disagree = _sign_disagreements(primary)
    rivers = list(getattr(event, "rivers", []))
    breaks = list(event.section_breaks_km or [])
    lines = [
        f"# Statistics: {event_key}",
        "",
        f"Reproduce with `python -m floodforest.cli stats --event {event_key}`. Seed {cfg.seed}.",
        "",
        f"## Verdict (8.10): {verdict}",
        "",
        "Wording follows the spec: forest is *associated with* reach; nothing here shows that "
        "forest protects or causes anything.",
        "",
        *([f"**The placebo shows a forest association ({checks['placebo_effects']}). Under 8.9 "
           "the method is producing artefacts and the main result cannot be trusted.**", ""]
          if not S.not_run(checks["placebo_null"]) and not checks["placebo_null"] else []),
        *pending_lines(checks),
        *([f"**Model A and Model B disagree in sign (8.3): {'; '.join(disagree)}.**", ""]
          if disagree else []),
        *([f"**{len(flips)} sensitivity estimates flip sign relative to the primary run "
           "(8.10); see the sign-flip table.**", ""] if flips else []),
        *(f"- {r}" for r in reasons),
        "",
        "| Check | Result |", "|---|---|",
        *(f"| {k} | {v} |" for k, v in checks.items()),
        "",
        *report_lines(flood_records(event, event_key)),
        "## Models",
        "",
        model_a_formula_line(primary),
        "",
        f"Model B (8.3): `(y_L - y_R) ~ {b_rhs_for(primary['_pairs'])}`, stations with both banks "
        "in the analysis set; not trimmed (14.8).",
        "",
        "Wall slope is 0 where no valley wall was found, with valley_wall_found in both models "
        "(user decision); a complete-case run is in the sensitivity table. Model B pairs are "
        "built from Model B's own variables, so a bank missing only a Model A control (for "
        "example curvature at the end of the line) still pairs.",
        "",
        "### Methods notes",
        "",
        f"- **Gap tolerance (7.2).** Samples are {cfg.grid_m:g} m apart and a run of N samples "
        f"spans (N - 1) x {cfg.grid_m:g} m. A non-impact gap of up to "
        f"{int(round(cfg.impact_gap_m / cfg.grid_m)) + 1} consecutive samples "
        f"(span <= {cfg.impact_gap_m:g} m) stays inside the run; a nodata gap of up to "
        f"{int(round(cfg.nodata_gap_m / cfg.grid_m)) + 1} samples (span <= {cfg.nodata_gap_m:g} m) "
        "does not end it, and a longer one ends it with censored = 1. A nodata sample neither "
        "resets nor adds to the non-impact count. The same convention sets the valley-wall run "
        f"({cfg.wall_run_m:g} m = {int(round(cfg.wall_run_m / cfg.grid_m)) + 1} samples) and the "
        f"channel-run merge ({cfg.channel_merge_gap_m:g} m, 14.18).",
        f"- **Channel width (14.18).** channel_width_m is the length of the chosen channel run "
        f"counted in whole pixels: last sample minus first sample plus one pixel ({cfg.grid_m:g} m). "
        f"Where no channel run is found it is {cfg.grid_m:g} m (one pixel, an upper bound).",
        "- **Censored model (14.7).** The Weibull AFT uses lifelines' robust (sandwich) standard "
        "errors, not clustered on blocks; it is a sensitivity conditional on reach > 0.",
        "- **Percentile intervals (8.4, 14.15).** The 1 km pairs percentile interval is computed "
        "for the primary run and the placebo's primary variant only; sensitivity runs report the "
        "14.15 interval alone.",
        "- **E-values (8.7, 14.9).** Each model uses its own SD(y): Model A the SD of y in its "
        "rows, Model B the SD of y_L - y_R in its pairs. The 14.9 prevalence of overtopped is "
        "taken over the run's full analysis set. 14.10 narrowness uses Model A's SD(y) for both "
        "models, as before.",
        *([standin_note] if standin_note else []),
        "",
        "Primary interval (14.15): wild cluster bootstrap-t, null imposed, CR1, Rademacher, "
        f"{cfg.wild_bootstrap_draws} draws, {cfg.bootstrap_block_km:g} km blocks; p(0) reported; "
        f"interval = b0 with p(b0) > {cfg.interval_alpha:g}. Logistic: cluster-robust SE with a t critical value "
        "(blocks - 1 df), approximate. Beside it: the 1 km pairs block bootstrap percentile "
        f"interval ({cfg.bootstrap_reps} reps), {S.PERCENTILE_LABEL}.",
        "",
        "## Common support (8.1, 14.8)",
        "",
        (f"Split: {sup.get('rule')}. Kept {sup.get('kept')} of {sup.get('rows')} rows "
         f"(dropped {sup.get('dropped')}), {sup.get('kept_segments')} segments; p bounds "
         f"{tuple(round(b, 3) for b in sup.get('bounds', ()))}. Too thin to interpret: "
         f"{sup.get('too_thin')}." if "error" not in sup else sup["error"]),
        f"Rows without every Model A value (curvature at the line end): "
        f"{primary['dropped_missing']}; these rows can still enter Model B.",
        "",
        "## Primary results",
        "",
        f"Model A: {primary['a_rows']} rows, {primary['a_stations']} stations, "
        f"{primary['a_segments']} segments, {primary['a_blocks']} blocks. "
        f"Model B: {primary['b_info']['stations']} stations, {primary['b_segments']} segments, "
        f"{primary['b_blocks']} blocks; "
        f"{primary['b_info']['stations_large_diff']} with |F_L - F_R| >= "
        f"{cfg.model_b_min_forest_diff:g} (minimum {cfg.model_b_min_stations}); "
        f"informative: {primary['b_info']['informative']}.",
        "",
        "| Model | Outcome | beta1 for F 0 -> 1 [95% interval], p(0) | Percentile interval (too narrow with few segments) | Per IQR of F | E-value (point / CI) | Reasonably narrow (14.10) |",
        "|---|---|---|---|---|---|---|",
    ]
    for model, part in (("A", A), ("B", B)):
        for y in ("y_height", "y_logdist"):
            r = part[y]
            if "error" in r:
                lines.append(f"| {model} | {OUT_NAMES[y]} | {_ci(r)} | | | | |")
                continue
            pc = r.get("percentile", {}).get("interval", (np.nan, np.nan))
            iqr = primary["iqr_F"]
            e = S.evalue_for(r, primary["sd"][y] if model == "A" else primary["sd_b"][y])
            nar = S.is_reasonably_narrow(r["interval"], primary["sd"][y], float(cfg.narrow_effect_sd))
            lines.append(
                f"| {model} | {OUT_NAMES[y]} | {_ci(r)} | [{_f(pc[0])}, {_f(pc[1])}] | "
                f"{_f(r['beta'] * iqr)} [{_f(r['interval'][0] * iqr)}, {_f(r['interval'][1] * iqr)}] "
                f"(IQR {_f(iqr)}) | {_f(e['evalue'])} / {_f(e['evalue_ci'])} | {nar} |")
    over = A["y_over"]
    if "error" in over:
        lines.append(f"| A | overtopped | {_ci(over)} | | | | |")
    else:
        e = S.evalue_logistic(over, primary["prevalence"], float(cfg.overtop_prevalence_for_rr))
        lines.append(f"| A | overtopped | {_ci(over)} | | | {_f(e['evalue'])} / {_f(e['evalue_ci'])} "
                     f"(prevalence {primary['prevalence']:.0%}) | "
                     f"{S.or_narrow(over['or_interval'], cfg.narrow_or_low, cfg.narrow_or_high)} |")
    lines += ["", "log1p distance: exp(beta1) - 1 is the proportional change in (1 + distance).",
              "",
              "E-value (8.7): an E-value of 2 means an unmeasured factor would need to be "
              "associated with both forest and flood reach by a risk ratio of 2 each to erase "
              "the result. The CI E-value is 1 when the interval includes the null.",
              "", "### Censored model, conditional on reach > 0 (14.7, sensitivity)", ""]
    for y, r in primary.get("aft", {}).items():
        if "error" in r:
            lines.append(f"- {y}: not estimable ({r['error'][:100]})")
        else:
            lines.append(f"- {y}: time ratio for F 0 -> 1 {_f(r['time_ratio'])} "
                         f"[{_f(r['ci_time_ratio'][0])}, {_f(r['ci_time_ratio'][1])}], p {_f(r['p'], 3)}; "
                         f"{r['rows']} rows, {r['censored']} censored; "
                         f"{r['zero_outcome_rows_left_out']} rows with a zero value left out. {r['method']}.")
    lines += ["", "## Sensitivity runs (8.8 and section 14)", "",
              "Model A reach height and Model B reach height, 14.15 interval. All runs, all outcomes: "
              "`stats_results.csv`. Forest plots: `fig_forest_plot_A_height.png`, "
              "`fig_forest_plot_B_height.png`.", "",
              "| Run | Group | Note | A rows | A height beta1 [interval], p(0) | A log1p distance beta1 [interval], p(0) | B stations / segments | B height beta1 [interval], p(0) |",
              "|---|---|---|---:|---|---|---:|---|"]
    for name, run in runs.items():
        group, note = notes[name]
        lines.append(f"| {name} | {group} | {note} | {run['a_rows']} | {_ci(run['A']['y_height'])} | "
                     f"{_ci(run['A']['y_logdist'])} | {run['b_info']['stations']} / "
                     f"{run['b_segments']} | {_ci(run['B']['y_height'])} |")
    for r in table[table["error"].notna() & (table["outcome"] == "all")].itertuples():
        lines.append(f"| {r.run} | {r.group} | {r.note} | | {r.error} | | | |")
    if breaks:
        lines += ["", f"Section breaks at km {', '.join(f'{b:g}' for b in breaks)} "
                      f"({' then '.join(rivers)}). The first-section-only run (14.16 rule 5) uses "
                      f"the {rivers[0] if rivers else 'first river'} above the first break. Other "
                      "rivers joining below a confluence carried their own flood water in the same "
                      "storm."]
    lines += ["", "### Sign flips across sensitivity runs (8.10)", "",
              "A run flips when its beta1 has the opposite sign to the primary run's for the same "
              "model and outcome (Model A and B reach height and log1p distance, Model A "
              "overtopped). Any flip makes the verdict inconclusive. All runs are compared, "
              "including the literal steep filter. Full list: `stats_sign_flips.csv`.", ""]
    if flips:
        lines += ["| Run | Model | Outcome | Primary beta1 | Run beta1 | Run p(0) |",
                  "|---|---|---|---:|---:|---:|"]
        lines += [f"| {f['run']} | {f['model']} | {OUT_NAMES[f['outcome']]} | "
                  f"{_f(f['primary_beta'])} | {_f(f['run_beta'])} | "
                  f"{_f(f['run_p0'], 3) if f['run_p0'] is not None else 'NA'} |" for f in flips]
    else:
        lines.append("No run flips sign.")
    nv = runs.get("river_not_visible rows as reach 0")
    lines += ["", f"(d) The real event's 14.20 sensitivity (river_not_visible rows included as reach 0) "
                  f"is reported: {'yes, row above' if nv else 'NO'}.", ""]
    if placebo_runs:
        lines += ["## Placebo through the same models (8.9)", "",
                  "Same Model A and Model B code and the same 14.15 interval as the real event.", "",
                  "| Variant | Note | A rows | A reach height | A log1p distance | A overtopped | B stations / segments | B reach height | B log1p distance |",
                  "|---|---|---:|---|---|---|---:|---|---|"]
        for name, run in placebo_runs.items():
            note = NULL_PUSH if "14.20" in name else ""
            lines.append(f"| {name} | {note} | {run['a_rows']} | {_ci(run['A']['y_height'])} | "
                         f"{_ci(run['A']['y_logdist'])} | {_ci(run['A']['y_over'])} | "
                         f"{run['b_info']['stations']} / {run['b_segments']} | {_ci(run['B']['y_height'])} | "
                         f"{_ci(run['B']['y_logdist'])} |")
        if placebo_drop:
            hi, lo = placebo_drop["High"], placebo_drop["Low"]
            lines += ["",
                      f"(c) The placebo drops {hi[0]:.0%} of High-forest rows for river_not_visible "
                      f"({hi[1]} of {hi[2]}), against {lo[0]:.0%} of Low-forest rows ({lo[1]} of {lo[2]}) "
                      f"(terciles of the {placebo_labels(event)['placebo_year']} forest over all "
                      "in-reach rows). The placebo therefore "
                      "cannot speak for the densest forest.",
                      "",
                      "The 14.20 assumed-zero variant adds those rows with reach 0; it "
                      f"{NULL_PUSH}, so it is not evidence that the placebo is null."]
    lines += event_limitations(event, event_key)
    lines += ["", "## Files", "", "- `stats_results.csv`, `stats_summary.json`, `stats_sign_flips.csv`",
              "- `fig_propensity.png`, `fig_model_b_forest_difference.png`, "
              "`fig_forest_plot_A_height.png`, `fig_forest_plot_B_height.png`", ""]
    return "\n".join(lines) + "\n"
