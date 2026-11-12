"""Descriptive collinearity diagnostic for Model A and Model B (no model is changed).

Rebuilds each event's primary Model A rows (analysis set, rows with every Model A value,
common-support trimmed) and Model B pairs exactly as stats_run.fit_run does, then writes
pairwise correlations and VIFs to outputs/diagnostics/. Fits nothing.
"""

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import patsy

from floodforest import stats as S
from floodforest.config import PROJECT_ROOT, load_config
from floodforest.stats_run import analysis_rows, b_rhs_for, rhs_for

R_LIMIT, VIF_LIMIT = 0.7, 5.0
cfg = load_config(PROJECT_ROOT / "config.yaml")
out_dir = PROJECT_ROOT / cfg.paths.outputs / "diagnostics"
out_dir.mkdir(parents=True, exist_ok=True)

# Model A predictors as they enter the model (8.2, 14.18, user decision on wall slope).
A_SCALAR = {
    "forest_frac": "forest_frac",
    "log_channel_width": "np.log(channel_width_m)",
    "gradient": "gradient",
    "log_valley_width": "np.log(valley_width_m)",
    "wall_slope_deg": "wall_slope_deg",
    "bank_slope_deg": "bank_slope_deg",
    "log10_upa": "log10_upa",
    "outer_bend": "outer_bend",
    "curvature": "curvature",
    "bank_elev_m": "bank_elev_m",
    "valley_wall_found": "valley_wall_found",
    "channel_resolved": "channel_resolved",
    "side_right": "C(side)[T.right]",
    "river_km": None,  # enters through the spline; its raw correlations are shown too
}


def term_vifs(design: pd.DataFrame, slices: dict) -> pd.DataFrame:
    """VIF per term from the design matrix (intercept dropped); GVIF for multi-column terms
    (Fox & Monette 1992), with GVIF^(1/(2*df)) to compare with sqrt of the VIF limit."""
    cols = [c for c in design.columns if c != "Intercept"]
    corr = np.corrcoef(design[cols].to_numpy(float), rowvar=False)
    det_all = np.linalg.det(corr)
    rows = []
    for term, sl in slices.items():
        if term == "Intercept":
            continue
        idx = [cols.index(c) for c in design.columns[sl] if c != "Intercept"]
        rest = [i for i in range(len(cols)) if i not in idx]
        gvif = (np.linalg.det(corr[np.ix_(idx, idx)]) * np.linalg.det(corr[np.ix_(rest, rest)])
                / det_all)
        df = len(idx)
        rows.append({"term": term, "df": df, "vif_or_gvif": float(gvif),
                     "gvif_1_2df": float(gvif ** (1 / (2 * df)))})
    return pd.DataFrame(rows)


def high_pairs(corr: pd.DataFrame, limit: float) -> list[dict]:
    out = []
    names = list(corr.columns)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            r = corr.loc[a, b]
            if np.isfinite(r) and abs(r) > limit:
                out.append({"a": a, "b": b, "r": round(float(r), 3)})
    return out


summary = {}
for event_key in ("lende2025", "melamchi2021"):
    event = getattr(cfg.events, event_key)
    root = PROJECT_ROOT / cfg.paths.outputs
    T = pd.read_csv(root / event_key / "transects.csv")
    L = pd.read_csv(root / event_key / "transects_labels.csv")
    base = analysis_rows(T, L, cfg)
    forest = f"forest_frac_{int(cfg.belt_width_m)}"
    df = S.model_frame(base[base["in_set"]], forest, "fab", None)
    sections = bool(event.section_breaks_km)

    # Model A rows: as fit_run (a_ok, then common-support trimming).
    a_all = df[df["a_ok"]]
    rhs_all = rhs_for(a_all, cfg, sections)
    cs = S.common_support(a_all, rhs_all.split(" + ", 1)[1], float(cfg.forest_high_threshold),
                          float(cfg.forest_high_min_share), float(cfg.overlap_min_row_fraction),
                          int(cfg.overlap_min_segments))
    a_df = a_all[cs["keep"]]
    rhs = rhs_for(a_df, cfg, sections)
    design = patsy.dmatrix(rhs, a_df, return_type="dataframe")

    scalar = pd.DataFrame(index=a_df.index)
    for name, col in A_SCALAR.items():
        if col is None:
            scalar[name] = a_df["river_km"].to_numpy(float)
        elif col in design.columns:
            scalar[name] = design[col].to_numpy(float)
    pearson = scalar.corr(method="pearson")
    spearman = scalar.corr(method="spearman")
    vif_full = term_vifs(design, design.design_info.term_name_slices)
    # VIF among the named predictors only (no spline, side or section): for comparison.
    named = [c for c in scalar.columns if c not in ("river_km", "side_right")]
    named_design = scalar[named].assign(Intercept=1.0)[["Intercept"] + named]
    vif_named = term_vifs(named_design, {c: slice(i, i + 1)
                                         for i, c in enumerate(named_design.columns)})

    # Model B pairs: as fit_run (b_ok rows, both banks, never trimmed).
    pairs = S.paired_frame(df[df["b_ok"]])
    b_terms = b_rhs_for(pairs).split(" + ")
    b_corr = pairs[b_terms].corr(method="pearson")
    b_spear = pairs[b_terms].corr(method="spearman")
    b_design = patsy.dmatrix(" + ".join(b_terms), pairs, return_type="dataframe")
    b_vif = term_vifs(b_design, b_design.design_info.term_name_slices)

    prefix = out_dir / event_key
    pearson.round(4).to_csv(f"{prefix}_model_a_pearson.csv")
    spearman.round(4).to_csv(f"{prefix}_model_a_spearman.csv")
    vif_full.round(4).to_csv(f"{prefix}_model_a_vif_fitted_design.csv", index=False)
    vif_named.round(4).to_csv(f"{prefix}_model_a_vif_named_predictors_only.csv", index=False)
    b_corr.round(4).to_csv(f"{prefix}_model_b_pearson.csv")
    b_spear.round(4).to_csv(f"{prefix}_model_b_spearman.csv")
    b_vif.round(4).to_csv(f"{prefix}_model_b_vif.csv", index=False)

    stats_json = json.loads((root / event_key / "stats" / "stats_summary.json")
                            .read_text(encoding="utf-8"))
    summary[event_key] = {
        "model_a_rows": int(len(a_df)), "model_a_rhs": rhs,
        "model_b_pairs": int(len(pairs)), "model_b_rhs": " + ".join(b_terms),
        "matches_stats_step": {
            "a_rows": stats_json["primary"]["a_rows"] == len(a_df),
            "rhs_a": stats_json["primary"]["rhs_a"] == rhs,
            "b_stations": stats_json["primary"]["b_info"]["stations"] == len(pairs)},
        "model_a_pearson_above_0.7": high_pairs(pearson, R_LIMIT),
        "model_a_spearman_above_0.7": high_pairs(spearman, R_LIMIT),
        "model_a_vif_above_5_fitted_design": vif_full[
            (vif_full.df == 1) & (vif_full.vif_or_gvif > VIF_LIMIT)].round(2).to_dict("records"),
        "model_a_multi_column_terms": vif_full[vif_full.df > 1].round(2).to_dict("records"),
        "model_a_vif_above_5_named_only": vif_named[
            vif_named.vif_or_gvif > VIF_LIMIT].round(2).to_dict("records"),
        "model_b_pearson_above_0.7": high_pairs(b_corr, R_LIMIT),
        "model_b_spearman_above_0.7": high_pairs(b_spear, R_LIMIT),
        "model_b_vif_above_5": b_vif[b_vif.vif_or_gvif > VIF_LIMIT].round(2).to_dict("records"),
    }

(out_dir / "collinearity_summary.json").write_text(json.dumps(summary, indent=2),
                                                   encoding="utf-8")
(out_dir / "README.md").write_text(
    "# Collinearity diagnostic (descriptive only)\n\n"
    "Pairwise correlations (Pearson and Spearman) and variance inflation factors for the "
    "Model A predictors on each event's primary Model A rows (analysis set, rows with every "
    "Model A value, common-support trimmed), and for the Model B between-bank differences on "
    "the Model B pairs. Predictors are taken as they enter the model (log channel width, log "
    "valley width). `_vif_fitted_design` uses the full fitted design (spline in river_km, "
    "side, river_section where used); multi-column terms get a generalised VIF and "
    "GVIF^(1/(2 df)). `_vif_named_predictors_only` leaves out the spline, side and section. "
    "Flags: |r| > 0.7, VIF > 5. No model, filter or result was changed; nothing here feeds "
    "the stats step. Script: `collinearity_diagnostic.py`.\n", encoding="utf-8")
shutil.copy(Path(__file__), out_dir / "collinearity_diagnostic.py")
print(json.dumps(summary, indent=1))
