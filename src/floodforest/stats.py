"""Models, uncertainty and E-value (spec section 8, clarifications 14.7 to
14.10, 14.15 and 14.16).

Primary interval (14.15): wild cluster bootstrap-t with the null imposed,
clustered on 3 km blocks along the river. Sensitivity interval (8.4): the
1 km pairs block bootstrap with a percentile interval, which is known to be
too narrow with few segments.

These functions take a prepared table, so the tests can exercise them on
synthetic rows without any satellite data.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import patsy
import statsmodels.formula.api as smf
from scipy import stats as scipy_stats

PERCENTILE_LABEL = "known to be too narrow with few segments"


def block_ids(river_km: np.ndarray, reach_start_km: float, block_km: float,
              final_block_min_km: float,
              reach_end_km: float | None = None) -> np.ndarray:
    """Cluster blocks along the river (14.15).

    block_id = floor((river_km - reach_start_km) / block_km). A final block
    shorter than final_block_min_km is merged into the previous one. Its
    length runs to reach_end_km, or to the last station if no end is given.
    """
    km = np.asarray(river_km, dtype=float)
    ids = np.floor((km - reach_start_km) / block_km).astype(int)
    end = float(km.max()) if reach_end_km is None else float(reach_end_km)
    last = int(ids.max())
    final_length = end - (reach_start_km + last * block_km)
    if last > 0 and final_length < final_block_min_km:
        ids[ids == last] = last - 1
    return ids


def winsorise_upper(values: np.ndarray, percentile: float) -> np.ndarray:
    """Cap values above the given percentile (14.19 f sensitivity). NaN stays NaN."""
    data = np.asarray(values, dtype=float)
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return data.copy()
    cap = float(np.percentile(finite, percentile))
    out = data.copy()
    out[np.isfinite(out) & (out > cap)] = cap
    return out


def enough_blocks(n_blocks: int, min_blocks: int) -> bool:
    """False means the event's verdict can be at most "inconclusive"."""
    return n_blocks >= min_blocks


def model_a_rhs(forest: str = "forest_frac", spline_df: int = 5,
                with_sections: bool = False, extra_terms: tuple = ()) -> str:
    """Right-hand side of Model A (8.2), with river_section when needed (14.16).

    river_section enters as a fixed effect only when the reach has at least
    one confirmed section break; Model B never includes it. extra_terms adds
    the indicators the run uses (channel_resolved, 14.18; valley_wall_found,
    user decision on missing wall slope).
    """
    rhs = (f"{forest} + np.log(channel_width_m) + gradient"
           " + np.log(valley_width_m) + wall_slope_deg + bank_slope_deg"
           " + log10_upa + outer_bend + curvature + bank_elev_m + C(side)"
           f" + cr(river_km, df={spline_df})")
    for term in extra_terms:
        rhs += f" + {term}"
    if with_sections:
        rhs += " + C(river_section)"
    return rhs


# Model B differences (8.3): valley width, wall slope, bank slope, outer bend, plus the
# valley-wall indicator that goes with the 0-filled wall slope (user decision).
MODEL_B_DIFFS = ("valley_width_m", "wall_slope_deg", "valley_wall_found", "bank_slope_deg",
                 "outer_bend")
MODEL_A_EXTRA = ("valley_wall_found", "channel_resolved")
MODEL_COLUMNS = ("forest_frac", "channel_width_m", "gradient", "valley_width_m",
                 "wall_slope_deg", "valley_wall_found", "bank_slope_deg", "log10_upa",
                 "outer_bend", "curvature", "bank_elev_m", "channel_resolved", "river_km",
                 "river_section", "side", "block_id", "segment_id", "station_id")
# Model B needs only its own variables (8.3): forest, the differenced controls, the
# outcomes it fits, and the identifiers used to pair and cluster.
MODEL_B_COLUMNS = ("forest_frac",) + MODEL_B_DIFFS + ("side", "block_id", "segment_id",
                                                      "station_id", "river_section")


def model_frame(table: pd.DataFrame, forest_col: str, dem: str = "fab",
                outcome_cols: dict | None = None) -> pd.DataFrame:
    """Model columns with standard names (8.2) for the rows given.

    Wall slope is 0 where no valley wall was found and valley_wall_found carries
    that (user decision). a_ok marks rows with every Model A value; b_ok marks rows
    with Model B's own variables (8.3), so a bank missing only a Model A control
    (for example curvature at the line end) can still pair in Model B. Rows usable
    by neither are dropped. attrs: dropped_missing (rows not usable by Model A) and
    dropped_missing_b (rows not usable by Model B).
    """
    outcome_cols = outcome_cols or {
        "y_height": f"reach_height_m_{dem}", "y_dist": f"reach_distance_m_{dem}",
        "y_over": f"overtopped_{dem}", "censored": f"censored_{dem}"}
    t = table
    wall_found = pd.to_numeric(t[f"valley_wall_found_{dem}"], errors="coerce")
    out = pd.DataFrame({
        "station_id": t["station_id"], "side": t["side"], "river_km": t["river_km"],
        "segment_id": t["segment_id"], "block_id": t["block_id"],
        "river_section": t["river_section"], "channel_resolved": t["channel_resolved"],
        "forest_frac": pd.to_numeric(t[forest_col], errors="coerce"),
        "channel_width_m": t["channel_width_m"], "gradient": t[f"gradient_{dem}"],
        "valley_width_m": t[f"valley_width_m_{dem}"],
        "wall_slope_deg": np.where(wall_found == 1, t[f"wall_slope_deg_{dem}"], 0.0),
        "valley_wall_found": wall_found,
        "bank_slope_deg": t[f"bank_slope_deg_{dem}"], "log10_upa": t["log10_upa"],
        "outer_bend": t["outer_bend"], "curvature": t["curvature"],
        "bank_elev_m": t[f"bank_elev_m_{dem}"],
    }, index=t.index)
    for name, col in outcome_cols.items():
        out[name] = pd.to_numeric(t[col], errors="coerce")
    out["y_logdist"] = np.log1p(out["y_dist"])
    a_ok = out[list(MODEL_COLUMNS) + ["y_height", "y_dist", "y_over"]].notna().all(axis=1)
    b_ok = out[list(MODEL_B_COLUMNS) + ["y_height", "y_dist"]].notna().all(axis=1)
    out["a_ok"], out["b_ok"] = a_ok, b_ok
    dropped_a, dropped_b = int((~a_ok).sum()), int((~b_ok).sum())
    out = out[a_ok | b_ok].copy()
    out.attrs["dropped_missing"] = dropped_a
    out.attrs["dropped_missing_b"] = dropped_b
    return out


def common_support(df: pd.DataFrame, controls_rhs: str, threshold: float,
                   min_share: float, min_row_fraction: float,
                   min_segments: int) -> dict:
    """Overlap trimming (8.1 with 14.8). Model A only; Model B is never trimmed.

    forest_high = F >= threshold, or F >= median(F) when either group has fewer
    than min_share of the rows. A logistic regression of forest_high on the
    controls gives p; rows with p inside the overlap of the two groups' ranges
    are kept. Too few rows or segments left means the comparison is too thin.
    """
    f = df["forest_frac"]
    split, rule = threshold, f"F >= {threshold:g}"
    high = f >= split
    if min(high.mean(), 1 - high.mean()) < min_share:
        split = float(f.median())
        rule = f"F >= median {split:.3f} (the {threshold:g} split left a group under {min_share:.0%})"
        high = f >= split
    if high.all() or not high.any():
        raise RuntimeError("Common support: one forest group is empty.")
    data = df.assign(forest_high=high.astype(int))
    model = smf.logit(f"forest_high ~ {controls_rhs}", data=data).fit(disp=0, maxiter=200)
    p = pd.Series(model.predict(data), index=df.index)
    lo = max(p[high].min(), p[~high].min())
    hi = min(p[high].max(), p[~high].max())
    keep = (p >= lo) & (p <= hi)
    kept_segments = int(df.loc[keep, "segment_id"].nunique())
    too_thin = bool(keep.mean() < min_row_fraction or kept_segments < min_segments)
    return {"rule": rule, "split": split, "p": p, "high": high, "keep": keep,
            "bounds": (float(lo), float(hi)), "rows": int(len(df)), "kept": int(keep.sum()),
            "dropped": int((~keep).sum()), "kept_segments": kept_segments,
            "too_thin": too_thin}


def paired_frame(df: pd.DataFrame, outcomes=("y_height", "y_logdist", "y_over")) -> pd.DataFrame:
    """Model B rows (8.3): left minus right at stations where both banks are in df."""
    left = df[df.side == "left"].set_index("station_id")
    right = df[df.side == "right"].set_index("station_id")
    both = left.index.intersection(right.index)
    left, right = left.loc[both], right.loc[both]
    out = pd.DataFrame({"station_id": both, "block_id": left["block_id"].to_numpy(),
                        "segment_id": left["segment_id"].to_numpy(),
                        "river_section": left["river_section"].to_numpy(),
                        "d_forest": (left.forest_frac - right.forest_frac).to_numpy()})
    for name in MODEL_B_DIFFS:
        out[f"d_{name}"] = (left[name] - right[name]).to_numpy()
    for name in outcomes:
        out[f"d_{name}"] = (left[name] - right[name]).to_numpy()
    return out


def sign_flips(primary: dict, runs: dict, keys) -> list[dict]:
    """8.10: sensitivity runs whose beta1 sign differs from the primary run's.

    primary and each run map (model, outcome) -> result dict with "beta" (or "error").
    Only pairs where both estimates exist are compared.
    """
    flips = []
    for name, run in runs.items():
        for key in keys:
            p, r = primary.get(key, {}), run.get(key, {})
            if "beta" not in p or "beta" not in r:
                continue
            if np.sign(p["beta"]) != np.sign(r["beta"]):
                flips.append({"run": name, "model": key[0], "outcome": key[1],
                              "primary_beta": p["beta"], "run_beta": r["beta"],
                              "run_p0": r.get("p0")})
    return flips


def model_b_rhs() -> str:
    return "d_forest + " + " + ".join(f"d_{n}" for n in MODEL_B_DIFFS)


def model_b_informative(pairs: pd.DataFrame, min_stations: int, min_diff: float) -> dict:
    """14.8: Model B is uninformative if fewer than min_stations have |F_L - F_R| >= min_diff."""
    n = int((pairs["d_forest"].abs() >= min_diff).sum())
    return {"stations": int(len(pairs)), "stations_large_diff": n,
            "informative": n >= min_stations,
            "abs_diff_quantiles": pairs["d_forest"].abs().quantile([0.25, 0.5, 0.75]).round(3).tolist()}


def aft_conditional(df: pd.DataFrame, outcome: str, rhs: str, term: str = "forest_frac") -> dict:
    """14.7 part 2: Weibull AFT on rows with reach > 0, right-censored where censored = 1.

    A sensitivity conditional on overtopping, not a headline result. Rows with a
    zero value of the outcome (a reach whose edge is not above z0) cannot enter a
    log-time model and are counted, not shifted (no constants are added to zeros).
    """
    from lifelines import WeibullAFTFitter

    rows = df[df["y_dist"] > 0]
    zero = int((rows[outcome] <= 0).sum())
    rows = rows[rows[outcome] > 0]
    x = patsy.dmatrix(rhs, rows, return_type="dataframe").drop(columns="Intercept")
    names = {c: f"x{i}" for i, c in enumerate(x.columns)}
    data = x.rename(columns=names)
    data = data.loc[:, data.std() > 0]
    data["T"] = rows.loc[data.index, outcome].to_numpy(float)
    data["E"] = 1 - rows.loc[data.index, "censored"].to_numpy(float)
    fit = WeibullAFTFitter().fit(data, duration_col="T", event_col="E", robust=True)
    s = fit.summary.loc[("lambda_", names[term])]
    return {"method": "Weibull AFT, right-censored, robust (not cluster) SE; conditional on reach > 0",
            "beta_log_time": float(s["coef"]), "time_ratio": float(np.exp(s["coef"])),
            "ci_time_ratio": (float(np.exp(s["coef lower 95%"])), float(np.exp(s["coef upper 95%"]))),
            "p": float(s["p"]), "rows": int(len(data)), "censored": int((data["E"] == 0).sum()),
            "zero_outcome_rows_left_out": zero}


def effect_sizes(df: pd.DataFrame, result: dict, forest_col: str = "forest_frac") -> dict:
    """8.6: beta for 0 -> 100 % forest and per interquartile range of F, with the interval."""
    q25, q75 = df[forest_col].quantile([0.25, 0.75])
    iqr = float(q75 - q25)
    lo, hi = result["interval"]
    return {"effect_0_to_100": result["beta"], "interval_0_to_100": (lo, hi), "iqr_F": iqr,
            "effect_per_iqr": result["beta"] * iqr, "interval_per_iqr": (lo * iqr, hi * iqr)}


def evalue_for(result: dict, sd_y: float) -> dict:
    """8.7 E-value for an OLS result, using the interval limit nearest the null."""
    lo, hi = result["interval"]
    nearest = 0.0 if lo <= 0 <= hi else (lo if abs(lo) < abs(hi) else hi)
    return evalue_continuous(result["beta"], sd_y, nearest)


def evalue_logistic(result: dict, prevalence: float, cutoff: float) -> dict:
    """14.9 E-value for overtopped, and for the OR limit nearest 1 (1 if the interval spans 1)."""
    lo, hi = result["or_interval"]
    point = evalue_binary(result["odds_ratio"], prevalence, cutoff)
    if lo <= 1 <= hi:
        ci = 1.0
    else:
        ci = evalue_binary(lo if lo > 1 else hi, prevalence, cutoff)
    return {"evalue": point, "evalue_ci": ci}


def or_narrow(or_interval: tuple[float, float], low: float, high: float) -> bool:
    """14.10: the whole odds-ratio interval lies within [low, high]."""
    return low <= or_interval[0] and or_interval[1] <= high


def decide_verdict(checks: dict) -> tuple[str, list[str]]:
    """The fixed decision rule (8.10) with the caps from 14.8, 14.13 and 14.15.

    checks holds booleans and short notes prepared by the runner:
      a_height_neg_sig, b_height_neg_sig, dist_same_sign, labels_same_sign,
      dems_same_sign, placebo_null, validation_ok, enough_blocks, overlap_ok,
      model_b_informative, validated, no_sign_flips, all_include_null, all_narrow.

    validated, validation_ok and placebo_null may instead be NOT_RUN when the event has no
    validation or no placebo yet. The verdict is then pending: a missing input is never read
    as a failed one. The other caps and unmet checks are still listed.
    """
    pending = [name for name, keys in PENDING_INPUTS
               if any(not_run(checks.get(k)) for k in keys)]
    caps, failed = _caps_and_unmet(checks)
    if pending:
        waiting = [f"{name} not yet run for this event; the verdict waits for it"
                   for name in pending]
        names = " and ".join(name.split(" (")[0] for name in pending)
        return f"Pending ({names} not yet run)", waiting + caps + failed
    if not caps and not failed:
        return "Consistent with a protective association", []
    if not caps and checks["all_include_null"] and checks["all_narrow"]:
        return "No evidence of association", []
    reasons = caps + failed
    if checks["all_include_null"] and not checks["all_narrow"]:
        reasons.append("intervals include the null but are not reasonably narrow (14.10)")
    return "Inconclusive", reasons


NOT_RUN = "not yet run"
PENDING_INPUTS = (("validation (14.13)", ("validated", "validation_ok")),
                  ("placebo (8.9)", ("placebo_null",)))


def not_run(value) -> bool:
    return isinstance(value, str) and value == NOT_RUN


def _caps_and_unmet(checks: dict) -> tuple[list[str], list[str]]:
    """The caps and the unmet protective checks; inputs not yet run count as neither."""
    def fails(key):
        return not not_run(checks[key]) and not checks[key]

    caps = []
    if not checks["enough_blocks"]:
        caps.append("fewer blocks than the minimum (14.15)")
    if not checks["overlap_ok"]:
        caps.append("overlap too thin (8.1)")
    if not checks["model_b_informative"]:
        caps.append("Model B uninformative (14.8)")
    if fails("validated"):
        caps.append("not validated (14.13)")
    if fails("placebo_null"):
        caps.append("placebo not null (8.9): the main result cannot be trusted")
    if not checks["no_sign_flips"]:
        caps.append("sign flips across sensitivity runs (8.10)")
    protective = [
        ("Model A reach height: beta < 0, interval excludes 0", "a_height_neg_sig"),
        ("Model B reach height: beta < 0, interval excludes 0", "b_height_neg_sig"),
        ("reach distance: same sign in Model A and B", "dist_same_sign"),
        ("same sign under optical-only and radar-only labels", "labels_same_sign"),
        ("same sign under both DEMs", "dems_same_sign"),
        ("placebo null", "placebo_null"),
        ("validation omission not materially higher on forested banks", "validation_ok"),
    ]
    return caps, [f"not met: {name}" for name, key in protective if fails(key)]


def fit_model_a(table: pd.DataFrame, outcome: str, formula_rhs: str):
    """OLS of one outcome on forest fraction plus controls.

    cov_type='nonrobust' skips the heavy default covariance calculation,
    which the block bootstrap replaces anyway.
    """
    return smf.ols(f"{outcome} ~ {formula_rhs}", data=table).fit(
        cov_type="nonrobust")


def fit_model_b(pairs: pd.DataFrame, outcome_diff: str, formula_rhs: str):
    """OLS of the left-minus-right outcome difference."""
    return smf.ols(f"{outcome_diff} ~ {formula_rhs}", data=pairs).fit(
        cov_type="nonrobust")


def _segment_groups(segment_ids: np.ndarray) -> list[np.ndarray]:
    """Row positions of each segment, in order of first appearance."""
    uniq = pd.unique(segment_ids)
    buckets: list[list[int]] = [[] for _ in range(len(uniq))]
    position = {s: i for i, s in enumerate(uniq)}
    for row, segment in enumerate(segment_ids):
        buckets[position[segment]].append(row)
    return [np.asarray(rows, dtype=np.intp) for rows in buckets]


def _bootstrap_ols_lstsq(table: pd.DataFrame, outcome: str, formula_rhs: str,
                         term: str, segment_col: str, reps: int,
                         seed: int) -> np.ndarray:
    """Block bootstrap of one OLS coefficient using one design matrix.

    The design matrix is built once. Each resample repeats whole segments and
    is solved with least squares, which is the same fit as refitting OLS.
    """
    y_df, x_df = patsy.dmatrices(
        f"{outcome} ~ {formula_rhs}", table, return_type="dataframe")
    if term not in x_df.columns:
        raise KeyError(
            f"Term {term!r} is not in the model. Columns: {list(x_df.columns)}")
    term_idx = list(x_df.columns).index(term)
    y = np.asarray(y_df, dtype=float).ravel()
    x = np.asarray(x_df, dtype=float)
    groups = _segment_groups(table.loc[x_df.index, segment_col].to_numpy())

    rng = np.random.default_rng(seed)
    n_segments = len(groups)
    draws = np.empty(reps)
    for b in range(reps):
        chosen = rng.integers(0, n_segments, size=n_segments)
        rows = np.concatenate([groups[i] for i in chosen])
        coef, *_ = np.linalg.lstsq(x[rows], y[rows], rcond=None)
        draws[b] = coef[term_idx]
    return draws


def block_bootstrap_beta(table: pd.DataFrame, outcome: str, formula_rhs: str,
                         term: str, segment_col: str = "segment_id",
                         reps: int = 2000, seed: int = 0,
                         fit=fit_model_a) -> np.ndarray:
    """Percentile bootstrap of one coefficient, resampling whole segments.

    All rows of a sampled segment move together, because neighbouring
    transects are not independent (8.4). The default of 2000 repetitions is
    the production setting. OLS refits are solved with numpy least squares on
    a single design matrix; other model types are refitted in full.
    """
    if fit in (fit_model_a, fit_model_b):
        return _bootstrap_ols_lstsq(
            table, outcome, formula_rhs, term, segment_col, reps, seed)

    rng = np.random.default_rng(seed)
    segments = table[segment_col].unique()
    draws = np.empty(reps)
    for b in range(reps):
        chosen = rng.choice(segments, size=len(segments), replace=True)
        parts = [table[table[segment_col] == s] for s in chosen]
        sample = pd.concat(parts, ignore_index=True)
        draws[b] = float(fit(sample, outcome, formula_rhs).params[term])
    return draws


def bootstrap_interval(draws: np.ndarray, level: float = 0.95) -> tuple[float, float]:
    tail = (1 - level) / 2 * 100
    return tuple(np.percentile(draws, [tail, 100 - tail]))


def pairs_percentile_interval(table: pd.DataFrame, outcome: str,
                              formula_rhs: str, term: str, reps: int,
                              seed: int, segment_col: str = "segment_id",
                              fit=fit_model_a) -> dict:
    """Sensitivity interval from 8.4, reported beside the primary one (14.15).

    1 km pairs block bootstrap with a percentile interval. It is labelled as
    too narrow with few segments, because it is.
    """
    draws = block_bootstrap_beta(table, outcome, formula_rhs, term,
                                 segment_col=segment_col, reps=reps,
                                 seed=seed, fit=fit)
    return {
        "method": "1 km pairs block bootstrap, percentile interval",
        "label": PERCENTILE_LABEL,
        "interval": bootstrap_interval(draws),
        "reps": reps,
    }


# ---------------------------------------------------------------------------
# Primary interval: wild cluster bootstrap-t with the null imposed (14.15)
# ---------------------------------------------------------------------------

def _cluster_indicator(clusters: np.ndarray) -> np.ndarray:
    """G x N matrix with a 1 where row n belongs to cluster g."""
    _, codes = np.unique(np.asarray(clusters), return_inverse=True)
    indicator = np.zeros((codes.max() + 1, len(codes)))
    indicator[codes, np.arange(len(codes))] = 1.0
    return indicator


def _cr1_factor(n_clusters: int, n_obs: int, n_params: int) -> float:
    """Small-sample factor of the CR1 cluster-robust variance."""
    return (n_clusters / (n_clusters - 1)) * ((n_obs - 1) / (n_obs - n_params))


class WildClusterBootstrap:
    """Wild cluster bootstrap-t for one OLS coefficient, null imposed.

    The design matrix and one set of Rademacher weights (one per cluster and
    draw) are fixed at construction, so p(b0) is a smooth function of b0 and
    the interval search is reproducible.
    """

    def __init__(self, table: pd.DataFrame, outcome: str, formula_rhs: str,
                 term: str, cluster_col: str, draws: int, seed: int):
        y_df, x_df = patsy.dmatrices(f"{outcome} ~ {formula_rhs}", table,
                                     return_type="dataframe")
        if term not in x_df.columns:
            raise KeyError(
                f"Term {term!r} is not in the model. Columns: {list(x_df.columns)}")
        self.term = term
        j = list(x_df.columns).index(term)
        self.y = np.asarray(y_df, dtype=float).ravel()
        x = np.asarray(x_df, dtype=float)
        self.f = x[:, j]
        x_rest = np.delete(x, j, axis=1)

        n_obs, n_params = x.shape
        self.indicator = _cluster_indicator(
            table.loc[x_df.index, cluster_col].to_numpy())
        self.n_clusters = self.indicator.shape[0]
        if self.n_clusters < 2:
            raise ValueError("The wild cluster bootstrap needs at least 2 clusters.")
        self.factor = _cr1_factor(self.n_clusters, n_obs, n_params)

        self.x = x
        self.proj = np.linalg.pinv(x)         # (X'X)^-1 X', K x N
        self.q = self.proj[j]                 # row that gives beta_term
        self.proj_rest = np.linalg.pinv(x_rest)
        self.x_rest = x_rest

        self.beta_hat = float(self.q @ self.y)
        self.se_cluster = float(self._se(self.y[:, None])[0])

        rng = np.random.default_rng(seed)
        weights = rng.choice([-1.0, 1.0], size=(self.n_clusters, draws))
        self.row_weights = self.indicator.T @ weights   # N x B, one sign per cluster

    def _se(self, y_mat: np.ndarray) -> np.ndarray:
        """CR1 standard error of the term for each column of y_mat."""
        resid = y_mat - self.x @ (self.proj @ y_mat)
        scores = self.indicator @ (self.q[:, None] * resid)  # G x B
        return np.sqrt(self.factor * np.sum(scores ** 2, axis=0))

    def p_value(self, b0: float) -> float:
        """Share of draws with |t*| >= |t_obs| when beta = b0 is imposed."""
        t_obs = (self.beta_hat - b0) / self.se_cluster
        restricted_y = self.y - b0 * self.f
        fitted = self.x_rest @ (self.proj_rest @ restricted_y)
        resid = restricted_y - fitted
        y_star = (b0 * self.f + fitted)[:, None] + self.row_weights * resid[:, None]
        beta_star = self.q @ y_star
        t_star = (beta_star - b0) / self._se(y_star)
        return float(np.mean(np.abs(t_star) >= abs(t_obs)))

    def interval(self, alpha: float = 0.05, tolerance_se: float = 0.01,
                 max_steps: int = 60) -> tuple[float, float]:
        """Set of b0 with p(b0) > alpha, by bisection outward from beta_hat."""
        tol = tolerance_se * self.se_cluster

        def edge(direction: int) -> float:
            inside = self.beta_hat
            step = self.se_cluster
            for _ in range(max_steps):
                outside = self.beta_hat + direction * step
                if self.p_value(outside) <= alpha:
                    break
                inside = outside
                step *= 2
            else:
                raise RuntimeError(
                    "Could not find the edge of the interval; p(b0) stayed above "
                    f"{alpha} out to {step:.3g} from the estimate.")
            while abs(outside - inside) > tol:
                mid = 0.5 * (inside + outside)
                if self.p_value(mid) > alpha:
                    inside = mid
                else:
                    outside = mid
            return 0.5 * (inside + outside)

        return edge(-1), edge(+1)


def wild_cluster_result(table: pd.DataFrame, outcome: str, formula_rhs: str,
                        term: str, cluster_col: str, draws: int, seed: int,
                        tolerance_se: float, alpha: float = 0.05) -> dict:
    """Primary result for one OLS coefficient: estimate, p(0) and interval.

    "Interval excludes 0" in 8.10 means p(0) <= alpha.
    """
    wb = WildClusterBootstrap(table, outcome, formula_rhs, term, cluster_col,
                              draws, seed)
    p0 = wb.p_value(0.0)
    return {
        "method": "wild cluster bootstrap-t, null imposed, Rademacher weights",
        "beta": wb.beta_hat,
        "se_cluster": wb.se_cluster,
        "p0": p0,
        "excludes_zero": p0 <= alpha,
        "interval": wb.interval(alpha=alpha, tolerance_se=tolerance_se),
        "n_clusters": wb.n_clusters,
        "draws": draws,
    }


def logistic_cluster_interval(table: pd.DataFrame, outcome: str,
                              formula_rhs: str, term: str, cluster_col: str,
                              level: float = 0.95) -> dict:
    """Approximate interval for the logistic model (14.15).

    Cluster-robust standard error on the blocks and a t critical value with
    (number of blocks - 1) degrees of freedom, reported as an odds ratio.
    """
    groups = pd.factorize(table[cluster_col])[0]
    model = smf.logit(f"{outcome} ~ {formula_rhs}", data=table).fit(
        disp=0, cov_type="cluster", cov_kwds={"groups": groups})
    n_clusters = int(groups.max() + 1)
    crit = scipy_stats.t.ppf(0.5 + level / 2, df=n_clusters - 1)
    beta = float(model.params[term])
    se = float(model.bse[term])
    return {
        "method": "cluster-robust SE with t critical value (approximate)",
        "beta": beta,
        "se_cluster": se,
        "odds_ratio": float(np.exp(beta)),
        "or_interval": (float(np.exp(beta - crit * se)),
                        float(np.exp(beta + crit * se))),
        "df": n_clusters - 1,
    }


def evalue_continuous(beta: float, sd_y: float,
                      ci_limit_nearest_null: float) -> dict:
    """E-value for a continuous outcome (8.7).

    An E-value of 2 means an unmeasured factor would need to be associated
    with both forest and flood reach by a risk ratio of 2 each to erase the
    result. The E-value for the confidence limit is 1 when the interval
    includes 0.
    """
    def _e(effect: float) -> float:
        d = effect / sd_y
        rr = np.exp(0.91 * abs(d))
        return float(rr + np.sqrt(rr * (rr - 1)))

    e_ci = 1.0 if (beta < 0) != (ci_limit_nearest_null < 0) or \
        ci_limit_nearest_null == 0 else _e(ci_limit_nearest_null)
    return {"evalue": _e(beta), "evalue_ci": e_ci}


def evalue_binary(odds_ratio: float, prevalence: float,
                  prevalence_cutoff: float = 0.15) -> float:
    """E-value for the yes/no outcome (14.9).

    Uses the odds ratio for forest going from 0 to 1. When overtopping is
    common (prevalence above the cutoff) the odds ratio is converted to a
    risk ratio with sqrt(OR); otherwise the odds ratio is used directly.
    """
    rr = np.sqrt(odds_ratio) if prevalence > prevalence_cutoff else odds_ratio
    rr_star = 1 / rr if rr < 1 else rr
    return float(rr_star + np.sqrt(rr_star * (rr_star - 1)))


def is_reasonably_narrow(ci: tuple[float, float], sd_y: float,
                         limit_sd: float, null: float = 0.0) -> bool:
    """True when the whole interval, scaled by SD(y), lies within +/- limit_sd.

    For the yes/no outcome pass sd_y=1 and the odds-ratio interval, with
    null=1 (14.10).
    """
    lo, hi = ci
    return abs(lo - null) / sd_y <= limit_sd and abs(hi - null) / sd_y <= limit_sd
