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
                with_sections: bool = False) -> str:
    """Right-hand side of Model A (8.2), with river_section when needed (14.16).

    river_section enters as a fixed effect only when the reach has at least
    one confirmed section break; Model B never includes it.
    """
    rhs = (f"{forest} + np.log(channel_width_m) + gradient"
           " + np.log(valley_width_m) + wall_slope_deg + bank_slope_deg"
           " + log10_upa + outer_bend + curvature + bank_elev_m + C(side)"
           f" + cr(river_km, df={spline_df})")
    if with_sections:
        rhs += " + C(river_section)"
    return rhs


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
