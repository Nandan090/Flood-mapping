"""Report figures. So far: the six example cross-sections (spec section 10, item 6)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

PERCENTILES = (5, 25, 50, 75, 95)


def choose_examples(table, column: str = "reach_height_m_fab",
                    censored_column: str = "censored_fab") -> list[tuple[str, int]]:
    """Fixed rule, not hand-picked (10.6).

    Analysis-set rows nearest the 5th, 25th, 50th, 75th and 95th percentile
    of reach height, then the censored row nearest the median height of the
    censored rows. Ties go to the lowest river_km, left bank first. No row is
    used twice. Returns (rule, table index) pairs.
    """
    pool = table[(table["analysis_set"] == 1) & np.isfinite(table[column])]
    if pool.empty:
        raise RuntimeError("The analysis set has no reach height, so no example was chosen.")
    chosen: list[tuple[str, int]] = []

    def nearest(candidates, target):
        ranked = candidates.assign(_gap=(candidates[column] - target).abs()).sort_values(
            ["_gap", "river_km", "side"], kind="mergesort")
        return int(ranked.index[0])

    for p in PERCENTILES:
        target = float(np.percentile(pool[column], p))
        free = pool[~pool.index.isin([i for _, i in chosen])]
        chosen.append((f"{p}th percentile of reach height ({target:.1f} m)", nearest(free, target)))
    censored = pool[pool[censored_column] == 1]
    free = censored[~censored.index.isin([i for _, i in chosen])]
    if free.empty:
        raise RuntimeError("No censored row is left in the analysis set for the sixth example.")
    target = float(np.median(censored[column]))
    chosen.append((f"censored, nearest the median censored height ({target:.1f} m)",
                   nearest(free, target)))
    return chosen


def example_cross_sections(table, layers: dict, forest: dict, grid: dict, cfg,
                           event_key: str) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    from .centerline import _sample
    from .config import PROJECT_ROOT
    from .profiles import station_profile

    out_dir = PROJECT_ROOT / cfg.paths.outputs / event_key / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    step = cfg.grid_m
    written = []
    for n, (rule, index) in enumerate(choose_examples(table), start=1):
        row = table.loc[index]
        pair = table[table.station_id == row.station_id]
        rows = {r.side: r for r in pair.itertuples(index=False)}
        points, offsets = station_profile(rows, step)
        z_fab = _sample(layers["z_fab"], grid, points)
        z_cop = _sample(layers["z_cop"], grid, points)
        channel = _sample(layers["channel"], grid, points)
        impact = _sample(layers["combined"], grid, points)
        trees = _sample(forest["forest"], grid, points)
        x = -offsets  # left bank on the left, looking downstream

        sign = 1.0 if row.side == "left" else -1.0
        edge = row.reach_edge_from_station_m_fab
        wall = (row.bank_top_from_station_m + row.valley_width_m_fab
                if row.valley_wall_found_fab == 1 else np.nan)
        reach_out = np.nanmax([300.0, (edge if np.isfinite(edge) else 0) + 150,
                               (wall if np.isfinite(wall) else 0) + 100])
        lim = min(reach_out, float(np.max(np.abs(x))))

        fig, (ax, strip) = plt.subplots(
            2, 1, figsize=(11, 5.6), sharex=True,
            gridspec_kw={"height_ratios": [4, 1.1], "hspace": 0.08})
        ax.plot(x, z_fab, color="black", lw=1.4, label="FABDEM elevation")
        ax.plot(x, z_cop, color="grey", lw=1.0, ls="--", label="Copernicus GLO-30 elevation")
        ax.axhline(row.z0_fab, color="#1f78b4", lw=1, ls=":", label="z0 (channel elevation)")
        for side_name, r in rows.items():
            s = 1.0 if side_name == "left" else -1.0
            ax.axvline(-s * r.bank_top_from_station_m, color="#1f78b4", lw=1, ls="-.",
                       label="bank top" if side_name == "left" else None)
        if np.isfinite(wall):
            ax.axvline(-sign * wall, color="#8c510a", lw=1.2, label=f"valley wall ({row.side} bank)")
        if np.isfinite(edge):
            ze = float(np.interp(-sign * edge, x[::-1], z_fab[::-1]))
            ax.plot([-sign * edge], [ze], "o", color="#e31a1c", ms=7,
                    label=f"reach edge ({row.side} bank)")
        visible = np.abs(x) <= lim
        zmin = np.nanmin(z_fab[visible])
        zmax = np.nanmax(z_fab[visible])
        ax.set_ylim(zmin - 10, zmax + 10)
        ax.set_ylabel("elevation (m)")
        ax.legend(loc="upper center", fontsize=8, ncol=3, framealpha=0.9)
        reason = row.censor_reason_fab if isinstance(row.censor_reason_fab, str) else "not censored"
        ax.set_title(
            f"Example {n}: {rule}\nstation {row.station_id}, river km {row.river_km:.2f}, "
            f"{row.side} bank. Reach distance {row.reach_distance_m_fab:.0f} m, reach height "
            f"{row.reach_height_m_fab:.1f} m, {reason}. Usable length {row.usable_length_m:.0f} m.",
            fontsize=10)

        codes = np.full((3, len(x)), 0)
        codes[0] = np.where(channel == 1, 1, 0)
        codes[1] = np.where(np.isnan(impact), 2, np.where(impact == 1, 3, 0))
        codes[2] = np.where(np.isnan(trees), 2, np.where(trees == 1, 4, 0))
        cmap = ListedColormap(["#f0f0f0", "#6baed6", "#bdbdbd", "#e31a1c", "#238b45"])
        order = np.argsort(x)
        edges = np.concatenate([[x[order][0] - step / 2], x[order] + step / 2])
        strip.pcolormesh(edges, [0, 1, 2, 3], codes[:, order], cmap=cmap, vmin=-0.5, vmax=4.5)
        strip.set_yticks([0.5, 1.5, 2.5])
        strip.set_yticklabels(["channel", "impact", "forest"], fontsize=8)
        strip.invert_yaxis()
        strip.set_xlim(-lim, lim)
        strip.set_xlabel("distance from the station along the transect (m); "
                         "left bank to the left, looking downstream")
        strip.legend(handles=[Patch(color="#6baed6", label="pre-event channel"),
                              Patch(color="#e31a1c", label="impacted (combined label)"),
                              Patch(color="#238b45", label="pre-event forest pixel"),
                              Patch(color="#bdbdbd", label="no data")],
                     loc="upper center", bbox_to_anchor=(0.5, -0.55), ncol=4, fontsize=8)
        path = out_dir / f"cross_section_{n}.png"
        fig.savefig(path, dpi=130, bbox_inches="tight")
        plt.close(fig)
        written.append(path)
    return written
