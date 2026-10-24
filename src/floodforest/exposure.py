"""Forest belt metrics (spec section 7.3, clarifications 14.3, 14.14, 14.21)."""

from __future__ import annotations

import numpy as np

NODATA = -9999.0


def forest_mask(trees_probability: np.ndarray, canopy_height_m: np.ndarray,
                hansen_loss_year: np.ndarray, pre_event_year: int,
                trees_probability_min: float = 0.5,
                canopy_height_min_m: float = 5.0,
                hansen_loss_from_year: int = 2021,
                hansen_newest_year: int | None = None,
                use_canopy_height: bool = True) -> np.ndarray:
    """Pixels that count as forest.

    A pixel is forest when the Dynamic World trees probability and (unless
    disabled) the canopy height both clear their thresholds, and Hansen shows
    no loss with year in [hansen_loss_from_year, pre_event_year].

    The Hansen range is empty when the pre-event year is before 2021, so the
    condition is skipped (14.3, the Melamchi case). When Hansen's newest year
    is earlier than the pre-event year, loss is applied through the newest
    available year and the caller must log the uncovered years.
    """
    forest = np.asarray(trees_probability, dtype=float) >= trees_probability_min
    if use_canopy_height:
        forest &= np.asarray(canopy_height_m, dtype=float) >= canopy_height_min_m

    loss_year = np.asarray(hansen_loss_year)
    loss_to = pre_event_year
    if hansen_newest_year is not None:
        loss_to = min(loss_to, hansen_newest_year)
    if loss_to >= hansen_loss_from_year:
        forest &= ~((loss_year >= hansen_loss_from_year) & (loss_year <= loss_to))
    return forest


def hansen_loss_window(pre_event_year: int, newest_year: int | None,
                       from_year: int = 2021) -> dict:
    """How section 14.3 applies Hansen loss for one pre-composite year.

    The condition is loss year in [from_year, pre_event_year]. When the
    pre-composite year is before from_year the range is empty and the
    condition is skipped. When Hansen's newest year stops short of the
    pre-composite year, loss is applied through that newest year and the
    uncovered years are reported. That shortfall is not an error.
    """
    if pre_event_year < from_year:
        return {
            "applied": False,
            "loss_through_year": None,
            "uncovered_years": [],
            "message": (
                f"Hansen loss condition skipped: the pre-composite year is "
                f"{pre_event_year}, which is before {from_year}, so the range "
                f"[{from_year}, {pre_event_year}] is empty."
            ),
        }
    if newest_year is None:
        raise ValueError(
            "Hansen's newest year is unknown, so the loss window cannot be applied."
        )
    if newest_year < from_year:
        raise ValueError(
            f"Hansen's newest year ({newest_year}) is before {from_year}, so "
            "there is no loss layer to apply."
        )
    if newest_year < pre_event_year:
        uncovered = list(range(newest_year + 1, pre_event_year + 1))
        return {
            "applied": True,
            "loss_through_year": newest_year,
            "uncovered_years": uncovered,
            "message": (
                f"Hansen loss applied through {newest_year}. Uncovered years "
                f"{uncovered[0]}–{uncovered[-1]} are left to Dynamic World."
            ),
        }
    return {
        "applied": True,
        "loss_through_year": pre_event_year,
        "uncovered_years": [],
        "message": (
            f"Hansen loss applied for years {from_year} through {pre_event_year}."
        ),
    }


def forest_fraction(forest: np.ndarray, valid: np.ndarray | None = None) -> float:
    """Forest share of the valid belt pixels. Raises if the belt is empty."""
    mask = np.asarray(forest, dtype=bool)
    if valid is not None:
        keep = np.asarray(valid, dtype=bool)
        if not keep.any():
            raise ValueError("Belt has no valid pixels; forest fraction is undefined.")
        mask = mask[keep]
    if mask.size == 0:
        raise ValueError("Belt is empty; forest fraction is undefined.")
    return float(mask.mean())


SOURCE_DYNAMIC_WORLD = 0
SOURCE_WORLDCOVER = 1


def forest_with_unknowns(trees_probability: np.ndarray, canopy_height_m: np.ndarray,
                         hansen_loss_year: np.ndarray, pre_event_year: int,
                         trees_probability_min: float, canopy_height_min_m: float,
                         hansen_loss_from_year: int, hansen_newest_year: int | None,
                         worldcover_class: np.ndarray | None = None,
                         worldcover_forest_class: int = 10) -> tuple[np.ndarray, np.ndarray]:
    """Forest rule of 7.3 as 1 forest, 0 non-forest, NaN unknown, and the
    source of each decided pixel (0 Dynamic World, 1 WorldCover, NaN none).

    A pixel below the Dynamic World threshold is non-forest whatever its
    canopy height. Where Dynamic World has no observation, WorldCover class
    10 takes the place of the tree probability test, with the same canopy
    height and Hansen conditions (14.23 a). A pixel is unknown when neither
    map has it, or when the tree test passes and the canopy height is missing.
    """
    trees = np.asarray(trees_probability, dtype=float)
    height = np.asarray(canopy_height_m, dtype=float)
    has_dw = np.isfinite(trees)
    tree_test = trees >= trees_probability_min
    source = np.where(has_dw, float(SOURCE_DYNAMIC_WORLD), np.nan)
    if worldcover_class is not None:
        wc = np.asarray(worldcover_class, dtype=float)
        use_wc = ~has_dw & np.isfinite(wc)
        tree_test = np.where(use_wc, wc == worldcover_forest_class, tree_test)
        source[use_wc] = SOURCE_WORLDCOVER
    decided = np.isfinite(source)
    proxy = np.where(tree_test, 1.0, 0.0)
    forest = forest_mask(proxy, np.nan_to_num(height, nan=0.0), hansen_loss_year,
                         pre_event_year, 0.5, canopy_height_min_m,
                         hansen_loss_from_year, hansen_newest_year)
    out = forest.astype(float)
    unknown = ~decided | (tree_test & ~np.isfinite(height))
    out[unknown] = np.nan
    source[np.isnan(out)] = np.nan
    return out, source


def belt_polygon(station_xy, along_unit, normal_unit, start_m: float, end_m: float,
                 length_m: float):
    """Belt rectangle: length_m along the river, centred on the transect, from
    start_m to end_m out along the normal. None when end_m <= start_m."""
    from shapely.geometry import Polygon
    if end_m <= start_m:
        return None
    origin = np.asarray(station_xy, dtype=float)
    half = np.asarray(along_unit, dtype=float) * length_m / 2.0
    near = origin + np.asarray(normal_unit, dtype=float) * start_m
    far = origin + np.asarray(normal_unit, dtype=float) * end_m
    return Polygon([near - half, far - half, far + half, near + half])


def belt_pixel_indices(polygons: list, shape: tuple[int, int], transform) -> list[np.ndarray]:
    """Flat grid indices of the pixels each belt covers, found with rasterstats.

    A pixel belongs to a belt when its centre is inside it. A missing belt
    (None) gets no pixels.
    """
    from rasterstats import zonal_stats
    index = np.arange(shape[0] * shape[1], dtype=np.float64).reshape(shape)
    out = [np.array([], dtype=np.int64) for _ in polygons]
    present = [i for i, poly in enumerate(polygons) if poly is not None]
    if not present:
        return out
    stats = zonal_stats([polygons[i] for i in present], index, affine=transform,
                        nodata=-1, stats=["count"], raster_out=True, all_touched=False)
    for i, item in zip(present, stats):
        mini = item["mini_raster_array"]
        out[i] = np.asarray(np.ma.compressed(mini), dtype=np.int64)
    return out


def belt_metrics(index: np.ndarray, layers: dict, width: int) -> dict:
    """Exposure columns for one belt (7.3, 14.14).

    `layers` holds flattened grids: forest, forest_source (0 Dynamic World,
    1 WorldCover), forest_dw (1/0/NaN), canopy_height,
    ndvi_pre, prior_loss (1/0/NaN) and worldcover_forest (1/0/NaN).
    """
    from skimage.measure import label as components
    nan = float("nan")
    forest = layers["forest"][index]
    valid = np.isfinite(forest)
    source = layers["forest_source"][index][valid]
    wc_share = float(np.mean(source == SOURCE_WORLDCOVER)) if valid.any() else nan
    if not valid.any():
        kind = "none"
    elif wc_share == 0:
        kind = "dynamic_world"
    elif wc_share == 1:
        kind = "worldcover"
    else:
        kind = "mixed"
    out = {"belt_pixels": int(valid.sum()), "forest_source": kind,
           "forest_source_worldcover_frac": wc_share, "forest_frac": nan,
           "forest_frac_dw_only": nan, "forest_definition_disagree_frac": nan,
           "forest_frac_worldcover": nan, "canopy_height_mean": nan,
           "ndvi_pre_mean": nan, "belt_continuity": nan, "prior_loss_frac": nan}
    if valid.any():
        is_forest = forest[valid] == 1
        out["forest_frac"] = forest_fraction(is_forest)
        out["canopy_height_mean"] = float(np.nanmean(layers["canopy_height"][index][valid]))
        out["prior_loss_frac"] = float(np.nanmean(layers["prior_loss"][index][valid]))
        dw = layers["forest_dw"][index]
        both = valid & np.isfinite(dw)
        if both.any():
            out["forest_definition_disagree_frac"] = float(np.mean(forest[both] != dw[both]))
        if is_forest.any():
            rows, cols = np.divmod(index[valid][is_forest], width)
            grid = np.zeros((rows.max() - rows.min() + 1, cols.max() - cols.min() + 1), bool)
            grid[rows - rows.min(), cols - cols.min()] = True
            sizes = np.bincount(components(grid, connectivity=2).ravel())[1:]
            out["belt_continuity"] = float(sizes.max() / sizes.sum())
    dw = layers["forest_dw"][index]
    if np.isfinite(dw).any():
        out["forest_frac_dw_only"] = float(np.nanmean(dw))
    wc = layers["worldcover_forest"][index]
    if np.isfinite(wc).any():
        out["forest_frac_worldcover"] = float(np.nanmean(wc))
    ndvi = layers["ndvi_pre"][index]
    if np.isfinite(ndvi).any():
        out["ndvi_pre_mean"] = float(np.nanmean(ndvi))
    return out


def forest_layers(event_key: str, cfg, provenance: dict, grid: dict) -> dict:
    """Pre-event exposure grids for the belts (7.3, 14.3, 14.14)."""
    from . import optical
    from .config import PROJECT_ROOT
    from .label import _read_named

    rasters = PROJECT_ROOT / cfg.paths.data / event_key / "rasters"
    trees = _read_named(rasters / "dynamic_world_trees_pre.tif", ("trees",),
                        cfg.crs, grid)[0]["trees"]
    height = _read_named(rasters / "canopy_height.tif", ("canopy_height_m",),
                         cfg.crs, grid)[0]["canopy_height_m"]
    loss = _read_named(rasters / "hansen_loss_year.tif", ("loss_year",),
                       cfg.crs, grid)[0]["loss_year"]
    hansen = provenance["hansen"]
    exported = provenance.get("worldcover_exported", {})
    if exported.get("year") != hansen["pre_composite_year"]:
        raise RuntimeError(
            f"worldcover.tif is the {exported.get('year')} map, not the pre-event year "
            f"{hansen['pre_composite_year']}, so it cannot stand in for Dynamic World.")
    wc = _read_named(rasters / "worldcover.tif", ("class",), cfg.crs, grid)[0]["class"]
    forest, source = forest_with_unknowns(
        trees, height, np.nan_to_num(loss, nan=0.0), hansen["pre_composite_year"],
        cfg.trees_probability_min, cfg.canopy_height_min_m, cfg.hansen_loss_from_year,
        hansen["newest_year"], wc, cfg.worldcover_forest_class)
    pre = _read_named(rasters / "s2_pre.tif", ("B4", "B8"), cfg.crs, grid)[0]

    forest_dw = (trees >= cfg.trees_probability_min).astype(float)
    forest_dw[~np.isfinite(trees)] = np.nan
    year = provenance["hansen"]["pre_composite_year"]
    loss_year = np.nan_to_num(loss, nan=0.0)
    prior_loss = ((loss_year > 0) & (loss_year <= year)).astype(float)
    prior_loss[~np.isfinite(forest)] = np.nan
    wc_forest = (wc == cfg.worldcover_forest_class).astype(float)
    wc_forest[~np.isfinite(wc)] = np.nan
    return {
        "forest": forest,
        "forest_source": source,
        "forest_dw": forest_dw,
        "canopy_height": height,
        "ndvi_pre": optical.ndvi(pre["B8"], pre["B4"]),
        "prior_loss": prior_loss,
        "worldcover_forest": wc_forest,
    }
