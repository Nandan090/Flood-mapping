"""The stats report prints the Model A formula actually fitted (after trimming)."""

import json

import pandas as pd
import pytest

from floodforest import stats as S
from floodforest.config import PROJECT_ROOT, load_config
from floodforest.stats_run import analysis_rows, model_a_formula_line, rhs_for

RHS = ("forest_frac + np.log(channel_width_m) + gradient + log10_upa + bank_elev_m"
       " + C(side) + cr(river_km, df=5) + valley_wall_found + channel_resolved")
OLD_LINE = (f"Model A (8.2, 14.16, 14.18): `y ~ {RHS}`. OLS for reach height and "
            "log1p(reach distance); logistic for overtopped.")


def test_formula_line_is_unchanged_when_nothing_is_left_out():
    assert model_a_formula_line({"rhs_a": RHS, "_rhs_a_fitted": RHS}) == OLD_LINE
    assert model_a_formula_line({"rhs_a": RHS}) == OLD_LINE      # older summaries


def test_formula_line_prints_the_fitted_formula_and_names_constant_terms():
    fitted = RHS.replace(" + valley_wall_found", "")
    line = model_a_formula_line({"rhs_a": RHS, "_rhs_a_fitted": fitted})
    assert line.startswith(f"Model A (8.2, 14.16, 14.18): `y ~ {fitted}`.")
    assert "`valley_wall_found` is in the formula before trimming but left out of the fit" in line


def _trimmed_rows(event_key):
    cfg = load_config()
    root = PROJECT_ROOT / cfg.paths.outputs / event_key
    if not (root / "stats" / "stats_summary.json").exists():
        pytest.skip(f"{event_key} stats outputs are not present")
    base = analysis_rows(pd.read_csv(root / "transects.csv"),
                         pd.read_csv(root / "transects_labels.csv"), cfg)
    df = S.model_frame(base[base["in_set"]], f"forest_frac_{int(cfg.belt_width_m)}", "fab", None)
    sections = bool(getattr(cfg.events, event_key).section_breaks_km)
    a_all = df[df["a_ok"]]
    cs = S.common_support(a_all, rhs_for(a_all, cfg, sections).split(" + ", 1)[1],
                          float(cfg.forest_high_threshold), float(cfg.forest_high_min_share),
                          float(cfg.overlap_min_row_fraction), int(cfg.overlap_min_segments))
    summary = json.loads((root / "stats" / "stats_summary.json").read_text(encoding="utf-8"))
    return cfg, root, a_all[cs["keep"]], sections, summary


@pytest.mark.filterwarnings("ignore")
def test_melamchi_report_formula_line_stays_the_same():
    cfg, root, a_df, sections, summary = _trimmed_rows("melamchi2021")
    primary = {"rhs_a": summary["primary"]["rhs_a"],
               "_rhs_a_fitted": rhs_for(a_df, cfg, sections)}
    report = (root / "stats" / "STATS_REPORT.md").read_text(encoding="utf-8")
    assert model_a_formula_line(primary) in report.splitlines()


@pytest.mark.filterwarnings("ignore")
def test_lende_formula_names_valley_wall_found():
    cfg, _, a_df, sections, summary = _trimmed_rows("lende2025")
    line = model_a_formula_line({"rhs_a": summary["primary"]["rhs_a"],
                                 "_rhs_a_fitted": rhs_for(a_df, cfg, sections)})
    assert "`valley_wall_found` is in the formula before trimming" in line
