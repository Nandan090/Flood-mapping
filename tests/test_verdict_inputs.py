"""Verdict inputs are per event: validation and placebo that have not been run are pending."""

from types import SimpleNamespace

import pandas as pd

from floodforest import stats as S
from floodforest.config import PROJECT_ROOT, load_config, validation_results_dir
from floodforest.stats_run import _validation_inputs, pending_lines


def _checks(**over):
    base = {k: True for k in ("a_height_neg_sig", "b_height_neg_sig", "dist_same_sign",
                              "labels_same_sign", "dems_same_sign", "placebo_null",
                              "validation_ok", "validated", "enough_blocks", "overlap_ok",
                              "model_b_informative", "no_sign_flips")}
    base.update(all_include_null=False, all_narrow=False)
    base.update(over)
    return base


NOT_RUN_VALIDATION = dict(validated=S.NOT_RUN, validation_ok=S.NOT_RUN)


def _cfg(results):
    return SimpleNamespace(paths=SimpleNamespace(results=str(results)),
                           validation_min_group_n=8, validation_recall_gap=0.10)


def test_validation_folder_is_per_event():
    cfg = load_config()
    assert validation_results_dir(cfg, "melamchi2021") == PROJECT_ROOT / "results"
    assert validation_results_dir(cfg, "lende2025") == PROJECT_ROOT / "results" / "lende2025"


def test_missing_inputs_give_a_pending_verdict_not_a_failure():
    verdict, reasons = S.decide_verdict(_checks(placebo_null=S.NOT_RUN, **NOT_RUN_VALIDATION))
    assert verdict == "Pending (validation and placebo not yet run)"
    assert reasons[:2] == [
        "validation (14.13) not yet run for this event; the verdict waits for it",
        "placebo (8.9) not yet run for this event; the verdict waits for it"]
    assert not any("cannot be trusted" in r or "not validated" in r or "placebo null" in r
                   or "validation omission" in r for r in reasons)
    assert S.decide_verdict(_checks(placebo_null=S.NOT_RUN))[0] == (
        "Pending (placebo not yet run)")
    assert S.decide_verdict(_checks(**NOT_RUN_VALIDATION))[0] == (
        "Pending (validation not yet run)")


def test_pending_still_lists_the_other_caps_and_unmet_checks():
    verdict, reasons = S.decide_verdict(_checks(placebo_null=S.NOT_RUN, no_sign_flips=False,
                                                b_height_neg_sig=False))
    assert verdict.startswith("Pending")
    assert "sign flips across sensitivity runs (8.10)" in reasons
    assert "not met: Model B reach height: beta < 0, interval excludes 0" in reasons


def test_inputs_that_were_run_keep_the_old_rule():
    """A failed placebo or validation is still a cap, as before (Melamchi's case)."""
    verdict, reasons = S.decide_verdict(_checks(placebo_null=False))
    assert verdict == "Inconclusive" and any("cannot be trusted" in r for r in reasons)
    verdict, reasons = S.decide_verdict(_checks(validated=False, validation_ok=False))
    assert verdict == "Inconclusive" and "not validated (14.13)" in reasons
    assert S.decide_verdict(_checks())[0] == "Consistent with a protective association"


def test_an_event_without_validation_files_is_not_run(tmp_path):
    (tmp_path / "validation_summary.csv").write_text("other event's file", encoding="utf-8")
    (tmp_path / "validation_per_point.csv").write_text("other event's file", encoding="utf-8")
    assert _validation_inputs(_cfg(tmp_path), "lende2025") == (S.NOT_RUN,) * 3


def test_an_event_reads_its_own_validation_files(tmp_path):
    folder = tmp_path / "lende2025"
    folder.mkdir()
    pd.DataFrame({"subset": ["main", "main"], "breakdown": ["forest_group"] * 2,
                  "category": ["Low", "High"], "n": [10, 10],
                  "share_under_detected": [0.1, 0.15]}).to_csv(
        folder / "validation_summary.csv", index=False)
    ok, note, validated = _validation_inputs(_cfg(tmp_path), "lende2025")
    assert ok is True or ok == True  # noqa: E712 (numpy bool)
    assert note.startswith("under-detected High 15.0% vs Low 10.0%")
    assert validated is False      # no per-point file


def test_melamchi_still_reads_the_shared_results_files():
    cfg = load_config()
    if not (PROJECT_ROOT / "results" / "validation_summary.csv").exists():
        import pytest
        pytest.skip("Melamchi validation results are not present")
    ok, note, validated = _validation_inputs(cfg, "melamchi2021")
    assert note.startswith("under-detected High 26.7% vs Low 10.0%")
    assert validated is True
    assert _validation_inputs(cfg, "lende2025") == (S.NOT_RUN,) * 3


def test_pending_note_in_the_report():
    assert pending_lines(_checks()) == []
    assert pending_lines(_checks(placebo_null=False)) == []
    lines = pending_lines(_checks(placebo_null=S.NOT_RUN, **NOT_RUN_VALIDATION))
    assert lines[0] == ("**Not yet run for this event: validation (14.13) and placebo (8.9). The "
                        "verdict is pending until they are run; their checks below read "
                        "'not yet run'.**")


def test_pending_note_names_only_what_is_missing():
    """Placebo run, validation not: the note names validation alone, in the singular."""
    only_validation = pending_lines(_checks(placebo_null=False, **NOT_RUN_VALIDATION))
    assert only_validation[0] == ("**Not yet run for this event: validation (14.13). The verdict "
                                  "is pending until it is run; its checks below read "
                                  "'not yet run'.**")
    assert "placebo" not in only_validation[0]
    only_placebo = pending_lines(_checks(placebo_null=S.NOT_RUN))
    assert only_placebo[0].startswith("**Not yet run for this event: placebo (8.9). The verdict "
                                      "is pending until it is run;")
    assert "validation" not in only_placebo[0]
