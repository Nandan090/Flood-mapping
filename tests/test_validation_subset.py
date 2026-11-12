"""The one-time follow-up subset of the validation sample (lende2025)."""

import pandas as pd
import pytest

from floodforest.config import load_config
from floodforest.validate import subset_draw

GROUPS = pd.Series({**{s: "low" for s in range(1, 11)}, **{s: "mid" for s in range(11, 21)},
                    **{s: "high" for s in range(21, 31)}})


def test_the_subset_tops_every_group_up_to_the_target():
    done = [1, 2, 3, 4, 11, 21]
    picked = subset_draw(GROUPS, done, 6, seed=43)
    assert picked == sorted(picked) and not set(picked) & set(done)
    every = done + picked
    assert [int(GROUPS[every].eq(g).sum()) for g in ("low", "mid", "high")] == [6, 6, 6]
    assert subset_draw(GROUPS, done, 6, seed=43) == picked       # same seed, same draw


def test_the_subset_refuses_what_it_cannot_balance():
    with pytest.raises(RuntimeError, match="more than"):
        subset_draw(GROUPS, [1, 2, 3], 2, seed=43)
    with pytest.raises(RuntimeError, match="only"):
        subset_draw(GROUPS, [], 11, seed=43)
    with pytest.raises(RuntimeError, match="not in the validation sample"):
        subset_draw(GROUPS, [99], 5, seed=43)


def test_lende_subset_settings_as_decided():
    sub = load_config().events.lende2025.validation.subset
    assert sub.done_before == [52, 53, 54, 55, 59, 60, 62, 66, 68, 69, 74]
    assert (sub.per_group, sub.seed) == (12, 43)
