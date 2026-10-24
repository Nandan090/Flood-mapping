import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def pytest_collection_modifyitems(config, items):
    """Skip slow tests unless the command already selects markers.

    A plain pytest run skips them. pytest -m slow runs only them, because
    passing -m sets a marker expression and this hook then does nothing.
    """
    if config.option.markexpr:
        return
    skip_slow = pytest.mark.skip(reason="slow test; run with: pytest -m slow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
