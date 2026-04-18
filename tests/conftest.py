"""
conftest.py — shared pytest fixtures for the test suite.

Fixtures defined here are automatically available to all test files in
this directory without an explicit import.

Fixtures:
    psic_geom       -- fresh psic() AdHocDiffractometer for each test
    reset_precision -- autouse; restores display precision after every test

For pure helper functions and constants used in parametrize() decorators
(which run at collection time, before fixtures are injected), see
helpers.py — import those explicitly with `from helpers import ...`.
"""

import pytest

from ad_hoc_diffractometer.display import set_precision
from ad_hoc_diffractometer.presets import psic

_DISPLAY_DEFAULT = 6  # matches display._DEFAULT_PRECISION


@pytest.fixture
def psic_geom():
    """Return a fresh psic() AdHocDiffractometer instance for each test."""
    return psic()


@pytest.fixture(autouse=True)
def reset_precision():
    """Restore package-level display precision to the default after every test."""
    yield
    set_precision(_DISPLAY_DEFAULT)
