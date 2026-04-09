"""
Unit tests for ad_hoc_diffractometer.display.

Covers:
  - get_precision()
  - set_precision()
  - fmt()

The reset_precision fixture in conftest.py restores the package-level
precision to its default after every test automatically.
"""

from contextlib import nullcontext as does_not_raise

import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import fmt
from ad_hoc_diffractometer import get_precision
from ad_hoc_diffractometer import set_precision

_DEFAULT = 6  # matches conftest._DISPLAY_DEFAULT


# ---------------------------------------------------------------------------
# get_precision() and set_precision()
# ---------------------------------------------------------------------------


def test_get_precision_default():
    assert get_precision() == _DEFAULT


@pytest.mark.parametrize(
    "digits, context",
    [
        pytest.param(0, does_not_raise(), id="zero-digits"),
        pytest.param(1, does_not_raise(), id="one-digit"),
        pytest.param(4, does_not_raise(), id="four-digits"),
        pytest.param(10, does_not_raise(), id="ten-digits"),
        pytest.param(
            -1,
            pytest.raises(ValueError, match="non-negative"),
            id="invalid-negative",
        ),
        pytest.param(
            1.5,
            pytest.raises(ValueError, match="non-negative integer"),
            id="invalid-float",
        ),
        pytest.param(
            "4",
            pytest.raises(ValueError, match="non-negative integer"),
            id="invalid-string",
        ),
        pytest.param(
            True,
            pytest.raises(ValueError, match="non-negative integer"),
            id="invalid-bool",
        ),
    ],
)
def test_set_precision(digits, context):
    with context:
        set_precision(digits)
        assert get_precision() == digits


def test_set_precision_affects_package():
    """set_precision updates the package-level default seen by ahd.get_precision()."""
    ahd.set_precision(3)
    assert ahd.get_precision() == 3


# ---------------------------------------------------------------------------
# fmt()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value, digits, expected, context",
    [
        pytest.param(4.785, 6, "4.785000", does_not_raise(), id="six-decimal-places"),
        pytest.param(4.785, 3, "4.785", does_not_raise(), id="three-decimal-places"),
        pytest.param(4.785, 0, "5", does_not_raise(), id="zero-decimal-places"),
        pytest.param(0.0, 4, "0.0000", does_not_raise(), id="zero-value"),
        pytest.param(-1.5, 2, "-1.50", does_not_raise(), id="negative-value"),
        pytest.param(1.0 / 3.0, 4, "0.3333", does_not_raise(), id="repeating-decimal"),
    ],
)
def test_fmt_explicit_digits(value, digits, expected, context):
    with context:
        assert fmt(value, digits) == expected


def test_fmt_uses_package_default():
    """fmt() with no digits argument uses the current package-level precision."""
    set_precision(3)
    assert fmt(1.23456) == "1.235"
