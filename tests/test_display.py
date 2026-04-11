# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
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
from ad_hoc_diffractometer import allclose
from ad_hoc_diffractometer import fmt
from ad_hoc_diffractometer import get_precision
from ad_hoc_diffractometer import precision_atol
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


# ---------------------------------------------------------------------------
# precision_atol() and allclose()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "digits, expected_atol, context",
    [
        pytest.param(6, 5e-7, does_not_raise(), id="6-digits"),
        pytest.param(4, 5e-5, does_not_raise(), id="4-digits"),
        pytest.param(3, 5e-4, does_not_raise(), id="3-digits"),
        pytest.param(0, 0.5, does_not_raise(), id="0-digits"),
    ],
)
def test_precision_atol(digits, expected_atol, context):
    with context:
        assert precision_atol(digits) == pytest.approx(expected_atol)


def test_precision_atol_uses_package_default():
    set_precision(4)
    assert precision_atol() == pytest.approx(5e-5)


@pytest.mark.parametrize(
    "a, b, digits, expected, context",
    [
        pytest.param(
            1.0000001, 1.0000002, 6, True, does_not_raise(), id="within-6-digit-tol"
        ),
        pytest.param(1.001, 1.002, 2, True, does_not_raise(), id="within-2-digit-tol"),
        pytest.param(
            1.001, 1.002, 4, False, does_not_raise(), id="outside-4-digit-tol"
        ),
        pytest.param(
            [1.0, 2.0], [1.0, 2.0], 6, True, does_not_raise(), id="equal-arrays"
        ),
        pytest.param(
            [1.0, 2.0], [1.0, 2.1], 6, False, does_not_raise(), id="unequal-arrays"
        ),
        pytest.param(0.0, 0.0, 6, True, does_not_raise(), id="both-zero"),
    ],
)
def test_allclose(a, b, digits, expected, context):
    with context:
        assert allclose(a, b, digits=digits) == expected


def test_allclose_explicit_atol():
    assert allclose(1.0, 1.05, atol=0.1) is True
    assert allclose(1.0, 1.15, atol=0.1) is False


def test_allclose_uses_package_default():
    set_precision(3)  # atol = 5e-4
    assert allclose(1.0000, 1.0004) is True  # within 5e-4
    assert allclose(1.0000, 1.0010) is False  # outside 5e-4
