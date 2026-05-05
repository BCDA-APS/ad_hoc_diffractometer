# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.axes.

Covers:
  - parse_axis()
  - axis_label()
  - axis_from_physical()
  - kappa_axis()
"""

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from helpers import STANDARD_BASIS as _STANDARD_BASIS

from ad_hoc_diffractometer.axes import axis_from_physical
from ad_hoc_diffractometer.axes import axis_label
from ad_hoc_diffractometer.axes import kappa_axis
from ad_hoc_diffractometer.axes import parse_axis
from ad_hoc_diffractometer.constants import XHAT
from ad_hoc_diffractometer.constants import YHAT
from ad_hoc_diffractometer.constants import ZHAT

# ---------------------------------------------------------------------------
# parse_axis()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "label, expected, context",
    [
        # Signed Cartesian labels
        pytest.param("+x", +XHAT, does_not_raise(), id="plus-x"),
        pytest.param("-x", -XHAT, does_not_raise(), id="minus-x"),
        pytest.param("+y", +YHAT, does_not_raise(), id="plus-y"),
        pytest.param("-y", -YHAT, does_not_raise(), id="minus-y"),
        pytest.param("+z", +ZHAT, does_not_raise(), id="plus-z"),
        pytest.param("-z", -ZHAT, does_not_raise(), id="minus-z"),
        # Case insensitive
        pytest.param("+X", +XHAT, does_not_raise(), id="plus-X-uppercase"),
        pytest.param("-Z", -ZHAT, does_not_raise(), id="minus-Z-uppercase"),
        # Physical direction names (require basis)
        pytest.param(
            "vertical",
            +XHAT,
            does_not_raise(),
            id="physical-vertical-no-sign",
        ),
        pytest.param(
            "+vertical",
            +XHAT,
            does_not_raise(),
            id="physical-plus-vertical",
        ),
        pytest.param(
            "-vertical",
            -XHAT,
            does_not_raise(),
            id="physical-minus-vertical",
        ),
        pytest.param(
            "-transverse",
            -ZHAT,
            does_not_raise(),
            id="physical-minus-transverse",
        ),
        pytest.param(
            "+longitudinal",
            +YHAT,
            does_not_raise(),
            id="physical-plus-longitudinal",
        ),
        # Failures
        pytest.param(
            "vertical",
            None,
            pytest.raises(ValueError, match=re.escape("basis dict is required")),
            id="invalid-physical-name-no-basis",
        ),
        pytest.param(
            "+bogus",
            None,
            pytest.raises(ValueError, match=re.escape("basis dict is required")),
            id="invalid-unknown-label-no-basis",
        ),
    ],
)
def test_parse_axis(label, expected, context):
    is_cartesian = label.strip().lower().lstrip("+-") in ("x", "y", "z")
    is_failure = expected is None
    if is_cartesian or is_failure:
        basis = None
    else:
        basis = _STANDARD_BASIS
    with context:
        result = parse_axis(label, basis=basis)
        np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize(
    "label, basis, expected, context",
    [
        pytest.param(
            "+x", None, +XHAT, does_not_raise(), id="cartesian-no-basis-needed"
        ),
        pytest.param(
            "vertical",
            _STANDARD_BASIS,
            +XHAT,
            does_not_raise(),
            id="physical-with-basis",
        ),
        pytest.param(
            "-transverse",
            _STANDARD_BASIS,
            -ZHAT,
            does_not_raise(),
            id="physical-minus-with-basis",
        ),
        pytest.param(
            "+bogus",
            _STANDARD_BASIS,
            None,
            pytest.raises(ValueError, match=re.escape("not found in basis dict")),
            id="invalid-direction-not-in-basis",
        ),
    ],
)
def test_parse_axis_with_basis(label, basis, expected, context):
    with context:
        result = parse_axis(label, basis=basis)
        np.testing.assert_array_equal(result, expected)


# ---------------------------------------------------------------------------
# axis_label()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "vector, expected_label, context",
    [
        pytest.param(+XHAT, "+x", does_not_raise(), id="plus-xhat"),
        pytest.param(-XHAT, "-x", does_not_raise(), id="minus-xhat"),
        pytest.param(+YHAT, "+y", does_not_raise(), id="plus-yhat"),
        pytest.param(-YHAT, "-y", does_not_raise(), id="minus-yhat"),
        pytest.param(+ZHAT, "+z", does_not_raise(), id="plus-zhat"),
        pytest.param(-ZHAT, "-z", does_not_raise(), id="minus-zhat"),
        pytest.param(
            np.array([0.5, 0.5, 0.0]),
            "[0.5, 0.5, 0]",
            does_not_raise(),
            id="non-standard-vector-numeric",
        ),
    ],
)
def test_axis_label(vector, expected_label, context):
    with context:
        assert axis_label(vector) == expected_label


# ---------------------------------------------------------------------------
# axis_from_physical()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "direction, sign, expected, context",
    [
        pytest.param("vertical", "+", +XHAT, does_not_raise(), id="plus-vertical"),
        pytest.param("vertical", "-", -XHAT, does_not_raise(), id="minus-vertical"),
        pytest.param("transverse", "-", -ZHAT, does_not_raise(), id="minus-transverse"),
        pytest.param(
            "longitudinal", "+", +YHAT, does_not_raise(), id="plus-longitudinal"
        ),
        pytest.param(
            "bogus",
            "+",
            None,
            pytest.raises(ValueError, match=re.escape("not found in basis dict")),
            id="invalid-unknown-direction",
        ),
    ],
)
def test_axis_from_physical(direction, sign, expected, context):
    with context:
        result = axis_from_physical(direction, sign, _STANDARD_BASIS)
        np.testing.assert_array_equal(result, expected)


# ---------------------------------------------------------------------------
# kappa_axis()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "alpha_deg, expected, context",
    [
        pytest.param(
            0.0,
            None,
            pytest.raises(ValueError, match=re.escape("must be in (0, 90)")),
            id="invalid-alpha-zero",
        ),
        pytest.param(
            90.0,
            None,
            pytest.raises(ValueError, match=re.escape("must be in (0, 90)")),
            id="invalid-alpha-90",
        ),
        pytest.param(
            -10.0,
            None,
            pytest.raises(ValueError, match=re.escape("must be in (0, 90)")),
            id="invalid-alpha-negative",
        ),
        pytest.param(
            50.0,
            np.cos(np.deg2rad(50)) * XHAT + np.sin(np.deg2rad(50)) * ZHAT,
            does_not_raise(),
            id="alpha-50-default",
        ),
        pytest.param(
            1.0,
            np.cos(np.deg2rad(1)) * XHAT + np.sin(np.deg2rad(1)) * ZHAT,
            does_not_raise(),
            id="alpha-1-near-vertical",
        ),
        pytest.param(
            89.0,
            np.cos(np.deg2rad(89)) * XHAT + np.sin(np.deg2rad(89)) * ZHAT,
            does_not_raise(),
            id="alpha-89-near-transverse",
        ),
        pytest.param(
            45.0,
            np.array([np.sqrt(2) / 2, 0.0, np.sqrt(2) / 2]),
            does_not_raise(),
            id="alpha-45-equal-components",
        ),
    ],
)
def test_kappa_axis(alpha_deg, expected, context):
    with context:
        result = kappa_axis(alpha_deg)
        np.testing.assert_allclose(result, expected, atol=1e-12)


@pytest.mark.parametrize(
    "alpha_deg, context",
    [
        pytest.param(50.0, does_not_raise(), id="kappa-axis-is-unit-vector-50"),
        pytest.param(30.0, does_not_raise(), id="kappa-axis-is-unit-vector-30"),
        pytest.param(70.0, does_not_raise(), id="kappa-axis-is-unit-vector-70"),
    ],
)
def test_kappa_axis_is_unit_vector(alpha_deg, context):
    with context:
        ax = kappa_axis(alpha_deg)
        np.testing.assert_allclose(np.linalg.norm(ax), 1.0, atol=1e-12)


def test_kappa_axis_bad_basis():
    with pytest.raises(ValueError, match=re.escape("missing: 'transverse'")):
        kappa_axis(50.0, basis={"vertical": XHAT, "longitudinal": YHAT})


def test_kappa_axis_with_explicit_basis():
    """``kappa_axis`` honors an explicit basis dict (covers the
    happy path that no longer goes through ``presets.py`` after
    issue #241).
    """
    basis = {"vertical": XHAT, "longitudinal": YHAT, "transverse": ZHAT}
    ax = kappa_axis(50.0, basis=basis)
    expected = np.cos(np.deg2rad(50.0)) * XHAT + np.sin(np.deg2rad(50.0)) * ZHAT
    np.testing.assert_allclose(ax, expected, atol=1e-12)
