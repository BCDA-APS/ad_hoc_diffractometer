"""
Unit tests for ad_hoc_diffractometer.stage.

Covers:
  - Stage construction
  - Stage.rotation_matrix()
  - Stage.limits (property, validation)
  - Stage.in_limits()
"""

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from helpers import Rx
from helpers import Ry
from helpers import Rz

from ad_hoc_diffractometer import XHAT
from ad_hoc_diffractometer import YHAT
from ad_hoc_diffractometer import ZHAT
from ad_hoc_diffractometer import Stage

# ---------------------------------------------------------------------------
# Stage construction
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, axis, parent, role, angle, context",
    [
        pytest.param(
            "mu", XHAT, None, "sample", 0.0, does_not_raise(), id="sample-floor-stage"
        ),
        pytest.param(
            "eta",
            -ZHAT,
            "mu",
            "sample",
            45.0,
            does_not_raise(),
            id="sample-stacked-stage",
        ),
        pytest.param(
            "delta",
            -ZHAT,
            "nu",
            "detector",
            0.0,
            does_not_raise(),
            id="detector-stacked-stage",
        ),
        pytest.param(
            "phi",
            -ZHAT,
            None,
            "sample",
            0.0,
            does_not_raise(),
            id="sample-floor-stage-neg-z",
        ),
    ],
)
def test_stage_construction(name, axis, parent, role, angle, context):
    with context:
        s = Stage(name, axis, parent=parent, role=role, angle=angle)
        assert s.name == name
        np.testing.assert_array_equal(s.axis, np.asarray(axis, dtype=float))
        assert s.parent == parent
        assert s.role == role
        assert s.angle == angle


# ---------------------------------------------------------------------------
# Stage.rotation_matrix()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "axis, angle_deg, expected_R, context",
    [
        pytest.param(
            XHAT, 0.0, np.eye(3), does_not_raise(), id="zero-angle-is-identity"
        ),
        pytest.param(XHAT, 90.0, Rx(90), does_not_raise(), id="90deg-about-x"),
        pytest.param(-ZHAT, 30.0, Rz(-30), does_not_raise(), id="left-handed-z-30deg"),
        pytest.param(YHAT, -45.0, Ry(-45), does_not_raise(), id="neg45deg-about-y"),
    ],
)
def test_stage_rotation_matrix(axis, angle_deg, expected_R, context):
    with context:
        s = Stage("test", axis, angle=angle_deg)
        np.testing.assert_allclose(s.rotation_matrix(), expected_R, atol=1e-10)


# ---------------------------------------------------------------------------
# Stage.limits — validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "limits, context",
    [
        pytest.param((-180.0, 180.0), does_not_raise(), id="default-limits"),
        pytest.param((0.0, 360.0), does_not_raise(), id="zero-to-360"),
        pytest.param((-90.0, 90.0), does_not_raise(), id="symmetric-90"),
        pytest.param(
            (0.0, 0.0),
            pytest.raises(ValueError, match=re.escape("must be less than max")),
            id="invalid-equal-limits",
        ),
        pytest.param(
            (10.0, -10.0),
            pytest.raises(ValueError, match=re.escape("must be less than max")),
            id="invalid-reversed-limits",
        ),
        pytest.param(
            (180.0, -180.0),
            pytest.raises(ValueError, match=re.escape("must be less than max")),
            id="invalid-reversed-full-range",
        ),
    ],
)
def test_stage_limits_validation(limits, context):
    with context:
        s = Stage("mu", XHAT, limits=limits)
        assert s.limits == limits


def test_stage_default_limits():
    """Default limits are (-180, 180) for all roles."""
    for role in ("sample", "detector"):
        s = Stage("test", XHAT, role=role)
        assert s.limits == (-180.0, 180.0)


def test_stage_limits_setter_updates():
    """Limits can be updated after construction."""
    s = Stage("mu", XHAT)
    s.limits = (-90.0, 90.0)
    assert s.limits == (-90.0, 90.0)
    assert s.in_limits(89.0) is True
    assert s.in_limits(91.0) is False


def test_stage_limits_setter_invalid_after_construction():
    """Setting invalid limits after construction raises ValueError."""
    s = Stage("mu", XHAT)
    with pytest.raises(ValueError, match=re.escape("must be less than max")):
        s.limits = (10.0, 5.0)


# ---------------------------------------------------------------------------
# Stage.in_limits()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "limits, angle, expected, context",
    [
        pytest.param((-180.0, 180.0), 0.0, True, does_not_raise(), id="zero-in-range"),
        pytest.param(
            (-180.0, 180.0), 180.0, True, does_not_raise(), id="at-max-boundary"
        ),
        pytest.param(
            (-180.0, 180.0), -180.0, True, does_not_raise(), id="at-min-boundary"
        ),
        pytest.param((-180.0, 180.0), 90.0, True, does_not_raise(), id="inside-range"),
        pytest.param(
            (-180.0, 180.0), 180.1, False, does_not_raise(), id="just-above-max"
        ),
        pytest.param(
            (-180.0, 180.0), -180.1, False, does_not_raise(), id="just-below-min"
        ),
        pytest.param((0.0, 90.0), 45.0, True, does_not_raise(), id="asymmetric-inside"),
        pytest.param(
            (0.0, 90.0), -1.0, False, does_not_raise(), id="asymmetric-below-min"
        ),
        pytest.param(
            (0.0, 90.0), 91.0, False, does_not_raise(), id="asymmetric-above-max"
        ),
        pytest.param(
            (-180.0, 180.0), 360.0, False, does_not_raise(), id="outside-full-range"
        ),
    ],
)
def test_stage_in_limits(limits, angle, expected, context):
    with context:
        s = Stage("mu", XHAT, limits=limits)
        assert s.in_limits(angle) == expected


def test_stage_repr_with_parent():
    """Stage.__repr__ includes the parent name when parent is set."""
    s = Stage("chi", XHAT, role="sample", parent="omega")
    r = repr(s)
    assert "chi" in r
    assert "omega" in r
