# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
Unit tests for ad_hoc_diffractometer.stage.

Covers:
  - Stage construction
  - Stage.rotation_matrix()
  - Stage.limits (property, validation)
  - Stage.in_limits()
  - Stage.to_dict() / from_dict(): name, axis, role, parent, angle, limits
"""

import json
import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from helpers import Rx
from helpers import Ry
from helpers import Rz
from helpers import fourcv

from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer.constants import XHAT
from ad_hoc_diffractometer.constants import YHAT
from ad_hoc_diffractometer.constants import ZHAT
from ad_hoc_diffractometer.stage import Stage

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


# ---------------------------------------------------------------------------
# Stage.to_dict() / from_dict()
# ---------------------------------------------------------------------------

_OMEGA = Stage(
    "omega",
    np.array([0.0, 0.0, -1.0]),
    role="sample",
    parent=None,
    limits=(-180.0, 180.0),
)
_OMEGA.angle = 20.97
_CHI = Stage(
    "chi",
    np.array([0.0, 1.0, 0.0]),
    role="sample",
    parent="omega",
    limits=(-90.0, 90.0),
)
_CHI.angle = 45.0


def test_stage_to_dict_structure():
    """to_dict() returns a JSON-serialisable dict with all required keys."""
    d = _OMEGA.to_dict()
    assert isinstance(d, dict)
    assert {"name", "axis", "role", "parent", "angle", "limits"} <= set(d.keys())
    assert json.dumps(d)  # must not raise


@pytest.mark.parametrize(
    "key, expected, context",
    [
        pytest.param("name", "omega", does_not_raise(), id="name"),
        pytest.param("role", "sample", does_not_raise(), id="role"),
        pytest.param("parent", None, does_not_raise(), id="parent-none"),
        pytest.param("angle", 20.97, does_not_raise(), id="angle"),
        pytest.param("limits", [-180.0, 180.0], does_not_raise(), id="limits"),
        pytest.param("axis", [0.0, 0.0, -1.0], does_not_raise(), id="axis"),
    ],
)
def test_stage_to_dict_values(key, expected, context):
    """to_dict() stores the correct value for each field."""
    with context:
        assert _OMEGA.to_dict()[key] == pytest.approx(expected)


def test_stage_to_dict_parent_name():
    """to_dict() stores the parent name string for a child stage."""
    assert _CHI.to_dict()["parent"] == "omega"


def test_stage_to_dict_axis_is_list():
    """to_dict() encodes the axis as a plain list of length 3."""
    axis = _OMEGA.to_dict()["axis"]
    assert isinstance(axis, list)
    assert len(axis) == 3


@pytest.mark.parametrize(
    "attr, accessor, context",
    [
        pytest.param("name", lambda s: s.name, does_not_raise(), id="name"),
        pytest.param("role", lambda s: s.role, does_not_raise(), id="role"),
        pytest.param("angle", lambda s: s.angle, does_not_raise(), id="angle"),
        pytest.param("limits", lambda s: s.limits, does_not_raise(), id="limits"),
    ],
)
def test_stage_from_dict_roundtrip(attr, accessor, context):
    """from_dict(to_dict()) recovers each scalar/string attribute."""
    with context:
        restored = Stage.from_dict(_OMEGA.to_dict())
        assert accessor(restored) == pytest.approx(accessor(_OMEGA))


def test_stage_from_dict_roundtrip_axis():
    """from_dict(to_dict()) recovers the axis vector."""
    np.testing.assert_allclose(Stage.from_dict(_OMEGA.to_dict()).axis, _OMEGA.axis)


def test_stage_from_dict_roundtrip_parent():
    """from_dict(to_dict()) recovers the parent name for a child stage."""
    assert Stage.from_dict(_CHI.to_dict()).parent == _CHI.parent


def test_geometry_to_dict_stage_delegation():
    """Each stage entry in geometry.to_dict() equals Stage.to_dict()."""
    g = fourcv()
    g.set_angle("omega", 20.97)
    geo_stages = {sd["name"]: sd for sd in g.to_dict()["stages"]}
    for name, stage_obj in g._stages.items():
        assert geo_stages[name] == stage_obj.to_dict()


def test_geometry_from_dict_stage_delegation():
    """Stages in AdHocDiffractometer.from_dict() match Stage.from_dict()."""
    g = fourcv()
    g.set_angle("chi", 45.0)
    g2 = AdHocDiffractometer.from_dict(g.to_dict())
    for name in g._stages:
        assert g2._stages[name].angle == pytest.approx(g._stages[name].angle)
        np.testing.assert_allclose(g2._stages[name].axis, g._stages[name].axis)
        assert g2._stages[name].role == g._stages[name].role
        assert g2._stages[name].limits == pytest.approx(g._stages[name].limits)
